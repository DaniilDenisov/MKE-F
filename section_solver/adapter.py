"""The sectionproperties boundary. Library x,y correspond to profile y,z."""
from dataclasses import dataclass
from importlib.metadata import version

import numpy as np
from sectionproperties.analysis import Section
from sectionproperties.pre import Material
from sectionproperties.pre.geometry import Geometry
from shapely.geometry import Polygon

from .profile import Profile, positive, validate_nu

GENERATOR = ("cytriangle", version("cytriangle"), version("sectionproperties"))
SOLVER = ("direct", version("sectionproperties"))


@dataclass
class Mesh:
    profile: Profile
    segments: int
    nodes: np.ndarray
    triangles: np.ndarray  # library ordering, including all six nodes

    @property
    def vertices(self):
        return self.nodes[self.triangles[:, :3]]

    @property
    def nbytes(self):
        return self.nodes.nbytes + self.triangles.nbytes


def generate(profile, segments, area):
    if not positive(area):
        raise ValueError("Maximum triangle area must be finite and positive")
    scale = max(profile.h,profile.b)
    if not 1e-3 <= scale <= 1e3:
        # The geometry library rounds contour coordinates at fixed precision.
        unit = Profile(profile.kind,profile.h/scale,profile.b/scale,profile.t/scale,profile.inner_radius/scale)
        mesh = generate(unit,segments,area/scale**2)
        return Mesh(profile,segments,mesh.nodes*scale,mesh.triangles)
    polygon = profile.polygon(segments)
    if profile.sharp:
        # Graded straight boundary facets resolve the integrable corner field.
        # Keep the exact corner and the entire material domain in the search.
        ring = np.array(polygon.interiors[0].coords)
        spacing = min(profile.t, np.sqrt(area))
        smallest = spacing * min(1., area/profile.t**2)
        inner = []
        for a,b in zip(ring,ring[1:]):
            length = np.linalg.norm(b-a)
            d = min(length/3,profile.t)
            offsets = []
            while d >= smallest:
                offsets.extend((d/length,1-d/length))
                d /= 2
            inner.append(a)
            inner.extend(a+u*(b-a) for u in sorted(offsets))
        polygon = Polygon(polygon.exterior, [inner])
    geometry = Geometry(polygon)
    raw = geometry.create_mesh(mesh_sizes=area).mesh
    return Mesh(profile, segments, np.asarray(raw["vertices"]), np.asarray(raw["triangles"]))


@dataclass
class Kernel:
    mesh: Mesh
    nu: float
    J_FE: float
    units: np.ndarray  # (element, local node, Vy/Vz/T, tauXY/tauXZ)

    def torsion_report(self, J_model=None):
        if J_model is not None and not positive(J_model):
            raise ValueError("J_model must be finite and positive")
        return dict(J_FE=self.J_FE, J_model=J_model,
                    relative_difference=None if J_model is None else (self.J_FE-J_model)/J_model)


def solve_fields(mesh, nu, check=lambda: None):
    nu = validate_nu(nu)
    p = mesh.profile
    scale = max(p.h,p.b)
    if not 1e-3 <= scale <= 1e3:
        unit = Profile(p.kind,p.h/scale,p.b/scale,p.t/scale,p.inner_radius/scale)
        normalized = Mesh(unit,mesh.segments,mesh.nodes/scale,mesh.triangles)
        field = solve_fields(normalized,nu,check)
        field.units[:,:,:2] /= scale**2
        field.units[:,:,2] /= scale**3
        return Kernel(mesh,nu,field.J_FE*scale**4,field.units)
    material = Material("profile", 1., nu, 1., 1., "gray")
    geometry = Geometry(mesh.profile.polygon(mesh.segments), material=material)
    geometry.mesh = dict(vertices=mesh.nodes, triangles=mesh.triangles,
                         triangle_attributes=np.zeros((len(mesh.triangles), 1), dtype=int))
    section = Section(geometry)
    section.calculate_geometric_properties()
    check()
    section.calculate_warping_properties(solver_type="direct")
    check()
    return _extract(section, mesh, nu, check)


def _arguments(section):
    cx, cy = section.get_c()
    ixx, iyy, ixy = section.get_eic()
    i11, i22 = section.get_eip()
    return dict(n=0, mxx=0, myy=0, m11=0, m22=0, vx=1, vy=1, mzz=1,
                ea=section.get_ea(), cx=cx, cy=cy, ixx=ixx, iyy=iyy, ixy=ixy,
                i11=i11, i22=i22, phi=section.get_phi(), j=section.get_ej(),
                nu=section.get_nu_eff(), delta_s=section.section_props.delta_s)


def _local_args(section, element):
    return {name: getattr(section.section_props, name)[element.node_ids]
            for name in ("omega", "psi_shear", "phi_shear")}


def _extract(section, mesh, nu, check):
    common = _arguments(section)
    units = np.empty((len(mesh.triangles), 6, 3, 2))
    for i, element in enumerate(section.elements):
        if i % 64 == 0:
            check()
        s = element.element_stress(**common, **_local_args(section, element))
        for field, start in enumerate((7, 9, 5)):
            units[i, :, field, :] = np.column_stack(s[start:start+2])
    if not np.isfinite(units).all() or not positive(common["j"]):
        raise ArithmeticError("Non-finite field or torsion constant")
    return Kernel(mesh, nu, common["j"], units)

