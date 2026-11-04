"""In-process refinement session. One owner per Engine, geometry in consistent units."""
from dataclasses import dataclass
from time import monotonic

import numpy as np
import psutil
from matplotlib.tri import Triangulation

from . import adapter
from .peak import Stopped, maximum
from .profile import positive, validate_nu
from .recovery import barycentric, forces, weights


@dataclass(frozen=True)
class Options:
    segments: int = 64
    area_fraction: float = 1/32
    levels: int = 5
    max_elements: int = 50000
    max_bytes: int = 1536*1024**2
    seconds: float = 120
    max_cells: int = 1000000
    search_tolerance: float = .001

    def __post_init__(self):
        for n in (self.segments,self.levels,self.max_elements,self.max_bytes,self.max_cells):
            if type(n) is not int or n < 1:
                raise ValueError("Positive integer limits required")
        if (self.segments < 8 or self.segments % 4 or self.levels > 20 or
                not positive(self.area_fraction) or not positive(self.seconds) or
                not 0 < self.search_tolerance <= .005):
            raise ValueError("Invalid mesh or search options")


def relative(a,b):
    return abs(a-b)/max(abs(a),abs(b)) if a or b else 0.


# Seven-point degree-five Dunavant quadrature on the coarser mesh.
QUAD = [([1/3]*3,.225)]
for a,b,w in ((.059715871789770,.470142064105115,.132394152788506),
              (.797426985353087,.101286507323456,.125939180544827)):
    QUAD.extend((np.roll([a,b,b],i),w) for i in range(3))


def field_changes(old, new, check):
    mesh = new.mesh
    scale = max(mesh.profile.h,mesh.profile.b)
    nodes = mesh.nodes/scale
    finder = Triangulation(*nodes.T, mesh.triangles[:,:3]).get_trifinder()
    vertices = old.mesh.vertices/scale
    units = np.array([scale**2,scale**2,scale**3])[:,None]
    areas = abs(np.linalg.det(vertices[:,1:]-vertices[:,:1]))/2
    errors, norms = np.zeros(3), np.zeros(3)
    for bary,w in QUAD:
        check()
        xy = np.einsum("n,enc->ec",bary,vertices)
        ids = finder(*xy.T)
        keep = ids >= 0  # compare on the intersection of polygonal material domains
        if not np.any(keep):
            raise ArithmeticError("No common material points")
        ids = ids[keep]
        local = barycentric(nodes[mesh.triangles[ids,:3]], xy[keep])
        a = np.einsum("n,enfc->efc", weights(bary), old.units[keep])*units
        b = np.einsum("en,enfc->efc", weights(local), new.units[ids])*units
        errors += w*np.einsum("e,efc->f",areas[keep],(a-b)**2)
        norms += w*np.einsum("e,efc->f",areas[keep],b*b)
    return np.sqrt(errors/norms).tolist()


def geometry_error(mesh):
    vertices = mesh.vertices
    areas = abs(np.linalg.det(vertices[:,1:]-vertices[:,:1]))/2
    # Exact second moments of affine coordinates over a triangle.
    iz,iy = areas @ (vertices.sum(axis=1)**2+(vertices**2).sum(axis=1))/12
    return float(max(abs(np.array([areas.sum(),iy,iz])/mesh.profile.properties-1)))


class Engine:
    def __init__(self):
        self.meshes = {}
        self.fields = {}
        self.comparisons = {}

    @property
    def resident_bytes(self):
        return (sum(m.nbytes for m in self.meshes.values()) +
                sum(k.units.nbytes for k in self.fields.values()) +
                512*(len(self.meshes)+len(self.fields)+len(self.comparisons)))

    def solve(self, profile, nu, load, options=Options(), *, J_model=None, cancelled=lambda: False):
        f = forces(load)
        def evaluate(k, check):
            return maximum(k,f,tolerance=options.search_tolerance,max_cells=options.max_cells,check=check)
        return self.refine(profile,nu,options,evaluate,J_model,cancelled)

    def refine(self, profile, nu, o, evaluate, J_model, cancelled):
        nu = validate_nu(nu)
        deadline = monotonic()+o.seconds
        process = psutil.Process()
        def check():
            if cancelled():
                raise Stopped("cancelled","Cancelled by caller")
            if monotonic() >= deadline or max(self.resident_bytes,process.memory_info().rss) > o.max_bytes:
                raise Stopped("resource_limit","Time or resident cache limit")
        history = []
        result = dict(status="nonconverged",fieldStatus="nonconverged",JStatus="nonconverged",
                      maximumStatus="nonconverged",history=history)
        previous = previous_key = None
        changes = []
        try:
            for level in range(o.levels):
                check()
                curved = profile.round or profile.inner_radius>0
                segments = o.segments * (2**level if curved else 1)
                fraction = o.area_fraction/4**level
                mesh_key = (profile,segments,fraction,level,adapter.GENERATOR)
                key = (mesh_key,nu,adapter.SOLVER)
                mesh = self.meshes.get(mesh_key)
                # area_fraction is relative to material area. Also budget boundary
                # segments and thin-wall triangles before invoking the generator.
                estimate = max(1/fraction, segments*(1+bool(profile.t)) if curved else 4,
                               max(profile.b,profile.h)/min(profile.b,profile.h),
                               max(profile.h,profile.b)/(profile.t or min(profile.h,profile.b)))
                if mesh is None:
                    if estimate > o.max_elements or process.memory_info().rss+estimate*12000 > o.max_bytes:
                        raise Stopped("resource_limit","Mesh generation budget")
                    mesh = adapter.generate(profile,segments,profile.properties[0]*fraction)
                count = len(mesh.triangles)
                if count > o.max_elements:
                    raise Stopped("resource_limit","Triangle limit")
                field = self.fields.get(key)
                # Reserve space for the mesh, sparse factorization and element objects
                # on a miss, or for quadrature and search scratch space on a hit.
                reserve = count*(1600 if field is not None else 40000)
                if max(process.memory_info().rss,self.resident_bytes)+reserve > o.max_bytes:
                    raise Stopped("resource_limit","Working memory estimate")
                self.meshes[mesh_key] = mesh
                if field is None:
                    field = adapter.solve_fields(mesh,nu,check)
                    self.fields[key] = field
                check()
                geom = geometry_error(mesh)
                point_changes = j_change = None
                if previous is not None:
                    pair = (previous_key,key)
                    if pair not in self.comparisons:
                        self.comparisons[pair] = field_changes(previous,field,check)
                    point_changes = self.comparisons[pair]
                    j_change = relative(previous.J_FE,field.J_FE)
                field_ok = point_changes is not None and max(point_changes)<=.005 and geom<=1e-4
                j_ok = j_change is not None and j_change<=.001 and geom<=1e-4
                peak = evaluate(field,check)
                value = peak.get("value")
                if value is not None and history and history[-1]["peak"] is not None:
                    changes.append(relative(value,history[-1]["peak"]))
                peak_status = peak["status"]
                analytic = peak.get("analytic",profile.sharp and value is not None)
                if peak_status == "ok" and not analytic:
                    if not field_ok or len(changes)<2 or max(changes[-2:])>.005:
                        peak_status = "nonconverged"
                history.append(dict(level=level, elements=count, segments=segments, geometry_error=geom,
                                    field_changes=point_changes,J_change=j_change,J_FE=field.J_FE,peak=value,
                                    peak_change=changes[-1] if changes else None))
                result.update(fieldStatus="ok" if field_ok else "nonconverged",
                              JStatus="ok" if j_ok else "nonconverged",maximumStatus=peak_status,
                              maximum=peak,kernel=field,**field.torsion_report(J_model))
                check()
                if field_ok and j_ok and peak_status in ("ok","unavailable"):
                    result["status"] = peak_status
                    break
                previous,previous_key = field,key
        except Stopped as stop:
            result.update(status=stop.status,reason=str(stop))
        except (ArithmeticError, RuntimeError) as err:
            result.update(status="nonconverged",reason=str(err))
        return result

