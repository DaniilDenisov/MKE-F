import numpy as np


def forces(values):
    f = np.asarray(values, dtype=float)
    if f.shape != (6,) or not np.isfinite(f).all():
        raise ValueError("Six finite profile-axis forces required")
    return f


def weights(bary):
    a,b,c = np.moveaxis(bary, -1, 0)
    return np.stack((a*(2*a-1), b*(2*b-1), c*(2*c-1), 4*a*b, 4*b*c, 4*c*a), axis=-1)


def barycentric(vertices, points):
    u, v = vertices[...,1,:]-vertices[...,0,:], vertices[...,2,:]-vertices[...,0,:]
    p = points-vertices[...,0,:]
    def cross(a, b):
        return a[...,0]*b[...,1]-a[...,1]*b[...,0]
    det = cross(u, v)
    b, c = cross(p,v)/det, cross(u,p)/det
    return np.stack((1-b-c,b,c), axis=-1)


def nodal(kernel, load, ids):
    n, vy, vz, torque, my, mz = load
    mesh = kernel.mesh
    area, iy, iz = mesh.profile.properties
    coords = mesh.nodes[mesh.triangles[ids][...,[0,1,2,5,3,4]]]
    sigma = n/area + my*coords[...,1]/iy - mz*coords[...,0]/iz
    tau = np.einsum("...nfk,f->...nk",kernel.units[ids],[vy,vz,torque])
    return np.concatenate((sigma[...,None], tau), axis=-1)


def stress_result(stress, coords, element):
    sigma, ty, tz = map(float, stress)
    value = np.hypot(sigma, np.sqrt(3)*np.hypot(ty,tz))
    if not np.isfinite(value):
        raise ArithmeticError("Stress overflow")
    return dict(status="ok", y=float(coords[0]), z=float(coords[1]), element=int(element),
                sigmaX=sigma, tauXY=ty, tauXZ=tz, value=float(value))


def point(kernel, load, y, z, element=None):
    """The lowest containing element id owns a shared edge unless a side is specified."""
    f = forces(load)
    mesh = kernel.mesh
    profile = mesh.profile
    if not np.isfinite([y,z]).all():
        raise ValueError("Finite coordinates required")
    if profile.round and profile.t and np.hypot(y,z) < profile.b/2-profile.t:
        return dict(status="unavailable", reason="Inside hole")
    if profile.inner_radius:
        offsets = np.abs([y,z])-np.array([profile.h,profile.b])/2+profile.t+profile.inner_radius
        if np.linalg.norm(np.maximum(offsets,0)) < profile.inner_radius:
            return dict(status="unavailable",reason="Inside hole")
    if profile.sharp and np.any(f[1:4]) and np.allclose(
            np.abs([y,z]), np.array([profile.h,profile.b])/2-profile.t,
            rtol=0, atol=1e-12*min(profile.h,profile.b)):
        return dict(status="unavailable", singularPeak=True)
    if element is not None and (not isinstance(element, (int, np.integer)) or
                                not 0 <= element < len(mesh.triangles)):
        raise ValueError("Invalid element id")
    ids = np.arange(len(mesh.triangles)) if element is None else np.array([element])
    bary = barycentric(mesh.nodes[mesh.triangles[ids,:3]], np.array([y,z]))
    hits = np.flatnonzero((bary >= -1e-12).all(axis=1))
    if not hits.size:
        return dict(status="unavailable", reason="Outside selected mesh material")
    i = hits[0]
    return stress_result(weights(bary[i]) @ nodal(kernel,f,ids[i]), [y,z], ids[i])

