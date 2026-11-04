import numpy as np

from .recovery import forces, nodal, stress_result, weights


class Stopped(Exception):
    def __init__(self, status, reason):
        self.status = status
        super().__init__(reason)


ROOT = np.eye(3)
MIDS = (ROOT + np.roll(ROOT,-1,axis=0))/2
SAMPLES = np.vstack((ROOT,MIDS,[1/3]*3))
SCALE = np.array([1,np.sqrt(3),np.sqrt(3)])


def search(vertices, values, tolerance=.001, max_cells=1000000, check=lambda: None):
    """Bound the stress norm on a single straight-sided quadratic mesh."""
    if not 0 < tolerance <= .005 or max_cells < 1:
        raise ValueError("Invalid search tolerance or cell limit")
    count = 0
    best = None

    def inspect(ids, corners):
        nonlocal best, count
        count += len(ids)
        check()
        if count > max_cells:
            raise Stopped("resource_limit","Maximum search cell limit")
        bary = SAMPLES @ corners
        samples = weights(bary) @ values[ids]
        vectors = samples*SCALE
        norms = np.hypot.reduce(vectors,axis=-1)
        if not np.isfinite(norms).all():
            raise ArithmeticError("Non-finite search norm")
        i,j = np.unravel_index(norms.argmax(),norms.shape)
        if best is None or norms[i,j] > best["value"]:
            best = stress_result(samples[i,j], bary[i,j] @ vertices[ids[i]], ids[i])
        controls = 2*vectors[:,3:6] - (vectors[:,:3]+np.roll(vectors[:,:3],-1,axis=1))/2
        return np.maximum(norms[:,:3].max(axis=1),np.hypot.reduce(controls,axis=-1).max(axis=1))

    bounds = np.empty(len(values))
    for start in range(0,len(values),256):
        ids = np.arange(start,min(start+256,len(values)))
        bounds[ids] = inspect(ids, np.broadcast_to(ROOT,(len(ids),3,3)))
    upper = best["value"]
    for i, bound in enumerate(bounds):
        stack = [(ROOT,bound)]
        while stack:
            corners, bound = stack.pop()
            if bound <= best["value"]*(1+tolerance):
                upper = max(upper,bound)
                continue
            a,b,c = corners
            ab,bc,ca = MIDS @ corners
            children = np.array([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]])
            stack.extend(zip(children, inspect(np.full(4,i),children)))
    check()
    upper = max(best["value"],upper)
    return dict(**best, upper=float(upper), relative_gap=upper/best["value"]-1 if best["value"] else 0,
                cells=count)


def maximum(kernel, load, **options):
    f = forces(load)
    p = kernel.mesh.profile
    if p.sharp and not np.any(f[1:4]):
        area,iy,iz = p.properties
        corners = np.array([[y*p.h/2,z*p.b/2] for y in (-1,1) for z in (-1,1)])
        normal = f[0]/area + f[4]*corners[:,1]/iy-f[5]*corners[:,0]/iz
        i = abs(normal).argmax()
        return dict(status="ok", value=abs(normal[i]), sigmaX=normal[i], tauXY=0., tauXZ=0.,
                    y=corners[i,0], z=corners[i,1], upper=abs(normal[i]), relative_gap=0., cells=0)
    peak = search(kernel.mesh.vertices, nodal(kernel,f,slice(None)), **options)
    if p.sharp:
        return dict(status="unavailable", singularPeak=True, mesh_peak=peak)
    return peak

