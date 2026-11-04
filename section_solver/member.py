from collections import deque
from dataclasses import replace

import numpy as np
from numpy.polynomial import polynomial as poly

from .engine import Options, relative
from .peak import Stopped, maximum
from .profile import positive


def coefficients(rows):
    if len(rows) != 6:
        raise ValueError("Six force polynomials required")
    c = np.zeros((6,4))
    for i,row in enumerate(rows):
        a = np.asarray(row,dtype=float)
        if a.ndim != 1 or not 1 <= a.size <= 4 or not np.isfinite(a).all():
            raise ValueError("Finite polynomial of degree at most three required")
        c[i,:a.size] = a
    return c


def planar(axial, shear, moment, angle=0):
    """Convert 2D local force polynomials once. 3D coefficients already use profile axes."""
    cs = {0:(1,0),90:(0,1),180:(-1,0),270:(0,-1)}
    if angle not in cs:
        raise ValueError("Unsupported profile angle")
    c,s = cs[angle]
    n,v,_,_,_,m = coefficients([axial,shear,[0],[0],[0],moment])
    return np.array([n,c*v,-s*v,np.zeros(4),s*m,c*m])


def stationary(rows):
    xs = {0.,1.}
    for row in rows:
        for root in poly.polyroots(poly.polyder(row)):
            if np.isreal(root) and 0 < root.real < 1:
                xs.add(float(root.real))
    return sorted(xs)


def truss(axial, area):
    if not positive(area):
        raise ValueError("Finite positive truss area required")
    row = coefficients([axial]+[[0]]*5)[0]
    samples = [dict(xi=x,sigmaX=poly.polyval(x,row)/area) for x in stationary([row])]
    if any(not np.isfinite(p["sigmaX"]) for p in samples):
        return dict(status="nonconverged",reason="Truss stress overflow")
    best = max(samples,key=lambda p: abs(p["sigmaX"]))
    return dict(status="ok",maximum=dict(**best,value=abs(best["sigmaX"]),tauXY=0.,tauXZ=0.),
                samples=samples)


def _length_peak(c, at, units, check, limit, tolerance):
    samples = {}
    def sample(x):
        if x not in samples:
            check()
            if len(samples)>=limit:
                raise Stopped("resource_limit","Length sampling limit")
            samples[x] = dict(at(x),xi=x)
        return samples[x]
    for x in sorted(set(stationary(c)).union(np.linspace(0,1,9))):
        sample(x)
    for _ in range(2):
        xs = sorted(samples)
        for a,b in zip(xs,xs[1:]):
            sample((a+b)/2)
    def interval(a,b):
        second = np.maximum(abs(2*c[:,2]+6*a*c[:,3]),abs(2*c[:,2]+6*b*c[:,3]))
        upper = max(samples[a]["upper"],samples[b]["upper"]) + (b-a)**2/8*(units @ second)
        return upper,a,b
    xs = sorted(samples)
    intervals = [interval(a,b) for a,b in zip(xs,xs[1:])]
    while True:
        check()
        best = max(samples.values(),key=lambda p:p["value"])
        upper,a,b = max(intervals)
        if upper <= best["value"]*(1+tolerance):
            return dict(best,upper=upper,relative_gap=upper/best["value"]-1 if best["value"] else 0.,
                        section_upper=best["upper"],samples=[samples[x] for x in sorted(samples)])
        x = (a+b)/2
        sample(x)
        intervals.remove((upper,a,b))
        intervals.extend([interval(a,x),interval(x,b)])


def member(engine, profile, nu, rows, options=Options(), *, max_samples=257,
           length_tolerance=.0025, cancelled=lambda: False, J_model=None, include_mesh_peak=False):
    c = coefficients(rows)
    if type(max_samples) is not int or max_samples < 2 or not 0 < length_tolerance <= .005:
        raise ValueError("Invalid length search options")
    singular = profile.sharp and c[1:4].any()
    if singular and not include_mesh_peak:
        return dict(status="unavailable",maximumStatus="unavailable",singularPeak=True)
    o = replace(options,search_tolerance=min(options.search_tolerance,length_tolerance/4))
    older = deque(maxlen=2)

    def evaluate(kernel,check):
        def on(k, f):
            result = maximum(k,f,tolerance=o.search_tolerance,max_cells=o.max_cells,check=check)
            return result.get("mesh_peak",result)
        cache = {}
        def at(x, k=kernel, store=cache):
            if x not in store:
                store[x] = on(k,poly.polyval(x,c.T))
            return store[x]
        if (singular or profile.inner_radius) and not c[:,2:].any():
            # With affine forces the stress norm at every material point is
            # convex along the span. One of the two ends attains the mesh peak.
            ends = [dict(at(x),xi=x) for x in (0.,1.)]
            best = max(ends,key=lambda p:p["value"])
            upper = max(p["upper"] for p in ends)
            peak = dict(best,section_upper=best["upper"],upper=upper,samples=ends,
                        relative_gap=upper/best["value"]-1 if best["value"] else 0.)
        else:
            units = np.zeros(6)
            for i in np.flatnonzero(c[:,2:].any(axis=1)):
                units[i] = on(kernel,np.eye(6)[i])["upper"]
            peak = _length_peak(c,at,units,check,max_samples,length_tolerance)
        if singular:
            return dict(status="unavailable",singularPeak=True,mesh_peak=peak)
        if not profile.sharp:
            # Adaptive grids may differ: compare the same physical sections on
            # both previous meshes, including newly inserted length samples.
            changes = []
            for row in peak["samples"]:
                values = [fn(row["xi"])["value"] for fn in older]+[row["value"]]
                changes.extend(relative(a,b) for a,b in zip(values,values[1:]))
            peak["sample_mesh_change"] = max(changes,default=None)
            if len(older)<2 or max(changes,default=1.)>.005:
                peak["status"] = "nonconverged"
        older.append(at)
        return peak
    return engine.refine(profile,nu,o,evaluate,J_model,cancelled)

