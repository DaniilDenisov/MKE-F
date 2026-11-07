"""Emit one isolated cold run per profile, then reuse its meshes for six unit loads."""
import json
from dataclasses import asdict
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
from time import perf_counter
from importlib.metadata import version

import numpy as np

from section_solver import Profile
from section_solver.engine import Engine, Options, QUAD
from section_solver.recovery import weights

OUT = Path(os.environ.get("SECTION_REPORT_DIR","output/section-solver"))
SECTIONS = {
    "circle": {"d":2}, "circular_tube":{"d":2,"t":.3}, "rectangle":{"h":3,"b":1},
    "rectangular_tube":{"h":3,"b":1,"t":.15}, "square_tube":{"b":1,"t":.15},
}
MIXED = [.1,1.,-.31,.17,.09,-.06]


def equilibrium_error(k):
    v = k.mesh.vertices
    area = abs(np.linalg.det(v[:,1:]-v[:,:1]))/2
    result = np.zeros((3,3))
    for bary,w in QUAD:
        y,z = np.einsum("n,enc->ce",bary,v)
        tau = np.einsum("n,enfc->efc",weights(bary),k.units)
        moment = y[:,None]*tau[:,:,1]-z[:,None]*tau[:,:,0]
        result += w*np.einsum("e,efc->fc",area,np.dstack((tau,moment)))
    scale = max(k.mesh.profile.b,k.mesh.profile.h)
    return abs((result-np.eye(3)) * [[1,1,1/scale],[1,1,1/scale],[scale,scale,1]]).max()


def run(kind):
    section = dict(type=kind,dimensions=SECTIONS[kind])
    p = Profile.from_section(section)
    e = Engine()
    options = Options(area_fraction=1/64) if kind == "rectangle" else Options()
    rows = []
    for load in [MIXED,*np.eye(6).tolist()]:
        start = perf_counter()
        r = e.solve(p,.3,load,options)
        k = r.pop("kernel")
        row = dict(**r,load=load,seconds=perf_counter()-start,equilibrium_error=equilibrium_error(k),
                   peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        rows.append(row)
        print(kind,load,row["status"],row["fieldStatus"],row["maximumStatus"],
              row["history"][-1],round(row["seconds"],2),flush=True)
        if len(rows)==1:
            # Interior integration points with physical area weights and the three
            # independent unit shear/torsion fields for independent analytical comparison.
            v = k.mesh.vertices
            area = abs(np.linalg.det(v[:,1:]-v[:,:1]))/2
            points = []
            if not p.sharp:
                for bary,w in QUAD:
                    xy = np.einsum("n,enc->ec",bary,v)
                    units = np.einsum("n,enfc->efc",weights(bary),k.units)
                    points.append(np.column_stack((xy,w*area,units.reshape(-1,6))))
                np.savetxt(OUT/(kind+".csv"),np.concatenate(points),delimiter=",")
    (OUT/(kind+".json")).write_text(json.dumps(dict(section=section,nu=.3,options=asdict(options),rows=rows),
                                               allow_nan=False))


if __name__ == "__main__":
    OUT.mkdir(parents=True,exist_ok=True)
    if len(sys.argv)>1:
        run(sys.argv[1])
    else:
        for kind in SECTIONS:
            subprocess.run([sys.executable,"-m","tests.integration.section_solver_acceptance.matrix",kind],check=True)
        meta = dict(python=platform.python_version(),versions={k:version(k) for k in
                    ("sectionproperties","cytriangle","numpy","scipy","shapely")})
        (OUT/"environment.json").write_text(json.dumps(meta,indent=2))
