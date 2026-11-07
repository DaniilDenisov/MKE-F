"""Validate matrix output against analytical fields and an independent grid search."""
import json
import os
from pathlib import Path

import numpy as np

from .reference import Reference
from .matrix import SECTIONS


def verify(directory,kind):
    case = json.loads((directory/(kind+".json")).read_text())
    ref = Reference(case["section"])
    l2 = None
    if not ref.sharp:
        data = np.loadtxt(directory/(kind+".csv"),delimiter=",")
        if ref.round and ref.t:
            data = data[np.hypot(*data[:,:2].T)>=ref.h/2-ref.t]
        y,z,w = data[:,:3].T
        expected = ref.fields(y,z,case["nu"])
        error = data[:,3:].reshape(-1,3,2)-expected
        l2 = np.sqrt(np.einsum("e,efc->f",w,error**2)/np.einsum("e,efc->f",w,expected**2))
        assert max(l2)<.005,(kind,l2)
        y,z = ref.grid()
        fields = ref.fields(y,z,case["nu"])
    for row in case["rows"]:
        n,vy,vz,t,my,mz = row["load"]
        h = row["history"][-1]
        assert row["fieldStatus"]==row["JStatus"]=="ok"
        assert h["geometry_error"]<=1e-4 and max(h["field_changes"])<=.005
        assert abs(row["J_FE"]/ref.J-1)<=.001 and row["equilibrium_error"]<=.001
        if ref.sharp and any((vy,vz,t)):
            assert row["maximumStatus"]=="unavailable" and row["maximum"]["singularPeak"]
            assert "value" not in row["maximum"]
            continue
        if ref.sharp:
            peak = abs(n)/ref.A+abs(my)*ref.b/(2*ref.Iy)+abs(mz)*ref.h/(2*ref.Iz)
        else:
            sigma = n/ref.A+my*z/ref.Iy-mz*y/ref.Iz
            ty,tz = np.einsum("efc,f->ce",fields,[vy,vz,t])
            peak = np.hypot(sigma,np.sqrt(3)*np.hypot(ty,tz)).max()
        row["peak_error"] = abs(row["maximum"]["value"]/peak-1)
        assert row["status"]=="ok" and row["peak_error"]<.01,(kind,row["load"],row["peak_error"])
        assert all(item["peak_change"]<=.005 for item in row["history"][-2:])
    print(kind,"passed",l2,flush=True)
    return dict(**case,refJ=ref.J,l2=None if l2 is None else l2.tolist())


if __name__=="__main__":
    directory = Path(os.environ.get("SECTION_REPORT_DIR","output/section-solver"))
    report = dict(environment=json.loads((directory/"environment.json").read_text()),
                  profiles=[verify(directory,kind) for kind in SECTIONS])
    (directory/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False))
