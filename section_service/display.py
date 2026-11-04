"""Version 1 binary field: header, physical nodes, Tri6 ids, then local samples.

Little endian. Header has four uint32s (version, nodes, elements, subdivisions).
Each straight triangle has 15 samples (i,j = 0..4, i+j<=4), with barycentric
coordinates (1-(i+j)/4,i/4,j/4). At each sample store four float64 components.
"""
import struct

import numpy as np

from section_solver.recovery import nodal, weights

DIVISIONS = 4
BARY = np.array([[(DIVISIONS-i-j)/DIVISIONS,i/DIVISIONS,j/DIVISIONS]
                 for i in range(DIVISIONS+1) for j in range(DIVISIONS-i+1)])
LIMIT = 24*1024**2
RETAINED_LIMIT = 48*1024**2


def export(kernel,load,limit=LIMIT):
    mesh = kernel.mesh
    nn,ne = len(mesh.nodes),len(mesh.triangles)
    size = 16+16*nn+24*ne+32*len(BARY)*ne
    if size>limit:
        return dict(status="resource_limit",bytes=size),None
    chunks = [struct.pack("<4I",1,nn,ne,DIVISIONS),mesh.nodes.astype("<f8",copy=False).tobytes(),
              mesh.triangles[:,[0,1,2,5,3,4]].astype("<u4").tobytes()]
    basis = weights(BARY)
    for start in range(0,ne,256):
        values = basis @ nodal(kernel,load,slice(start,start+256))
        sigma,ty,tz = values.transpose(2,0,1)
        vm = np.hypot(sigma,np.sqrt(3)*np.hypot(ty,tz))
        if not np.isfinite(vm).all():
            return dict(status="nonconverged",reason="Display stress overflow"),None
        chunks.append(np.dstack((values,vm)).astype("<f8",copy=False).tobytes())
    return dict(status="ok",bytes=size),b"".join(chunks)
