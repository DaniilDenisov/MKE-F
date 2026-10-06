# MKE-F postprocessor JSON format, version 2

Version 2 extends [version 1](schema-v1.md) for static uniform frame loads.
The exporter uses version 2 only when the model has element-load records;
other results retain version 1. The viewer accepts both. A version-1-only
viewer must reject version 2 rather than plot incorrect member diagrams.

All version 1 static fields remain, including the full assembled
`analysis.loadVector`. Version 2 requires these additional fields:

- `model.nodalLoads`: array of `{type: 10, nodeId, fx, fy, mz}` records.
  These are the original applied nodal actions, not equivalent element loads.
- `model.elementLoads`: array of `{type: 20, elementId, coordinateSystem, qx, qy}`.
  `coordinateSystem` is 1 (local) or 2 (global XY); intensities are force per
  unit actual member length. IDs are one-based and must identify frame 113.
  Components must be finite and not both zero. Multiple records add together.
- `analysis.elementResults[].equivalentLocalLoadVector`: six components in
  `[ux1, uy1, rz1, ux2, uy2, rz2]` order. This is the sum of consistent
  load vectors for the element; it is zero for unloaded elements.

Collections and vectors are arrays even when empty or containing one record.
Directions and coordinate conversion always use the original geometry.
The viewer draws original nodal actions and element distributions separately;
it does not additionally draw the equivalent actions from `loadVector`.

`localEndForces = Klocal * ulocal - equivalentLocalLoadVector` represents
forces exerted on the member at its ends. With the sum of local intensities
`qx, qy`, length `L`, and end-force vector `f`, the section diagrams are:

```text
N(x) = -f[0] - qx*x
V(x) = -f[1] - qy*x
M(x) = -f[2] + f[1]*x + qy*x*x/2,  0 <= x <= L
```

The viewer includes an interior stationary point of M in its samples and
maximum. `axialStrain`, `axialStress`, and `axialForce` are member means;
use N(x) for the spatially varying axial force. Bending stresses are not
exported. Displacement graphics retain finite-element interpolation, not
an analytical reconstruction of the internal loaded-beam deflection.

For both original and deformed geometry, distributed-load arrows keep their
original physical direction. SVG/PNG snapshots include the visible arrows
and sampled diagrams.

Support catalog extension: frame support types 5, 6, 7 restrain only ux, uy, thetaZ respectively. Trusses accept types 1–6, export canonical types 1, 3, 2, and reject 7. Existing types retain their original meaning. Older readers reject unknown support types.
