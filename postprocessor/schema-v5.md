# MKE-F postprocessor JSON format, version 5

Version 5 is emitted when at least one full-span linear distributed load is
present. Models without type 21 retain their v1–v4 version selection. A v5
dataset is static and contains the v2 `nodalLoads`, `elementLoads` and
`equivalentLocalLoadVector` fields. Older readers must reject version 5.

## Distributed loads

`model.elementLoads` may mix these records in any order:

- Uniform: `{type:20, elementId, coordinateSystem, qx, qy}` (unchanged).
- Linear: `{type:21, elementId, coordinateSystem, qx1, qy1, qx2, qy2}`.

Only the fields belonging to that load type are serialized. Elements must be
frames (113). `coordinateSystem` is 1 for local axes or 2 for global XY.
Components are finite forces per unit actual member length; all components
of a single record cannot be zero. Values are attached to element ends 1 and
2, in `nodeIds` order, and interpolate as `q(xi)=(1-xi)*q1+xi*q2`.
Zero at either end is valid, as are opposite signs, equal endpoints and
multiple records on one element. All records add together.

`equivalentLocalLoadVector` is the exact consistent vector in order
`[Fx1,Fy1,Mz1,Fx2,Fy2,Mz2]`. With local endpoint components a and b:

```text
[L(2ax+bx)/6, L(7ay+3by)/20, L²(3ay+2by)/60,
 L(ax+2bx)/6, L(3ay+7by)/20, -L²(2ay+3by)/60]
```

The v2 meaning of `localEndForces = Klocal*ulocal - flocal` is unchanged.
The mean axial force/stress fields retain their previous meanings. For
summed local q1 and gradient k=(q2-q1)/L, diagrams use:

```text
N(x) = -Fx1 - qx1*x - kx*x²/2
V(x) = -Fy1 - qy1*x - ky*x²/2
M(x) = -Mz1 + Fy1*x + qy1*x²/2 + ky*x³/6
```

Plot samples include every interior stationary point, including both roots
of V for M when present, and zeros of qx/qy for N/V.

## Optional constraints

MPC data follows v3 and is included only when MPCs are present. With end
releases, all v4 DOF ownership fields are required: `numberOfDOFs`,
`dofRegistry`, nonempty `releases`, `warnings`, and each element's
`globalDOFs`. Without releases, v5 retains the ordinary nodal DOF map and
does not require these fields. Releases and MPCs may coexist.

## Case input and editor

After the elements, add a repeatable section:

```text
eload_linear
1
21,1,2,0,-300,0,-900
```

This applies global Y intensity -300 at end 1 and -900 at end 2. Sections
may coexist with `eload_uniform`; serialization preserves record order.
The editor offers Uniform/Linear, copies uniform values to both ends, and
uses the endpoint mean when converting back. This preserves force but does
not in general preserve the moment of the original linear distribution.
Coordinate-system changes rotate both vectors. Geometry/property editing
and deletion follow the existing load invalidation and undo rules.

Examples: `CaseLinearFrame`, `CaseTriangleFrame`, `CaseLinearRelease`,
`CaseLinearMPC`. Verification: `test_linear_loads`, API integration tests,
and `scripts/test_linear_loads_browser.cjs` after exporting those fixtures.
