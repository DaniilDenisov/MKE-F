# MKE-F postprocessor JSON format, version 4

Version 4 is emitted for a model with at least one frame-end Mz release.
Models without releases retain the v1/v2/v3 selection and DOF ordering.
All existing static, modal and transient fields retain their meanings.
Element-load and MPC fields are included when those features are present;
a v4 dataset need not contain MPCs. Readers must reject unsupported versions.

## Model and DOF ownership

- `numberOfDOFs`: actual number of assembled coordinates, including private rotations.
- `dofPerNode`: possible nodal components (3 for frames), not an assembly size multiplier.
- `dofMap`: nodes × 3; positive entries identify shared nodal coordinates. A zero
  thetaZ entry means that all attached element ends release Mz. It is not a zero displacement.
- `releases`: nonempty array of `{elementId, end, component}`; end is 1 or 2,
  component is `"Mz"`. References must exist and each element/end occurs once.
- `elements[].globalDOFs`: six positive identifiers in global-component order
  `[ux1, uy1, thetaZ1, ux2, uy2, thetaZ2]`.
- `dofRegistry`: one record per global coordinate, ordered by `id` (1-based):
  `{id, kind, nodeId, component, elementId, end}`. Components are `ux`, `uy`, `thetaZ`.
  For `kind: "node"`, elementId and end are 0. For `kind: "elementEnd"`, component
  is thetaZ, elementId/end identify its sole owner and nodeId identifies its geometric node.
- `warnings`: array of diagnostic strings, including ignored thetaZ restraints.

Shared coordinates are numbered in node/component order, omitting rotations
disconnected by releases; private rotations follow in element/end order.
The registry, element maps, node maps and releases must agree. A private rotation
cannot be shared by another element. Every coordinate has exactly one owner.
An otherwise isolated node is not silently repaired by deleting its DOFs.

`supports` retains the original input masks. A mask on an absent thetaZ has no
effect; translations still apply. Viewers render the effective restraints and
show the warning. Nonzero nodal Mz, MPC terms and transient nodal monitors cannot
reference an absent thetaZ. Internal rotations cannot be targeted by case-file SPC/MPC.

## Results and forces

Static vectors and modal matrix rows use `numberOfDOFs`, not `nodes.length * 3`.
Transient `globalDOFIds` may include any registry coordinate; default export
includes all coordinates. `selectedGlobalDOFs`, `transientFields`, `timeStride`
and final-sample retention work as before. Single-item arrays and matrix nesting
are preserved. MPC multipliers retain equation ordering and are not DOF-filtered.

Static `localEndForces = Klocal * ulocal - flocal` uses actual end rotations and
includes the consistent element-load vector. Released components vanish to
numerical tolerance without being overwritten. Under dynamics the release
condition is the corresponding residual component `M*a + K*u - F = 0`, or
`K*phi - omega^2*M*phi = 0` for each mode. The elastic term alone need not vanish
with consistent mass. No artificial zeroing or static condensation is used.

The v3 definitions of `reactions`, `supportReactions`, `mpcForces` and
`mpcMultipliers` continue to apply when MPCs are present.
The global static moment balance counts each rotational coordinate once.

## Viewing

Frame interpolation uses `elements[].globalDOFs`, including private rotations.
Absent nodal thetaZ is displayed as “not defined”. Element cards show both end
rotations; modal card values use the stored eigenvector amplitude. Transient
charts offer element/end/thetaZ alongside nodal quantities. Geometry playback
requires displacement histories for every global coordinate; partial exports
remain available for charts and available element-end values.

## Case input

The repeatable section follows an element section and does not alter `elems_113`:

```text
releases
2
1,1,Mz
1,2,Mz
```

The editor writes one section sorted by element and end. Releases survive
geometry/property edits; deletion removes and renumbers associated records.
Changing a release preserves existing element loads. N/V releases are unsupported.

Examples: `CaseReleaseStatic`, `CaseReleaseModal`, `CaseReleaseTransient`.
Verification: `test_releases`, backend API tests, `scripts/test_releases_browser.cjs`.
