> Historical schema, no longer supported. Recalculate the original TXT model and use [v8](schema-v8.md).

# MKE-F spatial postprocessor JSON, version 7

This documents the existing spatial format; v1–v6 describe planar results.
The legacy exporter still selects v7 for spatial models by default.

## Model and coordinates

`format` is `mkef-postprocessor`, `version` is 7, and `model.dimension` is 3.
The supported homogeneous element families are 1133 (truss, three translations
per node) and 1436 (frame, three translations and three rotation components).
Physical order is ux, uy, uz, thetaX, thetaY, thetaZ, truncated to translations
for a truss. Node and element IDs are contiguous and one-based.

`model.dofRegistry` describes system coordinates. With releases, a frame's
nodal rotational coordinates can be an active basis rather than the three
physical rotations. Private element-end rotations have their own coordinates.
`dofMap` therefore does not always map directly to physical thetaX/Y/Z.

`nodeKinematics` is sparse COO `{rows, columns, entries}`. Each entry is
`[physicalRow, systemDOF, coefficient]`, with one-based indices. It maps system
displacements to physical nodal displacements. It is not a force transformation
for constrained nodal reactions.

Each element has `globalDOFs` and `localFromSystem`. Multiplying the latter by
the system vector restricted to `globalDOFs` gives its physical local end DOFs.
Order is all components of end 1 followed by all components of end 2.
`localBasis` contains the local axes as rows, mapping global vectors to local.
`referenceY` and `sectionAngle` document the frame orientation.

Geometry, material/section properties, the section catalog, supports, loads,
MPCs, releases, dependent/independent DOFs and warnings are included. V7 does
not export stiffness or mass matrices. The legacy browser reconstructs local
operators when recovering modal/transient element forces.

## Analysis

- Static results contain system displacements, loads and reactions, and
  `elementResults` with local displacements, end forces, equivalent local load
  vectors, mean axial strain/stress/force and optional section stress extrema.
- Frame `sectionForceCoefficients` is 6×4: rows N, Vy, Vz, Mx, My, Mz;
  columns are ascending powers of xi=x/L. Positive N is tension.
- Modal results contain `modeShapes` (system DOF × mode), `eigenvalues`,
  `frequenciesHz` and `angularFrequenciesRadPerSec`. Modal amplitude is arbitrary.
- Transient arrays are exported system DOF × exported time sample. The row order
  is `globalDOFIds`, not necessarily every system DOF. Optional arrays include
  displacements, velocities, accelerations, loads, reactions and spectra.
  `time`, `timeStep` and `sampling` describe the selected samples.
- `physicalDOFIds` identifies exported physical nodal rows for transient
  reactions and spectra. A mixed-coordinate amplitude spectrum cannot be
  reconstructed by transforming already unsigned system amplitudes; the solver
  exports the spectrum of transformed signed histories.
- Support reactions, MPC forces and MPC multipliers are distinct quantities.
  Private end rotations do not become physical nodal rotation coordinates.

Missing fields in a partial export mean unavailable, not zero. JSON arrays keep
their array shape even with one element, one DOF or one sample. Empty arrays are
arrays, and all exported numeric values must be finite.

The authoritative implementations are `createSpatialPostprocessorData.m`,
`prepareSpatialJsonData.m` and `postprocessor/js/spatial-data.js`. Fixed examples
are checked into `tests/fixtures/postprocessor-legacy/` and must not be regenerated
as part of a compatibility test.
