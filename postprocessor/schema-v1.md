# MKE-F postprocessor JSON format, version 1

The root object contains `format: "mkef-postprocessor"`, `version: 1`,
`metadata`, `model`, and `analysis`. Identifiers and global degrees of freedom
are one-based.

`metadata.units` provides display labels for length, force, moment, stress, and
time. The postprocessor does not convert numerical values. For an SI model
exported with `length: "m"`, `force: "N"`, `moment: "N*m"`, and
`stress: "Pa"`, coordinates and displacements are displayed in metres, forces
(including loads, reactions, and axial forces) in newtons, moments in
newton-metres, and stresses in pascals. The labels must match the consistent
unit system used by the source model.

`model` contains `dimension`, `dofPerNode`, `dofLabels`, `nodes`, `elements`,
`dofMap`, and `supports`. Element properties are named values. Solver matrices,
element matrices, and transformations are never exported.

Fields described as vectors, matrices, or collections are always JSON arrays,
including when they contain exactly one item. Matrix fields are arrays of row
arrays, including one-column and one-row matrices.

`analysis.type` is `static`, `modal`, or `transient`:

- Static data contains full one-dimensional displacement, load, and reaction
  vectors, the three-component equilibrium residual, and one result record per
  element.
- Modal data contains frequency vectors and `modeShapes[dof][mode]`.
  Mode-shape amplitudes and signs are arbitrary eigenvector conventions, not
  physical displacements in the declared length unit. The viewer copies each
  selected column, fixes its sign deterministically, and normalizes it for
  display without changing the JSON data. Its animation cycle is illustrative
  and does not represent the physical modal frequency.
- Transient data contains the exported `time` samples, `globalDOFIds`, selected
  histories in `[selected dof][time index]` orientation, optional spectrum data,
  and explicit sampling metadata. `globalDOFIds[row]` identifies the global DOF
  represented by a history or spectrum row. The optional history fields are
  `displacements`, `velocities`, `accelerations`, `loadHistory`, and `reactions`;
  `spectrumFrequencyHz` and `displacementAmplitudeSpectrum` are an optional pair.
  The viewer uses the exported `time` values directly, keeps one deformation
  scale for the full exported range, and never interpolates samples removed by
  export `timeStride`. Its frame-stride control only changes which exported
  samples are shown during stepping and playback.

`analysis.sampling` records `originalSampleCount`, `exportedSampleCount`, and
`timeStride`. A viewer must disclose decimation and omitted optional fields.
Transient `reactions` is the full dynamic residual `M*a + K*u - F`; only its
restrained rows are support reactions. SVG/PNG snapshots contain the currently
visible frame and chart, not a reconstructed or subsequently animated frame.

All numbers are finite JSON numbers. Complex, sparse, `NaN`, and infinite
values are invalid. Missing optional transient fields are omitted. Readers must
reject unknown major versions and may ignore additional fields in version 1.
