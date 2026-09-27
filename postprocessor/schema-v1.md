# MKE-F postprocessor JSON format, version 1

The root object contains `format: "mkef-postprocessor"`, `version: 1`,
`metadata`, `model`, and `analysis`. Identifiers and global degrees of freedom
are one-based.

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
- Transient data contains the exported `time` samples, `globalDOFIds`, selected
  histories in `[selected dof][time index]` orientation, optional spectrum data,
  and explicit sampling metadata. `globalDOFIds[row]` identifies the global DOF
  represented by a history or spectrum row.

All numbers are finite JSON numbers. Complex, sparse, `NaN`, and infinite
values are invalid. Missing optional transient fields are omitted. Readers must
reject unknown major versions and may ignore additional fields in version 1.
