# MKE-F postprocessor JSON, version 8: recovery operators

Status: the only supported result format for planar and spatial models.
Exporter, HTTP runner, backend and viewer use v8. Recalculate old JSON v1–v7
from the original TXT model; there is no downgrade mode or converter.

```octave
exportPostprocessorData(model, result, 'result.json', struct('formatVersion', 8));
```

V8 is the default. Explicit `formatVersion=8` is accepted; all
other version values are errors. The viewer rejects old versions before scene
creation and the backend does not publish results from an old solver.

V8 `operatorStorage` defaults to `auto`. It also accepts `inline` (operators
on each element) and `table` (force the shared layout for comparisons). This option
works with the default or explicit v8. Static results omit operators in every mode.

## Compatibility and provenance

`version` is 8; `model.dimension` is explicitly 2 or 3. Both coordinate and
analysis contracts are described below. The exporter adds:

```json
{
  "metadata": {
    "exporter": {"name": "mkef-octave", "contractVersion": 1}
  },
  "model": {
    "recovery": {
      "version": 1,
      "coordinates": "local-physical-end",
      "stiffnessFormulation": "euler-bernoulli-saint-venant",
      "massFormulation": "consistent-v1"
    }
  }
}
```

This is an excerpt, not a complete result file. The stiffness formulation is
`linear-bar` for 112/1133, `euler-bernoulli` for 113 or
`euler-bernoulli-saint-venant` for 1436. The mass model
is `consistent-v1`. Recovery contract and exporter contract versions describe
their respective semantics, not an application release or Git revision.
Unknown contracts, coordinate conventions and formulations are rejected.
Recovery contract 1 uses inline operators; contract 2 uses the shared table
described below. Both recovery contracts are supported within v8
because they describe the two useful storage layouts, not old result formats.

## Operators

Every element in a modal or transient result must either contain
`recoveryOperators: {K: [...], M: [...]}` (recovery contract 1), or reference
one K/M pair in the model table (contract 2). Both matrices are dense arrays of rows,
4×4 for 112, 6×6 for 113/1133 and 12×12 for 1436. They are copied directly from the solver's
`localStiffness` and `localMass`, before transformation to system coordinates.
They are not condensed at released ends. Matrix ordering matches the rows of
`localFromSystem`: end 1 components followed by end 2 components.

The matrices are independent of the number of modes or time samples.
No global matrix is exported.

### Shared layout (recovery contract 2)

`model.recovery.version` is 2. `model.recoveryOperators` is a nonempty array of
`{id, K, M}` records, with contiguous one-based IDs in array order. Each element
has `recoveryOperatorId` instead of `recoveryOperators`. No mixture of inline
operators and references is allowed, and every table entry must be referenced.

The exporter groups complete K/M pairs only when every coefficient matches
exactly. It hashes the unrounded IEEE double representations, verifies hash
matches with `isequal`, and assigns IDs in first-element occurrence order.
A hash collision between unequal matrices is an explicit export error with
an inline-storage fallback suggestion, never a reason to merge them. No
geometrical or numerical tolerance is used for grouping. Adjacent representable
numbers and different mass matrices remain distinct. Elements of one exported
model share the same family, local component ordering and formulation.

Local bases, connectivity, `globalDOFs` and `localFromSystem` remain on each
element. Rotation and releases do not prevent sharing an identical uncondensed
local pair; their different kinematic mappings are applied after lookup.
The browser validates each pair once and uses the same array objects for all
references. It does not add inline copies to elements, including on JSON download.

`auto` compares the compact-JSON operator payloads, including table IDs and
element references, and uses the table only when it is strictly smaller.
One-element and all-unique inputs keep inline storage; the all-unique path
skips table construction. Sparse repeats which do not pay for reference overhead
also keep inline storage. The comparison targets ordinary compact output,
not gzip size or pretty-print whitespace. Other model/history fields are identical.

For a system displacement vector u and acceleration vector a:

```text
u_local = localFromSystem * u[globalDOFs]
a_local = localFromSystem * a[globalDOFs]
elastic = K * u_local
inertia = M * a_local
```

Modal acceleration uses `a = -eigenvalue * modeShape` at the chosen phase.
For the frame, the torsional mass term uses the polar area inertia Iy+Iz,
not the Saint-Venant torsion constant J. Matrix entries carry the units implied
by these mappings; mixed translation/rotation blocks have different dimensions.

`elastic + inertia` retains the existing viewer's meaning as that sum, not a
complete dynamic force including external loads or damping. Current section
forces/stresses use the elastic recovery. Distributed element loads remain
restricted to static analysis.

Static v8 results use recovery contract 1 and omit operators, references and
the shared table: the exported end forces already
include the equivalent local load subtraction, and section force polynomials
are authoritative. Both planar and spatial static normal-stress extrema come
from Octave. The browser evaluates the field at a selected point, without
repeating the search over the member length.

## Accepted section properties and static extrema

Every element with a nonzero `sectionId` carries `sectionProperties`, copied
from the solver's assembled element record. The separate `model.sections`
catalog remains the drawing geometry and the input to the sectional FEM service.
No result consumer rebuilds A, I or J from that geometry.

- Both dimensions export `area`, `centroidY`, `centroidZ`, principal `Iy`/`Iz`,
  `momentOfInertia`, `etaMin` and `etaMax`. In 2D the last three reflect the
  element's section angle; the centroid and principal inertias describe the
  unrotated profile. In 3D properties use principal local section axes: the
  section angle is already included in the element basis and is not applied twice.
- 3D additionally exports `torsionConstant`, `torsionMethod`, `approximate`
  and `open`. Manual J is preserved exactly. J is never the polar mass inertia.
- Geometry helpers construct contours using the supplied centroid, and check
  the exported fiber bounds against those contours. They do not integrate
  area/inertia or run the torsion series. Element material properties must
  agree with the supplied section properties.
- Static `analysis.elementResults[].sectionStress` contains `minimum`,
  `maximum` and an increasing `positions` array including 0 and 1. Critical
  points contain `value`, normalized `xi`, and `eta` (2D) or `y,z` (3D).
  Positions include stationary points and envelope switches; consumers use
  them for plot sampling, coloring and exports. An unassigned truss also
  exports a constant axial summary, with just `value,xi` and `[0,1]`.
- Static frames without section geometry have no normal-stress summary;
  bending stresses cannot be inferred from A/I/J alone. Mean `axialStress`
  remains distinct from extreme normal stress due to axial force and bending.
- Validation checks finite values, positive properties, references, bounds,
  position order and the stress value at the supplied critical material points.
  It does not independently prove a global maximum: independent formula and
  dense-sampling comparisons belong to Docker tests.
- Modal/transient stress fields use the same accepted properties, but still
  find extrema in the browser for the current mode/frame. No stress histories
  or dynamic extrema tables are added to the JSON.

Missing required snapshots or static summaries are errors. Development v8
files exported before this contract change must be regenerated from TXT;
the viewer never silently reconstructs missing properties or static extrema.

## Spatial data (dimension 3)

Element families are 1133 (truss, `dofPerNode=3`) and 1436 (frame,
`dofPerNode=6`). Physical nodal order is ux, uy, uz, thetaX, thetaY, thetaZ,
with translations only for a truss. Node and element IDs are contiguous and
one-based. `numberOfDOFs`, `dofRegistry`, `dofMap`, `nodes`, `elements`,
`supports`, `nodalLoads`, `elementLoads`, `mpcs`, `releases`, `warnings`,
`dependentDOFs` and `independentDOFs` describe the assembled model.

Frame rotational system coordinates may be a mixed active basis; released
element-end rotations have private coordinates. `dofMap` columns 4–6 therefore
are basis slots, not necessarily physical thetaX/Y/Z. `nodeKinematics` is COO
`{rows, columns, entries}` with one-based `[physicalRow, systemDOF, coefficient]`
entries. It maps signed system displacements to physical nodal components; it
is not a transformation for constrained reaction forces. Each element carries
`globalDOFs`, `localFromSystem`, `localBasis` (axes as rows), `referenceY` and
`sectionAngle`. Private rotations do not become shared nodal components.

For assigned frame sections, only `rectangle`, `circle`, `square_tube`,
`rectangular_tube` and `circular_tube` are accepted. The same policy applies at
TXT import, editing, export and solving. Unassigned catalog entries do not
restrict a model; all catalog profiles remain available to trusses and 2D.

| Analysis | Required data | Optional data |
| --- | --- | --- |
| static | `displacements`, `loadVector`, `reactions`, `elementResults`, `equilibriumResidual`, `supportReactions`, `mpcForces`, `mpcMultipliers` | Section extrema depend on section geometry as described above. |
| modal | `modeShapes`, `eigenvalues`, `frequenciesHz`, `angularFrequenciesRadPerSec`, `supportReactions`, `mpcForces`, `mpcMultipliers` | No transient histories. |
| transient | `time`, positive `timeStep`, `globalDOFIds`, `physicalDOFIds`, `sampling` | Selected histories and spectra. |

Static frame force coefficients are 6×4, rows N, Vy, Vz, Mx, My, Mz, ascending
powers of xi=x/L; positive N is tension. Static element results contain local
displacements, end forces, equivalent local load, and mean axial quantities.
Modal matrices use system DOF × mode; mode amplitude is arbitrary.

Transient matrices use selected system DOF × saved sample. `globalDOFIds`
gives their row order. `physicalDOFIds` is always present, including when empty;
it lists physical rows whose required system coordinates are selected.
`sampling` contains `originalSampleCount`, `exportedSampleCount`, `timeStride`.
The last original sample is retained even if it falls between stride positions.
`timeStep` is the solver step, not the spacing of the exported frames.

For frames, static/modal results also require `nodalSupportReactions` and
`nodalMPCForces`, physical nodal row × mode (vectors in statics). Truss
static/modal system forces already use physical coordinates. Transient export
of `supportReactions` or `mpcForces` also exports its corresponding physical
array, for both families, in `physicalDOFIds` order. Selecting `spectrum`
exports `spectrumFrequencyHz`, `displacementAmplitudeSpectrum` (system rows)
and `nodalDisplacementAmplitudeSpectrum` (physical rows). The latter comes
from FFT of transformed signed histories, never from transforming unsigned
amplitudes. Missing physical representations are invalid, not reconstructed.

Unselected histories are unavailable, not zero. This includes acceleration,
reactions and spectra; partial exports still support their available charts.
The two coordinate representations and distinct support/MPC forces must not
be conflated. JSON retains array dimensions for singleton and empty arrays.

## Planar data (dimension 2)

Planar v8 uses planar geometry, support, load, MPC and analysis representations,
with the same recovery metadata, operators and storage rules as 3D.
Dispatch by `model.dimension`, not by version.

- `dofPerNode` is 2 for truss 112 or 3 for frame 113. `dofLabels` is
  `["ux","uy"]` or `["ux","uy","rz"]`. `numberOfDOFs` is mandatory.
- Every element has positive `length`, `globalDOFs` in physical-end order and
  `localFromSystem` of size 4×4 or 6×6. It is the actual solver transformation,
  not an MPC condensation. For frames, local order is ux,uy,rz at end 1, then
  ux,uy,rz at end 2. Released rotations reference their private global DOFs.
- `dofMap` maps nodal coordinates; zero rz denotes absent shared rotation.
  A released model carries `releases`, `dofRegistry` and `warnings`. Without
  releases, element `globalDOFs` must agree with the two nodal rows. With
  releases, registry ownership validates every private end coordinate.
- MPC fields retain planar master records and dependent/independent DOF
  partitions. Histories and modal vectors include the solver's expanded system
  coordinates; recovery must not apply the MPC transformation a second time.
- Modal results add `eigenvalues`, one per mode, consistent with squared
  `angularFrequenciesRadPerSec`. They determine modal accelerations.
- Static element results add `localDisplacements`, `equivalentLocalLoadVector`
  and `sectionForceCoefficients`: three rows N,V,M, four ascending powers of
  xi=x/L. Trusses have zero V/M rows. End forces include load subtraction.
  Coefficients carry force units in N/V rows and moment units in M rows;
  xi is dimensionless. Signs follow the node-1-side cut convention:
  N(0)=-f1x, V(0)=-f1y, M(0)=-m1z and dM/dx=-V.
- Static diagrams, section fibers and the member-Mises request builder use
  these coefficients directly. Planar catalog `sectionStress` supplies extrema
  and critical positions. The viewer uses accepted `sectionProperties` and
  checks their consistency with exported A/I and the drawing fiber bounds.
- Partial transient fields, `globalDOFIds`, time stride and spectra retain their
  existing semantics. Missing displacement/acceleration data is unavailable,
  never implicitly zero. Partial histories remain usable for DOF charts.

Absent optional load/MPC/section/release groups remain absent;
presence controls these capabilities, without version-number feature checks.
The 2D interface continues to provide
static diagrams, modal deformations and transient DOF histories; it has no
dynamic member-force panel. The shared recovery evaluator
already applies planar K/M for future consumers and numerical checks.

Missing selected displacement DOFs yield unavailable elastic/section results.
Missing acceleration data yields unavailable inertia. Operators do not supply
the missing history values. They remain mandatory even in a partial dynamic
export, so file interpretation does not depend on which field a UI opens first.

## Validation

- Check exact dimensions and finite real values of K and M. The exporter also
  checks their JSON round trip: losing a nonzero entry to zero is an error,
  with a suggestion to rescale model units.
- For each row, require `sum(abs(A[i][j]-A[j][i])) <= 1e-12*max(abs(A))`.
- Do not require positive definiteness of K: rigid modes make it singular.
- Validate IDs, DOF mappings, releases, bases and history dimensions.
- Require operators for modal/transient; reject them for static results.
- For contract 2, validate the table length, all IDs/references, absence of
  inline operators and absence of unreferenced entries; validate each K/M once.
- Never reconstruct a missing or invalid operator from geometry.

Validation establishes structural consistency, not authenticity or proof that
an arbitrary imported matrix matches material properties. Numerical agreement
with the actual solver is checked by tests using independently exported solver
references. No new matrix formula is implemented in the v8 reader.

## Verification and limits

Run `./scripts/testing/test_recovery_pipeline.ps1`. It mounts current source read-only
in Docker and writes to a new `output/recovery-runs/<run-id>/` directory.
It generates model inputs, runs Octave, checks fresh numerical references
and exercises the browser offline. Numerical input
cases are retained and exported afresh as v8. The pipeline also runs browser
unit tests and the backend suite with the real Octave runner.

`planar-recovery/` contains 23 fresh 2D cases in inline/table/auto layouts and
independent solver references. It covers both families, all analyses, rotated
geometry, uniform/linear/global loads, releases, MPC, one mode and partial
histories. Tests include a closed-form pinned-beam moment, 22 invalid contracts,
operator provenance, static coefficient consumers, section extrema, file opening,
JSON downloads, SVG and the member-Mises request payload. The browser checks
operator provenance and uses an independent stress reference.

`recovery/browser-measurements.json` records JSON size, parse/validation time
and browser heap measurements for static/modal/transient files with 1 and
2,000 elements. Repeated elements are synthetic storage/reader benchmarks,
not independently solved large structures. Heap totals include browser/test
overhead; timing is diagnostic and hardware-dependent. The 30-second test
ceiling is a smoke-test timeout, not the final performance acceptance target.
Resource acceptance and existing product limits are documented in
[result architecture](../reference/result-recovery-architecture.md#принятые-пределы).

`recovery-pool/measurements.json` records exporter wall time (one observation
per layout) and exact file bytes for inline/table/auto. Its scenarios include
a physically solved 48-node, 80-member two-storey frame with a release and MPC,
and labelled synthetic storage cases with 2,000 repeated, mixed or unique pairs.
The synthetic histories come from a one-element solve and must not be treated
as a 2,000-member physical solution.

`recovery-pool/browser-measurements.json` adds gzip size (level 9), median
JSON.parse/validation time and retained heap. It measures each layout in a
separate Chromium page, warms up once, then records five repetitions. CDP
`HeapProfiler.collectGarbage` and `Runtime.getHeapUsage` measure retained parsed
data with the input string already present. File I/O, text transfer and GC are
outside the parse/validation timer. This is not total renderer/GPU memory.

Octave 11.3 native `jsonencode` can replace small nonzero components with zero.
Result exports therefore use `encodeResultJson`: finite numeric doubles are
formatted with `%.17g`, while native encoding handles strings and booleans.
The prepared schema arrays retain their dimensions, including singleton rows.
Signed zero is normalized to positive zero, matching browser JSON downloads;
nonzero values, including subnormals, round-trip through JavaScript JSON.parse
without changing their IEEE double representation. No physical tolerance or
automatic removal of numerical noise is applied by the serializer.

`test_result_numbers.m` / `.cjs` check 887 bit representations, including tiny
operators, mode shapes and transient histories. Extreme magnitudes are synthetic
transport tests, not evidence of required accuracy in ordinary structural models.
Recovery references now use the original solver values: native Octave jsondecode
is not used as a numerical oracle. This is a transport guarantee, not a guarantee
against overflow, ill-conditioning or rounding inside the FE calculation.

`result-resources/browser-measurements.json` records static/modal/transient
results for a solved 200-member planar beam and an 80-member spatial frame,
plus a labelled synthetic 2,000-member planar display case. Parse/validation
and retained heap use three measured repetitions after warm-up; the input
string, file I/O, GC and GPU memory are outside these measurements. Animation
counts actual scene updates over one second, not requestAnimationFrame calls.
