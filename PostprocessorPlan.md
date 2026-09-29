# Unified HTML/SVG FEM Postprocessor Plan

## Summary

Create a single offline postprocessor for every analysis currently implemented
in MKE-F. GNU Octave remains the numerical backend and exports a small,
versioned result dataset. A standalone HTML application written in plain
JavaScript renders the model and results as SVG.

The implementation is organized as three completed layers: a retained static
scene, static diagrams plus modal shapes, and transient playback plus a
selected-DOF chart. Geometry is mounted once per dataset; layer switches,
scales, modes, and time steps update existing SVG attributes. Numerical values
are shown in HTML details and legend panels rather than as unbounded labels on
the structural canvas.

The supported analyses are:

- static (`analysisType = 'static'`);
- modal (`analysisType = 'modal'`);
- Newmark transient (`analysisType = 'transient'`).

The application must work by opening `postprocessor/index.html` directly from
disk. It must not require a web server, npm, a build step, a CDN, or an internet
connection. User data is loaded with a file picker or drag-and-drop through the
browser `FileReader` API, avoiding `fetch()` restrictions on `file://` pages.

The intended workflow is:

```octave
problem = StructFEProblem('ANSYSBeamStatic01.txt');
model = problem.GetAnalysisModel();
result = problem.RunStatic();

exportPostprocessorData(model, result, 'static-frame.json');
```

Then the user opens `postprocessor/index.html`, drops `static-frame.json` onto
the page, inspects the result, and exports an SVG or PNG for the documentation.

## Why HTML/SVG Is the Primary Renderer

HTML/SVG is preferred over an Octave plotting GUI because it provides:

- the same interface on Windows, Linux, and macOS;
- crisp vector geometry, labels, arrows, and diagrams;
- convenient zoom, pan, tooltips, layer switches, time sliders, and animation;
- direct SVG and PNG export;
- reuse inside the offline HTML reference;
- no dependence on Octave graphics toolkits, font caches, or a display server;
- clear separation between numerical analysis and presentation.

Octave plotting remains available only for backward compatibility. New
postprocessor features must not be implemented twice in Octave and JavaScript.

## Current Repository Contracts

The implementation must use the current production contracts rather than
introducing a second solver:

- every solver returns `analysisType`;
- the model contains `nodeCoordinates`, `dofMap`, `elementData`,
  `fixedBoundaryConditions`, and `forceBoundaryConditions`;
- each element contains its type, node numbers, properties, length, global DOF
  list, transformation, and element matrices;
- static results contain full displacements, load vector, reactions,
  equilibrium residual, and recovered element results;
- modal results contain frequencies and full mode-shape columns;
- transient results contain a common time grid and histories of displacement,
  velocity, acceleration, load, reaction, and displacement spectrum.

The exporter must copy only information required for postprocessing. It must
not serialize the global stiffness/mass matrices or the element matrices.

## Scope

The first complete version includes:

- original and deformed geometry;
- support symbols, node/element labels, nodal loads, and reactions;
- axial force and stress for element 112;
- `N`, `V`, and `M` diagrams for element 113 under the currently supported
  nodal loading model;
- modal shapes with frequency labels;
- transient histories, displacement spectrum, a time-step snapshot, and
  interactive animation;
- SVG and PNG export;
- reusable exported illustrations for the offline reference.

The first version does not include:

- new finite-element or analysis types;
- nonlinear, plastic, contact, or 3-D visualization;
- distributed loads or diagrams caused by them;
- result editing inside the browser;
- GIF/video export;
- a desktop wrapper such as Electron;
- a local HTTP service connecting a live Octave process to the browser;
- binary/compressed result files.

These exclusions keep the application offline and simple while leaving a
versioned extension path.

## Architecture

```text
Input file
   |
   v
StructFEProblem / pure solvers
   |
   +-- model struct
   +-- result struct
          |
          v
createPostprocessorData(model, result, options)
          |
          +-- validated Octave struct (testable without graphics)
          v
exportPostprocessorData(..., "result.json")
          |
          v
postprocessor/index.html
          |
          +-- FileReader JSON loading
          +-- schema validation
          +-- SVG structural renderer
          +-- SVG charts and diagrams
          +-- SVG/PNG download
```

Numerical solvers must not know that the web postprocessor exists. The browser
must not reproduce the FEM solve. It receives solved vectors and performs only
display interpolation, normalization, selection, and diagram construction.

## Proposed Repository Layout

```text
createPostprocessorData.m
exportPostprocessorData.m

postprocessor/
  index.html
  assets/
    postprocessor.css
  js/
    namespace.js
    data-model.js
    geometry.js
    svg-renderer.js
    charts.js
    static-view.js
    modal-view.js
    transient-view.js
    export.js
    app.js
  tests/
    index.html
    test-runner.js
    geometry-tests.js
    data-model-tests.js
  examples/
    README.md

tests/
  test_postprocessor_export.m

reference/examples/
  generate_postprocessor_datasets.m
```

The browser code must use ordinary deferred scripts in a fixed order, not ES
module imports. Chromium browsers commonly restrict module loading from
`file://`, while classic local scripts work without a server. Each file should
extend one namespace such as `window.MKEFPost` from an IIFE and avoid unrelated
global variables.

No minified or generated JavaScript is committed. Source files are loaded
directly by `index.html`.

## Public Octave API

### `createPostprocessorData`

```octave
data = createPostprocessorData(model, result)
data = createPostprocessorData(model, result, options)
```

This is the pure, testable conversion layer. It validates consistency between
the model and result and returns a JSON-ready struct without opening files or
figures.

### `exportPostprocessorData`

```octave
data = exportPostprocessorData(model, result, filename)
data = exportPostprocessorData(model, result, filename, options)
```

It delegates to `createPostprocessorData`, encodes the result with Octave's
`jsonencode`, writes UTF-8 JSON, and returns the same struct for inspection.
The writer must fail clearly on a bad path and must not leave a partially
written target file. Write to a sibling temporary file and replace the target
only after a successful close.

Supported exporter options:

- `title`: dataset title shown in the GUI;
- unit labels: `lengthUnit`, `forceUnit`, `momentUnit`, `stressUnit`,
  `timeUnit`; the solver remains unit-agnostic and the exporter does not convert
  values;
- `transientFields`: explicit subset of `displacements`, `velocities`,
  `accelerations`, `loadHistory`, `reactions`, and spectrum;
- `timeStride`: positive integer, default 1; decimation is never silent and is
  recorded in metadata;
- optional `selectedGlobalDOFs` for intentionally reduced time-history files;
- `prettyPrint` when supported without changing numerical content.

Static and modal exports reject transient-only options rather than silently
accepting a typo.

## Versioned JSON Contract

Top-level shape:

```json
{
  "format": "mkef-postprocessor",
  "version": 1,
  "metadata": {
    "title": "Static frame",
    "createdUtc": "2026-09-26T19:00:00Z",
    "generator": "MKE-F",
    "octaveVersion": "11.3.0",
    "units": {
      "length": "m",
      "force": "N",
      "moment": "N*m",
      "stress": "Pa",
      "time": "s"
    }
  },
  "model": {},
  "analysis": {}
}
```

Required model data:

```json
{
  "dimension": 2,
  "dofPerNode": 3,
  "dofLabels": ["ux", "uy", "rz"],
  "nodes": [
    {"id": 1, "x": 0.0, "y": 0.0}
  ],
  "elements": [
    {
      "id": 1,
      "type": 113,
      "nodeIds": [1, 2],
      "properties": {"area": 0.0001, "youngsModulus": 2e11,
                     "density": 7800, "momentOfInertia": 3.33e-9}
    }
  ],
  "dofMap": [[1, 2, 3]],
  "supports": [
    {"type": 1, "nodeId": 1}
  ]
}
```

Node and element identifiers remain one-based and explicit. JavaScript must not
silently treat them as zero-based array offsets.

The analysis object is a discriminated union selected by `type`.

Static analysis:

```json
{
  "type": "static",
  "displacements": [],
  "loadVector": [],
  "reactions": [],
  "equilibriumResidual": [0, 0, 0],
  "elementResults": []
}
```

Modal analysis:

```json
{
  "type": "modal",
  "frequenciesHz": [],
  "angularFrequenciesRadPerSec": [],
  "modeShapes": []
}
```

Transient analysis:

```json
{
  "type": "transient",
  "time": [],
  "displacements": [],
  "velocities": [],
  "accelerations": [],
  "loadHistory": [],
  "reactions": [],
  "spectrumFrequencyHz": [],
  "displacementAmplitudeSpectrum": [],
  "sampling": {
    "originalStepCount": 1000,
    "exportedStepCount": 501,
    "timeStride": 2
  }
}
```

Matrix orientation is part of the schema and must be documented and tested:

- `modeShapes[dof][mode]`;
- transient fields `[dof][timeIndex]`;
- `dofMap[nodeIndex][localDofIndex]`;
- `elementResults[elementIndex]`.

Every numeric value must be finite. `NaN`, `Inf`, sparse-matrix encoding, and
complex values are rejected before JSON generation. Missing optional fields are
omitted, not encoded with ambiguous empty placeholders.

The web app rejects unknown major versions with an actionable message. Additive
fields within version 1 may be ignored.

## Browser Loading and Validation

The opening screen contains:

- a file picker;
- a large drag-and-drop target;
- format/version help;
- an error panel that reports the JSON path of invalid data.

Loading uses `FileReader.readAsText`. Do not use `fetch()` for user files. JSON
is parsed as data only. Dataset strings are inserted through `textContent` and
never through `innerHTML`, preventing a result file from injecting markup or
script.

Browser validation checks:

- format and supported version;
- unique positive node and element IDs;
- valid element node references;
- supported element types 112 and 113;
- DOF-map dimensions and range;
- result-vector and matrix dimensions;
- matching time-history column counts;
- finite numeric values;
- valid support and load node references;
- availability of fields required by the selected view.

Validation errors must not leave a partially rendered old/new mixture. Parse
and validate into a new in-memory dataset, then replace the current application
state only on success.

## GUI Layout

Desktop layout:

```text
+--------------------------------------------------------------------+
| Open JSON | Analysis | View | Export SVG | Export PNG | Reset view |
+----------------------+---------------------------------------------+
| Layers and settings  |                                             |
|                      |              SVG viewport                   |
| [x] Original         |                                             |
| [x] Deformed         |                                             |
| [x] Supports         |                                             |
| [x] Loads            |                                             |
| [x] Reactions        |                                             |
| [ ] Node numbers     |                                             |
| [ ] Element numbers  |                                             |
| Scale: Auto / value  |                                             |
+----------------------+---------------------------------------------+
| Context panel: mode selector, time slider, chart, legend, values  |
+--------------------------------------------------------------------+
```

On desktop screens the toolbar and application shell fit within one browser
viewport. The structural canvas consumes the remaining flexible height so the
active legend, selection panel, and status stay visible without document-level
vertical scrolling; the settings column scrolls independently when necessary.
On narrow screens the settings panel moves above the viewport and normal page
scrolling is restored. Controls must remain usable without horizontal page
scrolling.

The viewport uses one responsive `<svg>` with a `viewBox`. It supports:

- wheel zoom centered at the pointer;
- pointer-drag pan;
- reset-to-fit;
- selection/highlight of a node or element;
- tooltips with IDs, coordinates, displacements, forces, and stress where
  available;
- keyboard-accessible controls and visible focus styles.

Model coordinates use mathematical positive Y upward. SVG positive Y is
downward, so the transformation must be isolated and tested. Text labels and
symbols must not be mirrored by the geometry transform.

## SVG Layer Model

Create stable groups in this order:

```text
viewport
  original-geometry
  diagrams
  deformed-geometry
  supports
  loads
  reactions
  nodes
  labels
  selection-overlay
```

Each element curve is a `<path>` with `data-element-id`. Nodes and result
symbols carry corresponding `data-*` identifiers. Updating a transient frame
changes existing path attributes instead of destroying and rebuilding the
entire SVG tree.

Colors are not the sole carrier of meaning: original/deformed geometry and
different result types also use line styles, symbols, and labels. The default
palette must remain readable on a white background and when printed.

## Deformed Geometry Algorithms

### Element 112

Interpolate the two global translational displacement components linearly
between the end nodes. This represents axial deformation and rigid rotation of
the straight truss member.

### Element 113

For a frame element:

1. Compute its direction cosines from the two exported node coordinates.
2. Transform global nodal DOFs to the local ordering
   `[u1, v1, theta1, u2, v2, theta2]`.
3. Interpolate axial displacement linearly.
4. Interpolate transverse displacement with cubic Hermite functions:

   ```text
   v(xi) = N1*v1 + N2*theta1 + N3*v2 + N4*theta2,
   0 <= xi <= 1
   ```

   `N2` and `N4` include the element length.
5. Add scaled displacement to the local undeformed axis.
6. Rotate and translate the sampled points back to global coordinates.

The default is 31 samples per frame element. The user may change this between
5 and 201. Geometry tests must check both endpoints, end slopes, horizontal and
inclined members, and sign conventions. Visual inspection alone is not enough.

### Automatic Scale

Automatic displacement scale targets a visible maximum deformation, for
example 10% of the model bounding-box diagonal. Degenerate one-dimensional
bounding boxes fall back to the longest element length.

- A zero field uses scale 1 and displays “zero deformation”.
- Manual scale accepts only finite positive values.
- Modal shapes are normalized before display scaling.
- A transient animation uses one scale computed from the selected range so the
  structure does not visually breathe between frames.
- Diagram, force-arrow, and reaction-arrow scales are independent and shown in
  the legend.

## Static View

The default static view shows:

- undeformed geometry as a muted dashed line;
- deformed geometry as a solid line, coloured by actual resultant displacement
  magnitude `|u| = sqrt(ux^2 + uy^2)` by default;
- a compact Viridis colour legend with zero, half-maximum, and maximum values
  in the declared length unit; a zero field uses one neutral value;
- the displacement magnification factor;
- support symbols derived from support types 1--4;
- nodal force and moment symbols from `loadVector`;
- reaction symbols only at restrained DOFs, with thicker purple dashed force
  arrows and thicker amber solid moment arcs;
- optional node and element labels.

Displacement colours always represent the unscaled physical result and do not
change with automatic or manual display magnification. Element 112 uses linear
displacement interpolation and element 113 uses the same cubic Hermite field as
the deformed curve. SVG gradients use at most nine stops per element so the
renderer retains one path per element and preserves the large-model target.

`Deflection |u|`, `None`, and the applicable force/stress or `N`/`V`/`M`
results are mutually exclusive selections. Selecting a node reports `|u|`,
labelled displacement components, rotation where present, and only nonzero
load/reaction components. Larger transparent node hit targets improve pointer
selection without adding visible labels.

For element 112, `axialForce` and `axialStress` are constant per element and can
be displayed as labels and a diverging color scale.

For element 113, derive `N`, `V`, and `M` diagrams from
`localEndForces = [N1, V1, M1, N2, V2, M2]`. With the current nodal-load-only
model, axial and shear force are constant within the element and bending moment
is linear. The internal-section sign convention must be stated in the GUI help
and verified on the cantilever example. Nodal forces acting on the element must
not be mislabeled as internal section forces.

Unsupported diagrams produce a clear message rather than an empty chart.

## Modal View

The modal view uses the same structural geometry renderer with
`modeShapes[:, modeIndex]` as the displacement field.

Because eigenvectors have arbitrary amplitude and sign:

- normalize by the maximum translational component;
- choose a deterministic sign by making the largest absolute translational
  component positive;
- fall back to the largest rotational component for rotation-only shapes;
- never modify the loaded dataset;
- show the mode number and frequency in Hz;
- provide bounded previous/next buttons and a select box;
- animate the selected shape as `q = cos(2*pi*t/T)` with a fixed illustrative
  two-second cycle at 1x, independent of the physical frequency;
- provide play/pause, 0.5x/1x/2x visualization speed, and a `[-1, 1]` phase
  factor slider;
- keep the peak deformation scale fixed throughout animation and show the
  phase factor, display multiplier, peak scale, and effective scale;
- export the currently visible phase with a caption and machine-readable SVG
  metadata;
- leave tiled comparison and automatic cycling through the mode list for a
  separate enhancement.

Stress and force diagrams are disabled for modal vectors because their absolute
amplitude is arbitrary.

## Transient View

The transient GUI provides:

- node and local-DOF selectors plus the resolved global DOF;
- quantity selector for displacement, velocity, acceleration, load, or
  reaction;
- SVG time-history chart with cursor and current value;
- displacement spectrum when exported;
- time slider and numeric step/time display;
- play, pause, step forward/back, speed, and frame-stride controls;
- a deformed snapshot/animation in the structural viewport.

Use the exported `time` array directly. Do not reconstruct it in JavaScript.
At `1x`, playback traverses the remaining exported time range in five seconds;
`0.5x` and `2x` change that illustrative viewing rate. The frame-stride control
applies to stepping and playback, always includes the final exported sample,
and does not fabricate samples removed by exporter `timeStride`. Changing speed
or frame stride while playing resumes from the currently visible frame.
For a displayed reaction, explain that transient `reactions` is the full
dynamic residual `M*a + K*u - F`; support reactions are its restrained rows.

The renderer calculates element curves only for the current frame. It must not
precompute SVG paths for every time step. Large charts may use a display-only
min/max bucket reduction while preserving the original values for cursor
readout and export.

Transient files can become large. The exporter therefore supports explicit
field selection, DOF selection, and `timeStride`. The GUI displays a visible
warning when a dataset was decimated.

SVG and PNG export synchronously clone the visible transient geometry and chart
before serialization/rasterization. Their caption and metadata identify the
sample number, exported time, history DOF/quantity/value, display scale, frame
stride, and sampling/decimation state.

## SVG and PNG Export

SVG export:

- clone the visible SVG;
- add the SVG namespace and explicit dimensions/viewBox;
- inline all styles needed by the clone;
- add metadata containing dataset title, analysis type, selected mode/time,
  and display scales;
- serialize to a Blob and trigger a download.

PNG export:

- serialize the prepared SVG;
- load it into an `Image` through an object URL;
- draw it to a canvas at the requested pixel size/device scale;
- export with `canvas.toBlob`;
- revoke object URLs after use.

Exported images must contain the legend, scale information, and result context
without including unrelated application controls.

## Offline and Security Requirements

- No network requests at runtime.
- No CDN, external fonts, analytics, service worker, or remote schema.
- No `eval`, `new Function`, or HTML insertion from result data.
- Use `textContent` and SVG DOM APIs for all user-controlled strings.
- Reject non-finite values and unreasonable structural dimensions before
  rendering.
- Place a configurable upper bound on JSON file size and show a helpful error;
  do not freeze the browser silently.
- The application must continue to work when browser networking is disabled.

## Backward Compatibility

- `solveStatic`, `solveModal`, and `solveTransient` remain pure and unchanged.
- `StructFEProblem` does not automatically open a browser or write JSON.
- Existing `plotting=false` default remains unchanged.
- Existing `FEMesh.Plot2DMesh` and `plotTransientResult` continue to work as
  lightweight legacy Octave plots.
- The web exporter is an explicit call, so existing scripts gain no filesystem
  side effects.
- GNU Octave 8.4 remains the minimum supported numerical runtime.
- The browser target is a current Chromium/Edge or Firefox with standard SVG,
  FileReader, Blob, and Canvas support.

## Six-Commit Implementation Sequence

Each commit must be independently reviewable and keep the existing Octave test
suite green. Local links introduced by a commit must resolve within that commit.

### Commit 1 — Add versioned postprocessor data export

Suggested message: `Add versioned postprocessor data export`

- Add `createPostprocessorData.m` and `exportPostprocessorData.m`.
- Define and document JSON schema version 1.
- Export a minimal model without stiffness/mass matrices.
- Implement static, modal, and transient discriminated result objects.
- Implement unit labels, field selection, DOF selection, and explicit transient
  stride metadata.
- Validate dimensions, finiteness, IDs, and analysis-specific fields.
- Add `tests/test_postprocessor_export.m` and register it in the normal test
  runner.
- Verify JSON round trips through `jsonencode`/`jsondecode` without numerical
  or orientation changes.

### Commit 2 — Add the offline SVG structural viewer

Suggested message: `Add offline SVG structural viewer`

- Add the application shell, local CSS, ordered classic scripts, FileReader
  loader, and version-aware validation.
- Add responsive desktop/narrow layout and error handling.
- Implement SVG fit, zoom, pan, reset, selection, labels, and tooltips.
- Implement original geometry plus deformed element 112/113 interpolation.
- Implement automatic/manual scaling and numerical browser tests.
- Confirm the page works from `file://` with networking disabled.

### Commit 3 — Add static results, diagrams, and image export

Suggested message: `Add static postprocessing and export`

- Add static layer controls, supports, loads, reactions, and values.
- Add 112 axial-force/stress display.
- Define and implement 113 `N`, `V`, and `M` diagrams.
- Add legends and independent geometry/arrow/diagram scales.
- Add standalone SVG and PNG downloads.
- Verify the single truss, cantilever, and `ANSYSBeamStatic01.txt` datasets.

### Commit 4 — Add modal result exploration

Suggested message: `Add modal shape exploration`

- Add mode selector and previous/next navigation.
- Implement deterministic sign, normalization, and frequency labels.
- Reuse the common SVG geometry renderer for sinusoidal single-mode playback.
- Add phase scrubbing and illustrative playback speed; stop playback when the
  tab becomes hidden.
- Export the visible frame with modal context and a visible caption.
- Test mode selection, normalization, fixed DOFs, rotation-only and zero modes,
  phase math, scale stability, and export snapshots.

### Commit 5 — Add transient charts and animation

Suggested message: `Add transient result animation`

- Add node/DOF and quantity selection.
- Add SVG history and spectrum charts.
- Add time cursor, slider, snapshot, and animation controls.
- Reuse SVG paths between frames and use a fixed range-wide scale.
- Display decimation/field-omission metadata.
- Test exact time/DOF samples, animation boundaries, pause/resume, and missing
  optional fields.

### Commit 6 — Integrate, document, and verify the postprocessor

Suggested message: `Integrate the HTML FEM postprocessor`

- Add `reference/examples/generate_postprocessor_datasets.m`.
- Add README instructions for export, opening the viewer, and saving images.
- Generate small example datasets for static, modal, and transient analysis.
- Export at least three consistent SVG/PNG illustrations and include only the
  images actually referenced by the offline documentation.
- Add a browser test page and a local-link/no-network check.
- Run all Octave tests, browser tests, desktop/narrow visual QA, keyboard QA,
  and offline QA.
- Document schema evolution and troubleshooting for large transient files.

## Testing Strategy

### Octave Tests

Tests must verify:

- the exporter does not mutate `model` or `result`;
- IDs, DOF labels, supports, and properties are mapped correctly;
- model matrices are absent from JSON;
- static vectors and element results preserve values;
- modal matrix orientation remains `[dof][mode]`;
- transient orientation remains `[dof][timeIndex]`;
- stride and selected-DOF export are explicit and correct;
- invalid/non-finite/inconsistent data fails with stable identifiers;
- failed writes do not damage an existing target file.

### Browser Numerical Tests

`postprocessor/tests/index.html` is a self-contained test harness that reports a
machine-readable final status in the DOM. It tests:

- schema acceptance/rejection;
- linear 112 interpolation;
- Hermite 113 endpoints and end slopes;
- inclined-element coordinate transforms;
- mathematical-Y to SVG-Y mapping;
- automatic scale and zero deformation;
- deterministic modal normalization;
- node/local-DOF mapping;
- exact transient frame/time selection;
- static diagram sign conventions;
- safe text rendering.

The application itself has no Node/npm runtime requirement. Node.js is not
currently available in `PATH`, so the first implementation must not depend on
an npm build or test runner. Browser tests can be opened normally and can also
be driven by a headless installed browser in CI.

### Graphics Smoke Tests

- Load one dataset of every analysis type.
- Assert that required SVG groups and finite path coordinates exist.
- Change every main selector and verify the context label/state.
- Export SVG and PNG and confirm nonempty output.
- Confirm static SVG and PNG exports include the active deflection colour
  legend, its numerical ticks and unit, and the same reaction styles as the
  interactive view.
- Do not use pixel-perfect golden-image comparison; fonts and antialiasing vary
  between browsers.

### Manual Visual QA

- Wide and narrow layouts have no unwanted horizontal scrolling.
- Shared nodes of neighboring elements remain connected.
- Frame curvature and end rotations have correct signs.
- Support/load/reaction symbols do not obscure the structure.
- Static reaction-force arrows and reaction-moment arcs remain distinguishable
  by both colour and line style where they overlap at a clamped node.
- Legends include units and all independent scale factors.
- A static cantilever colours from zero at the clamp to the maximum at the tip,
  and selecting the tip reports the analytical displacement.
- Modal shapes do not change sign between repeated loads.
- Transient axes and deformation scale remain fixed during animation.
- Keyboard navigation and visible focus work.
- The viewer works with browser networking disabled.

## Current Verification Environment

GNU Octave is installed at:

```text
C:\Users\user\Desktop\sft\octave-11.3.0-w64\mingw64\bin\octave-cli.exe
```

Numerical verification command:

```powershell
& 'C:\Users\user\Desktop\sft\octave-11.3.0-w64\mingw64\bin\octave-cli.exe' `
  --no-gui --quiet --norc --no-history `
  --eval "addpath(pwd); run_octave_tests;"
```

At plan creation, GNU Octave 11.3.0 passes all smoke and verification tests.
Node.js/npm is not in `PATH`; it is deliberately not required by the proposed
application. Microsoft Edge is available locally for interactive and headless
browser QA.

## Acceptance Criteria

The postprocessor is complete when:

- one versioned export API accepts all three existing result types;
- `postprocessor/index.html` opens directly from disk with no server/network;
- both 112 and 113 render with physically correct deformation interpolation;
- static results include supports, loads, reactions, stresses, and diagrams;
- static results default to an unscaled resultant-displacement colour field,
  provide focused numerical node inspection, and export its legend;
- modal results provide reliable mode navigation and frequency context;
- transient results provide histories, spectrum, snapshots, and animation;
- the current view exports to standalone SVG and PNG;
- user result data is handled safely as data, never executable markup;
- existing Octave APIs and tests remain compatible;
- exporter tests and browser numerical tests pass;
- at least three reproducible figures are included in the offline reference.

## Main Risks and Mitigations

- **`file://` restrictions:** use classic local scripts and FileReader; do not
  rely on module imports or fetch.
- **Large transient JSON:** explicit field/DOF selection and time stride,
  recorded in metadata; compute only the current animation frame.
- **Frame sign errors:** numerical endpoint/slope tests and analytical
  cantilever diagram tests.
- **Arbitrary modal sign:** deterministic sign canonicalization.
- **SVG Y-axis inversion:** isolate it in one transform and test labels
  separately from geometry.
- **Misleading scales:** keep deformation, arrows, and diagrams independent and
  print every scale in the legend/export.
- **Browser injection through JSON:** strict validation, `textContent`, and DOM
  APIs only; never use result values as HTML.
- **Fragile image tests:** assert numerical SVG geometry and state, not exact
  pixels.
- **Scope creep into a live desktop application:** keep the first version as an
  offline file viewer; a local server or WebView wrapper can be evaluated only
  after the format and renderer are stable.
