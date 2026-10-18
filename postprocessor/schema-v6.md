# MKE-F postprocessor JSON format, version 6

Version 6 is emitted when the model has a nonempty section catalog, for static,
modal and transient analyses. It composes with all v1–v5 features. Linear loads
still require static analysis; releases require the complete v4 DOF registry.
Older datasets retain their original version and behavior. Older readers reject v6.

## Sections and assignments

`model.sections` is an array of `{id, name, type, dimensions}` records. IDs are
unique positive integers, names are nonempty Unicode strings. Types and dimensions:

| Type | Dimensions | Geometry |
|---|---|---|
| `rectangle` | `b,h` | Solid rectangle |
| `circle` | `d` | Solid circle |
| `square_tube` | `b,t` | Square tube with outer side b and uniform wall t |
| `rectangular_tube` | `b,h,t` | Rectangular tube with outer width b, height h and uniform wall t |
| `circular_tube` | `d,t` | Circular tube with outer diameter d and uniform wall t |
| `channel` | `b,h,tw,tf` | Web at negative z, opening towards positive z |
| `tee` | `b,h,tw,tf` | Centered web, top flange at positive y |
| `ibeam` | `b,h,tw,tf` | Centered web and equal top/bottom flanges |

All dimensions use model length units and must be positive finite values.
For open profiles `tw < b`; `tf < h` for a tee, `2*tf < h` otherwise.
Tubes require `2*t < b` (square), `2*t < min(b,h)` (rectangular), or
`2*t < d` (circular). Square and rectangular tubes have sharp corners.
No fillets, flange slopes, torsion, eccentricity or out-of-plane DOFs.
Areas and centroidal second moments are calculated from nonoverlapping rectangles
and the parallel-axis theorem, or analytically for a circle or annulus.

Each element has `sectionId` (0 means manual A/I) and `sectionAngle` (0,90,180,270).
Geometry is authoritative: exported `properties.area` and, for frames,
`properties.momentOfInertia` must match it. Material E and density remain on elements.
Changing geometry affects stiffness and mass in every analysis type.

At 0°, x points from node 1 to node 2, y is the in-plane left normal to x, and
z points out of the drawing. Section axes form a right-handed frame. At angle θ,
`yθ = cosθ*n + sinθ*Z`, `zθ = -sinθ*n + cosθ*Z`.
The section origin is its centroid and lies on the element axis.
With `Iy=∫z²dA`, `Iz=∫y²dA`, the in-plane bending inertia is
`I=Iz*cos²θ+Iy*sin²θ`; a point's in-plane coordinate is `η=y*cosθ-z*sinθ`.

## Static stresses

`analysis.elementResults[].sectionStress` is an empty array for an unassigned
element, otherwise `{minimum, maximum, positions}`. Each extremum is
`{value, xi, eta}`: signed normal stress, normalized position x/L, and the
in-plane extreme-fiber coordinate. `positions` is an increasing array containing
0, 1, stationary points of both extreme-fiber stresses, and zeros of M. Repeated
positions within numerical tolerance are merged. Equal extrema select the first
position, then the lower fiber. Modal/transient results omit sectionStress.

Tension is positive. Frames use `σ=N/A-M*η/I`; trusses use `σ=N/A`.
For summed local distributed-load endpoints `(qx1,qy1),(qx2,qy2)`, L the length,
and f the recovered local end forces, polynomials in ξ=x/L are:

```text
N(ξ) = -f[0] - qx1*L*ξ - (qx2-qx1)*L*ξ²/2
M(ξ) = -f[2] + f[1]*L*ξ + qy1*L²*ξ²/2 + (qy2-qy1)*L²*ξ³/6
```

Use the two extreme η values to bound stresses across each section. Global extrema
use both endpoints and all real stationary points in [0,1], independently of plot
resolution. Legacy mean `axialStress` and `StressCalc` retain their meanings.
The browser reconstructs stresses from exported forces and loads; shared tests
compare its extrema with Octave. Missing sections are unavailable, never zero.

## Case text

Optional repeatable blocks after the elements, before or after other conditions:

```text
sections
1
{"id":1,"name":"Rectangle 100 x 200","type":"rectangle","dimensions":{"b":0.1,"h":0.2}}
element_sections
1
1,1,90
```

Each section is one JSON line, retaining names including commas and Unicode.
Assignments are `elementId,sectionId,angle`. References are resolved after reading
the file, before assembly. Duplicate section IDs, duplicate element assignments,
invalid dimensions/orientations and missing references are errors. Existing numeric
element records remain required; their A/I are replaced by assigned geometry.
The editor writes synchronized A/I for compatibility with older manual workflows.

## UI and verification

The editor's Sections panel supports create/copy/edit/delete and previews.
Select a member and use Assign to selected element in a catalog card, or select
a section and orientation in the element row. Each card lists assigned members. Editing updates all uses; deletion
of a referenced section is blocked; detaching keeps current A/I. These operations
preserve loads and support undo/redo. Rerun after edits; saved result JSON is an
immutable snapshot of its calculation, not a live model.

Static results offer σmin/σmax diagrams and a selected-member section at x/L.
The section color scale stays fixed over the selected member. Positive stress is
red, negative blue. Section axes can be toggled; SVG section export is separate
from the model's existing SVG/PNG export.

Example: `examples/cases/CaseSectionBeam.txt`. Run `test_sections` in Octave,
then `node scripts/test_sections.cjs` for numerical parity and format tests.
Browser acceptance: `node scripts/test_sections_browser.cjs` using Playwright.
The Docker recipe is `docker/browser-tests.Dockerfile`; mount the repository at
`/app` after Octave has generated the result fixtures in `output/`.
