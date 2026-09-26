# Student FEM Reference: Offline HTML Tutorial and Octave Labs

## Summary

Create `reference/` as a Russian-language, fully offline, multi-page HTML tutorial. It will explain the implemented FEM workflow at a technical, student-friendly level and connect each concept to small runnable Octave examples that display the actual matrices and results.

No solver behavior or production API will change. `ToDo.md` will be removed because its completed work will be documented in the reference. The existing input-file guide will be migrated into the site and then removed from the project root.

## Reference Site

Use `reference/index.html` as the entry point, with shared offline styling in `reference/assets/reference.css`. Every page will include contents, previous/next navigation, source-code links, and consistent callouts for:

- How to run the example
- What to inspect
- Expected invariant or analytical result
- Common mistakes
- Suggested experiments

The chapters will be:

1. **Model, nodes, elements, and DOFs**
   - Truss and frame degrees of freedom
   - Global DOF numbering
   - Model and result structures

2. **Element matrices**
   - Truss stiffness and mass
   - Euler-Bernoulli frame stiffness and mass
   - Units, symmetry, rank, and rigid-body modes

3. **Coordinate transformation and assembly**
   - Local/global displacement relationship
   - `K_global = T' * K_local * T`
   - DOF mapping and sparse triplet assembly

4. **Boundary conditions and nodal loads**
   - Free/fixed DOF reduction instead of matrix editing
   - Support types
   - Force and moment semantics
   - Static, harmonic, pulse, and step histories
   - Why `Mz` applies to frame nodes but not planar truss nodes

5. **Static analysis**
   - Reduced system solution
   - Expansion to the complete displacement vector
   - Reactions and global equilibrium
   - Bar and cantilever analytical checks

6. **Modal analysis**
   - Generalized eigenvalue problem
   - Natural frequencies and mode shapes
   - Restrained DOFs, mode scaling, and residual checks

7. **Newmark transient analysis**
   - `M*a + K*u = F(t)`
   - Average-acceleration method with `beta=1/4`, `gamma=1/2`
   - Initial acceleration, effective stiffness, time-grid convention
   - Impulse discretization, energy, convergence, and spectrum output

8. **Element-result recovery**
   - Gathering element DOFs
   - Transformation back to local coordinates
   - Truss strain, stress, and axial force
   - Frame end forces and sign conventions
   - Difference between reactions and element forces

9. **Verification and program architecture**
   - Pure solvers and immutable model data
   - Sparse assembly and cached element data
   - Parser validation and headless operation
   - Mapping every completed `ToDo.md` improvement to its implementation and regression test
   - Troubleshooting singular matrices and invalid models

10. **Input-file format**
    - Migrated and revised content from `Руководство по входному файлу.txt`
    - Exact section syntax, units, supported elements, supports, and loads
    - Repeated sections, mixed-mesh restrictions, examples, and error messages

Equations will use ordinary HTML, Unicode, `<sub>`, `<sup>`, tables, and `<pre>` blocks. No CDN, MathJax, web font, build generator, or internet connection will be required.

## Executable Octave Labs

Add assertion-backed functions under `reference/examples/`:

- `example_01_element_matrices.m`
- `example_02_sparse_assembly.m`
- `example_03_static_bar.m`
- `example_04_load_histories.m`
- `example_05_modal_frame.m`
- `example_06_newmark_sdof.m`
- `example_07_result_recovery.m`
- `run_reference_examples.m`

Each example will:

- Run from any working directory by locating the repository root itself.
- Accept an optional `verbose` argument.
- Print the important matrices and intermediate vectors when verbose.
- Return a result structure for further experimentation.
- Assert symmetry, equilibrium, analytical solutions, residuals, or other relevant invariants.
- Use the real production functions rather than reimplementing a separate educational solver.

The public educational entry point will be:

```octave
addpath(fullfile(pwd, "reference", "examples"));
run_reference_examples();
```

## Repository Integration

- Update `README.md` with a concise documentation link and the one-command example runner.
- Remove `ToDo.md` after its completed items are represented in the architecture/verification chapter.
- Remove `Руководство по входному файлу.txt` after its content is migrated and checked against the original.
- Keep filenames and code identifiers in English/ASCII while the explanatory prose remains Russian.

## Six-Commit Implementation Sequence

Each commit must be independently reviewable and must keep the existing Octave
test suite green. Documentation links introduced by a commit must resolve within
that commit; later commits must not be required to repair deliberately broken
navigation.

### Commit 1 -- Add the offline reference foundation

Suggested message: `Add offline FEM reference foundation`

- Create `reference/index.html`, `reference/assets/reference.css`, and the shared
  header, contents, previous/next navigation, callout, table, and code-block
  patterns.
- Add chapters 1--3: model/DOFs, element matrices, and transformation/assembly.
- Add `example_01_element_matrices.m` and `example_02_sparse_assembly.m` with
  optional `verbose` arguments, repository-root discovery, returned structures,
  and assertions.
- Link every discussed production function to its repository source using local
  relative links.
- Verify both examples directly in quiet and verbose modes, open all four HTML
  pages offline, and run the existing `run_octave_tests` suite.

### Commit 2 -- Document loads, constraints, and static analysis

Suggested message: `Add static-analysis tutorials and labs`

- Add chapters 4 and 5: boundary conditions/nodal loads and static analysis.
- Add `example_03_static_bar.m` and `example_04_load_histories.m`.
- Demonstrate free/fixed DOF partitioning, complete-vector expansion, reactions,
  equilibrium, `FL/(EA)`, frame moments, and static/pulse/step/harmonic loads.
- Assert analytical displacement, reaction equilibrium, and exact load-history
  samples; run the four accumulated examples and the existing test suite.

### Commit 3 -- Document modal and transient analysis

Suggested message: `Add modal and Newmark tutorials and labs`

- Add chapters 6 and 7: modal analysis and average-acceleration Newmark analysis.
- Add `example_05_modal_frame.m` and `example_06_newmark_sdof.m`.
- Cover the reduced generalized eigenproblem, mode expansion/scaling, residuals,
  initial acceleration, effective stiffness, time-grid convention, pulse
  discretization, convergence, energy interpretation, and spectrum output.
- Assert modal residuals and restrained DOFs, and compare the SDOF history with an
  analytical solution using explicit tolerances.
- Run the six accumulated examples and the existing test suite.

### Commit 4 -- Complete result recovery and the example runner

Suggested message: `Add result-recovery lab and reference runner`

- Add chapter 8: element-result recovery and sign conventions.
- Add `example_07_result_recovery.m` with truss stress/axial-force and frame
  end-force checks against reactions and analytical values.
- Add `run_reference_examples.m`; it must accept optional `verbose`, run from any
  working directory, return a result structure containing every example result,
  and fail immediately if an assertion fails.
- Run the public one-command entry point in both verbose and quiet modes, followed
  by the existing test suite.

### Commit 5 -- Migrate legacy project documentation into the site

Suggested message: `Migrate architecture and input documentation`

- Add chapter 9: architecture, verification, solver purity, sparse assembly,
  parser validation, headless operation, and troubleshooting.
- Map every completed `ToDo.md` item to the implementing source file and relevant
  regression test. Only then remove `ToDo.md` in this commit.
- Add chapter 10 by migrating and revising every section of
  `Руководство по входному файлу.txt`, including syntax, units, element/support/load
  tables, repeated sections, mixed-mesh restrictions, examples, and diagnostics.
- Compare the new chapter section-by-section with the original, then remove the
  root text guide in this commit.
- Finish index and previous/next navigation now that all ten chapters exist; run
  all examples and the existing test suite.

### Commit 6 -- Integrate and validate the complete reference

Suggested message: `Integrate and verify offline FEM reference`

- Add a quiet automated test for `run_reference_examples(false)` and register it
  in the repository's normal `run_octave_tests` path so the existing CI job runs
  it without a separate workflow.
- Add a reference-site test that checks the index, all chapter files, the shared
  stylesheet, local `href`/`src` targets, navigation consistency, and the absence
  of external runtime assets.
- Update `README.md` with the offline entry page and the one-command Octave runner.
- Perform the final copy-edit and terminology pass, with Russian prose and ASCII
  filenames/code identifiers.
- Run `run_octave_tests` under GNU Octave 11.3 and visually inspect every page at
  desktop and narrow viewport widths with networking disabled.
- Record the verification commands and results in the commit/PR description; do
  not add generated screenshots or reports to the repository.

## Verification and Acceptance

- Add a quiet automated test that runs all reference examples.
- Add a reference-site test that confirms required pages exist, local links resolve, and no external runtime assets are used.
- Run the complete existing Octave test suite to confirm no solver regression.
- Visually inspect the site in a browser at desktop and narrow viewport widths.
- Confirm GNU Octave 11.3 can execute every example headlessly.

The subproject is complete when a student can clone the repository, open `reference/index.html` without internet access, follow the full FEM workflow, inspect real matrices, and run all tutorial examples with one command.

## Assumptions

- The tutorial is primarily Russian and prioritizes implementation understanding over formal derivations.
- Examples use consistent SI units, although the solver itself remains unit-agnostic.
- HTML is the only prose source of truth; equivalent Markdown chapters will not be maintained.
- External references such as the Duke Newmark notes may be linked for further reading but are not required to use the tutorial.
- This work documents existing capabilities and does not introduce new element types or analysis methods.
