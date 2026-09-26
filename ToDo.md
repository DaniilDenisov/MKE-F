# MKE-F improvement plan

This plan intentionally splits the work into small, reviewable commits. The immediate
goal is to make the existing program run reliably in headless GNU Octave. Numerical
correctness fixes follow once a portable smoke-test loop is available.

GNU Octave 8.4 or newer is the only supported runtime. MATLAB compatibility is not
tested or guaranteed. Avoid gratuitous incompatibility, but do not complicate the
design or tests to preserve MATLAB support.

## Clarified behaviour

- `bcforce_stat` / load type `10` remains a static nodal load in `RunStatic`.
- In `RunTransient`, the same load type intentionally acts for one integration step
  as a rectangular pulse. It must **not** be silently changed into a persistent step
  load. For amplitude `F0` and time step `dt`, its discrete impulse is `F0*dt`.
- Harmonic load type `11` remains `F0*sin(2*pi*f*t)`.
- The solver remains a linear, small-displacement, undamped 2D truss/frame solver
  until damping or nonlinear behaviour is introduced explicitly.
- Headless GNU Octave execution is a supported use case; numerical analysis must
  not require a visible figure or an interactive desktop.

## Deferred provenance notes

This section records useful historical clues, but none of this research blocks the
Octave port or the numerical fixes. Book and ANSYS result recovery can be revisited
later if the original materials become readily available.

### Evidence present in this repository

- `CaseFig11.7p363 MarioPaz.txt` identifies Mario Paz, Figure 11.7, page 363.
- The input-file guide uses that case for modal and transient examples.
- Git commit `4b454b5` says that Euler-Bernoulli beam dynamics still required
  testing and records the known problem that transient analysis modified the global
  stiffness matrix.
- Git history says `ANSYSBeamStatic01.txt` was created for comparison with an
  article on the project website, but neither the expected ANSYS results nor an
  exact article citation are stored in the repository.
- No formal citations for the element matrices or Newmark method were found in the
  current files, deleted guides, or commit messages.

### Possible textbook reference

- Mario Paz and William Leigh, *Structural Dynamics: Theory and Computation*,
  5th ed., Springer/Kluwer, 2004, DOI: 10.1007/978-1-4615-0481-8.
- The fifth edition places "Dynamic Analysis of Plane Frames" on pages 353-379,
  so the case name `Fig11.7p363` is consistent with this edition. The repository
  does not contain enough information to prove the edition or to recover the
  published expected answers.
- Do not treat this case as an external regression oracle unless the original
  edition and published expected values become available.

## Commit plan

### Commit 1 -- Establish GNU Octave compatibility and CI

- Initial observations with portable GNU Octave 11.3.0 on Windows:
  - `octave-cli.exe` runs without installation.
  - The original files emitted invalid UTF-8 warnings; the source files have now
    been converted to UTF-8.
  - The original parser stopped at the abstract method prototypes in
    `FiniteElementStructural.m` with `external methods are only allowed in
    @-folders`; concrete fallback methods now resolve this incompatibility.
- [x] Choose and document GNU Octave 8.4 as the minimum target. Test 11.3.0
      locally and the Ubuntu 24.04 package in CI.
- [x] Add a headless smoke-test command:

      ```sh
      octave --no-gui --quiet --eval "addpath(pwd); run_octave_smoke_tests;"
      ```

- [x] Load and assemble every supplied case in Octave, then exercise static, modal,
      pulse-transient, and harmonic-transient paths.
- [x] Replace the unsupported abstract-method declarations in
      `FiniteElementStructural` with concrete base methods that throw a clear
      "must be overridden" error, while preserving the public API and subclass
      implementations.
- [x] Verify homogeneous arrays of beam and truss handle objects in `FEMesh` on
      Octave 11.3.0. They work without a cell-array workaround.
- [x] Make plotting and diagnostic printing optional so construction and all tests
      work without a graphics display.
- [x] Convert source, Markdown, and Russian documentation files to UTF-8
      without changing program behaviour.
- [x] Avoid APIs unavailable in supported Octave versions. For repeated linear
      solves, use Octave-supported `chol`/`lu` operations or a small helper.
- [x] Add a small smoke-test runner and CI job that executes it in Octave.
- [x] Preserve the original one-argument constructor and default interactive
      behaviour while making headless execution available.

Acceptance: all current example cases at least parse and start in headless Octave,
the portable test command passes, and Octave-specific workarounds are isolated and
documented.

### Commit 2 -- Add a self-contained Octave verification baseline

- [x] Expand the smoke tests into an Octave-only verification suite under `tests/`.
- [x] Add invariant tests for matrix dimensions, finite entries, global stiffness
      symmetry, and truss mass symmetry. Add beam mass symmetry with its known
      correction in Commit 3 so the main branch remains green.
- [x] Add elementary closed-form checks that do not need external books or ANSYS:
      an axial bar (`u=FL/EA`) and a cantilever beam (`v=PL^3/3EI`).
- [x] Check support reactions, axial translational mass, the intentional one-step
      pulse, and harmonic load samples.
- [x] Run every supplied input case as a no-crash regression test. Existing program
      output may be saved as a diagnostic snapshot, but not declared correct merely
      because the current code produced it.
- [x] Use explicit absolute and relative tolerances for numerical comparisons.
- [x] Run the complete verification suite in Octave CI.
- [x] Keep Mario Paz and ANSYS result recovery deferred; do not block subsequent
      commits on unavailable external material.

Acceptance: the Octave suite catches structural failures and verifies basic closed-
form behaviour without proprietary software or hard-to-find references.

## Post-Octave architecture direction

Begin this redesign only after Commit 1 runs the program in Octave and Commit 2
provides a portable safety net. Apply it incrementally; do not combine the whole
redesign with numerical formula corrections in one commit.

Use a hybrid architecture with a functional numerical core and lightweight model
data. The intended direction is:

```matlab
[K, M] = assembleModel(model);
staticResult  = solveStatic(K, model, loadCase);
modalResult   = solveModal(K, M, model);
dynamicResult = solveTransient(K, M, model, loadCase, options);
```

- Treat nodes, connectivity, properties, supports, and load cases as model data.
- Treat displacements, reactions, frequencies, modes, and time histories as result
  data returned by solvers, not mutable state left inside the model.
- Keep matrix assembly, analysis, result recovery, plotting, and file parsing as
  separate responsibilities.
- Numerical functions should receive all required inputs and avoid modifying their
  callers. Repeating analyses in any order must produce the same answers.
- Preserve `StructFEProblem` initially as a thin compatibility facade so existing
  example commands continue to work while it delegates to pure solver functions.
- Keep element-specific mathematics encapsulated. After seeing actual Octave
  behaviour, choose between lightweight value classes and simple type dispatch to
  pure functions such as `truss2DMatrices` and `beam2DMatrices`.
- Do not retain handle-class inheritance only for architectural purity. If the
  abstract hierarchy complicates Octave compatibility without providing useful
  shared behaviour, replace it with the simpler function/struct design.
- Conversely, do not rewrite working element classes merely for style. Make the
  choice using portability, clarity, tests, and assembly performance.

The implementation is distributed across Commits 4, 8, 9, and 10 below.

### Commit 3 -- Correct and test element matrices

- [x] Fix the missing symmetric `-22*L` entry in the beam consistent-mass matrix.
- [x] Express the beam stiffness coefficients directly as `EA/L`, `12EI/L^3`,
      `6EI/L^2`, and so on, avoiding unnecessary multiplication and division by
      `I`.
- [x] Reject zero-length elements and non-positive `A`, `E`, `rho`, or beam `I`.
- [x] Test truss and beam matrices in horizontal, vertical, and inclined
      orientations.
- [x] Test symmetry, rigid-body modes, coordinate-rotation invariance, and total
      translational mass.

Acceptance: element-level tests pass and the assembled beam mass matrix is
symmetric to numerical precision.

### Commit 4 -- Introduce the functional analysis core

- [x] Separate immutable model data from per-run analysis data and returned results.
- [x] Add pure or side-effect-free solver functions for static, modal, and transient
      analysis; let `StructFEProblem` delegate to them as a compatibility facade.
- [x] Stop modifying the stored unconstrained `K` and `M` during an analysis.
- [x] Build loads from a fresh zero vector/time-history for every run.
- [x] Ensure two identical consecutive calls give identical results.
- [x] Ensure static, modal, and transient analyses give the same result regardless
      of call order.
- [x] Return result structures with clearly named fields.
- [x] Move plotting and console formatting outside the numerical solver functions,
      while retaining them as optional convenience behaviour in the facade.

Acceptance: no analysis method leaves `K`, `M`, or a stale `F` in a state that
changes a later analysis, and the numerical core can be called without plotting or
printing.

### Commit 5 -- Replace matrix editing with free-DOF reduction

- [x] Convert support definitions into explicit `fixedDOFs` and `freeDOFs`.
- [x] Solve static and modal problems on `K(freeDOFs,freeDOFs)` and
      `M(freeDOFs,freeDOFs)`.
- [x] Remove artificial zero-frequency modes caused by unit masses on restrained
      DOFs.
- [x] Preserve loads on restrained DOFs when calculating reactions.
- [x] Check static equilibrium: applied loads plus reactions sum to zero.
- [x] Validate duplicate, invalid, or insufficient constraints and report useful
      errors for mechanisms/singular systems.

Acceptance: constrained displacements are exactly zero, modal results contain only
physical free-DOF modes, and reactions balance the applied load.

### Commit 6 -- Make nodal-load semantics explicit

- [x] Replace hard-coded writes to only DOFs 1 and 2 with indexed assembly over the
      element/model DOFs.
- [x] Treat the third frame-node component as a nodal moment `Mz`, not a `Fz`
      translation; update the input guide accordingly.
- [x] Document and test the intentional one-step rectangular pulse used by
      `RunTransient` for load type `10`.
- [x] Decide whether to add a distinct `bcforce_pulse` marker. If added, retain
      backward compatibility for existing case files.
- [x] Add a separate, explicit persistent step-load type rather than changing type
      `10` implicitly.
- [x] Test multiple loads at one node and loads at multiple nodes.

Acceptance: static, pulse, step, harmonic, and nodal-moment meanings are
unambiguous and each has a load-history test.

### Commit 7 -- Verify and improve Newmark transient analysis

- [x] Verify the average-acceleration Newmark equations (`beta=1/4`, `gamma=1/2`)
      against the chosen published reference.
- [x] Define the time convention precisely: column 1 is `t=0`, and load and state
      vectors are evaluated at consistent times.
- [x] Support explicit initial displacement and velocity, and calculate the initial
      acceleration from equilibrium.
- [x] Preallocate displacement, velocity, and acceleration histories with
      `tsNum+1` columns.
- [x] Factor the constant effective stiffness matrix once per run.
- [x] Correct the FFT length, frequency bins, and amplitude normalization using the
      actual number of response samples.
- [x] Return time histories and spectrum data independently of plotting.
- [x] Compare a simple SDOF pulse response with an analytical solution and check
      free-vibration energy conservation for the undamped case.

Acceptance: the time integrator passes analytical SDOF tests and the beam impulse
response is reproducible under time-step refinement.

### Commit 8 -- Repair result recovery

- [x] Replace or rewrite the obsolete `StressCalc` function; it currently calls the
      missing `ElemTransformCalc` and expects an old numeric element format.
- [x] Recover truss axial strain, stress, and axial force from each element's local
      displacement vector.
- [x] Recover frame local end forces and, if useful, axial/shear/moment diagrams.
      Local end forces are returned; diagrams are deferred until member loads are
      supported because nodal loads alone need no additional diagram sampling.
- [x] Add a one-bar stress test and cantilever end-force/reaction tests.

Acceptance: recovery uses the same DOF mapping and transformation conventions as
assembly and matches elementary closed-form solutions.

### Commit 9 -- Simplify and accelerate assembly

- [x] Based on the Octave compatibility results, decide whether elements remain
      lightweight value classes or become structs dispatched to pure element-kernel
      functions. Document the decision briefly. Plain structs and pure functions
      are used; the mutable handle element hierarchy has been removed.
- [x] If classes remain, avoid mutable handle state where value semantics suffice.
      Not applicable after selecting structs; no mutable element handles remain.
- [x] If structs/functions are selected, keep a single documented element interface
      for matrices, DOF mapping, and result recovery.
- [x] Give every element a single element-DOF vector and assemble with indexed
      matrix addition instead of individual scalar assignments.
- [x] Assemble sparse global matrices, preferably from triplet arrays.
- [x] Cache element length/transformation data where appropriate.
- [x] Keep this commit behaviour-preserving and compare all baseline results before
      and after the refactor.

Acceptance: reference results are unchanged within tolerance and larger meshes no
longer allocate dense `N x N` matrices.

### Commit 10 -- Harden input, documentation, and presentation

- [ ] Replace the fixed ten-marker parser loop with an EOF-driven parser.
- [ ] Check `fopen`, close files reliably on errors, and report line numbers for
      malformed input.
- [ ] Correctly append repeated boundary-condition/load sections.
- [ ] Validate node IDs, element connectivity, property counts, and requested plot
      DOFs.
- [ ] Decide explicitly whether mixed truss/frame meshes are supported; reject them
      clearly until heterogeneous assembly is implemented.
- [ ] Make mesh plotting and console printing optional and use equal plot axes.
- [ ] Update the README and Russian input guide with units, DOF meanings, load
      histories, solver assumptions, references, and reproducible examples.

Acceptance: malformed files fail with actionable messages, normal construction can
run headlessly, and all supported input semantics are documented.

## Suggested commit order

Do commits 1-2 first so every later change can be exercised in Octave against a
documented baseline. Commits 3-7 fix numerical correctness. Commit 8 adds
trustworthy post-processing. Commits 9-10 are mostly maintainability, performance,
and usability work and should not be mixed with formula corrections.
