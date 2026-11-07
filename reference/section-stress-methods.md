# Deferred static von Mises beam stress recovery

> Historical record of the implementation available at `f355300`. For the current Python calculation and Docker verification, see [section field diagnostics](EN/08-result-recovery.html#section-fem-title).

**Current application status (2026-11-01):** this kernel and its disabled UI/Worker
modules are retained for future work. The application does not load these scripts,
create a Mises Worker, or offer Mises coloring, ν settings or Mises CSV. The
sections below document the deferred implementation, not active viewer features.
Active coloring uses the existing normal-stress extrema:
`max |σx| = max(abs(σmin), abs(σmax))`, from axial force and bending only.
See the [current viewer guide](EN/08-result-recovery.html#normal-stress-title).

The retained implementation is `shared/mises.js`, designed for the planar and spatial
postprocessors. It evaluates specialized Saint-Venant solutions, not a sectional
finite element model. It neither changes the member stiffness nor adds shear
flexibility, restrained warping, transverse normal stresses, plasticity or local
stress concentrations to the analysis. Homogeneous isotropic material and the
existing small-displacement beam theory are assumed.

## Coordinates and superposition

Use the centroidal principal section coordinates `(y,z)` and the right-handed
member direction `x`. Resultants are `[N,Vy,Vz,T,My,Mz]` with `T=Mx`.

```
sigmaX = N/A + My*z/Iy - Mz*y/Iz
tau = Vy*tauVy + Vz*tauVz + T*tauT
vonMises = sqrt(sigmaX^2 + 3*(tauXY^2 + tauXZ^2))
```

The six stress contributions are combined at the SAME material point, then
maximized over the section. The implementation retains the components and
coordinates of that point. Independent component maxima are never combined.
For trusses only `abs(N)/A` is meaningful in the member model, regardless of the
catalog shape. Local nodal/connection stresses are outside this recovery.

The planar adapter transforms shear and moment into the original profile axes:
`Vy=V*cos(angle)`, `Vz=-V*sin(angle)`, `My=M*sin(angle)`, `Mz=M*cos(angle)`.
The spatial adapter already receives forces in the rolled section axes and must
not apply the section angle a second time. Uniform and linear member loads are
included in the force polynomials; end actions alone do not describe a loaded
member's interior.

## Circle and concentric annulus

For unit `Vy`, outer radius `R`, inner radius `r` (zero for a disk),
`I=pi*(R^4-r^4)/4`, `rho2=y^2+z^2`, define

```
a = (3+2*nu)/(8*(1+nu)*I)
b = (2*nu-1)/(8*(1+nu)*I)
c = -(1+2*nu)/(4*(1+nu)*I)
D = a*R^2*r^2/rho2^2                  (D=0 for a disk, including its centre)
tauXY = a*(R^2+r^2-y^2) + b*z^2 + D*(z^2-y^2)
tauXZ = c*y*z - 2*D*y*z
```

The second unit shear field follows by exchanging `y,z` and the components.
The dipole correction has zero divergence and ensures zero normal traction on
both concentric boundaries. The fields integrate to unit shear and zero torque.
This annular solution is not a disk field with the centre removed, nor a
thin-wall approximation. Unit torque gives `(-z,y)/(Iy+Iz)`.

## Rectangle

Let `y in [-a,a]`, `z in [-b,b]`, `Iz=4*b*a^3/3`. Unit `Vy` has the elementary
part `tauXY=(a^2-y^2)/(2*Iz)`, `tauXZ=0`. Add
`(dpsi/dz,-dpsi/dy)`, with `k=nu/((1+nu)*Iz)` and the explicit series

```
psi = -2*k*b^3/pi^3 * sum(n=1..infinity,
        (-1)^(n+1)/n^3 * (1-cosh(n*pi*y/b)/cosh(n*pi*a/b)) * sin(n*pi*z/b))
```

It satisfies `Laplacian(psi)=k*z` and `psi=0` on all four sides. Thus the
correction is divergence-free, leaves the resultant unchanged, satisfies the
free boundary and supplies the isotropic elasticity compatibility correction.
For `b>a`, the implementation uses the equivalent expansion in
`cos((2*n+1)*pi*y/(2*a))`, with coefficient
`k*4*(-1)^n/((2*n+1)*pi)/lambda^2` multiplying
`[b*sinh(lambda*z)/sinh(lambda*b)-z]`. This avoids a slow expansion across
the wrong aspect ratio. Polynomial parts are summed analytically. Swapping axes
gives the second unit shear field.

Unit-torque stresses come from the derivatives of the rectangular Prandtl
function. With long half-side `a` and short half-side `b`,

```
phi/(G*theta') = b^2-z^2 - 32*b^2/pi^3 * sum(n=0..infinity,
  (-1)^n/(2*n+1)^3 * cosh(lambda*y)/cosh(lambda*a) * cos(lambda*z))
lambda = (2*n+1)*pi/(2*b)
tauXY = dphi/dz; tauXZ = -dphi/dy
```

Normalize by the geometric Saint-Venant `J` to obtain unit torque. If the
rectangle's long axis is `z`, exchange derivatives with their proper signs.
Normal-to-boundary components and corner values are imposed analytically.
Hyperbolic ratios use scaled exponentials to avoid overflow.

The reference square-flexure factors at `nu=0.25` are approximately `0.940`
at the centre and `1.126` at the middle of a side, relative to `1.5*V/A`.
These and rectangular elasticity solutions are discussed in
[Timoshenko & Goodier, as compared in Gruttmann et al.](https://citeseerx.ist.psu.edu/document?doi=9ee6abc7763af187f3ada75007bd1e93cf2e1d19&repid=rep1&type=pdf)
and [Barretta & Barretta (2010)](https://doi.org/10.1016/j.euromechsol.2009.10.008).

## Closed rectangular tubes: nominal stresses

Admission policy: `t/min(b-t,h-t) <= 0.05`. This is a product limit for the
approximation, not a standard, accuracy certificate or buckling check.
Square sections use `h=b` even if a legacy record contains an extra `h` field.

On each straight mid-line wall, integrate
`dq/ds = -t*(Vy*y/Iz_mid + Vz*z/Iy_mid)`. The integration constant closes the
contour with `integral(q/t ds)=0` for shear through the shear centre. Mid-line
inertias make the shear resultants consistent within this thin-wall model.
Double symmetry ensures zero shear-induced torque. Add the Bredt flow of
magnitude `T/(2*Am)`, with sign matching the contour orientation. Convert `q/t`
to both Cartesian shear components using the wall tangent.

Normal stress uses the existing exact catalog area/inertias. Evaluate it on
both wall faces, with nominal shear constant through the thickness. At corner
overlaps each wall's limiting flow is continued separately; never add or average
the two wall vectors. The largest one-sided equivalent value is retained.
This is a nominal combination, not an elastic stress solution near sharp inner
corners or a prediction of a fillet/weld concentration.

## Evaluation and failure handling

All geometry is normalized by a characteristic length. Forces and moments scale
as `F/length^2` and `M/length^3`. A manual member `J` does not replace the
geometric torque normalization; the UI and CSV disclose this distinction.

Series sums double from 16 to at most 4096 terms. Two successive changes must
be at most `1e-5` of the relevant unit-field scale, with zero boundary components
handled analytically. Section sampling doubles from 4 to at most 128 divisions,
including the material interior and all boundaries (polar for disks/rings,
Cartesian for rectangles, both faces of each tube wall). Two successive peak
changes must be at most 0.5%. This is a stabilization test, not a rigorous bound
on the global peak or a guarantee on beam-model accuracy.

Along the member, seed endpoints, quarter points and roots of all force-polynomial
derivatives. Probe interval midpoints and refine locally when the interpolation
error exceeds 0.5% of the member peak. Require two successful peak/curve checks,
including a further global refinement confirmation. Limits are 18 passes and
4097 evaluation stations. An unresolved series, sectional maximum, length curve
or numeric overflow returns `nonconverged`; missing data or applicability failure
returns `unavailable`. Neither is a zero-stress result.

## Runtime contract and compatibility

`create(section,nu)` builds a unit-field kernel; `point(kernel,forces,y,z,wall)`
evaluates a material point; `maximum(kernel,forces)` searches a section;
`member({coeff,...},kernel)` produces length samples and a critical point.
Statuses, components, coordinates, method, nu source and manual-J notes travel
to the UI and CSV. Only static datasets are admitted by the UI.

Worker kernels retain up to eight geometries/nu combinations and bound their
point caches. The result session retains member curves. Camera, deformation
scale and selection changes reuse these curves. A new dataset terminates the
worker and resets all 2D assumptions. Hosted pages use a same-origin worker;
offline `file://` pages use a self-contained Blob worker. No new server endpoint,
JSON field, TXT record or Octave runtime calculation is introduced.

In 3D `nu=E/(2*G)-1`; invalid `nu` excludes frame stress recovery. In 2D the
explicit default assumption is `nu=0.3`, with common and per-member overrides
in the viewer only. Trusses need no nu. A green/yellow/purple stress color never
means that a material or connection has passed a strength assessment.

## Verification commands

The ordinary Octave verification suite includes `test_mises_exports`, which
generates genuine static fixtures for all five shapes in both dimensions.
To prepare just these additional fixtures, run from the repository root:

```
octave-cli --no-gui --quiet --eval "setup(); addpath('tests/octave'); test_mises_exports();"
node scripts/test_mises.cjs
node scripts/test_mises_adapters.cjs
node tests/browser/test_mises_browser.cjs  # verifies UI/scripts/Worker absence
node tests/node/test_normal_stress.cjs
node tests/browser/test_normal_stress_browser.cjs
node tests/browser/test_spatial_normal_stress_browser.cjs
```

The tests also consume existing section fixtures from `run_octave_tests`.
Browser tests use the project Playwright environment (`MKEF_PLAYWRIGHT_MODULE`
where needed). Active normal-stress tests cover both dimensions, all eight planar
shapes and orientations, thick tubes, manual trusses,
interior extrema, cache/reset, scale invariance, CSV, missing inputs, pure torsion
and operation without WebGL. The Mises browser test verifies absence of its UI,
scripts and Worker. Retained Mises mathematical checks include independent
resultant quadrature, boundary traction, differential equilibrium/compatibility,
published square-section factors, circular closed-form maxima and denser peak
sampling.
