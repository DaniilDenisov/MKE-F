# MKE-F postprocessor JSON format, version 3

Version 3 is emitted when a model contains homogeneous MPC equations. It
extends v1 for all three analyses and includes the v2 element-load fields
when static distributed loads are present. Without MPC, v1/v2 selection is
unchanged. Readers must reject versions they do not support.

## Model

- `mpcs`: nonempty array of `{id, depNode, depDOF, rhs, masters}`.
- `id`: consecutive equation number, starting at 1 in input order.
- `depNode`: existing node ID; `depDOF`: 1 = ux, 2 = uy, 3 = thetaZ.
- `rhs`: exactly zero.
- `masters`: nonempty array of `{node, dof, coefficient}`; coefficients are
  finite and at least one is nonzero. Duplicate master DOFs and self references
  are errors. Zero terms do not create graph edges.
- `dependentDOFs`: global DOFs in equation order.
- `independentDOFs`: ascending unconstrained, nondependent global DOFs.

The equation is `u(depNode,depDOF) = Σ coefficient*u(node,dof)`.
A dependent DOF cannot be fixed by an SPC or owned by another MPC. Master
DOFs may be fixed or dependent on other equations. Cycles are errors.
This directed form guarantees independent constraint rows without a dense
rank calculation. Structural mechanisms are checked after elimination.

## Results and sign convention

`C` contains one row per equation: dependent coefficient +1 and master
coefficients `-coefficient`. The original full residual is:

- Static: `r = K*u - F`.
- Transient: `r = M*a + K*u - F` (there is no damping in this solver).
- Modal, per stored eigenvector: `r = K*phi - omega^2*M*phi`.

`mpcMultipliers` contains λ in input equation order, with `mpcForces = C'*λ`
and `r = supportReactions + mpcForces` within numerical tolerance.
Positive λ acts positively on the dependent DOF; each master receives
`-coefficient*λ`. λ is work-conjugate to the dependent displacement: force
for ux/uy, moment for thetaZ. Coefficients relating translations to rotations
carry the corresponding units.

`supportReactions` is zero off SPC DOFs. On an SPC master it is the residual
minus the MPC force, not simply that residual component. Existing static and
transient `reactions` retains the full residual for API compatibility and
must not be displayed as support reactions for MPC models.

All three new fields are mandatory static vectors and modal matrices:

| Field | Static | Modal | Transient |
|---|---|---|---|
| supportReactions | global DOFs | global DOFs × modes | selected DOFs × samples |
| mpcForces | global DOFs | global DOFs × modes | selected DOFs × samples |
| mpcMultipliers | equations | equations × modes | equations × samples |

Transient fields can be omitted with `transientFields`. The default and API
runner include all three. `selectedGlobalDOFs` filters only nodal fields;
all equation multipliers are retained. `timeStride` applies equally to every
history and includes the final sample. Array nesting is preserved for one
equation, mode, selected DOF or sample.

Modal forces refer to the stored eigenvector amplitude. They are not physical
response amplitudes and the viewer's illustrative geometry scale does not
rescale their reported values.

## Case-file format and helper

The repeatable `mpc` section follows an element section:

```text
mpc
1
3,2,0,2,1,2,0.5,2,2,0.5
```

Each record is `depNode,depDOF,rhs,masterCount,node,dof,coefficient,...`.
The example imposes `uy(3) = 0.5*uy(1) + 0.5*uy(2)`.

The axis helper constructs `n·u_k = (1-xi)n·u_i + xi*n·u_j` using the
original geometry. The three nodes must be distinct, and k must be inside
the segment within perpendicular tolerance `1e-9*L`. The largest normal
component selects the dependent DOF (ux on a tie); the other nonzero
component is tried if necessary. The generated coefficients stay fixed
after geometry edits. Recreate the equation to follow a new axis.

End releases use [version 4](schema-v4.md), which also supports MPCs. Nonzero prescribed displacements are not implemented.
