# Three-tracer O-C FWL result v0.4

Status: **PASS DATA / WEAKLY IDENTIFIED QUADRATIC / DO NOT BUILD TIDAL TENSOR FROM THIS ESTIMATOR**

The source O-C dataset is unchanged from v0.2:

```text
N = 736
rows SHA256 = 78492d3031d0e3fd41c883a14117c79c66e5c165d072b2c0293749c12d9ebbf3
station failures = 0
```

The v0.4 calculation corrects the methodological error recorded in
`MPC_HORIZONS_OC_FWL_CORRECTION_V0_1.md`.

For each tangent coordinate, it estimates the quadratic term by the
Frisch-Waugh-Lovell profile

```math
Z=[1,t,S,C],\qquad
\tilde y=(I-P_Z)y,\qquad
\tilde q=(I-P_Z)t^2,
```

```math
\hat c_2={\tilde q^T\tilde y\over\tilde q^T\tilde q}.
```

## Core result

| Object | FWL d2(RA cos dec)/dt2 | FWL d2(Dec)/dt2 | FWL vector norm | surviving t^2 information |
|---|---:|---:|---:|---:|
| 2012 VP113 | -1.9235e-2 | +2.4493e-2 | 3.1143e-2 | 0.618% |
| 2015 TG387 | -1.6337e-2 | -1.5093e-2 | 2.2242e-2 | 0.135% |
| Sedna | +1.4262e-3 | -4.0600e-3 | 4.3032e-3 | 3.035% |

Units of the curvature columns are arcsec/year^2.

The apparently larger FWL coefficients are **not** evidence for a stronger
physical acceleration. They arise because the quadratic regressor is almost
entirely aliased with observing-station/catalog epoch structure. Once the
nuisance space is projected out, only a very small norm of the original
`t^2` design remains.

The corresponding iid-OLS diagnostic standard errors are of the same order
as the coefficients (and are not calibrated significances):

```text
VP113:  RA -0.01924 +/- 0.03465,  Dec +0.02449 +/- 0.02610
TG387:  RA -0.01634 +/- 0.02210,  Dec -0.01509 +/- 0.02681
Sedna:  RA +0.00143 +/- 0.00190,  Dec -0.00406 +/- 0.00185
```

Repeated same-night observations, catalog correlations, and the fact that
the current JPL orbit is itself estimated from overlapping astrometry make
those iid errors optimistic as inference errors.

## Gate

```text
MPC -> exact station/time Horizons O-C       PASS
raw data provenance                         PASS
correct station/catalog FWL                 PASS
quadratic identifiability                   FAIL / very weak
tensor construction from c2 coefficients   BLOCKED
PBH inference                               NOT PERMITTED
```

This does not reject a distant compact perturber. It rejects this particular
polynomial O-C compression as a sufficiently identified estimator.

The next layer must use a physical variational response matrix `B(t)` to a
common external STF tidal tensor and profile each object's six-dimensional
initial state in the same solve:

```math
r_j=A_j\,\delta x_j+B_j\,t+N_j\eta_j+\epsilon_j.
```

Only the component of `B_j` surviving projection through the joint nuisance
space `[A_j,N_j]` may contribute to the common tensor.
