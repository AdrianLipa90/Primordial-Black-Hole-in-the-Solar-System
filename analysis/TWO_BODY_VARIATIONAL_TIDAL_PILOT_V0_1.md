# Two-body variational common-tide pilot v0.1

Status: **FULL-RANK CHANNEL / POINT-MASS GATE FAIL / RELATIONAL_TENSOR_ONLY**

The frozen three-tracer input remained byte-identical:

```text
raw O-C rows      = 736
raw SHA256        = 78492d3031d0e3fd41c883a14117c79c66e5c165d072b2c0293749c12d9ebbf3
nightly equations = 552
nightly SHA256    = 802abdd923b6c0d8e22263dddb7a4d25ebfa932cb8904dd9b973ad5df6946844
```

The physical pilot used the linearized model

```math
r_j=A_j\delta x_j+B_j t+N_j\eta_j+\epsilon_j,
```

where (A_j) is the six-dimensional local orbital-state response,
(B_j) is the five-dimensional symmetric trace-free common-tide response,
and (N_j) contains station/catalog nuisance terms.

## Identifiability

After profiling ([A_j,N_j]), every object retained rank five in (B_{\perp,j}):

| tracer | nights | rank | condition |
|---|---:|---:|---:|
| 2012 VP113 | 24 | 5 | 223.7 |
| 2015 TG387 | 27 | 5 | 832.4 |
| Sedna | 225 | 5 | 153.5 |
| combined | 276 | 5 | **40.73** |

The combined singular values are

```text
84.2034, 18.9226, 12.6016, 4.02121, 2.06749  yr^-2
```

so complementary sky/orbit geometry materially improves the common-tensor
conditioning. The TNO channel is therefore not geometrically rank-dead.

However, each individual tracer loses almost all raw tide-column norm after
profiling local orbit and observational nuisance: retained fractions are
typically (2\times10^{-4}) to (1.4\times10^{-2}). This is a weak-signal
inverse problem.

## Fitted diagnostic tensor

The nuisance-profiled RMS changes only

```text
0.241391 arcsec -> 0.235950 arcsec
```

and the five-parameter tensor accounts for 4.46% of the profiled residual
sum of squares. No likelihood-to-sigma conversion is permitted here.

Its eigenvalues are

```text
-8.85045e-7, +2.90867e-7, +5.94178e-7 yr^-2.
```

GREMLIN v0.2 was run from the exact branch blob
`e6e0083f6c1cae4797642c00551c3ae12d747797`.

Result:

```text
finite                    PASS
symmetric                 PASS
trace-free compatible     PASS
mu > 0                    PASS
point-mass shape          FAIL   residual = 0.7525  (gate <= 0.20)
transverse degeneracy     FAIL   split/mu = 3.958   (gate <= 0.40)

GREMLIN status: RELATIONAL_TENSOR_ONLY
```

Therefore the displayed largest-eigenvector direction
(RA 213.98 deg, Dec +0.05 deg) is **not a source coordinate** and must not
be compared with the sealed historical PBH holdout.

## Boundary and next gate

This v0.1 response uses Sun-only two-body variational dynamics, geometric
geocentric projection, no planetary-gradient derivative, no light-time
derivative, equal nightly weights, and a JPL baseline estimated from
overlapping MPC data.

The next meaningful experiment is the same frozen inverse problem with
major-planet N-body dynamics and finite-difference/variational response
matrices. Only if the full-rank common tensor survives and moves toward the
GREMLIN point-mass eigenshape is a source-direction fit justified.
