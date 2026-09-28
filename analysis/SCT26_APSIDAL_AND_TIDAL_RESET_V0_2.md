# Orbital-relational PBH reset v0.2 — SCT26 + profiled tidal inversion

Status: **ACTIVE RESEARCH / ORBITAL-FIRST / NO PBH DETECTION CLAIM**

This branch supersedes the use of the older `varpi_p ~ 242 deg` relation as a fixed prior.

## 1. Why the reset is required

Siraj, Chyba & Tremaine (2026), *Measuring Apsidal Clustering* (ApJ 1003:76; arXiv:2604.25990), enlarged the stable/metastable distant-TNO sample to 25 objects and introduced a conditional-likelihood estimator designed to cancel uneven survey-footprint selection under the stated population assumptions.

For the SCT26 sample they report:

- simple circular statistic: `1.6 sigma`, direction `58 +/- 84 deg`;
- 2D conditional likelihood: `1.7 sigma`;
- 3D conditional likelihood: `1.9 sigma`;
- 3D preferred perihelion direction approximately `151 deg` with a broad asymmetric interval (`+34/-59 deg`).

Therefore apsidal clustering is currently a weak, directionally broad population-level hint, not a sufficiently rigid datum to define the sky direction of an unseen body.

Primary source:
- https://arxiv.org/abs/2604.25990

The older 2025 best-fit Planet-X result,
`m = 4.4 +/- 1.1 M_earth, a = 290 +/- 30 AU, e = 0.29 +/- 0.13, i = 6.8 +/- 5.0 deg`,
was inferred from the then-current 21-object sample and 300 N-body simulations:
- https://arxiv.org/abs/2410.18170

It remains a useful historical/model comparison, but it is not used as a prior in v0.2.

## 2. Independent reconstruction of the simple SCT26 statistic

The previous branch contained the 20 objects with `a > 170 AU`.  The 2026 sample uses the 25 stable/metastable objects with `a > 160 AU`: the former 20, `2021 LP43`, plus the four newly admitted objects `2017 OF201`, `2021 RR205`, `2023 KQ14`, and `2024 FX26`.

Using the tabulated/published longitudes of perihelion and the equal-weight circular statistic,

```math
z = {1\over N}\sum_i e^{i\varpi_i},\qquad
R=|z|,\qquad
\bar\varpi=\arg z,
```

the reconstructed SCT26 sample gives

```text
N                    = 25
simple mean varpi    = 57.5482 deg
resultant R          = 0.347288
circular std         = 83.3295 deg
```

which reproduces the published `58 +/- 84 deg` simple-statistic result to rounding.

For comparison:

```text
old 20-object (>170 AU) subset:
mean varpi = 61.6272 deg
R          = 0.537041

21-object (>160 AU pre-SCT26) subset:
mean varpi = 66.8391 deg
R          = 0.519788
```

The concentration loss in the 25-object sample is therefore real at the level of the raw circular statistic, before applying the 2026 selection-corrected conditional likelihood.

### Consequence

The former anti-apsidal relation near `242 deg` is reclassified:

```text
HISTORICAL_SAMPLE_CONDITIONAL
NOT_A_FIXED_RELATIONAL_AXIS
NOT_A_PRIOR_FOR_TIDAL_FIT
```

No candidate sky position is emitted from apsidal clustering in v0.2.

## 3. Primary observable: a common external tidal tensor

For a distant source with heliocentric position `R = R n` and a tracer at `r`,

```math
\delta\mathbf a =
GM_X\left[
{\mathbf R-\mathbf r\over |\mathbf R-\mathbf r|^3}
-
{\mathbf R\over R^3}
\right].
```

For `|r| << R`,

```math
\delta\mathbf a = \mathbf T\,\mathbf r + O(r^2/R^4),
```

with

```math
\mathbf T =
\mu(3\mathbf n\mathbf n^T-I),
\qquad
\mu={GM_X\over R^3}.
```

A point-mass tidal tensor is symmetric and trace free, with eigenvalue pattern

```math
(-\mu,-\mu,+2\mu).
```

This eigenvalue structure is a falsifiable signature.  A fitted generic tensor must not be silently projected onto a point-mass tensor and then called a source.

## 4. Profiled astrometric inversion

For object `j`, write the observed-minus-computed astrometric residual vector as

```math
\mathbf r_j =
A_j\,\delta\mathbf x_j +
B_j\,\mathbf t +
\boldsymbol\epsilon_j.
```

Here:

- `delta x_j` is the local six-dimensional correction to that object's osculating state;
- `t` contains the five independent components of a symmetric trace-free tensor;
- `A_j` is the local-state variational design matrix;
- `B_j` is the common-tide variational design matrix;
- `W_j=C_j^{-1}` is the astrometric precision matrix.

Parameterize

```math
\mathbf T =
\begin{pmatrix}
t_1&t_2&t_3\\
t_2&t_4&t_5\\
t_3&t_5&-t_1-t_4
\end{pmatrix}.
```

The nuisance orbital state is profiled out analytically.  Define

```math
W_{\perp,j}
=
W_j -
W_jA_j(A_j^TW_jA_j)^{-1}A_j^TW_j .
```

Then the common-tensor normal equations are

```math
F_T = \sum_j B_j^T W_{\perp,j} B_j,
\qquad
g_T = \sum_j B_j^T W_{\perp,j} r_j,
```

and

```math
\hat t = F_T^{-1}g_T,
\qquad
\mathrm{Cov}(\hat t)=F_T^{-1}.
```

This is the core relational operation: local orbit corrections are removed separately for every tracer, while only the common external field survives the cross-object sum.

Implementation:
- `analysis/profiled_tidal_tensor_fit_v02.py`

## 5. Model hierarchy

The inference order is frozen as:

```text
H0      : no common external tide
H_T     : unrestricted symmetric trace-free common tensor
H_point : T = mu (3 n n^T - I), mu > 0
H_orbit : finite-distance moving point mass, full nonlinear orbit
H_PBH   : compact-object identity from non-dynamical evidence
```

Promotion from `H_T` to `H_point` requires:

1. positive fitted `mu`;
2. small Frobenius point-mass shape residual;
3. near-degenerate transverse eigenvalues;
4. stability of the source axis under leave-one-object-out and leave-one-survey-out;
5. compatible amplitude/direction in independent data layers.

Gravity alone cannot distinguish a PBH from a dark planet with the same mass and trajectory.

## 6. External sensitivity sanity check

Gomes et al. (2023) forecast, using current-era major-planet ranging including Juno, Cassini and Mars orbiters, a median mass uncertainty

```math
\sigma_M \simeq
0.22 M_\oplus
\left({d_X\over400\,\mathrm{AU}}\right)^3.
```

This is a Fisher forecast, not a measured signal or significance.  It demonstrates that the common-tide layer is observationally meaningful and provides a scale check for any future fitted `mu`.

Source:
- https://arxiv.org/abs/2212.09594

## 7. Next data gate

The next scientific result must use real astrometry:

1. fetch MPC ADES observations for the frozen object sample;
2. use a modern baseline dynamical ephemeris;
3. compute real `O-C` residuals and per-observation covariance/weights;
4. propagate variational derivatives for local state and the five common tensor components;
5. fit `H0` and `H_T` with the profiled solver;
6. perform the point-mass eigenvalue gate;
7. calibrate improvement using object/survey/time-block resampling;
8. only then compare the inferred source region with any historical PBH coordinate.

MPC observations API:
- https://docs.minorplanetcenter.net/mpc-ops-docs/apis/get-obs/

No CMB information is used in this branch.

## Status

```text
SCT26 simple-statistic reproduction       PASS
old 242-deg direction as fixed prior      REJECTED
common-tidal-tensor formulation           PASS
profiled linear estimator                 IMPLEMENTED
real astrometric tensor fit               PENDING DATA INGEST
PBH identity                              NOT INFERRED
main                                      UNCHANGED
```
