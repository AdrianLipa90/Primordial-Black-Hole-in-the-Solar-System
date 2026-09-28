# PBH three-tracer STF geometry screen v0.3

Status: **SCREENING FAIL — POINT-MASS GATE**

Branch: `orbital-relational-tidal-v02-20260928`  
Runtime commit: `b39fec120c39c8f40cb950c0e92544519e7032f5`  
Railway deployment: `686f6190-8525-4650-8c93-1ce06e12698a`  
GREMLIN gate: `pbh-orbital-relational-v02-20260928`, content SHA `e6e0083f6c1cae4797642c00551c3ae12d747797`

## Input

The screen uses current MPC ADES observations and JPL Horizons topocentric predictions for:

- Sedna: 630 usable observations;
- 2012 VP113: 55;
- 2015 TG387 / Leleakuhonua: 51.

The historical PBH coordinate and the historical 576 AU scale are **not used** anywhere in the fit.

For each tangent coordinate the nuisance model is

```math
r_{j}(t)=b_{g(j)}+eta_1\tau+\frac12\beta_2\tau^2+\epsilon_j,
```

where `b_g` is a station-by-astrometric-catalog fixed effect.  The fitted angular curvature is then mapped through the actual mean line-of-sight geometry,

```math
\ddot\theta_\alpha \simeq e_\alpha^T T n,
\qquad
\ddot\theta_\delta \simeq e_\delta^T T n,
```

into a generic symmetric trace-free 3x3 tensor.

This is a **screening approximation**, not full six-dimensional variational orbit determination.

## GREMLIN point-mass gate

GREMLIN v0.2 requires, among other checks:

```text
point_mass_shape_residual <= 0.20
transverse_eigenvalue_split_over_mu <= 0.40
```

Observed screen:

| subset | point-mass shape residual | transverse split / mu | condition proxy | verdict |
|---|---:|---:|---:|---|
| ALL | 0.8027 | 4.6629 | 460.4 | FAIL |
| Gaia-only | 0.8310 | 5.1745 | 1328.5 | FAIL |

The leading eigen-axis happens to be similar in the two cuts:

```text
ALL       RA 322.94 deg, Dec +5.70 deg
Gaia-only RA 326.88 deg, Dec +7.33 deg
```

That axis is **not promoted to a source direction**.  The tensor itself fails the point-mass eigenshape, the geometry is poorly conditioned because the three tracers occupy a restricted sky region, and the Gaia-only system is even more ill-conditioned.

## Per-object curvature diagnostics

```text
ALL, arcsec/yr^2
Sedna       RA +0.001393 +/- 0.001986   Dec -0.002606 +/- 0.001915
2012 VP113  RA -0.020426 +/- 0.034723   Dec +0.022801 +/- 0.026301
2015 TG387  RA -0.014247 +/- 0.022378   Dec -0.017863 +/- 0.026761

Gaia-only, arcsec/yr^2
Sedna       RA -0.000119 +/- 0.001151   Dec -0.002683 +/- 0.001350
2012 VP113  RA -0.094191 +/- 0.044093   Dec +0.053323 +/- 0.031428
2015 TG387  RA -0.014328 +/- 0.019367   Dec -0.018189 +/- 0.029163
```

No Gaussian significance is assigned.  The quoted standard errors are conditional regression errors and do not incorporate all orbit-fit or catalog systematics.

## Interpretation

The current three-tracer residual field is **not consistent enough with the leading-order tide of one distant point mass to pass the GREMLIN shape gate**.

This does **not** rule out an unseen body.  The screen has only six tangent components for five tensor parameters, the objects are geometrically clustered, JPL orbit solutions are themselves fitted using MPC astrometry, and the local six-dimensional orbit state has not yet been variationally profiled.

The next valid step is therefore a frozen, high-quality, wide-sky TNO control sample selected independently of residual sign or the fitted axis.  Only after improving the tensor conditioning should the point-mass gate be rerun.

```text
three-tracer generic STF fit       PASS
ALL point-mass eigenshape gate     FAIL
Gaia-only point-mass gate          FAIL
source axis promotion              BLOCKED
PBH identity                       NOT INFERRED
historical-coordinate holdout      UNTOUCHED
main                               UNCHANGED
```
