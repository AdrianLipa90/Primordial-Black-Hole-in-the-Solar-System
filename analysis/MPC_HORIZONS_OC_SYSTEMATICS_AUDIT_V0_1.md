# MPC–Horizons O−C systematics audit v0.1

Status: **PASS DATA PATH / SYSTEMATICS-DOMINATED CURVATURE / NO SOURCE CLAIM**

## Data path

The v0.1 residual pipeline compares current MPC ADES astrometry with current JPL Horizons astrometric RA/Dec at the exact MPC observing station and exact observation time.

Frozen pilot:
- 2012 VP113: 55 usable MPC positions, 55 Horizons matches, 0 station failures.
- 541132 Leleakuhonua (2015 TG387): 51 usable MPC positions, 51 Horizons matches, 0 station failures.
- combined: 106/106 residual rows.
- canonical residual-row SHA256: `458ade608adc688e2856ce4c2d928d481495b8d315ce758040c4dc5d08b48bb2`.
- runtime revision: `parser-ra-dec-fixed-r2`.

The current Horizons solution is not independent of the MPC astrometry used to estimate the JPL orbit. This layer is therefore a residual/systematics diagnostic, not an independent detection experiment.

## Raw residual curvature

A quadratic tangent-plane profile

```math
r(t)=c_0+c_1(t-t_0)+c_2(t-t_0)^2
```

was used only as a diagnostic. The reported angular second derivative is `2 c2`.

Raw fits (arcsec/year^2):

| object | d2(RA cos dec)/dt2 | d2(Dec)/dt2 | vector norm |
|---|---:|---:|---:|
| 2012 VP113 | +2.9175e-3 | -8.0920e-3 | 8.60e-3 |
| 2015 TG387 | +4.7135e-3 | -4.2623e-3 | 6.36e-3 |

The matching sign pattern is superficially interesting but is not stable against observational nuisance structure.

## Catalog/station structure

Examples:
- VP113 USNOA2 (26 positions): mean residual about -0.101 arcsec in RA*cos(dec), +0.296 arcsec in Dec.
- VP113 Gaia3 (21 positions): mean about -0.016 arcsec in RA*cos(dec), -0.034 arcsec in Dec.
- TG387 station 950 / 2MASS (11 positions): mean about -0.623 arcsec in RA*cos(dec), +0.207 arcsec in Dec.
- TG387 Gaia2 (3 positions): residual scale is only about 0.019 arcsec radial RMS.

These offsets are large compared with the curvature left after nuisance profiling.

## Two-way fixed-effect profile

For each coordinate, alternating station and catalog means were profiled out before repeating the quadratic fit.

Result (arcsec/year^2):

| object | raw vector norm | station+catalog profiled vector norm | suppression |
|---|---:|---:|---:|
| 2012 VP113 | 8.60e-3 | 2.2646e-4 | ~38x |
| 2015 TG387 | 6.36e-3 | 2.9063e-5 | ~219x |

Profiled components:

```text
2012 VP113:
  d2(RA cos dec)/dt2 = -1.1787e-5 arcsec/yr^2
  d2(Dec)/dt2        = +2.2615e-4 arcsec/yr^2

2015 TG387:
  d2(RA cos dec)/dt2 = -2.6248e-5 arcsec/yr^2
  d2(Dec)/dt2        = -1.2480e-5 arcsec/yr^2
```

The common raw (+RA,-Dec) pattern therefore does not survive nuisance profiling.

## Leave-one-group-out leverage audit

The profiled quadratic is not uniformly identifiable after removal of groups that anchor the time baseline.

VP113:
- full profiled norm: `2.26e-4` arcsec/yr^2;
- typical leave-one-station-out solutions remain around `2e-4..3e-4`;
- removing station 645, which anchors the early arc, raises the norm to `4.22e-3`;
- removing Gaia3 similarly raises the norm to `7.51e-3`.

TG387:
- full profiled norm: `2.91e-5` arcsec/yr^2;
- typical leave-one-station-out solutions remain around `2e-5..3e-5`;
- removing early station 695 raises the norm to `7.14e-4`;
- removing Gaia3 raises it to `8.01e-4`.

This is temporal leverage / design-matrix fragility. It prevents interpreting the residual quadratic as a stable physical acceleration.

## Verdict

```text
exact-station exact-time O-C transport     PASS
raw same-sign curvature                    OBSERVED
survives station/catalog profiling         NO
leave-one-group stability                  FAIL / endpoint-leverage sensitive
common external-tide candidate             NOT ESTABLISHED
PBH inference                              NOT PERMITTED
```

The correct next step is to add a third long-baseline tracer (Sedna), repeat the identical nuisance-projected audit, and only then consider a 3D common-tensor diagnostic. A real common-tide fit still requires variational derivatives against the orbit state; polynomial O-C curvature is only a screening layer.
