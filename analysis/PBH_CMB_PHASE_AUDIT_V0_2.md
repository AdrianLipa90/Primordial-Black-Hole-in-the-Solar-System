# PBH CMB fixed-coordinate phase audit v0.2

Status: **candidate-only observational audit**. No PBH detection claim.

The pre-registered target is fixed at RA 107.576646 deg, Dec -25.907579 deg. The recovered Planck cutouts are 800×800, 0.025 deg/pixel, in thermodynamic K_CMB, and are centered on the rounded request (107 deg, -26 deg). All target statistics below use the exact WCS target pixel, not the cutout center.

## Recovered source maps

The NOEMA ZIP knowledge plane contained previously generated Planck cutouts for 100 GHz Release 3, 143 GHz Releases 1/2/3, and 217 GHz Release 3. Each recovered member was hash-checked against the selector index. The R2/R3 143-GHz maps are extremely stable: RMS difference is 1.88e-6 K, about 0.78% of the R3 field standard deviation.

## Fixed-coordinate amplitude test

Using a 0.2-deg aperture and a 0.5–1.0-deg local background annulus, robust target z values are:

| Band | robust z | empirical percentile |
| --- | ---: | ---: |
| 100 GHz | 0.065 | 57.6% |
| 143 GHz | 0.115 | 57.2% |
| 217 GHz | 0.466 | 70.8% |

The target is therefore not a local brightness outlier in any of the three bands.

## Phase-marker null test

For a fixed 192×192-pixel patch around the exact target and 250 deterministic control patches:

| Pair | target cross-phase coherence | percentile |
| --- | ---: | ---: |
| 100–143 | 0.684 | 33.6% |
| 100–217 | 0.480 | 4.0% |
| 143–217 | 0.716 | 16.8% |

The 100–217 value is unusually **low**, not unusually high. This does not look like a frequency-invariant common-sky phase marker.

The angular Fourier A2 marker is elevated at 100 GHz (88.0 percentile) and 143 GHz (89.2 percentile), but not at 217 GHz (58.0 percentile). That is an exploratory structural marker, not a detection statistic.

## Foreground control

The IRAS 100-micron map was reprojected against the Planck cutouts. Pixel correlations with IRAS are 0.709, 0.813 and 0.925 for 100, 143 and 217 GHz respectively. A simple linear IRAS template removes about 50.3%, 66.1% and 85.5% of the map variance respectively.

After the IRAS template subtraction the exact target robust z becomes:

- 100 GHz: -0.190
- 143 GHz: -0.190
- 217 GHz: -0.285

The strong increase of IRAS correlation toward 217 GHz and the disappearance of the small positive target residual after dust projection are consistent with strong Galactic foreground contamination in this low-latitude field.

## Live GREMLIN / PhaseNav 36D audit

Fresh NOEMA↔AUX guard: **ACTIVE**.

First pass routed through OWL + SPIDER + HOUND with PhaseNav/Terminal36D state consumed. The foreground-residual pass routed to OWL only. These are internal attention/routing quantities; the PhaseNav R values are not sigma values and are not astrophysical probabilities.

## Current verdict

**No significant fixed-coordinate CMB excess is found in this test.**

The result does not rule out a compact object at the orbital candidate coordinates; it specifically says that these recovered Planck frequency cutouts do not currently provide a clean CMB phase/amplitude signature for it. The next independent gate is WMAP and/or Planck component-separated half-mission maps, followed by a microlensing/astrometric test that does not depend on Galactic foreground emission.
