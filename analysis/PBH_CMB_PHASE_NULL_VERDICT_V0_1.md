# PBH CMB Phase Null Verdict v0.1

Date: 2026-09-26  
Status: CANDIDATE / FAIL-LOUD / NON-DETECTION / BRANCH ONLY

## Preregistered target

- RA (J2000): 107.576646 deg
- Dec (J2000): -25.907579 deg
- Galactic l: 238.1709218507047 deg
- Galactic b: -7.665782897020605 deg

The sky position was fixed before the CMB phase scan.

## Source stack

The centered analysis used 8 deg x 8 deg TAN cutouts from CDS HiPS2FITS for:

- Planck R3 SMICA CMB component map;
- Planck R3 HFI 100, 143 and 217 GHz;
- WMAP 9-year V and W bands.

All maps were resampled on the same 256 x 256 TAN grid. The comparison path additionally transported maps to a common 40 arcmin beam.

## GREMLIN measurement-derived T^36 binding

The central 4 deg region is transformed into an exact 36D phase cell using a fixed 6 radial x 6 angular Fourier-sector partition.

Each coordinate is the amplitude-weighted circular mean Fourier phase of one fixed sector.

No semantic meanings are assigned to the 36 coordinates.

The comparison algebra is pinned to:

- repository: AdrianLipa90/GREMLIN
- commit: 9fdc42e507d5fe54c13dda975c297173ba0e4eb6
- artifact: gremlin_mcp/holographic_orbital_memory.py
- operations: phase_rms_distance, moire_overlap, project_hologram

## Target result

Common-beam hologram coherence:

- all six maps: 0.8201877317640702
- frequency maps only: 0.8797559912026116

Selected target Moire overlaps:

- SMICA vs WMAP V: 0.5016568891854554
- SMICA vs WMAP W: 0.7198390848454271
- WMAP V vs WMAP W: 0.7260639242292846

SMICA local-center robust z:

- -1.1574638282688676

This is not a large local amplitude outlier.

## Matched Galactic-latitude pilot null

Seven deterministic control fields at the same Galactic latitude were evaluated with the same measurement-derived T^36 mapping and GREMLIN algebra.

Empirical upper-tail fractions:

| Metric | Target empirical upper tail |
| --- | ---: |
| abs(SMICA center robust z) | 0.375 |
| all-map hologram coherence | 0.500 |
| frequency-only hologram coherence | 0.500 |
| SMICA-WMAP V Moire overlap | 1.000 |
| SMICA-WMAP W Moire overlap | 0.750 |

Interpretation:

- the target is not an upper-tail outlier in hologram coherence;
- the target is not an upper-tail local SMICA amplitude anomaly;
- SMICA-WMAP V phase overlap is lower than every pilot control under the preregistered upper-tail test;
- SMICA-WMAP W overlap is also not unusually high.

Therefore the current phase marker does not distinguish the PBH target from matched low-Galactic-latitude controls.

## Fail-loud note

A larger 11-control scan was started on a separate Railway execution path. It reached 9/12 fields before that container was replaced during runtime reconfiguration. No incomplete result is promoted or combined with the completed pilot receipt.

## Verdict

```text
PBH_CMB_PHASE_MARKER_V0_1 = NOT_DISTINGUISHED_FROM_MATCHED_CONTROLS
PBH_DETECTION_CLAIM       = false
TARGET_POSITION_CHANGED   = false
POST_HOC_RECENTERING      = false
CANON_ALLOWED             = false
```

The result rejects promotion of the present CMB phase marker as PBH evidence. It does not test or falsify the independent orbital-dynamics prediction itself.

## Next falsification gate

The next useful gate is not to optimize this marker against the target. It is to hold the target fixed and add independent contamination/noise controls:

1. Planck half-mission or equivalent independent splits where source products permit it;
2. dust/IR foreground correlation at the identical target and matched controls;
3. a larger preregistered matched-latitude null sample;
4. only after those gates, test any new phase statistic with its direction fixed before looking at target rank.

