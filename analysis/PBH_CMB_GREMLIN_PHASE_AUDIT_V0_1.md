# PBH CMB GREMLIN / FERRET Phase Audit v0.1

Status: `CANDIDATE / EXPLORATORY / NO PBH DETECTION CLAIM`

Target fixed before this CMB pass:

- RA (J2000): `107.576646 deg`
- Dec (J2000): `-25.907579 deg`
- Galactic: `l = 238.1709218507047 deg`, `b = -7.665782897020605 deg`

## Provenance lane

FERRET retrieved fixed 8-degree TAN cutouts through the CDS HiPS2FITS service for:

- Planck PR3 CMB component map;
- Planck PR3 100 GHz;
- Planck PR3 143 GHz;
- Planck PR3 217 GHz;
- WMAP 9-year V;
- WMAP 9-year W.

All six returned FITS payloads and were hashed before analysis.

The phase transform is explicitly measurement-derived and candidate-only:

```text
central 4-degree image
 -> Fourier transform
 -> fixed k band 0.02..0.20 cycles/pixel
 -> 6 radial x 6 angular sectors
 -> amplitude-weighted circular phase per sector
 -> exact 36D phase cell
```

No semantic meaning is assigned to the 36 coordinates.

## Center test

After smoothing to a common 40 arcmin beam, robust center-vs-annulus z values are:

| map | robust z |
|---|---:|
| Planck 100 | -0.814 |
| Planck 143 | -0.627 |
| Planck 217 | +0.014 |
| Planck SMICA | -1.157 |
| WMAP V | -0.821 |
| WMAP W | -0.850 |

No map shows a strong centered intensity excursion at the preregistered coordinate.

## GREMLIN T36 marker

For the five frequency maps (Planck 100/143/217 + WMAP V/W), the common-beam T36 projection gives:

- holographic global coherence: `0.8797577579`;
- mean pairwise Moire overlap: `0.7600473809`;
- median pairwise Moire overlap: `0.7848278972`;
- minimum pairwise overlap: `0.5046256615`;
- maximum pairwise overlap: `0.9540580541`.

The strongest pair is Planck 143 vs Planck 217:

```text
Moire overlap = 0.9540580541
phase RMS distance = 0.3078449802 rad
```

The weakest frequency pair is Planck 217 vs WMAP W:

```text
Moire overlap = 0.5046256615
phase RMS distance = 1.2550264042 rad
```

## Matched-latitude null

A separate exploratory null used 20 control positions at the same Galactic latitude, with a 30-degree exclusion around the PBH longitude.

Primary descriptive statistic:

`median pairwise Moire overlap across the five frequency maps after common-beam smoothing`.

Observed:

```text
target median = 0.7848278972
control median = 0.8190411582
control mean = 0.8108895717
control range = 0.4787489767 .. 0.9274298383
empirical upper-rank fraction = 0.7142857143
```

Therefore the target is **not an unusually high T36 cross-frequency coherence location** in this matched-latitude pilot.

This empirical rank is not a formal p-value because the scalar summary statistic was selected after an earlier exploratory look at the target.

## Foreground gate

The dedicated foreground receipt compares SMICA with Planck 30 GHz and 353 GHz on 20 matched-latitude controls.

Target:

```text
SMICA vs 30 GHz image correlation = 0.3031611281758779
SMICA vs 30 GHz T36 Moire        = 0.5558586037940171
SMICA vs 353 GHz image correlation = 0.05988295207932019
SMICA vs 353 GHz T36 Moire         = 0.2812448168791079
```

Descriptive upper-rank fractions:

```text
|SMICA x 30 image| = 0.6666666666666666
SMICA x 30 T36    = 0.14285714285714285
|SMICA x 353 image| = 0.8095238095238095
SMICA x 353 T36    = 0.7619047619047619
```

The 30 GHz T36 overlap is elevated relative to this small control set, while the corresponding image correlation is not exceptional. This is a foreground diagnostic that should be tested against synchrotron/free-free structure before any PBH interpretation.

## Current verdict

```text
FITS provenance                    PASS
fixed-coordinate retrieval         PASS
cross-mission retrieval            PASS
centered intensity anomaly         NOT SEEN
T36 cross-frequency anomaly        NOT SEEN in matched-latitude pilot
foreground-independent PBH marker  NOT ESTABLISHED
PBH detection claim                NOT SUPPORTED by this CMB pass
```

The current result does not falsify a dark compact object at the orbital coordinate; it says only that this CMB/phase-marker pass does not provide a statistically distinctive PBH signature.

## Next falsification gate

1. add Planck 30/44/70 and 353/545/857 frequency channels;
2. fit the phase marker as a function of frequency;
3. test whether the elevated 30 GHz T36 relation follows a low-frequency Galactic foreground law;
4. add Planck half-mission / split-map component products for instrumental reproducibility;
5. only then define a preregistered phase statistic for a larger blind sky null.

Authority:

```text
canon_allowed = false
physical_PBH_detection = false
main_branch_mutation = false
```
