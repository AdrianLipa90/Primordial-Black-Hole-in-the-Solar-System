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


## Planck 30-857 GHz spectral T36 gate

The target-only spectral pass was repeated over the Planck frequency ladder after common smoothing to 60 arcmin.

SMICA-vs-frequency T36 Moire overlaps:

| GHz | T36 Moire |
|---:|---:|
| 30 | 0.277730 |
| 44 | 0.263341 |
| 70 | 0.527174 |
| 100 | 0.418707 |
| 143 | 0.277903 |
| 217 | 0.222267 |
| 353 | 0.357637 |
| 545 | 0.405610 |
| 857 | 0.419877 |

Band means:

```text
30-70 GHz   = 0.3560813173
100-217 GHz = 0.3062925981
353-857 GHz = 0.3943745598
```

This profile is non-monotonic. The 70 GHz channel is a local maximum, 217 GHz is a local minimum, and the high-frequency T36 relation rises again. That shape by itself is not a PBH discriminator.

## Matched-latitude spectral null

A follow-up control used 20 positions at the same Galactic latitude and the same 30-degree target-exclusion rule. The statistic was selected after the target-only spectral pass, so all ranks below are explicitly descriptive/post-hoc rather than formal p-values.

Target:

```text
T36(70 GHz)             = 0.5271735762
low 30-70 mean          = 0.3560813173
mid 100-217 mean        = 0.3062925981
high 353-857 mean       = 0.3943745598
70 GHz peak margin      = 0.1084664381
high minus mid          = 0.0880819617
low minus mid           = 0.0497887192
```

Matched-latitude control distributions:

```text
T36(70 GHz) control median      = 0.9295667053
low-band control median         = 0.7626887685
mid-band control median         = 0.6166737099
high-band control median        = 0.4307489616
70 GHz peak-margin median       = 0.0263674697
high-minus-mid median           = -0.1756026196
low-minus-mid median            = 0.0533736215
```

Descriptive ranks for the target:

```text
T36(70 GHz):       lower = 0.1429, upper = 0.9048
low-band mean:     lower = 0.0952, upper = 0.9524
mid-band mean:     lower = 0.1429, upper = 0.9048
high-band mean:    lower = 0.2857, upper = 0.7619
70 GHz margin:     lower = 0.7619, upper = 0.2857
high minus mid:    lower = 0.9524, upper = 0.0952
low minus mid:     lower = 0.5238, upper = 0.5238
```

The result does **not** support an unusually strong 70 GHz phase marker. The most distinctive scalar in this post-hoc pilot is instead the relatively high `high-minus-mid` contrast, produced primarily by a suppressed 100-217 GHz T36 relation at the target. Because that scalar was defined after inspecting the spectral profile, it must be preregistered and retested on a new control set before any inferential use.

Current interpretation is therefore narrower:

```text
unique PBH phase resonance          NOT ESTABLISHED
70 GHz-specific anomaly             NOT ESTABLISHED
broad cross-frequency coherence     BELOW matched-latitude typical values
mid-band suppression candidate      OBSERVED / POST-HOC
foreground-independent mechanism    NOT ESTABLISHED
```

The next clean gate is to freeze one scalar before looking at any further positions. A natural candidate is

`S = mean(T36_353,545,857) - mean(T36_100,143,217)`,

with the sign and aperture/beam fixed in advance, followed by a fresh blind longitude set not used in either existing control pass.


## Fresh-grid validation of the frozen high-minus-mid statistic

After the previous post-hoc spectral pass, one scalar was frozen before looking at a new sky grid:

```text
S = mean(M353, M545, M857) - mean(M100, M143, M217)
tail = upper
beam = 60 arcmin
field = 8 deg
analysis center = 4 deg
T36 layout = fixed radial6 x angular6
```

A new control grid was offset by 7.5 degrees from every longitude used in the earlier 15-degree grids. The same Galactic latitude and 30-degree exclusion around the PBH longitude were retained.

Target:

```text
S_target = 0.08808196166
```

Fresh controls:

```text
n = 20
mean = -0.2068328290
median = -0.2791850283
std(sample) = 0.2366903398
range = -0.7018729428 .. +0.2895963083
```

Validation result:

```text
standardized target vs control mean = +1.24599
upper add-one rank                  = 0.190476
```

The positive high-minus-mid contrast therefore **did not replicate as an extreme marker** on the fresh control grid. It remains a real measured property of the target cutout, but not a distinctive sky-location discriminator at the present resolution and statistic.

This closes the main exploratory T36 CMB branch with a negative discrimination result:

```text
centered intensity anomaly          NOT SEEN
broad T36 coherence anomaly         NOT SEEN
70 GHz resonance                    NOT ESTABLISHED
high-minus-mid contrast             DOES NOT VALIDATE AS EXTREME
CMB phase evidence for PBH          NOT ESTABLISHED
```

The orbital PBH candidate remains logically separate from this negative CMB marker result.


## Fresh-grid validation of frozen high-minus-mid statistic

After the post-hoc spectral pass, one scalar was frozen before opening a new control grid:

```text
S = mean(M353, M545, M857) - mean(M100, M143, M217)
tail = upper
common beam = 60 arcmin
central analysis field = 4 deg
T36 layout = 6 radial x 6 angular Fourier sectors
```

A fresh longitude grid was offset by 7.5 degrees from the prior 15-degree control grid, at the same Galactic latitude and with the same 30-degree target exclusion.

Result:

```text
target S                 = +0.0880819617
N controls               = 20
control mean             = -0.2068328290
control median           = -0.2791850283
control std (sample)     = 0.2366903398
control range            = -0.7018729428 .. +0.2895963083
standardized target      = +1.2459942001 sigma
upper add-one rank       = 0.1904761905
```

This new grid does not reproduce a strong tail event. The target is above the control mean, but it is not rare enough in this validation to support a detection claim. The proper status of `S` is therefore:

```text
directional contrast     REPRODUCED
strong anomaly           NOT REPRODUCED
formal significance      NOT ESTABLISHED
PBH interpretation       NOT ESTABLISHED
```

This validation is stronger than the original post-hoc rank because the scalar and tail were frozen before the fresh control grid was evaluated, but it is still an internal sky-control validation on the same public Planck products rather than an independent experiment.

## Source identity check

The CDS HiPS product used as the CMB component map, `CDS/P/PLANCK/R3/CMB`, is generated from the Planck R3 SMICA file `COM_CMB_IQU-smica_2048_R3.00_full.fits`. Therefore the internal label `planck_smica` is provenance-consistent and is retained.
