# PBH CMB Phase36D foreground gate v0.4

The exact pre-registered PBH coordinate was tested again after projecting IRAS 100 micron onto the recovered Planck 100/143/217 GHz cutouts and removing one linear dust template per Planck band.

A 128×128-pixel patch was used so that the IRAS-covered region contains a useful within-field null population. The marker is exactly 36-dimensional: twelve multiscale Fourier/rotation/amplitude descriptors for each of the three Planck frequencies. The 335 null patches overlap spatially, so their empirical percentiles are descriptive controls, not independent-sample p-values.

The target shrinkage-Mahalanobis distance is 18.645, at the 34.33 percentile of the within-field controls. RMS standardized displacement is 0.767 and the largest single standardized feature is 1.47 sigma-equivalent units relative to that empirical feature distribution. There is therefore no multivariate 36D outlier at the fixed target.

Low-spatial-frequency cross-phase coherence after the IRAS projection is 0.989 (100–143; 76.1 percentile), 0.913 (100–217; 47.2 percentile), and 0.939 (143–217; 48.4 percentile). Mid-band values lie between the 39.7 and 42.1 percentiles. None is an unusual common-sky phase lock at the target.

Fresh live GREMLIN/Terminal36D consumed this third-pass packet and routed it through OWL + SPIDER + HOUND. Its internal PhaseNav vector was [0.19787, 0.04144, 0.23425, 0.20042], with attention index 2. These are GREMLIN routing/attention quantities only.

**Verdict:** current recovered Planck frequency cutouts do not show a significant fixed-coordinate amplitude or 36D phase marker after foreground control. This does not reject the orbital compact-object/PBH hypothesis; it closes this particular CMB evidence lane unless an independent map (SMICA HM1/HM2 or WMAP) reproduces a signal.
