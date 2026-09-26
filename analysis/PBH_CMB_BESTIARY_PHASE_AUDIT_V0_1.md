# PBH CMB Bestiary Phase Audit v0.1

Status: **CANDIDATE / PROVENANCE-CLOSED / CURRENT MARKER SPECIFICITY FAIL**

Target is fixed before CMB inspection:

```text
RA  = 107.576646 deg
Dec = -25.907579 deg
l   = 238.1709218507047 deg
b   = -7.665782897020605 deg
```

## Data lane

The analysis uses centered CDS HiPS2FITS cutouts derived from Planck PR3 SMICA, Planck 100/143/217 GHz, and WMAP 9-year V/W. The repository branch stores source manifests, SHA-256 receipts, analysis code, and null-test receipts. No target coordinate was shifted to maximize a signal.

## Bestiary route

```text
FERRET
 -> HUMMINGBIRD
 -> BAT / HOUND / OWL / SPIDER
 -> BELZEBUB
 -> GREMLIN
```

FERRET acquires read-only astronomical sources. BAT extracts Fourier and phase descriptors. HOUND checks target-local anomaly. OWL audits provenance and epistemic status. SPIDER checks cross-mission relations. BELZEBUB applies the false-positive firewall. GREMLIN aggregates the 36D/Moire candidate.

## Target result

At a common 40 arcmin beam the target gives:

```text
|SMICA center z|                    = 1.203404
GREMLIN mean pairwise T36 overlap   = 0.821166
GREMLIN hologram coherence          = 0.911708
SMICA-WMAP V image correlation      = 0.771899
SMICA-WMAP W image correlation      = 0.796296
WMAP V-W image correlation          = 0.991076
```

The high WMAP V-W consistency is real as a cross-mission/map relation, but is not by itself a PBH marker.

## Matched-latitude null

A deterministic control grid was evaluated at the same Galactic latitude, with 20 longitudes and a 30 degree exclusion around the preregistered target longitude.

Add-one empirical upper-tail ranks:

```text
|SMICA center z|                    p_emp = 0.285714
GREMLIN mean pairwise T36 overlap   p_emp = 0.238095
GREMLIN hologram coherence          p_emp = 0.238095
SMICA-WMAP V correlation            p_emp = 0.904762
SMICA-WMAP W correlation            p_emp = 0.904762
```

These are coarse empirical ranks from this matched-latitude control set, **not formal survey p-values**.

## Bestiary verdict

BAT: the phase/harmonic descriptors are measurable and reproducible.

HOUND: the fixed target is not a strong local intensity outlier.

OWL: source/provenance lane is clean enough for candidate analysis; physical interpretation remains unpromoted.

SPIDER: Planck/WMAP relations exist, but comparable or stronger relations occur at controls.

BELZEBUB: reject any claim that the present CMB phase result is a PBH detection.

GREMLIN: the current 36D sector/Moire descriptor **fails the specificity gate**. Its target coherence is above much of the control set, but not rare enough to separate the PBH coordinate from matched-latitude sky structure.

This is a negative result for the **current marker definition**, not a proof against the orbital PBH hypothesis.

## Next gate

Do not tune the target location. Keep the coordinate fixed. Improve only the statistic:

1. enlarge the foreground-matched null set;
2. add dust/synchrotron templates and mask-sensitive controls;
3. test half-mission or independent Planck splits;
4. replace the coarse 6 x 6 Fourier-sector cell with a preregistered multiscale orbital/fractal representation;
5. require the new statistic to beat the same null without target-specific tuning.

Relevant receipts:

- `receipts/pbh_cmb_centered_phase_v0_2.json`
- `receipts/pbh_cmb_matched_latitude_pilot_null_v0_1.json`
- `receipts/pbh_cmb_matched_latitude_null_v0_2.json`
- `receipts/pbh_cmb_bestiary_handoff_v0_1.json`

No main-branch promotion is implied.
