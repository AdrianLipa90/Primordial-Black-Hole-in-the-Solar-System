# Legacy PBH FITS WCS audit v0.1

Status: **FAIL-LOUD CORRECTION**

The archived script

`docs/Studies/Metatime-main/simulations/PBH/primordialblackhole.py`

defines the predicted target as ICRS RA/Dec and then calls:

```python
target_x, target_y = wcs.all_world2pix(
    TARGET_COORD.ra.deg,
    TARGET_COORD.dec.deg,
    1,
)
```

The archived Planck cutouts do **not** have an equatorial WCS. Their headers use:

```text
CTYPE1 = GLON-TAN
CTYPE2 = GLAT-TAN
```

Therefore `all_world2pix` expects Galactic longitude/latitude, not the numerical ICRS RA/Dec values.

For the preregistered ICRS target

```text
RA  = 107.576646 deg
Dec = -25.907579 deg
```

the J2000 Galactic coordinate is

```text
l = 238.1709218507047 deg
b =  -7.665782897020605 deg
```

Using the recovered Planck cutout WCS, the correct target is approximately:

```text
FITS 1-based: x = 394.5953267620, y = 420.7185849002
NumPy 0-based: x = 393.5953267620, y = 419.7185849002
```

Passing the ICRS numbers directly as `(l,b)` instead projects to approximately:

```text
NumPy 0-based: x = -2643.86, y = 2690.98
```

which lies far outside the 800x800 cutout.

## Consequence

The legacy point-source proximity check in that script is not a valid test of the preregistered PBH coordinate.

This does **not** invalidate the recovered FITS pixels themselves. The current audit converts ICRS -> Galactic before applying the Galactic TAN WCS. The current fixed-coordinate aperture results are computed at the corrected location and remain candidate-only.

## Correct pattern

```python
gal = TARGET_COORD.galactic
target_x, target_y = wcs.all_world2pix(gal.l.deg, gal.b.deg, 1)
```

The correction is provenance-preserving: the archived script is left unchanged and this audit supersedes only its coordinate-frame interpretation.
