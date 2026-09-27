# Primordial-Black-Hole-in-the-Solar-System

## PBH FITS proxy staging

This branch is a provenance-first staging surface for the pre-registered PBH sky target:

- RA (J2000): `107.576646 deg`
- Dec (J2000): `-25.907579 deg`

It does not treat any CMB, IR or X-ray feature as a PBH detection by default.

Large/full-sky mission products stay at their official archives. This repository stores immutable source manifests, exact request parameters, small tangent-plane cutouts when appropriate, SHA-256/FITS-header receipts, and phase-marker analysis code.

GitHub Actions artifacts are used as temporary binary FITS proxy storage so full mission maps do not need to be committed to Git history.

The current fixed-coordinate lane requests a 12-degree field from official IRSA/LAMBDA sources for Planck PR3 100/143/217 GHz and WMAP 9-year V/W, then computes exploratory Fourier-phase markers. Markers remain candidate-only until they survive cross-frequency, cross-mission and foreground controls.
