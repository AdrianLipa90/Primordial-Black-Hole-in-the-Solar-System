# PBH CMB Standard-Lensing Scale Check v0.1

Status: **PHYSICAL SENSITIVITY CHECK / FAIL FOR DIRECT CMB LENSING INTERPRETATION**

This note checks only the standard weak/point-mass gravitational-lensing scale for the preregistered Solar-System PBH candidate.

## Candidate scale

The historical PBH material uses a distance near

```text
D_L ~ 576 AU
```

and contains two mass descriptions that are not numerically identical:

```text
~10 Earth masses
5e-5 solar masses
```

The second quantity is approximately **16.65 Earth masses**, not 10 Earth masses. Both are therefore carried as separate sensitivity cases rather than silently reconciled.

## Einstein angle

For a source effectively at cosmological distance relative to a Solar-System lens,

```text
theta_E ~= sqrt(4 G M / (c^2 D_L))
```

At 576 AU:

```text
M =  5.0 M_Earth   -> theta_E ~=  6.62 mas
M = 10.0 M_Earth   -> theta_E ~=  9.36 mas
M = 5e-5 M_Sun
  = 16.65 M_Earth  -> theta_E ~= 12.08 mas
```

For unresolved astrometric microlensing, the maximum centroid displacement is

```text
delta_theta_max = theta_E / (2 sqrt(2))
```

which gives approximately:

```text
5 M_Earth      -> 2.34 mas
10 M_Earth     -> 3.31 mas
16.65 M_Earth  -> 4.27 mas
```

## Consequence for CMB maps

The direct point-mass lensing scale is therefore of order **10 milliarcseconds**. Planck/WMAP temperature maps operate on arcminute and larger scales; the PBH CMB experiments in this branch additionally used a common 40 arcmin beam for cross-mission comparisons.

Thus the current degree/arcminute-scale CMB phase structures cannot be interpreted as a spatially resolved standard gravitational-lensing image of a 5--17 Earth-mass PBH at ~576 AU.

This does **not** rule out the orbital PBH hypothesis. It rules out one interpretation of the current CMB marker:

```text
degree-scale / arcminute-scale phase feature
!=
directly resolved standard PBH Einstein-scale lensing
```

Any claim that the CMB phase marker is physically produced by the candidate PBH therefore requires a separately specified and independently tested mechanism beyond the standard point-mass lensing response.

## Better matched observable

The same calculation points to a physically better search channel: **astrometric microlensing of background stars**, where a few-milliarcsecond centroid displacement is the relevant scale.

That search should use the preregistered orbital sky track rather than retuning the position after looking at the astrometric data.

Authority:

```text
standard_lensing_scale_check = true
pbh_existence_claim = not decided
nonstandard_mechanism_claim = not validated
canon_allowed = false
```
