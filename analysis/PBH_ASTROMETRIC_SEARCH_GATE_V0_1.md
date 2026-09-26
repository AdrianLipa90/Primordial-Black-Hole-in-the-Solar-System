# PBH Astrometric Search Gate v0.1

Status: **NEXT OBSERVATIONAL GATE / STANDARD PHYSICS / CANDIDATE-ONLY**

The CMB phase-marker lane has now been audited against matched-latitude and foreground controls. The current marker does not pass specificity. Standard lensing also shows that an Earth-mass-scale PBH at ~576 AU produces a milliarcsecond, not arcminute, Einstein scale.

## Fail-loud corrections to the historical model

The historical document states approximately:

```text
M = 5e-5 M_sun (~10 Earth masses)
D_L = 576 AU
theta_E ~ 0.31 mas
t_E ~ 1.19 h
```

These quantities are not mutually consistent under standard point-mass lensing.

At 576 AU:

```text
5e-5 M_sun = 16.65 M_Earth
theta_E(10 M_Earth)       ~=  9.36 mas
theta_E(5e-5 M_sun)       ~= 12.08 mas
max astrometric shift      ~=  3.31 -- 4.27 mas
physical Einstein radius   ~=  3910 -- 5045 km
```

An event time near one hour is possible for an effective transverse speed near 1 km/s, but then it is compatible with the ~10 mas Einstein angle, not with 0.31 mas.

## Orbital-state correction

The historical REBOUND excerpt used:

```text
r = 576 AU
vy = 0.001 AU/yr
```

while the same historical model also quotes roughly

```text
a = 430 AU
e = 0.78
```

The vis-viva speed at r~576 AU for a=430 AU is approximately:

```text
v ~= 1.01 km/s ~= 0.213 AU/yr
```

whereas 0.001 AU/yr is only ~4.74 m/s. The old state vector is therefore not dynamically consistent with the quoted orbital scale.

## Sky-track consequence

At 576 AU, the annual parallax amplitude is approximately:

```text
1 AU / 576 AU ~= 5.97 arcmin
```

and the orbital angular motion is of order:

```text
~1.3 arcmin/year
```

Therefore a microlensing search must use a time-dependent sky track. A fixed RA/Dec copied from the 2025 prediction is not sufficient for a current-epoch search.

## Recommended gate

The next physically matched observable is astrometric microlensing:

```text
orbital state + Earth ephemeris
 -> apparent lens track versus time
 -> background-star catalogue
 -> predicted theta_E and centroid-shift field
 -> preregistered search corridor
 -> null corridors with identical stellar-density selection
 -> GREMLIN T36 relational analysis only after standard astrometric residuals are formed
```

Do not tune the PBH trajectory against the stellar residuals. First reconstruct a dynamically consistent orbit from the original TNO fit or state vector, then freeze the track and test it.

Relevant correction receipt:

- `receipts/pbh_microlensing_orbital_scale_correction_v0_1.json`

No main-branch promotion is implied.
