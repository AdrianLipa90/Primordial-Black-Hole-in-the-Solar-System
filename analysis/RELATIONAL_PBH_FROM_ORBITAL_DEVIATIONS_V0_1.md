# Relational PBH inference from orbital deviations — v0.1

Status: **CANDIDATE / ORBITAL-FIRST / NO CMB INPUT**

This branch restarts the PBH problem from orbital structure only. The legacy 576 AU PBH state is not used as an input.

## Step 1 — stable/metastable distant TNOs

The input sample is the published stable/metastable TNO set with semimajor axis above 170 AU from Siraj, Chyba & Tremaine (2024). For these objects the eccentricity vector is reconstructed from the tabulated heliocentric elements.

The reproduced mean eccentricity-vector magnitude is

```text
|<e>| = 0.4452877847
```

which reproduces the paper's reported value of approximately 0.44.

## Step 2 — relational anti-apsidal direction

The eccentricity-weighted circular mean longitude of perihelion is

```text
<varpi>_e = 62.187746 deg
```

so the directly inferred anti-aligned latent apsidal direction is

```text
varpi_rel = 242.187746 deg
```

Sensitivity checks:

```text
equal-weight anti-varpi       = 241.627203 deg
survival-weighted anti-varpi  = 240.400591 deg
stable-only anti-varpi        = 246.153691 deg
```

Thus the relational axis is stable at the several-degree level under these weighting choices.

The full 3D antipode of the sample mean eccentricity vector is

```text
ecliptic lon = 243.806781 deg
ecliptic lat = +16.256911 deg
RA           = 244.829098 deg
Dec          = -4.922828 deg
```

This is **not** the instantaneous sky position of the perturber. It is an anti-apsidal relational axis.

## Step 3 — first dynamical-strength relation

Using the observed transition scale near `a_crit = 170 AU`, a leading-order balance between giant-planet-induced precession and the quadrupolar field of a distant perturber gives the scale relation

```text
mu_rel = m_p / [a_p^3 (1-e_p^2)^(3/2)]
       ~= 2.69087e-7 M_earth / AU^3
```

This produces, for example,

| e_p | 250 AU | 290 AU | 330 AU | 380 AU |
|---:|---:|---:|---:|---:|
| 0.00 | 4.20 | 6.56 | 9.67 | 14.77 |
| 0.20 | 3.95 | 6.17 | 9.10 | 13.89 |
| 0.29 | 3.69 | 5.75 | 8.48 | 12.94 |
| 0.40 | 3.24 | 5.05 | 7.44 | 11.37 |

Values are Earth masses and are **leading-order scale estimates**, not a substitute for a full N-body posterior.

## Current relational object

At v0.1 the PBH candidate is therefore not a point but a constrained manifold:

```text
apsidal direction: varpi_p ~ 242 deg
strength:          m_p/[a_p^3(1-e_p^2)^(3/2)] ~ 2.69e-7 Mearth/AU^3
orbital phase:     unresolved
node:              unresolved by current significant clustering
instantaneous RA/Dec: unresolved
```

The next gate is to use the observed `varpi(a)` relation and long-term stability to collapse this manifold in `(m,a,e)`, then only afterwards infer current position along the orbit.
