# Relational PBH Inversion from Orbital Deviations v0.1

Status: `STARTED / SOURCE-GROUNDED / NO PBH DETECTION CLAIM`

## Reset of the inference chain

The previous 2025 PBH calculation is not used as an input location for this branch.

That work assumed a PBH orbit first and then propagated a synthetic ensemble of ten TNOs. Therefore the old RA/Dec is treated only as a historical output to audit, not as evidence and not as a prior in the new inverse problem.

The new chain is:

```text
real observations
  -> baseline Solar-System dynamics
  -> post-fit orbital residuals
  -> common external tidal tensor
  -> relational source direction + M/r^3
  -> outer-TNO response
  -> mass/distance degeneracy breaking
  -> candidate trajectory
```

## Source-grounded sednoid layer

Huang & Gladman (2024) back-integrated Sedna, 2012 VP113 and Leleakuhonua, including orbital clones. They found a tight apsidal alignment around varpi ~ 200 deg approximately 4.5 Gyr ago under dynamics dominated by the four giant planets, with only weak modification from the Galactic tide.

This is important because it means the three canonical sednoids cannot simply be treated as unexplained present-day deviations caused by a still-existing object. Their present apsidal configuration has a viable no-current-perturber dynamical history.

Therefore the relational PBH inversion must begin with residuals that survive a fitted baseline model, not with present-day angle clustering alone.

## Relational observable

For a distant body whose position changes negligibly over a planetary ranging baseline, the leading perturbation is tidal:

```math
a_X(x) = G M_X / r_X^3 * [3 (x dot u) u - x].
```

Define

```math
T_X = mu_X (3 u u^T - I),  mu_X = G M_X / r_X^3.
```

This is the natural relational object.

Its eigensystem has one longitudinal eigenvector `u` and two transverse directions. Planetary residuals therefore constrain the **directional tidal tensor** first. They do not independently determine mass and distance.

The primary inverse unknown is thus

```text
(u, M/r^3)
```

rather than

```text
(M, r, a, e, i, Omega, omega, anomaly)
```

all at once.

## Data layers

### Layer A — inner Solar-System ephemerides
Cassini, Juno, Mars ranging and planetary astrometry.

Purpose: infer or bound the external tidal tensor.

### Layer B — Pluto / large-TNO astrometry
Long-baseline RA/Dec residuals.

Purpose: add leverage on direction and reject models that only repair one planetary subsystem.

### Layer C — detached TNO secular dynamics
Sedna, 2012 VP113, Leleakuhonua and later well-determined objects.

Purpose: test whether the same source tensor predicts their secular residual structure and, if possible, break the mass-distance degeneracy.

## Fail-closed rule

A PBH candidate is not emitted unless the same relational source improves independent orbital data layers without requiring incompatible source directions or tidal amplitudes.

CMB is excluded from this inference branch.
