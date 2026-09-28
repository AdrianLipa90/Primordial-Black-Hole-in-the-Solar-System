# O-C FWL correction note v0.1

Status: **METHODOLOGY CORRECTION / APPEND-ONLY**

## What was wrong

The v0.2 O-C diagnostic attempted to remove station and astrometric-catalog offsets by alternating demeaning of the dependent residual coordinate only, followed by a quadratic regression on the original time basis.

That is not the Frisch-Waugh-Lovell operation for the requested model.

For

```math
y = c_0 + c_1 t + c_2 t^2 + S\beta + C\gamma + \epsilon,
```

where `S` and `C` are station/catalog fixed effects, a correct profile of `c_2` must project **both** `y` and the quadratic regressor `t^2` through the same nuisance space

```math
Z=[1,t,S,C].
```

Then

```math
\tilde y=(I-P_Z)y,\qquad
\tilde q=(I-P_Z)t^2,
```

and

```math
\hat c_2=
{\tilde q^T\tilde y\over \tilde q^T\tilde q}.
```

Demeaning only `y` can bias the quadratic coefficient whenever survey/station membership is correlated with epoch, which is exactly the case for long TNO arcs.

## Invalidated derived quantities

The following outputs remain stored for provenance but are **not valid scientific results**:

- `PBH_MPC_HORIZONS_OC_RESIDUALS_V0_2.objects[*].station_catalog_profiled_curvature`;
- the leave-one-group ranges derived from that incorrect profiled coefficient;
- `PBH_OC_TO_STF_TIDAL_SCREEN_V0_1`, because its six input curvature components came from the invalid profiling step;
- any interpretation of its apparent `-2:+1:+1` signed eigenshape.

## Still valid

The underlying observational transport remains valid:

- MPC ADES acquisition;
- exact observing-station / exact-time Horizons matching;
- the 736 matched O-C rows in the three-tracer run;
- raw O-C values and station/catalog group offsets;
- their provenance hashes.

The correction changes the regression/statistical layer, not the source observations.

## Replacement

`runtime/mpc_horizons_oc_residuals_v03_fwl.py` implements the correct FWL profile

```text
Z = [intercept, linear time, station dummies, catalog dummies]
q = t^2
y_perp = (I - P_Z)y
q_perp = (I - P_Z)q
c2 = (q_perp'y_perp)/(q_perp'q_perp)
```

and additionally reports the surviving quadratic identifiability after nuisance projection.

No PBH/source inference is permitted until the v0.3 FWL receipt is available and a new tensor screen is built from it.
