#!/usr/bin/env python3
"""Profiled common-tidal-tensor least-squares solver.

Scientific boundary
-------------------
This module does not ingest synthetic astronomy data and does not claim a PBH.
It solves the linear algebra step after a real orbit/astrometry pipeline has
produced, for each independent tracer:

    residuals r
    local-state Jacobian A
    common-tide Jacobian B
    precision matrix W = C^{-1}

The six local state corrections are profiled out independently for each tracer.
Only the five shared components of a symmetric trace-free tidal tensor are
combined across tracers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class TracerBlock:
    tracer_id: str
    residual: np.ndarray
    local_design: np.ndarray
    tide_design: np.ndarray
    precision: np.ndarray


def vec5_to_tensor(t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=float).reshape(5)
    return np.array(
        [
            [t[0], t[1], t[2]],
            [t[1], t[3], t[4]],
            [t[2], t[4], -t[0] - t[3]],
        ],
        dtype=float,
    )


def tensor_to_vec5(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    if T.shape != (3, 3):
        raise ValueError("tensor must be 3x3")
    if not np.allclose(T, T.T, rtol=0.0, atol=1e-13):
        raise ValueError("tensor must be symmetric")
    if not np.isclose(np.trace(T), 0.0, rtol=0.0, atol=1e-13):
        raise ValueError("tensor must be trace free")
    return np.array([T[0, 0], T[0, 1], T[0, 2], T[1, 1], T[1, 2]])


def _validate_block(block: TracerBlock) -> tuple[np.ndarray, ...]:
    r = np.asarray(block.residual, dtype=float).reshape(-1)
    A = np.asarray(block.local_design, dtype=float)
    B = np.asarray(block.tide_design, dtype=float)
    W = np.asarray(block.precision, dtype=float)

    m = r.size
    if A.shape != (m, 6):
        raise ValueError(f"{block.tracer_id}: local_design must be (n_obs,6)")
    if B.shape != (m, 5):
        raise ValueError(f"{block.tracer_id}: tide_design must be (n_obs,5)")
    if W.shape != (m, m):
        raise ValueError(f"{block.tracer_id}: precision must be square n_obs x n_obs")
    if not all(np.all(np.isfinite(x)) for x in (r, A, B, W)):
        raise ValueError(f"{block.tracer_id}: non-finite input")
    if not np.allclose(W, W.T, rtol=1e-10, atol=1e-14):
        raise ValueError(f"{block.tracer_id}: precision matrix must be symmetric")
    return r, A, B, W


def profiled_precision(
    A: np.ndarray,
    W: np.ndarray,
    *,
    rcond: float = 1e-12,
) -> np.ndarray:
    """Return W projected away from the local-state tangent space."""
    normal_local = A.T @ W @ A
    inv_local = np.linalg.pinv(normal_local, rcond=rcond, hermitian=True)
    W_perp = W - W @ A @ inv_local @ A.T @ W
    return 0.5 * (W_perp + W_perp.T)


def fit_common_tidal_tensor(
    blocks: Iterable[TracerBlock],
    *,
    rcond: float = 1e-12,
) -> dict:
    blocks = list(blocks)
    if len(blocks) < 2:
        raise ValueError("at least two independent tracers are required")

    F = np.zeros((5, 5), dtype=float)
    g = np.zeros(5, dtype=float)
    profiled_chi2_null = 0.0
    cached = []

    for block in blocks:
        r, A, B, W = _validate_block(block)
        Wp = profiled_precision(A, W, rcond=rcond)
        F += B.T @ Wp @ B
        g += B.T @ Wp @ r
        profiled_chi2_null += float(r.T @ Wp @ r)
        cached.append((block.tracer_id, r, B, Wp))

    cov_t = np.linalg.pinv(F, rcond=rcond, hermitian=True)
    t_hat = cov_t @ g
    T_hat = vec5_to_tensor(t_hat)

    profiled_chi2_tide = 0.0
    per_tracer = {}
    for tracer_id, r, B, Wp in cached:
        rr = r - B @ t_hat
        c0 = float(r.T @ Wp @ r)
        c1 = float(rr.T @ Wp @ rr)
        profiled_chi2_tide += c1
        per_tracer[tracer_id] = {
            "chi2_null_profiled": c0,
            "chi2_tide_profiled": c1,
            "delta_chi2": c0 - c1,
        }

    eigvals, eigvecs = np.linalg.eigh(T_hat)
    axis = eigvecs[:, -1]
    for value in axis:
        if abs(value) > 1e-15:
            if value < 0:
                axis = -axis
            break

    Q = 3.0 * np.outer(axis, axis) - np.eye(3)
    mu = float(np.sum(T_hat * Q) / 6.0)
    point_tensor = mu * Q

    norm_T = float(np.linalg.norm(T_hat, ord="fro"))
    point_shape_residual = (
        float(np.linalg.norm(T_hat - point_tensor, ord="fro")) / norm_T
        if norm_T > 0.0
        else float("inf")
    )
    transverse_split = (
        float(abs(eigvals[1] - eigvals[0]) / max(abs(mu), 1e-300))
        if np.isfinite(mu)
        else float("inf")
    )

    return {
        "schema": "PBH_PROFILED_COMMON_TIDAL_FIT_V0_2",
        "status": "FIT_RESULT_NOT_IDENTITY_CLAIM",
        "n_tracers": len(blocks),
        "tensor_vec5": t_hat.tolist(),
        "tensor_3x3": T_hat.tolist(),
        "cov_tensor_vec5": cov_t.tolist(),
        "eigenvalues": eigvals.tolist(),
        "axis_sign_degenerate": axis.tolist(),
        "point_mass_projection": {
            "mu": mu,
            "tensor_3x3": point_tensor.tolist(),
            "frobenius_shape_residual": point_shape_residual,
            "transverse_eigenvalue_split_over_mu": transverse_split,
            "requires_mu_positive": True,
        },
        "fit": {
            "chi2_null_profiled": profiled_chi2_null,
            "chi2_tide_profiled": profiled_chi2_tide,
            "delta_chi2": profiled_chi2_null - profiled_chi2_tide,
            "per_tracer": per_tracer,
        },
        "interpretation_boundary": {
            "delta_chi2_is_not_sigma_without_null_calibration": True,
            "tensor_fit_does_not_identify_PBH": True,
            "mass_distance_degeneracy": "M/R^3 at leading tidal order",
        },
        "canon_allowed": False,
    }


def leave_one_out(blocks: Iterable[TracerBlock], *, rcond: float = 1e-12) -> dict:
    blocks = list(blocks)
    if len(blocks) < 3:
        raise ValueError("leave-one-out requires at least three tracers")
    out = {}
    for idx, block in enumerate(blocks):
        subset = blocks[:idx] + blocks[idx + 1 :]
        fit = fit_common_tidal_tensor(subset, rcond=rcond)
        out[block.tracer_id] = {
            "axis_sign_degenerate": fit["axis_sign_degenerate"],
            "mu": fit["point_mass_projection"]["mu"],
            "shape_residual": fit["point_mass_projection"]["frobenius_shape_residual"],
            "delta_chi2": fit["fit"]["delta_chi2"],
        }
    return out
