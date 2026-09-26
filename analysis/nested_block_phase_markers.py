#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np


def compact_bits(v: int) -> int:
    out = 0
    bit = 0
    dst = 0
    while v:
        if v & 1:
            out |= 1 << dst
        v >>= 2
        dst += 1
        bit += 1
    return out


def morton_xy(local_index: int) -> tuple[int, int]:
    # HEALPix NESTED face-local index interleaves x in even bits, y in odd bits.
    return compact_bits(local_index), compact_bits(local_index >> 1)


def block_to_grid(values: list[float]) -> np.ndarray:
    n = len(values)
    side = int(round(math.sqrt(n)))
    if side * side != n or side & (side - 1):
        raise ValueError("block length must be a square power of four")
    grid = np.full((side, side), np.nan, dtype=float)
    for local, value in enumerate(values):
        x, y = morton_xy(local)
        grid[y, x] = float(value)
    if not np.isfinite(grid).all():
        raise ValueError("block grid contains non-finite/missing cells")
    return grid


def phase_markers(grid: np.ndarray) -> dict:
    a = np.asarray(grid, dtype=float)
    centered = a - float(np.median(a))
    ny, nx = a.shape
    tapered = centered * np.outer(np.hanning(ny), np.hanning(nx))
    F = np.fft.fftshift(np.fft.fft2(tapered))
    amp = np.abs(F)
    power = amp * amp
    phase = np.angle(F)

    yy, xx = np.indices(a.shape)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    kx = (xx - cx) / nx
    ky = (yy - cy) / ny
    kr = np.hypot(kx, ky)
    theta = np.arctan2(ky, kx)

    # Side=8 WMAP blocks need a wider discrete annulus than side=32 Planck.
    lo = 0.5 / min(nx, ny)
    hi = 0.45
    band = (kr >= lo) & (kr <= hi)
    aw = amp[band]
    pw = power[band]
    if np.sum(aw) == 0 or np.sum(pw) == 0:
        raise ValueError("degenerate Fourier block")

    coherence = float(abs(np.sum(aw * np.exp(1j * phase[band]))) / np.sum(aw))
    hist, _ = np.histogram(phase[band], bins=24, range=(-np.pi, np.pi), weights=aw)
    p = hist / hist.sum()
    p = p[p > 0]
    entropy = float(-np.sum(p * np.log(p)) / np.log(24))

    result = {
        "shape": [int(ny), int(nx)],
        "mean": float(np.mean(a)),
        "std": float(np.std(a)),
        "phase_coherence": coherence,
        "phase_entropy": entropy,
    }
    for harmonic in (1, 2, 3, 6):
        result[f"A{harmonic}"] = float(
            abs(np.sum(pw * np.exp(1j * harmonic * theta[band]))) / np.sum(pw)
        )
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("receipt")
    ap.add_argument("--target-offset", type=int, required=True)
    ap.add_argument("--output", required=True)
    ns = ap.parse_args()

    src = json.loads(Path(ns.receipt).read_text(encoding="utf-8"))
    values = src["values"]
    grid = block_to_grid(values)

    tx, ty = morton_xy(ns.target_offset)
    result = {
        "schema": "PBH_NESTED_BLOCK_PHASE_MARKER_V0_1",
        "status": "ANALYZED",
        "source_receipt": str(ns.receipt),
        "source_config": src.get("config"),
        "target_offset": ns.target_offset,
        "target_xy": [tx, ty],
        "target_value": float(values[ns.target_offset]),
        "markers": phase_markers(grid),
        "interpretation": "LOCAL_NESTED_PARENT_DESCRIPTOR_ONLY",
        "authority": "CANDIDATE_ONLY",
        "canon_allowed": False,
    }
    out = Path(ns.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
