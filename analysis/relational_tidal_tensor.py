#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

M_EARTH_PER_M_SUN = 332946.0487


@dataclass(frozen=True)
class OrbitalRecord:
    name: str
    a_au: float
    sigma_a_au: float
    q_au: float
    i_deg: float
    Omega_deg: float
    varpi_deg: float

    @property
    def e(self) -> float:
        return 1.0 - self.q_au / self.a_au

    @property
    def omega_deg(self) -> float:
        return (self.varpi_deg - self.Omega_deg) % 360.0

    def perihelion_unit_vector(self) -> np.ndarray:
        O = math.radians(self.Omega_deg)
        w = math.radians(self.omega_deg)
        inc = math.radians(self.i_deg)
        v = np.array([
            math.cos(O) * math.cos(w) - math.sin(O) * math.sin(w) * math.cos(inc),
            math.sin(O) * math.cos(w) + math.cos(O) * math.sin(w) * math.cos(inc),
            math.sin(w) * math.sin(inc),
        ], dtype=float)
        return v / np.linalg.norm(v)


def unit_from_radec(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    return np.array([
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    ], dtype=float)


def tidal_tensor_shape(u: np.ndarray) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    u = u / np.linalg.norm(u)
    return 3.0 * np.outer(u, u) - np.eye(3)


def tidal_parameter_msun_per_au3(mass_msun: float, distance_au: float) -> float:
    return mass_msun / distance_au**3


def phase_mean_deg(values_deg: list[float]) -> tuple[float, float, float]:
    z = np.mean(np.exp(1j * np.radians(values_deg)))
    mean_deg = float(np.degrees(np.angle(z)) % 360.0)
    R = float(abs(z))
    circ_std_deg = float(np.degrees(np.sqrt(-2.0 * np.log(R)))) if R > 0 else math.inf
    return mean_deg, R, circ_std_deg


def main() -> int:
    src = json.loads(Path("data/real_orbital_constraints_v0_1.json").read_text())
    records = [OrbitalRecord(
        name=row["name"],
        a_au=row["a_au"],
        sigma_a_au=row["sigma_a_au"],
        q_au=row["q_au"],
        i_deg=row["i_deg"],
        Omega_deg=row["Omega_deg"],
        varpi_deg=row["varpi_deg"],
    ) for row in src["sednoids"]]

    peri = {r.name: r.perihelion_unit_vector() for r in records}
    pairwise = {}
    names = list(peri)
    for i,left in enumerate(names):
        for right in names[i+1:]:
            c = float(np.clip(np.dot(peri[left], peri[right]), -1.0, 1.0))
            pairwise[f"{left}__{right}"] = float(np.degrees(np.arccos(c)))

    mean_vec = np.sum(list(peri.values()), axis=0)
    mean_vec /= np.linalg.norm(mean_vec)
    lon = float(np.degrees(np.arctan2(mean_vec[1], mean_vec[0])) % 360.0)
    lat = float(np.degrees(np.arcsin(mean_vec[2])))
    mean_varpi, R, circ_std = phase_mean_deg([r.varpi_deg for r in records])

    old = src["old_work_audit"]
    xyz = np.asarray(old["old_output_xyz_au"], dtype=float)
    old_r = float(np.linalg.norm(xyz))
    m_old_msun = float(old["assumed_perturber"]["mass_Msun"])
    m_old_earth = m_old_msun * M_EARTH_PER_M_SUN

    u_old = unit_from_radec(old["old_output_ra_deg"], old["old_output_dec_deg"])
    Tshape = tidal_tensor_shape(u_old)
    mu_old = tidal_parameter_msun_per_au3(m_old_msun, old_r)
    mu_10earth = tidal_parameter_msun_per_au3(10.0 / M_EARTH_PER_M_SUN, old_r)

    # Published median Fisher sensitivity from Gomes et al. (2023):
    # sigma_M = 0.22 M_earth * (d/400 AU)^3
    sigma_M_earth = 0.22 * (old_r / 400.0)**3

    out = {
        "schema": "PBH_RELATIONAL_ORBITAL_INVERSION_BASELINE_V0_1",
        "status": "BASELINE_ONLY_NO_PBH_INFERENCE",
        "sednoids": {
            r.name: {
                "e": r.e,
                "omega_deg": r.omega_deg,
                "perihelion_unit_vector_ecliptic_like": peri[r.name].tolist(),
            } for r in records
        },
        "current_varpi_relation": {
            "circular_mean_deg": mean_varpi,
            "resultant_R": R,
            "circular_std_deg": circ_std,
            "pairwise_3d_perihelion_angles_deg": pairwise,
            "mean_perihelion_direction_lon_deg": lon,
            "mean_perihelion_direction_lat_deg": lat,
        },
        "source_result": {
            "four_giant_planets_backward_alignment_time_gyr_ago":
                src["source_interpretation"]["backward_alignment_time_gyr_ago"],
            "alignment_varpi_deg":
                src["source_interpretation"]["backward_alignment_varpi_deg"],
            "meaning":
                "These three sednoids do not by themselves provide a present-day residual requiring a current perturber."
        },
        "old_candidate_consistency_audit": {
            "xyz_norm_distance_au": old_r,
            "paper_stated_distance_au": 576.0,
            "assumed_mass_Msun": m_old_msun,
            "assumed_mass_Mearth": m_old_earth,
            "tidal_parameter_Msun_per_AU3": mu_old,
            "tidal_parameter_if_10_Mearth_Msun_per_AU3": mu_10earth,
            "tidal_shape_equatorial": Tshape.tolist(),
            "gomes2023_median_sigma_Mearth_at_candidate_distance": sigma_M_earth,
            "nominal_signal_to_median_sigma_if_10_Mearth": 10.0 / sigma_M_earth,
            "nominal_signal_to_median_sigma_if_assumed_5e-5_Msun": m_old_earth / sigma_M_earth,
            "warning":
                "The Fisher sensitivity is sky-position dependent and is not an observed detection significance."
        },
        "inverse_problem": {
            "measured_object":
                "tidal tensor T = (G M / r^3) (3 u u^T - I) from planetary ephemeris residuals",
            "directly_identifiable":
                ["sky axis u up to sign from tensor eigenvector", "tidal amplitude M/r^3"],
            "not_identifiable_from_tide_alone":
                ["mass M and distance r separately"],
            "degeneracy_breakers":
                ["outer-TNO secular response", "time dependence / proper motion", "direct imaging or microlensing"]
        },
        "authority": "CANDIDATE_ONLY",
        "canon_allowed": False
    }

    Path("receipts").mkdir(exist_ok=True)
    Path("receipts/orbital_relational_baseline_v0_1.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
