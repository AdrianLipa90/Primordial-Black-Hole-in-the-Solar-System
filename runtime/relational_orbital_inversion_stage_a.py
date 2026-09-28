#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import urllib.parse
import urllib.request

import numpy as np

OBJECTS = [
    "2000 CR105",
    "Sedna",
    "2004 VN112",
    "2010 GB174",
    "2012 VP113",
    "2013 FT28",
    "2013 RA109",
    "2013 SY99",
    "2013 UT15",
    "2014 SR349",
    "2015 RX245",
]

N_NULL = 5_000_000
NULL_CHUNK = 100_000
N_BOOT = 200_000
SEED = 20260928


def fetch(name: str) -> dict:
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api?" + urllib.parse.urlencode(
        {"sstr": name, "full-prec": "1", "cov": "mat", "cd-epoch": "1"}
    )
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "GREMLIN-PBH-relational-orbits/0.1 read-only"},
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())
    orbit = data["orbit"]
    elems = {e["name"]: e for e in orbit["elements"]}
    def val(k): return float(elems[k]["value"])
    def sig(k):
        x = elems[k].get("sigma")
        return None if x is None else float(x)
    return {
        "query": name,
        "fullname": data["object"].get("fullname"),
        "orbit_id": orbit.get("orbit_id"),
        "epoch": float(orbit["epoch"]),
        "rms": float(orbit["rms"]) if orbit.get("rms") is not None else None,
        "n_obs_used": int(orbit["n_obs_used"]) if orbit.get("n_obs_used") else None,
        "data_arc_days": float(orbit["data_arc"]) if orbit.get("data_arc") else None,
        "a": val("a"), "a_sigma": sig("a"),
        "e": val("e"), "e_sigma": sig("e"),
        "q": val("q"), "q_sigma": sig("q"),
        "i": val("i"), "i_sigma": sig("i"),
        "Omega": val("om"), "Omega_sigma": sig("om"),
        "omega": val("w"), "omega_sigma": sig("w"),
        "M": val("ma") if "ma" in elems else None,
        "M_sigma": sig("ma") if "ma" in elems else None,
        "source_url": url,
    }


def circ_mean_deg(values_deg: np.ndarray) -> tuple[float, float]:
    z = np.mean(np.exp(1j * np.deg2rad(values_deg)))
    return float(np.rad2deg(np.angle(z)) % 360.0), float(abs(z))


def node_order(Omega_deg: np.ndarray, inc_deg: np.ndarray) -> tuple[float, float]:
    weights = np.sin(np.deg2rad(inc_deg))
    z = np.sum(weights * np.exp(1j*np.deg2rad(Omega_deg))) / np.sum(weights)
    return float(np.rad2deg(np.angle(z)) % 360.0), float(abs(z))


def circular_delta_deg(a: float, b: float) -> float:
    return float(np.rad2deg(np.angle(np.exp(1j*np.deg2rad(a-b)))))


def robust_leave_one_out(varpi: np.ndarray) -> dict:
    best = None
    for k in range(len(varpi)):
        keep = np.arange(len(varpi)) != k
        mu, R = circ_mean_deg(varpi[keep])
        row = {"remove_index": int(k), "mu_deg": mu, "R": R}
        if best is None or row["R"] > best["R"]:
            best = row
    return best


def null_test(obs_Raps: float, obs_Rnode: float, inc_deg: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED)
    weights = np.sin(np.deg2rad(inc_deg))
    wsum = float(np.sum(weights))
    ge_aps = ge_node = ge_joint = 0
    total = 0
    obs_joint = obs_Raps * obs_Rnode
    for _ in range(0, N_NULL, NULL_CHUNK):
        m = min(NULL_CHUNK, N_NULL-total)
        varpi = rng.uniform(0.0, 2.0*np.pi, size=(m, len(inc_deg)))
        Omega = rng.uniform(0.0, 2.0*np.pi, size=(m, len(inc_deg)))
        Raps = np.abs(np.mean(np.exp(1j*varpi), axis=1))
        Rnode = np.abs(np.sum(weights[None,:]*np.exp(1j*Omega), axis=1)/wsum)
        joint = Raps * Rnode
        ge_aps += int(np.count_nonzero(Raps >= obs_Raps))
        ge_node += int(np.count_nonzero(Rnode >= obs_Rnode))
        ge_joint += int(np.count_nonzero(joint >= obs_joint))
        total += m
    return {
        "draws": total,
        "seed": SEED,
        "null": "uniform independent varpi and Omega; measured inclinations held fixed",
        "p_add_one_apsidal": (ge_aps+1)/(total+1),
        "p_add_one_node": (ge_node+1)/(total+1),
        "p_add_one_joint_product": (ge_joint+1)/(total+1),
        "observed_joint_product": obs_joint,
        "warning": "bias-unaware geometric null; not a survey-selection-corrected significance",
    }


def bootstrap_cluster(rows: list[dict], exclude_index: int) -> dict:
    keep = [r for j,r in enumerate(rows) if j != exclude_index]
    vp = np.array([(r["Omega"]+r["omega"])%360.0 for r in keep])
    Om = np.array([r["Omega"] for r in keep])
    inc = np.array([r["i"] for r in keep])
    rng = np.random.default_rng(SEED+1)
    n = len(keep)
    idx = rng.integers(0,n,size=(N_BOOT,n))

    z = np.mean(np.exp(1j*np.deg2rad(vp[idx])),axis=1)
    center = np.angle(z)
    mu0 = np.angle(np.mean(np.exp(1j*np.deg2rad(vp))))
    d = np.angle(np.exp(1j*(center-mu0)))
    unwrapped = mu0+d
    qvp = np.rad2deg(np.quantile(unwrapped,[0.025,0.5,0.975]))%360.0

    weights = np.sin(np.deg2rad(inc[idx]))
    zO = np.sum(weights*np.exp(1j*np.deg2rad(Om[idx])),axis=1)/np.sum(weights,axis=1)
    centerO = np.angle(zO)
    w0 = np.sin(np.deg2rad(inc))
    muO0 = np.angle(np.sum(w0*np.exp(1j*np.deg2rad(Om)))/np.sum(w0))
    dO = np.angle(np.exp(1j*(centerO-muO0)))
    unwrappedO = muO0+dO
    qO = np.rad2deg(np.quantile(unwrappedO,[0.025,0.5,0.975]))%360.0

    return {
        "draws": N_BOOT,
        "apsidal_cluster_deg_95_bootstrap": qvp.tolist(),
        "anti_apsidal_perturber_deg_95_bootstrap": ((qvp+180.0)%360.0).tolist(),
        "node_cluster_deg_95_bootstrap": qO.tolist(),
        "note": "object bootstrap only; does not include discovery-selection bias",
    }


def pairwise_graph(rows: list[dict]) -> dict:
    nodes=[]
    for r in rows:
        varpi=(r["Omega"]+r["omega"])%360.0
        nodes.append({**r,"varpi":varpi})
    edges=[]
    for i in range(len(nodes)):
        for j in range(i+1,len(nodes)):
            dv=abs(circular_delta_deg(nodes[i]["varpi"],nodes[j]["varpi"]))
            dO=abs(circular_delta_deg(nodes[i]["Omega"],nodes[j]["Omega"]))
            # Relation score only: 1 is identical in both angular channels.
            score=0.5*(math.cos(math.radians(dv))+1.0)*0.5*(math.cos(math.radians(dO))+1.0)
            edges.append({"i":i,"j":j,"d_varpi_deg":dv,"d_Omega_deg":dO,"relation_score":score})
    strength=np.zeros(len(nodes))
    for e in edges:
        strength[e["i"]]+=e["relation_score"]
        strength[e["j"]]+=e["relation_score"]
    order=np.argsort(-strength)
    return {
        "nodes":[{"index":i,"name":n["fullname"],"varpi_deg":n["varpi"],"Omega_deg":n["Omega"],"i_deg":n["i"]} for i,n in enumerate(nodes)],
        "edges":edges,
        "node_strength":[{"index":int(i),"name":nodes[int(i)]["fullname"],"strength":float(strength[i])} for i in order],
    }


def main() -> int:
    rows=[fetch(name) for name in OBJECTS]
    varpi=np.array([(r["Omega"]+r["omega"])%360.0 for r in rows])
    Omega=np.array([r["Omega"] for r in rows])
    inc=np.array([r["i"] for r in rows])

    mu_aps,R_aps=circ_mean_deg(varpi)
    mu_node,R_node=node_order(Omega,inc)

    loo=robust_leave_one_out(varpi)
    k=loo["remove_index"]
    mask=np.arange(len(rows))!=k
    mu_aps_rob,R_aps_rob=circ_mean_deg(varpi[mask])
    mu_node_rob,R_node_rob=node_order(Omega[mask],inc[mask])

    receipt={
        "schema":"PBH_RELATIONAL_ORBITAL_INVERSION_STAGE_A_V0_1",
        "status":"ANALYZED",
        "source":{
            "authority":"NASA_JPL_SBDB",
            "sample_origin":"Brown_Batygin_2021_Table_1",
            "selection":"150<a<1000 AU and q>42 AU; 11-object historical sample, current JPL solutions",
        },
        "objects":[{**r,"varpi":float((r["Omega"]+r["omega"])%360.0)} for r in rows],
        "relations":{
            "all_11":{
                "apsidal_order":{"angle_deg":mu_aps,"R":R_aps},
                "node_order":{"angle_deg":mu_node,"R":R_node},
                "anti_apsidal_perturber_phase_deg":float((mu_aps+180.0)%360.0),
            },
            "leave_one_out_robust":{
                "algorithm":"remove exactly one object whose removal maximizes apsidal resultant R",
                "removed_index":int(k),
                "removed_object":rows[k]["fullname"],
                "apsidal_order":{"angle_deg":mu_aps_rob,"R":R_aps_rob},
                "node_order":{"angle_deg":mu_node_rob,"R":R_node_rob},
                "anti_apsidal_perturber_phase_deg":float((mu_aps_rob+180.0)%360.0),
            },
        },
        "bootstrap":bootstrap_cluster(rows,k),
        "geometric_null":null_test(R_aps,R_node,inc),
        "relational_graph":pairwise_graph(rows),
        "inference_boundary":{
            "what_stage_A_can_identify":[
                "coherent apsidal phase",
                "coherent nodal/pole phase",
                "anti-apsidal orientation anchor for an external perturber",
            ],
            "what_stage_A_cannot_identify":[
                "mass",
                "semimajor_axis",
                "eccentricity",
                "current true anomaly or sky position",
                "whether the perturber is a PBH rather than a planet",
            ],
        },
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    print("PBH_RELATIONAL_STAGE_A="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
