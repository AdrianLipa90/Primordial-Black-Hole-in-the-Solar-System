#!/usr/bin/env python3
"""Candidate-only bridge from six descriptive CMB phase markers to GREMLIN T^36.

This uses the exact GREMLIN v0.1 circular coherence/similarity definitions, but
the 6->36 embedding is explicitly local to this PBH experiment and has no
canonical semantic-axis authority.
"""
from __future__ import annotations
import cmath
import json
import math
from pathlib import Path
from typing import Sequence

DIM=36
TAU=2.0*math.pi
FEATURE_ORDER=("phase_coherence","phase_entropy","A1","A2","A3","A6")


def circular_delta(a: float,b: float)->float:
    return (a-b+math.pi)%TAU-math.pi


def phase_coherence(a: Sequence[float],b: Sequence[float])->float:
    z=sum(cmath.exp(1j*circular_delta(x,y)) for x,y in zip(a,b))/DIM
    return abs(z)


def phase_similarity(a: Sequence[float],b: Sequence[float])->float:
    z=sum(cmath.exp(1j*circular_delta(x,y)) for x,y in zip(a,b))/DIM
    return 0.5*(1.0+z.real)


def embed_marker6(marker: dict[str,float])->tuple[float,...]:
    base=[TAU*(float(marker[k])%1.0) for k in FEATURE_ORDER]
    return tuple(base*6)


def circular_mean(states: list[tuple[float,...]])->tuple[float,...]:
    out=[]
    for i in range(DIM):
        z=sum(cmath.exp(1j*s[i]) for s in states)/len(states)
        out.append(cmath.phase(z)%TAU)
    return tuple(out)


def torus_distance(a: Sequence[float],b: Sequence[float])->float:
    return math.sqrt(sum(circular_delta(x,y)**2 for x,y in zip(a,b)))


def analyze(markers: dict[str,dict[str,float]])->dict:
    states={name:embed_marker6(m) for name,m in markers.items()}
    names=sorted(states)
    edges=[]
    for i,left in enumerate(names):
        for right in names[i+1:]:
            coh=phase_coherence(states[left],states[right])
            edges.append({
                "left":left,"right":right,
                "coherence":coh,
                "similarity":phase_similarity(states[left],states[right]),
                "spider_edge_at_0_80":coh>=0.80,
            })
    aggregate=circular_mean(list(states.values()))
    return {
        "schema":"PBH_GREMLIN_T36_MARKER_BRIDGE_V0_1",
        "status":"CANDIDATE_ONLY",
        "feature_order":list(FEATURE_ORDER),
        "embedding":"theta_i = 2pi*marker_i, six-marker block repeated six times",
        "semantic_axis_assignment":False,
        "canonical_embedding":False,
        "gremlin_contract":"GREMLIN_BESTIARY_PHASENAV_PHASE_GATES_V0_1",
        "pairwise":edges,
        "aggregate":{
            "mean_torus_disagreement":sum(torus_distance(aggregate,s) for s in states.values())/len(states),
            "coherence_to_aggregate":{n:phase_coherence(aggregate,s) for n,s in states.items()},
            "similarity_to_aggregate":{n:phase_similarity(aggregate,s) for n,s in states.items()},
        },
        "physical_detection_claim":False,
        "canon_allowed":False,
    }


if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("markers_json")
    ap.add_argument("--output",required=True)
    ns=ap.parse_args()
    payload=json.loads(Path(ns.markers_json).read_text())
    result=analyze(payload)
    Path(ns.output).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
