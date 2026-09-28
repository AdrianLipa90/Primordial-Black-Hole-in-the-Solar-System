#!/usr/bin/env python3
from __future__ import annotations
import json, math, urllib.parse, urllib.request

M_EARTH_PER_M_SUN = 332946.0487

# Nominal giant-planet mass ratios Mplanet/Msun and semimajor axes [AU].
GIANTS = {
    "Jupiter": {"mu": 0.0009543, "a_au": 5.2044},
    "Saturn": {"mu": 0.0002857, "a_au": 9.5826},
    "Uranus": {"mu": 0.00004365, "a_au": 19.2184},
    "Neptune": {"mu": 0.00005149, "a_au": 30.1104},
}

OBJECTS=[
"2000 CR105","Sedna","2004 VN112","2010 GB174","2012 VP113",
"2013 FT28","2013 RA109","2013 SY99","2013 UT15","2014 SR349","2015 RX245"
]

def fetch(name):
    url="https://ssd-api.jpl.nasa.gov/sbdb.api?"+urllib.parse.urlencode(
        {"sstr":name,"full-prec":"1","cd-epoch":"1"}
    )
    req=urllib.request.Request(url,headers={"User-Agent":"GREMLIN-PBH-secular-strength/0.1 read-only"})
    with urllib.request.urlopen(req,timeout=30) as resp:
        d=json.loads(resp.read())
    o=d["orbit"]; el={e["name"]:e for e in o["elements"]}
    return {
        "name":d["object"].get("fullname"),
        "epoch":float(o["epoch"]),
        "a":float(el["a"]["value"]),
        "e":float(el["e"]["value"]),
        "q":float(el["q"]["value"]),
        "i":float(el["i"]["value"]),
        "Omega":float(el["om"]["value"]),
        "omega":float(el["w"]["value"]),
        "varpi":(float(el["om"]["value"])+float(el["w"]["value"]))%360.0,
    }

def kcand(m_earth,a_p,e_p):
    mu=m_earth/M_EARTH_PER_M_SUN
    return mu/(a_p**3*(1.0-e_p**2)**1.5)

def mass_for_K(K,a_p,e_p):
    return K*M_EARTH_PER_M_SUN*a_p**3*(1.0-e_p**2)**1.5

def main():
    rows=[fetch(x) for x in OBJECTS]
    epochs=sorted({r["epoch"] for r in rows})
    if len(epochs)!=1:
        raise RuntimeError(f"common epoch failed: {epochs}")

    S=sum(v["mu"]*v["a_au"]**2 for v in GIANTS.values())
    for r in rows:
        r["K_req_au^-3"]=S/(r["a"]**5*(1.0-r["e"]**2)**2)

    high=[r for r in rows if r["a"]>250.0]
    low=[r for r in rows if r["a"]<=250.0]

    K_low=max(r["K_req_au^-3"] for r in high)
    K_high=min(r["K_req_au^-3"] for r in low)
    lower_driver=max(high,key=lambda r:r["K_req_au^-3"])
    upper_driver=min(low,key=lambda r:r["K_req_au^-3"])

    if not K_low < K_high:
        raise RuntimeError("no secular-strength bracket from selected cohorts")

    grid=[]
    for a_p in (250,300,350,380,430,500,600,700):
        for e_p in (0.0,0.2,0.4,0.6,0.78):
            q_p=a_p*(1.0-e_p)
            grid.append({
                "a_p_au":a_p,"e_p":e_p,"q_p_au":q_p,
                "mass_min_Mearth":mass_for_K(K_low,a_p,e_p),
                "mass_max_Mearth":mass_for_K(K_high,a_p,e_p),
                "mass_mid_geom_Mearth":math.sqrt(
                    mass_for_K(K_low,a_p,e_p)*mass_for_K(K_high,a_p,e_p)
                ),
            })

    old_mass_earth=5e-5*M_EARTH_PER_M_SUN
    old={"mass_Msun":5e-5,"mass_Mearth":old_mass_earth,"a_p_au":430.0,"e_p":0.78}
    old["q_p_au"]=old["a_p_au"]*(1-old["e_p"])
    old["K_au^-3"]=kcand(old["mass_Mearth"],old["a_p_au"],old["e_p"])
    old["K_over_bracket_low"]=old["K_au^-3"]/K_low
    old["K_over_bracket_high"]=old["K_au^-3"]/K_high

    reference={
        "label":"Brown_Batygin_2021_central_values_comparison_only",
        "mass_Mearth":6.2,"a_p_au":380.0,"q_p_au":300.0
    }
    reference["e_p"]=1.0-reference["q_p_au"]/reference["a_p_au"]
    reference["K_au^-3"]=kcand(reference["mass_Mearth"],reference["a_p_au"],reference["e_p"])

    receipt={
        "schema":"PBH_RELATIONAL_SECULAR_STRENGTH_STAGE_B1_V0_1",
        "status":"ANALYZED",
        "epoch_jd":epochs[0],
        "giant_quadrupole_moment_proxy_S_au2":S,
        "definition":{
            "K":"(m_p/Msun)/(a_p^3*(1-e_p^2)^(3/2))",
            "K_req_i":"S/(a_i^5*(1-e_i^2)^2)",
            "interpretation":"leading-order coplanar quadrupole precession-balance proxy; common prefactors cancel; inclination/resonant/octupole effects omitted",
        },
        "cohorts":{
            "outer":{"criterion":"a>250 AU","n":len(high)},
            "inner_boundary":{"criterion":"a<=250 AU within the 11-object detached sample","n":len(low)},
        },
        "objects":rows,
        "relational_strength_bracket":{
            "K_min_au^-3":K_low,
            "K_max_au^-3":K_high,
            "lower_boundary_object":lower_driver["name"],
            "upper_boundary_object":upper_driver["name"],
            "ratio_width":K_high/K_low,
        },
        "mass_orbit_ridge":grid,
        "old_2025_candidate_check":old,
        "external_reference_not_used_for_fit":reference,
        "inference_boundary":{
            "supported":"a one-parameter secular torque-strength ridge K, conditional on the 250 AU transition model",
            "not_supported":[
                "unique mass","unique semimajor axis","unique eccentricity",
                "current anomaly or sky location","PBH identity"
            ]
        },
        "authority":"CANDIDATE_ONLY",
        "canon_allowed":False,
    }
    print("PBH_STAGE_B1="+json.dumps(receipt,sort_keys=True,separators=(",",":")),flush=True)
    return 0

if __name__=="__main__": raise SystemExit(main())
