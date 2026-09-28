#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from datetime import datetime

import requests

CATALOG = "J/AJ/152/80/obs"
URL = (
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
    "?-source=J%2FAJ%2F152%2F80%2Fobs"
    "&-out.all"
    "&-out.max=unlimited"
)


def parse_date(text: str) -> float:
    # YYYY:MM:DD.dddddd -> decimal year, sufficient for catalog span summaries.
    y, m, d = text.strip().split(":")
    year = int(y)
    month = int(m)
    day = float(d)
    return year + ((month - 1) + (day - 1) / 31.0) / 12.0


def main() -> int:
    r = requests.get(
        URL,
        timeout=180,
        headers={"User-Agent": "GREMLIN-FERRET-PBH-ORBITAL/0.1 read-only"},
    )
    r.raise_for_status()
    payload = r.content
    text = payload.decode("utf-8", "replace")
    lines = text.splitlines()

    data_lines = [line for line in lines if line and not line.startswith("#")]
    if not data_lines:
        raise RuntimeError("VizieR response contains no table rows")

    # ASU-TSV normally includes header, units, and separator lines after comments.
    # Detect the header by required catalog labels.
    header_idx = None
    for i, line in enumerate(data_lines):
        if "Planet" in line and "Date" in line and "RAJ2000" in line and "DEJ2000" in line:
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError("could not find VizieR table header")

    header = data_lines[header_idx].split("\t")
    rows = []
    for line in data_lines[header_idx + 1:]:
        parts = line.split("\t")
        if len(parts) != len(header):
            continue
        row = dict(zip(header, parts))
        if row.get("Planet") == "Planet" or row.get("Date") in ("Date", "---", ""):
            continue
        if not row.get("Date", "").strip()[:4].isdigit():
            continue
        rows.append(row)

    if len(rows) < 40000:
        raise RuntimeError(f"unexpectedly few parsed rows: {len(rows)}")

    by_object = Counter()
    by_observatory = Counter()
    outliers = Counter()
    date_span = {}
    uncertainty = defaultdict(list)

    for row in rows:
        key = (row.get("Planet", "").strip() + " " + row.get("Prov", "").strip()).strip()
        by_object[key] += 1
        by_observatory[row.get("Obs", "").strip()] += 1
        outliers[row.get("Out", "").strip()] += 1
        yr = parse_date(row["Date"])
        if key not in date_span:
            date_span[key] = [yr, yr]
        else:
            date_span[key][0] = min(date_span[key][0], yr)
            date_span[key][1] = max(date_span[key][1], yr)
        for field in ("e_RAJ2000", "e_DEJ2000"):
            raw = row.get(field, "").strip()
            try:
                uncertainty[field].append(float(raw))
            except Exception:
                pass

    pluto_keys = [k for k in by_object if k.strip() in ("134340", "Pluto") or "Pluto" in k]
    # MPC number 134340 is expected for Pluto, but preserve observed labels.
    top_objects = by_object.most_common(30)

    receipt = {
        "schema": "PBH_HOLMAN_PAYNE_CDS_INGEST_V0_1",
        "status": "INGESTED",
        "catalog": CATALOG,
        "source_url": URL,
        "http_status": r.status_code,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "parsed_rows": len(rows),
        "unique_objects": len(by_object),
        "outlier_flag_counts": dict(outliers),
        "pluto_candidate_labels": pluto_keys,
        "top_objects": top_objects,
        "top_observatories": by_observatory.most_common(30),
        "selected_object_spans": {
            k: date_span[k] for k, _ in top_objects[:30]
        },
        "uncertainty_arcsec_summary": {
            field: {
                "n": len(vals),
                "min": min(vals) if vals else None,
                "max": max(vals) if vals else None,
                "mean": (sum(vals) / len(vals)) if vals else None,
            }
            for field, vals in uncertainty.items()
        },
        "first_data_rows": rows[:5],
        "authority": "OBSERVATIONAL_DATA_ONLY",
        "canon_allowed": False,
    }

    compact = {
        k: receipt[k]
        for k in (
            "schema", "status", "bytes", "sha256", "parsed_rows",
            "unique_objects", "outlier_flag_counts", "pluto_candidate_labels",
            "top_objects", "top_observatories",
            "uncertainty_arcsec_summary",
        )
    }
    print(
        "PBH_ORBITAL_CATALOG_SUMMARY="
        + json.dumps(compact, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    print(
        "PBH_ORBITAL_CATALOG_RECEIPT="
        + json.dumps(receipt, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
