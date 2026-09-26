#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://datasets-server.huggingface.co/rows"


def fetch_rows(dataset: str, config: str, split: str, offset: int, length: int) -> list[dict]:
    out = []
    cursor = offset
    remaining = length
    while remaining:
        take = min(100, remaining)
        query = urllib.parse.urlencode({
            "dataset": dataset,
            "config": config,
            "split": split,
            "offset": cursor,
            "length": take,
        })
        request = urllib.request.Request(
            API + "?" + query,
            headers={"User-Agent": "PBH-CMB-row-proxy/0.1 read-only"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = json.load(response)
        rows = payload.get("rows") or []
        if len(rows) != take:
            raise RuntimeError(f"short row slice at offset {cursor}: expected {take}, got {len(rows)}")
        expected = list(range(cursor, cursor + take))
        actual = [int(item["row_idx"]) for item in rows]
        if actual != expected:
            raise RuntimeError(f"row order mismatch at offset {cursor}")
        out.extend(rows)
        cursor += take
        remaining -= take
        time.sleep(0.05)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--split", default="train")
    ap.add_argument("--offset", type=int, required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--field", required=True)
    ap.add_argument("--output", required=True)
    ns = ap.parse_args()

    rows = fetch_rows(ns.dataset, ns.config, ns.split, ns.offset, ns.length)
    values = []
    for item in rows:
        row = item.get("row") or {}
        if ns.field not in row:
            raise RuntimeError(f"field {ns.field!r} missing at row {item.get('row_idx')}")
        values.append(row[ns.field])

    receipt = {
        "schema": "PBH_HF_ROW_SLICE_RECEIPT_V0_1",
        "status": "FETCHED",
        "dataset": ns.dataset,
        "config": ns.config,
        "split": ns.split,
        "offset": ns.offset,
        "length": ns.length,
        "field": ns.field,
        "row_indices": [int(item["row_idx"]) for item in rows],
        "values": values,
        "source_api": API,
        "authority": "SECONDARY_SOURCE_FAITHFUL_PROXY",
        "canon_allowed": False,
    }
    out = Path(ns.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
