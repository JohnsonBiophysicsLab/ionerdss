import csv
import json
import math
import time
from pathlib import Path

import requests

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
OUTDIR = Path("pdb_benchmark")
OUTDIR.mkdir(exist_ok=True)

def post_search(payload):
    r = requests.post(SEARCH_URL, json=payload, timeout=120)
    r.raise_for_status()
    return r.json()

def make_query(n_chains, homo=True, max_resolution=3.5):
    entity_count_rule = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_assembly_info.polymer_entity_count",
            "operator": "equals" if homo else "greater",
            "value": 1 if homo else 1
        }
    }

    # for heteromers, "greater than 1"
    if not homo:
        entity_count_rule["parameters"]["operator"] = "greater"
        entity_count_rule["parameters"]["value"] = 1

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.selected_polymer_entity_types",
                        "operator": "exact_match",
                        "value": "Protein (only)"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "in",
                        "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_assembly_info.polymer_entity_instance_count_protein",
                        "operator": "equals",
                        "value": n_chains
                    }
                },
                entity_count_rule,
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": max_resolution
                    }
                }
            ]
        },
        "return_type": "assembly",
        "request_options": {
            "return_all_hits": True,
            "results_verbosity": "minimal"
        }
    }

def fetch_ids(n_chains, homo):
    payload = make_query(n_chains=n_chains, homo=homo)
    data = post_search(payload)
    ids = [x["identifier"] for x in data.get("result_set", [])]
    return ids

def classify_label(n_chains, homo):
    prefix = "homo" if homo else "hetero"
    names = {
        2: "dimer",
        3: "trimer",
        4: "tetramer",
        5: "pentamer",
        6: "hexamer",
        7: "heptamer",
        8: "octamer",
    }
    return prefix + names.get(n_chains, f"_{n_chains}mer")

def build_candidate_table(chain_counts=(2, 3, 4, 5, 6, 8)):
    rows = []
    for n in chain_counts:
        for homo in (True, False):
            ids = fetch_ids(n, homo)
            label = classify_label(n, homo)
            for assembly_id in ids:
                pdb_id, asm = assembly_id.split("-")
                rows.append({
                    "pdb_id": pdb_id.lower(),
                    "assembly_id": asm,
                    "class_label": label,
                    "n_protein_chains": n,
                    "homo": homo,
                })
            time.sleep(0.2)
    return rows

def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

if __name__ == "__main__":
    rows = build_candidate_table(chain_counts=(2, 3, 4, 5, 6, 8))
    write_csv(rows, OUTDIR / "assembly_candidates.csv")
    print(f"Wrote {len(rows)} candidates to {OUTDIR / 'assembly_candidates.csv'}")