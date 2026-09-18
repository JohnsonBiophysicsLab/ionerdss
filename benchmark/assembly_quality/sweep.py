#!/usr/bin/env python3
"""Compare hyperparameter settings by the assembly their NERDSS runs actually produce.

Each named setting in the config file is built once with the PDB pipeline, then simulated
once per seed. Every run is scored with `assembly_metrics`: whether any complex is
mis-assembled, whether NERDSS left a stretched bond, and how completely the design's own
bonds formed. One CSV row per run.

    python benchmark/assembly_quality/sweep.py \
        --source 6bno --configs benchmark/assembly_quality/configs_6bno.json \
        --nerdss_dir ~/Workspace/NERDSS --workspace_root /tmp/6bno_sweep \
        --seeds 1 2 3 --iterations 100000 --output results_6bno.csv

The config file maps a name to hyperparameter overrides, everything else staying at the
pipeline defaults. Names become directory names under --workspace_root, so keep them to
characters a path accepts:

    {"default": {}, "cutoff_1.0 + limit_3.0": {"interface_detect_distance_cutoff": 1.0,
                                               "nerdss_overlap_sep_limit": 3.0}}
"""
import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.main import PDBModelBuilder

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assembly_metrics import Model, analyse_run, closure_report  # noqa: E402

FIELDS = [
    "config", "hyperparameters", "seed", "iterations",
    "n_sites", "sites", "n_reactions", "n_self_binding_reactions",
    "mis_assembled_snapshots", "snapshots", "first_mis_assembled_step", "largest_mis_assembled",
    "worst_min_com_separation_nm", "self_binding_bonds_final",
    "stretched_bond_snapshots", "first_stretched_bond_step",
    "largest_complex_final", "complex_sizes_final", "lattice_bond_completeness_final",
    "closure_translation_error_nm", "closure_rotation_error_deg", "closure_closes",
    "wall_s", "nerdss_returncode", "error",
]


def build_model(source, name, overrides, workspace_root):
    """Build one setting's model and return (workspace, topology fields)."""
    workspace = Path(workspace_root) / "builds" / name
    hyperparams = PDBModelHyperparameters(**{"generate_visualizations": False, **overrides})
    PDBModelBuilder(source, hyperparams=hyperparams).build_system(workspace_path=str(workspace))

    model = Model.from_nerdss_dir(workspace / "nerdss_files")
    sites = sorted(f"{mol}({site})" if len(model.sites) > 1 else site
                   for mol, mol_sites in model.sites.items() for site in mol_sites)
    self_binding = sum(1 for (first, second) in model.radii if first == second)
    return workspace, {
        "n_sites": len(sites),
        "sites": " ".join(sites),
        "n_reactions": len(model.radii),
        "n_self_binding_reactions": self_binding,
    }


def simulate(nerdss_bin, workspace, sim_dir, seed, iterations, snapshots, timeout):
    """Run NERDSS on a built model and return (returncode, wall seconds)."""
    sim_dir = Path(sim_dir)
    if sim_dir.exists():
        shutil.rmtree(sim_dir)
    shutil.copytree(Path(workspace) / "nerdss_files", sim_dir)

    parms = (sim_dir / "parms.inp").read_text()
    parms = re.sub(r"nItr = [\d.eE+]+", f"nItr = {iterations}", parms)
    # Trajectories and restarts are only needed at the end; the snapshots carry the analysis.
    for key in ("trajWrite", "restartWrite", "checkPoint", "pdbWrite"):
        parms = re.sub(rf"{key} = [\d.eE+]+", f"{key} = {iterations}", parms)
    interval = max(1, iterations // snapshots)
    parms = re.sub(r"timeWrite = [\d.eE+]+", f"timeWrite = {interval}", parms)
    parms = parms.replace("end parameters", f"    bondedComplexWrite = {interval}\nend parameters", 1)
    (sim_dir / "parms.inp").write_text(parms)

    started = time.time()
    with open(sim_dir / "output.log", "w") as log:
        completed = subprocess.run(
            [str(nerdss_bin), "-f", "parms.inp", "-s", str(seed)],
            cwd=sim_dir, stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
        )
    return completed.returncode, round(time.time() - started, 1)


def run_one(args, name, overrides, workspace, topology, seed):
    row = {"config": name, "hyperparameters": json.dumps(overrides, sort_keys=True),
           "seed": seed, "iterations": args.iterations, **topology}
    sim_dir = Path(args.workspace_root) / "sims" / f"{name}_n{args.iterations}_s{seed}"
    try:
        returncode, wall = simulate(
            Path(args.nerdss_dir) / "bin" / "nerdss", workspace, sim_dir,
            seed, args.iterations, args.snapshots, args.timeout,
        )
        system_json = next((Path(workspace) / "outputs" / "systems").glob("*_system.json"))
        row.update(analyse_run(sim_dir, system_json, clash_nm=args.clash_nm))
        row["nerdss_returncode"] = returncode
        row["wall_s"] = wall
        closure = closure_report(sim_dir, system_json)
        if closure:
            row["closure_translation_error_nm"] = closure.get("translation_error_nm")
            row["closure_rotation_error_deg"] = closure.get("rotation_error_deg")
            row["closure_closes"] = closure.get("closes")
    except Exception as exc:  # one failed run must not stop the sweep
        row["error"] = f"{type(exc).__name__}: {exc}"[:300]
    row["complex_sizes_final"] = " ".join(str(size) for size in row.get("complex_sizes_final", [])[:10])
    return {field: row.get(field, "") for field in FIELDS}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True, help="PDB ID or path to a structure file")
    parser.add_argument("--configs", required=True, help="JSON file mapping a name to hyperparameter overrides")
    parser.add_argument("--nerdss_dir", required=True, help="NERDSS repository directory (expects bin/nerdss)")
    parser.add_argument("--workspace_root", required=True, help="Directory for builds and simulations")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3], help="NERDSS random seeds, one run each")
    parser.add_argument("--iterations", type=int, default=100000, help="NERDSS iterations per run")
    parser.add_argument("--snapshots", type=int, default=50, help="Bonded-complex snapshots per run")
    parser.add_argument("--clash_nm", type=float, default=3.0,
                        help="Subunit centres closer than this count as steric overlap")
    parser.add_argument("--timeout", type=float, default=3600.0, help="Per-run wall-clock limit in seconds")
    parser.add_argument("--jobs", type=int, default=4, help="Simulations to run at once")
    args = parser.parse_args()

    configs = json.loads(Path(args.configs).read_text())
    jobs, failures = [], []
    for name, overrides in configs.items():
        try:
            workspace, topology = build_model(args.source, name, overrides, args.workspace_root)
        except Exception as exc:
            print(f"{name}: build failed: {type(exc).__name__}: {exc}")
            failures.append({field: "" for field in FIELDS} | {
                "config": name, "hyperparameters": json.dumps(overrides, sort_keys=True),
                "error": f"build failed: {type(exc).__name__}: {exc}"[:300],
            })
            continue
        print(f"{name}: {topology['n_sites']} sites ({topology['sites']}), "
              f"{topology['n_reactions']} reactions, {topology['n_self_binding_reactions']} self-binding")
        jobs.extend((name, overrides, workspace, topology, seed) for seed in args.seeds)

    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(lambda job: run_one(args, *job), jobs))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(failures + rows)

    print(f"\nwrote {len(rows)} runs to {output}")
    for name in configs:
        runs = [row for row in rows if row["config"] == name]
        if not runs:
            continue
        bad = sum(1 for row in runs if row["mis_assembled_snapshots"])
        stretched = sum(1 for row in runs if row["stretched_bond_snapshots"])
        print(f"  {name}: {bad}/{len(runs)} runs mis-assembled, {stretched}/{len(runs)} with a stretched bond")


if __name__ == "__main__":
    main()
