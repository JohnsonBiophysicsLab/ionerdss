#!/usr/bin/env python3
"""
This script runs the ioNERDSS validation suite on a list of PDB IDs.
It builds the coarse-grained model, runs the NERDSS simulation, and computes the RMSD between the designed and observed structures.

Usage:
    python run_validation_benchmark.py --pdb_ids 1dlh 5l93 8y7s --nerdss_dir /path/to/nerdss --output benchmark_results.csv

    python benchmark/run_validation_benchmark.py --pdb_list_file benchmark/pdb_ids_validation.csv --nerdss_dir ~/Workspace/Reaction_ode/nerdss_development
"""
import argparse
import csv
import logging
import traceback
from pathlib import Path

from ionerdss.model.pdb import PDBModelBuilder
from ionerdss.model import pdb

FAST_VALIDATION_ITERATIONS = 100000


def _run_validation_attempt(
    system,
    workspace_manager,
    titration_rates,
    box_size: float,
    iterations: int,
    nerdss_dir: str,
    sim_dir_name: str,
):
    artifacts = pdb.validation.setup_simulation(
        system,
        workspace_manager=workspace_manager,
        box_nm=(box_size, box_size, box_size),
        initial_molecule_count=1,
        titration_on_rate=titration_rates,
        parms_overrides={
            "nItr": iterations,
            "timeWrite": 1000,
            "trajWrite": iterations,
            "restartWrite": iterations,
            "checkPoint": iterations,
            "pdbWrite": iterations,
        },
    )

    print(f"  -> Running NERDSS validation simulation for {iterations} iterations...")
    sim_result = pdb.validation.run_simulation(
        artifacts,
        nerdss_dir=nerdss_dir,
        sim_dir_name=sim_dir_name,
    )
    return artifacts, sim_result


def _dump_nerdss_log(simulation_dir: Path, pdb_id: str):
    try:
        log_src = simulation_dir / "output.log"
        if log_src.exists():
            import shutil

            log_dst_dir = Path("benchmark_output/logs")
            log_dst_dir.mkdir(parents=True, exist_ok=True)
            log_dst = log_dst_dir / f"{pdb_id}_nerdss_crash.log"
            shutil.copy(log_src, log_dst)
            print(f"     Dumping NERDSS stdout to {log_dst}")
    except Exception as ex:
        print(f"     Could not dump NERDSS log for {pdb_id}: {ex}")


def main():
    parser = argparse.ArgumentParser(description="Run ioNERDSS validation suite on a list of PDB IDs.")
    parser.add_argument("--pdb_ids", nargs="+", help="List of PDB IDs to benchmark")
    parser.add_argument("--pdb_list_file", type=str, help="Path to a text or CSV file containing PDB IDs separated by commas or newlines")
    parser.add_argument("--nerdss_dir", required=True, type=str, help="Path to the compiled NERDSS binary directory")
    parser.add_argument("--output", default="benchmark_output/benchmark_results.csv", type=str, help="Output CSV file path")
    parser.add_argument("--iterations", default=1000000, type=int, help="Number of NERDSS iterations for the long rerun after the initial 100000-step probe")
    parser.add_argument("--box_size", default=50.0, type=float, help="Box dimensions in nm")
    
    args = parser.parse_args()
    
    import re
    all_pdb_ids = []
    if args.pdb_ids:
        all_pdb_ids.extend(args.pdb_ids)
    if args.pdb_list_file:
        list_path = Path(args.pdb_list_file)
        if list_path.exists():
            content = list_path.read_text()
            content = re.sub(r'#.*', '', content) # Ignore comments
            file_pdbs = [p.strip() for p in re.split(r'[,\n]', content) if p.strip()]
            all_pdb_ids.extend(file_pdbs)
        else:
            parser.error(f"File not found: {args.pdb_list_file}")

    if not all_pdb_ids:
        parser.error("You must provide either --pdb_ids or --pdb_list_file containing at least one PDB ID.")
        
    # Remove duplicates but preserve order
    seen = set()
    all_pdb_ids = [x for x in all_pdb_ids if not (x in seen or seen.add(x))]
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Initialize CSV if it doesn't exist
    if not output_path.exists():
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["PDB ID", "Number of chains in PDB", "Number of chain types in PDB", "Status", "RMSD"])
    
    for count, pdb_id in enumerate(all_pdb_ids, 1):
        print(f"\n[{count}/{len(all_pdb_ids)}] Testing PDB: {pdb_id}")
        
        chains_count = 0
        chain_types_count = 0
        status = "Crashed"
        rmsd = None
        
        try:
            # 1. Full coarse graining
            builder = PDBModelBuilder(source=pdb_id)
            system = builder.build_system(
                workspace_path=f"benchmark/trials/{pdb_id}",
                generate_visualizations=False,
                generate_nerdss_files=False,
                logger_level=logging.WARNING,
                #interface_detect_distance_cutoff=1.5,
                #interface_detect_n_residue_cutoff=6,
            )
            
            # Extract basic metrics
            chains_count = len(system.molecule_instances)
            rep_instances = pdb.validation.get_representative_instances(system)
            chain_types_count = len(rep_instances)
            
            print(f"  -> Extracted {chains_count} chains across {chain_types_count} unique types.")
            
            # Create a range of titration rates, scaled for each molecular species
            base_rate = 0.0 # 0.25e-3
            titration_rates = {}
            for idx, mol_name in enumerate(rep_instances.keys()):
                # Linearly increment the rate by 50% for each unique species
                titration_rates[mol_name] = base_rate * (1.0 + (idx * 0.5))
                
            artifacts = None
            sim_result = None

            # 2. Run a short validation simulation first.
            artifacts, sim_result = _run_validation_attempt(
                system=system,
                workspace_manager=builder.workspace_manager,
                titration_rates=titration_rates,
                box_size=args.box_size,
                iterations=FAST_VALIDATION_ITERATIONS,
                nerdss_dir=args.nerdss_dir,
                sim_dir_name="validation_output_fast",
            )

            # 3. If the target never appears in the histogram, rerun longer.
            if not sim_result.full_assembly_found:
                print(
                    "  -> Target composition did not appear within "
                    f"{FAST_VALIDATION_ITERATIONS} iterations; rerunning with {args.iterations} iterations."
                )
                artifacts, sim_result = _run_validation_attempt(
                    system=system,
                    workspace_manager=builder.workspace_manager,
                    titration_rates=titration_rates,
                    box_size=args.box_size,
                    iterations=args.iterations,
                    nerdss_dir=args.nerdss_dir,
                    sim_dir_name="validation_output_full",
                )
            
            # 4. Check results and compute RMSD if successful
            if sim_result.full_assembly_found and sim_result.observed_coordinates:
                print("  -> Full assembly found! Aligning structure...")
                alignment = pdb.validation.align_structure(
                    artifacts.designed_coordinates,
                    sim_result.observed_coordinates,
                    backend='kabsch',
                )
                rmsd = alignment.rmsd
                status = "Success"
                print(f"  -> Validation Complete! RMSD: {rmsd:.4f} nm")
            else:
                print("  -> Simulation ran but did not yield a full matching assembly.")
                if sim_result.warning_message:
                    print(f"     Warning: {sim_result.warning_message}")
                    if "[Errno 2] No such file" in sim_result.warning_message or "invalid literal format" in sim_result.warning_message or "invalid literal for int" in sim_result.warning_message:
                        status = "Crashed"
                    else:
                        status = "Failed_Assembly"
                        print("     Error: target composition did not appear in the histogram after both validation runs.")
                else:
                    status = "Failed_Assembly"
                
                _dump_nerdss_log(Path(sim_result.simulation_dir), pdb_id)
                
        except Exception as e:
            status = "Crashed"
            print(f"  -> Error encountered benchmarking {pdb_id}: {e}")
            if "DATA/restart.dat or any RESTART snapshot" in str(e):
                print("     Error: target composition appeared in the histogram, but no restart snapshot contained the full assembly.")
            traceback.print_exc()
            
        # Write to CSV
        with open(output_path, 'a', newline='') as f:
            writer = csv.writer(f)
            rmsd_val = f"{rmsd:.4f} nm" if rmsd is not None else ""
            writer.writerow([pdb_id, chains_count, chain_types_count, status, rmsd_val])
            
if __name__ == "__main__":
    main()
