#!/usr/bin/env python3
"""
This script runs the ioNERDSS validation suite on a list of PDB IDs.
It builds the coarse-grained model, runs the NERDSS simulation, and computes the RMSD between the designed and observed structures.

Usage:
    python run_validation_benchmark.py --pdb_ids 1dlh 5l93 8y7s --nerdss_path /path/to/nerdss --output benchmark_results.csv

    python benchmark/run_validation_benchmark.py --pdb_list_file benchmark/pdb_ids_validation.csv --nerdss_path ~/Workspace/Reaction_ode/nerdss_development


| Status | Meaning | When it is assigned |
|---|---|---|
| `Success` | Validation succeeded | The target assembly was found and RMSD was computed successfully |
| `FP` | Too few protein chains | After coarse-graining, the structure has fewer than 2 protein chains, so assembly is guaranteed to fail |
| `DC` | Disconnected graph | The designed assembly graph is disconnected, so it cannot form a single target `N`-mer |
| `IC` | Interface at centre of mass | A molecule type's interface sites coincide or all lie within `--interface_com_proximity_threshold` of it, so NERDSS cannot define the binding angles (it would exit with "Cannot resolve phi angle"); the model is not simulated |
| `NC` | NERDSS crash | NERDSS itself crashed, segfaulted, aborted, or otherwise died during simulation/export |
| `UA` | Underassembly | The target assembly was not found, and the largest observed assembly size is smaller than the target assembly size |
| `OA` | Overassembly | The target assembly was not found, and the largest observed assembly size is greater than or equal to the target assembly size |
| `Failed_Assembly` | Old fallback failure label | Used only if the run failed assembly validation but does not cleanly fit `UA` or `OA` because the needed size metadata is unavailable |
| `Crashed` | Old fallback crash label | Used only if the run crashed but there is not enough evidence to specifically classify it as `NC`, `FP`, or `DC` |

"""
import argparse
import csv
import logging
import traceback
from pathlib import Path
from typing import Optional

from ionerdss.model.pdb import PDBModelBuilder
from ionerdss.model import pdb
from ionerdss.nerdss_simulation import resolve_nerdss_executable
from ionerdss.model.pdb.structure_validation import DEFAULT_INTERFACE_COM_PROXIMITY_THRESHOLD_NM
from ionerdss.model.pdb.structure_validation import get_box_fit_message
from ionerdss.model.pdb.structure_validation import get_degenerate_site_layouts
from ionerdss.model.pdb.structure_validation import get_disconnected_design_message
from ionerdss.model.pdb.structure_validation import get_free_interface_capacity
from ionerdss.model.pdb.structure_validation import get_interface_com_proximity_message
from ionerdss.model.pdb.structure_validation import get_near_com_interface_sites

FAST_VALIDATION_ITERATIONS = 100000


def _classify_failed_assembly(sim_result, target_assembly_size: int) -> str:
    """Return UA when the largest observed assembly is smaller than target, otherwise OA."""
    if sim_result.largest_observed_assembly_size < target_assembly_size:
        return "UA"
    return "OA"


def _contains_nerdss_crash_signature(message: Optional[str]) -> bool:
    """Return True when the text clearly indicates a NERDSS runtime crash."""
    if not message:
        return False

    normalized = message.lower()
    crash_signatures = (
        "segmentation fault",
        "segfault",
        "core dumped",
        "signal 11",
        "sigsegv",
        "abort trap",
        "stack trace",
    )
    return any(signature in normalized for signature in crash_signatures)


def _detect_nerdss_crash(sim_result) -> bool:
    """Return True when the simulation result includes clear evidence of a NERDSS crash."""
    if _contains_nerdss_crash_signature(getattr(sim_result, "warning_message", None)):
        return True

    simulation_dir = getattr(sim_result, "simulation_dir", None)
    if simulation_dir is None:
        return False

    log_path = Path(simulation_dir) / "output.log"
    if not log_path.exists():
        return False

    try:
        return _contains_nerdss_crash_signature(log_path.read_text(encoding="utf-8", errors="replace"))
    except OSError:
        return False


def _status_for_failed_validation(sim_result, target_assembly_size: int) -> str:
    """Prefer specific failure codes and fall back to the older ambiguous status when needed."""
    if _detect_nerdss_crash(sim_result):
        return "NC"

    largest_observed_assembly_size = getattr(sim_result, "largest_observed_assembly_size", None)
    if largest_observed_assembly_size is not None:
        return _classify_failed_assembly(sim_result, target_assembly_size)

    return "Failed_Assembly"


def _interface_com_proximity_threshold(builder) -> float:
    """Return the threshold the builder ran with, or the module default before it ran."""
    hyperparams = getattr(builder, "hyperparams", None)
    threshold = getattr(hyperparams, "interface_com_proximity_threshold", None)
    if threshold is None:
        return DEFAULT_INTERFACE_COM_PROXIMITY_THRESHOLD_NM
    return float(threshold)


def _site_geometry_status(system, threshold_nm: float) -> tuple[Optional[str], Optional[str]]:
    """Return ("IC", message) when NERDSS could not define the model's binding angles.

    Only a degenerate site layout -- a multi-interface molecule type whose sites
    coincide or all sit within the threshold of its centre of mass -- earns the
    status: those are the models NERDSS aborts on at the first
    association. A single site near the centre of mass is reported in the message but
    is still simulated, since NERDSS tolerates it.
    """
    message = get_interface_com_proximity_message(
        system,
        prefix="Validation preflight warning",
        threshold_nm=threshold_nm,
    )
    if message is None:
        return None, None
    if get_degenerate_site_layouts(system, threshold_nm=threshold_nm):
        return "IC", message
    return None, message


def _closest_site_com_distance(system) -> Optional[float]:
    """Return the smallest distance from a reacting interface site to its molecule's COM."""
    # A threshold of infinity lists every reacting site, sorted by distance.
    sites = get_near_com_interface_sites(system, threshold_nm=float("inf"))
    if not sites:
        return None
    return sites[0].distance_nm


def _partial_chain_counts(builder) -> tuple[int, int]:
    """Return best-effort chain and chain-type counts from partial builder state."""
    coarse_summary = getattr(builder, "coarse_summary", None) or {}
    group_summary = getattr(builder, "group_summary", None) or {}

    chains_count = int(coarse_summary.get("num_chains", 0) or 0)
    chain_types_count = int(group_summary.get("num_groups", 0) or 0)
    return chains_count, chain_types_count


def _status_from_partial_builder(builder) -> Optional[str]:
    """Infer FP/DC from partial pipeline state after a build failure."""
    chains_count, _chain_types_count = _partial_chain_counts(builder)

    # Only claim FP when coarse-graining actually ran and counted the chains. A builder
    # that died earlier -- a failed download, a full disk -- also reports zero chains, and
    # calling that "too few protein chains" turns an infrastructure failure into a
    # scientific result.
    coarse_summary = getattr(builder, "coarse_summary", None)
    coarse_graining_completed = bool(coarse_summary) and "num_chains" in coarse_summary

    if coarse_graining_completed and chains_count < 2:
        return "FP"
    if not coarse_graining_completed:
        return None

    system = getattr(builder, "system", None)
    if system is not None:
        disconnected_design_message = get_disconnected_design_message(
            system,
            prefix="Validation preflight warning",
        )
        if disconnected_design_message is not None:
            return "DC"

        site_geometry_status, _message = _site_geometry_status(
            system, _interface_com_proximity_threshold(builder)
        )
        if site_geometry_status is not None:
            return site_geometry_status

    return None


def _run_validation_attempt(
    system,
    workspace_manager,
    titration_rates,
    box_size: float,
    iterations: int,
    nerdss_path: str,
    sim_dir_name: str,
    copies: int = 1,
):
    # Supplying `copies` times the deposited stoichiometry is what makes over-assembly
    # reachable at all: with one copy the largest possible assembly IS the target, so
    # OA can never be observed. Scale the box by copies**(1/3) so the added subunits
    # raise the particle count without also raising the concentration.
    scaled_box = box_size * (float(copies) ** (1.0 / 3.0))
    artifacts = pdb.validation.setup_simulation(
        system,
        workspace_manager=workspace_manager,
        box_nm=(scaled_box, scaled_box, scaled_box),
        initial_molecule_count=copies,
        titration_on_rate=titration_rates,
        parms_overrides={
            "nItr": iterations,
            "timeWrite": 1000,
            "trajWrite": iterations,
            "restartWrite": iterations,
            "checkPoint": iterations,
            "pdbWrite": iterations,
            # NERDSS leaves bondedComplexWrite at -1 unless it is asked for, so DATA/COMPLEXES
            # stays empty and validation silently falls back to the terminal restart snapshot --
            # which only sees the final state and misses an assembly that formed and then grew.
            # 100 snapshots per run keeps the JSON output bounded at any nItr.
            "bondedComplexWrite": max(1, iterations // 100),
        },
    )

    print(f"  -> Running NERDSS validation simulation for {iterations} iterations...")
    sim_result = pdb.validation.run_simulation(
        artifacts,
        nerdss_path=nerdss_path,
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
    parser.add_argument("--nerdss_path", "--nerdss_dir", dest="nerdss_path", type=str, default=None,
                        help="The NERDSS executable, or a directory containing nerdss or nerdss_mpi directly "
                             "or in bin/ (e.g. a NERDSS checkout); defaults to nerdss on PATH. "
                             "--nerdss_dir is the old name")
    parser.add_argument("--output", default="benchmark_output/benchmark_results.csv", type=str, help="Output CSV file path")
    parser.add_argument("--iterations", default=1000000, type=int, help="Number of NERDSS iterations for the long rerun after the initial 100000-step probe")
    parser.add_argument("--box_size", default=50.0, type=float,
                        help="Box edge in nm for a single copy of the assembly; the actual box is "
                             "scaled by copies**(1/3) so concentration stays fixed")
    parser.add_argument("--geometric_regularization", default="off", type=str,
                        choices=["off", "auto"],
                        help="Snap a detected cyclic assembly onto exact Cn geometry before "
                             "simulating. Cyclic designs only close when the generator "
                             "transform is exactly n-fold; 'auto' enforces that")
    parser.add_argument("--copies", default=1, type=int,
                        help="Copies of the deposited stoichiometry to supply. >1 makes "
                             "over-assembly (OA) observable; the target composition is unchanged")
    parser.add_argument("--interface_com_proximity_threshold", default=None, type=float,
                        help="Distance in nm below which a reacting interface site counts as "
                             "sitting on its molecule's centre of mass (the IC preflight). "
                             f"Defaults to the hyperparameter value "
                             f"({DEFAULT_INTERFACE_COM_PROXIMITY_THRESHOLD_NM} nm)")
    parser.add_argument("--interface_site_placement", default="centroid", type=str,
                        choices=["centroid", "auto"],
                        help="Where interface sites are placed: 'centroid' of the contacting "
                             "residues, or 'auto' to move sites the IC preflight flags onto the "
                             "chain surface facing the partner so NERDSS can define the angles")
    
    args = parser.parse_args()
    # Resolve once up front so a wrong path fails here, not after the first model is built.
    args.nerdss_path = str(resolve_nerdss_executable(args.nerdss_path))
    print(f"Using NERDSS executable {args.nerdss_path}")

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
            writer.writerow(["PDB ID", "Number of chains in PDB", "Number of chain types in PDB", "Status", "RMSD", "Free interface slots", "Symmetry", "Closest site-COM distance (nm)"])
    
    for count, pdb_id in enumerate(all_pdb_ids, 1):
        print(f"\n[{count}/{len(all_pdb_ids)}] Testing PDB: {pdb_id}")
        
        chains_count = 0
        chain_types_count = 0
        status = "Crashed"
        rmsd = None
        free_interface_slots = ""
        symmetry = ""
        closest_site_com_distance = ""
        
        builder = None

        try:
            # 1. Full coarse graining
            builder_kwargs = {}
            if args.interface_com_proximity_threshold is not None:
                builder_kwargs["interface_com_proximity_threshold"] = args.interface_com_proximity_threshold
            builder = PDBModelBuilder(source=pdb_id)
            system = builder.build_system(
                workspace_path=f"benchmark/trials/{pdb_id}",
                generate_visualizations=False,
                generate_nerdss_files=False,
                logger_level=logging.WARNING,
                geometric_regularization=args.geometric_regularization,
                interface_site_placement=args.interface_site_placement,
                #interface_detect_distance_cutoff=1.5,
                #interface_detect_n_residue_cutoff=6,
                **builder_kwargs,
            )
            
            # Extract basic metrics
            chains_count = len(system.molecule_instances)
            rep_instances = pdb.validation.get_representative_instances(system)
            chain_types_count = len(rep_instances)
            
            print(f"  -> Extracted {chains_count} chains across {chain_types_count} unique types.")

            if chains_count < 2:
                status = "FP"
                print("  -> Too few protein chains remain after coarse-graining; marking as FP.")
            
            # Unused binding capacity predicts over-assembly; record it so the run can be
            # cross-tabulated against the outcome actually observed.
            detection = getattr(system, "symmetry_detection", None)
            symmetry = getattr(detection, "group", "") if detection is not None else ""

            free_interface_slots = sum(
                sum(interfaces.values())
                for interfaces in get_free_interface_capacity(system).values()
            )
            if free_interface_slots:
                print(f"  -> {free_interface_slots} free interface slot(s); over-assembly possible.")

            disconnected_design_message = get_disconnected_design_message(
                system,
                prefix="Validation preflight warning",
            )
            if disconnected_design_message is not None:
                status = "DC"
                print(f"  -> {disconnected_design_message}")

            # Interface sites on the centre of mass leave NERDSS unable to define the
            # binding angles; it aborts at the first association, which would otherwise
            # be scored as an empty run (under-assembly). Record how close the closest
            # site is so the threshold can be re-examined from the CSV.
            distance = _closest_site_com_distance(system)
            closest_site_com_distance = f"{distance:.4f}" if distance is not None else ""
            site_geometry_status, site_geometry_message = _site_geometry_status(
                system, _interface_com_proximity_threshold(builder)
            )
            if site_geometry_message is not None:
                print(f"  -> {site_geometry_message}")
            if site_geometry_status is not None and status not in {"FP", "DC"}:
                status = site_geometry_status
                print("  -> NERDSS cannot define the binding angles for this model; not simulating it.")

            box_edge = args.box_size * (float(args.copies) ** (1.0 / 3.0))
            box_fit_message = get_box_fit_message(
                system, (box_edge, box_edge, box_edge), prefix="Validation preflight warning"
            )
            if box_fit_message is not None:
                print(f"  -> {box_fit_message}")
            
            if status not in {"FP", "DC", "IC"}:
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
                    nerdss_path=args.nerdss_path,
                    sim_dir_name="validation_output_fast",
                    copies=args.copies,
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
                        nerdss_path=args.nerdss_path,
                        sim_dir_name="validation_output_full",
                        copies=args.copies,
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
                    target_assembly_size = sum(artifacts.target_counts.values())
                    status = _status_for_failed_validation(sim_result, target_assembly_size)
                    if sim_result.warning_message:
                        print(f"     Warning: {sim_result.warning_message}")
                        if "[Errno 2] No such file" in sim_result.warning_message or "invalid literal format" in sim_result.warning_message or "invalid literal for int" in sim_result.warning_message:
                            status = "Crashed"
                        elif status in {"UA", "OA"}:
                            print(
                                "     Error: target assembly was not found after both validation runs. "
                                f"Largest observed assembly size was {sim_result.largest_observed_assembly_size} "
                                f"vs target size {target_assembly_size}."
                            )
                    
                    _dump_nerdss_log(Path(sim_result.simulation_dir), pdb_id)
                
        except Exception as e:
            if builder is not None:
                partial_chains_count, partial_chain_types_count = _partial_chain_counts(builder)
                if partial_chains_count:
                    chains_count = partial_chains_count
                if partial_chain_types_count:
                    chain_types_count = partial_chain_types_count

                partial_status = _status_from_partial_builder(builder)
                if partial_status is not None:
                    status = partial_status
                else:
                    status = "NC" if _contains_nerdss_crash_signature(str(e)) else "Crashed"
            else:
                status = "NC" if _contains_nerdss_crash_signature(str(e)) else "Crashed"
            print(f"  -> Error encountered benchmarking {pdb_id}: {e}")
            if "DATA/restart.dat or any RESTART snapshot" in str(e):
                print("     Error: target composition appeared in the histogram, but no restart snapshot contained the full assembly.")
            traceback.print_exc()
            
        # Write to CSV
        with open(output_path, 'a', newline='') as f:
            writer = csv.writer(f)
            rmsd_val = f"{rmsd:.4f} nm" if rmsd is not None else ""
            writer.writerow([pdb_id, chains_count, chain_types_count, status, rmsd_val,
                             free_interface_slots, symmetry, closest_site_com_distance])
            
if __name__ == "__main__":
    main()
