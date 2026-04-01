import ionerdss as ion
from ionerdss import build_system_from_pdb

pdb_id = "7pg9"

print("--- Test 1: Triggering warning ---")
try:
    system = build_system_from_pdb(
        source=pdb_id,
        workspace_path=f"test_overlap_dir",
        interface_detect_distance_cutoff=0.6,
        interface_detect_n_residue_cutoff=3,
        chain_grouping_seq_threshold=0.5,
        nerdss_overlap_sep_limit=3.0,  # This should be too large
        ode_enabled=False
    )
except Exception as e:
    import traceback
    traceback.print_exc()

print("\n--- Test 2: Bypassing safety check ---")
try:
    system2 = build_system_from_pdb(
        source=pdb_id,
        workspace_path=f"test_overlap_dir2",
        interface_detect_distance_cutoff=0.6,
        interface_detect_n_residue_cutoff=3,
        chain_grouping_seq_threshold=0.5,
        nerdss_overlap_sep_limit=3.0,
        disable_overlap_sep_limit_check=True,
        ode_enabled=False
    )
except Exception as e:
    import traceback
    traceback.print_exc()
