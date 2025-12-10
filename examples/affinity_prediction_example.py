"""Example demonstrating affinity prediction in coarse-graining.

This example shows how to use the predict_affinity option in the coarse_grain method
to predict binding energies between interacting chains using ProAffinity-GNN.

The implementation uses the already-downloaded PDB file, avoiding redundant downloads.
"""

import sys
sys.path.append("..")

from ionerdss.model.pdb_model import PDBModel

# Example 1: Using default fixed binding energy (-16 RT)
print("=" * 80)
print("Example 1: Coarse-graining without affinity prediction (default)")
print("=" * 80)

model1 = PDBModel(pdb_id="8erq", save_dir="./8erq_default")
model1.coarse_grain(
    distance_cutoff=0.35,
    residue_cutoff=3,
    predict_affinity=False,  # Use default fixed energy
    standard_output=True
)

print("\n" + "=" * 80)
print("Example 2: Coarse-graining with ProAffinity prediction")
print("=" * 80)

# Example 2: Using ProAffinity to predict binding energies
model2 = PDBModel(pdb_id="8erq", save_dir="./8erq_predicted")
model2.coarse_grain(
    distance_cutoff=0.35,
    residue_cutoff=3,
    predict_affinity=True,  # Enable affinity prediction
    adfr_path='/home/local/WIN/msang2/ADFRsuite-1.0/bin/prepare_receptor',  # Path to ADFR
    standard_output=True
)

print("\n" + "=" * 80)
print("Comparison Summary")
print("=" * 80)
print(f"Model 1 (fixed energy): {len(model1.all_interfaces)} chains with interfaces")
print(f"Model 2 (predicted energy): {len(model2.all_interfaces)} chains with interfaces")
print("\nInterface energies comparison:")
for i, chain in enumerate(model1.all_chains):
    if model1.all_interface_energies[i]:
        print(f"\nChain {chain.id}:")
        print(f"  Fixed energies: {[f'{e:.2f}' for e in model1.all_interface_energies[i]]}")
        print(f"  Predicted energies: {[f'{e:.2f}' for e in model2.all_interface_energies[i]]}")
