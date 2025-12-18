"""
Example: Using ODE Auto-Pipeline with ionerdss

This script demonstrates how to automatically calculate ODE solutions
for molecular assembly before running NERDSS simulations.

The ODE pipeline predicts concentration time courses based on reaction
kinetics, which can be compared with particle-based NERDSS results.
"""

from ionerdss.model import pdb
import subprocess
import os

# PDB ID for the example system
pdb_id = "6bno"

# Use local file (or can use PDB ID to download)
cif_path = "workspace_6BNO/structures/downloaded/6BNO.cif"

# Create model builder
model = pdb.PDBModelBuilder(source=cif_path)

# Configure hyperparameters with ODE pipeline enabled
model.set_hyperparameters(
    # Interface detection parameters
    interface_detect_distance_cutoff=1.0,
    ring_regularization_mode="off",
    
    # Enable ODE pipeline
    ode_enabled=True,
    ode_time_span=(0.0, 10.0),  # Simulation time in seconds
    ode_solver_method="BDF",     # Good for stiff systems
    ode_plot=True,               # Generate concentration plots
    ode_save_csv=True,           # Save results to CSV
    
    # Optional: Set custom initial concentrations
    # ode_initial_concentrations={'C1': 1.0, 'C2': 0.0}  # Start with monomer at 1.0 μM
)

# Build system with ODE calculation enabled
# This will:
# 1. Parse the PDB structure
# 2. Detect interfaces
# 3. Group chains
# 4. Build templates
# 5. Generate NERDSS files
# 6. Calculate ODE solution (NEW!)
# 7. Save all results
# Hyperparameters are automatically used from builder!
system = model.build_system(workspace_path="6bno_dir")

print("\n" + "="*60)
print("ODE Pipeline Completed!")
print("="*60)
print(f"Check the 'ode_results' directory in workspace for:")
print("  - ode_solution.csv: Time series data")
print("  - ode_solution.png: Concentration plots")
print("="*60)

# Optional: Run NERDSS simulation for comparison
print("\nRunning NERDSS simulation...")
nerdss_dir = "6bno_dir/nerdss_files"
nerdss_cmd = "~/Workspace/Reaction_ode/nerdss_development/bin/nerdss -f parms.inp"

# Change to the directory and run the command
subprocess.run(nerdss_cmd, shell=True, cwd=nerdss_dir, executable='/bin/bash')

print("\nWorkflow complete! You can now compare:")
print("  - ODE predictions: 6bno_dir/ode_results/ode_solution.csv")
print("  - NERDSS results: 6bno_dir/nerdss_files/...")
