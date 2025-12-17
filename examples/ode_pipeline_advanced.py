"""
Advanced Example: Direct ODE Pipeline Usage

This script demonstrates how to use the ODE pipeline functions directly
for more customized analysis and visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from ionerdss import ParseComplexes, ODEPipelineConfig, calculate_ode_solution
from ionerdss.model import pdb

# Build the system (without automatic ODE calculation)
model = pdb.PDBModelBuilder(source="6bno")
model.set_hyperparameters(
    interface_detect_distance_cutoff=1.0,
    ring_regularization_mode="off",
    generate_nerdss_files=True,
    ode_enabled=False  # We'll run ODE manually
)

system = model.build_system(workspace_path="6bno_advanced")

# Generate complex reaction system using ParseComplexes
print("Generating complex reaction network...")
complex_list, complex_reaction_system = ParseComplexes(system)

print(f"\nFound {len(complex_list)} complex species")
print(f"Generated {len(complex_reaction_system.reactions)} reactions")

# Print reaction network
print("\n" + "="*60)
print("Reaction Network:")
print("="*60)
for i, reaction in enumerate(complex_reaction_system.reactions):
    print(f"{i+1}. {reaction.expression} (rate = {reaction.rate})")

# Configure ODE calculation with custom settings
ode_config = ODEPipelineConfig(
    t_span=(0.0, 20.0),           # Longer time span
    solver_method="BDF",
    atol=1e-6,                     # Tighter tolerance
    plot=False,                    # We'll make custom plots
    save_csv=False,
    initial_concentrations=None    # Default: monomer at 1.0
)

# Calculate ODE solution
print("\nSolving ODE system...")
time, concentrations, species_names = calculate_ode_solution(
    complex_reaction_system,
    config=ode_config
)

print(f"Solved for {len(time)} time points")
print(f"Species: {species_names}")

# Create custom visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: All species
for i, species in enumerate(species_names):
    ax1.plot(time, concentrations[:, i], label=species, linewidth=2, alpha=0.7)

ax1.set_xlabel('Time (s)', fontsize=12)
ax1.set_ylabel(r'Concentration $(\mu\mathrm{M})$', fontsize=12)
ax1.set_title('All Species Concentrations', fontsize=14)
ax1.legend(loc='best', fontsize=8)
ax1.grid(True, alpha=0.3)

# Plot 2: Selected species or aggregated view
# Example: Plot monomer vs sum of all multimers
monomer_conc = concentrations[:, 0]  # First species is typically monomer
multimer_conc = concentrations[:, 1:].sum(axis=1)  # Sum all other species

ax2.plot(time, monomer_conc, label='Monomer', linewidth=3)
ax2.plot(time, multimer_conc, label='All Multimers', linewidth=3)
ax2.set_xlabel('Time (s)', fontsize=12)
ax2.set_ylabel(r'Concentration $(\mu\mathrm{M})$', fontsize=12)
ax2.set_title('Monomer vs Multimers', fontsize=14)
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('6bno_advanced/ode_custom_analysis.png', dpi=300, bbox_inches='tight')
print(f"\nCustom plots saved to: 6bno_advanced/ode_custom_analysis.png")

# Calculate some metrics
print("\n" + "="*60)
print("Analysis Metrics:")
print("="*60)
print(f"Initial monomer concentration: {monomer_conc[0]:.4f} μM")
print(f"Final monomer concentration: {monomer_conc[-1]:.4f} μM")
print(f"Monomer depletion: {(monomer_conc[0] - monomer_conc[-1])/monomer_conc[0]*100:.2f}%")
print(f"Final total multimer concentration: {multimer_conc[-1]:.4f} μM")

# Find dominant species at equilibrium
final_concentrations = concentrations[-1, :]
dominant_idx = np.argmax(final_concentrations)
print(f"\nDominant species at t={time[-1]:.1f}s:")
print(f"  {species_names[dominant_idx]}: {final_concentrations[dominant_idx]:.4f} μM")

plt.show()
