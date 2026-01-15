#!/usr/bin/env python3
"""
example_platonic_solid.py

A complete example demonstrating how to:
1. Generate Platonic Solid models using PlatonicSolidsModel.
2. Inspect the generated System and ReactionRules.
3. Export the models to NERDSS format (.mol and parms.inp).

Usage:
    python example_platonic_solid.py
"""

import os
from ionerdss.model.PlatonicSolids import PlatonicSolidsModel

def main():
    # Define output directory
    output_dir = "nerdss_output"
    
    # 1. Create a Cube
    print("--- Generating Cube ---")
    cube_system, cube_reactions = PlatonicSolidsModel.create_solid(
        solid_type="cube", 
        radius=10.0,  # nm
        sigma=1.0     # nm
    )
    
    print(f"Generated System for Cube:")
    print(f"  Molecule Types: {len(cube_system.molecule_types)}")
    print(f"  Molecule Instances: {len(cube_system.molecule_instances)}")
    print(f"  Interface Types: {len(cube_system.interface_types)}")
    print(f"  Reactions Generated: {len(cube_reactions)}")
    
    # Export Cube
    cube_out = os.path.join(output_dir, "cube_sim")
    print(f"Exporting Cube to '{cube_out}'...")
    PlatonicSolidsModel.export_nerdss(cube_system, cube_out, cube_reactions)
    print("Export complete.\n")

    # 2. Create a Dodecahedron (more complex)
    print("--- Generating Dodecahedron ---")
    dode_system, dode_reactions = PlatonicSolidsModel.create_solid(
        solid_type="dode",
        radius=15.0,
        sigma=1.5
    )
    
    print(f"Generated System for Dodecahedron:")
    print(f"  Molecule Types: {len(dode_system.molecule_types)}")
    print(f"  Molecule Instances: {len(dode_system.molecule_instances)}")
    print(f"  Interface Types: {len(dode_system.interface_types)}")
    print(f"  Reactions Generated: {len(dode_reactions)}")

    # Export Dodecahedron
    dode_out = os.path.join(output_dir, "dode_sim")
    print(f"Exporting Dodecahedron to '{dode_out}'...")
    PlatonicSolidsModel.export_nerdss(dode_system, dode_out, dode_reactions)
    print("Export complete.\n")

    print(f"All examples finished. Check the '{output_dir}' directory for outputs.")

if __name__ == "__main__":
    main()
