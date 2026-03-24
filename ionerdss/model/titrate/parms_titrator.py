#!/usr/bin/env python3.10
"""
This script reads the parms.inp file, finds the molecules, 
and generates a new parms.inp file with titration reactions added.
"""

import os
import re
import random
import argparse

def parse_mol_file(mol_file_path):
    """
    Parses a NERDSS .mol file to extract binding interfaces.
    Returns a list of interface names.
    """
    interfaces = []
    try:
        with open(mol_file_path, 'r') as f:
            lines = f.readlines()
            
        in_coords_section = False
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                # Check if we are entering the coords section
                if "coordinates" in line.lower() or "coords" in line.lower():
                    in_coords_section = True
                continue
                
            # If we hit another section header like # bonds, stop parsing coords
            if line.startswith('#') and in_coords_section:
                in_coords_section = False
                continue
                
            if in_coords_section or len(line.split()) >= 4:
                # Typically, coordinate lines have 4 columns: Name X Y Z
                parts = line.split()
                if len(parts) >= 1:
                    name = parts[0]
                    # Exclude COM and REF
                    if name.lower() not in ['com', 'ref'] and not name.startswith('#'):
                        # Additional safety check to avoid grabbing non-coordinate lines
                        try:
                            if len(parts) >= 4:
                                float(parts[1])
                                float(parts[2])
                                float(parts[3])
                                interfaces.append(name)
                        except ValueError:
                            pass
                            
    except FileNotFoundError:
        print(f"Warning: {mol_file_path} not found.")
    
    return interfaces

def generate_titration_reactions(parms_file, output_file, rate_min=0.00001, rate_max=0.0001):
    """
    Reads parms.inp, extracts molecules, finds their .mol files, 
    generates titration reactions, and writes to parms_titrate.inp.
    """
    base_dir = os.path.dirname(os.path.abspath(parms_file))
    if not base_dir:
        base_dir = "."
        
    molecules = []
    
    # Read parms.inp
    with open(parms_file, 'r') as f:
        parms_lines = f.readlines()
        
    # Extract molecules
    in_molecules_section = False
    for line in parms_lines:
        line_strip = line.strip()
        if line_strip == "start molecules":
            in_molecules_section = True
            continue
        elif line_strip == "end molecules":
            in_molecules_section = False
            break
            
        if in_molecules_section and line_strip and not line_strip.startswith('#'):
            # Format is typically "A0 : 29"
            mol_name = line_strip.split(':')[0].strip()
            molecules.append(mol_name)
            
    print(f"Found molecules: {molecules}")
    
    # Generate titration reactions text
    titration_reactions_text = ["\n    # Titration Reactions\n"]
    
    for mol in molecules:
        mol_file = os.path.join(base_dir, f"{mol}.mol")
        interfaces = parse_mol_file(mol_file)
        
        # Format the reaction
        # 0 -> A0(a0b01, a0c1, a0g1, a0h1, a0i1)
        if interfaces:
            interfaces_str = ", ".join(interfaces)
            reaction_str = f"    0 -> {mol}({interfaces_str})"
        else:
            reaction_str = f"    0 -> {mol}()"
            
        titration_reactions_text.append(reaction_str)
        
        # Generate random rate
        rate = random.uniform(rate_min, rate_max)
        titration_reactions_text.append(f"    onRate3Dka = {rate:.6f} # M/s\n")
        
    # Write new file
    with open(output_file, 'w') as f:
        for line in parms_lines:
            f.write(line)
            
            if line.strip() == "start reactions":
                # Write titration reactions immediately after "start reactions"
                for titr_line in titration_reactions_text:
                    f.write(titr_line + "\n")
                    
    print(f"Successfully wrote {output_file} with {len(molecules)} titration reactions.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add titration reactions to NERDSS parms.inp")
    parser.add_argument("input_file", help="Path to input parms.inp")
    parser.add_argument("-o", "--output", default="parms_titrate.inp", help="Output file name (default: parms_titrate.inp)")
    parser.add_argument("--min-rate", type=float, default=0.00001, help="Minimum onRate3Dka (default: 0.00001)")
    parser.add_argument("--max-rate", type=float, default=0.0001, help="Maximum onRate3Dka (default: 0.0001)")
    
    args = parser.parse_args()
    
    # Handle output path to be in same dir as input if only filename is provided
    output_path = args.output
    if not os.path.dirname(output_path):
        output_path = os.path.join(os.path.dirname(os.path.abspath(args.input_file)), args.output)
        
    generate_titration_reactions(args.input_file, output_path, args.min_rate, args.max_rate)
