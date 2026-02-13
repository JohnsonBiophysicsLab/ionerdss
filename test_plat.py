"""
platonic solids

"""


from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.parser import PDBParser

from ionerdss import build_system_from_plat

import subprocess
import os

pdb_id = "8y7s"

# Use local file
system, reactions = build_system_from_plat("dode", radius=10, sigma=3, output_nerdss=True, output_dir="./dode_dir")



##
# Step 2: Build system (now includes ODE calculation)
# Build NERDSS system and output NERDSS parameter files

## Run NERDSS simulation
print("\nNow running NERDSS simulation...\n")
# Change directory and run NERDSS simulation
nerdss_dir = f"./dode_dir/nerdss_files"
nerdss_cmd = "~/Workspace/Reaction_ode/nerdss_development/bin/nerdss -f parms.inp"

# Change to the directory and run the command
subprocess.run(nerdss_cmd, shell=True, cwd=nerdss_dir, executable='/bin/bash')
