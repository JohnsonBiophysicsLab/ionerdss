"""
core.py

High-level orchestration class PDBModel; imports helper functions from submodules.
"""

# model/pdb/core.py

import os
from . import io
from .coarse_grain import coarse_grain_structure
from .homolog_detection import identify_homologous_chains
from .homolog_alignment import regularize_model
from .reaction import build_binding_reactions
from ..components import Model  # assuming you have a base Model class

class PDBModel(Model):
    """Main driver for converting a PDB file into NERDSS molecule types and reactions."""

    def __init__(self, pdb_file=None, pdb_id=None, save_dir=None):
        super().__init__(save_dir)
        self.pdb_file = pdb_file
        self.pdb_id = pdb_id
        self.save_dir = os.path.abspath(os.path.expanduser(save_dir or os.getcwd()))
        os.makedirs(self.save_dir, exist_ok=True)

        # Step 1: Download if necessary
        if not pdb_file and pdb_id:
            self.pdb_file = io.download_pdb(pdb_id, self.save_dir)
        elif not self.pdb_file:
            raise ValueError("Must provide either pdb_file or pdb_id.")

        # Step 2: Parse structure
        self.structure = io.parse_structure(self.pdb_file)

        # Initialize metadata
        self.chains_map = {}
        self.chains_group = []

    def run_pipeline(self, options=None):
        """Runs the full model generation pipeline."""

        # 1. Coarse-grain the structure
        cg_result = coarse_grain_structure(self.structure, options=options)

        # 2. Identify homologous chains
        self.chains_map, self.chains_group = identify_homologous_chains(
            self.pdb_file, self.structure
        )

        # 3. Regularize molecules (alignment, interface generation)
        model_data = regularize_structure(
            cg_result, self.chains_map, self.chains_group
        )

        # 4. Generate reactions
        reactions = build_reactions(model_data)

        # 5. Rescale energies
        rescale_reaction_energies(reactions)

        # 6. Save model files
        self.save_model(self.pdb_id + "_model.json", model_data, reactions)

        # 7. Optionally plot or write CIFs
        if options and options.get("plot"):
            plot_structure(model_data, self.save_dir)
        if options and options.get("save_cif"):
            save_structure_outputs(model_data, self.save_dir)

        print("Pipeline completed.")
