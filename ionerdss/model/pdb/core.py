"""
core.py

High-level orchestration class PDBModel; imports helper functions from submodules.
"""

import os
from . import io
from .coarse_grain import coarse_grain_structure
from .detect_repeats import identify_repeated_chains
from .repeated_chain_alignment import regularize_model
from .reaction import build_binding_reactions
from .capsid_sphere_pipeline import run_capsid_pipeline  # new module for spherical gag pipeline
from ..components import Model  # assuming you have a base Model class

class PDBModel(Model):
    """Main driver for converting a PDB file into NERDSS molecule types and reactions."""

    def __init__(self, pdb_file=None, pdb_id=None, save_dir=None):
        """
        Parameters
        ----------
        pdb_file : str or None
            Path to the PDB or CIF file.
        pdb_id : str or None
            RCSB PDB ID to fetch if pdb_file is not provided.
        save_dir : str or None
            Directory to save generated model files.
        """
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
        """Runs the full model generation pipeline.

        Parameters
        ----------
        options : dict or None
            Optional dictionary of settings. Can include:
            - 'plot': bool
            - 'save_cif': bool
            - 'is_on_sphere': bool (if True, run spherical capsid pipeline)
        """
        options = options or {}



        # 1. Coarse-grain the structure
        cg_result = coarse_grain_structure(self.structure, options=options)

        # 2. Identify repeated chains
        self.chains_map, self.chains_group = identify_repeated_chains(
            self.pdb_file, self.structure
        )

        # 3. Regularize molecules (alignment, interface generation)
        model_data = regularize_model(
            cg_result, self.chains_map, self.chains_group
        )

        # 4. Generate reactions
        reactions = build_binding_reactions(model_data)

        # 5. Rescale energies (optional)
        # rescale_reaction_energies(reactions)  # Uncomment if needed

        # If modeling a spherical capsid, use the alternate pipeline
        if options.get("is_on_sphere"):
            run_spherical_pipeline(self.pdb_file, self.save_dir, options)
            print("Capsid pipeline completed.")
            return

        # 6. Save model files
        self.save_model(self.pdb_id + "_model.json", model_data, reactions)

        # 7. Optionally plot or write CIFs
        if options.get("plot"):
            plot_structure(model_data, self.save_dir)
        if options.get("save_cif"):
            save_structure_outputs(model_data, self.save_dir)

        print("Pipeline completed.")
