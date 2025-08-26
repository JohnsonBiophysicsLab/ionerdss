"""
core.py

High-level orchestration class PDBModel; imports helper functions from submodules.
"""

import os
import logging

from ionerdss.model.graph_based.complexes.graphize import (
    networkx_graph_to_string,
    build_simple_graph,
)
from ionerdss.model.graph_based.complexes.subcomplexes import (
    get_unique_fully_connected_subgraphs,
)
from ionerdss.model.graph_based.reactions.dimer import (
    find_all_dimer_reactions,
    get_broken_edges,
)
from ionerdss.model.graph_based.reactions.transformation import (         
    find_all_transformable_subgraph_pairs,
)

from . import io
from .coarse_grain import coarse_grain_structure
from .detect_repeats import identify_repeated_chains
from .regularize_repeats import regularize_repeated_chains
from .hyperparameters import PDBModelHyperparameters
from .visualize import plot_coarse_grain_model, save_coarse_grained_structure
from ..components import Model  # assuming you have a base Model class

# get module level logger
logger = logging.getLogger("ionerdss.model.pdb")       # module-level logger

class PDBModel(Model):
    """Main driver for converting a PDB file into NERDSS molecule types and reactions."""

    def __init__(self, pdb_file=None, pdb_id=None, save_dir=None,
                 auto_run=True,
                 options=None,
                 **kwargs):
        """
        Parameters
        ----------
        pdb_file : str or None
            Path to the PDB or CIF file.
        pdb_id : str or None
            RCSB PDB ID to fetch if pdb_file is not provided.
        save_dir : str or None
            Directory to save generated model files.
        options :
            Provide a dictionary or a PDBModelHyperparameters instance for pipeline
            hyperparameters:
        auto_run : bool (default : True)
            automaitcally run pipleline if set to true
        """
        # -------------------------------------------------------------------
        # Initialize parent
        # -------------------------------------------------------------------
        super().__init__(save_dir)
        self.pdb_file = pdb_file
        self.pdb_id = pdb_id
        self.save_dir = os.path.abspath(
            os.path.expanduser(save_dir or os.getcwd()))
        os.makedirs(self.save_dir, exist_ok=True)

        # -------------------------------------------------------------------
        # Configure logging once in your main script (NOT in every module)
        # -------------------------------------------------------------------
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        )

        # Step 0: Setup hyperparameters
        # if provided a hyperparameter class, directly use that
        if isinstance(options, PDBModelHyperparameters):
            self.params = options
        # Validate dict/None path
        elif options is not None and not isinstance(options, dict):
            raise TypeError(
            "hyperparameters only supports None, dict, or PDBModelHyperparameters"
        )
        else:
            # merge dict-style options and keyword-style options
            opts = {**(options or {}), **kwargs}
            self.params = PDBModelHyperparameters(options=opts)

        # set verbose mode if enabled in hyperparameters
        if self.params.logger_level == "DEBUG":
            logger.setLevel(logging.DEBUG)     # enable verbose per-run
        else:
            logger.setLevel(logging.INFO)      # disable verbose per-run

        # Step 1: Download if necessary
        if not pdb_file and pdb_id:
            self.pdb_file = io.download_pdb(pdb_id, self.save_dir)
        elif not self.pdb_file:
            raise ValueError("Must provide either pdb_file or pdb_id.")

        # Step 2: Parse structure
        self.structure = io.parse_structure(self.pdb_file)

        # Initialize metadata
        self.chains_map = {}
        self.chains_groups = []

        # Step 3: Auto run if set to True
        if auto_run:
            self.run_pipeline()

    def run_pipeline(self,
                     do_plot = False):
        """Runs the full model generation pipeline.

        Parameters
        ----------
        options : dict or None
            Optional dictionary of settings. Can include:
            - 'plot': bool
            - 'save_cif': bool
            - 'is_on_sphere': bool (if True, run spherical capsid pipeline)
        """

        logger.info("Starting Coarse-graining the structure...")

        # 1. Coarse-grain the structure
        cg_model = coarse_grain_structure(self.structure,
                                          params=self.params)

        logger.debug("Generated cg_model :")
        logger.debug(cg_model)

        # returned cg model is a dictionary with the following k,v pairs
        # {
        # 'chains': chains,
        # 'COMs': coms,
        # 'radii': chain_radii,
        # 'interfaces': interfaces,
        # 'interface_coords': interface_coords,
        # 'interface_residues': interface_residues,
        # 'interface_energies': interface_energies,
        # }

        # 2. Identify repeated chains
        self.chains_map, self.chains_groups = identify_repeated_chains(
            self.pdb_file, self.structure,
            params=self.params
        )

        # 3. Regularize molecules (alignment, interface generation)
        regularized_model_data = regularize_repeated_chains(
            cg_model, self.chains_map, self.chains_groups,
            params=self.params
        )
        updated_cg_model = regularized_model_data["updated_cg_model"]

        # logging
        logger.debug("Generated regularized model data : ")
        logger.debug(regularized_model_data)
        logger.debug(".....................")

        for it in regularized_model_data["interface_templates"]:
            logger.debug("Interface template %s", it.name)  # Expect nonzero Coords object
            logger.debug("IT Coordinate %s", it.coords)

        # Generate a regularized model while keeping the original
        # cg_model intact via deep copy
        #regularize_cg_model = cg_model.copy()
        #regularize_cg_model["COMs"] = [
        #    mol.coord for mol in regularized_model_data["molecules"]]
        #regularize_cg_model["interface_coords"] = [
        #    mol.interface_list for mol in regularized_model_data["molecules"]]

        # Draw plot if prompted
        logger.debug("Interface coordinates : ")
        logger.debug(updated_cg_model["interface_coords"])

        if do_plot:
            plot_coarse_grain_model(updated_cg_model)
        save_coarse_grained_structure(
            updated_cg_model, self.save_dir, self.pdb_file)

        # 4. Generate reactions
        
        # Suppose you have cg_model (from your coarse graining),
        # and a chains_map like {"A":"A","B":"A",...}
        G = build_simple_graph(cg_model, chains_map=self.chains_map)

        # If you want 5l93-like edge labels:
        # annotate_edges_by_cycles(G)

        # debug output
        logger.debug(networkx_graph_to_string(G))

        # get all subspecies
        species = get_unique_fully_connected_subgraphs(G)
        logger.debug(species)

        reactions = find_all_dimer_reactions(species, use_multiprocessing=False)

        reactions_list = [
            {
                "product": list(r[2].nodes),
                "part1": list(r[0]),
                "part2": list(r[1]),
                "bonds_broken": get_broken_edges(r[2], r[0], r[1])
            }
            for r in reactions
        ]

        logger.debug(reactions_list)
        
        transformations = find_all_transformable_subgraph_pairs(G, species)

        transformations_list = [
            {
                "monomer_1_nodes": list(t1.nodes),
                "monomer_2_nodes": list(t2.nodes),
                "diff": list(set(t1.edges(data="type")) ^ set(t2.edges(data="type")))
            }
            for t1, t2 in transformations
        ]

        logger.debug(transformations_list)


        # 5. Rescale energies (optional)
        # rescale_reaction_energies(reactions)  # Uncomment if needed

        # If modeling a spherical capsid, use the alternate pipeline
        # if options.get("is_on_sphere"):
        #    run_spherical_pipeline(self.pdb_file, self.save_dir, options)
        #    print("Capsid pipeline completed.")
        #    return

        # 6. Save model files
        # self.save_model(self.pdb_id + "_model.json", regularized_model_data, reactions)

        # 7. Optionally plot or write CIFs
        # if options.get("plot"):
        #    plot_structure(regularized_model_data, self.save_dir)
        # if options.get("save_cif"):
        #    save_structure_outputs(regularized_model_data, self.save_dir)

        logger.info("Pipeline completed.")
