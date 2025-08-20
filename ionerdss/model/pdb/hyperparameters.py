"""
Hyperparameters for PDB Models.

This module defines the `PDBModelHyperparameters` class, which centralizes 
all tunable hyperparameters used in PDB model processing. The class 
provides default values for parameters but allows overriding them 
through an options dictionary. 

Hyperparameter groups include:
    - **Structure coarse-graining**: distance and residue cutoffs, energy tables.
    - **Chain identification**: thresholds for RMSD and sequence similarity, 
      with optional custom aligners and matching modes.
    - **Chain regularization**: intra- and inter-chain distance thresholds 
      and angular thresholds for detecting repeated or symmetric interfaces.
    - **Output control**: toggling standard vs. verbose output.

Example:
    >>> from model.pdb.hyperparameters import PDBModelHyperparameters
    >>> hp = PDBModelHyperparameters({"distance_cutoff": 0.8, "verbose_mode": True})
    >>> hp.distance_cutoff
    0.8
    >>> hp.verbose_mode
    True
    >>> hp.seq_threshold  # uses default if not overridden
    0.5
"""

from ionerdss.model.pdb.energy_table import get_default_energy_table

class PDBModelHyperparameters():
    """
    Store hyperparameters and default values for PDB Model
    """

    def __init__(self, options=None):
        # Override default with emtpy dictionary
        if options is None:
            options = {}
        # Coarse-grain the structure
        self.distance_cutoff = options.get("distance_cutoff", 0.6)
        self.residue_cutoff = options.get("residue_cutoff", 3)
        self.energy_table = options.get("energy_table", None)
        # get default energy table if is set to None
        if self.energy_table is None:
            self.energy_table = get_default_energy_table()
        # Identify repeated chains
        self.rmsd_threshold = options.get("rmsd_threshold", 2.0)
        self.seq_threshold = options.get("seq_threshold", 0.5)
        self.custom_aligner = options.get("custom_aligner", None)
        self.matching_mode = options.get("matching_mode", "default")
        # Regularize repeated chains
        self.dist_threshold_intra = options.get("dist_threshold_intra", 3.5)
        self.dist_threshold_inter = options.get("dist_threshold_inter", 3.5)
        self.angle_threshold = options.get("angle_threshold", 25.0)
        # Output control
        self.standard_output = options.get("standard_output", False)
        self.verbose_mode = options.get("verbose_mode", False)
