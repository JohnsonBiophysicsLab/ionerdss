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
        self.logger_level = options.get("logger_level", "INFO")

    def __repr__(self):
        """
        Developer-oriented representation. Should be unambiguous and ideally valid Python.
        """
        attrs = (
            f"distance_cutoff={self.distance_cutoff}, "
            f"residue_cutoff={self.residue_cutoff}, "
            f"energy_table={self.energy_table!r}, "
            f"rmsd_threshold={self.rmsd_threshold}, "
            f"seq_threshold={self.seq_threshold}, "
            f"custom_aligner={self.custom_aligner!r}, "
            f"matching_mode={self.matching_mode!r}, "
            f"dist_threshold_intra={self.dist_threshold_intra}, "
            f"dist_threshold_inter={self.dist_threshold_inter}, "
            f"angle_threshold={self.angle_threshold}, "
            f"standard_output={self.standard_output}, "
            f"logger_level={self.logger_level}"
        )
        return f"{self.__class__.__name__}({attrs})"

    def __str__(self):
        """
        User-friendly string for pretty printing (multi-line).
        """
        return (
            f"PDBModelHyperparameters:\n"
            f"  distance_cutoff      = {self.distance_cutoff}\n"
            f"  residue_cutoff       = {self.residue_cutoff}\n"
            f"  energy_table         = {self.energy_table}\n"
            f"  rmsd_threshold       = {self.rmsd_threshold}\n"
            f"  seq_threshold        = {self.seq_threshold}\n"
            f"  custom_aligner       = {self.custom_aligner}\n"
            f"  matching_mode        = {self.matching_mode}\n"
            f"  dist_threshold_intra = {self.dist_threshold_intra}\n"
            f"  dist_threshold_inter = {self.dist_threshold_inter}\n"
            f"  angle_threshold      = {self.angle_threshold}\n"
            f"  standard_output      = {self.standard_output}\n"
            f"  logger_level         = {self.logger_level}"
        )
