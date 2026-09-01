"""
ionerdss.model.pdb.hyperparameters

Simulation hyperparameters for PDB to NERDSS conversion pipeline.

This module defines the PDBModelHyperparameters class that contains all
configurable parameters for the molecular model building process, including
distance cutoffs, thresholds, and algorithmic choices.

"""

from dataclasses import dataclass, field, fields
from typing import Optional, Literal
from Bio.Align import PairwiseAligner

from ionerdss.model.components.units import Units


@dataclass
class PDBModelHyperparameters:
    """Hyperparameters for PDB to NERDSS parameter pipeline.

    Contains all configurable parameters that control the molecular model
    building process, from interface detection to template generation.

    """

    # Core detection parameters
    interface_detect_distance_cutoff: float = field(
        default=0.9,
        metadata={"description": "Contact search radius per atom pair for interface detection", "unit": "nm"}
    )
    interface_detect_n_residue_cutoff: int = field(
        default=3,
        metadata={"description": "Minimum number of contacting residues (on each chain) to accept an interface", "unit": "residues"}
    )
    min_chain_length: int = field(
        default=4,
        metadata={"description": "Minimum number of residues for a chain to be included (filters out small molecules)", "unit": "residues"}
    )
    include_modified_residues: bool = field(
        default=True,
        metadata={"description": "Keep non-standard residues that belong to a polymer entity (modified amino acids such as D-amino acids or selenomethionine, and terminal caps). These are stored as HETATM and are invisible to BioPython's is_aa, but they often form the binding interfaces. Set False for the older standard-amino-acids-only behaviour."}
    )
    
    # Interface type assignment parameters (for template building)
    interface_type_assignment_distance_threshold: float = field(
        default=2.0,
        metadata={"description": "Distance threshold (Angstroms) for assigning interfaces to the same type during template building. Interfaces within this distance are merged into one type.", "unit": "Å"}
    )
    interface_type_assignment_angle_threshold: float = field(
        default=0.2,
        metadata={"description": "Angle threshold (radians) for assigning interfaces to the same type during template building. ~11 degrees. More restrictive than the previous 0.5 radians default.", "unit": "radians"}
    )
    split_spatial_interface_patches: bool = field(
        default=False,
        metadata={"description": "Split interfaces into spatial patches for faster template building"}
    )

    # Chain grouping parameters
    chain_grouping_rmsd_threshold: float = field(
        default=2.0,
        metadata={"description": "RMSD threshold for structure superposition to determine repeated chains", "unit": "Å"}
    )
    chain_grouping_seq_threshold: float = field(
        default=0.5,
        metadata={"description": "Sequence identity threshold for sequence alignment to determine repeated chains (0.5 = 50%)"}
    )
    chain_grouping_custom_aligner: Optional[PairwiseAligner] = field(
        default=None,
        metadata={"description": "Custom Bio.Align.PairwiseAligner for sequence alignment (None uses default settings)"}
    )
    chain_grouping_matching_mode: Literal["default", "sequence", "structure", "sequence_structure"] = field(
        default="default",
        metadata={"description": "Mode for determining repeated chains: 'default' (mmCIF header with sequence fallback), 'sequence' (sequence-based), 'structure' (structure-based), 'sequence_structure' (both must pass; separates quasi-equivalent conformers of one sequence)"}
    )

    # Steric clash detection
    steric_clash_mode: Literal["off", "auto", "custom"] = field(
        default="off",
        metadata={"description": "Mode for detecting steric clashes: 'off' (disabled), 'auto' (automatic Cα clash detection), 'custom' (user-provided lists)"}
    )

    # Template building parameters
    signature_precision: int = field(
        default=6,
        metadata={"description": "Number of decimal places for geometric signature normalization to avoid floating-point errors", "unit": "decimal places"}
    )
    homodimer_distance_threshold: float = field(
        default=0.5,
        metadata={"description": "Distance threshold for homodimer detection", "unit": "nm"}
    )
    homodimer_angle_threshold: float = field(
        default=0.5,
        metadata={"description": "Angle threshold for homodimer detection", "unit": "radians"}
    )

    # Enhanced homotypic detection parameters
    homotypic_detection: Literal["auto", "signature", "off"] = field(
        default="auto",
        metadata={"description": "Mode for homotypic binding detection: 'auto', 'signature', or 'off'"}
    )
    homotypic_detection_residue_similarity_threshold: float = field(
        default=0.7,
        metadata={"description": "Residue similarity threshold for homotypic detection (0.7 = 70% similarity)"}
    )
    homotypic_detection_interface_radius: float = field(
        default=8.0,
        metadata={"description": "Interface detection radius for homotypic binding", "unit": "Å"}
    )

    # Ring regularizer parameters
    is_on_sphere: bool = field(
        default=False,
        metadata={"description": "Project all molecules onto best-fit spheres centered on the most abundant element's center"}
    )

    # Chain regularizer parameters
    template_regularization_strength: float = field(
        default=0.0,
        metadata={"description": "Regularization strength for template fitting"}
    )

    # Visualizer options
    generate_visualizations: bool = field(
        default=True,
        metadata={"description": "Generate visualization outputs"}
    )

    # NERDSS file options
    generate_nerdss_files: bool = field(
        default=True,
        metadata={"description": "Generate NERDSS simulation files"}
    )
    nerdss_water_box: list = field(
        default_factory=lambda: [500.0, 500.0, 500.0],
        metadata={"description": "Water box dimensions for NERDSS simulation", "unit": "nm"}
    )

    nerdss_total_molecule_count: int = field(
        default=75,
        metadata={"description": "Total target number of molecules for NERDSS simulation, distributed according to stoichiometry", "unit": "molecules"}
    )
    nerdss_time_step: Optional[float] = field(
        default=None,
        metadata={"description": "Time step for NERDSS simulation (us). If None, calculated automatically based on stability criteria.", "unit": "us"}
    )
    
    nerdss_n_itr: int = field(
        default=int(1e5),
        metadata={"description": "Number of iterations for NERDSS simulation", "unit": "steps"}
    )
    
    nerdss_overlap_sep_limit: float = field(
        default=0.1,
        metadata={"description": "Minimum allowed separation distance between molecule centers of mass to prevent steric overlap. Increase this if flexible ring assemblies infinitely spiral.", "unit": "nm"}
    )
    
    disable_overlap_sep_limit_check: bool = field(
        default=False,
        metadata={"description": "Disable the safety check that automatically caps nerdss_overlap_sep_limit to 0.9 * minimum chain COM distance"}
    )
    
    # ProAffinity binding energy prediction options
    predict_affinity: bool = field(
        default=False,
        metadata={"description": "Enable ProAffinity-GNN binding affinity prediction"}
    )
    adfr_path: Optional[str] = field(
        default=None,
        metadata={"description": "Path to ADFR prepare_receptor tool (optional, will auto-detect if not provided)"}
    )
    proaffinity_backend: str = field(
        default="auto",
        metadata={"description": "Where to run ProAffinity: 'auto' (sidecar if one is configured, else in-process), 'in_process', or 'sidecar'. A sidecar keeps ProAffinity's numpy 1.x pin out of this environment."}
    )
    proaffinity_python: Optional[str] = field(
        default=None,
        metadata={"description": "Interpreter of the environment ProAffinity is installed in, used when running as a sidecar. Defaults to $IONERDSS_PROAFFINITY_PYTHON."}
    )
    
    # PDB file format options
    pdb_file_format: str = field(
        default="bioassembly1",
        metadata={"description": "Format for PDB file download: 'pdb', 'cif', 'mmcif', 'bioassembly1', 'bioassembly2', etc. (case-insensitive)"}
    )
    
    # ODE pipeline options
    ode_enabled: bool = field(
        default=False,
        metadata={"description": "Enable ODE pipeline for kinetic modeling"}
    )
    max_complex_size_ode: int = field(
        default=12,
        metadata={"description": "Maximum complex size (number of molecules) for ODE generation. ODE will be skipped if assembly exceeds this."}
    )
    ode_time_span: Optional[tuple] = field(
        default=None,
        metadata={"description": "Time span for ODE solving (start, end). If None, will be [0, nerdss_time_step * nerdss_n_itr].", "unit": "seconds"}
    )
    ode_solver_method: str = field(
        default="BDF",
        metadata={"description": "Solver method for stiff ODE systems (e.g., 'BDF', 'LSODA')"}
    )
    ode_atol: float = field(
        default=1e-4,
        metadata={"description": "Absolute tolerance for ODE solver"}
    )
    ode_plot: bool = field(
        default=True,
        metadata={"description": "Generate plots from ODE results"}
    )
    ode_save_csv: bool = field(
        default=True,
        metadata={"description": "Save ODE results to CSV file"}
    )
    ode_initial_concentrations: Optional[dict] = field(
        default=None,
        metadata={"description": "Custom initial concentrations for ODE (dict of species: concentration in uM). If not set, use the concentration equivalent to count / volume for nerdss simulation."}
    )

    ode_plot_species_indices: Optional[list] = field(
        default=None,
        metadata={"description": "List of species indices to plot specifically (None plots all or top N)"}
    )
    ode_plot_sample_points: int = field(
        default=1000,
        metadata={"description": "Number of time points to sample for plotting"}
    )
    ode_species_labels: Optional[dict] = field(
        default=None,
        metadata={"description": "Custom labels for species indices in plots"}
    )
    
    # Kinetic parameters
    default_on_rate_3d_ka: float = field(
        default=120.0,
        metadata={"description": "Default 3D association rate (ka) for diffusion-limited reactions", "unit": "nm^3/us"}
    )

    # Transition matrix output options
    count_transition: bool = field(
        default=False,
        metadata={"description": "Enable transition matrix tracking during NERDSS simulation"}
    )
    transition_matrix_size: Optional[int] = field(
        default=None,
        metadata={"description": "Size of transition matrix, i.e. the largest number of copies of a single molecule type NERDSS can track inside one complex. Must be >= the number of copies of that molecule type in the simulation, otherwise NERDSS writes past the end of the matrix and crashes. If None, it is set automatically to the molecule count of each type."}
    )
    transition_write: Optional[int] = field(
        default=None,
        metadata={"description": "Interval to write transition matrix (defaults to nItr/10)"}
    )
    
    # units
    units: Units = field(
        default_factory=Units,
        metadata={"description": "Unit system for the model (internal use, not serialized)"}
    )

    

    def __post_init__(self):
        """Initialize default aligner if not provided."""
        if self.chain_grouping_custom_aligner is None:
            self.chain_grouping_custom_aligner = self._create_default_aligner()

    def _create_default_aligner(self) -> PairwiseAligner:
        """Create default PairwiseAligner with standard parameters.

        Returns:
            Configured PairwiseAligner with global alignment mode and
            standard scoring parameters.
        """
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.5
        return aligner

    def to_dict(self) -> dict:
        """Convert hyperparameters to dictionary representation.

        Returns:
            Dictionary containing all hyperparameter values.
            Custom aligner is serialized as its parameter dictionary.
        """
        result = {}

        # Get all dataclass fields
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            # Handle special cases
            if field_name == 'chain_grouping_custom_aligner':
                if field_value is not None:
                    # Serialize aligner parameters
                    result[field_name] = {
                        'mode': getattr(field_value, 'mode', 'global'),
                        'match_score': getattr(field_value, 'match_score', 1.0),
                        'mismatch_score': getattr(field_value, 'mismatch_score', 0.0),
                        'open_gap_score': getattr(field_value, 'open_gap_score', -0.5),
                        'extend_gap_score': getattr(field_value, 'extend_gap_score', -0.5),
                    }
                else:
                    result[field_name] = None
            elif field_name == 'units':
                # Skip units field for JSON serialization
                continue
            else:
                # Regular field - just copy the value
                result[field_name] = field_value

        return result

    @classmethod
    def from_dict(cls, data: dict) -> "PDBModelHyperparameters":
        """Create hyperparameters from dictionary.

        Args:
            data: Dictionary containing hyperparameter values.

        Returns:
            New PDBModelHyperparameters instance.
        """
        if not data:
            return cls()

        # Get all valid field names from the dataclass
        valid_fields = {f.name for f in fields(cls)}

        # Filter out unknown keys and prepare data
        filtered_data = {}
        for key, value in data.items():
            if key in valid_fields:
                if key == 'chain_grouping_custom_aligner' and value is not None:
                    # Reconstruct aligner from parameters
                    aligner = PairwiseAligner()
                    if isinstance(value, dict):
                        for param, param_value in value.items():
                            if hasattr(aligner, param):
                                setattr(aligner, param, param_value)
                    filtered_data[key] = aligner
                elif key == 'ode_time_span' and isinstance(value, list):
                    # Convert list back to tuple (JSON serialization converts tuples to lists)
                    filtered_data[key] = tuple(value)
                else:
                    filtered_data[key] = value

        return cls(**filtered_data)

    def validate(self) -> list:
        """Validate hyperparameter values.

        Returns:
            List of validation error messages. Empty list if all valid.
        """
        errors = []
        
        # Validate enhanced homotypic detection parameters
        if self.homotypic_detection not in ["auto", "signature", "off"]:
            errors.append("homotypic_detection must be 'auto', 'signature', or 'off'")

        if not (0.0 <= self.homotypic_detection_residue_similarity_threshold <= 1.0):
            errors.append("residue_similarity_threshold must be between 0.0 and 1.0")

        if self.homotypic_detection_interface_radius <= 0:
            errors.append("interface_radius must be positive")

        # Validate distance_cutoff
        if self.interface_detect_distance_cutoff <= 0:
            errors.append("distance_cutoff must be positive")

        # Validate residue_cutoff
        if self.interface_detect_n_residue_cutoff < 1:
            errors.append("residue_cutoff must be at least 1")

        # Validate rmsd_threshold
        if self.chain_grouping_rmsd_threshold < 0:
            errors.append("rmsd_threshold must be non-negative")

        # Validate seq_threshold
        if not (0 <= self.chain_grouping_seq_threshold <= 1):
            errors.append("seq_threshold must be between 0 and 1")

        # Validate signature_precision
        if self.signature_precision < 0:
            errors.append("signature_precision must be non-negative")

        # Validate homodimer thresholds
        if self.homodimer_distance_threshold < 0:
            errors.append("homodimer_distance_threshold must be non-negative")

        if self.homodimer_angle_threshold < 0:
            errors.append("homodimer_angle_threshold must be non-negative")



        return errors

    def __str__(self) -> str:
        """Return string representation of hyperparameters."""
        return (f"PDBModelHyperparameters("
                f"homotypic_detection='{self.homotypic_detection}', "
                f"residue_similarity_threshold={self.homotypic_detection_residue_similarity_threshold}, "
                f"distance_cutoff={self.interface_detect_distance_cutoff}, "
                f"residue_cutoff={self.interface_detect_n_residue_cutoff}, "
                f"matching_mode='{self.chain_grouping_matching_mode}', "
                f"steric_clash_mode='{self.steric_clash_mode}')")

    def __repr__(self) -> str:
        """Return detailed representation of hyperparameters."""
        return self.__str__()
