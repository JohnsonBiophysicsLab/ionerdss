"""
ionerdss.model.pdb.hyperparameters

Simulation hyperparameters for PDB to NERDSS conversion pipeline.

This module defines the PDBModelHyperparameters class that contains all
configurable parameters for the molecular model building process, including
distance cutoffs, thresholds, and algorithmic choices.

## Hyperparameters Reference

| Hyperparameter | Definition | Default Value |
|----------------|------------|---------------|
| **Core Detection Parameters** |
| `distance_cutoff` | Contact search radius per atom pair for interface detection | 0.6 nm |
| `residue_cutoff` | Minimum number of contacting residues (on each chain) to accept an interface | 3 residues |
| **Chain Grouping Parameters** |
| `rmsd_threshold` | RMSD threshold for structure superposition to determine repeated chains | 2.0 nm |
| `seq_threshold` | Sequence identity threshold for sequence alignment to determine repeated chains | 0.5 (50%) |
| `custom_aligner` | Custom Bio.Align.PairwiseAligner for sequence alignment (None uses default settings) | None |
| `matching_mode` | Mode for determining repeated chains: "default" (mmCIF header with sequence fallback), "sequence" (sequence-based), "structure" (structure-based) | "default" |
| **Steric Clash Detection** |
| `steric_clash_mode` | Mode for detecting steric clashes: "off" (disabled), "auto" (automatic Cα clash detection), "custom" (user-provided lists) | "off" |
| **Template Building Parameters** |
| `signature_precision` | Number of decimal places for geometric signature normalization to avoid floating-point errors | 6 decimal places |
| `homodimer_distance_threshold` | Distance threshold for homodimer detection | 0.1 Å |
| `homodimer_angle_threshold` | Angle threshold for homodimer detection | 0.1 radians |
| **Ring Regularization Parameters** |
| `ring_regularization_mode` | Ring structure regularization mode: "off" (disabled), "separate" (individual ring fitting), "uniform" (single fit for all rings) | "uniform" |
| `ring_geometry` | Target geometry for ring regularization: "cylinder" or "sphere" | "cylinder" |
| `min_ring_size` | Minimum number of subunits required to form a ring | 3 subunits |

## Usage Examples

### Basic Usage
```python
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

# Use default parameters
params = PDBModelHyperparameters()

# Customize specific parameters
params = PDBModelHyperparameters(
    distance_cutoff=0.8,      # Looser contact detection
    residue_cutoff=5,         # Require more contacts
    matching_mode="sequence"  # Force sequence-based grouping
)
```

### Configuration Management
```python
# Save configuration
config_dict = params.to_dict()

# Load configuration
params = PDBModelHyperparameters.from_dict(config_dict)

# Validate parameters
errors = params.validate()
if errors:
    print("Configuration errors:", errors)
```

### Common Parameter Sets

**High-Resolution Structures:**
```python
high_res_params = PDBModelHyperparameters(
    distance_cutoff=0.5,      # Tight contacts
    residue_cutoff=5,         # Substantial interfaces
    rmsd_threshold=1.0,       # Strict structural similarity
    seq_threshold=0.9         # High sequence identity
)
```

**Low-Resolution Structures:**
```python
low_res_params = PDBModelHyperparameters(
    distance_cutoff=1.2,      # Loose contacts
    residue_cutoff=3,         # Minimal interfaces
    rmsd_threshold=5.0,       # Permissive structural similarity
    seq_threshold=0.3         # Low sequence identity
)
```

**Ring Structure Processing:**
```python
ring_params = PDBModelHyperparameters(
    ring_regularization_mode="separate",  # Individual ring fitting
    ring_geometry="sphere",               # Spherical geometry
    min_ring_size=4                       # Require 4+ subunits
)
```

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
    interface_detect_distance_cutoff: float = 0.6  # nm
    interface_detect_n_residue_cutoff: int = 3

    # Chain grouping parameters
    chain_grouping_rmsd_threshold: float = 2.0  # A
    chain_grouping_seq_threshold: float = 0.5
    chain_grouping_custom_aligner: Optional[PairwiseAligner] = field(default=None)
    chain_grouping_matching_mode: Literal["default", "sequence", "structure"] = "default"

    # Steric clash detection
    steric_clash_mode: Literal["off", "auto", "custom"] = "off"

    # Template building parameters
    signature_precision: int = 6
    homodimer_distance_threshold: float = 0.5  # nm
    homodimer_angle_threshold: float = 0.5  # radians

    # Enhanced homotypic detection parameters
    homotypic_detection: Literal["auto", "signature", "off"] = "auto"
    homotypic_detection_residue_similarity_threshold: float = 0.7  # 70% similarity
    homotypic_detection_interface_radius: float = 8.0  # A

    # Ring regularizer parameters
    ring_regularization_mode: str = "uniform"  # "off", "separate", "uniform"
    ring_geometry: str = "cylinder"  # "cylinder", "sphere"
    min_ring_size: int = 3

    # Chain regularizer parameters
    template_regularization_strength: float = 0.0

    # Visualizer options
    generate_visualizations: bool = True

    # NERDSS file options
    generate_nerdss_files: bool = True
    
    # ProAffinity binding energy prediction options
    predict_affinity: bool = False  # Enable ProAffinity-GNN prediction
    adfr_path: Optional[str] = None  # Path to ADFR prepare_receptor tool
    
    # ODE pipeline options
    ode_enabled: bool = False
    ode_time_span: tuple = (0.0, 10.0)  # (start, end) in seconds
    ode_solver_method: str = "BDF"  # Solver method for stiff systems
    ode_atol: float = 1e-4  # Absolute tolerance
    ode_plot: bool = True  # Generate plots
    ode_save_csv: bool = True  # Save results to CSV
    ode_initial_concentrations: Optional[dict] = None  # Custom initial concentrations
    
    # Transition matrix output options
    count_transition: bool = False  # Enable transition matrix tracking
    transition_matrix_size: int = 500  # Size of transition matrix
    transition_write: Optional[int] = None  # Interval to write transition matrix (defaults to nItr/10)
    
    # units
    units = Units()

    

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
            if field_name == 'custom_aligner':
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
                if key == 'custom_aligner' and value is not None:
                    # Reconstruct aligner from parameters
                    aligner = PairwiseAligner()
                    if isinstance(value, dict):
                        for param, param_value in value.items():
                            if hasattr(aligner, param):
                                setattr(aligner, param, param_value)
                    filtered_data[key] = aligner
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
