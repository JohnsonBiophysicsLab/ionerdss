"""
ionerdss.model.pdb.coarse_graining

Interface detection and coarse-grained representation generation.

This module implements the coarse-graining pipeline that converts all-atom
protein structures into simplified representations suitable for NERDSS
modeling. It performs interface detection via KD-tree spatial queries and
generates coarse-grained molecular representations. It converts all-atom protein
structures into simplified representations suitable for NERDSS molecular
simulations by detecting binding interfaces through spatial proximity analysis
and generating coarse-grained molecular models.

## Key Concepts

### Coarse-Graining
**Coarse-graining** in this code is referring to the process of simplifying
complex molecular structures by reducing atomic-level detail while preserving 
essential geometric and interaction properties. In ionerdss, this involves:

- Converting protein chains into single particles with effective radii
- Identifying binding interfaces as discrete interaction sites
- Representing molecular assemblies as networks of interacting particles

### Interface Detection
**Interface detection** identifies regions where protein chains come into close
contact and can potentially bind. The algorithm uses:

- **Spatial proximity**: Cα atoms within a distance cutoff
- **Contact persistence**: Minimum number of contacting residues
- **Geometric characterization**: Interface coordinates and contact maps
- **Residue composition**: Amino acid types and positions at interfaces

### Partner Mapping
**Partner mapping** creates explicit relationships between binding interfaces,
enabling the construction of molecular interaction networks for simulation setup.

## Classes and Data Structures

### `ResidueInfo`

Represents detailed information about a residue at an interface.

```python
@dataclass
class ResidueInfo:
    residue_id: int           # Residue sequence number
    residue_name: str         # Three-letter amino acid code (e.g., "ALA")
    position: np.ndarray      # Cα coordinate position (Å)
    chain_id: str            # Chain identifier
```

### `InterfaceString`

Represents a detected binding interface between two protein chains with detailed
residue composition information.

```python
@dataclass
class InterfaceString:
    chain_i: str                        # First chain ID
    chain_j: str                        # Second chain ID  
    coord_i: np.ndarray                 # Interface coordinate on chain i (Å)
    coord_j: np.ndarray                 # Interface coordinate on chain j (Å)
    residues_i: Set[int]                # Contacting residue IDs on chain i
    residues_j: Set[int]                # Contacting residue IDs on chain j
    residue_details_i: List[ResidueInfo] # Detailed residue info for chain i
    residue_details_j: List[ResidueInfo] # Detailed residue info for chain j
    energy: float = -1.0                # Binding energy (default -1.0)
```

### `CoarseGrainedChain`

Coarse-grained representation of a protein chain containing essential geometric
and topological information. Interface references are managed separately to avoid
JSON serialization issues.

```python
@dataclass  
class CoarseGrainedChain:
    chain_id: str                      # Chain identifier
    com: np.ndarray                    # Center of mass (Å)
    radius: float                      # Effective radius (Å)
    sequence: str                      # Amino acid sequence
    bbox_min: np.ndarray              # Bounding box minimum (Å)
    bbox_max: np.ndarray              # Bounding box maximum (Å)
```

### `CoarseGrainer`

Main engine that orchestrates the coarse-graining pipeline.

```python
class CoarseGrainer:
    def __init__(self, parser: PDBParser, hyperparams: PDBModelHyperparameters)
```

**Attributes:**
- `parser`: PDB parser containing structure data
- `hyperparams`: Configuration parameters
- `chains`: Dictionary of coarse-grained chain representations
- `interfaces`: List of all detected interfaces
- `partner_map`: Mapping of interface partnerships

## Enhanced Interface Detection

The interface detection now captures detailed amino acid composition:

```python
# For each contacting residue, store:
residue_info = ResidueInfo(
    residue_id=101,           # Sequence number
    residue_name="ALA",       # Alanine
    position=[x, y, z],       # Cα coordinates
    chain_id="A"              # Chain identifier
)
```

This enables downstream residue-based homotypic interface validation.

## JSON Serialization

All data structures are designed to be JSON serializable to avoid system export issues.
InterfaceString objects are stored separately and accessed through helper methods.

"""

from typing import Any, Dict, List, Tuple, Set, Optional
from dataclasses import dataclass
import logging

import numpy as np
from scipy.spatial import KDTree
from scipy.cluster.hierarchy import fcluster, linkage

from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser

logger = logging.getLogger(__name__)

@dataclass
class ResidueInfo:
    """Detailed information about a residue at an interface.

    Attributes:
        residue_id: Residue sequence number.
        residue_name: Three-letter amino acid code (e.g., "ALA").
        position: Cα coordinate position in Angstroms.
        chain_id: Chain identifier.
    """
    residue_id: int
    residue_name: str
    position: np.ndarray
    chain_id: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert ResidueInfo to JSON-serializable dictionary."""
        return {
            'residue_id': self.residue_id,
            'residue_name': self.residue_name,
            'position': self.position.tolist() if isinstance(self.position, np.ndarray) else self.position,
            'chain_id': self.chain_id
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ResidueInfo':
        """Create ResidueInfo from dictionary."""
        return cls(
            residue_id=data['residue_id'],
            residue_name=data['residue_name'],
            position=np.array(data['position']),
            chain_id=data['chain_id']
        )

    def __str__(self) -> str:
        """Return human-readable representation."""
        return f"{self.residue_name}{self.residue_id}"

    def __repr__(self) -> str:
        """Return detailed representation."""
        return f"ResidueInfo({self.residue_name}{self.residue_id}, chain={self.chain_id})"


@dataclass
class InterfaceString:
    """Detected interface between two chains with detailed residue information.

    Separate from interface instance and use simple string to
    avoid assigning complicated crosslink at first

    Represents a validated binding interface with geometric and
    topological information for downstream processing.

    Attributes:
        chain_i: First chain ID.
        chain_j: Second chain ID.
        coord_i: Interface coordinate on chain i (Angstroms).
        coord_j: Interface coordinate on chain j (Angstroms).
        residues_i: Set of contacting residue IDs on chain i.
        residues_j: Set of contacting residue IDs on chain j.
        residue_details_i: Detailed residue information for chain i.
        residue_details_j: Detailed residue information for chain j.
        energy: Optional binding energy (default -1.0).
    """
    chain_i: str
    chain_j: str
    coord_i: np.ndarray
    coord_j: np.ndarray
    residues_i: Set[int]
    residues_j: Set[int]
    residue_details_i: List[ResidueInfo]
    residue_details_j: List[ResidueInfo]
    is_hht: bool = False
    hht_orientation: str = "ji"
    energy: float = -1.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert InterfaceString to JSON-serializable dictionary."""
        return {
            'chain_i': self.chain_i,
            'chain_j': self.chain_j,
            'coord_i': self.coord_i.tolist() if isinstance(self.coord_i, np.ndarray) else self.coord_i,
            'coord_j': self.coord_j.tolist() if isinstance(self.coord_j, np.ndarray) else self.coord_j,
            'residues_i': list(self.residues_i),
            'residues_j': list(self.residues_j),
            'residue_details_i': [residue.to_dict() for residue in self.residue_details_i],
            'residue_details_j': [residue.to_dict() for residue in self.residue_details_j],
            'energy': self.energy
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'InterfaceString':
        """Create InterfaceString from dictionary."""
        return cls(
            chain_i=data['chain_i'],
            chain_j=data['chain_j'],
            coord_i=np.array(data['coord_i']),
            coord_j=np.array(data['coord_j']),
            residues_i=set(data['residues_i']),
            residues_j=set(data['residues_j']),
            residue_details_i=[ResidueInfo.from_dict(r) for r in data['residue_details_i']],
            residue_details_j=[ResidueInfo.from_dict(r) for r in data['residue_details_j']],
            energy=data.get('energy', -1.0)
        )

    def get_residue_composition_i(self) -> Dict[str, int]:
        """Get amino acid composition for chain i side of interface.

        Returns:
            Dictionary mapping amino acid names to counts.
        """
        composition = {}
        for residue in self.residue_details_i:
            composition[residue.residue_name] = composition.get(residue.residue_name, 0) + 1
        return composition

    def get_residue_composition_j(self) -> Dict[str, int]:
        """Get amino acid composition for chain j side of interface.

        Returns:
            Dictionary mapping amino acid names to counts.
        """
        composition = {}
        for residue in self.residue_details_j:
            composition[residue.residue_name] = composition.get(residue.residue_name, 0) + 1
        return composition

    def get_residue_sequence_i(self) -> str:
        """Get ordered residue sequence for chain i side of interface.

        Returns:
            String of single-letter amino acid codes.
        """
        # Sort by residue ID to get proper sequence order
        sorted_residues = sorted(self.residue_details_i, key=lambda r: r.residue_id)
        return ''.join(self._three_to_one(r.residue_name) for r in sorted_residues)

    def get_residue_sequence_j(self) -> str:
        """Get ordered residue sequence for chain j side of interface.

        Returns:
            String of single-letter amino acid codes.
        """
        # Sort by residue ID to get proper sequence order
        sorted_residues = sorted(self.residue_details_j, key=lambda r: r.residue_id)
        return ''.join(self._three_to_one(r.residue_name) for r in sorted_residues)

    def _three_to_one(self, three_letter: str) -> str:
        """Convert three-letter amino acid code to one-letter code."""
        conversion = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        return conversion.get(three_letter.upper(), 'X')

    def calculate_kon(self, default_ka: float = None) -> float:
        """Calculate association rate constant.
        
        Returns fixed diffusion-limited association rate based on:
        kon = 120 nm³/μs
        
        Args:
            default_ka: Default association rate in nm³/μs. If None, fetches 
                       from PDBModelHyperparameters default.
        
        Returns:
            float: Association rate constant in nm³/μs
        """
        if default_ka is None:
            # Fetch default from hyperparameters class
            # We instantiate to get the default value defined in field()
            default_ka = PDBModelHyperparameters().default_on_rate_3d_ka
            
        return default_ka  # nm³/μs (diffusion-limited)

    def calculate_koff(self, temperature: float = 298.0, kon: float = None) -> float:
        """Calculate dissociation rate constant from binding energy and on-rate.
        
        Uses thermodynamic relationship:
        koff = (kon * C0) * exp(ΔG/RT)
        
        where:
        - kon is the association rate (default 120 nm³/μs)
        - C0 is standard concentration (~0.6022 nm⁻³)
        - ΔG is the binding free energy
        
        Args:
            temperature: Temperature in Kelvin (default: 298K)
            kon: Association rate in nm³/μs (optional, uses calculate_kon() default if None)
        
        Returns:
            float: Dissociation rate constant in s⁻¹
        """
        import math
        
        R = 0.008314  # Gas constant in kJ/(mol·K)
        
        # Get kon if not provided
        if kon is None:
            kon = self.calculate_kon()
            
        # Standard concentration C0 calculated from 1 Molar
        # 1 M = 1 mol/L = 6.022e23 particles / 1e24 nm³ = 0.602214 nm⁻³
        C0 = 0.602214076
        
        # Calculate prefactor A = kon * C0
        # kon is in nm³/μs, C0 is in nm⁻³
        # A_us = kon * C0 (in μs⁻¹)
        # A_s = A_us * 1e6 (in s⁻¹)
        # For kon=120, A_s ≈ 7.22e7 s⁻¹
        
        prefactor = kon * C0 * 1e6
        
        # Use binding energy (ΔG in kJ/mol)
        if self.energy == -1.0:
            # Use default energy if not predicted by ProAffinity
            delta_G = -16 * R * temperature  # -16RT in kJ/mol
        else:
            delta_G = self.energy  # kJ/mol from ProAffinity
        
        # Calculate koff using: koff = (kon * C0) * exp(ΔG/RT)
        
        koff = prefactor * math.exp(delta_G / (R * temperature))
        
        return koff  # s⁻¹

    def get_rates(self, temperature: float = 298.0) -> tuple:
        """Get both association and dissociation rate constants.
        
        Args:
            temperature: Temperature in Kelvin (default: 298K)
        
        Returns:
            tuple: (kon, koff) where kon is in nm³/μs and koff is in s⁻¹
        """
        return (self.calculate_kon(), self.calculate_koff(temperature))


@dataclass
class CoarseGrainedChain:
    """Coarse-grained representation of a protein chain.

    Contains the essential geometric and topological information
    needed for molecular simulation setup. Interface references are
    managed separately to avoid JSON serialization issues.

    Attributes:
        chain_id: Chain identifier.
        com: Center of mass in Angstroms.
        radius: Coarse-grained radius in Angstroms.
        sequence: Amino acid sequence.
        bbox_min: Minimum coordinates of bounding box.
        bbox_max: Maximum coordinates of bounding box.
    """
    chain_id: str
    com: np.ndarray
    radius: float
    sequence: str
    bbox_min: np.ndarray
    bbox_max: np.ndarray
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert CoarseGrainedChain to JSON-serializable dictionary."""
        return {
            'chain_id': self.chain_id,
            'com': self.com.tolist() if isinstance(self.com, np.ndarray) else self.com,
            'radius': self.radius,
            'sequence': self.sequence,
            'bbox_min': self.bbox_min.tolist() if isinstance(self.bbox_min, np.ndarray) else self.bbox_min,
            'bbox_max': self.bbox_max.tolist() if isinstance(self.bbox_max, np.ndarray) else self.bbox_max
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CoarseGrainedChain':
        """Create CoarseGrainedChain from dictionary."""
        return cls(
            chain_id=data['chain_id'],
            com=np.array(data['com']),
            radius=data['radius'],
            sequence=data['sequence'],
            bbox_min=np.array(data['bbox_min']),
            bbox_max=np.array(data['bbox_max'])
        )


class CoarseGrainer:
    """Coarse-graining engine for protein structure analysis.

    Implements the interface detection pipeline using KD-tree spatial
    queries to identify binding interfaces between protein chains.
    Generates coarse-grained representations suitable for molecular
    simulation with detailed residue composition information.

    Attributes:
        parser: PDB parser containing structure data.
        hyperparams: Configuration parameters for coarse-graining.
        chains: Dictionary of coarse-grained chain representations.
        interfaces: List of all detected interfaces.
        partner_map: Mapping of interface partnerships.
    """

    def __init__(self, parser: PDBParser, hyperparams: PDBModelHyperparameters):
        """Initialize coarse-grainer with parser and parameters.

        Args:
            parser: Initialized PDB parser with structure data.
            hyperparams: Hyperparameters for coarse-graining process.
        """
        self.parser = parser
        self.hyperparams = hyperparams
        self.chains: Dict[str, CoarseGrainedChain] = {}
        self.interfaces: List[InterfaceString] = []
        self.partner_map: Dict[Tuple[str, int], Tuple[str, int]] = {}

        # Run coarse-graining pipeline
        self._run_coarse_graining()

    def _run_coarse_graining(self) -> None:
        """Execute the complete coarse-graining pipeline."""
        # Initialize coarse-grained chains
        self._initialize_chains()

        # Detect interfaces between all chain pairs
        self._detect_all_interfaces()
        
        # Run ProAffinity batch prediction if enabled
        if self.hyperparams.predict_affinity:
            self._predict_interface_energies()

        # Build partner mapping
        self._build_partner_mapping()

    def _initialize_chains(self) -> None:
        """Initialize coarse-grained chain representations."""
        min_length = getattr(self.hyperparams, 'min_chain_length', 4)
        
        for chain_id in self.parser.get_chain_ids():
            chain_data = self.parser.get_chain_data(chain_id)
            
            # Filter out short chains (small molecules)
            sequence = chain_data.get('sequence', '')
            if len(sequence) < min_length:
                if hasattr(self, 'parser') and hasattr(self.parser, 'workspace_manager'):
                    if self.parser.workspace_manager:
                        self.parser.workspace_manager.logger.info(
                            "Skipping chain %s: only %d residues (min_chain_length=%d)",
                            chain_id, len(sequence), min_length
                        )
                continue

            self.chains[chain_id] = CoarseGrainedChain(
                chain_id=chain_id,
                com=chain_data['com'].copy(),
                radius=chain_data['radius'],
                sequence=chain_data['sequence'],
                bbox_min=chain_data['bbox_min'].copy(),
                bbox_max=chain_data['bbox_max'].copy()
            )

    def _detect_all_interfaces(self) -> None:
        """Detect interfaces between all valid chain pairs."""
        chain_ids = list(self.chains.keys())

        for i, chain_i in enumerate(chain_ids):
            for j, chain_j in enumerate(chain_ids[i+1:], i+1):
                # Check if chains can potentially interact via bounding boxes
                if self._can_chains_interact(chain_i, chain_j):
                    interfaces = self._detect_interface(chain_i, chain_j)
                    if interfaces:
                        self.interfaces.extend(interfaces)
                        # Don't store InterfaceString objects in chains to avoid serialization issues
                        # The interfaces are already stored in self.interfaces
                        # Access interfaces by chain using get_interfaces_for_chain()

    def _can_chains_interact(self, chain_i: str, chain_j: str) -> bool:
        """Check if two chains can potentially interact via bounding box test.

        Args:
            chain_i: First chain ID.
            chain_j: Second chain ID.

        Returns:
            True if chains' bounding boxes are within interaction distance.
        """
        chain_data_i = self.chains[chain_i]
        chain_data_j = self.chains[chain_j]

        # Convert distance cutoff to Angstroms
        r_cut_angstrom = self.parser.convert_distance_to_angstrom(
            self.hyperparams.interface_detect_distance_cutoff
        )

        # Check separation along each axis
        for d in range(3):  # x, y, z
            # Check if bounding boxes are separated by more than cutoff
            if (chain_data_j.bbox_min[d] - chain_data_i.bbox_max[d] > r_cut_angstrom or
                    chain_data_i.bbox_min[d] - chain_data_j.bbox_max[d] > r_cut_angstrom):
                return False

        return True

    def _detect_interface(self, chain_i: str, chain_j: str) -> List[InterfaceString]:
        """Detect interfaces between two specific chains with detailed residue information.
        Splits distinct spatial patches into separate interfaces.

        Args:
            chain_i: First chain ID.
            chain_j: Second chain ID.

        Returns:
            List of Interface objects found (empty list if none).
        """
        # Get chain data
        data_i = self.parser.get_chain_data(chain_i)
        data_j = self.parser.get_chain_data(chain_j)

        # Extract Cα coordinates and residue info
        coords_i = data_i['ca_coords']
        coords_j = data_j['ca_coords']
        residues_i = data_i['residues']
        residues_j = data_j['residues']

        if len(coords_i) == 0 or len(coords_j) == 0:
            return []

        # Convert distance cutoff to Angstroms
        r_cut_angstrom = self.parser.convert_distance_to_angstrom(
            self.hyperparams.interface_detect_distance_cutoff
        )

        # Build KD-tree for chain j
        tree_j = KDTree(coords_j)

        # Query neighbors for each residue in chain i
        neighbor_lists = tree_j.query_ball_point(coords_i, r=r_cut_angstrom)

        # Find hitting indices
        hit_i_indices = np.where([len(nbrs) > 0 for nbrs in neighbor_lists])[0]
        hit_j_indices = np.array(sorted({j_idx for nbrs in neighbor_lists for j_idx in nbrs}))
        
        n_i_hits = len(hit_i_indices)
        n_j_hits = len(hit_j_indices)

        if n_i_hits < self.hyperparams.interface_detect_n_residue_cutoff or n_j_hits < self.hyperparams.interface_detect_n_residue_cutoff:
            return []
            
        detected_interfaces = []

        # Optionally disable spatial patch splitting and treat the whole contacting
        # chain-pair surface as a single interface.
        if not getattr(self.hyperparams, "split_spatial_interface_patches", False):
            interface_coord_i = coords_i[hit_i_indices].mean(axis=0)
            interface_coord_j = coords_j[hit_j_indices].mean(axis=0)

            contacting_residues_i = {residues_i[idx]['id'] for idx in hit_i_indices}
            contacting_residues_j = {residues_j[idx]['id'] for idx in hit_j_indices}

            residue_details_i = []
            for idx in hit_i_indices:
                residue_data = residues_i[idx]
                residue_details_i.append(ResidueInfo(
                    residue_id=residue_data['id'],
                    residue_name=residue_data['name'],
                    position=coords_i[idx].copy(),
                    chain_id=chain_i
                ))

            residue_details_j = []
            for idx in hit_j_indices:
                residue_data = residues_j[idx]
                residue_details_j.append(ResidueInfo(
                    residue_id=residue_data['id'],
                    residue_name=residue_data['name'],
                    position=coords_j[idx].copy(),
                    chain_id=chain_j
                ))

            detected_interfaces.append(InterfaceString(
                chain_i=chain_i,
                chain_j=chain_j,
                coord_i=interface_coord_i,
                coord_j=interface_coord_j,
                residues_i=contacting_residues_i,
                residues_j=contacting_residues_j,
                residue_details_i=residue_details_i,
                residue_details_j=residue_details_j,
                energy=-1.0
            ))

            if self.parser.workspace_manager:
                self.parser.workspace_manager.logger.debug(
                    f"Found single merged interface {chain_i}-{chain_j} with "
                    f"{len(hit_i_indices)}/{len(hit_j_indices)} residues "
                    f"(spatial patch splitting disabled)"
                )

            return detected_interfaces

        # Otherwise keep the existing clustering-based splitting
        points_i = coords_i[hit_i_indices]
        points_j = coords_j[hit_j_indices]

        all_points = np.vstack((points_i, points_j))

        separation_threshold = max(15.0, r_cut_angstrom * 2.0)
        if len(all_points) > 1:
            z = linkage(all_points, method='single')
            labels = fcluster(z, separation_threshold, criterion='distance')
        else:
            labels = np.array([1])

        num_clusters = labels.max()

        labels_i = labels[:len(hit_i_indices)]
        labels_j = labels[len(hit_i_indices):]

        for cluster_id in range(1, num_clusters + 1):
            cluster_mask_i_local = (labels_i == cluster_id)
            cluster_mask_j_local = (labels_j == cluster_id)

            cluster_indices_i = hit_i_indices[cluster_mask_i_local]
            cluster_indices_j = hit_j_indices[cluster_mask_j_local]

            if len(cluster_indices_i) < self.hyperparams.interface_detect_n_residue_cutoff or \
            len(cluster_indices_j) < self.hyperparams.interface_detect_n_residue_cutoff:
                continue

            interface_coord_i = coords_i[cluster_indices_i].mean(axis=0)
            interface_coord_j = coords_j[cluster_indices_j].mean(axis=0)

            contacting_residues_i = {residues_i[idx]['id'] for idx in cluster_indices_i}
            contacting_residues_j = {residues_j[idx]['id'] for idx in cluster_indices_j}

            residue_details_i = []
            for idx in cluster_indices_i:
                residue_data = residues_i[idx]
                residue_details_i.append(ResidueInfo(
                    residue_id=residue_data['id'],
                    residue_name=residue_data['name'],
                    position=coords_i[idx].copy(),
                    chain_id=chain_i
                ))

            residue_details_j = []
            for idx in cluster_indices_j:
                residue_data = residues_j[idx]
                residue_details_j.append(ResidueInfo(
                    residue_id=residue_data['id'],
                    residue_name=residue_data['name'],
                    position=coords_j[idx].copy(),
                    chain_id=chain_j
                ))

            detected_interfaces.append(InterfaceString(
                chain_i=chain_i,
                chain_j=chain_j,
                coord_i=interface_coord_i,
                coord_j=interface_coord_j,
                residues_i=contacting_residues_i,
                residues_j=contacting_residues_j,
                residue_details_i=residue_details_i,
                residue_details_j=residue_details_j,
                energy=-1.0
            ))

            if self.parser.workspace_manager:
                self.parser.workspace_manager.logger.debug(
                    f"Cluster {cluster_id}: Found interface {chain_i}-{chain_j} "
                    f"with {len(cluster_indices_i)}/{len(cluster_indices_j)} residues"
                )

        return detected_interfaces

    def _build_partner_mapping(self) -> None:
        """Build explicit binding partner index mapping.

        Creates the π_intf mapping: (i,k) → (j,m) where i,k refers to
        the k-th partner on chain i and j,m refers to the m-th partner on chain j.
        """
        # Count interfaces per chain to assign indices
        chain_partner_counts: Dict[str, int] = {
            chain_id: 0 for chain_id in self.chains}

        for interface in self.interfaces:
            # Assign partner indices
            i_partner_idx = chain_partner_counts[interface.chain_i]
            j_partner_idx = chain_partner_counts[interface.chain_j]

            # Create bidirectional mapping
            self.partner_map[(interface.chain_i, i_partner_idx)] = (
                interface.chain_j, j_partner_idx)
            self.partner_map[(interface.chain_j, j_partner_idx)] = (
                interface.chain_i, i_partner_idx)

            # Increment counters
            chain_partner_counts[interface.chain_i] += 1
            chain_partner_counts[interface.chain_j] += 1

    def _predict_interface_energies(self) -> None:
        """Predict binding energies for all interfaces using ProAffinity-GNN.
        
        Runs batch prediction for all detected interfaces and updates
        their energy values. Falls back to default energy if prediction fails.
        """
        if not self.interfaces:
            return
        
        # Prepare batch prediction data
        affinity_prediction_pairs = []
        for interface in self.interfaces:
            affinity_prediction_pairs.append({
                'pdb_file': str(self.parser.filepath),
                'chains': f"{interface.chain_i},{interface.chain_j}",
                'interface': interface  # Store reference to update later
            })
        logger.warning("NOTE: Using ProAffinity-GNN without FASTA sequences for binding energy prediction. ")
        logger.info("\n" + "="*80)
        logger.info("NOTE: Using ProAffinity-GNN for binding energy prediction")
        logger.info("="*80)
        logger.info("This is an easy-to-use version that skips sequence alignment with")
        logger.info("canonical FASTA sequences. For better accuracy and advanced options,")
        logger.info("please visit: https://github.com/legendzzy/ProAffinity-GNN")
        logger.info("="*80 + "\n")
        logger.info(f"Predicting energies for {len(affinity_prediction_pairs)} interfaces...")
        logger.info("="*80 + "\n")
        
        try:
            # Import here to avoid dependency issues if not used
            from ..proaffinity_predictor import predict_proaffinity_binding_energy_batch
            
            # Run batch predictions
            binding_energies = predict_proaffinity_binding_energy_batch(
                predictions_list=affinity_prediction_pairs,
                adfr_path=self.hyperparams.adfr_path,
                verbose=True
            )
            
            # Update interface energies
            for pair_info, binding_energy in zip(affinity_prediction_pairs, binding_energies):
                interface = pair_info['interface']
                if np.isnan(binding_energy):
                    # Fall back to default
                    R = 0.008314  # kJ/(mol·K)
                    T = 298.0
                    interface.energy = -16 * R * T  # -16RT
                    logger.error(f"ProAffinity prediction failed for {interface.chain_i}-{interface.chain_j}, using default energy")
                else:
                    interface.energy = binding_energy
                    logger.info(f"Predicted energy for {interface.chain_i}-{interface.chain_j}: {binding_energy:.2f} kJ/mol")
                    
        except Exception as e:
            logger.error(f"Batch affinity prediction failed: {e}")
            logger.error("Using default energies for all interfaces")
            # Set default energies
            R = 0.008314
            T = 298.0
            for interface in self.interfaces:
                interface.energy = -16 * R * T

    def get_coarse_grained_chains(self) -> Dict[str, CoarseGrainedChain]:
        """Get all coarse-grained chain representations.

        Returns:
            Dictionary mapping chain IDs to CoarseGrainedChain objects.
        """
        return self.chains.copy()

    def get_interfaces(self) -> List[InterfaceString]:
        """Get all detected interfaces.

        Returns:
            List of Interface objects with detailed residue information.
        """
        return self.interfaces.copy()

    def get_interfaces_for_chain(self, chain_id: str) -> List[InterfaceString]:
        """Get all interfaces involving a specific chain.
        
        Args:
            chain_id: Chain ID to get interfaces for.
            
        Returns:
            List of interfaces involving the specified chain.
        """
        return [interface for interface in self.interfaces 
                if interface.chain_i == chain_id or interface.chain_j == chain_id]

    def get_partner_mapping(self) -> Dict[Tuple[str, int], Tuple[str, int]]:
        """Get the binding partner index mapping.

        Returns:
            Dictionary mapping (chain_id, partner_index) to (partner_chain_id, partner_index).
        """
        return self.partner_map.copy()

    def get_summary(self) -> Dict[str, any]:
        """Get summary statistics of coarse-graining results.

        Returns:
            Dictionary containing summary information including residue composition.
        """
        # Calculate residue composition statistics
        total_interface_residues = 0
        amino_acid_counts = {}

        for interface in self.interfaces:
            total_interface_residues += len(interface.residue_details_i) + len(interface.residue_details_j)

            # Count amino acids
            for residue in interface.residue_details_i + interface.residue_details_j:
                amino_acid_counts[residue.residue_name] = amino_acid_counts.get(residue.residue_name, 0) + 1

        return {
            "num_chains": len(self.chains),
            "num_interfaces": len(self.interfaces),
            "chains": list(self.chains.keys()),
            "interface_pairs": [(intf.chain_i, intf.chain_j) for intf in self.interfaces],
            "total_interface_residues": total_interface_residues,
            "amino_acid_composition": amino_acid_counts,
            "average_interface_size": total_interface_residues / (2 * len(self.interfaces)) if self.interfaces else 0
        }

    def print_interface_details(self) -> None:
        """Print detailed information about all detected interfaces."""
        print(f"Detected {len(self.interfaces)} interfaces:")
        print("=" * 80)

        for i, interface in enumerate(self.interfaces):
            print(f"\nInterface {i+1}: {interface.chain_i} ↔ {interface.chain_j}")
            print(f"Distance: {np.linalg.norm(interface.coord_i - interface.coord_j):.2f} Å")
            print(f"Contacts: {len(interface.residue_details_i)} ↔ {len(interface.residue_details_j)} residues")

            print(f"\nChain {interface.chain_i} residues:")
            for residue in sorted(interface.residue_details_i, key=lambda r: r.residue_id):
                print(f"  {residue}")

            print(f"\nChain {interface.chain_j} residues:")
            for residue in sorted(interface.residue_details_j, key=lambda r: r.residue_id):
                print(f"  {residue}")

            print(f"\nAmino acid composition:")
            comp_i = interface.get_residue_composition_i()
            comp_j = interface.get_residue_composition_j()
            print(f"  Chain {interface.chain_i}: {comp_i}")
            print(f"  Chain {interface.chain_j}: {comp_j}")

            print("-" * 40)
