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

### Partner Mapping
**Partner mapping** creates explicit relationships between binding interfaces,
enabling the construction of molecular interaction networks for simulation setup.

## Classes and Data Structures

### `InterfaceString`

Represents a detected binding interface between two protein chains.

```python
@dataclass
class InterfaceString:
    chain_i: str              # First chain ID
    chain_j: str              # Second chain ID  
    coord_i: np.ndarray       # Interface coordinate on chain i (Å)
    coord_j: np.ndarray       # Interface coordinate on chain j (Å)
    residues_i: Set[int]      # Contacting residue IDs on chain i
    residues_j: Set[int]      # Contacting residue IDs on chain j
    energy: float = -1.0      # Binding energy (default -1.0)
```

### `CoarseGrainedChain`

Coarse-grained representation of a protein chain containing essential geometric
and topological information.

```python
@dataclass  
class CoarseGrainedChain:
    chain_id: str                      # Chain identifier
    com: np.ndarray                    # Center of mass (Å)
    radius: float                      # Effective radius (Å)
    interfaces: List[InterfaceString]  # Associated interfaces
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

## Algorithm Overview

### 1. Chain Initialization
```python
# Extract essential properties from parsed structure
for chain_id in parser.get_chain_ids():
    chain_data = parser.get_chain_data(chain_id)
    coarse_chain = CoarseGrainedChain(
        chain_id=chain_id,
        com=chain_data['com'],
        radius=chain_data['radius'],
        # ... other properties
    )
```

### 2. Interface Detection Pipeline
```python
# Detect interfaces between all chain pairs
for chain_i, chain_j in all_chain_pairs:
    if can_chains_interact(chain_i, chain_j):      # Bounding box pre-filter
        interface = detect_interface(chain_i, chain_j)  # Detailed detection
        if interface:
            store_interface(interface)
```

### 3. Partner Mapping Construction
```python
# Build explicit partner relationships
partner_map[(chain_i, interface_idx)] = (chain_j, partner_idx)
```

## Interface Detection Pipeline

### Step 1: Bounding Box Pre-filtering

**Purpose:** Eliminate distant chain pairs to reduce computational cost.

**Algorithm:**
```python
def can_chains_interact(chain_i, chain_j):
    r_cut = distance_cutoff_in_angstrom
    
    # Check separation along each axis
    for dimension in [x, y, z]:
        gap = min(chain_j.bbox_min[dim] - chain_i.bbox_max[dim],
                  chain_i.bbox_min[dim] - chain_j.bbox_max[dim])
        if gap > r_cut:
            return False  # Chains too far apart
    return True
```

### Step 2: Detailed Interface Detection

**Purpose:** Identify actual binding interfaces using atomic-level proximity.

**Algorithm:**
```python
def detect_interface(chain_i, chain_j):
    # Build KD-tree for efficient spatial queries
    tree_j = KDTree(chain_j.ca_coords)
    
    # Find neighbors for each residue in chain_i
    neighbor_lists = tree_j.query_ball_point(
        chain_i.ca_coords, 
        r=distance_cutoff
    )
    
    # Count contacting residues
    contacting_i = [idx for idx, neighbors in enumerate(neighbor_lists) 
                   if len(neighbors) > 0]
    contacting_j = set().union(*neighbor_lists)
    
    # Apply residue cutoff
    if len(contacting_i) >= residue_cutoff and len(contacting_j) >= residue_cutoff:
        return create_interface(contacting_i, contacting_j)
    
    return None
```

### Step 3: Interface Coordinate Calculation

**Purpose:** Determine representative coordinates for each interface.

**Algorithm:**
```python
# Calculate interface coordinates as centroids of contacting residues
interface_coord_i = mean(ca_coords_i[contacting_residues_i])
interface_coord_j = mean(ca_coords_j[contacting_residues_j])
```

### Step 4: Partner Index Mapping

**Purpose:** Create explicit partnership relationships for simulation setup.

**Algorithm:**
```python
# Assign sequential indices to interfaces on each chain
chain_interface_counts = defaultdict(int)

for interface in all_interfaces:
    i_idx = chain_interface_counts[interface.chain_i]
    j_idx = chain_interface_counts[interface.chain_j]
    
    # Create bidirectional mapping
    partner_map[(interface.chain_i, i_idx)] = (interface.chain_j, j_idx)
    partner_map[(interface.chain_j, j_idx)] = (interface.chain_i, i_idx)
    
    chain_interface_counts[interface.chain_i] += 1
    chain_interface_counts[interface.chain_j] += 1
```

## Usage Examples

### Basic Usage

```python
from ionerdss.model.pdb.coarse_graining import CoarseGrainer
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.parser import PDBParser

# Initialize components
parser = PDBParser("structure.pdb")
hyperparams = PDBModelHyperparameters(
    distance_cutoff=0.8,  # 8 Å cutoff in nm
    residue_cutoff=5      # Minimum 5 contacting residues
)

# Run coarse-graining
grainer = CoarseGrainer(parser, hyperparams)

# Access results
chains = grainer.get_coarse_grained_chains()
interfaces = grainer.get_interfaces()
partner_map = grainer.get_partner_mapping()
```

### Analyzing Results

```python
# Get summary statistics
summary = grainer.get_summary()
print(f"Processed {summary['num_chains']} chains")
print(f"Detected {summary['num_interfaces']} interfaces")
print(f"Interface pairs: {summary['interface_pairs']}")

# Examine specific chains
for chain_id, chain in chains.items():
    print(f"Chain {chain_id}:")
    print(f"  COM: {chain.com}")
    print(f"  Radius: {chain.radius:.2f} Å")
    print(f"  Interfaces: {len(chain.interfaces)}")
    print(f"  Sequence length: {len(chain.sequence)}")
```

### Interface Analysis

```python
# Analyze detected interfaces
for i, interface in enumerate(interfaces):
    print(f"Interface {i+1}: {interface.chain_i} ↔ {interface.chain_j}")
    print(f"  Distance: {np.linalg.norm(interface.coord_i - interface.coord_j):.2f} Å")
    print(f"  Contacts: {len(interface.residues_i)} ↔ {len(interface.residues_j)}")
    print(f"  Energy: {interface.energy}")
```

### Partner Mapping Usage

```python
# Explore partner relationships
for (chain, idx), (partner_chain, partner_idx) in partner_map.items():
    print(f"Chain {chain} interface {idx} ↔ Chain {partner_chain} interface {partner_idx}")

# Find partners for specific chain
chain_a_partners = [(partner_chain, partner_idx) 
                   for (chain, idx), (partner_chain, partner_idx) in partner_map.items()
                   if chain == "A"]
print(f"Chain A has {len(chain_a_partners)} binding partners")
```

## Configuration Parameters

### Distance Cutoff (`distance_cutoff`)

**Definition:** Maximum distance between Cα atoms to consider them in contact.

**Units:** Nanometers (converted to Angstroms internally)

**Typical Values:**
- `0.5 nm` (5 Å): Very tight contacts only
- `0.8 nm` (8 Å): Standard protein-protein interfaces  
- `1.2 nm` (12 Å): Loose contacts, includes water-mediated

### Residue Cutoff (`residue_cutoff`)

**Definition:** Minimum number of contacting residues required on each chain to form a valid interface.

**Units:** Number of residues

**Typical Values:**
- `3-5`: Small interfaces, individual contacts
- `5-10`: Standard protein-protein interfaces
- `10+`: Large, extensive interfaces only

## Performance Considerations

### Computational Complexity

| Operation | Time Complexity | Space Complexity | Notes |
|-----------|----------------|------------------|-------|
| Bounding box check | O(n²) | O(n) | n = number of chains |
| Interface detection | O(n² × m × log m) | O(n × m) | m = average residues per chain |
| Partner mapping | O(k) | O(k) | k = number of interfaces |
| **Total** | **O(n² × m × log m)** | **O(n × m)** | Dominated by KD-tree queries |

### Optimization Strategies

#### 1. Bounding Box Pre-filtering
```python
# Eliminates ~90% of chain pairs in typical structures
if not can_chains_interact(chain_i, chain_j):
    continue  # Skip expensive interface detection
```

#### 2. KD-tree Spatial Indexing
```python
# O(log m) neighbor queries instead of O(m) brute force
tree = KDTree(coordinates)
neighbors = tree.query_ball_point(query_points, radius)
```

#### 3. Early Termination
```python
# Stop as soon as residue cutoff is violated
if len(contacting_residues) < residue_cutoff:
    return None  # No need to continue processing
```

### Memory Usage

**Typical Memory Requirements:**
- Small protein (100 residues): ~1 MB
- Medium protein (500 residues): ~25 MB  
- Large complex (2000 residues): ~400 MB
- Viral capsid (10000+ residues): ~10+ GB

### Downstream Usage

**Chain Grouping:**
```python
# Uses coarse-grained chains for structural comparison
chains = coarse_grainer.get_coarse_grained_chains()
for chain_id, chain in chains.items():
    # Compare COM positions, radii, interface patterns
    group_chains_by_similarity(chain)
```

**Template Building:**
```python
# Uses interfaces for template creation
interfaces = coarse_grainer.get_interfaces()
for interface in interfaces:
    # Create interface templates from geometric data
    template = create_interface_template(interface)
```

**System Building:**
```python
# Uses partner mapping for instance creation
partner_map = coarse_grainer.get_partner_mapping()
for (chain, idx), (partner, p_idx) in partner_map.items():
    # Create molecular instances with proper connectivity
    create_molecule_instance(chain, partner, connectivity_info)
```

"""

from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass

import numpy as np
from scipy.spatial import KDTree

from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser


@dataclass
class InterfaceString:
    """Detected interface between two chains.

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
        energy: Optional binding energy (default -1.0).
    """
    chain_i: str
    chain_j: str
    coord_i: np.ndarray
    coord_j: np.ndarray
    residues_i: Set[int]
    residues_j: Set[int]
    energy: float = -1.0


@dataclass
class CoarseGrainedChain:
    """Coarse-grained representation of a protein chain.

    Contains the essential geometric and topological information
    needed for molecular simulation setup.

    Attributes:
        chain_id: Chain identifier.
        com: Center of mass in Angstroms.
        radius: Coarse-grained radius in Angstroms.
        interfaces: List of interfaces involving this chain.
        sequence: Amino acid sequence.
        bbox_min: Minimum coordinates of bounding box.
        bbox_max: Maximum coordinates of bounding box.
    """
    chain_id: str
    com: np.ndarray
    radius: float
    interfaces: List[InterfaceString]
    sequence: str
    bbox_min: np.ndarray
    bbox_max: np.ndarray


class CoarseGrainer:
    """Coarse-graining engine for protein structure analysis.

    Implements the interface detection pipeline using KD-tree spatial
    queries to identify binding interfaces between protein chains.
    Generates coarse-grained representations suitable for molecular
    simulation.

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

        # Build partner mapping
        self._build_partner_mapping()

    def _initialize_chains(self) -> None:
        """Initialize coarse-grained chain representations."""
        for chain_id in self.parser.get_chain_ids():
            chain_data = self.parser.get_chain_data(chain_id)

            self.chains[chain_id] = CoarseGrainedChain(
                chain_id=chain_id,
                com=chain_data['com'].copy(),
                radius=chain_data['radius'],
                interfaces=[],
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
                    interface = self._detect_interface(chain_i, chain_j)
                    if interface:
                        self.interfaces.append(interface)
                        self.chains[chain_i].interfaces.append(interface)
                        self.chains[chain_j].interfaces.append(interface)

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
            self.hyperparams.distance_cutoff
        )

        # Check separation along each axis
        for d in range(3):  # x, y, z
            # Check if bounding boxes are separated by more than cutoff
            if (chain_data_j.bbox_min[d] - chain_data_i.bbox_max[d] > r_cut_angstrom or
                    chain_data_i.bbox_min[d] - chain_data_j.bbox_max[d] > r_cut_angstrom):
                return False

        return True

    def _detect_interface(self, chain_i: str, chain_j: str) -> Optional[InterfaceString]:
        """Detect interface between two specific chains.

        Args:
            chain_i: First chain ID.
            chain_j: Second chain ID.

        Returns:
            Interface object if valid interface found, None otherwise.
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
            return None

        # Convert distance cutoff to Angstroms
        r_cut_angstrom = self.parser.convert_distance_to_angstrom(
            self.hyperparams.distance_cutoff
        )

        # Build KD-tree for chain j
        tree_j = KDTree(coords_j)

        # Query neighbors for each residue in chain i
        neighbor_lists = tree_j.query_ball_point(coords_i, r=r_cut_angstrom)

        # Find residues in chain i that have neighbors in chain j
        hit_i_mask = np.array([len(nbrs) > 0 for nbrs in neighbor_lists])
        n_i_hits = int(hit_i_mask.sum())

        # Find residues in chain j that were matched
        hit_j_indices = {j_idx for nbrs in neighbor_lists for j_idx in nbrs}
        n_j_hits = len(hit_j_indices)

        # Check if both chains meet the residue cutoff
        if n_i_hits >= self.hyperparams.residue_cutoff and n_j_hits >= self.hyperparams.residue_cutoff:
            # Calculate interface coordinates as mean of contacting Cα positions
            interface_coord_i = coords_i[hit_i_mask].mean(axis=0)
            interface_coord_j = coords_j[list(hit_j_indices)].mean(axis=0)

            # Get contacting residue IDs
            contacting_residues_i = {residues_i[idx]['id'] for idx in np.where(hit_i_mask)[
                0]}
            contacting_residues_j = {
                residues_j[idx]['id'] for idx in hit_j_indices}

            return InterfaceString(
                chain_i=chain_i,
                chain_j=chain_j,
                coord_i=interface_coord_i,
                coord_j=interface_coord_j,
                residues_i=contacting_residues_i,
                residues_j=contacting_residues_j,
                energy=-1.0  # Default energy
            )

        return None

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

    def get_coarse_grained_chains(self) -> Dict[str, CoarseGrainedChain]:
        """Get all coarse-grained chain representations.

        Returns:
            Dictionary mapping chain IDs to CoarseGrainedChain objects.
        """
        return self.chains.copy()

    def get_interfaces(self) -> List[InterfaceString]:
        """Get all detected interfaces.

        Returns:
            List of Interface objects.
        """
        return self.interfaces.copy()

    def get_partner_mapping(self) -> Dict[Tuple[str, int], Tuple[str, int]]:
        """Get the binding partner index mapping.

        Returns:
            Dictionary mapping (chain_id, partner_index) to (partner_chain_id, partner_index).
        """
        return self.partner_map.copy()

    def get_summary(self) -> Dict[str, any]:
        """Get summary statistics of coarse-graining results.

        Returns:
            Dictionary containing summary information.
        """
        return {
            "num_chains": len(self.chains),
            "num_interfaces": len(self.interfaces),
            "chains": list(self.chains.keys()),
            "interface_pairs": [(intf.chain_i, intf.chain_j) for intf in self.interfaces]
        }
