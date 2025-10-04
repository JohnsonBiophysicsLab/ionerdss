"""
ionerdss.model.pdb.chain_grouping

@author yying7@jh.edu

Repeated chain detection and classification.

This module implements multiple strategies for detecting repeated (symmetry-related
or highly similar) protein chains, including header-based, sequence-based, and
structure-based grouping methods. This is essential for identifying biological
assemblies, symmetry relationships, and reducing computational complexity by
treating similar chains as equivalent entities.

## Key Concepts

### Chain Groups
A **chain group** represents a collection of protein chains that are considered
equivalent based on structural, sequence, or metadata similarity. Each group has:
- A **representative chain**: The first chain encountered in the group
- **Member chains**: All chains belonging to the group
- **Grouping method**: The strategy used to create the group

### Grouping Strategies
The module supports multiple strategies with automatic fallback:
1. **Header-based**: Uses mmCIF entity information (fastest, most reliable)
2. **Sequence-based**: Compares amino acid sequences using alignment scores
3. **Structure-based**: Performs 3D structural superposition and RMSD calculation

## Classes

### `ChainGroup`

Represents a group of similar chains.

```python
class ChainGroup:
    def __init__(self, representative: str, members: List[str], method: str)
```

**Attributes:**
- `representative`: Chain ID of the group representative
- `members`: Sorted list of all chain IDs in the group
- `grouping_method`: Method used to create this group ("header", "sequence",
"structure", "singleton")

**Methods:**
- `__len__()`: Returns number of chains in group
- `__contains__(chain_id)`: Checks if chain is in group

### `ChainGrouper`

Main engine for detecting and grouping chains.

```python
class ChainGrouper:
    def __init__(self, parser: PDBParser, coarse_grainer: CoarseGrainer, 
                 hyperparams: PDBModelHyperparameters)
```

**Key Methods:**
- `get_groups()`: Returns all detected chain groups
- `get_group_for_chain(chain_id)`: Finds group containing specific chain
- `get_representative(chain_id)`: Gets representative for a chain
- `get_summary()`: Returns grouping statistics

## Grouping Strategies

### 1. Header-Based Grouping (Default Primary)

Uses mmCIF entity information to group chains that belong to the same biological entity.

**Example:**
```python
# mmCIF entity information:
# Entity 1: chains A, B (same protein)
# Entity 2: chains C, D (different protein)
# Result: Group 1 = [A, B], Group 2 = [C, D]
```

### 2. Sequence-Based Grouping (Default Fallback)

Compares amino acid sequences using pairwise alignment scores.

**Algorithm:**
1. Extract sequences for all chains
2. Perform pairwise sequence alignments
3. Calculate identity score: `alignment_score / max(len(seq1), len(seq2))`
4. Group chains with identity ≥ threshold (default: 0.8)

**Example:**
```python
# Chain A: "ACGTMKL..."
# Chain B: "ACGTMKL..." (identical)
# Chain C: "TGCAMKL..." (different)
# Result: Group 1 = [A, B], Group 2 = [C]
```

### 3. Structure-Based Grouping (Explicit Mode)

Performs 3D structural superposition using Cα atoms.

**Algorithm:**
1. Extract Cα coordinates for all chains
2. Perform pairwise structural superposition
3. Calculate RMSD after optimal alignment
4. Group chains with RMSD ≤ threshold (default: 2.0 Å)

**Example:**
```python
# Chain A: RMSD = 0.5 Å vs Chain B → Similar
# Chain A: RMSD = 5.2 Å vs Chain C → Different
# Result: Group 1 = [A, B], Group 2 = [C]
```

## Usage Examples

### Basic Usage

```python
from ionerdss.model.pdb.chain_grouping import ChainGrouper
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

# Initialize with default settings
hyperparams = PDBModelHyperparameters(matching_mode="default")
grouper = ChainGrouper(parser, coarse_grainer, hyperparams)

# Get all groups
groups = grouper.get_groups()
for group in groups:
    print(f"Group {group.representative}: {group.members}")
    print(f"Method: {group.grouping_method}, Size: {len(group)}")
```

### Sequence-Based Grouping

```python
# Force sequence-based grouping with custom threshold
hyperparams = PDBModelHyperparameters(
    matching_mode="sequence",
    seq_threshold=0.9  # 90% sequence identity required
)
grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
```

### Structure-Based Grouping

```python
# Use structural similarity with tight RMSD threshold
hyperparams = PDBModelHyperparameters(
    matching_mode="structure",
    rmsd_threshold=1.5  # 1.5 Å RMSD threshold
)
grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
```

### Finding Specific Chain Information

```python
# Find which group a chain belongs to
group = grouper.get_group_for_chain("A")
if group:
    print(f"Chain A belongs to group with representative {group.representative}")
    print(f"Other members: {[m for m in group.members if m != 'A']}")

# Get representative for a chain
rep = grouper.get_representative("B")
print(f"Chain B is represented by chain {rep}")
```

### Getting Summary Statistics

```python
summary = grouper.get_summary()
print(f"Found {summary['num_groups']} groups using {summary['grouping_method']} method")

for group_info in summary['groups']:
    print(f"Group {group_info['representative']}: "
          f"{group_info['size']} chains, method: {group_info['method']}")
```

## Configuration

### Hyperparameters

Configure grouping behavior through `PDBModelHyperparameters`:

```python
hyperparams = PDBModelHyperparameters(
    matching_mode="default",      # "default", "sequence", "structure"
    seq_threshold=0.8,           # Sequence identity threshold (0.0-1.0)
    rmsd_threshold=2.0,          # RMSD threshold in Angstroms
    custom_aligner=None          # Custom sequence aligner (optional)
)
```

### Matching Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `"default"` | Header-based with sequence fallback | General use, recommended |
| `"sequence"` | Force sequence-based grouping | When header info unreliable |
| `"structure"` | Force structure-based grouping | Structural similarity focus |

### Thresholds

**Sequence Threshold (`seq_threshold`):**
- Range: 0.0 - 1.0
- Default: 0.8 (80% identity)
- Higher values = stricter grouping
- Lower values = more permissive grouping

**RMSD Threshold (`rmsd_threshold`):**
- Range: > 0.0 Angstroms
- Default: 2.0 Å
- Lower values = stricter structural similarity
- Higher values = more permissive structural similarity

## Algorithm Details

### Sequence Alignment Scoring

```python
# Identity calculation
identity = alignment_score / max(len(sequence1), len(sequence2))

# Grouping decision
if identity >= seq_threshold:
    # Add to same group
```

### Structural Superposition

```python
# For chains with same number of Cα atoms
superimposer = Superimposer()
superimposer.set_atoms(coords1, coords2)
rmsd = superimposer.rms

# Grouping decision
if rmsd <= rmsd_threshold:
    # Add to same group
```

### Short Chain Handling

For chains with < 3 Cα atoms, structure-based grouping uses distance-based comparison:

```python
# Calculate mean distance between corresponding atoms
distances = np.linalg.norm(coords1 - coords2, axis=1)
mean_distance = np.mean(distances)

# Use mean distance as similarity measure
similar = mean_distance <= rmsd_threshold
```

### Fallback Strategy

```python
if matching_mode == "default":
    success = group_by_header()
    if not success:
        group_by_sequence()  # Automatic fallback
```

## Performance Considerations

### Computational Complexity

| Method | Time Complexity | Space Complexity | Notes |
|--------|----------------|------------------|-------|
| Header-based | O(n) | O(n) | Fastest, linear scan |
| Sequence-based | O(n² × L) | O(n × L) | L = average sequence length |
| Structure-based | O(n² × m) | O(n × m) | m = average structure size |

## Troubleshooting

### Common Issues

#### 1. No Groups Created
```python
# Check if chains exist
chain_ids = parser.get_chain_ids()
print(f"Available chains: {chain_ids}")

# Check grouping results
summary = grouper.get_summary()
print(f"Groups created: {summary['num_groups']}")
```

#### 2. All Chains in Separate Groups
```python
# Lower thresholds for more permissive grouping
hyperparams.seq_threshold = 0.5    # Lower sequence threshold
hyperparams.rmsd_threshold = 5.0   # Higher RMSD threshold
```

#### 3. Alignment Errors
```python
# Check sequence quality
for chain_id in chain_ids:
    chain_data = parser.get_chain_data(chain_id)
    seq = chain_data['sequence']
    print(f"Chain {chain_id}: length={len(seq)}, sequence='{seq[:50]}...'")
```

#### 4. Superposition Failures
```python
# Check coordinate availability
for chain_id in chain_ids:
    chain_data = parser.get_chain_data(chain_id)
    coords = chain_data['ca_coords']
    print(f"Chain {chain_id}: {len(coords)} Cα atoms")
```

Chain groups are used downstream for:
- **Template creation**: One template per group
- **Interface detection**: Between group representatives
- **System building**: Instances based on groups
- **Simulation setup**: Reduced complexity through grouping

"""

from typing import Dict, List, Optional
import numpy as np
from Bio.PDB.Superimposer import Superimposer

from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser
from .coarse_graining import CoarseGrainer


class ChainGroup:
    """Represents a group of repeated/similar chains.

    Attributes:
        representative: Chain ID of the group representative (first encountered).
        members: List of all chain IDs in the group.
        grouping_method: Method used to create this group.
    """

    def __init__(self, representative: str, members: List[str], method: str):
        self.representative = representative
        self.members = sorted(members)  # Ensure deterministic ordering
        self.grouping_method = method

    def __len__(self) -> int:
        return len(self.members)

    def __contains__(self, chain_id: str) -> bool:
        return chain_id in self.members


class ChainGrouper:
    """Chain grouping engine for detecting repeated chains.

    Implements multiple strategies for grouping chains based on similarity:
    - Header-based: Uses mmCIF entity information
    - Sequence-based: Uses sequence alignment similarity
    - Structure-based: Uses structural superposition RMSD

    Attributes:
        parser: PDB parser with structure data.
        coarse_grainer: Coarse-grainer with chain information.
        hyperparams: Configuration parameters.
        groups: List of detected chain groups.
        chain_to_group: Mapping from chain ID to group representative.
    """

    def __init__(self, parser: PDBParser, coarse_grainer: CoarseGrainer,
                 hyperparams: PDBModelHyperparameters):
        """Initialize chain grouper.

        Args:
            parser: PDB parser with structure data.
            coarse_grainer: Coarse-grainer with processed chains.
            hyperparams: Hyperparameters including grouping method.
        """
        self.parser = parser
        self.coarse_grainer = coarse_grainer
        self.hyperparams = hyperparams
        self.groups: List[ChainGroup] = []
        self.chain_to_group: Dict[str, str] = {}

        # Run grouping based on specified method
        self._run_grouping()

    def _run_grouping(self) -> None:
        """Execute chain grouping based on specified method."""
        if self.hyperparams.matching_mode == "default":
            success = self._group_by_header()
            if not success:
                self._group_by_sequence()
        elif self.hyperparams.matching_mode == "sequence":
            self._group_by_sequence()
        elif self.hyperparams.matching_mode == "structure":
            self._group_by_structure()
        else:
            raise ValueError(
                f"Unknown matching_mode: {self.hyperparams.matching_mode}")

        # Ensure all chains are assigned to groups
        self._ensure_all_chains_grouped()

        # Sort groups by representative for determinism
        self.groups.sort(key=lambda g: g.representative)

    def _group_by_header(self) -> bool:
        """Group chains using mmCIF header entity information.

        Returns:
            True if successful, False if header information unavailable.
        """
        strand_ids = self.parser.get_strand_ids()
        if not strand_ids:
            return False

        # Group chains by entity ID
        entity_groups: Dict[str, List[str]] = {}
        chain_ids = self.parser.get_chain_ids()

        for entity_id, chains in strand_ids.items():
            # Only include chains that are in our processed chain list
            valid_chains = [c for c in chains if c in chain_ids]
            if len(valid_chains) > 0:
                entity_groups[entity_id] = valid_chains

        # Convert to ChainGroup objects
        for entity_id, chains in entity_groups.items():
            if chains:  # Skip empty groups
                representative = chains[0]  # First chain as representative
                group = ChainGroup(representative, chains, "header")
                self.groups.append(group)

                # Update mapping
                for chain_id in chains:
                    self.chain_to_group[chain_id] = representative

        return len(self.groups) > 0

    def _group_by_sequence(self) -> None:
        """Group chains based on sequence similarity using score-based approach."""
        # Get chain data
        chain_ids = self.parser.get_chain_ids()
        chains = [self.parser.get_chain_data(
            chain_id) for chain_id in chain_ids]

        # Extract sequences
        sequences = {}
        for chain_id in chain_ids:
            chain_data = self.parser.get_chain_data(chain_id)
            sequences[chain_id] = chain_data['sequence']

        # Use the score-based grouping logic from the sample
        aligner = self.hyperparams.custom_aligner
        chains_groups = []
        visited = set()

        for i in range(len(chains)):
            ci = chain_ids[i]
            if ci in visited:
                continue
            group = [ci]
            seq_i = sequences[ci]

            for j in range(i + 1, len(chains)):
                cj = chain_ids[j]
                if cj in visited:
                    continue
                seq_j = sequences[cj]
                if not seq_i or not seq_j:
                    continue

                try:
                    # Get alignment score directly
                    score = aligner.align(seq_i, seq_j).score
                    identity = score / max(len(seq_i), len(seq_j))
                    if identity >= self.hyperparams.seq_threshold:
                        group.append(cj)
                        visited.add(cj)
                except Exception as e:
                    # Handle any alignment errors gracefully
                    print(
                        f"Warning: Alignment failed for chains {ci} and {cj}: {e}")
                    continue

            visited.update(group)
            chains_groups.append(group)

        # Convert to ChainGroup objects
        for group_members in chains_groups:
            if group_members:  # Skip empty groups
                # First chain as representative
                representative = group_members[0]
                group = ChainGroup(representative, group_members, "sequence")
                self.groups.append(group)

                # Update mapping
                for chain_id in group_members:
                    self.chain_to_group[chain_id] = representative

    def _group_by_structure(self) -> None:
        """Group chains based on structural similarity."""
        # Get chain data
        chain_ids = self.parser.get_chain_ids()
        chains = [self.parser.get_chain_data(
            chain_id) for chain_id in chain_ids]

        # Use similar logic for structure-based grouping
        chains_groups = []
        visited = set()

        for i, ci in enumerate(chain_ids):
            if ci in visited:
                continue
            group = [ci]
            coords_i = chains[i]['ca_coords']

            for j, cj in enumerate(chain_ids):
                if cj in visited:
                    continue
                coords_j = chains[j]['ca_coords']

                # Check if structures are similar
                if self._are_structures_similar_coords(coords_i, coords_j):
                    group.append(cj)
                    visited.add(cj)

            visited.update(group)
            chains_groups.append(group)

        # Convert to ChainGroup objects
        for group_members in chains_groups:
            if group_members:  # Skip empty groups
                # First chain as representative
                representative = group_members[0]
                group = ChainGroup(representative, group_members, "structure")
                self.groups.append(group)

                # Update mapping
                for chain_id in group_members:
                    self.chain_to_group[chain_id] = representative

    def _are_structures_similar_coords(self, coords_i: np.ndarray, coords_j: np.ndarray) -> bool:
        """Check if two coordinate sets represent similar structures.

        Args:
            coords_i: Cα coordinates of first chain.
            coords_j: Cα coordinates of second chain.

        Returns:
            True if structures are similar below RMSD threshold.
        """
        try:
            # Need same number of Cα atoms for superposition
            if len(coords_i) != len(coords_j) or len(coords_i) == 0:
                return False

            # Handle very short chains
            if len(coords_i) < 3:
                # For very short chains, use simple distance comparison
                distances = np.linalg.norm(coords_i - coords_j, axis=1)
                mean_distance = np.mean(distances)
                return mean_distance <= self.hyperparams.rmsd_threshold

            # Perform structural superposition
            sup = Superimposer()
            sup.set_atoms(coords_i, coords_j)
            rmsd = sup.rms
            return rmsd <= self.hyperparams.rmsd_threshold

        except Exception as e:
            print(f"Warning: Structure comparison failed: {e}")
            return False

    def _ensure_all_chains_grouped(self) -> None:
        """Ensure all chains are assigned to groups (create singleton groups if needed)."""
        all_chain_ids = set(self.parser.get_chain_ids())
        grouped_chains = set(self.chain_to_group.keys())

        # Create singleton groups for ungrouped chains
        for chain_id in all_chain_ids - grouped_chains:
            group = ChainGroup(chain_id, [chain_id], "singleton")
            self.groups.append(group)
            self.chain_to_group[chain_id] = chain_id

    def get_groups(self) -> List[ChainGroup]:
        """Get all chain groups.

        Returns:
            List of ChainGroup objects.
        """
        return self.groups.copy()

    def get_group_for_chain(self, chain_id: str) -> Optional[ChainGroup]:
        """Get the group containing a specific chain.

        Args:
            chain_id: Chain identifier.

        Returns:
            ChainGroup containing the chain, or None if not found.
        """
        representative = self.chain_to_group.get(chain_id)
        if representative:
            for group in self.groups:
                if group.representative == representative:
                    return group
        return None

    def get_representative(self, chain_id: str) -> Optional[str]:
        """Get the group representative for a chain.

        Args:
            chain_id: Chain identifier.

        Returns:
            Representative chain ID, or None if not found.
        """
        return self.chain_to_group.get(chain_id)

    def get_summary(self) -> Dict[str, any]:
        """Get summary of grouping results.

        Returns:
            Dictionary with grouping statistics.
        """
        return {
            "num_groups": len(self.groups),
            "grouping_method": self.hyperparams.matching_mode,
            "groups": [
                {
                    "representative": group.representative,
                    "members": group.members,
                    "size": len(group),
                    "method": group.grouping_method
                }
                for group in self.groups
            ]
        }
