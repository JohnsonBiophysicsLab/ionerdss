"""
ionerdss.model.pdb.template_builder

Molecular and interface template generation with proper signature-based deduplication.

This module builds molecular templates and interface templates from grouped
chains, handles geometric signature calculation, template deduplication based
on signatures, and regularization across symmetry mates.


## Key Concepts

### Template-Based Modeling Rule

**Molecular Templates**: Reusable molecular components that capture the essential
geometric and interaction properties of protein chains, allowing multiple
instances to be created from a single template.

**Interface Templates**: Geometric and energetic descriptions of protein-protein
interaction sites that can be reused across multiple molecular instances.

**Signature-Based Deduplication**: Multiple interface types can exist between the
same molecule types if they have different geometric signatures, enabling complex
multi-site interactions.

### Template Hierarchy

```
System Level:
├── Molecule Types (Templates)
│   ├── Radius, diffusion constants
│   ├── Interface binding sites
│   └── Geometric properties
└── Interface Types (Templates)
    ├── Geometric signature
    ├── Binding energy
    ├── Local coordinates
    └── Partner relationships
```

### Coordinate System Management

**Input**: Coordinates in Angstroms (from structural data)
**Processing**: Automatic conversion to nanometers for NERDSS compatibility
**Local Coordinates**: Interface positions relative to molecule center of mass
**Absolute Coordinates**: Interface positions in global coordinate system

## Geometric Signatures

### Signature Components

The `GeometricSignature` captures the essential geometric relationship of an interface:

```python
@dataclass
class GeometricSignature:
    d_i: float      # COM-to-interface distance for chain i
    d_j: float      # COM-to-interface distance for chain j  
    theta_i: float  # Angle between COM-to-interface and COM-to-COM vectors (chain i)
    theta_j: float  # Angle between COM-to-interface and COM-to-COM vectors (chain j)
```

### Signature Calculation

**Distance Components**:
```python
d_i = ||interface_coord_i - COM_i||
d_j = ||interface_coord_j - COM_j||
```

**Angular Components**:
```python
offset_i = interface_coord_i - COM_i
com_vector_ij = COM_j - COM_i

theta_i = arccos(dot(offset_i, com_vector_ij) / (||offset_i|| * ||com_vector_ij||))
```

**Geometric Interpretation**:
- **d_i, d_j**: How far the interface is from each molecule's center
- **theta_i, theta_j**: The orientation of the interface relative to the intermolecular vector

### Signature Applications

**Similarity Detection**:
```python
signature1.is_similar_to(signature2, 
                        distance_threshold=5.0,  # Angstroms
                        angle_threshold=0.5)     # Radians (~30°)
```

**Homotypic Detection**:
```python
# Detect symmetric homodimer interfaces
is_homotypic = signature.is_homotypic(
    distance_threshold=1.0,  # Angstroms
    angle_threshold=0.2      # Radians (~11°)
)
```

## Template Generation Process

### 1. Molecular Template Creation

```python
def _build_molecule_template(self, group: ChainGroup) -> None:
```

**Process**:
1. **Representative Selection**: Use group representative as template basis
2. **Name Generation**: Create unique template names with conflict resolution
3. **Radius Calculation**: Convert from Angstroms to nanometers
4. **Diffusion Constants**: Calculate from molecular radius using Stokes-Einstein relation
5. **Metadata Storage**: Store grouping information and original chain names

**Template Properties**:
```python
molecule_template = MoleculeType(
    name="ProteinA",                    # Unique template name
    radius_nm=1.5,                      # Radius in nanometers
    diffusion_constants_calculated=True  # Auto-calculated from radius
)

# Metadata for traceability
molecule_template.signature = {
    'group_representative': 'A',
    'group_members': ['A', 'C', 'E'],
    'grouping_method': 'sequence_similarity',
    'original_chain_names': ['A', 'C', 'E']
}
```

### 2. Interface Template Creation

```python
def _build_all_interface_templates(self) -> None:
```

**Workflow**:
1. **Interface Processing**: Iterate through all detected interfaces
2. **Group Resolution**: Map chains to their representative groups
3. **Template Mapping**: Find corresponding molecule templates
4. **Signature Calculation**: Compute geometric signature for each interface
5. **Deduplication**: Check for existing similar interface types
6. **Template Creation**: Create new templates or reuse existing ones

### 3. Template Naming Strategy

**Molecule Templates**:
```python
# Preferred: Use representative chain name
"A" → "A"

# Conflict resolution: Add descriptive suffix
"A" (taken) → "A_group"

# Final fallback: Numeric suffix
"A_group" (taken) → "A_1"
```

**Interface Templates**:
```python
# Homodimer homotypic
"A" + "A" + index → "A_A_1"

# Heterotypic (bidirectional)
"A" + "B" + index → "A_B_1", "B_A_1"

# Multiple interfaces between same types
"A" + "B" + different_signature → "A_B_2", "B_A_2"
```

## Deduplication Strategy

### Signature-Based Matching

**Matching Criteria**:
```python
def _find_matching_interface_type(self, template_i: str, template_j: str,
                                  signature: GeometricSignature) -> Optional[str]:
```

**Process**:
1. **Template Pair Check**: Find interfaces between same molecule types
2. **Signature Comparison**: Use relaxed thresholds for matching
3. **Bidirectional Matching**: Handle both A→B and B→A orientations

**Matching Thresholds**:
```python
distance_threshold = 5.0   # 5 Angstroms tolerance
angle_threshold = 0.5      # ~30 degrees tolerance
```

### Multiple Interface Types

**Same Molecule Pair, Different Signatures**:
```python
# Interface 1: Close to molecule centers
signature_1 = GeometricSignature(d_i=3.0, d_j=3.0, theta_i=0.2, theta_j=0.2)
→ Creates: A_B_1, B_A_1

# Interface 2: Far from molecule centers  
signature_2 = GeometricSignature(d_i=8.0, d_j=8.0, theta_i=1.4, theta_j=1.4)
→ Creates: A_B_2, B_A_2
```

## Interface Types

### Homodimer Homotypic Interfaces

**Characteristics**:
- Same molecule type on both sides
- Symmetric geometric signature
- Single shared interface template

**Creation**:
```python
# Symmetric signature detection
is_homotypic = (abs(d_i - d_j) < threshold and 
                abs(theta_i - theta_j) < threshold)

if is_homotypic:
    # Create single shared template
    interface_template = InterfaceType(
        this_mol_type_name="A",
        partner_mol_type_name="A",
        interface_index=1,
        # ... coordinates and energy
    )
```

## Usage Examples

### Basic Template Building

```python
from ionerdss.model.pdb.template_builder import TemplateBuilder

# Build templates from processed components
builder = TemplateBuilder(
    parser=parser,
    coarse_grainer=coarse_grainer,
    chain_grouper=chain_grouper,
    hyperparams=hyperparams,
    workspace_manager=workspace_manager
)

# Get generated templates
molecule_templates = builder.get_molecule_templates()
interface_templates = builder.get_interface_templates()

print(f"Created {len(molecule_templates)} molecule templates:")
for name, template in molecule_templates.items():
    print(f"  {name}: radius={template.radius_nm:.2f} nm")

print(f"Created {len(interface_templates)} interface templates:")
for name, template in interface_templates.items():
    print(f"  {name}: {template.this_mol_type_name} ↔ {template.partner_mol_type_name}")
```

### Template Analysis

```python
# Get comprehensive summary
summary = builder.get_summary()

print("Template Building Summary:")
print(f"  Molecule Templates: {summary['num_molecule_templates']}")
print(f"  Interface Templates: {summary['num_interface_templates']}")

# Analyze interface type distribution
interface_counts = summary['interface_type_counts_by_molecule_pair']
for mol_pair, count in interface_counts.items():
    print(f"  {mol_pair[0]} ↔ {mol_pair[1]}: {count} interface types")

# Chain to template mapping
chain_mapping = summary['chain_name_mapping']
print("Chain → Template Mapping:")
for chain, template in chain_mapping.items():
    print(f"  Chain {chain} → Template {template}")
```


"""

from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass

import numpy as np
from Bio.PDB.Superimposer import Superimposer

from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.units import Units
from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser
from .coarse_graining import CoarseGrainer, InterfaceString, CoarseGrainedChain
from .chain_grouping import ChainGrouper, ChainGroup
from .file_manager import WorkspaceManager


@dataclass
class GeometricSignature:
    """Geometric signature for interface template deduplication.

    Attributes:
        d_i: COM-to-interface distance for chain i.
        d_j: COM-to-interface distance for chain j.
        theta_i: Angle between COM-to-interface vector and COM-to-COM vector for chain i.
        theta_j: Angle between COM-to-interface vector and COM-to-COM vector for chain j.
    """
    d_i: float
    d_j: float
    theta_i: float
    theta_j: float

    def normalize(self, precision: int = 6) -> Tuple[float, float, float, float]:
        """Normalize signature to avoid floating-point errors.

        Args:
            precision: Number of decimal places for rounding.

        Returns:
            Normalized signature tuple.
        """
        return (
            round(self.d_i, precision),
            round(self.d_j, precision),
            round(self.theta_i, precision),
            round(self.theta_j, precision)
        )

    def is_similar_to(self, other: 'GeometricSignature',
                      distance_threshold: float, angle_threshold: float) -> bool:
        """Check if this signature is similar to another signature.

        Args:
            other: Other geometric signature to compare with.
            distance_threshold: Distance threshold for similarity.
            angle_threshold: Angle threshold for similarity.

        Returns:
            True if signatures are similar within thresholds.
        """
        return (abs(self.d_i - other.d_i) < distance_threshold and
                abs(self.d_j - other.d_j) < distance_threshold and
                abs(self.theta_i - other.theta_i) < angle_threshold and
                abs(self.theta_j - other.theta_j) < angle_threshold)

    def is_homotypic(self, distance_threshold: float, angle_threshold: float) -> bool:
        """Check if signature represents a homotypic interaction.

        Args:
            distance_threshold: Distance threshold for homodimer detection.
            angle_threshold: Angle threshold for homotypic detection.

        Returns:
            True if signature is nearly symmetric (homotypic).
        """
        return (abs(self.d_i - self.d_j) < distance_threshold and
                abs(self.theta_i - self.theta_j) < angle_threshold)


class TemplateBuilder:
    """Builder for molecular and interface templates.

    Generates molecular templates from chain groups and creates interface
    templates with geometric signatures for proper deduplication. Multiple
    interface types can exist between the same molecule types if they have
    different geometric signatures.

    Attributes:
        parser: PDB parser with structure data.
        coarse_grainer: Coarse-grainer with interface data.
        chain_grouper: Chain grouper with group information.
        hyperparams: Configuration parameters.
        units: Unit system for the model.
        workspace_manager: Workspace manager for logging (optional).
        molecule_templates: Dictionary of molecular templates by name.
        interface_templates: Dictionary of interface templates by name.
        interface_signatures: Dictionary mapping interface names to signatures.
        group_to_template: Mapping from group representative to template name.
        used_template_names: Set of already used template names.
        interface_type_counters: Counter for interface types between molecule pairs.
    """

    def __init__(self, parser: PDBParser, coarse_grainer: CoarseGrainer,
                 chain_grouper: ChainGrouper, hyperparams: PDBModelHyperparameters,
                 units: Optional[Units] = None, workspace_manager: Optional[WorkspaceManager] = None):
        """Initialize template builder.

        Args:
            parser: PDB parser with structure data.
            coarse_grainer: Coarse-grainer with processed data.
            chain_grouper: Chain grouper with group information.
            hyperparams: Configuration parameters.
            units: Unit system (defaults to standard units).
            workspace_manager: Workspace manager for logging (optional).
        """
        self.parser = parser
        self.coarse_grainer = coarse_grainer
        self.chain_grouper = chain_grouper
        self.hyperparams = hyperparams
        self.units = units or Units()
        self.interface_to_type_mapping = {}

        # Get workspace manager from parser if not provided
        if workspace_manager is None and hasattr(parser, 'workspace_manager'):
            self.workspace_manager = parser.workspace_manager
        else:
            self.workspace_manager = workspace_manager

        # Template storage
        self.molecule_templates: Dict[str, MoleculeType] = {}
        self.interface_templates: Dict[str, InterfaceType] = {}
        self.interface_signatures: Dict[str, GeometricSignature] = {}

        # Mapping from group representative to template name
        self.group_to_template: Dict[str, str] = {}

        # Track used template names to avoid conflicts
        self.used_template_names: Set[str] = set()

        # Counter for interface types between molecule pairs
        self.interface_type_counters: Dict[Tuple[str, str], int] = {}

        # Build templates
        self._build_templates()

        # Regularize across groups
        self._regularize_templates()

        # Detect steric clashes if enabled
        if self.hyperparams.steric_clash_mode == "auto":
            self._detect_steric_clashes()

    def _generate_template_name(self, group: ChainGroup) -> str:
        """Generate a template name based on the chain group.

        Uses the representative chain name if unique, otherwise creates
        a descriptive name based on the group members.

        Args:
            group: Chain group to generate name for.

        Returns:
            Unique template name.
        """
        # Start with the representative chain name
        representative_name = group.representative

        # Check if the representative name is already used
        if representative_name not in self.used_template_names:
            self.used_template_names.add(representative_name)
            return representative_name

        # If representative name is taken, try variations
        base_name = representative_name

        # For groups with multiple members, try adding suffix
        if len(group.members) > 1:
            # Try adding "_group" suffix
            candidate = f"{base_name}_group"
            if candidate not in self.used_template_names:
                self.used_template_names.add(candidate)
                return candidate

        # If still conflicts, add numeric suffix
        counter = 1
        while True:
            candidate = f"{base_name}_{counter}"
            if candidate not in self.used_template_names:
                self.used_template_names.add(candidate)
                return candidate
            counter += 1

    def _build_templates(self) -> None:
        """Build molecular and interface templates from chain groups."""
        groups = self.chain_grouper.get_groups()

        # Sort groups by representative name for deterministic ordering
        groups_sorted = sorted(groups, key=lambda g: g.representative)

        for group in groups_sorted:
            # Create molecular template from group representative
            self._build_molecule_template(group)

        # Build interface templates - process all interfaces
        self._build_all_interface_templates()

    def _build_molecule_template(self, group: ChainGroup) -> None:
        """Build molecular template from chain group.

        Args:
            group: Chain group to build template from.
        """
        representative_id = group.representative
        chain_data = self.coarse_grainer.get_coarse_grained_chains()[
            representative_id]

        # Generate template name based on chain names
        template_name = self._generate_template_name(group)

        # Store mapping from group representative to template name
        self.group_to_template[representative_id] = template_name

        # Convert radius from Angstroms to nanometers
        radius_nm = chain_data.radius / 10.0

        # Create molecule type with chain-based name
        molecule_template = MoleculeType(
            name=template_name,
            radius_nm=radius_nm
        )

        # Set diffusion constants from radius
        molecule_template.set_diffusion_constants_from_radius()

        self.molecule_templates[template_name] = molecule_template

        # Add group information to template for reference
        molecule_template.signature = {
            'group_representative': representative_id,
            'group_members': group.members.copy(),
            'grouping_method': group.grouping_method,
            'original_chain_names': group.members.copy()
        }

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Created molecule template: %s for group %s",
                template_name, group.members
            )

    def _build_all_interface_templates(self) -> None:
        """Build interface templates for all detected interfaces."""
        interfaces = self.coarse_grainer.get_interfaces()

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Processing %d interfaces for template creation", len(
                    interfaces)
            )

        for interface in interfaces:
            # Get groups for both chains
            group_i = self.chain_grouper.get_group_for_chain(interface.chain_i)
            group_j = self.chain_grouper.get_group_for_chain(interface.chain_j)

            if not group_i or not group_j:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Missing group for interface %s <-> %s",
                        interface.chain_i, interface.chain_j
                    )
                continue

            # Get template names
            template_i = self.group_to_template.get(group_i.representative)
            template_j = self.group_to_template.get(group_j.representative)

            if not template_i or not template_j:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Missing template for interface %s <-> %s",
                        interface.chain_i, interface.chain_j
                    )
                continue

            # Calculate geometric signature
            signature = self._calculate_geometric_signature(interface)

            # Find or create interface type for this signature
            self._process_interface_with_signature(
                interface, template_i, template_j, signature)

    def _process_interface_with_signature(self, interface: InterfaceString,
                                      template_i: str, template_j: str,
                                      signature: GeometricSignature) -> None:
        """Process an interface and assign it to the appropriate interface type based on signature.

        Args:
            interface: Interface object.
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature for this interface.
        """
        # Check if we already have a similar interface type
        matching_interface_name = self._find_matching_interface_type(
            template_i, template_j, signature)

        if matching_interface_name:
            # Use existing interface type - update the interface object
            interface.interface_type = matching_interface_name

            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Interface %s <-> %s assigned to existing type: %s",
                    interface.chain_i, interface.chain_j, matching_interface_name
                )
        else:
            # Create new interface type(s)
            new_interface_names = self._create_new_interface_type(
                interface, template_i, template_j, signature)

            # Update the interface object with the appropriate new type
            if len(new_interface_names) == 1:
                # Homodimeric homotypic case - single type
                interface.interface_type = new_interface_names[0]
            else:
                # heterotypic case - need to determine which type this interface gets
                interface_type_name = self._assign_interface_to_heterotypic_type(
                    interface, new_interface_names, template_i, template_j)
                interface.interface_type = interface_type_name
    
    def _assign_interface_to_heterotypic_type(self, interface: InterfaceString, 
                                        interface_type_names: List[str],
                                        template_i: str, template_j: str) -> str:
        """Assign an interface to one of the heterotypic interface types.
        
        Args:
            interface: Interface object to assign.
            interface_type_names: List of available interface type names.
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            
        Returns:
            Selected interface type name.
        """
        if len(interface_type_names) != 2:
            return interface_type_names[0]
        
        # Get the two interface templates
        type_1_template = self.interface_templates[interface_type_names[0]]
        type_2_template = self.interface_templates[interface_type_names[1]]
        
        # Determine which type this interface should use based on the template's "this_side"
        # and the interface's chain information
        
        # Check which template corresponds to which side
        if (type_1_template.this_mol_type_name == template_i and 
            hasattr(type_1_template, 'signature') and 
            type_1_template.signature.get('this_side') == 'i'):
            return interface_type_names[0]
        elif (type_1_template.this_mol_type_name == template_j and 
            hasattr(type_1_template, 'signature') and 
            type_1_template.signature.get('this_side') == 'j'):
            return interface_type_names[0]
        else:
            return interface_type_names[1]
            
    def _store_interface_mapping(self, interface: InterfaceString, interface_type_name: str) -> None:
        """Store mapping from interface to interface type for later use in system building.

        Args:
            interface: Interface object.
            interface_type_name: Name of the interface type.
        """
        # Store mapping for system builder to use
        if not hasattr(self, 'interface_to_type_mapping'):
            self.interface_to_type_mapping = {}

        # For heterotypic interfaces, we need to store mappings for both sides
        # Create unique keys for both sides of the interface
        interface_key_i = f"{interface.chain_i}_{interface.chain_j}_{interface.coord_i[0]:.3f}_{interface.coord_i[1]:.3f}_{interface.coord_i[2]:.3f}"
        interface_key_j = f"{interface.chain_j}_{interface.chain_i}_{interface.coord_j[0]:.3f}_{interface.coord_j[1]:.3f}_{interface.coord_j[2]:.3f}"
        
        # Store the mapping
        self.interface_to_type_mapping[interface_key_i] = interface_type_name
        
        # For heterotypic interfaces, also store the reverse mapping if it's a different type
        if interface_type_name in self.interface_templates:
            interface_template = self.interface_templates[interface_type_name]
            if hasattr(interface_template, 'partner_interface_type') and interface_template.partner_interface_type:
                partner_type_name = interface_template.partner_interface_type.name if hasattr(interface_template.partner_interface_type, 'name') else str(interface_template.partner_interface_type)
                self.interface_to_type_mapping[interface_key_j] = partner_type_name
                
    def _create_new_interface_type(self, interface: InterfaceString,
                                   template_i: str, template_j: str,
                                   signature: GeometricSignature) -> List[str]:
        """Create new interface type(s) for the given interface.

        Args:
            interface: Interface object.
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature for this interface.

        Returns:
            List of created interface type names.
        """
        # Get next interface index for this template pair
        template_pair = tuple(sorted([template_i, template_j]))
        interface_index = self.interface_type_counters.get(
            template_pair, 0) + 1
        self.interface_type_counters[template_pair] = interface_index

        # Check if this is a homodimeric homotypic interaction
        is_homodimeric_homotypic = (template_i == template_j and
                        signature.is_homotypic(
                            self.hyperparams.homodimer_distance_threshold,
                            self.hyperparams.homodimer_angle_threshold
                        ))

        created_names = []

        if is_homodimeric_homotypic:
            # Create single shared interface template
            interface_name = self._create_homotypic_interface_template(
                interface, template_i, signature, interface_index
            )
            created_names.append(interface_name)
        else:
            # Create separate interface templates for each side
            interface_names = self._create_heterotypic_interface_templates(
                interface, template_i, template_j, signature, interface_index
            )
            created_names.extend(interface_names)

        return created_names

    def _find_matching_interface_type(self, template_i: str, template_j: str,
                                      signature: GeometricSignature) -> Optional[str]:
        """Find existing interface type that matches the given signature.

        Args:
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature to match.

        Returns:
            Name of matching interface type, or None if no match found.
        """
        # Look for existing interface types between these templates
        for interface_name, existing_signature in self.interface_signatures.items():
            interface_template = self.interface_templates[interface_name]

            # Check if this interface is between the same template types (in either direction)
            templates_match = (
                (interface_template.this_mol_type_name == template_i and
                 interface_template.partner_mol_type_name == template_j) or
                (interface_template.this_mol_type_name == template_j and
                 interface_template.partner_mol_type_name == template_i)
            )

            if templates_match:
                # Use more relaxed thresholds for signature matching
                distance_threshold = 5.0  # 5 Angstroms tolerance
                angle_threshold = 0.5     # ~30 degrees tolerance

                # Check if signatures are similar
                if signature.is_similar_to(existing_signature, distance_threshold, angle_threshold):
                    if self.workspace_manager:
                        self.workspace_manager.logger.debug(
                            "Found matching interface type %s for signature d_i=%.2f, d_j=%.2f, theta_i=%.3f, theta_j=%.3f",
                            interface_name, signature.d_i, signature.d_j, signature.theta_i, signature.theta_j
                        )
                    return interface_name

        return None

    def _create_homotypic_interface_template(self, interface: InterfaceString,
                                             template_name: str, signature: GeometricSignature,
                                             interface_index: int) -> str:
        """Create shared interface template for homodimeric homotypic interaction.

        Returns:
            Name of created interface template.
        """
        # Generate interface name using index
        interface_name = f"{template_name}_{template_name}_{interface_index}"

        # Convert coordinates to nanometers and calculate local coordinates
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_i]
        com_nm = self.parser.convert_coords_to_nm(chain_i_data.com)
        intf_coord_nm = self.parser.convert_coords_to_nm(interface.coord_i)
        local_coord_nm = intf_coord_nm - com_nm

        # Create interface template
        interface_template = InterfaceType(
            this_mol_type_name=template_name,
            partner_mol_type_name=template_name,
            interface_index=interface_index,
            absolute_coord=intf_coord_nm,
            local_coord=local_coord_nm,
            energy=interface.energy
        )

        # Add metadata about original chains
        interface_template.signature = {
            'original_chain_i': interface.chain_i,
            'original_chain_j': interface.chain_j,
            'interaction_type': 'homodimer',
            'contacting_residues_i': list(interface.residues_i),
            'contacting_residues_j': list(interface.residues_j),
            'geometric_signature': signature.normalize(self.hyperparams.signature_precision)
        }

        # Store template and signature
        self.interface_templates[interface_name] = interface_template
        self.interface_signatures[interface_name] = signature

        # Add to molecule template's interface map
        mol_template = self.molecule_templates[template_name]
        mol_template.interfaces_neighbors_map[interface_name] = template_name

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Created homodimer interface type: %s (d_i=%.2f, d_j=%.2f, theta_i=%.3f, theta_j=%.3f)",
                interface_name, signature.d_i, signature.d_j, signature.theta_i, signature.theta_j
            )

        return interface_name

    def _create_heterotypic_interface_templates(self, interface: InterfaceString,
                                                template_i: str, template_j: str,
                                                signature: GeometricSignature,
                                                interface_index: int) -> List[str]:
        """Create separate interface templates for heterotypicic interaction.

        Returns:
            List of created interface template names.
        """
        created_names = []

        # Create interface template for side i
        interface_name_i = f"{template_i}_{template_j}_{interface_index}"
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_i]
        com_i_nm = self.parser.convert_coords_to_nm(chain_i_data.com)
        intf_i_nm = self.parser.convert_coords_to_nm(interface.coord_i)
        local_i_nm = intf_i_nm - com_i_nm

        interface_template_i = InterfaceType(
            this_mol_type_name=template_i,
            partner_mol_type_name=template_j,
            interface_index=interface_index,
            absolute_coord=intf_i_nm,
            local_coord=local_i_nm,
            energy=interface.energy
        )

        # Create interface template for side j
        interface_name_j = f"{template_j}_{template_i}_{interface_index}"
        chain_j_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_j]
        com_j_nm = self.parser.convert_coords_to_nm(chain_j_data.com)
        intf_j_nm = self.parser.convert_coords_to_nm(interface.coord_j)
        local_j_nm = intf_j_nm - com_j_nm

        interface_template_j = InterfaceType(
            this_mol_type_name=template_j,
            partner_mol_type_name=template_i,
            interface_index=interface_index,
            absolute_coord=intf_j_nm,
            local_coord=local_j_nm,
            energy=interface.energy
        )

        # Add metadata
        interface_template_i.signature = {
            'original_chain_i': interface.chain_i,
            'original_chain_j': interface.chain_j,
            'interaction_type': 'heterotypic',
            'this_side': 'i',
            'contacting_residues': list(interface.residues_i),
            'geometric_signature': signature.normalize(self.hyperparams.signature_precision)
        }

        interface_template_j.signature = {
            'original_chain_i': interface.chain_i,
            'original_chain_j': interface.chain_j,
            'interaction_type': 'heterotypic',
            'this_side': 'j',
            'contacting_residues': list(interface.residues_j),
            'geometric_signature': signature.normalize(self.hyperparams.signature_precision)
        }

        # Set up cross-references
        interface_template_i.partner_interface_type = interface_template_j
        interface_template_j.partner_interface_type = interface_template_i
        interface_template_i.this_mol_type = self.molecule_templates[template_i]
        interface_template_i.partner_mol_type = self.molecule_templates[template_j]
        interface_template_j.this_mol_type = self.molecule_templates[template_j]
        interface_template_j.partner_mol_type = self.molecule_templates[template_i]

        # Store templates and signatures
        self.interface_templates[interface_name_i] = interface_template_i
        self.interface_templates[interface_name_j] = interface_template_j
        self.interface_signatures[interface_name_i] = signature
        self.interface_signatures[interface_name_j] = signature

        # Add to molecule templates' interface maps
        self.molecule_templates[template_i].interfaces_neighbors_map[interface_name_i] = template_j
        self.molecule_templates[template_j].interfaces_neighbors_map[interface_name_j] = template_i

        created_names = [interface_name_i, interface_name_j]

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Created heterotypic interface types: %s, %s (d_i=%.2f, d_j=%.2f, theta_i=%.3f, theta_j=%.3f)",
                interface_name_i, interface_name_j, signature.d_i, signature.d_j, signature.theta_i, signature.theta_j
            )

        return created_names

    def get_interface_type_for_interface(self, interface: InterfaceString) -> Optional[str]:
        """Get the interface type name for a specific interface.

        Args:
            interface: Interface object.

        Returns:
            Interface type name if found, None otherwise.
        """
        if not hasattr(self, 'interface_to_type_mapping'):
            return None

        # Create the same key used during mapping
        interface_key = f"{interface.chain_i}_{interface.chain_j}_{interface.coord_i[0]:.3f}_{interface.coord_i[1]:.3f}_{interface.coord_i[2]:.3f}"
        return self.interface_to_type_mapping.get(interface_key)

    def _calculate_geometric_signature(self, interface: InterfaceString) -> GeometricSignature:
        """Calculate geometric signature for an interface.

        Args:
            interface: Interface to calculate signature for.

        Returns:
            GeometricSignature object.
        """
        # Get chain data
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_i]
        chain_j_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_j]

        # COMs and interface coordinates
        com_i = chain_i_data.com
        com_j = chain_j_data.com
        intf_i = interface.coord_i
        intf_j = interface.coord_j

        # Calculate offset vectors
        offset_i = intf_i - com_i
        offset_j = intf_j - com_j
        com_vector_ij = com_j - com_i
        com_vector_ji = com_i - com_j

        # Calculate distances
        d_i = np.linalg.norm(offset_i)
        d_j = np.linalg.norm(offset_j)

        # Calculate angles (dot product normalized)
        if d_i > 0 and np.linalg.norm(com_vector_ij) > 0:
            theta_i = np.dot(offset_i, com_vector_ij) / \
                (d_i * np.linalg.norm(com_vector_ij))
            theta_i = np.arccos(np.clip(theta_i, -1.0, 1.0)
                                )  # Ensure valid range
        else:
            theta_i = 0.0

        if d_j > 0 and np.linalg.norm(com_vector_ji) > 0:
            theta_j = np.dot(offset_j, com_vector_ji) / \
                (d_j * np.linalg.norm(com_vector_ji))
            theta_j = np.arccos(np.clip(theta_j, -1.0, 1.0)
                                )  # Ensure valid range
        else:
            theta_j = 0.0

        return GeometricSignature(d_i, d_j, theta_i, theta_j)

    def _regularize_templates(self) -> None:
        """Regularize template geometry across repeated chains."""
        # For each group, propagate reference geometry to all members
        for group in self.chain_grouper.get_groups():
            if len(group.members) > 1:
                self._regularize_group(group)

    def _regularize_group(self, group: ChainGroup) -> None:
        """Regularize geometry within a chain group.

        Args:
            group: Chain group to regularize.
        """
        reference_id = group.representative
        reference_data = self.coarse_grainer.get_coarse_grained_chains()[
            reference_id]

        # Define reference normal vector (arbitrary choice: (0,0,1))
        reference_normal = np.array([0.0, 0.0, 1.0])

        # For each non-reference member, compute rigid transform and apply
        for member_id in group.members[1:]:  # Skip reference (first member)
            member_data = self.coarse_grainer.get_coarse_grained_chains()[
                member_id]

            # Compute rigid transform from reference to member
            transform = self._compute_rigid_transform(
                reference_data, member_data)

            # Apply transform to template geometry
            # (This would involve updating the actual chain data with transformed coordinates)
            # For now, we'll store the transform information
            setattr(member_data, 'transform_from_reference', transform)

    def _compute_rigid_transform(self, reference_data: CoarseGrainedChain,
                                 member_data: CoarseGrainedChain) -> np.ndarray:
        """Compute rigid transform between two chains.

        Args:
            reference_data: Reference chain data.
            member_data: Member chain data to transform.

        Returns:
            4x4 transformation matrix.
        """
        # Get Cα coordinates for both chains
        ref_coords = self.parser.get_chain_data(
            reference_data.chain_id)['ca_coords']
        mem_coords = self.parser.get_chain_data(
            member_data.chain_id)['ca_coords']

        if len(ref_coords) == len(mem_coords) and len(ref_coords) > 0:
            # Use Superimposer to compute transformation
            sup = Superimposer()
            try:
                sup.set_atoms(ref_coords, mem_coords)
                # Extract rotation matrix and translation vector
                rotation = sup.rotran[0]
                translation = sup.rotran[1]

                # Build 4x4 transformation matrix
                transform = np.eye(4)
                transform[:3, :3] = rotation
                transform[:3, 3] = translation
                return transform
            except Exception:
                # Fallback to identity transform
                return np.eye(4)

        return np.eye(4)

    def _detect_steric_clashes(self) -> None:
        """Detect steric clashes between interface templates."""
        # This is a simplified implementation
        # Full implementation would require detailed Cα clash checking

        for template_name, interface_template in self.interface_templates.items():
            # For each interface, check for potential clashes with other interfaces
            # on the same molecule type
            mol_type = interface_template.this_mol_type
            if mol_type:
                other_interfaces = [
                    name for name, intf in self.interface_templates.items()
                    if (intf.this_mol_type == mol_type and name != template_name)
                ]

                # Simple distance-based clash detection
                # (Real implementation would use detailed atomic coordinates)
                for other_name in other_interfaces:
                    other_intf = self.interface_templates[other_name]
                    distance = np.linalg.norm(
                        interface_template.local_coord - other_intf.local_coord
                    )

                    # If interfaces are very close, mark as mutually exclusive
                    if distance < 0.5:  # 0.5 nm threshold (adjustable)
                        interface_template.required_free.append(other_name)
                        other_intf.required_free.append(template_name)

    def get_molecule_templates(self) -> Dict[str, MoleculeType]:
        """Get all molecular templates.

        Returns:
            Dictionary mapping template names to MoleculeType objects.
        """
        return self.molecule_templates.copy()

    def get_interface_templates(self) -> Dict[str, InterfaceType]:
        """Get all interface templates.

        Returns:
            Dictionary mapping interface names to InterfaceType objects.
        """
        return self.interface_templates.copy()

    def get_template_name_for_group(self, group_representative: str) -> Optional[str]:
        """Get the template name for a group representative.

        Args:
            group_representative: Chain ID of group representative.

        Returns:
            Template name if found, None otherwise.
        """
        return self.group_to_template.get(group_representative)

    def get_chain_name_mapping(self) -> Dict[str, str]:
        """Get mapping from original chain names to template names.

        Returns:
            Dictionary mapping chain IDs to their template names.
        """
        chain_to_template = {}
        for group in self.chain_grouper.get_groups():
            template_name = self.group_to_template.get(group.representative)
            if template_name:
                for chain_id in group.members:
                    chain_to_template[chain_id] = template_name
        return chain_to_template

    def get_summary(self) -> Dict[str, any]:
        """Get summary of template building results.

        Returns:
            Dictionary with template statistics and naming information.
        """
        # Count interface types by molecule pair
        interface_type_counts = {}
        for interface_name, interface_template in self.interface_templates.items():
            mol_pair = tuple(sorted([interface_template.this_mol_type_name,
                                     interface_template.partner_mol_type_name]))
            interface_type_counts[mol_pair] = interface_type_counts.get(
                mol_pair, 0) + 1

        return {
            "num_molecule_templates": len(self.molecule_templates),
            "num_interface_templates": len(self.interface_templates),
            "molecule_templates": list(self.molecule_templates.keys()),
            "interface_templates": list(self.interface_templates.keys()),
            "interface_type_counts_by_molecule_pair": interface_type_counts,
            "group_to_template_mapping": self.group_to_template.copy(),
            "chain_name_mapping": self.get_chain_name_mapping(),
            "template_naming_strategy": "signature_based"
        }
