"""
ionerdss.model.pdb.system_builder

Final system assembly and validation with visualization support.

This module assembles the complete ionerdss System from processed templates
and creates molecule/interface instances for the final simulation system.

## Key Concepts

### System Assembly Pipeline

The SystemBuilder follows a structured assembly process that transforms parsed
structural data into a complete simulation-ready system:

```
PDB/mmCIF → Parser → CoarseGrainer → ChainGrouper → TemplateBuilder →
SystemBuilder → Complete System
```

### Instance Creation Strategy

**Molecule Instances**: Created from coarse-grained chains using the maximum
binding sites selection strategy to capture full connectivity potential.

**Interface Instances**: Generated bidirectionally for each detected interface,
ensuring proper partner relationships and cross-references.

**Cross-Reference Network**: Establishes comprehensive relationships between all
system components for efficient navigation and validation.

### Coordinate System Management

**Input**: Coordinates in Angstroms (from structural data)
**Processing**: Automatic conversion to nanometers for NERDSS compatibility
**Output**: System with coordinates in nanometers and proper unit handling

## System Assembly Process

### 1. Molecule Instance Creation

```python
def _create_molecule_instances(self) -> List[MoleculeInstance]:
    pass
```

**Process**:
1. **Chain Selection**: Iterate through all coarse-grained chains
2. **Group Mapping**: Map each chain to its representative group
3. **Template Resolution**: Find corresponding molecule template
4. **Coordinate Conversion**: Convert COM coordinates from Å to nm
5. **Instance Creation**: Create MoleculeInstance with proper attributes

**Example Chain Processing**:
```python
# For chain "A" with template "ProteinA"
molecule_instance = MoleculeInstance(
    name="A_ProteinA",           # Unique identifier
    molecule_type=protein_a_type, # Reference to template
    com=np.array([1.0, 2.0, 3.0]), # COM in nm
    norm=np.array([0.0, 0.0, 1.0]) # Default normal vector
)
```

### 2. Interface Instance Creation

```python
def _create_interface_instances(self) -> List[InterfaceInstance]:
    pass
```

**Bidirectional Creation Strategy**:
```python
# For interface A ↔ B, create two instances:
instance_i = InterfaceInstance(
    this_mol_name="A_ProteinA",
    partner_mol_name="B_ProteinB", 
    interface_type=interface_template,
    absolute_coord=coord_i_nm
)

instance_j = InterfaceInstance(
    this_mol_name="B_ProteinB",
    partner_mol_name="A_ProteinA",
    interface_type=partner_template,
    absolute_coord=coord_j_nm
)

# Pre-link partners
instance_i.partner_interface = instance_j
instance_j.partner_interface = instance_i
```

**Interface Type Resolution**:
1. **Primary Check**: Look for pre-assigned interface type on the interface object
2. **Fallback Lookup**: Query template builder for interface type assignment
3. **Template Retrieval**: Get interface template from template builder registry
4. **Partner Template**: Handle complementary interface types for heterotypic interactions

### 3. Cross-Reference Establishment

**Reference Network Creation**:

**Molecule ↔ Interface References**:
```python
# Set this_mol references for interface instances
for interface_instance in self.interface_instances:
    mol_instance = mol_instances_by_name[interface_instance.this_mol_name]
    interface_instance.this_mol = mol_instance
```

**Interface ↔ Interface Partner References**:
```python
# Bidirectional interface linking
interface_instance.partner_interface = partner_interface
partner_interface.partner_interface = interface_instance
```

**Molecule Interfaces Map**:
```python
# Build interfaces_neighbors_map: InterfaceInstance → partner MoleculeInstance
mol_instance.interfaces_neighbors_map[interface_instance] = partner_mol_instance
```

### 4. System Object Creation

```python
def _create_system(self) -> None:
    pass
```

**Registry Population**:
```python
# Add templates to registries
for molecule_type in molecule_templates.values():
    self.system.molecule_types.add(molecule_type)

for interface_type in interface_templates.values():
    self.system.interface_types.add(interface_type)

# Add instances to registries  
for molecule_instance in self.molecule_instances:
    self.system.molecule_instances.add(molecule_instance)

for interface_instance in self.interface_instances:
    self.system.interface_instances.add(interface_instance)

# Rebuild cross-references in system context
self.system._rebuild_cross_references()
```

## Component Integration

### Ring Regularization Integration

```python
# Optional ring structure regularization
if hasattr(self.hyperparams, 'ring_regularization_mode'):
    ring_regularizer = RingRegularizer(
        system=self.system,
        workspace_manager=self.workspace_manager,
        mode=getattr(self.hyperparams, 'ring_regularization_mode', 'off'),
        geometry=getattr(self.hyperparams, 'ring_geometry', 'cylinder')
    )
    ring_regularizer.regularize()
```

**Integration Benefits**:
- **Automatic Detection**: Ring regularization is applied automatically
if enabled in hyperparameters
- **Coordinate Correction**: Regularized coordinates are updated in place
within the system
- **Validation Integration**: Ring regularization results are included
in system validation

### Workspace Integration

```python
# Comprehensive workspace integration
workspace/
├── structures/
│   └── downloaded/          # Original PDB/mmCIF files
├── processed/
│   ├── coarse_grained/     # Coarse-graining results
│   ├── templates/          # Molecular templates
│   └── system/             # Final system data
├── visualizations/         # Generated plots and images
├── nerdss_files/          # NERDSS simulation files
└── logs/
    └── pipeline.log        # Complete processing log
```

## Cross-Reference Management

### Reference Types and Relationships

**1. Molecule Instance References**:
```python
molecule_instance.interfaces_neighbors_map = {
    interface_instance_1: partner_molecule_1,
    interface_instance_2: partner_molecule_2,
    # ... more interface → partner molecule mappings
}
```

**2. Interface Instance References**:
```python
interface_instance.this_mol = owning_molecule_instance
interface_instance.partner_interface = complementary_interface_instance
interface_instance.interface_type = interface_template
```

**3. Template References**:
```python
interface_instance.interface_type → InterfaceType template
molecule_instance.molecule_type → MoleculeType template
```

### Cross-Reference Validation

**Consistency Checks**:
```python
# Bidirectional consistency
assert interface_a.partner_interface.partner_interface == interface_a

# Molecule-interface consistency  
assert interface_instance.this_mol.name == interface_instance.this_mol_name

# Template consistency
assert interface_instance.interface_type in system.interface_types
```

## Usage Examples

### Basic System Building

```python
from ionerdss.model.pdb.system_builder import SystemBuilder

# Build system from processed components
builder = SystemBuilder(
    parser=parser,
    coarse_grainer=coarse_grainer, 
    chain_grouper=chain_grouper,
    template_builder=template_builder,
    hyperparams=hyperparams,
    workspace_path="/path/to/workspace",
    pdb_id="1ABC",
    workspace_manager=workspace_manager
)

# Get the assembled system
system = builder.get_system()

print(f"System contains:")
print(f"  Molecule types: {len(system.molecule_types)}")
print(f"  Interface types: {len(system.interface_types)}")
print(f"  Molecule instances: {len(system.molecule_instances)}")
print(f"  Interface instances: {len(system.interface_instances)}")
```

### System Validation

```python
# Validate the assembled system
validation_results = builder.validate_system()

if validation_results["errors"]:
    print("System validation errors:")
    for error in validation_results["errors"]:
        print(f"  - {error}")
else:
    print("System validation passed!")

if validation_results["warnings"]:
    print("System validation warnings:")
    for warning in validation_results["warnings"]:
        print(f"  - {warning}")
```

### System Summary and Statistics

```python
# Get comprehensive system summary
summary = builder.get_summary()

print("System Summary:")
print(f"  PDB ID: {summary.get('pdb_id', 'Unknown')}")
print(f"  Molecule Types: {summary['molecule_types']}")
print(f"  Interface Types: {summary['interface_types']}")
print(f"  Total Instances: {summary['molecule_instances']}")
print(f"  Total Interfaces: {summary['interface_instances']}")

# Hyperparameter information
hyperparams = summary['hyperparameters']
print(f"  Distance Cutoff: {hyperparams['distance_cutoff']} nm")
print(f"  Residue Cutoff: {hyperparams['residue_cutoff']}")

# Validation status
validation = summary['validation']
print(f"  Validation Errors: {len(validation['errors'])}")
print(f"  Validation Warnings: {len(validation['warnings'])}")
```

### Complete Pipeline Integration

```python
from ionerdss.model.pdb.parser import PDBParser
from ionerdss.model.pdb.coarse_graining import CoarseGrainer
from ionerdss.model.pdb.chain_grouping import ChainGrouper
from ionerdss.model.pdb.template_builder import TemplateBuilder
from ionerdss.model.pdb.system_builder import SystemBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.file_manager import WorkspaceManager

# Complete pipeline
with WorkspaceManager("/workspace", "1ABC") as workspace:
    # Configure parameters
    hyperparams = PDBModelHyperparameters(
        distance_cutoff=0.6,
        residue_cutoff=3,
        ring_regularization_mode="separate",
        ring_geometry="cylinder"
    )
    
    # Parse structure
    parser = PDBParser("1ABC", fetch_from_pdb=True, workspace_manager=workspace)
    
    # Coarse-grain
    coarse_grainer = CoarseGrainer(parser, hyperparams, workspace_manager=workspace)
    
    # Group chains
    chain_grouper = ChainGrouper(coarse_grainer, hyperparams, workspace_manager=workspace)
    
    # Build templates
    template_builder = TemplateBuilder(chain_grouper, hyperparams, workspace_manager=workspace)
    
    # Assemble system
    builder = SystemBuilder(
        parser=parser,
        coarse_grainer=coarse_grainer,
        chain_grouper=chain_grouper,
        template_builder=template_builder,
        hyperparams=hyperparams,
        workspace_path=workspace.workspace_path,
        pdb_id="1ABC",
        workspace_manager=workspace
    )
    
    # Get final system
    system = builder.get_system()
```

## Advanced Features

### Visualization Generation

```python
# Generate all visualizations
viz_outputs = builder.generate_visualizations()

print("Generated visualizations:")
for viz_type, viz_path in viz_outputs.items():
    print(f"  {viz_type}: {viz_path}")

# Typical outputs:
# structure_overview: /workspace/visualizations/structure_overview.png
# chain_grouping: /workspace/visualizations/chain_grouping.png  
# interfaces: /workspace/visualizations/interfaces.png
# templates: /workspace/visualizations/templates.png
```

### NERDSS Export Integration

```python
# Export NERDSS simulation files
nerdss_outputs = builder.export_nerdss_files(
    molecule_counts={"ProteinA": 50, "ProteinB": 25},
    box_nm=(200.0, 200.0, 200.0),
    parms_overrides={
        "nItr": 5e5,
        "timestep": 0.1,
        "onRate3Dka": 500.0
    }
)

print("Generated NERDSS files:")
for file_type, file_path in nerdss_outputs.items():
    print(f"  {file_type}: {file_path}")

# Outputs:
# ProteinA_mol: /workspace/nerdss_files/ProteinA.mol
# ProteinB_mol: /workspace/nerdss_files/ProteinB.mol
# parms: /workspace/nerdss_files/parms.inp
```

### Ring Regularization Control

```python
# Configure ring regularization in hyperparameters
hyperparams = PDBModelHyperparameters(
    ring_regularization_mode="uniform",  # "off", "separate", "uniform"
    ring_geometry="sphere",              # "cylinder", "sphere"
    min_ring_size=4                      # Minimum ring size to consider
)

# Ring regularization is applied automatically during system building
builder = SystemBuilder(..., hyperparams=hyperparams, ...)

# Check if ring regularization was applied
system = builder.get_system()
# Ring regularization results are integrated into the system coordinates
```

## Validation and Export

### System Validation Workflow

```python
# Comprehensive validation
validation = builder.validate_system()

# Check different validation categories
molecular_errors = [e for e in validation["errors"] if "molecule" in e.lower()]
interface_errors = [e for e in validation["errors"] if "interface" in e.lower()]
template_errors = [e for e in validation["errors"] if "template" in e.lower()]

print(f"Molecular validation errors: {len(molecular_errors)}")
print(f"Interface validation errors: {len(interface_errors)}")
print(f"Template validation errors: {len(template_errors)}")

# Validation includes:
# - Cross-reference consistency
# - Template completeness
# - Coordinate validity
# - Instance relationships
# - Registry integrity
```

### Export Options

**1. NERDSS Simulation Files**:
```python
# Complete NERDSS export
nerdss_files = builder.export_nerdss_files(
    molecule_counts={"ProteinA": 100},
    box_nm=(500.0, 500.0, 500.0)
)
```

**2. System Serialization**:
```python
# Get system for further processing
system = builder.get_system()

# System can be serialized, analyzed, or passed to other tools
# All cross-references and relationships are preserved
```

**3. Visualization Outputs**:
```python
# Generate publication-ready visualizations
visualizations = builder.generate_visualizations()

# Includes structure overviews, interface maps, template diagrams
```


"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import numpy as np

from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance
from ionerdss.model.components.units import Units
from .nerdss_exporter import NERDSSExporter
from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser
from .coarse_graining import CoarseGrainer
from .chain_grouping import ChainGrouper
from .template_builder import TemplateBuilder
from .file_manager import WorkspaceManager
from .visualizer import PDBVisualizer
from .ring_regularizer import RingRegularizer


class SystemBuilder:
    """Builder for complete ionerdss System objects.

    Assembles the final simulation system from processed molecular templates,
    creates instances for each chain, and populates all registries with
    proper cross-references.

    Attributes:
        parser: PDB parser with structure data.
        coarse_grainer: Coarse-grainer with processed data.
        chain_grouper: Chain grouper with group information.
        template_builder: Template builder with molecular templates.
        hyperparams: Configuration parameters.
        units: Unit system for the model.
        workspace_manager: Workspace manager for file organization.
        system: Assembled System object.
        molecule_instances: List of created molecule instances.
        interface_instances: List of created interface instances.
    """

    def __init__(self, parser: PDBParser, coarse_grainer: CoarseGrainer,
                 chain_grouper: ChainGrouper, template_builder: TemplateBuilder,
                 hyperparams: PDBModelHyperparameters, workspace_path: str,
                 pdb_id: Optional[str] = None, units: Optional[Units] = None,
                 workspace_manager: Optional[WorkspaceManager] = None):
        """Initialize system builder.

        Args:
            parser: PDB parser with structure data.
            coarse_grainer: Coarse-grainer with processed data.
            chain_grouper: Chain grouper with group information.
            template_builder: Template builder with molecular templates.
            hyperparams: Configuration parameters.
            workspace_path: Workspace directory path.
            pdb_id: PDB identifier (optional).
            units: Unit system (defaults to standard units).
            workspace_manager: Workspace manager for file organization.
        """
        self.parser = parser
        self.coarse_grainer = coarse_grainer
        self.chain_grouper = chain_grouper
        self.template_builder = template_builder
        self.hyperparams = hyperparams
        self.units = units or Units()
        self.workspace_manager = workspace_manager
        self.workspace_path = workspace_path
        self.pdb_id = pdb_id

        # Initialize instance storage
        self.molecule_instances: List[MoleculeInstance] = []
        self.interface_instances: List[InterfaceInstance] = []

        # Build complete system
        self._build_system()

    def _build_system(self) -> None:
        """Build the complete system with proper cross-references."""
        if self.workspace_manager:
            self.workspace_manager.logger.info("Building complete system...")

        # Step 1: Create molecule instances
        self.molecule_instances = self._create_molecule_instances()

        # Step 2: Create interface instances
        self.interface_instances = self._create_interface_instances()

        # Step 3: Establish cross-references between instances
        self._establish_cross_references()

        # Step 4: Create the final system
        self._create_system()

        # Step 5: Ring regularization (if enabled)
        if hasattr(self.hyperparams, 'ring_regularization_mode'):
            ring_regularizer = RingRegularizer(
                system=self.system,
                workspace_manager=self.workspace_manager,
                mode=getattr(self.hyperparams,
                             'ring_regularization_mode', 'off'),
                geometry=getattr(self.hyperparams, 'ring_geometry', 'cylinder')
            )
            ring_regularizer.regularize()

        if self.workspace_manager:
            self.workspace_manager.logger.info("System building completed")

    def _create_molecule_instances(self) -> List[MoleculeInstance]:
        """Create molecule instances from chains."""
        instances = []
        chains = self.coarse_grainer.get_coarse_grained_chains()

        if self.workspace_manager:
            self.workspace_manager.logger.debug("Creating molecule instances for %d chains",
                                                len(chains))

        for chain_id, chain_data in chains.items():
            # Get the group and template for this chain
            group = self.chain_grouper.get_group_for_chain(chain_id)
            if not group:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "No group found for chain %s", chain_id)
                continue

            template_name = self.template_builder.get_template_name_for_group(
                group.representative)
            if not template_name:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning("No template found for group %s",
                                                          group.representative)
                continue

            molecule_type = self.template_builder.molecule_templates.get(
                template_name)
            if not molecule_type:
                if self.workspace_manager:
                    self.workspace_manager.logger.error("Molecule type %s not found",
                                                        template_name)
                continue

            # Convert coordinates to nanometers
            com_nm = self.parser.convert_coords_to_nm(chain_data.com)

            # Create arbitrary normal vector (could be computed from structure)
            norm = np.array([0.0, 0.0, 1.0])

            # Create molecule instance
            molecule_instance = MoleculeInstance(
                name=f"{chain_id}_{template_name}",
                molecule_type=molecule_type,
                com=com_nm,
                norm=norm
            )

            instances.append(molecule_instance)

            if self.workspace_manager:
                self.workspace_manager.logger.debug("Created molecule instance: %s",
                                                    molecule_instance.name)

        return instances

    def _create_interface_instances(self) -> List[InterfaceInstance]:
        """Create interface instances from processed interfaces."""
        instances = []
        interfaces = self.coarse_grainer.get_interfaces()

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                f"Processing {len(interfaces)} interfaces for instance creation")

        for i, interface in enumerate(interfaces):
            # Check if the interface has been assigned a type during template building
            if hasattr(interface, 'interface_type') and interface.interface_type:
                interface_type_name = interface.interface_type
                if self.workspace_manager:
                    self.workspace_manager.logger.debug(
                        f"Interface {i} has assigned type: {interface_type_name}")
            else:
                # Fallback: ask template builder for the type
                interface_type_name = self.template_builder.get_interface_type_for_interface(
                    interface)
                if self.workspace_manager:
                    self.workspace_manager.logger.debug(
                        f"Interface {i} fallback type lookup: {interface_type_name}")

            if not interface_type_name:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "No interface type found for interface %s <-> %s, skipping",
                        interface.chain_i, interface.chain_j
                    )
                continue

            # Get the interface template
            interface_template = self.template_builder.interface_templates.get(
                interface_type_name)
            if not interface_template:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Interface template %s not found, skipping interface %s <-> %s",
                        interface_type_name, interface.chain_i, interface.chain_j
                    )
                continue

            # Get template names for the chains
            group_i = self.chain_grouper.get_group_for_chain(interface.chain_i)
            group_j = self.chain_grouper.get_group_for_chain(interface.chain_j)
            template_i = self.template_builder.group_to_template.get(
                group_i.representative) if group_i else None
            template_j = self.template_builder.group_to_template.get(
                group_j.representative) if group_j else None

            if not template_i or not template_j:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Missing template for interface %s <-> %s",
                        interface.chain_i, interface.chain_j
                    )
                continue

            # Always create interface instance for side i
            instance_i = InterfaceInstance(
                absolute_coord=self.parser.convert_coords_to_nm(
                    interface.coord_i),
                interface_type=interface_template,
                this_mol_name=f"{interface.chain_i}_{template_i}",
                partner_mol_name=f"{interface.chain_j}_{template_j}",
                interface_index=getattr(
                    interface_template, 'interface_index', 1),
                residues=list(interface.residues_i) if hasattr(
                    interface, 'residues_i') else [],
                energy=interface.energy
            )
            instances.append(instance_i)

            # Always create interface instance for side j (bidirectional)
            # Get the partner interface type
            partner_type = getattr(
                interface_template, 'partner_interface_type', None)

            if isinstance(partner_type, str):
                partner_template = self.template_builder.interface_templates.get(
                    partner_type)
            else:
                partner_template = partner_type

            # If no explicit partner template, use the same template (for homodimeric cases)
            if not partner_template:
                partner_template = interface_template

            instance_j = InterfaceInstance(
                absolute_coord=self.parser.convert_coords_to_nm(
                    interface.coord_j),
                interface_type=partner_template,
                this_mol_name=f"{interface.chain_j}_{template_j}",
                partner_mol_name=f"{interface.chain_i}_{template_i}",
                interface_index=getattr(
                    partner_template, 'interface_index', 1),
                residues=list(interface.residues_j) if hasattr(
                    interface, 'residues_j') else [],
                energy=interface.energy
            )
            instances.append(instance_j)

            # Pre-link the partner interfaces
            instance_i.partner_interface = instance_j
            instance_j.partner_interface = instance_i

            if self.workspace_manager:
                self.workspace_manager.logger.debug(
                    "Created bidirectional interface instances: %s <-> %s",
                    instance_i.get_name(), instance_j.get_name()
                )

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                f"Created {len(instances)} interface instances")

        return instances

    def _establish_cross_references(self) -> None:
        """Establish all cross-references between instances."""
        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Establishing cross-references between instances")

        # Create lookup maps
        mol_instances_by_name = {
            mol.name: mol for mol in self.molecule_instances}
        # Create lookup map for interface instances by their identifying properties
        interface_instances_by_key = {}
        for intf in self.interface_instances:
            key = (intf.this_mol_name, intf.partner_mol_name,
                   intf.interface_index)
            interface_instances_by_key[key] = intf
        print("-----------")
        print(interface_instances_by_key.keys())

        # Set this_mol references for interface instances
        for interface_instance in self.interface_instances:
            mol_instance = mol_instances_by_name.get(
                interface_instance.this_mol_name)
            if mol_instance:
                interface_instance.this_mol = mol_instance
            else:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Could not find molecule instance %s for interface",
                        interface_instance.this_mol_name
                    )

        # Establish partner_interface cross-references
        for interface_instance in self.interface_instances:
            # Look for the complementary interface instance
            partner_key = (interface_instance.partner_mol_name,
                           interface_instance.this_mol_name,
                           interface_instance.interface_index)

            partner_interface = interface_instances_by_key.get(partner_key)
            if partner_interface:
                # Set up bidirectional references
                interface_instance.partner_interface = partner_interface
                partner_interface.partner_interface = interface_instance

                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "Linked partner interfaces: %s <-> %s",
                        interface_instance.get_name(), partner_interface.get_name()
                    )
            else:
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "No partner interface found for %s (looking for key: %s)",
                        interface_instance.get_name(), partner_key
                    )

        # Build interfaces_neighbors_map for molecule instances
        # This maps InterfaceInstance objects to their partner MoleculeInstance objects
        for mol_instance in self.molecule_instances:
            # Find all interface instances belonging to this molecule
            mol_interfaces = [
                intf for intf in self.interface_instances
                if intf.this_mol_name == mol_instance.name
            ]

            # Build the interfaces_neighbors_map: InterfaceInstance -> partner MoleculeInstance
            for interface_instance in mol_interfaces:
                # Get the partner molecule instance
                partner_mol_instance = mol_instances_by_name.get(
                    interface_instance.partner_mol_name)
                if partner_mol_instance:
                    # Map this InterfaceInstance to its partner MoleculeInstance
                    mol_instance.interfaces_neighbors_map[interface_instance] = partner_mol_instance
                else:
                    if self.workspace_manager:
                        self.workspace_manager.logger.warning(
                            "Could not find partner molecule instance %s for interface %s",
                            interface_instance.partner_mol_name, interface_instance.get_name()
                        )

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Cross-references established successfully")

    def _create_system(self) -> None:
        """Create the final system object."""
        if self.workspace_manager:
            self.workspace_manager.logger.debug("Creating final system object")

        # Create system with correct constructor arguments
        self.system = System(
            workspace_path=self.workspace_path,
            pdb_id=self.pdb_id,
            units=self.units
        )

        # Add molecule types to registry
        molecule_templates = self.template_builder.get_molecule_templates()
        for template_name, molecule_type in molecule_templates.items():
            self.system.molecule_types.add(molecule_type)

        # Add interface types to registry
        interface_templates = self.template_builder.get_interface_templates()
        for interface_name, interface_type in interface_templates.items():
            self.system.interface_types.add(interface_type)

        # Add molecule instances to registry
        for molecule_instance in self.molecule_instances:
            self.system.molecule_instances.add(molecule_instance)

        # Add interface instances to registry
        for interface_instance in self.interface_instances:
            self.system.interface_instances.add(interface_instance)

        # Rebuild cross-references in the system
        self.system._rebuild_cross_references()

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "System created with %d molecule types, %d interface types, %d molecule instances, %d interface instances",
                len(self.system.molecule_types),
                len(self.system.interface_types),
                len(self.system.molecule_instances),
                len(self.system.interface_instances)
            )

    def generate_visualizations(self) -> Dict[str, any]:
        """Generate all visualizations for the built system.

        Returns:
            Dictionary mapping visualization types to output file paths.
        """
        if self.workspace_manager:
            self.workspace_manager.logger.info("Generating visualizations...")
            visualizer = PDBVisualizer(self.workspace_manager)
            viz_outputs = visualizer.visualize_all(
                self.parser, self.coarse_grainer, self.chain_grouper, self.template_builder
            )

            # Log each generated visualization
            for viz_type, viz_path in viz_outputs.items():
                self.workspace_manager.logger.info(
                    "Generated %s: %s", viz_type, viz_path)

            return viz_outputs
        else:
            return {}

    def get_system(self) -> System:
        """Get the assembled system.

        Returns:
            Complete System object ready for simulation.
        """
        return self.system

    def validate_system(self) -> Dict[str, list]:
        """Validate the assembled system.

        Returns:
            Dictionary with validation results from the system.
        """
        return self.system.validate_system()

    def get_summary(self) -> Dict[str, any]:
        """Get summary of the assembled system.

        Returns:
            Dictionary with system statistics and validation results.
        """
        summary = self.system.get_summary()
        validation = self.validate_system()

        summary.update({
            "validation": validation,
            "hyperparameters": self.hyperparams.to_dict()
        })

        return summary

    def export_nerdss_files(self, molecule_counts: Optional[Dict[str, int]] = None,
                            box_nm: Tuple[float, float, float] = (
                                100.0, 100.0, 100.0),
                            parms_overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Path]:
        """Export NERDSS simulation files.

        Args:
            molecule_counts: Number of molecules per type. Defaults to 10 each.
            box_nm: Simulation box size in nm. Default (100, 100, 100).
            parms_overrides: Additional parameters for parms.inp.

        Returns:
            Dictionary mapping file types to output paths.
        """
        exporter = NERDSSExporter(self.system, self.workspace_manager)
        return exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=box_nm,
            parms_overrides=parms_overrides
        )
