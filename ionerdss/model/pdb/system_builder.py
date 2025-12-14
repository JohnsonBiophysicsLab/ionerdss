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

from .template_builder import _enforce_identical_local_geometry_after_com

def _kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """
    Calculate best rotation R that aligns P to Q (Q ~ R @ P).
    P, Q: (N, 3) arrays.
    """
    # Center coordinates
    P_cent = P - np.mean(P, axis=0)
    Q_cent = Q - np.mean(Q, axis=0)
    
    # Compute covariance matrix
    H = np.dot(P_cent.T, Q_cent)
    
    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)
    
    # Calculate Rotation Matrix
    R = np.dot(Vt.T, U.T)
    
    # Special reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)
        
    return R

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

            # Calculate rotation and transform ref1/ref2
            ref1 = molecule_type.ref1_local
            ref2 = molecule_type.ref2_local
            
            # If this chain is not the representative, we need to find the rotation
            # that maps the representative to this chain instance
            if chain_id != group.representative:
                try:
                    coords_rep = self.parser.get_chain_data(group.representative)['ca_coords']
                    coords_curr = self.parser.get_chain_data(chain_id)['ca_coords']
                    
                    if len(coords_rep) == len(coords_curr) and len(coords_rep) > 2:
                        # Use Kabsch algorithm to find rotation
                        # P = coords_rep (template/source)
                        # Q = coords_curr (instance/target)
                        rot = _kabsch_rotation(coords_rep, coords_curr)
                        
                        # Apply rotation to reference vectors
                        ref1 = rot @ ref1
                        ref2 = rot @ ref2
                    else:
                        if self.workspace_manager:
                            self.workspace_manager.logger.warning(
                                "Cannot align chain %s to representative %s for orientation (different lengths or too short). Using default orientation.",
                                chain_id, group.representative
                            )
                except Exception as e:
                    if self.workspace_manager:
                        self.workspace_manager.logger.warning(
                            "Failed to calculate orientation for chain %s: %s. Using default orientation.",
                            chain_id, str(e)
                        )

            # Create arbitrary normal vector (could be computed from structure)
            norm = np.array([0.0, 0.0, 1.0])

            # Create molecule instance
            molecule_instance = MoleculeInstance(
                name=f"{chain_id}_{template_name}",
                molecule_type=molecule_type,
                com=com_nm,
                norm=norm,
                ref1=ref1,
                ref2=ref2
            )

            instances.append(molecule_instance)

            if self.workspace_manager:
                self.workspace_manager.logger.debug("Created molecule instance: %s",
                                                    molecule_instance.name)

        return instances

    def _create_interface_instances(self) -> List[InterfaceInstance]:
        """
        Create interface instances with:
        - de-duplication of undirected edges
        - robust partner template resolution
        - deterministic f/b assignment for homodimeric-heterotypic pairs:
            representative chain side -> ..._1f, non-representative side -> ..._1b
        """
        instances: List[InterfaceInstance] = []
        interfaces = self.coarse_grainer.get_interfaces()

        log = self.workspace_manager.logger if self.workspace_manager else None
        if log:
            log.info("Processing %d interfaces for instance creation", len(interfaces))
            log.info("Beginning interface instance creation with de-duplication")

        # ---------- helpers ----------
        def _nm(coord):
            return self.parser.convert_coords_to_nm(coord)

        def _round_sig(v: np.ndarray, nd=3) -> tuple:
            return tuple(np.round(v, nd))

        def _iface_name(obj) -> Optional[str]:
            if obj is None:
                return None
            if hasattr(obj, "get_name"):
                return obj.get_name()
            if isinstance(obj, str):
                return obj
            return None

        def _lookup_iface(name: Optional[str]):
            return self.template_builder.interface_templates.get(name) if name else None

        def _parse_iface_name(name: str):
            # A_B_1, A_A_1f, A_A_1b
            try:
                p = name.split("_")
                if len(p) < 3:
                    return (None, None, None, None)
                fam_i, fam_j = p[0], p[1]
                last = p[2]
                if len(last) >= 2 and last[-1] in ("f", "b") and last[:-1].isdigit():
                    return (fam_i, fam_j, last[:-1], last[-1])
                if last.isdigit():
                    return (fam_i, fam_j, last, None)
            except Exception:
                pass
            return (None, None, None, None)

        def _compose_name(f1, f2, idx, suf):
            return f"{f1}_{f2}_{idx}{suf}" if suf in ("f", "b") else f"{f1}_{f2}_{idx}"

        def _infer_partner_name(primary: Optional[str]) -> Optional[str]:
            if not primary:
                return None
            fi, fj, idx, suf = _parse_iface_name(primary)
            if fi is None:
                return None
            if fi != fj:
                # heterodimeric: swap families, keep idx/suf
                return _compose_name(fj, fi, idx, suf)
            # homodimeric
            if suf in ("f", "b"):  # heterotypic
                return _compose_name(fi, fj, idx, "b" if suf == "f" else "f")
            # homotypic
            return _compose_name(fi, fj, idx, None)

        def _resolve_partner_template(primary_template):
            # 1) explicit partner
            partner = getattr(primary_template, "partner_interface_type", None)
            if partner is not None:
                if hasattr(partner, "get_name"):
                    return partner
                cand = _lookup_iface(_iface_name(partner))
                if cand is not None:
                    return cand
            # 2) name inference
            primary_name = _iface_name(primary_template)
            inferred = _infer_partner_name(primary_name)
            cand = _lookup_iface(inferred)
            if cand is not None:
                return cand
            # 3) homodimeric-homotypic safe fallback to same
            fi, fj, idx, suf = _parse_iface_name(primary_name or "")
            if fi is not None and fi == fj and suf is None:
                return primary_template
            # 4) last resort: same
            if log:
                log.warning("Could not find partner template for %s; using same template", primary_name)
            return primary_template

        def _lookup_hht_meta_from_catalog(iface_template):
            """
            Try to pull canonical HHT info from template_builder.hht_catalog.

            Key shape in template_builder:
              (template_name, ordered_signature_tuple) -> {
                  'canon_order': 'ij' | 'ji',
                  'f': name_f,
                  'b': name_b,
                  'index': int,
              }

            We assume the interface template exposes:
              - get_name()
              - signature (a dict) with some ordered tuple we can use
            """
            if not hasattr(self.template_builder, "hht_catalog"):
                return None

            # get template name
            tname = None
            if hasattr(iface_template, "get_name"):
                tname = iface_template.get_name()
            else:
                tname = getattr(iface_template, "name", None)

            if not tname:
                return None

            sig = getattr(iface_template, "signature", None)
            if not sig:
                return None

            # try to find an ordered signature tuple in the signature dict
            # adjust the key name here to *your* actual signature layout
            ordered = (
                sig.get("ordered_signature")
                or sig.get("ordered_signature_tuple")
                or sig.get("ordered")
            )

            # last resort: if signature itself is already a tuple (rare)
            if ordered is None and isinstance(sig, (tuple, list)):
                ordered = tuple(sig)

            if ordered is None:
                return None

            key = (tname, tuple(ordered))
            return self.template_builder.hht_catalog.get(key)

        def _get_hht_pair(primary_template):
            """
            If this is a homodimeric-heterotypic family, return (f_template, b_template, fam, idx).
            Otherwise return (None, None, None, None).
            """
            name = _iface_name(primary_template)
            fi, fj, idx, suf = _parse_iface_name(name or "")
            if fi is None or fi != fj:
                return (None, None, None, None)
            # Needs f/b suffix on at least one side
            if suf not in ("f", "b"):
                return (None, None, None, None)

            fam = fi
            # Find both f and b templates via names; fall back to partner link if needed
            f_name = _compose_name(fam, fam, idx, "f")
            b_name = _compose_name(fam, fam, idx, "b")
            f_t = _lookup_iface(f_name)
            b_t = _lookup_iface(b_name)
            # If one is missing, try from explicit partner connection
            if f_t is None or b_t is None:
                partner_t = _resolve_partner_template(primary_template)
                # ensure both
                for nm in (f_name, b_name):
                    if _lookup_iface(nm) is None and _iface_name(partner_t) == nm:
                        if nm.endswith("f"):
                            f_t = partner_t
                        else:
                            b_t = partner_t
            if f_t is None or b_t is None:
                return (None, None, None, None)
            return (f_t, b_t, fam, idx)

        def _enforce_hht_orientation(primary_template, chain_i: str, chain_j: str):
            """
            Deterministically choose which side is 'f' and which is 'b' for
            homodimeric-heterotypic (HHT) interfaces.

            Priority:
            1. If template_builder.hht_catalog has an entry for this interface
               template + ordered signature, obey it.
               - 'canon_order' == 'ij'  -> chain_i gets 'f', chain_j gets 'b'
               - 'canon_order' == 'ji'  -> chain_j gets 'f', chain_i gets 'b'
            2. Else, fall back to your current representative-based logic.
            """
            # detect HHT pair the old way
            f_t, b_t, fam, idx = _get_hht_pair(primary_template)
            if f_t is None:
                # not an HHT case -> old path
                partner = _resolve_partner_template(primary_template)
                return (primary_template, partner)

            # 1) try catalog
            hht_meta = _lookup_hht_meta_from_catalog(primary_template)
            if hht_meta is not None:
                canon_order = hht_meta.get("canon_order")  # 'ij' or 'ji'
                name_f = hht_meta.get("f")
                name_b = hht_meta.get("b")

                # resolve to actual templates, because catalog stores names
                tmpl_f = _lookup_iface(name_f) or f_t
                tmpl_b = _lookup_iface(name_b) or b_t

                if canon_order == "ij":
                    # original order i -> f, j -> b
                    return (tmpl_f, tmpl_b)
                elif canon_order == "ji":
                    # original order j -> f, i -> b
                    return (tmpl_b, tmpl_f)
                else:
                    # unknown string: just fall back
                    if log:
                        log.warning(
                            "HHT catalog entry for %s has unknown canon_order=%s; falling back",
                            primary_template.get_name() if hasattr(primary_template, "get_name") else str(primary_template),
                            canon_order,
                        )

            # 2) catalog not available or incomplete -> use your representative rule

            g_i = self.chain_grouper.get_group_for_chain(chain_i)
            g_j = self.chain_grouper.get_group_for_chain(chain_j)
            rep = None
            if g_i and g_i.representative:
                rep = g_i.representative
            elif g_j and g_j.representative:
                rep = g_j.representative

            if rep is None:
                partner = _resolve_partner_template(primary_template)
                if log:
                    log.warning(
                        "HHT naming detected but representative not found; using default assignment"
                    )
                return (primary_template, partner)

            # rep side -> f, other side -> b
            if chain_i == rep and chain_j != rep:
                return (f_t, b_t)
            if chain_j == rep and chain_i != rep:
                return (b_t, f_t)

            # tie / pathological -> alphabetical to stay deterministic
            if chain_i <= chain_j:
                return (f_t, b_t)
            else:
                return (b_t, f_t)


        # de-dup key over undirected edges using rounded nm coords
        seen_edges = set()

        for i, interface in enumerate(interfaces):
            # Resolve interface type name
            if hasattr(interface, "interface_type") and interface.interface_type:
                iface_type_name = interface.interface_type
                if log:
                    log.info("Interface %d has assigned type: %s", i, iface_type_name)
            else:
                iface_type_name = self.template_builder.get_interface_type_for_interface(interface)
                if log:
                    log.info("Interface %d fallback type lookup: %s", i, iface_type_name)

            if not iface_type_name:
                if log:
                    log.warning("No interface type for %s <-> %s; skipping",
                                getattr(interface, "chain_i", "?"), getattr(interface, "chain_j", "?"))
                continue

            iface_template = _lookup_iface(iface_type_name)
            if not iface_template:
                if log:
                    log.warning("Interface template %s not found; skipping %s <-> %s",
                                iface_type_name, interface.chain_i, interface.chain_j)
                continue

            # Resolve molecule templates for the two chains
            gi = self.chain_grouper.get_group_for_chain(interface.chain_i)
            gj = self.chain_grouper.get_group_for_chain(interface.chain_j)
            ti = self.template_builder.group_to_template.get(gi.representative) if gi else None
            tj = self.template_builder.group_to_template.get(gj.representative) if gj else None
            if not ti or not tj:
                if log:
                    log.warning("Missing template for interface %s <-> %s", interface.chain_i, interface.chain_j)
                continue

            # Prepare de-dup edge key
            ci_nm = _nm(interface.coord_i)
            cj_nm = _nm(interface.coord_j)
            ci_sig = _round_sig(ci_nm, 3)
            cj_sig = _round_sig(cj_nm, 3)
            a, b = sorted([interface.chain_i, interface.chain_j])
            pos_pair = tuple(sorted([ci_sig, cj_sig]))
            edge_key = (a, b, pos_pair)
            if edge_key in seen_edges:
                if log:
                    log.debug("Skipping duplicate interface edge %s<->%s (key=%s)",
                            interface.chain_i, interface.chain_j, edge_key)
                continue
            seen_edges.add(edge_key)

            # =========================
            # Deterministic HHT mapping
            # =========================
            # If HHT, orient by representative => (f on rep, b on non-rep)
            tmpl_i, tmpl_j = _enforce_hht_orientation(iface_template, interface.chain_i, interface.chain_j)

            # Build the two instances
            inst_i = InterfaceInstance(
                absolute_coord=ci_nm,
                interface_type=tmpl_i,
                this_mol_name=f"{interface.chain_i}_{ti}",
                partner_mol_name=f"{interface.chain_j}_{tj}",
                interface_index=getattr(tmpl_i, "interface_index", 1),
                residues=list(getattr(interface, "residues_i", [])),
                energy=getattr(interface, "energy", None),
            )
            inst_j = InterfaceInstance(
                absolute_coord=cj_nm,
                interface_type=tmpl_j,
                this_mol_name=f"{interface.chain_j}_{tj}",
                partner_mol_name=f"{interface.chain_i}_{ti}",
                interface_index=getattr(tmpl_j, "interface_index", 1),
                residues=list(getattr(interface, "residues_j", [])),
                energy=getattr(interface, "energy", None),
            )

            # Pre-link partners
            inst_i.partner_interface = inst_j
            inst_j.partner_interface = inst_i

            instances.extend([inst_i, inst_j])

            if log:
                try:
                    log.info("Created bidirectional interface instances: %s (%s) <-> %s (%s)",
                            inst_i.get_name(), inst_i.interface_type.get_name(), inst_j.get_name(), inst_j.interface_type.get_name())
                except Exception:
                    log.info("Created bidirectional interface instances: (%s -> %s) and (%s -> %s)",
                            inst_i.this_mol_name, inst_i.partner_mol_name,
                            inst_j.this_mol_name, inst_j.partner_mol_name)

        if log:
            log.info("Created %d interface instances", len(instances))
        return instances



    def _is_homodimeric_heterotypic_interface(self, interface_instance: InterfaceInstance) -> bool:
        """Check if an interface instance represents a hom het interface.

        Args:
            interface_instance: Interface instance to check.

        Returns:
            True if this is a hom het interface.
        """
        # Extract molecule type from names (remove instance part)
        this_mol_type = interface_instance.this_mol_name.split('_')[-1]
        partner_mol_type = interface_instance.partner_mol_name.split('_')[-1]

        # Failed homotypic: same molecule type but different interface indices
        return (this_mol_type == partner_mol_type and
                hasattr(interface_instance.interface_type, 'signature') and
                interface_instance.interface_type.signature.get('interaction_type') == 'failed_homotypic')

    def _find_failed_homotypic_partner(self, interface_instance: InterfaceInstance,
                                       interface_instances_by_key: Dict) -> Optional[InterfaceInstance]:
        """Find the complementary partner for a failed homotypic interface.

        For failed homotypic interfaces:
        - A_A_1 partners with A_A_2
        - A_A_3 partners with A_A_4

        Args:
            interface_instance: Interface instance to find partner for.
            interface_instances_by_key: Lookup dictionary of interface instances.

        Returns:
            Partner interface instance or None.
        """
        current_index = interface_instance.interface_index

        # Determine the complementary index
        if current_index % 2 == 1:  # Odd index (1, 3, 5...)
            partner_index = current_index + 1  # Look for next even index
        else:  # Even index (2, 4, 6...)
            partner_index = current_index - 1  # Look for previous odd index

        # Look for the complementary interface
        # For failed homotypic, both this_mol_name and partner_mol_name should be the same
        # but we need to find the one with the complementary index
        for key, partner_candidate in interface_instances_by_key.items():
            candidate_this_mol, candidate_partner_mol, candidate_index = key

            # Check if this is the complementary interface:
            # 1. Same molecule types (since it's failed homotypic)
            # 2. Complementary index
            # 3. Reverse direction (this becomes partner, partner becomes this)
            if (candidate_this_mol == interface_instance.partner_mol_name and
                candidate_partner_mol == interface_instance.this_mol_name and
                    candidate_index == partner_index):
                return partner_candidate

        return None

    def _get_expected_partner_key(self, interface_instance: InterfaceInstance) -> Tuple:
        """Get the expected partner key for debugging purposes.

        Args:
            interface_instance: Interface instance.

        Returns:
            Expected partner key tuple.
        """
        if self._is_homodimeric_heterotypic_interface(interface_instance):
            current_index = interface_instance.interface_index
            if current_index % 2 == 1:
                partner_index = current_index + 1
            else:
                partner_index = current_index - 1

            return (interface_instance.partner_mol_name,
                    interface_instance.this_mol_name,
                    partner_index)
        else:
            return (interface_instance.partner_mol_name,
                    interface_instance.this_mol_name,
                    interface_instance.interface_index)

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

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Interface instances keys: %s", list(
                    interface_instances_by_key.keys())
            )

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

        # Establish partner_interface cross-references with special handling for failed homotypic
        for interface_instance in self.interface_instances:
            if interface_instance.partner_interface:
                # Already linked during creation
                continue

            # Determine if this is a failed homotypic case
            is_hom_het = self._is_homodimeric_heterotypic_interface(
                interface_instance)
            
            print("-------------")
            print(f"DEBUG: Interface {interface_instance.get_name()} is hom het!")

            if is_hom_het:
                # For failed homotypic, find the complementary interface
                partner_interface = self._find_failed_homotypic_partner(
                    interface_instance, interface_instances_by_key
                )
            else:
                # For regular heterotypic, look for reverse direction
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
                    expected_key = self._get_expected_partner_key(
                        interface_instance)
                    self.workspace_manager.logger.info(
                        "No partner interface found for %s (looking for key: %s)",
                        interface_instance.get_name(), expected_key
                    )

        # Build interfaces_neighbors_map for molecule instances
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
        for mol_instance in self.system.molecule_instances:
            print("======================")
            print(mol_instance.name)
            print(mol_instance.com)
            print(type(mol_instance.interfaces_neighbors_map))
            for (intf, neighbor) in mol_instance.interfaces_neighbors_map.items():
                print(f"{intf.get_name()}:{neighbor.name}, {intf.absolute_coord}")
        return exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=box_nm,
            parms_overrides=parms_overrides
        )
