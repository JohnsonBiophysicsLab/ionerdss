# PDB Module

## Overview

The `ionerdss.model.pdb` module provides a comprehensive pipeline for processing Protein Data Bank (PDB) structures and converting them into coarse-grained molecular models suitable for NERDSS (Numerical Evaluation of Reaction-Diffusion Spatial Stochasticity) simulations. This module transforms atomic-resolution protein structures into simplified representations that capture essential geometric and interaction properties while enabling efficient large-scale simulations.

## Table of Contents

- [Architecture Overview](#architecture-overview)
- [Pipeline Components](#pipeline-components)
- [Key Features](#key-features)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Advanced Features](#advanced-features)
- [Integration Examples](#integration-examples)

## Architecture Overview

The PDB module follows a modular, pipeline-based architecture where each component performs a specific transformation step:

```
PDB/mmCIF File → Parser → CoarseGrainer → ChainGrouper → TemplateBuilder → SystemBuilder → NERDSS Files
                    ↓         ↓             ↓              ↓              ↓
                Structure  Interfaces    Groups        Templates      System
                   Data      &COMs                                  Assembly
```

### Core Philosophy

**Modular Design**: Each component has a single responsibility and can be used independently or as part of the complete pipeline.

**Data-Driven Processing**: All transformations preserve provenance and provide detailed logging for reproducibility and debugging.

**Flexible Configuration**: Extensive hyperparameter system allows fine-tuning for different types of molecular systems.

**Workspace Management**: Organized file structure with automatic cleanup and comprehensive logging.

## Pipeline Components

### 1. PDB Parser (`parser.py`)
**Purpose**: Download, parse, and extract structural data from PDB/mmCIF files.

**Key Functions**:
- Automatic PDB download from RCSB database
- mmCIF and PDB format support
- Chain extraction and validation
- Coordinate system management
- Missing atom handling

**Output**: Structured chain data with atomic coordinates, sequences, and metadata.

### 2. Coarse Grainer (`coarse_graining.py`)
**Purpose**: Convert atomic structures to coarse-grained representations and detect protein-protein interfaces.

**Key Functions**:
- Center-of-mass calculation for protein chains
- Radius estimation using convex hull or Cα-based methods
- Interface detection via distance-based criteria
- Binding site identification
- Energy estimation for interactions

**Output**: Coarse-grained chains with centers of mass, radii, and detected interfaces.

### 3. Chain Grouper (`chain_grouping.py`)
**Purpose**: Group similar protein chains to reduce system complexity and identify symmetries.

**Key Functions**:
- Sequence similarity analysis
- Structural similarity comparison
- Symmetry detection
- Template reduction strategies
- Group validation and optimization

**Output**: Chain groups with representative chains and similarity metrics.

### 4. Template Builder (`template_builder.py`)
**Purpose**: Generate reusable molecular and interface templates with geometric signatures.

**Key Functions**:
- Molecular template creation from chain groups
- Interface template generation with geometric signatures
- Template deduplication based on similarity
- Cross-reference establishment
- Steric clash detection

**Output**: Molecular and interface templates ready for system assembly.

### 5. System Builder (`system_builder.py`)
**Purpose**: Assemble complete molecular systems from templates and instances.

**Key Functions**:
- Molecule instance creation
- Interface instance generation
- Cross-reference network establishment
- System validation
- Ring regularization (optional)

**Output**: Complete molecular system with all components and relationships.

### 6. Visualizer (`visualizer.py`)
**Purpose**: Generate comprehensive visualizations for validation and analysis.

**Key Functions**:
- 3D structure plots
- Interface connectivity diagrams
- Template property analysis
- PyMOL script generation
- Quality assessment plots

**Output**: Publication-ready plots, interactive visualizations, and analysis reports.

### Supporting Components

**Hyperparameters (`hyperparameters.py`)**: Centralized configuration management for all pipeline parameters.

**File Manager (`file_manager.py`)**: Workspace organization, logging, and file lifecycle management.

**Units (`units.py`)**: Unit system management and coordinate conversions.

## Key Features

### Automated PDB Processing
```python
from ionerdss.model.pdb.file_manager import WorkspaceManager
from ionerdss.model.pdb.parser import PDBParser

# Automatic download and processing
with WorkspaceManager("/workspace", "1ABC") as workspace:
    parser = PDBParser("1ABC", fetch_from_pdb=True, workspace_manager=workspace)
    # Structure automatically downloaded and parsed
```

### Intelligent Chain Grouping
```python
# Automatic detection of symmetric chains
chain_grouper = ChainGrouper(coarse_grainer, hyperparams)
groups = chain_grouper.get_groups()

for group in groups:
    print(f"Group {group.representative}: {group.members}")
    print(f"  Similarity method: {group.grouping_method}")
    print(f"  Average similarity: {group.average_similarity:.3f}")
```

### Multi-Interface Support
```python
# Detect multiple binding modes between same molecule types
template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, hyperparams)
interface_templates = template_builder.get_interface_templates()

# Example: A_B_1, A_B_2 for two different binding modes between A and B
for name, template in interface_templates.items():
    print(f"{name}: {template.this_mol_type_name} ↔ {template.partner_mol_type_name}")
```

### Comprehensive Validation
```python
system_builder = SystemBuilder(...)
validation_results = system_builder.validate_system()

if validation_results["errors"]:
    print("System validation failed:")
    for error in validation_results["errors"]:
        print(f"  - {error}")
else:
    print("System validation passed!")
```

### Rich Visualization Suite
```python
visualizer = PDBVisualizer(workspace)
outputs = visualizer.visualize_all(parser, coarse_grainer, chain_grouper, template_builder)

# Generates: structure plots, interface diagrams, template analysis, PyMOL scripts
for viz_type, file_path in outputs.items():
    print(f"{viz_type}: {file_path}")
```

## Quick Start

### Basic Pipeline Execution

```python
from ionerdss.model.pdb.file_manager import WorkspaceManager
from ionerdss.model.pdb.parser import PDBParser
from ionerdss.model.pdb.coarse_graining import CoarseGrainer
from ionerdss.model.pdb.chain_grouping import ChainGrouper
from ionerdss.model.pdb.template_builder import TemplateBuilder
from ionerdss.model.pdb.system_builder import SystemBuilder
from ionerdss.model.pdb.visualizer import PDBVisualizer
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

# Configure parameters
hyperparams = PDBModelHyperparameters()
hyperparams.distance_cutoff = 0.6  # nm
hyperparams.residue_cutoff = 3

# Process structure
with WorkspaceManager("/workspace", "1ABC") as workspace:
    # Parse structure
    parser = PDBParser("1ABC", fetch_from_pdb=True, workspace_manager=workspace)
    
    # Coarse-grain
    coarse_grainer = CoarseGrainer(parser, hyperparams)
    
    # Group chains
    chain_grouper = ChainGrouper(coarse_grainer, hyperparams)
    
    # Build templates
    template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                     hyperparams, workspace_manager=workspace)
    
    # Assemble system
    system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                 template_builder, hyperparams, 
                                 str(workspace.workspace_path), "1ABC", workspace)
    
    # Generate visualizations
    visualizer = PDBVisualizer(workspace)
    viz_outputs = visualizer.visualize_all(parser, coarse_grainer, 
                                          chain_grouper, template_builder)
    
    # Export NERDSS files
    nerdss_outputs = system_builder.export_nerdss_files(
        molecule_counts={"ProteinA": 50, "ProteinB": 25},
        box_nm=(200.0, 200.0, 200.0)
    )
    
    print("Pipeline completed successfully!")
    print(f"Generated {len(viz_outputs)} visualizations")
    print(f"Generated {len(nerdss_outputs)} NERDSS files")
```

### Batch Processing

```python
def process_pdb_batch(pdb_ids, workspace_base):
    """Process multiple PDB structures."""
    results = {}
    
    for pdb_id in pdb_ids:
        try:
            workspace_path = workspace_base / pdb_id
            with WorkspaceManager(workspace_path, pdb_id) as workspace:
                # Run complete pipeline
                parser = PDBParser(pdb_id, fetch_from_pdb=True, workspace_manager=workspace)
                coarse_grainer = CoarseGrainer(parser, hyperparams)
                chain_grouper = ChainGrouper(coarse_grainer, hyperparams)
                template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                                 hyperparams, workspace_manager=workspace)
                system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                             template_builder, hyperparams, 
                                             str(workspace_path), pdb_id, workspace)
                
                # Collect results
                results[pdb_id] = {
                    "system_summary": system_builder.get_summary(),
                    "validation": system_builder.validate_system(),
                    "workspace": str(workspace_path)
                }
                
        except Exception as e:
            print(f"Failed to process {pdb_id}: {e}")
            results[pdb_id] = {"error": str(e)}
    
    return results

# Usage
pdb_list = ["1ABC", "2DEF", "3GHI"]
results = process_pdb_batch(pdb_list, Path("/batch_workspace"))
```

## Detailed Usage

### Advanced Chain Grouping

```python
# Custom grouping parameters
hyperparams = PDBModelHyperparameters()
# Configure sequence similarity threshold
# Configure structure similarity threshold
# Configure minimum group size

chain_grouper = ChainGrouper(coarse_grainer, hyperparams)

# Analyze grouping results
groups = chain_grouper.get_groups()
summary = chain_grouper.get_summary()

print(f"Grouping Summary:")
print(f"  Original chains: {summary['num_chains']}")
print(f"  Final groups: {summary['num_groups']}")
print(f"  Reduction ratio: {summary['reduction_ratio']:.2f}")

# Detailed group analysis
for group in groups:
    print(f"Group {group.representative}:")
    print(f"  Members: {group.members}")
    print(f"  Method: {group.grouping_method}")
    if hasattr(group, 'similarity_matrix'):
        print(f"  Avg similarity: {group.similarity_matrix.mean():.3f}")
```

### Interface Analysis

```python
# Detailed interface analysis
interfaces = coarse_grainer.get_interfaces()

print(f"Detected {len(interfaces)} interfaces:")
for i, interface in enumerate(interfaces):
    print(f"Interface {i+1}: {interface.chain_i} ↔ {interface.chain_j}")
    print(f"  Distance: {np.linalg.norm(interface.coord_i - interface.coord_j):.2f} Å")
    print(f"  Energy: {interface.energy:.2f}")
    print(f"  Residues i: {len(interface.residues_i)}")
    print(f"  Residues j: {len(interface.residues_j)}")

# Interface statistics
summary = coarse_grainer.get_summary()
print(f"Interface Statistics:")
print(f"  Total interfaces: {summary['num_interfaces']}")
print(f"  Average energy: {summary['avg_interface_energy']:.2f}")
print(f"  Energy range: {summary['energy_range']}")
```

### Template Customization

```python
# Custom template building
template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                 hyperparams, workspace_manager=workspace)

# Analyze templates
mol_templates = template_builder.get_molecule_templates()
intf_templates = template_builder.get_interface_templates()

print(f"Molecular Templates:")
for name, template in mol_templates.items():
    print(f"  {name}: radius={template.radius_nm:.3f} nm")
    print(f"    D_trans={template.D_t_nm2_us:.2e} nm²/μs")
    print(f"    D_rot={template.D_r_rad2_us:.2e} rad²/μs")

print(f"Interface Templates:")
for name, template in intf_templates.items():
    print(f"  {name}: {template.this_mol_type_name} ↔ {template.partner_mol_type_name}")
    print(f"    Energy: {template.energy:.2f}")
    print(f"    Index: {template.interface_index}")
```

### System Analysis

```python
# Comprehensive system analysis
system = system_builder.get_system()
summary = system_builder.get_summary()

print(f"System Summary:")
print(f"  Molecule types: {len(system.molecule_types)}")
print(f"  Interface types: {len(system.interface_types)}")
print(f"  Molecule instances: {len(system.molecule_instances)}")
print(f"  Interface instances: {len(system.interface_instances)}")

# Validation
validation = system_builder.validate_system()
print(f"Validation Results:")
print(f"  Errors: {len(validation['errors'])}")
print(f"  Warnings: {len(validation['warnings'])}")

if validation['errors']:
    for error in validation['errors']:
        print(f"    ERROR: {error}")

if validation['warnings']:
    for warning in validation['warnings']:
        print(f"    WARNING: {warning}")
```

## Configuration

### Hyperparameter Configuration

```python
# Create custom hyperparameters
hyperparams = PDBModelHyperparameters()

# Coarse-graining parameters
hyperparams.distance_cutoff = 0.6  # nm - interface detection distance
hyperparams.residue_cutoff = 3     # minimum residues for interface

# Chain grouping parameters (if available)
# hyperparams.sequence_similarity_threshold = 0.9
# hyperparams.structure_similarity_threshold = 2.0  # Å RMSD

# Template building parameters
hyperparams.signature_precision = 6  # decimal places for geometric signatures
hyperparams.homodimer_distance_threshold = 1.0  # Å
hyperparams.homodimer_angle_threshold = 0.2     # radians

# Advanced features
hyperparams.ring_regularization_mode = "off"  # "off", "separate", "uniform"
hyperparams.steric_clash_mode = "off"         # "off", "auto"

# Export configuration
config = hyperparams.to_dict()
with open("hyperparams.json", "w") as f:
    json.dump(config, f, indent=2)
```

### Workspace Configuration

```python
# Custom workspace setup
workspace_config = {
    "base_path": "/custom/workspace",
    "pdb_id": "1ABC",
    "create_subdirs": True,
    "cleanup_temp": True,
    "log_level": "INFO"
}

with WorkspaceManager(**workspace_config) as workspace:
    # Pipeline execution
    pass
```

## Output Files

### Workspace Structure

```
workspace/
├── logs/
│   └── pipeline.log                    # Comprehensive processing log
├── structures/
│   └── downloaded/
│       └── 1ABC.cif                   # Original structure file
├── visualizations/
│   ├── basic_coarse_grained_structure.png
│   ├── interface_connections.png
│   ├── chain_groups.png
│   ├── template_overview.png
│   ├── 1ABC_coarse_grained.cif       # Coarse-grained structure
│   ├── 1ABC_visualization.pml        # PyMOL script
│   └── visualization_summary.txt      # Analysis report
├── nerdss_files/
│   ├── ProteinA.mol                   # Molecule definition
│   ├── ProteinB.mol
│   ├── parms.inp                      # Simulation parameters
│   └── system.inp                     # System configuration
└── outputs/
    └── reports/
        └── 1ABC_summary.txt           # Pipeline summary
```

### Key Output Files

**Visualization Files**:
- `basic_coarse_grained_structure.png`: 3D plot of molecular centers and interfaces
- `interface_connections.png`: Interface connectivity diagram
- `chain_groups.png`: Color-coded chain groups
- `template_overview.png`: Template property analysis dashboard
- `1ABC_visualization.pml`: PyMOL script for interactive visualization

**NERDSS Simulation Files**:
- `*.mol`: Molecular template definitions with binding sites
- `parms.inp`: Simulation parameters (timestep, iterations, etc.)
- `system.inp`: System configuration (box size, molecule counts)

**Analysis Reports**:
- `visualization_summary.txt`: Comprehensive analysis report
- `1ABC_summary.txt`: Pipeline execution summary
- `pipeline.log`: Detailed processing log with timestamps

## Advanced Features

### Ring Regularization

```python
# Enable ring regularization for cyclic structures
hyperparams = PDBModelHyperparameters()
hyperparams.ring_regularization_mode = "uniform"  # or "separate"
hyperparams.ring_geometry = "cylinder"            # or "sphere"

# Ring regularization automatically applied during system building
system_builder = SystemBuilder(...)
# Regularized coordinates integrated into final system
```

### Steric Clash Detection

```python
# Enable automatic steric clash detection
hyperparams.steric_clash_mode = "auto"

# Clashes automatically detected and marked as mutually exclusive
template_builder = TemplateBuilder(...)
interface_templates = template_builder.get_interface_templates()

for name, template in interface_templates.items():
    if template.required_free:
        print(f"{name} conflicts with: {template.required_free}")
```

### Custom Visualization

```python
# Generate custom visualizations
visualizer = PDBVisualizer(workspace)

# Individual visualization types
basic_plot = visualizer.plot_basic_coarse_grained_structure(coarse_grainer, figsize=(16, 12))
interface_plot = visualizer.plot_interface_connections(coarse_grainer)
groups_plot = visualizer.plot_chain_groups(coarse_grainer, chain_grouper)

# Custom analysis plots
quality_metrics = visualizer.plot_quality_metrics(coarse_grainer, chain_grouper)
```

### Batch Analysis

```python
# Analyze multiple structures
def analyze_pdb_set(pdb_ids, output_dir):
    """Comparative analysis of multiple PDB structures."""
    results = {}
    
    for pdb_id in pdb_ids:
        with WorkspaceManager(output_dir / pdb_id, pdb_id) as workspace:
            # Run pipeline
            parser = PDBParser(pdb_id, fetch_from_pdb=True, workspace_manager=workspace)
            coarse_grainer = CoarseGrainer(parser, hyperparams)
            chain_grouper = ChainGrouper(coarse_grainer, hyperparams)
            template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                             hyperparams, workspace_manager=workspace)
            system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                         template_builder, hyperparams, 
                                         str(workspace.workspace_path), pdb_id, workspace)
            
            # Collect metrics
            results[pdb_id] = {
                "chains": len(parser.get_chain_ids()),
                "groups": len(chain_grouper.get_groups()),
                "interfaces": len(coarse_grainer.get_interfaces()),
                "templates": len(template_builder.get_molecule_templates()),
                "validation": system_builder.validate_system()
            }
    
    # Generate comparative report
    generate_comparative_report(results, output_dir / "analysis_report.txt")
    return results
```

## Integration Examples

### Integration with NERDSS

```python
# Complete pipeline to NERDSS simulation
def pdb_to_nerdss_simulation(pdb_id, workspace_path, simulation_params):
    """Convert PDB structure to ready-to-run NERDSS simulation."""
    
    with WorkspaceManager(workspace_path, pdb_id) as workspace:
        # Process structure
        parser = PDBParser(pdb_id, fetch_from_pdb=True, workspace_manager=workspace)
        coarse_grainer = CoarseGrainer(parser, hyperparams)
        chain_grouper = ChainGrouper(coarse_grainer, hyperparams)
        template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                         hyperparams, workspace_manager=workspace)
        system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                     template_builder, hyperparams, 
                                     str(workspace_path), pdb_id, workspace)
        
        # Export NERDSS files
        nerdss_files = system_builder.export_nerdss_files(
            molecule_counts=simulation_params["molecule_counts"],
            box_nm=simulation_params["box_size"],
            parms_overrides=simulation_params["parameters"]
        )
        
        # Generate run script
        generate_nerdss_run_script(nerdss_files, workspace_path / "run_simulation.sh")
        
        return nerdss_files

# Usage
simulation_config = {
    "molecule_counts": {"ProteinA": 100, "ProteinB": 50},
    "box_size": (500.0, 500.0, 500.0),  # nm
    "parameters": {
        "nItr": 1e6,
        "timestep": 0.1,
        "onRate3Dka": 1000.0
    }
}

nerdss_files = pdb_to_nerdss_simulation("1ABC", "/simulation_workspace", simulation_config)
```

### Integration with Molecular Viewers

```python
# Generate files for molecular visualization
def create_visualization_package(pdb_id, workspace_path):
    """Create comprehensive visualization package."""
    
    with WorkspaceManager(workspace_path, pdb_id) as workspace:
        # Process structure
        parser = PDBParser(pdb_id, fetch_from_pdb=True, workspace_manager=workspace)
        coarse_grainer = CoarseGrainer(parser, hyperparams)
        chain_grouper = ChainGrouper(coarse_grainer, hyperparams)
        template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                         hyperparams, workspace_manager=workspace)
        
        # Generate visualizations
        visualizer = PDBVisualizer(workspace)
        viz_outputs = visualizer.visualize_all(parser, coarse_grainer, 
                                              chain_grouper, template_builder)
        
        # Create viewer-specific files
        viewer_files = {
            "pymol": viz_outputs.get("pymol"),
            "chimera": create_chimera_script(coarse_grainer, workspace_path),
            "vmd": create_vmd_script(coarse_grainer, workspace_path),
            "coarse_grained_pdb": viz_outputs.get("cg_structure")
        }
        
        return viewer_files
```

---

*The PDB module provides a complete, flexible, and extensible framework for converting protein structures into coarse-grained models suitable for large-scale molecular simulations, with comprehensive validation, visualization, and analysis capabilities.*