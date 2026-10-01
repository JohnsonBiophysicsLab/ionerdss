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

**Workspace Management**: Organized file structure with comprehensive logging.

## Pipeline Components

### 1. PDB Parser (`parser.py`)
**Purpose**: Download, parse, and extract structural data from PDB/mmCIF files.

**Key Functions**:
- Automatic PDB download from RCSB database (`PDBModelBuilder` fetches biological assembly 1 by default; see `pdb_file_format`)
- mmCIF and PDB format support
- Chain extraction and validation
- Coordinate system management
- Missing atom handling
- Modified residues of polymer chains kept by default (`include_modified_residues`)

**Output**: Structured chain data with atomic coordinates, sequences, and metadata.

### 2. Coarse Grainer (`coarse_graining.py`)
**Purpose**: Convert atomic structures to coarse-grained representations and detect protein-protein interfaces.

**Key Functions**:
- Center-of-mass calculation for protein chains
- Radius estimation as the RMS distance of the chain's atoms from its center of mass
- Interface detection from Cα–Cα contacts within `interface_detect_distance_cutoff` (default 0.9 nm), requiring at least `interface_detect_n_residue_cutoff` (default 2) contacting residues on each chain
- Binding site identification
- Optional binding-energy prediction with ProAffinity-GNN (`predict_affinity=True`, see [docs/Proaffinity.md](../../../docs/Proaffinity.md)); otherwise each interface gets the default ΔG = −16 RT

**Output**: Coarse-grained chains with centers of mass, radii, and detected interfaces.

### 3. Chain Grouper (`chain_grouping.py`)
**Purpose**: Group similar protein chains to reduce system complexity and identify symmetries.

**Key Functions**:
- mmCIF entity (header) grouping with a sequence fallback (`chain_grouping_matching_mode="default"`)
- Sequence similarity analysis (`"sequence"`; `chain_grouping_seq_threshold`, default 0.5)
- Structural similarity comparison (`"structure"`; Cα RMSD `chain_grouping_rmsd_threshold`, default 2.0 Å)
- Combined sequence and structure matching (`"sequence_structure"`), which keeps quasi-equivalent conformers of one sequence in separate groups
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
- Geometric regularization (optional): cyclic (Cn) rings snapped to exact n-fold symmetry (`geometric_regularization="auto"`) and sphere projection (`is_on_sphere=True`)

**Output**: Complete molecular system with all components and relationships.

### 6. Visualizer (`visualizer.py`)
**Purpose**: Generate comprehensive visualizations for validation and analysis.

**Key Functions**:
- 3D structure plots
- Interface connectivity diagrams
- Template property analysis
- PyMOL script generation
- Coarse-grained structure export (mmCIF) and a text summary report

**Output**: Publication-ready plots, interactive visualizations, and analysis reports.

### Supporting Components

**Model Builder (`main.py`)**: `PDBModelBuilder` runs the whole pipeline, writes the system JSON and reports, and exports NERDSS files; `ionerdss.build_system_from_pdb` wraps it.

**NERDSS Exporter (`nerdss_exporter.py`)**: Writes one `.mol` file per molecule type and `parms.inp` into `nerdss_files/`.

**Structure Validation (`structure_validation.py`)**: Preflight warnings for the built system (disconnected designs, interface sites too close to a molecule's center of mass) and the one-copy-per-type structure validation setup.

**Hyperparameters (`hyperparameters.py`)**: Centralized configuration management for all pipeline parameters.

**File Manager (`file_manager.py`)**: Workspace organization, logging, and file lifecycle management.

**Units (`ionerdss/model/components/units.py`)**: Unit system management and coordinate conversions.

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
chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
groups = chain_grouper.get_groups()

for group in groups:
    print(f"Group {group.representative}: {group.members}")
    print(f"  Similarity method: {group.grouping_method}")
    print(f"  Size: {len(group)}")
```

### Multi-Interface Support
```python
# Detect multiple binding modes between same molecule types
template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, hyperparams)
interface_templates = template_builder.get_interface_templates()

# Example: AB1, AB2 for two different binding modes between A and B
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
hyperparams.interface_detect_distance_cutoff = 0.9  # nm (default)
hyperparams.interface_detect_n_residue_cutoff = 2   # residues per chain (default)

# Process structure
with WorkspaceManager("/workspace", "1ABC") as workspace:
    # Parse structure
    parser = PDBParser("1ABC", fetch_from_pdb=True, workspace_manager=workspace)
    
    # Coarse-grain
    coarse_grainer = CoarseGrainer(parser, hyperparams)
    
    # Group chains
    chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
    
    # Build templates
    template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                     hyperparams, workspace_manager=workspace)
    
    # Assemble system
    system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                 template_builder, hyperparams, 
                                 str(workspace.workspace_path), "1ABC",
                                 workspace_manager=workspace)
    
    # Generate visualizations
    visualizer = PDBVisualizer(workspace)
    viz_outputs = visualizer.visualize_all(parser, coarse_grainer, 
                                          chain_grouper, template_builder)
    
    # Export NERDSS files (keys are molecule type names, i.e. representative chain IDs such as "A")
    nerdss_outputs = system_builder.export_nerdss_files(
        molecule_counts={"A": 50, "B": 25},
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
                chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
                template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                                 hyperparams, workspace_manager=workspace)
                system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                             template_builder, hyperparams, 
                                             str(workspace_path), pdb_id,
                                             workspace_manager=workspace)
                
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
hyperparams = PDBModelHyperparameters(
    chain_grouping_matching_mode="sequence_structure",  # "default", "sequence", "structure", "sequence_structure"
    chain_grouping_seq_threshold=0.9,   # sequence identity threshold
    chain_grouping_rmsd_threshold=2.0,  # Å, Cα RMSD threshold after superposition
)

chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)

# Analyze grouping results
groups = chain_grouper.get_groups()
summary = chain_grouper.get_summary()

print(f"Grouping Summary:")
print(f"  Original chains: {coarse_grainer.get_summary()['num_chains']}")
print(f"  Final groups: {summary['num_groups']}")
print(f"  Matching mode: {summary['grouping_method']}")

# Detailed group analysis
for group in groups:
    print(f"Group {group.representative}:")
    print(f"  Members: {group.members}")
    print(f"  Method: {group.grouping_method}")
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
print(f"  Average interface size: {summary['average_interface_size']:.1f} residues per side")
print(f"  Total interface residues: {summary['total_interface_residues']}")
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
hyperparams.interface_detect_distance_cutoff = 0.9  # nm - interface detection distance (default)
hyperparams.interface_detect_n_residue_cutoff = 2   # minimum contacting residues per chain (default)

# Chain grouping parameters
hyperparams.chain_grouping_matching_mode = "default"  # "default", "sequence", "structure", "sequence_structure"
hyperparams.chain_grouping_seq_threshold = 0.5        # sequence identity
hyperparams.chain_grouping_rmsd_threshold = 2.0       # Å RMSD

# Template building parameters
hyperparams.signature_precision = 6  # decimal places for geometric signatures
hyperparams.homodimer_distance_threshold = 0.5  # nm
hyperparams.homodimer_angle_threshold = 0.2     # radians

# Advanced features
hyperparams.geometric_regularization = "off"  # "off", "auto"
hyperparams.steric_clash_mode = "off"         # "off", "auto"

# Export configuration
config = hyperparams.to_dict()
with open("hyperparams.json", "w") as f:
    json.dump(config, f, indent=2)
```

### Workspace Configuration

```python
# Custom workspace setup: logs/, structures/, outputs/ and temp/ are created automatically
with WorkspaceManager("/custom/workspace", pdb_id="1ABC") as workspace:
    # Pipeline execution
    pass
# On exit, outputs/reports/1ABC_summary.txt lists the workspace contents
```

## Output Files

### Workspace Structure

```
workspace/
├── logs/
│   └── pipeline.log                    # Comprehensive processing log
├── structures/
│   ├── downloaded/
│   │   └── 1abc-assembly1.cif         # Original structure file (1abc.cif with file_format="mmcif")
│   └── processed/
├── visualizations/
│   ├── basic_coarse_grained_structure.png
│   ├── interface_connections.png
│   ├── chain_groups.png
│   ├── template_overview.png
│   ├── 1ABC_coarse_grained.cif       # Coarse-grained structure
│   ├── 1ABC_visualization.pml        # PyMOL script
│   └── visualization_summary.txt      # Analysis report (from PDBVisualizer.generate_summary_report)
├── nerdss_files/
│   ├── A.mol                          # Molecule definition, one per molecule type
│   ├── B.mol
│   └── parms.inp                      # Simulation parameters, box size, molecule counts and reactions
├── outputs/
│   ├── systems/
│   │   └── 1ABC_system.json           # Serialized System (PDBModelBuilder)
│   └── reports/
│       ├── 1ABC_validation.txt        # System validation report (PDBModelBuilder)
│       ├── 1ABC_detailed_summary.txt  # Per-step pipeline summary (PDBModelBuilder)
│       └── 1ABC_summary.txt           # Workspace summary (written when a WorkspaceManager context exits)
├── ode_results/                       # ODE pipeline output (only with ode_enabled=True)
└── temp/
```

### Key Output Files

**Visualization Files**:
- `basic_coarse_grained_structure.png`: 3D plot of molecular centers and interfaces
- `interface_connections.png`: Interface connectivity diagram
- `chain_groups.png`: Color-coded chain groups
- `template_overview.png`: Template property analysis dashboard
- `1ABC_visualization.pml`: PyMOL script for interactive visualization

**NERDSS Simulation Files**:
- `*.mol`: Molecular template definitions with binding sites, named after the molecule type
- `parms.inp`: Simulation parameters (timestep, iterations, etc.), water box size, molecule counts and reactions

**Analysis Reports**:
- `visualization_summary.txt`: Comprehensive analysis report
- `1ABC_validation.txt` and `1ABC_detailed_summary.txt`: System validation and per-step summary written by `PDBModelBuilder`
- `1ABC_summary.txt`: Listing of the workspace contents
- `pipeline.log`: Detailed processing log with timestamps

**System File**:
- `1ABC_system.json`: The built System (molecule/interface types and instances) serialized to JSON

## Advanced Features

### Ring Regularization

```python
# Enable ring regularization for cyclic structures
hyperparams = PDBModelHyperparameters()
hyperparams.geometric_regularization = "auto"  # detect Cn rings and snap them to exact n-fold symmetry
hyperparams.symmetry_fold_tolerance = 0.15     # accept n-fold within this fraction of the assembly radius
# hyperparams.is_on_sphere = True              # separately: project molecules onto concentric spheres

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

# Template analysis and text report
template_plot = visualizer.plot_template_overview(template_builder)
summary_report = visualizer.generate_summary_report(coarse_grainer, chain_grouper, template_builder)
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
            chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
            template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                             hyperparams, workspace_manager=workspace)
            system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                         template_builder, hyperparams, 
                                         str(workspace.workspace_path), pdb_id,
                                         workspace_manager=workspace)
            
            # Collect metrics
            results[pdb_id] = {
                "chains": len(parser.get_chain_ids()),
                "groups": len(chain_grouper.get_groups()),
                "interfaces": len(coarse_grainer.get_interfaces()),
                "templates": len(template_builder.get_molecule_templates()),
                "validation": system_builder.validate_system()
            }
    
    # Generate comparative report (your own helper; not part of ionerdss)
    generate_comparative_report(results, output_dir / "analysis_report.txt")
    return results
```

## Integration Examples

### Tutorial: Build a 6BNO-Based NERDSS Model

This walkthrough uses the current high-level PDB pipeline to:

- download and parse `6BNO`
- tune interface detection and chain grouping
- export `nerdss_files/`
- launch a NERDSS simulation with `parms.inp`

```python
from pathlib import Path
import subprocess

from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.parser import PDBParser

# Choose an output workspace for all generated files
save_folder = Path("6bno_tutorial")

# Optional: inspect the deposited structure before building the coarse-grained model.
# Because "6bno" is a 4-character PDB ID, PDBParser will download it automatically.
parser = PDBParser("6bno")
print(f"PDB ID: {parser.get_pdb_id()}")
print(f"Detected chains: {parser.get_chain_ids()}")

# Configure the interface-detection and grouping settings from this tutorial.
hyperparams = PDBModelHyperparameters(
    interface_detect_distance_cutoff=1.0,
    interface_detect_n_residue_cutoff=2,
    nerdss_overlap_sep_limit=3.0,
    chain_grouping_seq_threshold=0.5,
    generate_nerdss_files=True,
    generate_visualizations=True,
    nerdss_water_box=[500.0, 500.0, 500.0],
    nerdss_total_molecule_count=75,
)

# Build the ionerdss System and export NERDSS input files into save_folder/nerdss_files.
builder = PDBModelBuilder(source="6bno", hyperparams=hyperparams)
system = builder.build_system(workspace_path=str(save_folder))

print("Molecule types:", [mol.name for mol in system.molecule_types])
print("NERDSS files written to:", save_folder / "nerdss_files")

# Run NERDSS with the generated parms.inp file.
nerdss_dir = save_folder / "nerdss_files"
nerdss_cmd = "/path/to/NERDSS/bin/nerdss -f parms.inp"  # placeholder: point this at your NERDSS build

subprocess.run(
    nerdss_cmd,
    shell=True,
    cwd=nerdss_dir,
    executable="/bin/bash",
    check=True,
)
```

What this configuration changes:

- `interface_detect_distance_cutoff=1.0`: considers Cα–Cα pairs within 1.0 nm as potential contacts when identifying interfaces. This merges all of the long-pitch contacts into a single `aa2f`/`aa2b` pair, giving the correct four binding sites; the 0.9 nm default splits the C–E and F–H contacts off into a spurious self-binding site.
- `interface_detect_n_residue_cutoff=2`: keeps interfaces that have at least two contacting residues on each side (also the default).
- `nerdss_overlap_sep_limit=3.0`: keeps subunit centres at least 3.0 nm apart during the NERDSS run. Without it, the filament mis-assembles and subunits collapse on top of one another. Values from 2.5 to 3.75 nm work for 6BNO; the pipeline caps this at 0.9 × the minimum chain COM distance (about 3.79 nm here).
- `chain_grouping_seq_threshold=0.5`: groups repeated chains when their sequence identity is at least 50% (also the default).

Where to look after the build finishes:

- `6bno_tutorial/logs/pipeline.log`: step-by-step pipeline log
- `6bno_tutorial/nerdss_files/parms.inp`: NERDSS simulation parameters
- `6bno_tutorial/nerdss_files/A.mol`: the actin molecule template with its four binding sites
- `6bno_tutorial/outputs/reports/`: validation and summary reports
- `6bno_tutorial/visualizations/`: coarse-grained plots and PyMOL helper files

If you want a shorter version, the public wrapper exposes the same hyperparameters directly:

```python
from ionerdss import build_system_from_pdb

system = build_system_from_pdb(
    source="6bno",
    workspace_path="6bno_tutorial",
    interface_detect_distance_cutoff=1.0,
    interface_detect_n_residue_cutoff=2,
    nerdss_overlap_sep_limit=3.0,
    chain_grouping_seq_threshold=0.5,
    generate_nerdss_files=True,
)
```

### Integration with NERDSS

```python
# Complete pipeline to NERDSS simulation
def pdb_to_nerdss_simulation(pdb_id, workspace_path, simulation_params):
    """Convert PDB structure to ready-to-run NERDSS simulation."""
    
    with WorkspaceManager(workspace_path, pdb_id) as workspace:
        # Process structure
        parser = PDBParser(pdb_id, fetch_from_pdb=True, workspace_manager=workspace)
        coarse_grainer = CoarseGrainer(parser, hyperparams)
        chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
        template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                         hyperparams, workspace_manager=workspace)
        system_builder = SystemBuilder(parser, coarse_grainer, chain_grouper, 
                                     template_builder, hyperparams, 
                                     str(workspace_path), pdb_id,
                                     workspace_manager=workspace)
        
        # Export NERDSS files; passing hyperparams applies nerdss_time_step,
        # default_on_rate_3d_ka, nerdss_overlap_sep_limit, etc.
        nerdss_files = system_builder.export_nerdss_files(
            molecule_counts=simulation_params["molecule_counts"],
            box_nm=simulation_params["box_size"],
            parms_overrides={**simulation_params["parameters"], "hyperparams": hyperparams}
        )
        
        # Generate run script (your own helper; not part of ionerdss)
        generate_nerdss_run_script(nerdss_files, Path(workspace_path) / "run_simulation.sh")
        
        return nerdss_files

# Usage
hyperparams.nerdss_time_step = 0.1          # μs; None (default) calculates a stable time step
hyperparams.default_on_rate_3d_ka = 1000.0  # nm³/μs, written as onRate3Dka for every reaction
simulation_config = {
    "molecule_counts": {"A": 100, "B": 50},  # keys are molecule type names
    "box_size": (500.0, 500.0, 500.0),  # nm
    "parameters": {  # written to the parameters block of parms.inp
        "nItr": 1e6,
    }
}

nerdss_files = pdb_to_nerdss_simulation("1ABC", "/simulation_workspace", simulation_config)
```

### Tutorial: Build a Synthetic Dodecahedron and Run NERDSS

The platonic-solid generator is useful when you want a designed reference assembly
instead of a structure inferred from a deposited PDB file.

```python
from pathlib import Path
import subprocess

from ionerdss import build_system_from_plat

save_folder = Path("dode_tutorial")

sys, rxn = build_system_from_plat(
    solid_type="dode",
    radius=8,
    sigma=1,
    output_nerdss=True,
    output_dir=str(save_folder),
)

print(f"Built {len(sys.molecule_types.molecule_types)} molecule type(s)")
print(f"Built {len(rxn)} reaction rule(s)")

print("\nNow running NERDSS simulation...\n")

nerdss_dir = save_folder / "nerdss_files"
nerdss_cmd = "/path/to/NERDSS/bin/nerdss -f parms.inp"  # placeholder: point this at your NERDSS build

subprocess.run(
    nerdss_cmd,
    shell=True,
    cwd=nerdss_dir,
    executable="/bin/bash",
    check=True,
)
```

In this case, `build_system_from_plat(..., output_nerdss=True)` writes the same
`nerdss_files/` directory structure, so the simulation launch step is identical to
the `6BNO` tutorial.

### Integration with Molecular Viewers

```python
# Generate files for molecular visualization
def create_visualization_package(pdb_id, workspace_path):
    """Create comprehensive visualization package."""
    
    with WorkspaceManager(workspace_path, pdb_id) as workspace:
        # Process structure
        parser = PDBParser(pdb_id, fetch_from_pdb=True, workspace_manager=workspace)
        coarse_grainer = CoarseGrainer(parser, hyperparams)
        chain_grouper = ChainGrouper(parser, coarse_grainer, hyperparams)
        template_builder = TemplateBuilder(parser, coarse_grainer, chain_grouper, 
                                         hyperparams, workspace_manager=workspace)
        
        # Generate visualizations
        visualizer = PDBVisualizer(workspace)
        viz_outputs = visualizer.visualize_all(parser, coarse_grainer, 
                                              chain_grouper, template_builder)
        
        # Create viewer-specific files (the Chimera/VMD helpers are your own; not part of ionerdss)
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
