"""
ionerdss.model.pdb.visualizer

Visualization utilities for inspecting coarse-grained molecular assemblies.

This module provides comprehensive visualization capabilities for PDB processing
results, including 3D plots of coarse-grained models, interface visualization,
and PyMOL script generation. All outputs are saved to the workspace directory.

## Key Features

### Comprehensive Visualization Suite

**Multi-Format Output**: Generates PNG plots, CIF structure files, PyMOL 
scripts, and text reports for complete analysis coverage.

**3D Molecular Visualization**: Creates interactive 3D plots showing molecular
centers of mass, interface positions, and connectivity patterns.

**Template Analysis**: Visualizes molecular templates with their properties
including radii, diffusion constants, and interface distributions.

**Chain Grouping Display**: Color-coded visualization of chain groups to
validate grouping algorithms and symmetry detection.

### Automated Workspace Integration

**Organized Output Structure**:
```
workspace/
└── visualizations/
    ├── basic_coarse_grained_structure.png
    ├── interface_connections.png
    ├── chain_groups.png
    ├── template_overview.png
    ├── {PDB_ID}_coarse_grained.cif
    ├── {PDB_ID}_visualization.pml
    ├── {PDB_ID}_comparison.png
    └── visualization_summary.txt
```

**Consistent Naming**: All outputs use PDB ID and descriptive names
for easy identification and organization.

**Logging Integration**: Comprehensive logging of all visualization
operations with progress tracking and error reporting.

## Visualization Types

### 1. Basic Coarse-Grained Structure

```python
visualizer.plot_basic_coarse_grained_structure(coarse_grainer)
```

**Features**:
- **Chain Centers of Mass**: Red spheres showing molecular centers
- **Interface Positions**: Blue points marking binding sites
- **COM-Interface Connections**: Gray lines connecting centers to their interfaces
- **Chain Labels**: Text labels identifying each molecular chain
- **Equal Aspect Ratio**: Proper 3D perspective for accurate spatial relationships

**Use Cases**:
- Initial structure validation
- Spatial arrangement assessment
- Interface positioning verification

**Visual Elements**:
```python
# Chain COMs: Large red spheres (s=200)
ax.scatter(coms[:, 0], coms[:, 1], coms[:, 2], 
          c='red', s=200, alpha=0.8, label='Chain COMs')

# Interface positions: Small blue spheres (s=50)
ax.scatter(interface_coords[:, 0], interface_coords[:, 1], interface_coords[:, 2],
          c='blue', s=50, alpha=0.7, label='Interfaces')

# Connection lines: Gray dashed lines
ax.plot([com_x, intf_x], [com_y, intf_y], [com_z, intf_z],
        'gray', alpha=0.5, linewidth=1)
```

### 2. Interface Connections

```python
visualizer.plot_interface_connections(coarse_grainer)
```

**Features**:
- **Interface Pairs**: Connected blue spheres showing binding partners
- **Connection Lines**: Green lines highlighting protein-protein interactions
- **Chain Centers**: Red spheres for spatial reference
- **Connection Count**: Legend showing total number of interfaces

**Applications**:
- Interaction network analysis
- Binding site validation
- Assembly pathway visualization

**Color Scheme**:
- **Red**: Molecular centers of mass
- **Blue**: Interface positions
- **Green**: Active protein-protein connections

### 3. Chain Groups Visualization

```python
visualizer.plot_chain_groups(coarse_grainer, chain_grouper)
```

**Features**:
- **Color-Coded Groups**: Each chain group gets a unique color from matplotlib's Set1 colormap
- **Group Legends**: Labels showing group representative and member count
- **Interface Overlay**: Black interface points for context
- **Symmetry Validation**: Visual confirmation of grouping algorithms

**Group Information Display**:
```python
# Each group gets distinct color
colors = plt.cm.Set1(np.linspace(0, 1, len(groups)))

# Legend format: "Group A (3 chains)"
label=f'Group {group.representative} ({len(group.members)} chains)'
```

**Benefits**:
- **Symmetry Verification**: Confirm that similar chains are grouped together
- **Grouping Quality**: Assess spatial clustering of group members
- **Template Validation**: Ensure template assignments make geometric sense

### 4. Template Overview Dashboard

```python
visualizer.plot_template_overview(template_builder)
```

**Four-Panel Analysis**:

**Panel 1: Molecule Template Radii**
- Bar chart of template radii in nanometers
- Identifies size differences between molecular types
- Validates radius calculations

**Panel 2: Diffusion Constants**
- Dual bar chart showing translational and rotational diffusion
- Compares mobility between different molecule types
- Validates Stokes-Einstein calculations

**Panel 3: Interface Count per Template**
- Bar chart of interface types per molecular template
- Shows binding site complexity
- Identifies highly connected molecules

**Panel 4: Interface Energy Distribution**
- Histogram of binding energies
- Shows energy landscape of interactions
- Identifies strong vs. weak binding sites

### 5. Coarse-Grained Structure File

```python
visualizer.save_coarse_grained_structure(coarse_grainer)
```

**CIF Format Output**:
```cif
# Coarse-grained structure CIF file
data_coarse_grained
_audit_conform_dict.text 'Coarse-grained model generated by ionerdss'
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.type_symbol
ATOM      1  COM  MOL A    0.000   0.000   0.000  1.00  0.00  C
ATOM      2  INT  INT A   10.000   0.000   0.000  1.00  0.00  O
```

**Atom Types**:
- **COM**: Centers of mass (Carbon atoms for visualization)
- **INT**: Interface positions (Oxygen atoms for distinction)

**Applications**:
- Import into molecular viewers (ChimeraX, VMD, PyMOL)
- Further analysis with structural biology tools
- Comparison with original atomic structures

### 6. PyMOL Visualization Script

```python
visualizer.generate_pymol_script(parser, coarse_grainer, template_builder)
```

**Generated Script Features**:

**Structure Loading**:
```python
# Load original structure
load /path/to/original.pdb, original
hide everything, original
show cartoon, original
spectrum chain, original
set cartoon_transparency, 0.7, original

# Load coarse-grained structure
load /path/to/coarse_grained.cif, coarse_grained
show spheres, name COM
show spheres, name INT
color red, name COM
color blue, name INT
```

**Interactive Elements**:
```python
# Create pseudoatoms for better control
pseudoatom com_A, pos=[0.000, 0.000, 0.000], color=red, label=ProteinA

# Create distance measurements
distance interface1, int_A_1, int_B_1
set dash_width, 4, interface1
color green, interface1
```

**Multiple Views**:
- **Comparison view**: Original + coarse-grained overlay
- **Side view**: Profile perspective
- **Top view**: Plan perspective
- **High-resolution images**: 1200x1200 at 300 DPI


## Usage Examples

### Basic Visualization Generation

```python
from ionerdss.model.pdb.visualizer import PDBVisualizer

# Initialize visualizer
visualizer = PDBVisualizer(workspace_manager)

# Generate all visualizations
outputs = visualizer.visualize_all(
    parser=parser,
    coarse_grainer=coarse_grainer,
    chain_grouper=chain_grouper,
    template_builder=template_builder
)

# Print generated files
for viz_type, file_path in outputs.items():
    print(f"{viz_type}: {file_path}")

# Output:
# basic_cg: /workspace/visualizations/basic_coarse_grained_structure.png
# interfaces: /workspace/visualizations/interface_connections.png
# groups: /workspace/visualizations/chain_groups.png
# templates: /workspace/visualizations/template_overview.png
# pymol: /workspace/visualizations/1ABC_visualization.pml
# cg_structure: /workspace/visualizations/1ABC_coarse_grained.cif
```

### Individual Visualization Types

```python
# Create specific visualizations
basic_plot = visualizer.plot_basic_coarse_grained_structure(coarse_grainer)
interface_plot = visualizer.plot_interface_connections(coarse_grainer)
groups_plot = visualizer.plot_chain_groups(coarse_grainer, chain_grouper)
template_plot = visualizer.plot_template_overview(template_builder)

# Generate structure files
cif_file = visualizer.save_coarse_grained_structure(coarse_grainer)
pymol_script = visualizer.generate_pymol_script(parser, coarse_grainer, template_builder)

# Create summary report
summary_report = visualizer.generate_summary_report(
    coarse_grainer, chain_grouper, template_builder
)
```

"""

from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from .file_manager import WorkspaceManager
from .coarse_graining import CoarseGrainer
from .chain_grouping import ChainGrouper
from .template_builder import TemplateBuilder


class PDBVisualizer:
    """Visualization engine for PDB processing pipeline results.

    Generates various plots and visualization files for coarse-grained
    molecular models, saving all outputs to the workspace directory.

    Attributes:
        workspace_manager: Workspace manager for file organization.
        logger: Logger instance for this visualizer.
    """

    def __init__(self, workspace_manager: WorkspaceManager):
        """Initialize visualizer with workspace manager.

        Args:
            workspace_manager: Workspace manager for file organization.
        """
        self.workspace_manager = workspace_manager
        self.logger = workspace_manager.logger

        # Create visualization subdirectory
        self.viz_dir = workspace_manager.workspace_path / 'visualizations'
        self.viz_dir.mkdir(exist_ok=True)

        self.logger.info("Initialized PDB visualizer")

    def visualize_all(self, parser, coarse_grainer: CoarseGrainer,
                      chain_grouper: ChainGrouper, template_builder: TemplateBuilder) -> Dict[str, Path]:
        """Generate all visualization outputs.

        Args:
            parser: PDB parser with structure data.
            coarse_grainer: Coarse-grainer with processed data.
            chain_grouper: Chain grouper with group information.
            template_builder: Template builder with molecular templates.

        Returns:
            Dictionary mapping visualization types to output file paths.
        """
        self.logger.info("Generating all visualizations...")

        outputs = {}

        try:
            # 1. Basic coarse-grained structure plot
            outputs['basic_cg'] = self.plot_basic_coarse_grained_structure(
                coarse_grainer)

            # 2. Interface connections plot
            outputs['interfaces'] = self.plot_interface_connections(
                coarse_grainer)

            # 3. Chain groups visualization
            outputs['groups'] = self.plot_chain_groups(
                coarse_grainer, chain_grouper)

            # 4. Template overview
            outputs['templates'] = self.plot_template_overview(
                template_builder)

            # 5. PyMOL visualization script
            outputs['pymol'] = self.generate_pymol_script(
                parser, coarse_grainer, template_builder)

            # 6. Coarse-grained structure file
            outputs['cg_structure'] = self.save_coarse_grained_structure(
                coarse_grainer)

            self.logger.info(
                "Generated %d visualization outputs", len(outputs))

        except Exception as e:
            self.logger.error("Error generating visualizations: %s", str(e))
            raise

        return outputs

    def plot_basic_coarse_grained_structure(self, coarse_grainer: CoarseGrainer,
                                            figsize: Tuple[int, int] = (12, 10)) -> Path:
        """Plot basic coarse-grained structure with COMs and interfaces.

        Args:
            coarse_grainer: Coarse-grainer with processed data.
            figsize: Figure size tuple.

        Returns:
            Path to saved plot file.
        """
        self.logger.info("Plotting basic coarse-grained structure...")

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        chains = coarse_grainer.get_coarse_grained_chains()
        interfaces = coarse_grainer.get_interfaces()

        # Plot chain COMs
        coms = []
        chain_ids = []
        for chain_id, chain_data in chains.items():
            com_nm = chain_data.com / 10.0  # Convert to nm
            coms.append(com_nm)
            chain_ids.append(chain_id)

        if coms:
            coms = np.array(coms)
            ax.scatter(coms[:, 0], coms[:, 1], coms[:, 2],
                       c='red', s=200, alpha=0.8, label='Chain COMs')

            # Add chain labels
            for i, chain_id in enumerate(chain_ids):
                ax.text(coms[i, 0], coms[i, 1], coms[i, 2],
                        f'  {chain_id}', fontsize=10, fontweight='bold')

        # Plot interfaces
        interface_coords = []
        for interface in interfaces:
            coord_i_nm = interface.coord_i / 10.0  # Convert to nm
            coord_j_nm = interface.coord_j / 10.0  # Convert to nm
            interface_coords.extend([coord_i_nm, coord_j_nm])

        if interface_coords:
            interface_coords = np.array(interface_coords)
            ax.scatter(interface_coords[:, 0], interface_coords[:, 1], interface_coords[:, 2],
                       c='blue', s=50, alpha=0.7, label='Interfaces')

        # Plot connections between COMs and their interfaces
        for interface in interfaces:
            chain_i_data = chains[interface.chain_i]
            chain_j_data = chains[interface.chain_j]

            com_i_nm = chain_i_data.com / 10.0
            com_j_nm = chain_j_data.com / 10.0
            intf_i_nm = interface.coord_i / 10.0
            intf_j_nm = interface.coord_j / 10.0

            # Draw lines from COMs to interfaces
            ax.plot([com_i_nm[0], intf_i_nm[0]],
                    [com_i_nm[1], intf_i_nm[1]],
                    [com_i_nm[2], intf_i_nm[2]],
                    'gray', alpha=0.5, linewidth=1)

            ax.plot([com_j_nm[0], intf_j_nm[0]],
                    [com_j_nm[1], intf_j_nm[1]],
                    [com_j_nm[2], intf_j_nm[2]],
                    'gray', alpha=0.5, linewidth=1)

        ax.set_xlabel('X (nm)')
        ax.set_ylabel('Y (nm)')
        ax.set_zlabel('Z (nm)')
        ax.set_title(
            f'Coarse-Grained Structure - {self.workspace_manager.pdb_id}')
        ax.legend()

        # Set equal aspect ratio
        self._set_equal_aspect_3d(ax, coms if len(
            coms) > 0 else np.array([[0, 0, 0]]))

        # Save plot
        output_path = self.viz_dir / 'basic_coarse_grained_structure.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        self.logger.info(
            "Saved basic coarse-grained structure plot to: %s", output_path)
        return output_path

    def plot_interface_connections(self, coarse_grainer: CoarseGrainer,
                                   figsize: Tuple[int, int] = (12, 10)) -> Path:
        """Plot coarse-grained structure with interface connections highlighted.

        Args:
            coarse_grainer: Coarse-grainer with processed data.
            figsize: Figure size tuple.

        Returns:
            Path to saved plot file.
        """
        self.logger.info("Plotting interface connections...")

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        chains = coarse_grainer.get_coarse_grained_chains()
        interfaces = coarse_grainer.get_interfaces()

        # Plot chain COMs
        coms = []
        chain_ids = []
        for chain_id, chain_data in chains.items():
            com_nm = chain_data.com / 10.0  # Convert to nm
            coms.append(com_nm)
            chain_ids.append(chain_id)

        if coms:
            coms = np.array(coms)
            ax.scatter(coms[:, 0], coms[:, 1], coms[:, 2],
                       c='red', s=200, alpha=0.8, label='Chain COMs')

            # Add chain labels
            for i, chain_id in enumerate(chain_ids):
                ax.text(coms[i, 0], coms[i, 1], coms[i, 2],
                        f'  {chain_id}', fontsize=10, fontweight='bold')

        # Plot interface connections
        connection_count = 0
        for interface in interfaces:
            intf_i_nm = interface.coord_i / 10.0
            intf_j_nm = interface.coord_j / 10.0

            # Plot interface points
            ax.scatter([intf_i_nm[0]], [intf_i_nm[1]], [intf_i_nm[2]],
                       c='blue', s=80, alpha=0.8)
            ax.scatter([intf_j_nm[0]], [intf_j_nm[1]], [intf_j_nm[2]],
                       c='blue', s=80, alpha=0.8)

            # Draw connection between interface points
            ax.plot([intf_i_nm[0], intf_j_nm[0]],
                    [intf_i_nm[1], intf_j_nm[1]],
                    [intf_i_nm[2], intf_j_nm[2]],
                    'green', alpha=0.7, linewidth=2)

            connection_count += 1

        # Add legend entry for connections
        if connection_count > 0:
            ax.plot([], [], [], 'green', linewidth=2,
                    label=f'Interface Connections ({connection_count})')

        ax.set_xlabel('X (nm)')
        ax.set_ylabel('Y (nm)')
        ax.set_zlabel('Z (nm)')
        ax.set_title(
            f'Interface Connections - {self.workspace_manager.pdb_id}')
        ax.legend()

        # Set equal aspect ratio
        self._set_equal_aspect_3d(ax, coms if len(
            coms) > 0 else np.array([[0, 0, 0]]))

        # Save plot
        output_path = self.viz_dir / 'interface_connections.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        self.logger.info(
            "Saved interface connections plot to: %s", output_path)
        return output_path

    def plot_chain_groups(self, coarse_grainer: CoarseGrainer, chain_grouper: ChainGrouper,
                          figsize: Tuple[int, int] = (12, 10)) -> Path:
        """Plot coarse-grained structure colored by chain groups.

        Args:
            coarse_grainer: Coarse-grainer with processed data.
            chain_grouper: Chain grouper with group information.
            figsize: Figure size tuple.

        Returns:
            Path to saved plot file.
        """
        self.logger.info("Plotting chain groups...")

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        chains = coarse_grainer.get_coarse_grained_chains()
        groups = chain_grouper.get_groups()

        # Define colors for groups
        colors = plt.cm.Set1(np.linspace(0, 1, len(groups)))

        # Plot each group with different colors
        for group_idx, group in enumerate(groups):
            color = colors[group_idx]

            # Plot COMs for this group
            group_coms = []
            for chain_id in group.members:
                if chain_id in chains:
                    com_nm = chains[chain_id].com / 10.0
                    group_coms.append(com_nm)

            if group_coms:
                group_coms = np.array(group_coms)
                ax.scatter(group_coms[:, 0], group_coms[:, 1], group_coms[:, 2],
                           c=[color], s=200, alpha=0.8,
                           label=f'Group {group.representative} ({len(group.members)} chains)')

                # Add chain labels
                for i, chain_id in enumerate(group.members):
                    if chain_id in chains:
                        com_nm = chains[chain_id].com / 10.0
                        ax.text(com_nm[0], com_nm[1], com_nm[2],
                                f'  {chain_id}', fontsize=9, color=color, fontweight='bold')

        # Plot interfaces
        interfaces = coarse_grainer.get_interfaces()
        interface_coords = []
        for interface in interfaces:
            coord_i_nm = interface.coord_i / 10.0
            coord_j_nm = interface.coord_j / 10.0
            interface_coords.extend([coord_i_nm, coord_j_nm])

        if interface_coords:
            interface_coords = np.array(interface_coords)
            ax.scatter(interface_coords[:, 0], interface_coords[:, 1], interface_coords[:, 2],
                       c='black', s=30, alpha=0.5, label='Interfaces')

        ax.set_xlabel('X (nm)')
        ax.set_ylabel('Y (nm)')
        ax.set_zlabel('Z (nm)')
        ax.set_title(f'Chain Groups - {self.workspace_manager.pdb_id}')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        # Set equal aspect ratio
        all_coms = []
        for chain_data in chains.values():
            all_coms.append(chain_data.com / 10.0)

        if all_coms:
            self._set_equal_aspect_3d(ax, np.array(all_coms))

        # Save plot
        output_path = self.viz_dir / 'chain_groups.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        self.logger.info("Saved chain groups plot to: %s", output_path)
        return output_path

    def plot_template_overview(self, template_builder: TemplateBuilder,
                               figsize: Tuple[int, int] = (14, 10)) -> Path:
        """Plot overview of molecular templates and their properties.

        Args:
            template_builder: Template builder with molecular templates.
            figsize: Figure size tuple.

        Returns:
            Path to saved plot file.
        """
        self.logger.info("Plotting template overview...")

        molecule_templates = template_builder.get_molecule_templates()
        interface_templates = template_builder.get_interface_templates()

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)

        # 1. Molecule template sizes
        template_names = list(molecule_templates.keys())
        template_radii = [mol.radius_nm for mol in molecule_templates.values()]

        ax1.bar(template_names, template_radii, color='skyblue', alpha=0.7)
        ax1.set_title('Molecule Template Radii')
        ax1.set_ylabel('Radius (nm)')
        ax1.tick_params(axis='x', rotation=45)

        # 2. Diffusion constants
        d_trans = [mol.D_t_nm2_us for mol in molecule_templates.values()]
        d_rot = [mol.D_r_rad2_us for mol in molecule_templates.values()]

        x = np.arange(len(template_names))
        width = 0.35

        ax2.bar(x - width/2, d_trans, width,
                label='D_trans (nm²/μs)', alpha=0.7)
        ax2.bar(x + width/2, d_rot, width, label='D_rot (rad²/μs)', alpha=0.7)
        ax2.set_title('Diffusion Constants')
        ax2.set_ylabel('Diffusion Constant')
        ax2.set_xticks(x)
        ax2.set_xticklabels(template_names)
        ax2.legend()
        ax2.tick_params(axis='x', rotation=45)

        # 3. Interface count per template
        interface_counts = {}
        for template_name in template_names:
            count = sum(1 for intf_name in interface_templates.keys()
                        if intf_name.startswith(template_name + '_'))
            interface_counts[template_name] = count

        ax3.bar(interface_counts.keys(), interface_counts.values(),
                color='lightcoral', alpha=0.7)
        ax3.set_title('Interface Count per Template')
        ax3.set_ylabel('Number of Interfaces')
        ax3.tick_params(axis='x', rotation=45)

        # 4. Interface energies histogram
        energies = [intf.energy for intf in interface_templates.values()]
        if energies:
            ax4.hist(energies, bins=10, color='lightgreen',
                     alpha=0.7, edgecolor='black')
            ax4.set_title('Interface Energy Distribution')
            ax4.set_xlabel('Energy')
            ax4.set_ylabel('Count')
        else:
            ax4.text(0.5, 0.5, 'No interface data', ha='center', va='center',
                     transform=ax4.transAxes)
            ax4.set_title('Interface Energy Distribution')

        plt.tight_layout()

        # Save plot
        output_path = self.viz_dir / 'template_overview.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        self.logger.info("Saved template overview plot to: %s", output_path)
        return output_path

    def save_coarse_grained_structure(self, coarse_grainer: CoarseGrainer) -> Path:
        """Save coarse-grained structure to CIF file.

        Args:
            coarse_grainer: Coarse-grainer with processed data.

        Returns:
            Path to saved CIF file.
        """
        self.logger.info("Saving coarse-grained structure to CIF...")

        output_path = self.viz_dir / \
            f'{self.workspace_manager.pdb_id}_coarse_grained.cif'

        chains = coarse_grainer.get_coarse_grained_chains()
        interfaces = coarse_grainer.get_interfaces()

        with open(output_path, 'w', encoding='utf-8') as cif_file:
            atom_id = 1

            # Write CIF header
            cif_file.write("# Coarse-grained structure CIF file\n")
            cif_file.write("data_coarse_grained\n")
            cif_file.write(
                "_audit_conform_dict.text 'Coarse-grained model generated by ionerdss'\n")
            cif_file.write("loop_\n")
            cif_file.write("_atom_site.group_PDB\n")
            cif_file.write("_atom_site.id\n")
            cif_file.write("_atom_site.label_atom_id\n")
            cif_file.write("_atom_site.label_comp_id\n")
            cif_file.write("_atom_site.label_asym_id\n")
            cif_file.write("_atom_site.Cartn_x\n")
            cif_file.write("_atom_site.Cartn_y\n")
            cif_file.write("_atom_site.Cartn_z\n")
            cif_file.write("_atom_site.occupancy\n")
            cif_file.write("_atom_site.B_iso_or_equiv\n")
            cif_file.write("_atom_site.type_symbol\n")

            # Write COM atoms for each chain
            for chain_id, chain_data in chains.items():
                com = chain_data.com
                cif_file.write(
                    f"ATOM  {atom_id:5d}  COM  MOL {chain_id}  "
                    f"{com[0]:8.3f} {com[1]:8.3f} {com[2]:8.3f}  1.00  0.00  C\n"
                )
                atom_id += 1

            # Write interface atoms
            interface_id = 1
            for interface in interfaces:
                # Interface on chain i
                coord_i = interface.coord_i
                cif_file.write(
                    f"ATOM  {atom_id:5d}  INT  INT {interface.chain_i}  "
                    f"{coord_i[0]:8.3f} {coord_i[1]:8.3f} {coord_i[2]:8.3f}  1.00  0.00  O\n"
                )
                atom_id += 1

                # Interface on chain j
                coord_j = interface.coord_j
                cif_file.write(
                    f"ATOM  {atom_id:5d}  INT  INT {interface.chain_j}  "
                    f"{coord_j[0]:8.3f} {coord_j[1]:8.3f} {coord_j[2]:8.3f}  1.00  0.00  O\n"
                )
                atom_id += 1

                interface_id += 1

        self.logger.info("Saved coarse-grained structure to: %s", output_path)
        return output_path

    def generate_pymol_script(self, parser, coarse_grainer: CoarseGrainer,
                              template_builder: TemplateBuilder) -> Path:
        """Generate PyMOL script for visualization.

        Args:
            parser: PDB parser with original structure.
            coarse_grainer: Coarse-grainer with processed data.
            template_builder: Template builder with molecular templates.

        Returns:
            Path to saved PyMOL script.
        """
        self.logger.info("Generating PyMOL visualization script...")

        output_path = self.viz_dir / \
            f'{self.workspace_manager.pdb_id}_visualization.pml'

        # Get paths to structure files
        original_structure = parser.filepath
        cg_structure = self.viz_dir / \
            f'{self.workspace_manager.pdb_id}_coarse_grained.cif'

        chains = coarse_grainer.get_coarse_grained_chains()
        interfaces = coarse_grainer.get_interfaces()
        templates = template_builder.get_molecule_templates()

        with open(output_path, 'w', encoding='utf-8') as pml_file:
            pml_file.write(
                "# PyMOL script for coarse-grained structure visualization\n")
            pml_file.write(
                f"# Generated for PDB: {self.workspace_manager.pdb_id}\n\n")

            # Load original structure
            pml_file.write(f"# Load original structure\n")
            pml_file.write(f"load {original_structure}, original\n")
            pml_file.write("hide everything, original\n")
            pml_file.write("show cartoon, original\n")
            pml_file.write("spectrum chain, original\n")
            pml_file.write("set cartoon_transparency, 0.7, original\n\n")

            # Load coarse-grained structure if it exists
            if cg_structure.exists():
                pml_file.write(f"# Load coarse-grained structure\n")
                pml_file.write(f"load {cg_structure}, coarse_grained\n")
                pml_file.write("hide everything, coarse_grained\n")
                pml_file.write("show spheres, name COM\n")
                pml_file.write("show spheres, name INT\n")
                pml_file.write("set sphere_scale, 1.0\n")
                pml_file.write("color red, name COM\n")
                pml_file.write("color blue, name INT\n\n")

            # Create pseudoatoms for better visualization
            pml_file.write("# Create pseudoatoms for COMs and interfaces\n")
            atom_index = 1

            for chain_id, chain_data in chains.items():
                com = chain_data.com
                template_name = template_builder.get_template_name_for_group(
                    chain_id)

                # Create COM pseudoatom
                pml_file.write(
                    f"pseudoatom com_{chain_id}, pos=[{com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}], "
                    f"color=red, label={template_name or chain_id}\n"
                )

            # Create interface pseudoatoms and connections
            pml_file.write("\n# Create interface connections\n")
            for i, interface in enumerate(interfaces):
                coord_i = interface.coord_i
                coord_j = interface.coord_j

                # Create interface pseudoatoms
                pml_file.write(
                    f"pseudoatom int_{interface.chain_i}_{i}, "
                    f"pos=[{coord_i[0]:.3f}, {coord_i[1]:.3f}, {coord_i[2]:.3f}], color=blue\n"
                )
                pml_file.write(
                    f"pseudoatom int_{interface.chain_j}_{i}, "
                    f"pos=[{coord_j[0]:.3f}, {coord_j[1]:.3f}, {coord_j[2]:.3f}], color=blue\n"
                )

                # Connect COM to interfaces
                pml_file.write(
                    f"distance line{atom_index}, com_{interface.chain_i}, int_{interface.chain_i}_{i}\n"
                )
                pml_file.write(f"set dash_width, 2, line{atom_index}\n")
                pml_file.write(f"set dash_gap, 0.3, line{atom_index}\n")
                atom_index += 1

                pml_file.write(
                    f"distance line{atom_index}, com_{interface.chain_j}, int_{interface.chain_j}_{i}\n"
                )
                pml_file.write(f"set dash_width, 2, line{atom_index}\n")
                pml_file.write(f"set dash_gap, 0.3, line{atom_index}\n")
                atom_index += 1

                # Connect interfaces
                pml_file.write(
                    f"distance interface{i}, int_{interface.chain_i}_{i}, int_{interface.chain_j}_{i}\n"
                )
                pml_file.write(f"set dash_width, 4, interface{i}\n")
                pml_file.write(f"set dash_gap, 0.1, interface{i}\n")
                pml_file.write(f"color green, interface{i}\n")

            # Final visualization settings
            pml_file.write("\n# Final visualization settings\n")
            pml_file.write("set sphere_transparency, 0.2\n")
            pml_file.write("bg_color white\n")
            pml_file.write("zoom all\n")

            # Save images
            #comparison_image = self.viz_dir / \
            #    f'{self.workspace_manager.pdb_id}_comparison.png'
            #pml_file.write(f"\n# Save comparison image\n")
            #pml_file.write(f"png {comparison_image}, 1200, 1200, 300, 1\n")

            # Create views for different perspectives
            #pml_file.write("\n# Create different views\n")
            #pml_file.write("# View 1: Side view\n")
            #pml_file.write(
            #    "set_view (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -50.0, 0.0, 0.0, 0.0, 40.0, 60.0, -20.0)\n")
            #side_view_image = self.viz_dir / \
            #    f'{self.workspace_manager.pdb_id}_side_view.png'
            #pml_file.write(f"png {side_view_image}, 1200, 1200, 300, 1\n")

            #pml_file.write("\n# View 2: Top view\n")
            #pml_file.write(
            #    "set_view (1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -50.0, 0.0, 0.0, 0.0, 40.0, 60.0, -20.0)\n")
            #top_view_image = self.viz_dir / \
            #    f'{self.workspace_manager.pdb_id}_top_view.png'
            #pml_file.write(f"png {top_view_image}, 1200, 1200, 300, 1\n")

        self.logger.info("Generated PyMOL script: %s", output_path)
        self.logger.info(
            "Run 'pymol %s' to visualize the structure", output_path)

        return output_path

    def _set_equal_aspect_3d(self, ax: Axes3D, coords: np.ndarray):
        """Set equal aspect ratio for 3D plot.

        Args:
            ax: 3D matplotlib axis.
            coords: Array of coordinates to determine limits.
        """
        if len(coords) == 0:
            return

        # Find the range of coordinates
        max_range = np.array([coords[:, 0].max() - coords[:, 0].min(),
                             coords[:, 1].max() - coords[:, 1].min(),
                             coords[:, 2].max() - coords[:, 2].min()]).max() / 2.0

        # Find the center
        mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
        mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
        mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5

        # Set the limits
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

    def generate_summary_report(self, coarse_grainer: CoarseGrainer,
                                chain_grouper: ChainGrouper,
                                template_builder: TemplateBuilder) -> Path:
        """Generate a comprehensive visualization summary report.

        Args:
            coarse_grainer: Coarse-grainer with processed data.
            chain_grouper: Chain grouper with group information.
            template_builder: Template builder with molecular templates.

        Returns:
            Path to saved summary report.
        """
        self.logger.info("Generating visualization summary report...")

        output_path = self.viz_dir / 'visualization_summary.txt'

        chains = coarse_grainer.get_coarse_grained_chains()
        interfaces = coarse_grainer.get_interfaces()
        groups = chain_grouper.get_groups()
        templates = template_builder.get_molecule_templates()

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(f"Visualization Summary Report\n")
            f.write(f"{'=' * 50}\n\n")
            f.write(f"PDB ID: {self.workspace_manager.pdb_id}\n")
            f.write(f"Generated visualizations in: {self.viz_dir}\n\n")

            f.write(f"Structure Overview:\n")
            f.write(f"  - Chains: {len(chains)}\n")
            f.write(f"  - Interfaces: {len(interfaces)}\n")
            f.write(f"  - Chain Groups: {len(groups)}\n")
            f.write(f"  - Molecule Templates: {len(templates)}\n\n")

            f.write(f"Chain Details:\n")
            for chain_id, chain_data in chains.items():
                f.write(f"  {chain_id}: radius={chain_data.radius:.2f}Å, "
                        f"sequence_length={len(chain_data.sequence)}\n")

            f.write(f"\nInterface Details:\n")
            for i, interface in enumerate(interfaces):
                f.write(f"  Interface {i+1}: {interface.chain_i} <-> {interface.chain_j}, "
                        f"energy={interface.energy}, "
                        f"residues_i={len(interface.residues_i)}, "
                        f"residues_j={len(interface.residues_j)}\n")

            f.write(f"\nChain Groups:\n")
            for group in groups:
                f.write(f"  Group {group.representative}: {group.members} "
                        f"(method: {group.grouping_method})\n")

            f.write(f"\nMolecule Templates:\n")
            for name, template in templates.items():
                f.write(f"  {name}: radius={template.radius_nm:.3f}nm, "
                        f"D_trans={template.D_t_nm2_us:.2e}nm²/μs, "
                        f"D_rot={template.D_r_rad2_us:.2e}rad²/μs\n")

            f.write(f"\nGenerated Files:\n")
            viz_files = list(self.viz_dir.glob('*'))
            for viz_file in sorted(viz_files):
                if viz_file.is_file():
                    f.write(f"  - {viz_file.name}\n")

        self.logger.info(
            "Generated visualization summary report: %s", output_path)
        return output_path
