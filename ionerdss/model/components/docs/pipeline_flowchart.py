import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
import numpy as np


# Create figure with even larger size
fig, ax = plt.subplots(1, 1, figsize=(22, 26))

# Define colors for different types of operations
colors = {
    'input': '#FFE5B4',      # Peach - Input/Output
    'process': '#B4E5FF',    # Light Blue - Processing
    'decision': '#FFB4B4',   # Light Red - Decision
    'data': '#B4FFB4',       # Light Green - Data structures
    'output': '#E5B4FF',     # Light Purple - Final output
    'validation': '#FFFFB4'   # Light Yellow - Validation
}

# Positions remain the same
positions = {
    'pdb_file': (11, 25),
    'parse_structure': (11, 23.5),
    'extract_chains': (11, 22.5),
    'calculate_properties': (11, 21.5),
    'bounding_box_filter': (11, 20),
    'interface_detection': (11, 18.5),
    'kdtree_query': (11, 17.5),
    'residue_cutoff_check': (11, 16.5),
    'interface_decision': (11, 15),
    'record_interface': (17, 15),
    'skip_interface': (5, 15),
    'grouping_mode_decision': (11, 13),
    'header_grouping': (2, 11),
    'sequence_grouping': (11, 11),
    'structure_grouping': (20, 11),
    'header_available': (2, 9.5),
    'fallback_sequence': (6, 9.5),
    'build_mol_templates': (11, 8),
    'calculate_signatures': (11, 7),
    'homodimer_check': (11, 5.5),
    'create_shared_template': (5, 4),
    'create_separate_templates': (17, 4),
    'regularize_geometry': (11, 2.5),
    'steric_clash_check': (11, 1),
    'detect_clashes': (17, 1),
    'skip_clashes': (5, 1),
    'create_instances': (11, -0.5),
    'populate_registries': (11, -1.5),
    'rebuild_references': (11, -2.5),
    'validate_system': (11, -3.5),
    'final_system': (11, -5),
}

# Box dimensions
box_dims = {
    'pdb_file': (2.8, 0.7),
    'parse_structure': (2.8, 0.7),
    'extract_chains': (3.0, 0.7),
    'calculate_properties': (3.2, 0.7),
    'bounding_box_filter': (2.8, 0.7),
    'interface_detection': (3.0, 0.7),
    'kdtree_query': (2.8, 0.7),
    'residue_cutoff_check': (2.8, 0.7),
    'interface_decision': (3.5, 1.2),  # Larger for diamond
    'record_interface': (2.8, 0.7),
    'skip_interface': (2.0, 0.7),
    'grouping_mode_decision': (3.0, 1.2),  # Larger for diamond
    'header_grouping': (2.5, 0.7),
    'sequence_grouping': (2.5, 0.7),
    'structure_grouping': (2.5, 0.7),
    'header_available': (2.5, 1.2),  # Larger for diamond
    'fallback_sequence': (2.2, 0.7),
    'build_mol_templates': (3.5, 0.7),
    'calculate_signatures': (3.0, 0.7),
    'homodimer_check': (3.2, 1.2),  # Larger for diamond
    'create_shared_template': (2.8, 0.7),
    'create_separate_templates': (3.0, 0.7),
    'regularize_geometry': (3.0, 0.7),
    'steric_clash_check': (3.2, 1.2),  # Larger for diamond
    'detect_clashes': (2.5, 0.7),
    'skip_clashes': (2.2, 0.7),
    'create_instances': (3.2, 0.7),
    'populate_registries': (2.8, 0.7),
    'rebuild_references': (3.0, 0.7),
    'validate_system': (2.8, 0.7),
    'final_system': (3.0, 0.7),
}

def create_diamond_shape(center_x, center_y, width, height):
    """Create a proper diamond shape using four points."""
    # Calculate the four points of the diamond
    top = (center_x, center_y + height/2)
    right = (center_x + width/2, center_y)
    bottom = (center_x, center_y - height/2)
    left = (center_x - width/2, center_y)
    
    # Return the points in order
    return [top, right, bottom, left]

def get_box_edges(pos, dims, is_diamond=False):
    """Calculate the edges of a box or diamond for connection points."""
    x, y = pos
    w, h = dims
    
    if is_diamond:
        # For diamond shapes, use the actual diamond points
        return {
            'top': (x, y + h/2),
            'bottom': (x, y - h/2),
            'left': (x - w/2, y),
            'right': (x + w/2, y),
        }
    else:
        # Regular rectangular edges
        return {
            'top': (x, y + h/2),
            'bottom': (x, y - h/2),
            'left': (x - w/2, y),
            'right': (x + w/2, y),
        }

def is_decision_box(name):
    """Check if a box should be a decision diamond."""
    decision_names = [
        'interface_decision', 
        'grouping_mode_decision', 
        'header_available', 
        'homodimer_check', 
        'steric_clash_check'
    ]
    return name in decision_names

def find_best_connection_points(start_name, end_name):
    """Find the best edge points to connect two boxes."""
    start_pos = positions[start_name]
    end_pos = positions[end_name]
    start_dims = box_dims[start_name]
    end_dims = box_dims[end_name]
    
    # Check if either is a decision (diamond)
    start_is_diamond = is_decision_box(start_name)
    end_is_diamond = is_decision_box(end_name)
    
    start_edges = get_box_edges(start_pos, start_dims, start_is_diamond)
    end_edges = get_box_edges(end_pos, end_dims, end_is_diamond)
    
    # Calculate direction
    dx = end_pos[0] - start_pos[0]
    dy = end_pos[1] - start_pos[1]
    
    # Choose connection points based on relative positions
    if abs(dx) > abs(dy):  # More horizontal movement
        if dx > 0:  # End is to the right
            return start_edges['right'], end_edges['left']
        else:  # End is to the left
            return start_edges['left'], end_edges['right']
    else:  # More vertical movement
        if dy > 0:  # End is above
            return start_edges['top'], end_edges['bottom']
        else:  # End is below
            return start_edges['bottom'], end_edges['top']

def create_flowchart_box(ax, name, text, box_type):
    """Create a flowchart box with appropriate styling."""
    pos = positions[name]
    dims = box_dims[name]
    x, y = pos
    width, height = dims
    color = colors[box_type]
    
    if box_type == 'decision':
        # Create proper diamond shape using Polygon
        diamond_points = create_diamond_shape(x, y, width, height)
        diamond = Polygon(diamond_points, 
                         facecolor=color, 
                         edgecolor='black', 
                         linewidth=2,
                         zorder=2)
        ax.add_patch(diamond)
    else:
        # Rectangle for other types
        box = FancyBboxPatch((x-width/2, y-height/2), width, height,
                           boxstyle="round,pad=0.1", 
                           facecolor=color, edgecolor='black', linewidth=1.5,
                           zorder=2)
        ax.add_patch(box)
    
    # Add text with appropriate sizing
    font_size = 13 if box_type == 'decision' else 12
    ax.text(x, y, text, ha='center', va='center', fontsize=font_size, fontweight='bold',
            zorder=3)

def create_connecting_arrow(ax, start_name, end_name, label=None, color='black', style='->'):
    """Create an arrow between box edges with proper connection."""
    try:
        start_point, end_point = find_best_connection_points(start_name, end_name)
        
        # Create arrow with proper styling
        arrow = FancyArrowPatch(start_point, end_point,
                              arrowstyle=style, 
                              color=color, linewidth=2.5, 
                              alpha=0.8, zorder=1,
                              mutation_scale=20)
        ax.add_patch(arrow)
        
        # Add label if provided
        if label:
            mid_x = (start_point[0] + end_point[0]) / 2
            mid_y = (start_point[1] + end_point[1]) / 2
            
            # Smart label positioning to avoid overlaps
            dx = end_point[0] - start_point[0]
            dy = end_point[1] - start_point[1]
            
            if abs(dx) > abs(dy):  # Horizontal arrow
                offset_y = 0.4 if dy >= 0 else -0.4
                offset_x = 0
            else:  # Vertical arrow
                offset_x = 1.0 if dx >= 0 else -1.0
                offset_y = 0
            
            ax.text(mid_x + offset_x, mid_y + offset_y, label, 
                   fontsize=11, color=color, weight='bold',
                   bbox=dict(boxstyle="round,pad=0.2", facecolor='white', 
                            alpha=0.9, edgecolor=color, linewidth=1),
                   ha='center', va='center', zorder=4)
    except Exception as e:
        print(f"Error creating arrow from {start_name} to {end_name}: {e}")

# Create all flowchart elements

# Input stage
create_flowchart_box(ax, 'pdb_file', 'PDB/mmCIF File\nInput', 'input')

# Parser stage
create_flowchart_box(ax, 'parse_structure', 'Parse Structure\n(BioPython)', 'process')
create_flowchart_box(ax, 'extract_chains', 'Extract Valid Chains\n(with amino acids)', 'process')
create_flowchart_box(ax, 'calculate_properties', 'Calculate Chain Properties\n(COM, radius, bbox)', 'process')

# Coarse graining stage
create_flowchart_box(ax, 'bounding_box_filter', 'Bounding Box\nPre-filter', 'process')
create_flowchart_box(ax, 'interface_detection', 'For Each Chain Pair\nInterface Detection', 'process')
create_flowchart_box(ax, 'kdtree_query', 'KD-Tree Query\n(distance_cutoff)', 'process')
create_flowchart_box(ax, 'residue_cutoff_check', 'Count Contacting\nResidues', 'process')
create_flowchart_box(ax, 'interface_decision', 'Both chains ≥\nresidue_cutoff?', 'decision')
create_flowchart_box(ax, 'record_interface', 'Record Interface\n(coords, residues)', 'data')
create_flowchart_box(ax, 'skip_interface', 'Skip\nInterface', 'process')

# Chain grouping stage
create_flowchart_box(ax, 'grouping_mode_decision', 'Grouping Mode?', 'decision')
create_flowchart_box(ax, 'header_grouping', 'Header-based\nGrouping', 'process')
create_flowchart_box(ax, 'sequence_grouping', 'Sequence-based\nGrouping', 'process')
create_flowchart_box(ax, 'structure_grouping', 'Structure-based\nGrouping', 'process')
create_flowchart_box(ax, 'header_available', 'Header\nAvailable?', 'decision')
create_flowchart_box(ax, 'fallback_sequence', 'Fallback to\nSequence', 'process')

# Template building stage
create_flowchart_box(ax, 'build_mol_templates', 'Build Molecule Templates\n(from group representatives)', 'process')
create_flowchart_box(ax, 'calculate_signatures', 'Calculate Geometric\nSignatures', 'process')
create_flowchart_box(ax, 'homodimer_check', 'Homodimeric\nInteraction?', 'decision')
create_flowchart_box(ax, 'create_shared_template', 'Create Shared\nInterface Template', 'data')
create_flowchart_box(ax, 'create_separate_templates', 'Create Separate\nInterface Templates', 'data')

# Regularization stage
create_flowchart_box(ax, 'regularize_geometry', 'Regularize Geometry\nAcross Groups', 'process')
create_flowchart_box(ax, 'steric_clash_check', 'Steric Clash\nMode = auto?', 'decision')
create_flowchart_box(ax, 'detect_clashes', 'Detect Steric\nClashes', 'validation')
create_flowchart_box(ax, 'skip_clashes', 'Skip Clash\nDetection', 'process')

# System building stage
create_flowchart_box(ax, 'create_instances', 'Create Molecule &\nInterface Instances', 'process')
create_flowchart_box(ax, 'populate_registries', 'Populate All\nRegistries', 'process')
create_flowchart_box(ax, 'rebuild_references', 'Rebuild Cross-\nReferences', 'process')
create_flowchart_box(ax, 'validate_system', 'Validate System\nIntegrity', 'validation')

# Output
create_flowchart_box(ax, 'final_system', 'Complete ionerdss\nSystem Object', 'output')

# Create all arrows with proper edge connections (same as before)
create_connecting_arrow(ax, 'pdb_file', 'parse_structure')
create_connecting_arrow(ax, 'parse_structure', 'extract_chains')
create_connecting_arrow(ax, 'extract_chains', 'calculate_properties')
create_connecting_arrow(ax, 'calculate_properties', 'bounding_box_filter')
create_connecting_arrow(ax, 'bounding_box_filter', 'interface_detection')
create_connecting_arrow(ax, 'interface_detection', 'kdtree_query')
create_connecting_arrow(ax, 'kdtree_query', 'residue_cutoff_check')
create_connecting_arrow(ax, 'residue_cutoff_check', 'interface_decision')

# Decision branches for interface detection
create_connecting_arrow(ax, 'interface_decision', 'record_interface', 'Yes', color='green')
create_connecting_arrow(ax, 'interface_decision', 'skip_interface', 'No', color='red')

# Convergence back to main flow
create_connecting_arrow(ax, 'record_interface', 'grouping_mode_decision')
create_connecting_arrow(ax, 'skip_interface', 'grouping_mode_decision')

# Grouping mode branches
create_connecting_arrow(ax, 'grouping_mode_decision', 'header_grouping', 'default', color='blue')
create_connecting_arrow(ax, 'grouping_mode_decision', 'sequence_grouping', 'sequence', color='blue')
create_connecting_arrow(ax, 'grouping_mode_decision', 'structure_grouping', 'structure', color='blue')

# Header grouping sub-flow
create_connecting_arrow(ax, 'header_grouping', 'header_available')
create_connecting_arrow(ax, 'header_available', 'fallback_sequence', 'No', color='red')
create_connecting_arrow(ax, 'header_available', 'build_mol_templates', 'Yes', color='green')
create_connecting_arrow(ax, 'fallback_sequence', 'sequence_grouping')

# Other grouping paths
create_connecting_arrow(ax, 'sequence_grouping', 'build_mol_templates')
create_connecting_arrow(ax, 'structure_grouping', 'build_mol_templates')

# Template building flow
create_connecting_arrow(ax, 'build_mol_templates', 'calculate_signatures')
create_connecting_arrow(ax, 'calculate_signatures', 'homodimer_check')
create_connecting_arrow(ax, 'homodimer_check', 'create_shared_template', 'Yes', color='green')
create_connecting_arrow(ax, 'homodimer_check', 'create_separate_templates', 'No', color='red')

# Convergence to regularization
create_connecting_arrow(ax, 'create_shared_template', 'regularize_geometry')
create_connecting_arrow(ax, 'create_separate_templates', 'regularize_geometry')

# Regularization flow
create_connecting_arrow(ax, 'regularize_geometry', 'steric_clash_check')
create_connecting_arrow(ax, 'steric_clash_check', 'detect_clashes', 'Yes', color='green')
create_connecting_arrow(ax, 'steric_clash_check', 'skip_clashes', 'No', color='red')

# Convergence to system building
create_connecting_arrow(ax, 'detect_clashes', 'create_instances')
create_connecting_arrow(ax, 'skip_clashes', 'create_instances')

# Final system building flow
create_connecting_arrow(ax, 'create_instances', 'populate_registries')
create_connecting_arrow(ax, 'populate_registries', 'rebuild_references')
create_connecting_arrow(ax, 'rebuild_references', 'validate_system')
create_connecting_arrow(ax, 'validate_system', 'final_system')

# Add stage labels, information boxes, and legend (same as before)
stage_labels = [
    (0.5, 23, 'PARSING\nSTAGE'),
    (0.5, 17, 'COARSE\nGRAINING\nSTAGE'),
    (0.5, 11.5, 'CHAIN\nGROUPING\nSTAGE'),
    (0.5, 6.5, 'TEMPLATE\nBUILDING\nSTAGE'),
    (0.5, 2, 'REGULARIZATION\nSTAGE'),
    (0.5, -2, 'SYSTEM\nBUILDING\nSTAGE'),
]

for x, y, label in stage_labels:
    ax.text(x, y, label, fontsize=16, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.8, edgecolor='black'),
            zorder=2)

# Add information boxes
hyperparams_text = """Key Hyperparameters:
• distance_cutoff (0.9 nm)
• residue_cutoff (2)
• matching_mode (default/sequence/structure)
• rmsd_threshold (2.0 Å)
• seq_threshold (0.5)
• steric_clash_mode (off/auto/custom)"""

ax.text(19.5, 20, hyperparams_text, fontsize=11, ha='left', va='top',
        bbox=dict(boxstyle="round,pad=0.5", facecolor='lightyellow', alpha=0.9, edgecolor='black'),
        zorder=2)

data_flow_text = """Data Structures Created:
• CoarseGrainedChain objects
• Interface objects with coordinates
• ChainGroup objects
• MoleculeType templates
• InterfaceType templates
• MoleculeInstance & InterfaceInstance
• Complete System with all registries"""

ax.text(19.5, 10, data_flow_text, fontsize=11, ha='left', va='top',
        bbox=dict(boxstyle="round,pad=0.5", facecolor='lightblue', alpha=0.9, edgecolor='black'),
        zorder=2)

# Add legend
legend_elements = [
    patches.Patch(color=colors['input'], label='Input/Output'),
    patches.Patch(color=colors['process'], label='Processing Step'),
    patches.Patch(color=colors['decision'], label='Decision Point (Diamond)'),
    patches.Patch(color=colors['data'], label='Data Creation'),
    patches.Patch(color=colors['validation'], label='Validation'),
]

ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=12)

# Set axis properties
ax.set_xlim(-1, 24)
ax.set_ylim(-6, 26.5)
ax.set_aspect('equal')
ax.axis('off')

# Add title
plt.title('PDB to NERDSS Parameter Pipeline Flowchart\nionerdss.model.pdb Module Architecture', 
          fontsize=19, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig("/Users/yueying/Workspace/ionerdss/ionerdss/model/components/docs/pipeline_flowchart.png")
