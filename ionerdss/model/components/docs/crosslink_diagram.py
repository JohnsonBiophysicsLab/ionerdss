"""
draw the class hierachy and cross-link diagram with python
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Create figure and axis with more space
fig, ax = plt.subplots(1, 1, figsize=(20, 16))

# Define colors for different categories
colors = {
    'system': '#FF6B6B',      # Red - Top level
    'registry': '#4ECDC4',     # Teal - Registries
    'type': '#45B7D1',        # Blue - Type definitions
    'instance': '#96CEB4',     # Green - Runtime instances
    'config': '#FECA57',      # Yellow - Configuration
    'utility': '#DDA0DD',     # Purple - Utility classes
    'attribute': '#F0F0F0'    # Light gray - Attributes
}

# Define positions with much better spacing
positions = {
    # Top level - System (centered at top)
    'System': (10, 14),

    # Configuration components (spread out more)
    'Units': (3, 12),
    'workspace_path': (3, 12.4),
    'pdb_id': (3, 12.8),

    # Registries (well-spaced horizontally)
    'MoleculeTypeRegistry': (2, 10),
    'InterfaceTypeRegistry': (7, 10),
    'MoleculeInstanceRegistry': (13, 8.5),
    'InterfaceInstanceRegistry': (18, 8.5),

    # Type definitions
    'MoleculeType': (2, 8),
    'InterfaceType': (7, 8),

    # Instance definitions
    'MoleculeInstance': (13, 6.5),
    'InterfaceInstance': (18, 6.5),

    # Utility classes (left side, more spaced)
    'ReactionRule': (1, 6),
    'ReactionGeometrySet': (5, 6),

    # Attributes spread across bottom with much more space
    'mol_name': (1, 3),
    'radius_nm': (3.5, 3),
    'diffusion': (6, 3),
    'coordinates': (8.5, 3),
    'energy': (11, 3),
    'binding_state': (13.5, 3),
    'interfaces_map': (16, 3),
    'com_norm': (18.5, 3),
}

# Box dimensions
box_dims = {
    'System': (2.5, 0.8),
    'Units': (1.6, 0.2),
    'workspace_path': (1.6, 0.2),
    'pdb_id': (1.8, 0.2),
    'MoleculeTypeRegistry': (2.2, 0.6),
    'InterfaceTypeRegistry': (2.2, 0.6),
    'MoleculeInstanceRegistry': (2.4, 0.6),
    'InterfaceInstanceRegistry': (2.4, 0.6),
    'MoleculeType': (2.0, 0.5),
    'InterfaceType': (2.0, 0.5),
    'MoleculeInstance': (2.2, 0.5),
    'InterfaceInstance': (2.2, 0.5),
    'ReactionRule': (1.8, 0.5),
    'ReactionGeometrySet': (2.0, 0.5),
    # Attributes
    'mol_name': (1.2, 0.6),
    'radius_nm': (1.2, 0.6),
    'diffusion': (1.2, 0.6),
    'coordinates': (1.2, 0.6),
    'energy': (1.2, 0.6),
    'binding_state': (1.2, 0.6),
    'interfaces_map': (1.2, 0.6),
    'com_norm': (1.2, 0.6),
}


def get_box_edges(pos, dims):
    """Calculate the edges of a box for connection points."""
    x, y = pos
    w, h = dims
    return {
        'top': (x, y + h/2),
        'bottom': (x, y - h/2),
        'left': (x - w/2, y),
        'right': (x + w/2, y),
        'top_left': (x - w/2, y + h/2),
        'top_right': (x + w/2, y + h/2),
        'bottom_left': (x - w/2, y - h/2),
        'bottom_right': (x + w/2, y - h/2)
    }


def find_best_connection_points(start_pos, start_dims, end_pos, end_dims):
    """Find the best edge points to connect two boxes."""
    start_edges = get_box_edges(start_pos, start_dims)
    end_edges = get_box_edges(end_pos, end_dims)

    # Determine the general direction
    dx = end_pos[0] - start_pos[0]
    dy = end_pos[1] - start_pos[1]

    # Choose connection points based on relative positions
    if abs(dx) > abs(dy):  # More horizontal
        if dx > 0:  # End is to the right
            start_point = start_edges['right']
            end_point = end_edges['left']
        else:  # End is to the left
            start_point = start_edges['left']
            end_point = end_edges['right']
    else:  # More vertical
        if dy > 0:  # End is above
            start_point = start_edges['top']
            end_point = end_edges['bottom']
        else:  # End is below
            start_point = start_edges['bottom']
            end_point = end_edges['top']

    return start_point, end_point


def create_box(ax_, pos, text_, color, dims, zorder=2):
    """Create a box with specified dimensions."""
    x, y = pos
    width, height = dims
    box = FancyBboxPatch((x-width/2, y-height/2), width, height,
                         boxstyle="round,pad=0.1",
                         facecolor=color, edgecolor='black', linewidth=1.5,
                         zorder=zorder)
    ax_.add_patch(box)
    ax_.text(x, y, text_, ha='center', va='center',
             fontsize=9, fontweight='bold', zorder=zorder+1)


def create_labeled_arrow(ax_, start_name, end_name, label, style='->',
                         color='black', linestyle='-', alpha=0.8, label_offset=(0, 0)):
    """Create an arrow between box edges with a label positioned very close to the arrow.

    Args:
        label_offset: Tuple (x_offset, y_offset) for label displacement from arrow midpoint.
                     Default (0, 0) places label directly on the arrow.
    """
    start_pos = positions[start_name]
    end_pos = positions[end_name]
    start_dims = box_dims[start_name]
    end_dims = box_dims[end_name]

    start_point, end_point = find_best_connection_points(
        start_pos, start_dims, end_pos, end_dims)

    # Create arrow
    arrow = FancyArrowPatch(start_point, end_point,
                            arrowstyle=style,
                            color=color, linestyle=linestyle, alpha=alpha,
                            linewidth=2, zorder=1)
    ax_.add_patch(arrow)

    # Position label at arrow midpoint with minimal offset
    mid_x = (start_point[0] + end_point[0]) / 2
    mid_y = (start_point[1] + end_point[1]) / 2

    # Apply custom offset (default is 0,0 for very close positioning)
    label_x = mid_x + label_offset[0]
    label_y = mid_y + label_offset[1]

    # Add label with background
    ax_.text(label_x, label_y, label, ha='center', va='center',
             fontsize=8, color=color, weight='bold',
             bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9,
                       edgecolor=color, linewidth=1),
             zorder=4)


# Create all boxes
create_box(ax, positions['System'], 'System',
           colors['system'], box_dims['System'], zorder=3)

# Configuration
create_box(ax, positions['Units'], 'Units',
           colors['config'], box_dims['Units'], zorder=3)
create_box(ax, positions['workspace_path'], 'workspace_path',
           colors['config'], box_dims['workspace_path'], zorder=3)
create_box(ax, positions['pdb_id'], 'pdb_id\n(Optional)',
           colors['config'], box_dims['pdb_id'], zorder=3)

# Registries
create_box(ax, positions['MoleculeTypeRegistry'], 'MoleculeType\nRegistry',
           colors['registry'], box_dims['MoleculeTypeRegistry'], zorder=3)
create_box(ax, positions['InterfaceTypeRegistry'], 'InterfaceType\nRegistry',
           colors['registry'], box_dims['InterfaceTypeRegistry'], zorder=3)
create_box(ax, positions['MoleculeInstanceRegistry'], 'MoleculeInstance\nRegistry',
           colors['registry'], box_dims['MoleculeInstanceRegistry'], zorder=3)
create_box(ax, positions['InterfaceInstanceRegistry'], 'InterfaceInstance\nRegistry',
           colors['registry'], box_dims['InterfaceInstanceRegistry'], zorder=3)

# Types
create_box(ax, positions['MoleculeType'], 'MoleculeType',
           colors['type'], box_dims['MoleculeType'], zorder=3)
create_box(ax, positions['InterfaceType'], 'InterfaceType',
           colors['type'], box_dims['InterfaceType'], zorder=3)

# Instances
create_box(ax, positions['MoleculeInstance'], 'MoleculeInstance',
           colors['instance'], box_dims['MoleculeInstance'], zorder=3)
create_box(ax, positions['InterfaceInstance'], 'InterfaceInstance',
           colors['instance'], box_dims['InterfaceInstance'], zorder=3)

# Utility classes
create_box(ax, positions['ReactionRule'], 'ReactionRule',
           colors['utility'], box_dims['ReactionRule'], zorder=3)
create_box(ax, positions['ReactionGeometrySet'], 'ReactionGeometry\nSet',
           colors['utility'], box_dims['ReactionGeometrySet'], zorder=3)

# Attributes
attributes = {
    'mol_name': 'name',
    'radius_nm': 'radius_nm',
    'diffusion': 'D_trans\nD_rot',
    'coordinates': 'coordinates',
    'energy': 'energy',
    'binding_state': 'binding\nstate',
    'interfaces_map': 'interfaces\nmap',
    'com_norm': 'com/norm'
}

for key, text in attributes.items():
    create_box(ax, positions[key], text,
               colors['attribute'], box_dims[key], zorder=3)

# Create hierarchical relationships with clear labels
create_labeled_arrow(ax, 'System', 'Units', 'units', color='red')
create_labeled_arrow(ax, 'System', 'workspace_path', 'workspace', color='red')
create_labeled_arrow(ax, 'System', 'pdb_id', 'pdb_id', color='red')
create_labeled_arrow(ax, 'System', 'MoleculeTypeRegistry',
                     'molecule_types', color='red')
create_labeled_arrow(ax, 'System', 'InterfaceTypeRegistry',
                     'interface_types', color='red')
create_labeled_arrow(ax, 'System', 'MoleculeInstanceRegistry',
                     'molecule_instances', color='red')
create_labeled_arrow(ax, 'System', 'InterfaceInstanceRegistry',
                     'interface_instances', color='red')

# Registry to type relationships
create_labeled_arrow(ax, 'MoleculeTypeRegistry',
                     'MoleculeType', 'contains', color='blue')
create_labeled_arrow(ax, 'InterfaceTypeRegistry',
                     'InterfaceType', 'contains', color='blue')
create_labeled_arrow(ax, 'MoleculeInstanceRegistry',
                     'MoleculeInstance', 'contains', color='blue')
create_labeled_arrow(ax, 'InterfaceInstanceRegistry',
                     'InterfaceInstance', 'contains', color='blue')

# Cross-references with very specific labels
# Add intermediate points for complex references to avoid overlaps
# InterfaceType references to MoleculeType
ax.annotate('', xy=get_box_edges(positions['MoleculeType'], box_dims['MoleculeType'])['right'],
            xytext=get_box_edges(
                positions['InterfaceType'], box_dims['InterfaceType'])['left'],
            arrowprops=dict(arrowstyle='<-', color='purple', linestyle='--', lw=2, alpha=0.8))
ax.text(4.5, 8.3, 'this_mol_type', ha='center', va='center', fontsize=8, color='purple', weight='bold',
        bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9, edgecolor='purple'))

# Instance to Type references
ax.annotate('', xy=get_box_edges(positions['MoleculeType'], box_dims['MoleculeType'])['bottom'],
            xytext=get_box_edges(
                positions['MoleculeInstance'], box_dims['MoleculeInstance'])['top'],
            arrowprops=dict(arrowstyle='<-', color='green', linestyle='--', lw=2, alpha=0.8))
ax.text(7.5, 7.25, 'molecule_type', ha='center', va='center', fontsize=8, color='green', weight='bold',
        bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9, edgecolor='green'))

ax.annotate('', xy=get_box_edges(positions['InterfaceType'], box_dims['InterfaceType'])['right'],
            xytext=get_box_edges(
                positions['InterfaceInstance'], box_dims['InterfaceInstance'])['left'],
            arrowprops=dict(arrowstyle='<-', color='green', linestyle='--', lw=2, alpha=0.8))
ax.text(12.5, 7.2, 'interface_type', ha='center', va='center', fontsize=8, color='green', weight='bold',
        bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.9, edgecolor='green'))

# InterfaceInstance to MoleculeInstance
create_labeled_arrow(ax, 'InterfaceInstance', 'MoleculeInstance', 'this_mol',
                     color='green', linestyle='--', style='<-')

# Reaction references
create_labeled_arrow(ax, 'ReactionRule',
                     'ReactionGeometrySet', 'geometry', color='purple')

# Some key attribute relationships (select few to avoid clutter)
create_labeled_arrow(ax, 'MoleculeType', 'mol_name',
                     'name', color='gray', alpha=0.6)
create_labeled_arrow(ax, 'MoleculeType', 'radius_nm',
                     'radius', color='gray', alpha=0.6)
create_labeled_arrow(ax, 'InterfaceType', 'coordinates',
                     'coords', color='gray', alpha=0.6)
create_labeled_arrow(ax, 'InterfaceType', 'energy',
                     'energy', color='gray', alpha=0.6)

# Position legend in empty space (top left)
legend_elements = [
    patches.Patch(color=colors['system'], label='System (Top Level)'),
    patches.Patch(color=colors['config'], label='Configuration'),
    patches.Patch(color=colors['registry'], label='Registries'),
    patches.Patch(color=colors['type'], label='Type Definitions'),
    patches.Patch(color=colors['instance'], label='Runtime Instances'),
    patches.Patch(color=colors['utility'], label='Utility Classes'),
    patches.Patch(color=colors['attribute'], label='Attributes'),
]

ax.legend(handles=legend_elements, loc='upper left',
          bbox_to_anchor=(0, 0.98), fontsize=10)

# Position relationship legend in bottom right empty space
relationship_text = """Relationship Types:
━━━ Containment (red/blue)
┅┅┅ References (purple/green)
← Points to/Uses"""

ax.text(0.75, 0.25, relationship_text, transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor='lightyellow',
                  alpha=0.9, edgecolor='black'),
        verticalalignment='top')

# Set axis properties with more space
ax.set_xlim(-1, 21)
ax.set_ylim(1.5, 15.5)
ax.set_aspect('equal')
ax.axis('off')

# Add title
plt.title('ionerdss System Architecture: Object Hierarchy and Cross-References',
          fontsize=16, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig(
    "./crosslink_diagram.png", dpi=300)
