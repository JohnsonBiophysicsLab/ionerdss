"""
visualize.py

Visualization utilities for inspecting coarse-grained molecular assemblies
detected from PDB/mmCIF files.

This module renders:
- Molecule centers (COMs) as spheres
- Interfaces as smaller colored markers
- Interface connections (optional) as lines/arcs
- Labels for molecules and interfaces

Requires `matplotlib` and `mpl_toolkits.mplot3d` for 3D rendering.

Functions
---------
- plot_coarse_grain_model(model, show_reactions=False)
- plot_interfaces(ax, coords, color, label=None)
"""

import matplotlib.pyplot as plt
import numpy as np


def plot_interfaces(ax, coords, color="red", label=None, size=40):
    """
    Plot spherical markers at interface positions.

    Parameters
    ----------
    ax : Axes3D
        The matplotlib 3D axis to draw on.
    coords : list of (3,) ndarray
        Coordinates of interface positions.
    color : str
        Color for marker.
    label : str or None
        Optional label for the group of interfaces.
    size : int
        Marker size.
    """
    coords = np.asarray(coords)
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
               c=color, s=size, label=label)

def plot_coarse_grain_model(model, show_reactions=False, figsize=(10, 8), elev=30, azim=30):
    """
    Render a 3D plot of the coarse-grained model.

    Parameters
    ----------
    model : dict
        The regularized model output from `regularize_model()`, containing:
        - 'COMs', 'interface_coords', 'interfaces', 'chain_ids', etc.
    show_reactions : bool
        If True, draw lines between each interface pair.
    figsize : tuple
        Size of the matplotlib figure.
    elev : float
        Elevation angle for 3D view.
    azim : float
        Azimuthal angle for 3D view.
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    coms = np.array(model["COMs"])
    ax.scatter(coms[:, 0], coms[:, 1], coms[:, 2],
               c="black", s=80, label="COMs")

    # Plot interfaces for each chain
    for i, iface_coords in enumerate(model["interface_coords"]):
        plot_interfaces(ax, iface_coords, color="red")

    # Optionally draw edges between interfaces
    if show_reactions:
        for i, neighbors in enumerate(model["interfaces"]):
            for j, neighbor in enumerate(neighbors):
                if i < neighbor:
                    coord1 = model["interface_coords"][i][j]
                    # Find back-link index
                    try:
                        k = model["interfaces"][neighbor].index(i)
                        coord2 = model["interface_coords"][neighbor][k]
                        x, y, z = zip(coord1, coord2)
                        ax.plot(x, y, z, "gray", linewidth=1)
                    except ValueError:
                        continue

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=elev, azim=azim)
    ax.legend()
    plt.tight_layout()
    plt.show()
