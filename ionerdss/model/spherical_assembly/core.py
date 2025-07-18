"""
capsid_sphere_pipeline.py

Modularized driver for spherical reconstruction pipeline.
Generalized from HIV Gag to support arbitrary structure mapped to sphere.

TODO:
(1) Check for distance between center of masses and regularize
(2) Check if the subunits are rotationally aligned and rotate them to the correct position if not
(3) Figure out curvature from center of mass and calculate the radius of fitted sphere
    - RMSD based gradient descent fitting?

"""
import os
import numpy as np
import pandas as pd
from .geometry import calculate_gradient, calculate_rmsd


def run_spherical_reconstruction_pipeline(pdb_file, save_dir, options=None):
    """
    Reconstructs and outputs a spherical structure by projecting subunits onto a defined sphere.

    Parameters
    ----------
    pdb_file : str
        Path to the input file containing center and interface positions per subunit.
    save_dir : str
        Directory to save output files.
    options : dict or None
        Optional parameters:
            - 'num_units': int, number of repeating units (default=18)
            - 'points_per_unit': int, number of points per unit (default=6)
            - 'radius': float, target radius for projection (nm)
            - 'fit_sphere': bool, whether to compute best-fit sphere (default=True)
    """
    # Parse input options
    options = options or {}
    num_units = options.get("num_units", 18)
    points_per_unit = options.get("points_per_unit", 6)
    radius = options.get("radius", 25.0)
    fit_sphere = options.get("fit_sphere", True)

    # Step 1: Load and normalize positions
    df = pd.read_csv(pdb_file, header=None, sep=r'\s+')
    df.columns = ['x', 'y', 'z']
    df = df / 10.0  # Convert from Å to nm
    positions_vec = df.to_numpy()

    # Step 2: Determine or fit sphere
    center_vec = np.array([positions_vec[points_per_unit * i] for i in range(num_units)])

    if fit_sphere:
        xyz_r = [0, 0, 0, 70]
        force = calculate_gradient(center_vec, xyz_r)
        rmsd_old = calculate_rmsd(center_vec, xyz_r)
        while np.linalg.norm(force) > 0.01:
            step = 1.0
            while True:
                trial = np.array(xyz_r) - step * force
                rmsd_new = calculate_rmsd(center_vec, trial)
                if rmsd_new < rmsd_old:
                    break
                step *= 0.8
            xyz_r = trial
            rmsd_old = rmsd_new
            force = calculate_gradient(center_vec, xyz_r)
        x0, y0, z0, r0 = xyz_r
        print("Fitted sphere center and radius:", xyz_r)
    else:
        x0, y0, z0 = options.get("center", [0.0, 0.0, 0.0])
        r0 = radius

    # Step 3: Normalize positions to origin + surface
    shift_vec = np.array([x0, y0, z0])
    positions_vec -= shift_vec
    center_vec -= shift_vec
    for i in range(num_units):
        shift = (center_vec[i] / np.linalg.norm(center_vec[i]) * r0) - center_vec[i]
        center_vec[i] += shift
        for j in range(points_per_unit):
            positions_vec[points_per_unit * i + j] += shift

    # Step 4: Output the constructed positions
    os.makedirs(save_dir, exist_ok=True)
    coord_out = os.path.join(save_dir, f"coordR{int(radius)}.txt")
    with open(coord_out, 'w', encoding="utf-8") as f:
        for i in range(num_units):
            unit = positions_vec[points_per_unit * i : points_per_unit * (i + 1)]
            f.write(f"UNIT{i:02d}\n")
            for xyz in unit:
                f.write(f"{xyz[0]:.8f} {xyz[1]:.8f} {xyz[2]:.8f}\n")

    print("Spherical reconstruction coordinates written to:", coord_out)

    return
