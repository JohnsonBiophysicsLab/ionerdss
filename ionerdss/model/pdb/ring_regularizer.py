"""
ionerdss.model.pdb.ring_regularizer

Sphere projection regularization for molecular assemblies.

This module implements a regularization strategy where all molecules in the system
are projected onto concentric spheres. The center of these spheres is determined
by the best-fit sphere of the most abundant molecular species in the system (by stoichiometry).

This approach assumes that the biological assembly is roughly spherical or
can be approximated as layers of spherical shells (e.g., viral capsids).

ALGORITHM:
1. Identify the "reference species": the molecule type with the highest count (stoichiometry).
2. Fit a sphere to the Center of Mass (COM) positions of all instances of the reference species.
   This defines the global center O_0 and the reference radius R_ref.
3. Project all instances of the reference species onto the sphere (O_0, R_ref) along the
   vector from O_0 to their current COM.
4. For each other species S_i:
   a. Calculate the average distance R_i of all instances of S_i from the global center O_0.
   b. Project all instances of S_i onto the sphere (O_0, R_i) along the vector from O_0 to their current COM.

This preserves the relative angular positions of molecules while enforcing perfect spherical symmetry
in terms of radial distance.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
from scipy.optimize import minimize

from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance
from .file_manager import WorkspaceManager


class RingRegularizer:
    """Regularizes system geometry by projecting onto concentric spheres.

    Attributes:
        system: ionerdss System to process
        workspace_manager: Workspace manager for logging
    """

    def __init__(self, system: System, workspace_manager: Optional[WorkspaceManager] = None):
        """Initialize regularizer.

        Args:
            system: ionerdss System to process
            workspace_manager: Workspace manager for logging
        """
        self.system = system
        self.workspace_manager = workspace_manager

    def regularize(self) -> bool:
        """Perform sphere projection regularization.

        Returns:
            True if regularization was performed successfully
        """
        if self.workspace_manager:
            self.workspace_manager.logger.info("Starting Sphere Projection Regularization")

        # Step 1: Identify reference species (highest stoichiometry)
        ref_species_name, counts = self._get_stoichiometry()
        if not ref_species_name:
            if self.workspace_manager:
                self.workspace_manager.logger.warning("No molecules found in system, skipping regularization")
            return False

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Reference species for sphere fitting: %s (count=%d)",
                ref_species_name, counts[ref_species_name]
            )

        # Step 2: Fit sphere to reference species
        ref_molecules = [m for m in self.system.molecule_instances 
                         if m.molecule_type.name == ref_species_name]
        
        # Collect positions
        ref_positions = np.array([m.com for m in ref_molecules])
        
        # Fit sphere
        center, ref_radius, error = self._fit_sphere(ref_positions)
        
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Fitted reference sphere: center=%s, radius=%.3f nm, error=%.4f",
                np.round(center, 3), ref_radius, error
            )

        # Step 3 & 4: Project all species onto their respective shells
        self._project_all_species(center, counts.keys())

        if self.workspace_manager:
            self.workspace_manager.logger.info("Sphere projection completed")

        return True

    def _get_stoichiometry(self) -> Tuple[Optional[str], Dict[str, int]]:
        """Calculate stoichiometry and find the most abundant species.

        Returns:
            Tuple of (most_abundant_species_name, counts_dict)
        """
        counts = {}
        for mol in self.system.molecule_instances:
            m_type = mol.molecule_type.name
            counts[m_type] = counts.get(m_type, 0) + 1
            
        if not counts:
            return None, {}
            
        # Find key with max value
        most_abundant = max(counts, key=counts.get)
        return most_abundant, counts

    def _fit_sphere(self, points: np.ndarray) -> Tuple[np.ndarray, float, float]:
        """Fit a sphere to a set of 3D points.

        Args:
            points: (N, 3) array of points

        Returns:
            Tuple of (center, radius, mean_error)
        """
        # Initial guess: centroid
        centroid = np.mean(points, axis=0)
        
        # Helper to calculate distances from a center
        def _calc_R(c):
            return np.linalg.norm(points - c, axis=1)

        # Objective: minimize variance of radius (which makes points equidistant)
        def _loss(c):
            R = _calc_R(c)
            return np.var(R)

        # Optimization
        res = minimize(_loss, centroid, method='Nelder-Mead', tol=1e-6)
        center = res.x
        
        # Radius is the mean distance to the optimized center
        radii = np.linalg.norm(points - center, axis=1)
        radius = np.mean(radii)
        error = np.std(radii) # Standard deviation is a good measure of error/sphericity
        
        return center, radius, error

    def _project_all_species(self, center_O0: np.ndarray, species_names: List[str]) -> None:
        """Project all species onto concentric spheres centered at O0.

        Args:
            center_O0: The global system center derived from reference species
            species_names: List of all molecule type names to process
        """
        for s_name in species_names:
            # Get all instances of this species
            molecules = [m for m in self.system.molecule_instances 
                         if m.molecule_type.name == s_name]
            
            if not molecules:
                continue

            # Calculate average radius for this species
            positions = np.array([m.com for m in molecules])
            distances = np.linalg.norm(positions - center_O0, axis=1)
            avg_radius = np.mean(distances)
            
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Projecting species %s to sphere radius %.3f nm",
                    s_name, avg_radius
                )

            # Project each molecule
            for i, mol in enumerate(molecules):
                # vector from center to current position
                vec = mol.com - center_O0
                current_dist = np.linalg.norm(vec)
                
                if current_dist < 1e-6:
                    # Point is practically at the center, cannot project distinct direction
                    # Keep it at center or warn?
                    # Ideally this shouldn't happen for shell structures
                     if self.workspace_manager:
                        self.workspace_manager.logger.warning(
                            "Molecule %s is at the system center, cannot project to sphere.", mol.name
                        )
                else:
                    # Normalize direction and scale to avg_radius
                    direction = vec / current_dist
                    new_pos = center_O0 + direction * avg_radius
                    
                    # Update molecule position
                    mol.com = new_pos
