"""
ionerdss.model.pdb.nerdss_exporter
===================================

Export ionerdss System to NERDSS simulation files.

This module converts a complete ionerdss System object into the file format
required by NERDSS simulations, generating `.mol` files for each molecule type
and a `parms.inp` file with reaction parameters and simulation settings.

**CRITICAL UNIT CONVENTIONS:**
    ⚠️ **This module expects ALL input coordinates in nanometers (nm)!** ⚠️
    
    - **Input (from System)**: Nanometers (nm) - already converted by PDB parser
    - **Output (.mol files)**: Nanometers (nm) - NERDSS format
    - **Output (parms.inp)**: Nanometers (nm) - NERDSS format
    - **Angles**: Radians - converted from geometric calculations
    - **Energies**: kT (dimensionless) - binding free energies

**Coordinate Systems:**
    - **Global frame**: All molecule COM and interface coordinates in simulation box (nm)
    - **Local frame**: Interface positions relative to molecule COM (nm)
    - **Reference vectors**: Unit vectors defining molecule orientation in global frame

**Key Methods:**
    - ``export()``: Main export method - generates all NERDSS files
    - ``_write_mol_file()``: Creates .mol file for each molecule type
    - ``_write_parms_file()``: Creates parms.inp with all reactions
    - ``_generate_reaction_angles()``: Calculates binding angles (θ, φ, ω)

**File Outputs:**
    - ``<MoleculeType>.mol``: Molecule geometry and binding sites (nm)
    - ``parms.inp``: Simulation parameters and reaction definitions

**Homotypic Parameter Sharing:**
    **Strategy**: Calculate normal vectors automatically, then share both sigma and
    angles between all homotypic interfaces since they use the same normal vector.

    **Cache Key**: Use representative site for both sigma and angles since normal
    vectors are now consistent within homotypic groups.

**Angle Calculation:**
    Uses reference frame from each molecule instance to calculate NERDSS-compatible
    binding angles. See ``_generate_reaction_angles()`` for mathematical details.

.. warning::
    **Unit Mismatch Will Cause Simulation Failure!**
    If coordinates are accidentally in Ångström instead of nm, NERDSS will
    produce incorrect results or crash. Always verify units before export.

See Also:
    - ``ionerdss.model.pdb.parser``: Converts PDB (Å) → System (nm)
    - NERDSS documentation: https://nerdss.github.io/
"""

import re
from typing import Dict, Any, Iterable, Optional, List, Tuple
from dataclasses import dataclass
from pathlib import Path
import math

import numpy as np

from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.utils.vectors import convert_to_unit, get_magnitude
from ionerdss.utils.angles import signed_angle_arccos
from ionerdss.model.pdb import interface_naming
from .file_manager import WorkspaceManager


#------------ exporter class ---------------

class NERDSSExporter:
    """Exporter for converting ionerdss System to NERDSS simulation files.

    Generates .mol files for each molecule type and parms.inp with simulation
    parameters and reaction definitions. Automatically calculates normal vectors
    based on interface geometry.

    Attributes:
        system: ionerdss System to export.
        workspace_manager: Workspace manager for file organization.
        output_dir: Directory for NERDSS output files.
        interface_to_site_map: Mapping from interface names to site labels.
        reaction_metadata: Metadata about reactions for rate calculation.
        homotypic_interface_map: Maps site labels to representative labels.
        calculated_normals: Automatically calculated normal vectors per molecule type.
        reaction_params_cache: Cache for calculated reaction parameters.
    """

    # --------------- constructor ------------

    def __init__(self, system: System, workspace_manager: Optional[WorkspaceManager] = None):
        """Initialize NERDSS exporter."""
        self.system = system
        self.workspace_manager = workspace_manager

        # Store mapping from interface names to site labels
        self.interface_to_site_map: Dict[str, str] = {}

        # Store reaction metadata for rate calculation
        self.reaction_metadata: List[Dict[str, Any]] = []

        # Cache for homotypic interface mapping: site_label -> representative_site_label
        self.homotypic_interface_map: Dict[str, str] = {}

        # Store automatically calculated normal vectors: mol_type -> {site_label: normal_vector}
        self.calculated_normals: Dict[str, Dict[str, np.ndarray]] = {}

        # Cache for reaction parameters: (representative_site1, representative_site2) -> (sigma, angles)
        self.reaction_params_cache: Dict[Tuple[str, str],
                                         Tuple[float, Tuple[float, float, float, float, float]]] = {}

        # Precalculated geometry map: (mol1, type1, mol2, type2) -> (sigma, angles)
        # Allows manually injecting exact parameters (e.g. from Platonic Solids model)
        # preventing re-measurement from structures.
        self.precalculated_geometry: Dict[Tuple[str, str, str, str], 
                                          Tuple[float, Tuple[float, float, float, float, float]]] = {}

        # Precalculated rates map: (mol1, site1, mol2, site2) -> (ka, kb)
        # Allows manually injecting kinetic parameters
        self.precalculated_rates: Dict[Tuple[str, str, str, str], Tuple[float, float]] = {}

        # Create NERDSS output directory in workspace
        if workspace_manager:
            self.output_dir = workspace_manager.workspace_path / 'nerdss_files'
            self.output_dir.mkdir(exist_ok=True)
            self.workspace_manager.logger.info(
                "Created NERDSS export directory: %s", self.output_dir)
        else:
            self.output_dir = Path("nerdss_files")
            self.output_dir.mkdir(exist_ok=True)
            
    # ------------- helpers -----------------

    def _circular_mean_signed(self, vals: List[float]) -> float:
        """Mean on [-π, π] using atan2(⟨sin⟩, ⟨cos⟩)."""
        if not vals:
            return 0.0
        s = np.mean(np.sin(vals))
        c = np.mean(np.cos(vals))
        return float(np.arctan2(s, c))

    def _circular_var_signed(self, vals: List[float]) -> float:
        """Circular variance on [-π, π]; 0 → tight, 1 → uniform."""
        if not vals:
            return 1.0
        s = np.mean(np.sin(vals))
        c = np.mean(np.cos(vals))
        R = np.hypot(s, c)
        return float(1.0 - R)

    def _enumerate_binding_instance_pairs(self, mol1_name: str, site1: str,
                                      mol2_name: str, site2: str):
        """
        Return only (m1, m2, intf1, intf2) where intf1.type ∈ types(site1) AND intf2.type ∈ types(site2),
        and intf1↔intf2 are mutual partners (actual binding pair).
        """
        site1_types = self._site_to_interface_types(site1)
        site2_types = self._site_to_interface_types(site2)
        pairs = []

        for m1 in self.system.molecule_instances:
            if not (m1.molecule_type and m1.molecule_type.name == mol1_name):
                continue
            for intf1, m2 in m1.interfaces_neighbors_map.items():
                if not (m2 and m2.molecule_type and m2.molecule_type.name == mol2_name):
                    continue
                if self._iface_type_name(intf1) not in site1_types:
                    continue

                # find partner interface that (a) points back to m1 and (b) is one of the site2 types
                for intf2, back in m2.interfaces_neighbors_map.items():
                    if back is m1 and self._iface_type_name(intf2) in site2_types:
                        pairs.append((m1, m2, intf1, intf2))
                        break  # one partner per intf1 is enough
        return pairs

    def _enumerate_exact_type_pairs(self, mol1_name: str, type1: str,
                                mol2_name: str, type2: str):
        """
        Find all (m1, m2, intf1, intf2) where:
        - m1 is mol1_name and has an interface intf1 with type == type1
        - m2 is mol2_name and has an interface intf2 with type == type2
        - intf1's partner is m2 and intf2's partner is m1 (actual bound pair)
        """
        out = []
        for m1 in self.system.molecule_instances:
            if not (m1.molecule_type and m1.molecule_type.name == mol1_name):
                continue
            for intf1, m2 in m1.interfaces_neighbors_map.items():
                if not (intf1 and m2 and m2.molecule_type and m2.molecule_type.name == mol2_name):
                    continue
                if self._iface_type_name(intf1) != type1:
                    continue
                # find partner interface on m2 that points back and has the exact requested type
                for intf2, back in m2.interfaces_neighbors_map.items():
                    if back is m1 and self._iface_type_name(intf2) == type2:
                        out.append((m1, m2, intf1, intf2))
                        break
        return out


    def _mean_params_from_pairs(self, pairs: List[Tuple[Any, Any, Any, Any]],
                            mol1: str, site1: str, mol2: str, site2: str
                        ) -> Tuple[float, Tuple[float,float,float,float,float]]:
        """
        Compute averaged σ and angles across all provided binding pairs.
        θ's get arithmetic mean/std; φ/ω use circular mean/std.
        Also prints a per-pair debug line including local and rotated normals.
        """
        sigmas = []
        th1s, th2s = [], []
        ph1s, ph2s, ws = [], [], []

        self.workspace_manager.logger.info("\n=== DEBUG: Averaging pairs for "
            f"{mol1}({site1}) + {mol2}({site2}) ===")

        header = ("idx  |  sigma   "
                "|  n1(local)              n2(local)              "
                "|  n1f(global)              n2f(global)             "
                "|  theta1   theta2   phi1      phi2      omega")
        self.workspace_manager.logger.info(header)
        self.workspace_manager.logger.info("-" * len(header))

        for k, (m1, m2, intf1, intf2) in enumerate(pairs, start=1):
            com1, com2 = m1.com, m2.com
            p1, p2 = intf1.absolute_coord, intf2.absolute_coord

            # --- normals ---
            # local normals are what will be written to parms.inp (NERDSS expects local)
            n1_local = self._get_species_normal_vector(mol1, site1, m1)
            n2_local = self._get_species_normal_vector(mol2, site2, m2)
            # rotated/global normals used in angle calc (n1f, n2f)
            n1f = self._get_species_normal_vector_with_rotation(mol1, site1, m1)
            n2f = self._get_species_normal_vector_with_rotation(mol2, site2, m2)

            # --- angles ---
            sigma, (th1, th2, ph1, ph2, w) = self._generate_reaction_angles(
                p1, p2, com1, com2, mol1, mol2, site1, site2
            )

            sigmas.append(sigma)
            th1s.append(th1); th2s.append(th2)
            ph1s.append(ph1); ph2s.append(ph2); ws.append(w)

            def vfmt(v):
                return f"[{v[0]: .6f},{v[1]: .6f},{v[2]: .6f}]"

            self.workspace_manager.logger.info(f"{k:>3d}  | {sigma: .6f} "
                f"| {vfmt(n1_local)} {vfmt(n2_local)} "
                f"| {vfmt(n1f)} {vfmt(n2f)} "
                f"| {th1: .6f} {th2: .6f} {ph1: .6f} {ph2: .6f} {w: .6f}")

        # --- summary stats ---
        if sigmas:
            sigma_mean = float(np.mean(sigmas))
            sigma_std  = float(np.std(sigmas))
        else:
            sigma_mean = 1.0; sigma_std = 0.0

        if th1s:
            theta1_mean = float(np.mean(th1s)); theta1_std = float(np.std(th1s))
        else:
            theta1_mean = 0.0; theta1_std = 0.0
        if th2s:
            theta2_mean = float(np.mean(th2s)); theta2_std = float(np.std(th2s))
        else:
            theta2_mean = 0.0; theta2_std = 0.0

        phi1_mean, phi1_std = self._circular_mean_std(ph1s) if ph1s else (0.0, 0.0)
        phi2_mean, phi2_std = self._circular_mean_std(ph2s) if ph2s else (0.0, 0.0)
        w_mean, w_std       = self._circular_mean_std(ws)   if ws   else (0.0, 0.0)

        # Return the means that will be used downstream
        return sigma_mean, (theta1_mean, theta2_mean, phi1_mean, phi2_mean, w_mean)

    
    def _site_to_interface_types(self, site_label: str) -> set[str]:
        """
        Return the exact InterfaceType names that this site_label represents.
        Built from interface_to_site_map (the mapping you create in _write_mol_file).
        """
        types = set()
        for iface_type_name, mapped_site in self.interface_to_site_map.items():
            if mapped_site == site_label:
                types.add(iface_type_name)
        return types

    def _iface_type_name(self, interface_instance) -> str:
        # Consistent accessor
        return interface_instance.interface_type.get_name() if interface_instance.interface_type else ""

    # ---------- main export function --------


    def export_all(self, molecule_counts: Optional[Dict[str, int]] = None,
                   box_nm: Tuple[float, float, float] = (100.0, 100.0, 100.0),
                   parms_overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Path]:
        """Export complete NERDSS simulation setup."""
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Exporting NERDSS simulation files...")

        output_files = {}

        # Set default molecule counts
        if molecule_counts is None:
            molecule_counts = {}
            for mol_type in self.system.molecule_types:
                molecule_counts[mol_type.name] = 10

        # Clear mappings and caches for fresh export
        self.interface_to_site_map.clear()
        self.reaction_metadata.clear()
        self.homotypic_interface_map.clear()
        self.homotypic_interface_map.clear()
        self.calculated_normals.clear()
        self.reaction_params_cache.clear()
        # Note: We DO NOT clear self.precalculated_geometry as it is user-provided configuration

        # Export .mol files for each molecule type (this builds the mapping)
        for mol_type in self.system.molecule_types:
            mol_file_path = self._write_mol_file(mol_type, parms_overrides.get('hyperparams') if parms_overrides else None)
            output_files[f"{mol_type.name}_mol"] = mol_file_path

        # after the loop that calls _write_mol_file(...) for all mol types
        # NOTE: legacy 1↔2 homodimer validation removed.
        # With f/b scheme, validation happens implicitly when generating reactions below.


        # Calculate normal vectors for each molecule type
        self._calculate_normal_vectors()

        # Generate reactions using the stored mapping
        reactions = self._generate_reactions()

        # Calculate bond lengths and angles for reactions
        sigma_list, angles_list = self._calculate_reaction_parameters(
            reactions)

        # Export parms.inp
        parms_path = self._write_parms_file(
            reactions=reactions,
            molecule_counts=molecule_counts,
            box_nm=box_nm,
            sigma_list=sigma_list,
            angles_list=angles_list,
            parms_overrides=parms_overrides
        )
        output_files['parms'] = parms_path

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Exported %d NERDSS files", len(output_files))
            for file_type, file_path in output_files.items():
                self.workspace_manager.logger.info(
                    "  %s: %s", file_type, file_path)
                
        # DEBUG: Print representative instance information
        for mol_type in self.system.molecule_types:
            self._debug_representative_instance(mol_type.name)

        return output_files

    def _calculate_normal_vectors(self):
        """Calculate normal vectors for each molecule instance based on its specific interface geometry."""
        # Clear the storage - now we store per instance, not per type
        self.calculated_normals.clear()

        # Calculate normal vectors for each molecule instance individually
        for mol_instance in self.system.molecule_instances:
            if not mol_instance.molecule_type:
                continue

            mol_type_name = mol_instance.molecule_type.name
            # Use instance ID as unique identifier
            instance_id = id(mol_instance)

            # Initialize storage for this instance
            if mol_type_name not in self.calculated_normals:
                self.calculated_normals[mol_type_name] = {}
            self.calculated_normals[mol_type_name][instance_id] = {}

            # Group interfaces by their type for this specific instance
            homotypic_groups = self._group_interfaces_by_type(mol_instance)

            for interface_type_name, interface_group in homotypic_groups.items():
                # Get site labels and coordinates for this specific instance
                site_labels = []
                interface_coords = []

                for interface_data in interface_group:
                    # Find corresponding site label
                    for key, site_label in self.interface_to_site_map.items():
                        if (key.startswith(interface_type_name) or key == interface_type_name):
                            if site_label not in site_labels:
                                site_labels.append(site_label)
                                interface_coords.append(
                                    interface_data['coord'])
                                break

                if not site_labels:
                    continue

                # Calculate normal vector for this instance's interfaces
                if len(site_labels) == 1:
                    # Singleton interface
                    normal_vector = self._calculate_singleton_normal_vector(
                        interface_coords[0], mol_type_name, interface_type_name
                    )
                else:
                    # Multiple interfaces - no rotation correction needed since we're using actual instance
                    normal_vector = self._calculate_group_normal_vector(
                        interface_coords, mol_type_name, interface_type_name
                    )

                # Assign normal vector to all sites in this group for this instance
                for site_label in site_labels:
                    self.calculated_normals[mol_type_name][instance_id][site_label] = normal_vector

                if self.workspace_manager:
                    group_type = "singleton" if len(
                        site_labels) == 1 else "group"
                    self.workspace_manager.logger.info(
                        "Calculated normal vector for instance %s of %s %s interface %s (sites %s): %s",
                        instance_id, mol_type_name, group_type, interface_type_name,
                        site_labels, normal_vector
                    )

    def _get_species_normal_vector(self, mol_name: str, site_label: str,
                                   mol_instance=None) -> np.ndarray:
        """Get calculated normal vector for a specific molecule instance and site.

        Args:
            mol_name: Molecule type name (e.g., 'A').
            site_label: Site label (e.g., 'a1', 'a2').
            mol_instance: Specific molecule instance (if None, tries to find from context).

        Returns:
            Calculated normal vector for this specific instance.
        """
        if mol_instance is None:
            # Try to find the instance from current context (fallback behavior)
            mol_instance = self._find_molecule_instance_for_site(
                mol_name, site_label)

        if mol_instance is not None:
            instance_id = id(mol_instance)

            if (mol_name in self.calculated_normals and
                instance_id in self.calculated_normals[mol_name] and
                    site_label in self.calculated_normals[mol_name][instance_id]):
                return self.calculated_normals[mol_name][instance_id][site_label]

        # Fallback: try to find any instance with this site
        if mol_name in self.calculated_normals:
            for instance_id, sites in self.calculated_normals[mol_name].items():
                if site_label in sites:
                    if self.workspace_manager:
                        self.workspace_manager.logger.info(
                            "Using fallback normal vector for %s.%s from instance %s",
                            mol_name, site_label, instance_id
                        )
                    return sites[site_label]

        # Final fallback to default
        if self.workspace_manager:
            self.workspace_manager.logger.warning(
                "No calculated normal for %s.%s (instance %s), using default [1,0,0]",
                mol_name, site_label, id(
                    mol_instance) if mol_instance else "unknown"
            )
        return np.array([1.0, 0.0, 0.0])

    def _get_species_normal_vector_with_rotation(self, mol_name: str, site_label: str, 
                                            mol_instance=None) -> np.ndarray:
        """Get normal vector for a specific molecule instance with rotation correction."""
        if mol_instance is None:
            return self._get_species_normal_vector(mol_name, site_label, mol_instance)
        
        # Get the base normal vector from representative structure
        base_normal = self._get_base_normal_vector_for_site(mol_name, site_label)
        
        # Check if this instance is the representative
        representative_instance = self._get_representative_instance(mol_name)
        
        if representative_instance is None or mol_instance == representative_instance:
            # This is the representative instance, no rotation needed
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Using unrotated normal for representative instance %s.%s (id=%s): %s",
                    mol_name, site_label, id(mol_instance), base_normal
                )
            return base_normal
        
        # This is not the representative, need to apply rotation
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Calculating rotation for non-representative instance %s.%s (id=%s)",
                mol_name, site_label, id(mol_instance)
            )
        
        # Remove the extra site_label parameter
        rotation_matrix = self._calculate_rotation_from_representative(
            mol_name, representative_instance, mol_instance
        )
        
        if rotation_matrix is not None:
            # Apply rotation to the normal vector
            rotated_normal = np.dot(rotation_matrix, base_normal)
            
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Applied rotation to normal for instance %s.%s (id=%s): %s -> %s",
                    mol_name, site_label, id(mol_instance), base_normal, rotated_normal
                )
            
            return rotated_normal
        else:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Could not calculate rotation for instance %s.%s (id=%s), using base normal",
                    mol_name, site_label, id(mol_instance)
                )
            return base_normal

    def _calculate_rotation_from_representative(self, mol_name: str, representative_instance, target_instance) -> Optional[np.ndarray]:
        """
        Proper rotation mapping rep → target using LOCAL interface vectors
        matched by exact type name. Works with 1, 2, or ≥3 matches (see _rotation_from_vectors).
        """
        ref_vecs, tgt_vecs = [], []
        
        def dump_map(inst):
            rows = []
            for intf, partner in inst.interfaces_neighbors_map.items():
                rows.append((intf.interface_type.get_name(), partner.name, tuple(intf.absolute_coord)))
            rows.sort()
            if self.workspace_manager:
                self.workspace_manager.logger.info("Interface map for %s: %s", mol_name, rows)
            # Count by family and f/b:
            from collections import Counter
            fam = [name[:-1] if name[-1] in ("f","b") else name for (name,_,_) in rows]
            ends = [name[-1] if name[-1] in ("f","b") else "-" for (name,_,_) in rows]
            if self.workspace_manager:
                self.workspace_manager.logger.info("By family:", Counter(fam))
                self.workspace_manager.logger.info("Ends f/b:", Counter(ends))
            
        dump_map(representative_instance)
        
        for ref_intf, _ in representative_instance.interfaces_neighbors_map.items():
            tname = ref_intf.interface_type.get_name()
            ref_local = ref_intf.absolute_coord - representative_instance.com
            match = None
            for tgt_intf, _ in target_instance.interfaces_neighbors_map.items():
                if tgt_intf.interface_type.get_name() == tname:
                    match = tgt_intf
                    break
            if match is None:
                continue
            tgt_local = match.absolute_coord - target_instance.com
            ref_vecs.append(ref_local); tgt_vecs.append(tgt_local)

        if not ref_vecs:
            if self.workspace_manager:
                self.workspace_manager.logger.warning("No matched interfaces to define rotation for %s", mol_name)
            return None

        R = self._calculate_kabsch_rotation(ref_vecs, tgt_vecs)
        if R is not None and self._is_valid_rotation_matrix(R):
            return R
        if self.workspace_manager:
            self.workspace_manager.logger.warning("Computed rotation invalid; falling back to single-vector alignment")
        # last fallback: align first vector only
        return self._calculate_single_vector_rotation(ref_vecs[0], tgt_vecs[0])


    def _calculate_single_vector_rotation(self, ref_vec: np.ndarray, target_vec: np.ndarray) -> np.ndarray:
        """
        Kept for compatibility; now delegates to the triad builder to guarantee det=+1.
        """
        def _orthonormal_triad(v1: np.ndarray) -> np.ndarray:
            v1 = np.asarray(v1, float)
            n = np.linalg.norm(v1)
            e1 = v1 / n if n > 1e-12 else np.array([1.0, 0.0, 0.0])
            cand = np.array([1.0, 0.0, 0.0]) if abs(e1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            v2p = cand - np.dot(cand, e1) * e1
            e2 = v2p / (np.linalg.norm(v2p) + 1e-12)
            e3 = np.cross(e1, e2); e3 /= (np.linalg.norm(e3) + 1e-12)
            if np.linalg.det(np.column_stack([e1, e2, e3])) < 0:
                e3 = -e3
            return np.column_stack([e1, e2, e3])

        B_ref = _orthonormal_triad(ref_vec)
        B_tgt = _orthonormal_triad(target_vec)
        R = B_tgt @ B_ref.T
        if np.linalg.det(R) < 0:
            B_tgt[:, 2] *= -1
            R = B_tgt @ B_ref.T
        return R


    def _calculate_axis_alignment_rotation(self, axis_ref: np.ndarray, axis_target: np.ndarray) -> np.ndarray:
        """Calculate rotation to align two axes (collinear interfaces case)."""
        
        # This is the same as single vector rotation
        return self._calculate_single_vector_rotation(axis_ref, axis_target)

    def _rodrigues_rotation(self, axis: np.ndarray, angle: float) -> np.ndarray:
        """Calculate rotation matrix using Rodrigues' rotation formula."""
        
        axis = axis / np.linalg.norm(axis)  # Ensure unit vector
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        
        # Rodrigues' formula
        K = np.array([[0, -axis[2], axis[1]],
                    [axis[2], 0, -axis[0]],
                    [-axis[1], axis[0], 0]])
        
        R = np.eye(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)
        
        return R

    def _calculate_kabsch_rotation(self, ref_vecs: list[np.ndarray], tgt_vecs: list[np.ndarray]) -> Optional[np.ndarray]:
        """
        Return proper rotation R (det=+1) mapping representative → target given
        matched LOCAL vectors (each defined as interface.absolute_coord - COM).

        - ≥3 vectors: Kabsch on centered sets, then fix handedness.
        - 2 vectors : build right-handed bases with Gram–Schmidt (twist fixed by v2⊥).
        - 1 vector  : build bases from v and a deterministic pseudo-perp (twist fixed deterministically).
        """
        ref_vecs = [np.asarray(v, float) for v in ref_vecs]
        tgt_vecs = [np.asarray(v, float) for v in tgt_vecs]
        n = min(len(ref_vecs), len(tgt_vecs))

        def _unit(x):
            nrm = np.linalg.norm(x)
            return x / nrm if nrm > 1e-12 else np.array([1.0, 0.0, 0.0], float)

        def _basis_from(v1, v2=None):
            e1 = _unit(v1)
            if v2 is None:
                # deterministic pseudo-perp
                cand = np.array([1.0, 0.0, 0.0]) if abs(e1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
                v2p = cand - e1 * np.dot(cand, e1)
            else:
                v2p = v2 - e1 * np.dot(v2, e1)
                if np.linalg.norm(v2p) < 1e-12:
                    # nearly collinear: pick a stable pseudo-perp
                    cand = np.array([1.0, 0.0, 0.0]) if abs(e1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
                    v2p = cand - e1 * np.dot(cand, e1)
            e2 = _unit(v2p)
            e3 = np.cross(e1, e2)
            e3n = np.linalg.norm(e3)
            if e3n < 1e-12:
                # last-ditch: nudge e2
                e2 = _unit(e2 + 1e-6*np.array([0.3,0.5,0.2]))
                e3 = np.cross(e1, e2)
                e3n = np.linalg.norm(e3)
            e3 /= max(e3n, 1e-12)
            B = np.column_stack([e1, e2, e3])
            # enforce right-handed basis
            if np.linalg.det(B) < 0:
                B[:, 2] *= -1
            return B

        if n >= 3:
            ref = np.stack(ref_vecs[:n], 0); tgt = np.stack(tgt_vecs[:n], 0)
            ref_c = ref - ref.mean(0); tgt_c = tgt - tgt.mean(0)
            H = ref_c.T @ tgt_c
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            return R

        if n == 2:
            B_ref = _basis_from(ref_vecs[0], ref_vecs[1])
            B_tgt = _basis_from(tgt_vecs[0], tgt_vecs[1])
            R = B_tgt @ B_ref.T
            if np.linalg.det(R) < 0:
                B_tgt[:, 2] *= -1
                R = B_tgt @ B_ref.T
            return R

        if n == 1:
            B_ref = _basis_from(ref_vecs[0], None)
            B_tgt = _basis_from(tgt_vecs[0], None)
            R = B_tgt @ B_ref.T
            if np.linalg.det(R) < 0:
                B_tgt[:, 2] *= -1
                R = B_tgt @ B_ref.T
            return R

        return None


    def _get_base_normal_vector_for_site(self, mol_name: str, site_label: str) -> np.ndarray:
        """Get the base normal vector for a site from the representative structure.

        Args:
            mol_name: Molecule type name.
            site_label: Site label.

        Returns:
            Base normal vector from representative structure.
        """
        # Get representative instance
        representative_instance = self._get_representative_instance(mol_name)
        if representative_instance:
            instance_id = id(representative_instance)

            if (mol_name in self.calculated_normals and
                instance_id in self.calculated_normals[mol_name] and
                    site_label in self.calculated_normals[mol_name][instance_id]):
                return self.calculated_normals[mol_name][instance_id][site_label]

        # Fallback to default
        return np.array([1.0, 0.0, 0.0])

    def _are_complementary_interfaces(self, interface_type1, interface_type2) -> bool:
        """Check if two interface types are complementary (for homodimeric heterotypic cases).

        Args:
            interface_type1: First interface type.
            interface_type2: Second interface type.

        Returns:
            True if they are complementary interfaces.
        """
        # For failed homotypic cases, check if they are complementary pairs
        # A_A_1 is complementary to A_A_2, A_A_3 is complementary to A_A_4, etc.

        try:
            return interface_naming.are_complementary_homodimeric_heterotypic(interface_type1.get_name(), interface_type2.get_name())
        except Exception:
            return False

    def _is_valid_rotation_matrix(self, R: np.ndarray, tolerance: float = 1e-6) -> bool:
        """Check if matrix is a valid rotation matrix."""
        try:
            if R.shape != (3, 3):
                return False

            # Check if R^T * R = I (orthogonal)
            should_be_identity = np.dot(R.T, R)
            identity = np.eye(3)
            if not np.allclose(should_be_identity, identity, atol=tolerance):
                return False

            # Check if det(R) = 1 (proper rotation, not reflection)
            det_R = np.linalg.det(R)
            if not np.isclose(det_R, 1.0, atol=tolerance):
                return False

            return True

        except Exception:
            return False

    def _find_molecule_instance_for_site(self, mol_name: str, site_label: str):
        """Find molecule instance that has the given site label.

        Args:
            mol_name: Molecule type name.
            site_label: Site label to find.

        Returns:
            Molecule instance or None if not found.
        """
        for mol_instance in self.system.molecule_instances:
            if (mol_instance.molecule_type and
                    mol_instance.molecule_type.name == mol_name):

                # Check if this instance has the site label
                instance_id = id(mol_instance)
                if (mol_name in self.calculated_normals and
                    instance_id in self.calculated_normals[mol_name] and
                        site_label in self.calculated_normals[mol_name][instance_id]):
                    return mol_instance

        return None

    def _calculate_group_normal_vector(self, interface_coords: List[np.ndarray],
                                       mol_name: str, interface_type_name: str) -> np.ndarray:
        """Calculate normal vector for multiple interfaces of the same type."""
        n_interfaces = len(interface_coords)

        if n_interfaces == 0:
            return np.array([1.0, 0.0, 0.0])  # Default

        elif n_interfaces == 1:
            # This shouldn't happen here since singletons are handled separately
            # But include for safety
            return self._calculate_singleton_normal_vector(interface_coords[0], mol_name, interface_type_name)

        elif n_interfaces == 2:
            # Two interfaces: use cross product
            v1 = interface_coords[0]
            v2 = interface_coords[1]

            cross_prod = np.cross(v1, v2)
            if np.linalg.norm(cross_prod) > 1e-10:
                normal = cross_prod / np.linalg.norm(cross_prod)
            else:
                # Vectors are parallel, use perpendicular direction
                normal = self._get_perpendicular_vector(v1)

            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Two interfaces %s: cross product %s", interface_type_name, normal
                )
            return normal

        else:
            # Three or more interfaces: fit plane
            return self._fit_plane_normal(interface_coords, mol_name, interface_type_name)

    def _calculate_singleton_normal_vector(self, interface_coord: np.ndarray,
                                           mol_name: str, interface_type_name: str) -> np.ndarray:
        """Calculate normal vector for a singleton interface using simple default method.

        Uses the simple default method:
        1. Try [1, 0, 0] as default
        2. If approximately collinear with COM-to-site vector, use [0, 0, -1]

        Args:
            interface_coord: Interface coordinate relative to COM.
            mol_name: Molecule type name.
            interface_type_name: Interface type name.

        Returns:
            Normal vector for the singleton interface.
        """
        # Default normal vector - Changed back to Z-axis for parms.inp consistency
        default_normal = np.array([0.0, 0.0, 1.0])

        # COM-to-site vector (interface_coord is already relative to COM in local coordinates)
        site_vector = interface_coord
        site_vector_norm = np.linalg.norm(site_vector)

        if site_vector_norm < 1e-10:
            # Interface is at COM, use default
            normal = default_normal
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Singleton interface %s.%s at COM, using default normal %s",
                    mol_name, interface_type_name, normal
                )
        else:
            # Normalize site vector
            site_unit = site_vector / site_vector_norm

            # Check if default normal is approximately collinear with site vector
            dot_product = abs(np.dot(default_normal, site_unit))
            collinearity_threshold = 0.9  # cos(~25 degrees)

            if dot_product > collinearity_threshold:
                # Too collinear, use alternative normal
                normal = np.array([1.0, 0.0, 0.0])
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "Singleton interface %s.%s collinear with default (dot=%.3f), using alternative normal %s",
                        mol_name, interface_type_name, dot_product, normal
                    )
            else:
                # Use default normal
                normal = default_normal
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "Singleton interface %s.%s using default normal %s (dot=%.3f)",
                        mol_name, interface_type_name, normal, dot_product
                    )

        return normal

    def _fit_plane_normal(self, interface_coords: List[np.ndarray],
                          mol_name: str, interface_type_name: str) -> np.ndarray:
        """Fit plane to interface coordinates and return normal vector."""
        coords_array = np.array(interface_coords)  # Shape: (n, 3)

        # Center the coordinates
        centroid = np.mean(coords_array, axis=0)
        centered_coords = coords_array - centroid

        # Use SVD to find the best-fit plane
        try:
            U, s, Vt = np.linalg.svd(centered_coords, full_matrices=False)

            # The normal vector is the last column of V (last row of Vt)
            normal = Vt[-1, :]

            # Calculate fitting residual
            residuals = np.abs(np.dot(centered_coords, normal))
            max_residual = np.max(residuals)
            rms_residual = np.sqrt(np.mean(residuals**2))

            # Check if fit is reasonable
            max_coord_range = np.max(np.ptp(coords_array, axis=0))
            relative_error = rms_residual / max_coord_range if max_coord_range > 0 else 0

            if relative_error > 0.1:  # 10% relative error threshold
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Poor plane fit for %s.%s: RMS residual=%.3f, max=%.3f, relative=%.1f%%",
                        mol_name, interface_type_name, rms_residual, max_residual, relative_error*100
                    )

            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Plane fit for %s.%s: normal=%s, RMS residual=%.3f",
                    mol_name, interface_type_name, normal, rms_residual
                )

            return normal / np.linalg.norm(normal)

        except np.linalg.LinAlgError:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "SVD failed for %s.%s, using default normal", mol_name, interface_type_name
                )
            return np.array([0.0, 0.0, 1.0])

    def _get_perpendicular_vector(self, v: np.ndarray) -> np.ndarray:
        """Get a vector perpendicular to the input vector."""
        v = v / \
            np.linalg.norm(v) if np.linalg.norm(
                v) > 1e-10 else np.array([1.0, 0.0, 0.0])

        # Find the component with smallest absolute value
        min_idx = np.argmin(np.abs(v))

        # Create perpendicular vector
        perp = np.zeros(3)
        perp[min_idx] = 1.0

        # Make it orthogonal using Gram-Schmidt
        perp = perp - np.dot(perp, v) * v
        return perp / np.linalg.norm(perp)

    def _get_representative_instance(self, mol_type_name: str):
        """Get representative instance with maximum interfaces for a molecule type."""
        representative_instance = None
        max_interfaces = 0

        for mol_instance in self.system.molecule_instances:
            if mol_instance.molecule_type and mol_instance.molecule_type.name == mol_type_name:
                interface_count = len(mol_instance.interfaces_neighbors_map)
                if interface_count > max_interfaces:
                    max_interfaces = interface_count
                    representative_instance = mol_instance

        return representative_instance

    def _group_interfaces_by_type(self, mol_instance):
        """Group interfaces by their type name."""
        interface_groups = {}

        for interface_instance, partner_instance in mol_instance.interfaces_neighbors_map.items():
            if interface_instance.interface_type:
                type_name = interface_instance.interface_type.get_name()

                if type_name not in interface_groups:
                    interface_groups[type_name] = []

                interface_groups[type_name].append({
                    'instance': interface_instance,
                    'coord': interface_instance.interface_type.local_coord,
                    'partner': partner_instance.molecule_type.name if (partner_instance and partner_instance.molecule_type) else "unknown",
                    'type_name': type_name
                })

        return interface_groups

    def _validate_site_labels(self, all_site_labels: List[str]) -> bool:
        """Validate that all site labels are unique and unambiguous.

        Args:
            all_site_labels: List of all generated site labels.

        Returns:
            True if all labels are unique.
        """
        unique_labels = set(all_site_labels)

        if len(unique_labels) != len(all_site_labels):
            # Find duplicates
            seen = set()
            duplicates = set()
            for label in all_site_labels:
                if label in seen:
                    duplicates.add(label)
                seen.add(label)

            if self.workspace_manager:
                self.workspace_manager.logger.error(
                    "Duplicate site labels detected: %s", list(duplicates)
                )
            return False

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "All %d site labels are unique: %s",
                len(unique_labels), sorted(unique_labels)
            )

        return True

    def _write_mol_file(self, mol_type: MoleculeType, hyperparams=None, dry_run: bool = False) -> Optional[Path]:
        mol_file_path = self.output_dir / f"{mol_type.name}.mol"

        # Get the representative instance for this molecule type
        rep_inst = self._get_representative_instance(mol_type.name)
        if rep_inst is None:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "No representative instance found for molecule type %s", mol_type.name
                )
            return mol_file_path
        # Initialize dictionaries to collect interface data
        per_type_local: dict[str, np.ndarray] = {}
        per_type_partner_ids: dict[str, int] = {}

        # Get ALL interface types for this molecule type
        mol_interface_types = [it for it in self.system.interface_types if it.this_mol_type_name == mol_type.name]
        
        # Get all instances of this molecule type
        all_instances = [inst for inst in self.system.molecule_instances if inst.molecule_type and inst.molecule_type.name == mol_type.name]

        # Priority 1: Get coordinates from representative instance (Ensures consistent reference frame)
        if rep_inst:
            for iface, partner in rep_inst.interfaces_neighbors_map.items():
                if iface.interface_type:
                    tname = iface.interface_type.get_name()
                    if tname in [it.get_name() for it in mol_interface_types]:
                        # Only update if not already set (though rep_inst should be primary)
                        if tname not in per_type_local:
                            per_type_local[tname] = iface.absolute_coord - rep_inst.com
                            per_type_partner_ids[tname] = id(partner) if partner else -1

        # Priority 2: Fallback to any instance (Warning: May have inconsistent frame if not aligned)
        # Only for types missing from representative
        missing_types = [it.get_name() for it in mol_interface_types if it.get_name() not in per_type_local]
        
        if missing_types:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Representative instance %s missing interfaces: %s. Falling back to mixed-frame instances (RISKY).",
                    id(rep_inst), missing_types
                )
            
            for tname in missing_types:
                found = False
                for inst in all_instances:
                    for iface, partner in inst.interfaces_neighbors_map.items():
                        if iface.interface_type and iface.interface_type.get_name() == tname:
                            per_type_local[tname] = iface.absolute_coord - inst.com
                            per_type_partner_ids[tname] = id(partner) if partner else -1
                            found = True
                            break
                    if found:
                        break

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Using representative instance (id=%s) with %d interfaces for molecule type %s",
                id(rep_inst), len(per_type_local), mol_type.name
            )


        # Build stable list of interface records with site labels
        interfaces = []
        for interface_type_name, local_coord in sorted(per_type_local.items()):
            site_label = self._get_base_site_label(mol_type.name, interface_type_name)
            # Store mapping for *every* concrete type we found
            self.interface_to_site_map[interface_type_name] = site_label
            interfaces.append({
                "type_name": interface_type_name,
                "site_label": site_label,
                "coord": local_coord,
                "partner": per_type_partner_ids.get(interface_type_name, 0),
            })
            if self.workspace_manager:
                self.workspace_manager.logger.debug(
                    "Interface %s -> site %s at actual coord %s",
                    interface_type_name, site_label, local_coord
                )

        # Sort interfaces by site label for consistent output
        interfaces.sort(key=lambda x: x["site_label"])
        
        if dry_run:
             return None

        # --- existing file writing logic below unchanged ---
        with open(mol_file_path, 'w', encoding='utf-8') as f:
            f.write(f"Name = {mol_type.name}\n\n")
            f.write("checkOverlap = true\n")
            
            # Write transition matrix parameters if hyperparameters provided
            if hyperparams:
                count_trans = str(hyperparams.count_transition).lower()
                f.write(f"countTransition = {count_trans}\n")
                f.write(f"transitionMatrixSize = {hyperparams.transition_matrix_size}\n")
            
            D_t = mol_type.D_t_nm2_us; D_r = mol_type.D_r_rad2_us
            f.write("# translational diffusion constants\n")
            f.write(f"D = [{D_t:.6g}, {D_t:.6g}, {D_t:.6g}]\n\n")
            f.write("# rotational diffusion constants\n")
            f.write(f"Dr = [{D_r:.6g}, {D_r:.6g}, {D_r:.6g}]\n\n")
            f.write("# Coordinates\n")
            f.write(f"COM   {self._format_vec([0.0, 0.0, 0.0])}\n")
            for rec in interfaces:
                f.write(f"{rec['site_label']}  {self._format_vec(rec['coord'])}\n")
            f.write("\n# bonds\n")
            f.write(f"bonds = {len(interfaces)}\n")
            for rec in interfaces:
                f.write(f"com {rec['site_label']}\n")

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Wrote .mol file with %d interfaces using representative coordinates: %s",
                len(interfaces), mol_file_path
            )

        return mol_file_path


    def _get_base_site_label(self, mol_name: str, interface_type_name: str) -> str:
        """Get base site label from molecule name and interface type name.

        Format rules:
        - If both mol names are single character: A_A_1 -> aa1
        - If any mol name >= 2 characters: AH_Q_1 -> ah_q1, YDF_UU_2 -> ydf_uu_2
        - No special rule for hmodimeric heterotypic case: A_A_1f -> aa1f

        Args:
            mol_name: Molecule type name.
            interface_type_name: Interface type name (e.g., "A_A_1", "AH_Q_2").

        Returns:
            Formatted site label.
        """
        # Parse interface type name using proper parser
        try:
            parsed = interface_naming.parse_interface_name(interface_type_name)
            mol1_name = parsed.this_mol
            mol2_name = parsed.partner_mol
            index = str(parsed.index)
            if parsed.tag:
                index += parsed.tag
        except Exception as e:
            # Fallback for unexpected format
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Failed to parse interface type: %s, error: %s, using fallback",
                    interface_type_name, str(e)
                )
            initial = mol_name[0].lower() if mol_name else "x"
            return f"{initial}1"

        # Convert to lowercase
        mol1_lower = mol1_name.lower()
        mol2_lower = mol2_name.lower()

        # Apply formatting rules
        if len(mol1_name) == 1 and len(mol2_name) == 1:
            # Homodimeric labels (mol1 == mol2): use format like aa0ac11
            site_label = f"{mol1_lower}{mol2_lower}{index}"
        else:
            # Heterodimeric labels (mol1 != mol2): use format like aa0ab01 (no underscore)
            site_label = f"{mol1_lower}{mol2_lower}{index}"

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Generated site label: %s -> %s (mol1=%s, mol2=%s, index=%s)",
                interface_type_name, site_label, mol1_name, mol2_name, index
            )

        return site_label

    def _get_required_free_sites(self, mol_name: str, type_name: str) -> List[str]:
        """Get steric exclusion sites (required free) for a given interface type."""
        required_sites = []
        # Find interface type definition
        interface_type = None
        for it in self.system.interface_types:
            if it.this_mol_type_name == mol_name and it.get_name() == type_name:
                interface_type = it
                break
        
        if interface_type and hasattr(interface_type, 'required_free'):
            for req_type_name in interface_type.required_free:
                # Map required interface type to site label
                if req_type_name in self.interface_to_site_map:
                    site = self.interface_to_site_map[req_type_name]
                    if site not in required_sites:
                        required_sites.append(site)
        
        return sorted(required_sites)

    def _generate_reactions(self) -> List[str]:
        """Generate BNGL reaction strings with proper handling of failed homotypic interfaces."""
        reactions = []
        processed_pairs = set()

        # Group interface types by molecule pair and index
        interface_pairs = {}
        for interface_type in self.system.interface_types:
            mol1 = interface_type.this_mol_type_name
            mol2 = interface_type.partner_mol_type_name
            index = interface_type.interface_index

            pair_key = (mol1, mol2, index)
            if pair_key not in interface_pairs:
                interface_pairs[pair_key] = []
            interface_pairs[pair_key].append(interface_type)

        for pair_key, interface_types in interface_pairs.items():
            mol1, mol2, index = pair_key

            # Create canonical key for duplicate checking
            canonical_key = tuple(sorted([mol1, mol2]) + [index])

            if canonical_key in processed_pairs:
                continue
            processed_pairs.add(canonical_key)

            # Check if this is a failed homotypic case
            is_homodimer = (mol1 == mol2)
        
            if is_homodimer:
                # --- NEW: build f/b pairs at this (mol,index) ---
                # Gather concrete type names present for this (mol,mol,index)
                concrete = [t.get_name() for t in interface_types]
                # Keep only validly parsed names (robustness)
                parsed = []
                for tname in concrete:
                    try:
                        p = interface_naming.parse_interface_name(tname)
                        parsed.append((tname, p))
                    except Exception:
                        continue
                # Group by index; for our single "index" in key, this is already aligned
                # Find any f/b combo(s)
                fb_pairs: list[tuple[str,str]] = []
                names_set = set(t for t, _ in parsed)
                for tname, p in parsed:
                    if p.tag == 'f':
                        # Use make_interface_name to match the new format without underscores
                        candidate_b = interface_naming.make_interface_name(p.this_mol, p.partner_mol, p.index, 'b')
                        if candidate_b in names_set:
                            fb_pairs.append((tname, candidate_b))
                # For each f/b pair, map type → site and create reactions
                # For each f/b pair, map type → site and create reactions
                for t_f, t_b in fb_pairs:
                    sites_f = [s for (k, s) in self.interface_to_site_map.items() if k == t_f]
                    sites_b = [s for (k, s) in self.interface_to_site_map.items() if k == t_b]
                    
                    # Get required free sites
                    rf1 = self._get_required_free_sites(mol1, t_f)
                    rf2 = self._get_required_free_sites(mol2, t_b)
                    
                    # Format requirement strings
                    req1 = "," + ",".join(rf1) if rf1 else ""
                    req2 = "," + ",".join(rf2) if rf2 else ""

                    for s1 in sites_f:
                        for s2 in sites_b:
                            # Add steric exclusions
                            reaction = f"{mol1}({s1}{req1}) + {mol2}({s2}{req2}) <-> {mol1}({s1}!1{req1}).{mol2}({s2}!1{req2})"
                            reactions.append(reaction)
                            self.reaction_metadata.append({
                                'reaction': reaction,
                                'is_cross_reaction': False,
                                'mol1': mol1, 'mol2': mol2,
                                'site1': s1, 'site2': s2,
                                'interaction_type': 'hom_het'
                            })
                
                # Handle homodimeric homotypic (self-binding, tag=None)
                # These interfaces bind to themselves: A(aa1) + A(aa1) <-> A(aa1!1).A(aa1!1)
                homotypic_types = [tname for tname, p in parsed if p.tag is None]
                for type_name in homotypic_types:
                    # Get the site label for this interface type
                    sites = [s for (k, s) in self.interface_to_site_map.items() if k == type_name]
                    
                    # Get required free sites (same for both)
                    rf = self._get_required_free_sites(mol1, type_name)
                    req = "," + ",".join(rf) if rf else ""

                    for site in sites:
                        # Self-binding reaction: same site on both sides
                        reaction = f"{mol1}({site}{req}) + {mol2}({site}{req}) <-> {mol1}({site}!1{req}).{mol2}({site}!1{req})"
                        reactions.append(reaction)
                        self.reaction_metadata.append({
                            'reaction': reaction,
                            'is_cross_reaction': False,
                            'mol1': mol1, 'mol2': mol2,
                            'site1': site, 'site2': site,
                            'interaction_type': 'hom_hom'
                        })


            else:
                # Handle true heterodimeric cases
                # Need to find BOTH interface types (type_name and partner_type_name)
                type_name = interface_types[0].get_name()
                # Construct partner interface name  
                partner_type_name = interface_naming.make_interface_name(mol2, mol1, index, None)
                mol1_sites = []
                mol2_sites = []

                # Find sites for mol1 - use exact match since no underscores
                for key, site_label in self.interface_to_site_map.items():
                    if key == type_name:
                        if site_label not in mol1_sites:
                            mol1_sites.append(site_label)

                # Find sites for mol2
                for key, site_label in self.interface_to_site_map.items():
                    if key == partner_type_name:
                        if site_label not in mol2_sites:
                            mol2_sites.append(site_label)
                
                
                # Get required free sites
                rf1 = self._get_required_free_sites(mol1, type_name)
                rf2 = self._get_required_free_sites(mol2, partner_type_name)
                
                req1 = "," + ",".join(rf1) if rf1 else ""
                req2 = "," + ",".join(rf2) if rf2 else ""

                # Generate all combinations for heterotypic
                for site1 in mol1_sites:
                    for site2 in mol2_sites:
                        reaction = f"{mol1}({site1}{req1}) + {mol2}({site2}{req2}) <-> {mol1}({site1}!1{req1}).{mol2}({site2}!1{req2})"
                        reactions.append(reaction)

                        self.reaction_metadata.append({
                            'site1': site1, 'site2': site2,
                            'interaction_type': 'het'
                        })

        # ADDED: Include reactions from precalculated_geometry if not present
        # This allows PlatonicSolids explicit reactions (e.g. cross interactions) to be included
        existing_reactions = set(reactions)
        for (mol1, iface1, mol2, iface2) in self.precalculated_geometry.keys():
            # Map interface types to site labels
            s1_list = [s for (k, s) in self.interface_to_site_map.items() if k == iface1]
            s2_list = [s for (k, s) in self.interface_to_site_map.items() if k == iface2]
            
            # If not found, maybe the iface name IS the site label (if simple)
            if not s1_list: s1_list = [iface1]
            if not s2_list: s2_list = [iface2]
            
            for s1 in s1_list:
                for s2 in s2_list:
                    reaction = f"{mol1}({s1}) + {mol2}({s2}) <-> {mol1}({s1}!1).{mol2}({s2}!1)"
                    reaction_rev = f"{mol2}({s2}) + {mol1}({s1}) <-> {mol2}({s2}!1).{mol1}({s1}!1)"
                    
                    if reaction not in existing_reactions and reaction_rev not in existing_reactions:
                        reactions.append(reaction)
                        existing_reactions.add(reaction)
                        self.reaction_metadata.append({
                            'reaction': reaction,
                            'is_cross_reaction': (mol1 != mol2),
                            'mol1': mol1, 'mol2': mol2,
                            'site1': s1, 'site2': s2,
                            'interaction_type': 'explicit'
                        })

        return reactions

    def _calculate_reaction_parameters(self, reactions: List[str]) -> Tuple[List[float], List[Tuple[float, float, float, float, float]]]:
        """Calculate bond lengths and angles for reactions; print all contributing pairs + summary."""

        def _circ_mean_std(vals: List[float]) -> Tuple[float, float]:
            """Circular mean and circular std (sqrt(-2 ln R)) on [-π, π]."""
            if not vals:
                return 0.0, 0.0
            s = float(np.mean(np.sin(vals)))
            c = float(np.mean(np.cos(vals)))
            mean = float(np.arctan2(s, c))
            R = float(np.hypot(s, c))
            R = np.clip(R, 1e-12, 1.0)
            std = float(np.sqrt(-2.0 * np.log(R)))
            return mean, std

        sigma_list: List[float] = []
        angles_list: List[Tuple[float, float, float, float, float]] = []

        reaction_re = re.compile(
            r"^\s*([A-Za-z0-9_]+)\(([^)]+)\)\s*\+\s*([A-Za-z0-9_]+)\(([^)]+)\)"
        )

        for reaction in reactions:
            match = reaction_re.match(reaction)
            if not match:
                if self.workspace_manager:
                     self.workspace_manager.logger.warning(f"Reaction regex failed to match: {reaction}")
                sigma_list.append(1.0)
                angles_list.append((0.0, 0.0, 0.0, 0.0, 0.0))
                continue

            mol1, site1_full, mol2, site2_full = match.groups()
            
            # The site string might contain required_free interfaces (comma-separated)
            # The binding site is always the first one
            site1 = site1_full.split(',')[0].strip()
            site2 = site2_full.split(',')[0].strip()

            # Resolve sites to exact interface type names
            type1 = self._site_to_single_interface_type(site1)
            type2 = self._site_to_single_interface_type(site2)

            if not type1 or not type2:
                sigma, angles = 1.0, (0.0, 0.0, 0.0, 0.0, 0.0)
                sigma_list.append(sigma); angles_list.append(angles)
                continue

            cache_key = (mol1, site1, type1, mol2, site2, type2)

            # Cache hit?
            if cache_key in self.reaction_params_cache:
                sigma, angles = self.reaction_params_cache[cache_key]
                sigma_list.append(sigma); angles_list.append(angles)
                if self.workspace_manager:
                    self.workspace_manager.logger.info("Using cached params for %s: sigma=%.6f", cache_key, sigma)
                continue

            # Check precalculated geometry (user overrides)
            precalc_key = (mol1, type1, mol2, type2)
            if precalc_key in self.precalculated_geometry:
                sigma, angles = self.precalculated_geometry[precalc_key]
                self.reaction_params_cache[cache_key] = (sigma, angles)
                sigma_list.append(sigma); angles_list.append(angles)
                if self.workspace_manager:
                    self.workspace_manager.logger.info("Using precalculated geometry for %s: %s", precalc_key, angles)
                continue
            # Also check reverse key just in case
            precalc_key_rev = (mol2, type2, mol1, type1)
            if precalc_key_rev in self.precalculated_geometry:
                # If reverse, we might need to swap angles theta1/theta2 etc? 
                # Reaction geometry is directional: theta1 is angle on mol1.
                # If we swap mol1/mol2, we must swap theta1<->theta2 and phi1<->phi2.
                # Omega remains same? Omega is torsional.
                sigma, (th1, th2, ph1, ph2, om) = self.precalculated_geometry[precalc_key_rev]
                angles = (th2, th1, ph2, ph1, om) # Swapped
                self.reaction_params_cache[cache_key] = (sigma, angles)
                sigma_list.append(sigma); angles_list.append(angles)
                if self.workspace_manager:
                    self.workspace_manager.logger.info("Using precalculated geometry (reversed) for %s: %s", precalc_key_rev, angles)
                continue

            # Enumerate ONLY exact-type bound pairs
            pairs = self._enumerate_exact_type_pairs(mol1, type1, mol2, type2)

            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Found %d bound instance pairs for %s(%s:%s) + %s(%s:%s)",
                    len(pairs), mol1, site1, type1, mol2, site2, type2
                )

            if not pairs:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "No exact-type bound pairs for %s(%s:%s)+%s(%s:%s); using default.",
                        mol1, site1, type1, mol2, site2, type2
                    )
                sigma, angles = 1.0, (0.0, 0.0, 0.0, 0.0, 0.0)
                self.reaction_params_cache[cache_key] = (sigma, angles)
                sigma_list.append(sigma); angles_list.append(angles)
                continue

            # Compute per-pair params, then average
            sigmas = []
            angles_acc = []  # list of (theta1, theta2, phi1, phi2, omega)

            # track a "representative" pair to anchor torsion mode
            repr_pair_idx = 0
            repr_inst = self._get_representative_instance(mol1)
            for idx, (m1i, m2i, intf1, intf2) in enumerate(pairs):
                if m1i is repr_inst:        # prefer the one that contains the representative of mol1
                    repr_pair_idx = idx
                com1 = m1i.com; com2 = m2i.com
                p1   = intf1.absolute_coord; p2 = intf2.absolute_coord
                sigma_i, angles_i = self._generate_reaction_angles(
                    p1, p2, com1, com2, mol1, mol2, site1, site2
                )
                sigmas.append(sigma_i)
                angles_acc.append(angles_i)

            # arithmetic for sigma/thetas
            sigma_calc = float(np.mean(sigmas))
            
            # send out a warning for small sigma
            
            if sigma_calc < 0.5:
                self.workspace_manager.logger.warning("WARNING: small sigma values : %f; consider increasing binding radius threshold", sigma_calc)
            #sigma = max(sigma_calc, 0.5)
            sigma = sigma_calc # <- placeholder, disable arbitrary clamping
            theta1_vals = [a[0] for a in angles_acc]
            theta2_vals = [a[1] for a in angles_acc]
            theta1 = float(np.mean(theta1_vals))
            theta2 = float(np.mean(theta2_vals))

            # torsions: circular mode selection if broad
            phi1_vals = [a[2] for a in angles_acc]
            phi2_vals = [a[3] for a in angles_acc]
            omega_vals= [a[4] for a in angles_acc]

            phi1_ref  = angles_acc[repr_pair_idx][2]
            phi2_ref  = angles_acc[repr_pair_idx][3]
            omega_ref = angles_acc[repr_pair_idx][4]

            phi1_mu, R1, used1 = self._circ_mode_mean(phi1_vals,  phi1_ref, var_threshold=0.5)
            phi2_mu, R2, used2 = self._circ_mode_mean(phi2_vals,  phi2_ref, var_threshold=0.5)
            omg_mu , Rw, usedw = self._circ_mode_mean(omega_vals, omega_ref, var_threshold=0.5)

            angles = (theta1, theta2, phi1_mu, phi2_mu, omg_mu)

            # Diagnostics / warnings
            if self.workspace_manager:
                def _cvar(R): return (1.0 - R)
                self.workspace_manager.logger.info(
                    ("Angle dispersion for %s(%s)+%s(%s): "
                    "σ_std=%.4g | Var(θ1)=%.4g Var(θ2)=%.4g | CVar(φ1)=%.4g%s CVar(φ2)=%.4g%s CVar(ω)=%.4g%s  [N=%d]"),
                    mol1, site1, mol2, site2,
                    float(np.std(sigmas)),
                    float(np.var(theta1_vals)), float(np.var(theta2_vals)),
                    _cvar(R1), " [MODE]" if used1 else "",
                    _cvar(R2), " [MODE]" if used2 else "",
                    _cvar(Rw), " [MODE]" if usedw else "",
                    len(sigmas)
                )

            # Cache & store
            self.reaction_params_cache[cache_key] = (sigma_calc, angles)
            sigma_list.append(sigma_calc); angles_list.append(angles)

        return sigma_list, angles_list


    def _get_coms_interfaces(self, mol1_name: str, site1: str, mol2_name: str, site2: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
        """Get center of mass and interface site coordinates for two specific molecule instances.

        For homotypic reactions (A + A), ensures we get two different instances.
        """
        # For homotypic reactions, we need to find two different instances
        is_homotypic = (mol1_name == mol2_name)

        if is_homotypic:
            return self._get_coms_interfaces_homotypic(mol1_name, site1, site2)
        else:
            return self._get_coms_interfaces_heterotypic(mol1_name, site1, mol2_name, site2)

    def _get_coms_interfaces_homotypic(self, mol_name: str, site1: str, site2: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
        """Get coordinates for homotypic reaction using representative for site1."""
        
        # Always use representative instance for site1 (mol1)
        representative_instance = self._get_representative_instance(mol_name)
        if not representative_instance:
            return None, None, None, None
        
        # Find the interface on representative that matches site1
        repr_interface = None
        repr_coord = None
        
        for interface, partner in representative_instance.interfaces_neighbors_map.items():
            if self._interface_matches_site(interface, site1):
                repr_interface = interface
                repr_coord = interface.absolute_coord
                
                # Find the partner instance that has site2
                partner_instance = partner
                partner_coord = None
                
                # Find partner's interface that matches site2 and connects back to repr
                for partner_interface, partner_neighbor in partner_instance.interfaces_neighbors_map.items():
                    if (partner_neighbor == representative_instance and
                        self._interface_matches_site(partner_interface, site2)):
                        partner_coord = partner_interface.absolute_coord
                        break
                
                if partner_coord is not None:
                    if self.workspace_manager:
                        self.workspace_manager.logger.info(
                            "Homotypic reaction %s(%s) + %s(%s): repr=%s (site1) -> partner=%s (site2)",
                            mol_name, site1, mol_name, site2,
                            id(representative_instance), id(partner_instance)
                        )
                    
                    return (representative_instance.com, partner_instance.com, 
                            repr_coord, partner_coord)
        
        if self.workspace_manager:
            self.workspace_manager.logger.warning(
                "Could not find binding pair for %s(%s) + %s(%s) using representative",
                mol_name, site1, mol_name, site2
            )
        
        return None, None, None, None

    def _get_coms_interfaces_heterotypic(self, mol1_name: str, site1: str, mol2_name: str, site2: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
        """Get coordinates for heterotypic reaction (A + B)."""
        # Find instances that are actually binding to each other
        for mol_instance in self.system.molecule_instances:
            if mol_instance.molecule_type and mol_instance.molecule_type.name == mol1_name:
                for this_interface, neighbor_instance in mol_instance.interfaces_neighbors_map.items():
                    if (neighbor_instance.molecule_type and
                        neighbor_instance.molecule_type.name == mol2_name and
                            self._interface_matches_site(this_interface, site1)):

                        # Find the corresponding interface on the partner
                        for partner_interface, partner_neighbor in neighbor_instance.interfaces_neighbors_map.items():
                            if (partner_neighbor == mol_instance and
                                    self._interface_matches_site(partner_interface, site2)):

                                return (mol_instance.com, neighbor_instance.com,
                                        this_interface.absolute_coord, partner_interface.absolute_coord)

        return None, None, None, None

    def _interface_matches_site(self, interface_instance, site_label: str) -> bool:
        """Check if an interface instance corresponds to a site label."""
        # Get the interface type name
        interface_type_name = interface_instance.interface_type.get_name()

        # Check if this interface type maps to the site label
        for key, mapped_site in self.interface_to_site_map.items():
            if mapped_site == site_label and key == interface_type_name:
                return True

        return False

    def _find_instance_from_coordinates(self, mol_name: str, com: np.ndarray,
                                        intf_coord: np.ndarray) -> Optional:
        """Find molecule instance based on COM and interface coordinates.

        Args:
            mol_name: Molecule type name.
            com: Center of mass coordinates.
            intf_coord: Interface coordinates.

        Returns:
            Matching molecule instance or None.
        """
        tolerance = 1e-6  # Coordinate matching tolerance

        for mol_instance in self.system.molecule_instances:
            if (mol_instance.molecule_type and
                    mol_instance.molecule_type.name == mol_name):

                # Check COM match
                if np.linalg.norm(mol_instance.com - com) < tolerance:
                    # Check if any interface coordinate matches
                    for interface, _ in mol_instance.interfaces_neighbors_map.items():
                        if np.linalg.norm(interface.absolute_coord - intf_coord) < tolerance:
                            return mol_instance

        return None

    def _generate_reaction_angles(self, intf1: np.ndarray, intf2: np.ndarray,
                                com1: np.ndarray, com2: np.ndarray,
                                mol1_name: str, mol2_name: str,
                                site1: str, site2: str,
                                tol: float = 1e-8) -> Tuple[float, Tuple[float, float, float, float, float]]:
        """
        Infer NERDSS-style reaction parameters (σ, θ₁, θ₂, φ₁, φ₂, ω) using 
        reference vectors (ref1) from MoleculeInstance.
        
        Logic adopted from ref_angles.py.
        """
        
        # DEBUG: Log all inputs
        if self.workspace_manager:
            self.workspace_manager.logger.debug("="*80)
            self.workspace_manager.logger.debug(f"ANGLE CALC for {mol1_name}({site1}) + {mol2_name}({site2})")
            self.workspace_manager.logger.debug(f"  COM1:   {com1}")
            self.workspace_manager.logger.debug(f"  COM2:   {com2}")
            self.workspace_manager.logger.debug(f"  intf1:  {intf1}")
            self.workspace_manager.logger.debug(f"  intf2:  {intf2}")
        
        mol1_instance = self._find_instance_from_coordinates(mol1_name, com1, intf1)
        mol2_instance = self._find_instance_from_coordinates(mol2_name, com2, intf2)
        
        if mol1_instance is None or mol2_instance is None:
             if self.workspace_manager:
                 self.workspace_manager.logger.error("Could not find instances for angle calculation")
             return 0.0, (0.0, 0.0, 0.0, 0.0, 0.0)

        # helpers
        def mag(x): 
            return np.linalg.norm(x)
        def unit(x):
            m = mag(x)
            return x / m if m > 1e-12 else np.zeros_like(x)

        # 1. Determine normal points (ref points) from ref1
        # Use the same normal vector logic as parms.inp to ensure consistent reference frames
        n1 = self._n_global_from_local_x(mol1_name, site1, mol1_instance)
        # Check if we need to normalize (though _n_global usually returns unit vectors)
        n1 = unit(n1)
        if self.workspace_manager:
            self.workspace_manager.logger.debug(f"  n1 (norm): {n1}")
        
        n2 = self._n_global_from_local_x(mol2_name, site2, mol2_instance)
        n2 = unit(n2)
        if self.workspace_manager:
            self.workspace_manager.logger.debug(f"  n2 (norm): {n2}")

        # 2. Basic vectors
        v1 = intf1 - com1
        v2 = intf2 - com2
        sigma1 = intf1 - intf2
        sigma2 = intf2 - intf1
        sigma_magnitude = mag(sigma1)
        
        if self.workspace_manager:
            self.workspace_manager.logger.debug(f"  v1: {v1}, |v1|={mag(v1):.6f}")
            self.workspace_manager.logger.debug(f"  v2: {v2}, |v2|={mag(v2):.6f}")
            self.workspace_manager.logger.debug(f"  sigma: {sigma_magnitude:.6f}")
        
        if sigma_magnitude < tol:
            return 0.0, (0.0, 0.0, 0.0, 0.0, 0.0)

        # 3. Calculate theta
        # theta = acos( (v . sigma) / (|v| |sigma|) )
        theta1 = math.acos(np.clip(np.dot(v1, sigma1) / (mag(v1) * sigma_magnitude), -1.0, 1.0))
        theta2 = math.acos(np.clip(np.dot(v2, sigma2) / (mag(v2) * sigma_magnitude), -1.0, 1.0))

        # 4. Calculate phi
        # t1 = unit(cross(v, sigma))
        # t2 = unit(cross(v, n))
        # phi = acos( t1 . t2 )
        # For linear molecules: if molecule has only 1 interface, phi is undefined (set to NaN)
        
        # Check if molecules have only 1 interface (linear molecule case)
        # Get the molecule instances to check their interface count
        mol1_instance = self._find_instance_from_coordinates(mol1_name, com1, intf1)
        mol2_instance = self._find_instance_from_coordinates(mol2_name, com2, intf2)
        
        # Count interfaces for each molecule type (excluding reference vectors)
        mol1_interface_count = len(mol1_instance.molecule_type.interfaces_neighbors_map) if mol1_instance and mol1_instance.molecule_type else 0
        mol2_interface_count = len(mol2_instance.molecule_type.interfaces_neighbors_map) if mol2_instance and mol2_instance.molecule_type else 0
        
        # Calculate phi1
        if mol1_interface_count == 1:
            # Molecule has only 1 interface - phi1 is physically meaningless for linear molecules
            phi1 = float('nan')
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    f"Linear molecule detected for {mol1_name}({site1}): only 1 interface, setting phi1=NaN")
        else:
            t1_1 = unit(np.cross(v1, sigma1))
            t2_1 = unit(np.cross(v1, n1))
            phi1 = math.acos(np.clip(np.dot(t1_1, t2_1), -1.0, 1.0))
        
        # Calculate phi2
        if mol2_interface_count == 1:
            # Molecule has only 1 interface - phi2 is physically meaningless for linear molecules
            phi2 = float('nan')
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    f"Linear molecule detected for {mol2_name}({site2}): only 1 interface, setting phi2=NaN")
        else:
            t1_2 = unit(np.cross(v2, sigma2))
            t2_2 = unit(np.cross(v2, n2))
            phi2 = math.acos(np.clip(np.dot(t1_2, t2_2), -1.0, 1.0))

        # 5. Determine sign of phi
        # Project n and sigma onto plane perpendicular to v
        # Skip this for linear molecules where phi is NaN
        v1_uni = unit(v1)
        v2_uni = unit(v2)
        
        # Only calculate sign for non-NaN phi values
        if not np.isnan(phi1):
            n1_proj = n1 - v1_uni * np.dot(v1_uni, n1)
            sigma1_proj = sigma1 - v1_uni * np.dot(v1_uni, sigma1)
            phi1_dir = unit(np.cross(sigma1_proj, n1_proj))
            
            # Determine sign of phi - using full 3D vector comparison (robust for arbitrary orientations)
            # Check if v_uni and phi_dir are parallel (dot ≈ 1) or anti-parallel (dot ≈ -1)
            tol_sign = 1e-6
            dot_v1_phi1 = np.dot(v1_uni, phi1_dir)
            if abs(dot_v1_phi1 - 1.0) < tol_sign:  # parallel
                phi1 = -phi1
            elif abs(dot_v1_phi1 + 1.0) < tol_sign:  # anti-parallel
                phi1 = phi1
            else:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        f"Phi1 sign ambiguous: dot(v1,phi1_dir)={dot_v1_phi1:.6f}")
        
        if not np.isnan(phi2):
            n2_proj = n2 - v2_uni * np.dot(v2_uni, n2)
            sigma2_proj = sigma2 - v2_uni * np.dot(v2_uni, sigma2)
            phi2_dir = unit(np.cross(sigma2_proj, n2_proj))
            
            tol_sign = 1e-6
            dot_v2_phi2 = np.dot(v2_uni, phi2_dir)
            if abs(dot_v2_phi2 - 1.0) < tol_sign:  # parallel
                phi2 = -phi2
            elif abs(dot_v2_phi2 + 1.0) < tol_sign:  # anti-parallel
                phi2 = phi2
            else:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        f"Phi2 sign ambiguous: dot(v2,phi2_dir)={dot_v2_phi2:.6f}")

        # 6. Calculate omega
        # a1 = cross(sigma1, v1)
        # a2 = cross(sigma1, v2)
        # omega = acos(a1 . a2)
        a1 = unit(np.cross(sigma1, v1))
        a2 = unit(np.cross(sigma1, v2))
        
        # fallback for degenerate omega logic not explicitly in ref_angles.py but in previous code
        # retaining basic logic from ref_angles.py which assumes non-degenerate
        omega = math.acos(np.clip(np.dot(a1, a2), -1.0, 1.0))

        # 7. Determine sign of omega
        sigma1_uni = unit(sigma1)
        # Project v1, v2 onto plane perpendicular to sigma1
        v1_proj_om = v1 - sigma1_uni * np.dot(sigma1_uni, v1)
        v2_proj_om = v2 - sigma1_uni * np.dot(sigma1_uni, v2)
        
        omega_dir = unit(np.cross(v1_proj_om, v2_proj_om))
        
        # Determine sign of omega - using full 3D vector comparison (robust for arbitrary orientations)
        dot_sigma_omega = np.dot(sigma1_uni, omega_dir)
        if abs(dot_sigma_omega - 1.0) < tol_sign:  # parallel
            omega = -omega
        elif abs(dot_sigma_omega + 1.0) < tol_sign:  # anti-parallel
            omega = omega
        else:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    f"Omega sign ambiguous: dot(sigma1,omega_dir)={dot_sigma_omega:.6f}")
            
        if self.workspace_manager:
              self.workspace_manager.logger.debug(
                  "Calculated angles (ref-based) for %s(%s) + %s(%s): σ=%.4f | θ₁=%.4f θ₂=%.4f φ₁=%.4f φ₂=%.4f ω=%.4f",
                  mol1_name, site1, mol2_name, site2,
                  sigma_magnitude, theta1, theta2, phi1, phi2, omega
              )

        return sigma_magnitude, (theta1, theta2, phi1, phi2, omega)



    
    def _calculate_auto_time_step(self, reactions: List[str], molecule_counts: Dict[str, int],
                                  box_nm: Tuple[float, float, float], sigma_list: List[float]) -> Optional[float]:
        """Calculate automatic time step based on stability criteria.
        
        Formula: Δt = (1 / 56(DA+DB)) * [ ((3/(4πρ)) + σ³)^(1/3) - σ ]²
        """
        V = box_nm[0] * box_nm[1] * box_nm[2]
        if V <= 0: return None
        
        diff_map = {m.name: m.D_t_nm2_us for m in self.system.molecule_types}
        min_dt = float('inf')
        
        reaction_re = re.compile(
            r"^\s*([A-Za-z0-9_]+)\(([A-Za-z0-9_]+)\)\s*\+\s*([A-Za-z0-9_]+)\(([A-Za-z0-9_]+)\)"
        )
        
        found_interaction = False
        
        for i, rxn_str in enumerate(reactions):
            match = reaction_re.match(rxn_str)
            if not match: continue
            
            mol1, _, mol2, _ = match.groups()
            sigma = sigma_list[i]
            
            D1 = diff_map.get(mol1, 0.0)
            D2 = diff_map.get(mol2, 0.0)
            D_sum = D1 + D2
            
            if D_sum <= 1e-12: continue
            
            # Check both directions for density dependence
            for mol_target in [mol1, mol2]:
                N = molecule_counts.get(mol_target, 0)
                if N <= 0: continue
                
                rho = N / V
                if rho <= 0: continue
                
                # Formula: dt = (1 / 56(D1+D2)) * [ ((3/(4 pi rho)) + sigma^3)^(1/3) - sigma ]^2
                term1 = 1.0 / (56.0 * D_sum)
                
                r_avg_term = 3.0 / (4.0 * math.pi * rho)
                bracket_inner = r_avg_term + (sigma**3)
                bracket = (bracket_inner**(1.0/3.0)) - sigma
                
                if bracket < 0: bracket = 0.0
                
                dt = term1 * (bracket**2)
                
                if dt < min_dt:
                    min_dt = dt
                    found_interaction = True
                    
        return min_dt if found_interaction else None

    def calculate_simulation_timestep(self, molecule_counts: Dict[str, int], 
                                    box_nm: Tuple[float, float, float]) -> Optional[float]:
        """Calculate the required time step based on system geometry and physics.
        
        This method runs a partial export process (without writing files) to determine
        the stable time step for the system.
        """
        # Clear state
        self.interface_to_site_map.clear()
        self.reaction_metadata.clear()
        self.homotypic_interface_map.clear()
        self.calculated_normals.clear()
        self.reaction_params_cache.clear()

        # Build mappings via dry-run of mol file generation
        for mol_type in self.system.molecule_types:
            self._write_mol_file(mol_type, dry_run=True)
            
        # Calculate normals
        self._calculate_normal_vectors()
        
        # Generate reactions
        reactions = self._generate_reactions()
        
        # Calculate parameters (sigma needed for timestep)
        sigma_list, _ = self._calculate_reaction_parameters(reactions)
        
        # Calculate timestep
        return self._calculate_auto_time_step(reactions, molecule_counts, box_nm, sigma_list)

    def _write_parms_file(self, reactions: List[str], molecule_counts: Dict[str, int],
                          box_nm: Tuple[float, float, float], sigma_list: List[float],
                          angles_list: List[Tuple[float, float, float, float, float]],
                          parms_overrides: Optional[Dict[str, Any]] = None) -> Path:
        """Write parms.inp file with calculated normal vectors."""
        parms_path = self.output_dir / "parms.inp"

        # Default parameters
        # NOTE: onRate3Dka and offRatekb are now calculated per-reaction based on interface energies
        # The values below are only used as fallback defaults if energy data is unavailable
        params = {
            'nItr': 1e5,
            'timestep': 0.5,
            'timeWrite': 1e3,
            'trajWrite': 1e5,
            'restartWrite': 1e5,
            'checkPoint': 1e5,
            'pdbWrite': 1e5,
            'overlapSepLimit': 2.0,
            'scaleMaxDisplace': 100.0,
        }
        
        # Extract default_ka from hyperparams provided in overrides
        default_ka_val = 120.0
        if parms_overrides and 'hyperparams' in parms_overrides:
            hp = parms_overrides['hyperparams']
            if hasattr(hp, 'default_on_rate_3d_ka'):
                default_ka_val = hp.default_on_rate_3d_ka

        # Add transitionWrite from hyperparams if provided
        if parms_overrides and 'hyperparams' in parms_overrides:
            hyperparams = parms_overrides['hyperparams']
            if hasattr(hyperparams, 'transition_write') and hyperparams.transition_write is not None:
                params['transitionWrite'] = hyperparams.transition_write

        # Apply overrides
        if parms_overrides:
            # Create copy to avoid modifying original or injecting objects
            safe_overrides = parms_overrides.copy()
            if 'hyperparams' in safe_overrides:
                del safe_overrides['hyperparams']
            params.update(safe_overrides)

        # --- Time Step Logic ---
        # 1. Start with default or overridden value (from params)
        final_dt = params['timestep']
        
        # 2. Check for auto-calculation if not explicitly forced by hyperparams
        hyperparams_dt = None
        if parms_overrides and 'hyperparams' in parms_overrides:
            hp = parms_overrides['hyperparams']
            if hasattr(hp, 'nerdss_time_step'):
                hyperparams_dt = hp.nerdss_time_step

        if hyperparams_dt is not None:
             # Force using the hyperparameter value
             final_dt = hyperparams_dt
             if self.workspace_manager:
                 self.workspace_manager.logger.info(
                     f"Using forced NERDSS time step from hyperparameters: {final_dt} us")
        else:
             # Auto-calculate if no forced hyperparameter
             auto_dt = self._calculate_auto_time_step(reactions, molecule_counts, box_nm, sigma_list)
             if auto_dt is not None:
                 final_dt = auto_dt
                 if self.workspace_manager:
                     self.workspace_manager.logger.info(
                         f"Using automatically calculated NERDSS time step: {final_dt:.6f} us")
             else:
                 if self.workspace_manager:
                     self.workspace_manager.logger.info(
                         f"Could not auto-calculate time step (no dynamic species?). Using default/override: {final_dt} us")

        params['timestep'] = final_dt
        # -----------------------

        # Regex to parse reactions
        reaction_re = re.compile(
            r"^\s*([A-Za-z0-9_]+)\(([A-Za-z0-9_]+)\)\s*\+\s*([A-Za-z0-9_]+)\(([A-Za-z0-9_]+)\)"
        )

        with open(parms_path, 'w', encoding='utf-8') as f:
            # Parameters section
            f.write("start parameters\n")
            for key, val in params.items():
                f.write(f"    {key} = {val}\n")
            f.write("end parameters\n\n")

            # Boundaries section
            f.write("start boundaries\n")
            f.write(
                f"    WaterBox = [{box_nm[0]}, {box_nm[1]}, {box_nm[2]}] #nm\n")
            f.write("end boundaries\n\n")

            # Molecules section
            f.write("start molecules\n")
            
            # Filter molecule_counts to only include molecules that:
            # 1. Have corresponding molecule types in the system
            # 2. Have at least one instance (so a .mol file was created)
            mol_type_names = {mol_type.name for mol_type in self.system.molecule_types}
            for mol_name, count in molecule_counts.items():
                if mol_name not in mol_type_names:
                    if self.workspace_manager:
                        self.workspace_manager.logger.warning(
                            "Skipping molecule '%s' in parms.inp - no corresponding molecule type found (may have been renamed)",
                            mol_name
                        )
                    continue
                
                # Check if this molecule has a representative instance (i.e., .mol file was created)
                if self._get_representative_instance(mol_name) is None:
                    if self.workspace_manager:
                        self.workspace_manager.logger.warning(
                            "Skipping molecule '%s' in parms.inp - no instances found (no .mol file created)",
                            mol_name
                        )
                    continue
                    
                f.write(f"    {mol_name} : {count}\n")
            
            f.write("end molecules\n\n")

            # Reactions section
            f.write("start reactions\n")
            f.write("    \n")
            f.write("    # Binding reactions\n")

            for i, reaction in enumerate(reactions):
                f.write(f"    {reaction}\n")

                # Parse reaction to get molecule types and sites
                match = reaction_re.match(reaction)
                if match:
                    mol1, site1, mol2, site2 = match.groups()

                    # Get calculated normal vectors for each species and site
                    norm1_local = self._local_x_with_degeneracy(mol1, site1)
                    norm2_local = self._local_x_with_degeneracy(mol2, site2)

                else:
                    # Fallback
                    norm1_local = np.array([0.0, 0.0, 1.0])
                    norm1_local = np.array([0.0, 0.0, 1.0])
                    norm2_local = np.array([0.0, 0.0, 1.0])

                # Determine Rates
                ka_val = default_ka_val # Default from hyperparams or fallback (nm^3/us)
                kb_val = 1000.0 # Default fallback (s^-1)
                
                # Check precalculated rates
                rate_key = (mol1, site1, mol2, site2)
                rate_key_rev = (mol2, site2, mol1, site1)
                
                if rate_key in self.precalculated_rates:
                    ka_val, kb_val = self.precalculated_rates[rate_key]
                elif rate_key_rev in self.precalculated_rates:
                    ka_val, kb_val = self.precalculated_rates[rate_key_rev]
                else:
                    # Calculate from energy if available
                    # Standard Concentration C0 ~ 0.6022 nm^-3 (1 Molar)
                    C0 = 0.602214076
                    
                    # Look up energy from interface types
                    # Resolve sites to types
                    t1_name = self._site_to_single_interface_type(site1)
                    # t2_name = self._site_to_single_interface_type(site2) # assumed symmetric energy usually
                    
                    # Find InterfaceType object to get energy
                    delta_G_raw = -1.0 # Default binding energy (kJ/mol if ProAffinity, else flag -1.0)
                    
                    found_energy = False
                    for it in self.system.interface_types:
                        if it.get_name() == t1_name:
                            if it.energy is not None:
                                delta_G_raw = it.energy
                                found_energy = True
                            break
                    
                    # Calculate kb
                    # kb = ka * C0 * exp(delta_G/RT)
                    # ka is in nm^3/us. C0 is nm^-3. Product is 1/us.
                    # We output kb in s^-1, so multiply by 1e6.
                    
                    prefactor_us = ka_val * C0

                    # Convert energy to RT units or handle default
                    R_kJ = 0.008314
                    T = 298.0
                    
                    if delta_G_raw == -1.0:
                        # Flag for default strong binding: -16RT
                        # delta_G / RT = -16
                        exponent = -16.0
                    else:
                        # Assumed explicitly set in kJ/mol (e.g. from ProAffinity)
                        # Normalize by RT
                        exponent = delta_G_raw / (R_kJ * T)

                    exp_term = math.exp(exponent)
                    kb_us = prefactor_us * exp_term # units: us^-1
                    
                    kb_val = kb_us * 1e6 # Convert to s^-1

                f.write(f"    onRate3Dka = {ka_val}\n")
                f.write(f"    offRatekb = {kb_val}\n")


                # Write calculated normal vectors
                f.write(
                    f"    norm1 = [{norm1_local[0]}, {norm1_local[1]}, {norm1_local[2]}]\n")
                f.write(
                    f"    norm2 = [{norm2_local[0]}, {norm2_local[1]}, {norm2_local[2]}]\n")

                # Use calculated sigma and angles
                sigma = sigma_list[i] if i < len(sigma_list) else 1.0
                f.write(f"    sigma = {sigma}\n")

                angles = angles_list[i] if i < len(
                    angles_list) else (0.0, 0.0, 0.0, 0.0, 0.0)
                f.write(
                    f"    assocAngles = [{angles[0]},{angles[1]},{angles[2]},{angles[3]},{angles[4]}]\n")
                f.write("    \n")

            f.write("end reactions\n")

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Wrote parms.inp file: %s", parms_path)

        return parms_path

    def _format_vec(self, vec: Iterable[float], precision: int = 7) -> str:
        """Format vector for output files."""
        return "   ".join(f"{float(v):.{precision}f}" for v in vec)

    #######
    #debug
    #######
    
    def _debug_representative_instance(self, mol_name: str):
        """Debug print all information about the representative instance."""
        
        representative = self._get_representative_instance(mol_name)
        
        if not representative:
            self.workspace_manager.logger.debug(f"DEBUG REPRESENTATIVE: No representative found for {mol_name}")
            return
        self.workspace_manager.logger.debug(f"DEBUG REPRESENTATIVE INSTANCE for {mol_name}:")
        self.workspace_manager.logger.debug(f"=" * 60)
        self.workspace_manager.logger.debug(f"Instance ID: {id(representative)}")
        self.workspace_manager.logger.debug(f"COM (absolute): {representative.com}")
        self.workspace_manager.logger.debug(f"Number of interfaces: {len(representative.interfaces_neighbors_map)}")
        
        self.workspace_manager.logger.debug("INTERFACES AND BINDING PARTNERS:")
        for i, (interface, partner) in enumerate(representative.interfaces_neighbors_map.items(), 1):
            interface_type_name = interface.interface_type.get_name()
            interface_absolute_coord = interface.absolute_coord
            interface_local_coord = interface_absolute_coord - representative.com
            
            # Find the site label for this interface
            site_label = "UNKNOWN"
            for key, label in self.interface_to_site_map.items():
                if key.startswith(interface_type_name) or key == interface_type_name:
                    site_label = label
                    break
            
            self.workspace_manager.logger.debug(f"  Interface {i}:")
            self.workspace_manager.logger.debug(f"    Type: {interface_type_name}")
            self.workspace_manager.logger.debug(f"    Site label: {site_label}")
            self.workspace_manager.logger.debug(f"    Absolute coord: {interface_absolute_coord}")
            self.workspace_manager.logger.debug(f"    Local coord (relative to COM): {interface_local_coord}")
            if partner:
                self.workspace_manager.logger.debug(f"    Partner molecule ID: {id(partner)}")
                self.workspace_manager.logger.debug(f"    Partner molecule COM: {partner.com}")
                
                # Find the partner's interface that connects back
                partner_interface = None
                for p_interface, p_neighbor in partner.interfaces_neighbors_map.items():
                    if p_neighbor == representative:
                        partner_interface = p_interface
                        break
            else:
                 self.workspace_manager.logger.debug("    Partner molecule: None (Unbound)")
                 partner_interface = None
            
            if partner_interface:
                partner_type_name = partner_interface.interface_type.get_name()
                partner_absolute_coord = partner_interface.absolute_coord
                partner_local_coord = partner_absolute_coord - partner.com
                
                # Find partner site label
                partner_site_label = "UNKNOWN"
                for key, label in self.interface_to_site_map.items():
                    if key.startswith(partner_type_name) or key == partner_type_name:
                        partner_site_label = label
                        break
                
                self.workspace_manager.logger.debug(f"    Partner interface type: {partner_type_name}")
                self.workspace_manager.logger.debug(f"    Partner site label: {partner_site_label}")
                self.workspace_manager.logger.debug(f"    Partner interface absolute coord: {partner_absolute_coord}")
                self.workspace_manager.logger.debug(f"    Partner interface local coord: {partner_local_coord}")
                self.workspace_manager.logger.debug(f"    Bond length: {np.linalg.norm(interface_absolute_coord - partner_absolute_coord):.6f}")
            else:
                self.workspace_manager.logger.debug(f"    (No connected partner interface found)")
        
        self.workspace_manager.logger.debug("INTERFACE-TO-SITE MAPPING:")
        self.workspace_manager.logger.debug("Interface type -> Site label:")
        for key, site_label in self.interface_to_site_map.items():
            self.workspace_manager.logger.debug(f"  {key} -> {site_label}")
        
        self.workspace_manager.logger.debug("EXPECTED REACTIONS (based on interface types):")
        interface_types = [intf.interface_type.get_name()
                           for intf in representative.interfaces_neighbors_map.keys()]
        # f/b complementary preview (homodimeric heterotypic)
        for t in interface_types:
            try:
                p = interface_naming.parse_interface_name(t)
                if p.this_mol == p.partner_mol and p.tag == 'f':
                    partner = f"{p.this_mol}_{p.partner_mol}_{p.index}b"
                    if partner in interface_types:
                        s1 = self._get_site_label_for_interface_type(t)
                        s2 = self._get_site_label_for_interface_type(partner)
                        self.workspace_manager.logger.info(f"  {mol_name}({s1}) + {mol_name}({s2}) <-> {mol_name}({s1}!1).{mol_name}({s2}!1)")
            except Exception:
                continue
        
        self.workspace_manager.logger.info("=" * 60)

    def _get_site_label_for_interface_type(self, interface_type_name: str) -> str:
        """Get site label for a given interface type name."""
        for key, site_label in self.interface_to_site_map.items():
            if key == interface_type_name:
                return site_label
        return "UNKNOWN"
    
    def _iface_type_name(self, interface_instance) -> str:
        return interface_instance.interface_type.get_name() if interface_instance and interface_instance.interface_type else ""

    def _site_to_single_interface_type(self, site_label: str) -> str:
        """
        Resolve a site label to exactly ONE interface type name.
        If multiple or none are found, emit a clear warning and choose deterministically.
        """
        hits = [t for t, s in self.interface_to_site_map.items() if s == site_label]
        if not hits:
            if self.workspace_manager:
                self.workspace_manager.logger.error("Site %s maps to no interface types.", site_label)
            return ""  # caller will guard
        uniq = sorted(set(hits))
        if len(uniq) > 1:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Site %s maps to multiple interface types %s; using the first: %s",
                    site_label, uniq, uniq[0]
                )
        return uniq[0]

    def _rep_local_site_vector(self, mol_name: str, site_label: str) -> Optional[np.ndarray]:
        """Return the representative instance's LOCAL coordinate of this site (p_local)."""
        rep = self._get_representative_instance(mol_name)
        if rep is None:
            return None
        # Resolve site -> exact interface type name
        iface_type = self._site_to_single_interface_type(site_label)
        if not iface_type:
            return None
        # Find the matching interface on the representative
        for intf, _ in rep.interfaces_neighbors_map.items():
            if self._iface_type_name(intf) == iface_type:
                return (intf.absolute_coord - rep.com)  # local vector in representative frame
        return None

    def _local_x_with_degeneracy(self, mol_name: str, site_label: str, thr: float = 0.99) -> np.ndarray:
        """
        Choose the local base normal for a site in the REPRESENTATIVE frame.
        Default [0,0,1] (Z-axis); if collinear with the site's local vector (at rep), use [1,0,0].
        
        Reverted to 0.99 to match user manual reference frame behavior.
        """
        n_local = np.array([0.0, 0.0, 1.0], dtype=float)
        p_local = self._rep_local_site_vector(mol_name, site_label)
        if p_local is None:
            return n_local
        nl = np.linalg.norm(p_local)
        if nl > 1e-12:
            vhat = p_local / nl
            if abs(float(np.dot(n_local, vhat))) > thr:
                return np.array([1.0, 0.0, 0.0], dtype=float)
        return n_local

    def _n_global_from_local_x(self, mol_name: str, site_label: str, inst) -> np.ndarray:
        """
        Take the base local normal (chosen in representative frame) and rotate it
        to the *instance's* global frame using the existing R(repr->inst).
        """
        rep = self._get_representative_instance(mol_name)
        if rep is None or inst is None or inst is rep:
            # No rotation needed (rep itself)
            return self._local_x_with_degeneracy(mol_name, site_label)

        R = self._calculate_rotation_from_representative(mol_name, rep, inst)
        if R is None:
            # Fallback: no rotation available
            return self._local_x_with_degeneracy(mol_name, site_label)

        n_local = self._local_x_with_degeneracy(mol_name, site_label)
        return R @ n_local
    
    def _circular_mean_std(self, vals: List[float]) -> Tuple[float, float]:
        """
        Circular mean and 'circular std' (sqrt(-2 ln R)) on [-π, π].
        If vals empty, returns (0, 0).
        """
        if not vals:
            return 0.0, 0.0
        s = float(np.mean(np.sin(vals)))
        c = float(np.mean(np.cos(vals)))
        mean = float(np.arctan2(s, c))
        R = float(np.hypot(s, c))
        # Clamp R into (0,1] to avoid log issues
        R = np.clip(R, 1e-12, 1.0)
        std = float(np.sqrt(-2.0 * np.log(R)))
        return mean, std


    def _circ_mean_R(self, angles: list[float]) -> tuple[float, float]:
        """Return (circular_mean, resultant_length R in [0,1])."""
        if not angles:
            return 0.0, 0.0
        s = np.mean(np.sin(angles)); c = np.mean(np.cos(angles))
        mu = float(np.arctan2(s, c))
        R = float(np.hypot(s, c))
        return mu, R

    def _circ_mode_mean(self, angles: list[float], ref_angle: float, var_threshold: float = 0.5) -> tuple[float, float, bool]:
        """
        If the circular resultant R < (1 - var_threshold), treat as bimodal and
        choose the cluster closest to ref_angle. Returns (mean, R_selected, used_mode).
        """
        mu_all, R_all = self._circ_mean_R(angles)
        # circular variance ≈ 1 - R; high if R small
        if (1.0 - R_all) <= var_threshold:
            return mu_all, R_all, False

        # --- simple 2-means on the unit circle with seeded antipodal centers ---
        ref = np.array([np.cos(ref_angle), np.sin(ref_angle)], float)
        ctrs = np.stack([ref, -ref], axis=0)

        pts = np.stack([np.cos(angles), np.sin(angles)], axis=1)
        assign = None

        for _ in range(10):
            # cosine distance maximization == dot maximization
            dots = pts @ ctrs.T                         # N x 2
            new_assign = np.argmax(dots, axis=1)        # 0 or 1
            if assign is not None and np.all(new_assign == assign):
                break
            assign = new_assign
            for k in (0, 1):
                sel = pts[assign == k]
                if len(sel) > 0:
                    v = sel.mean(axis=0)
                    n = np.linalg.norm(v)
                    ctrs[k] = v / n if n > 1e-12 else ctrs[k]

        # choose cluster whose center is closer to ref on the circle
        k_ref = np.argmax(ctrs @ ref)
        sel = pts[assign == k_ref]
        if len(sel) == 0:
            # fallback to all
            return mu_all, R_all, False

        mu = float(np.arctan2(sel[:,1].mean(), sel[:,0].mean()))
        R  = float(np.linalg.norm(sel.mean(axis=0)))
        return mu, R, True
