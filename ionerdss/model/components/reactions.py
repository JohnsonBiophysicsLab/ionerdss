"""
ionerdss.model.components.reactions

Reaction definitions and rules for molecular binding simulations.

This module provides classes for defining molecular binding reactions,
including reaction angles and complete reaction rules. These components 
define how molecules can bind to each other during simulation, including 
geometric constraints and kinetic parameters.

Classes:
    ReactionAngleSet: Defines the angular constraints for binding reactions
        using the NERDSS software angle convention.
    ReactionRule: Complete definition of a binding reaction including
        reactant interfaces, kinetics, and geometric constraints.

Key Features:
    - NERDSS-compatible angle definitions
    - Serialization support for configuration files
    - Geometric constraint specification
    - Kinetic parameter storage
    - Direct interface type references

Example:
    ```python
    # Create interface types (assuming they exist)
    interface_a_b = InterfaceType(...)  # A->B interface
    interface_b_a = InterfaceType(...)  # B->A interface
    
    # Define reaction angles
    angles = ReactionAngleSet(
        theta1=0.0, theta2=0.0, phi1=0.0, phi2=0.0, omega=0.0
    )
    
    # Create reaction rule
    reaction = ReactionRule(
        expr="A + B -> A.B",
        reactant_interfaces=(interface_a_b, interface_b_a),
        bond_length_nm=1.5,
        angles=angles,
        ka=1e6,  # association rate
        kb=1e-3  # dissociation rate
    )
    
    # Access molecule types through interfaces
    mol_type_a = reaction.reactant_interfaces[0].this_mol_type
    mol_type_b = reaction.reactant_interfaces[1].this_mol_type
    ```

See Also:
    ionerdss.model.components.types: Molecule and interface type definitions
    ionerdss.simulation.kinetics: Kinetic simulation engine
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

# Import InterfaceType for type hints
from ionerdss.model.components.types import InterfaceType, MoleculeType

from ionerdss.utils.bond_geometry import compute_bond_angles_and_length_auto


@dataclass
class ReactionGeometrySet:
    """Defines angular and bond length constraints for molecular binding reactions.

    Uses the NERDSS software angle convention to specify the relative
    orientations required for successful binding between two molecules.

    Attributes:
        theta1: First polar angle in radians.
        theta2: Second polar angle in radians.
        phi1: First azimuthal angle in radians.
        phi2: Second azimuthal angle in radians.
        omega: Dihedral angle in radians.
    """
    # NERDSS software reaction angles
    theta1: float
    theta2: float
    phi1: float
    phi2: float
    omega: float

    # bond distance
    sigma_nm: float

    # normal vector
    norm1: np.ndarray
    norm2: np.ndarray

    def __init__(self, theta1: float, theta2: float,
                 phi1: float, phi2: float,
                 omega: float, sigma_nm: float,
                 norm1: np.ndarray, norm2: np.ndarray):
        """
        Simple assign bond angles and bond lengths
        """
        self.theta1 = theta1
        self.theta2 = theta2
        self.phi1 = phi1
        self.phi2 = phi2
        self.omega = omega
        self.sigma_nm = sigma_nm
        self.norm1 = norm1
        self.norm2 = norm2

    def as_list(self) -> List[float]:
        """Convert angles to a list of floats.

        Returns:
            List containing [theta1, theta2, phi1, phi2, omega, sigma].
        """
        return [self.theta1, self.theta2, self.phi1, self.phi2, self.omega, self.sigma_nm]

    def from_absolute_coords(self,
                             com1: np.ndarray,
                             com2: np.ndarray,
                             bind_site1: np.ndarray,
                             bind_site2: np.ndarray) -> "ReactionGeometrySet":
        """Calculate reaction angles from absolute coordinates and normal vectors.

        Computes the NERDSS angle set from the absolute positions and
        orientation vectors of two binding sites.

        Args:
            coord1: Absolute coordinates of first binding site.
            coord2: Absolute coordinates of second binding site.
            norm1: Normal vector of first binding site.
            norm2: Normal vector of second binding site.

        Returns:
            New ReactionAngleSet with calculated angles.

        Note:
            This is a placeholder implementation. The actual calculation
            would require the specific NERDSS angle convention formulas.
        """
        # actual NERDSS angle calculations
        # This would involve complex geometric calculations based on:
        # - Vector between binding sites: coord2 - coord1
        # - Normal vectors of both sites
        # - NERDSS-specific angle definitions
        theta1, theta2, phi1, phi2, omega, sigma, normal_point1, normal_point2 = compute_bond_angles_and_length_auto(com1, com2,
                                                                                                                     bind_site1, bind_site2)

        # For now, return current angles (this should be replaced with actual calculation)
        return ReactionGeometrySet(
            theta1,
            theta2,
            phi1,
            phi2,
            omega,
            sigma,
            normal_point1,
            normal_point2
        )


@dataclass
class ReactionRule:
    """Complete definition of a molecular binding reaction.

    Defines all aspects of a binding reaction including the reactant interfaces,
    geometric constraints, kinetic parameters, and requirements for the reaction
    to occur. Molecule types are accessed through the interface types.

    Attributes:
        expr: Human-readable reaction expression (e.g., "A + B -> A.B").
        reactant_interfaces: Tuple of two InterfaceType objects defining the 
            binding partners. Molecule types can be accessed via 
            interface.this_mol_type.
        required_free: Tuple of lists specifying interface names that must be
            unbound for each reactant. Defaults to empty lists.
        geometry: see ReactionGeometrySet above
        ka: Association rate constant. Defaults to 0.0.
        kb: Dissociation rate constant. Defaults to 0.0.
    """
    expr: str
    reactant_interfaces: Tuple[InterfaceType, InterfaceType]
    # list of interface names that are required to be free for the reaction to happen
    required_free: Tuple[List[str], List[str]] = None
    geometry: Optional[ReactionGeometrySet] = None
    ka: float = 0.0
    kb: float = 0.0

    def __post_init__(self):
        """Initialize default values after dataclass creation."""
        if self.required_free is None:
            self.required_free = ([], [])
        self.update_expr()

    def update_expr(self):
        """
        Update the BNGL reaction string with interface information
        """
        # Get molecule and interface names
        mol1_name = self.reactant_interfaces[0].this_mol_type_name
        mol2_name = self.reactant_interfaces[1].this_mol_type_name

        interface1_name = self.reactant_interfaces[0].get_name()
        interface2_name = self.reactant_interfaces[1].get_name()

        # Build reactant side with interface states and required_free constraints
        reactant1 = self.build_molecule_expression(
            mol1_name,
            interface1_name,
            "free",
            self.required_free[0] if self.required_free else []
        )

        reactant2 = self.build_molecule_expression(
            mol2_name,
            interface2_name,
            "free",
            self.required_free[1] if self.required_free else []
        )

        # Build product side with bound interfaces
        product1 = self.build_molecule_expression(
            mol1_name,
            interface1_name,
            "bound",
            self.required_free[0] if self.required_free else []
        )

        product2 = self.build_molecule_expression(
            mol2_name,
            interface2_name,
            "bound",
            self.required_free[1] if self.required_free else []
        )

        # Construct the full BNGL reaction expression
        self.expr = f"{reactant1} + {reactant2} <-> {product1}.{product2}"

    def build_molecule_expression(self, mol_name: str, binding_interface: str,
                                   binding_state: str, required_free_interfaces: List[str]) -> str:
        """
        Build a BNGL molecule expression with interface states.

        Args:
            mol_name: Name of the molecule
            binding_interface: Name of the primary binding interface
            binding_state: State of the binding interface ("free" or "bound")
            required_free_interfaces: List of interfaces that must be free

        Returns:
            BNGL molecule expression string
        """
        interfaces = []

        # Add the primary binding interface
        if binding_state == "free":
            interfaces.append(f"{binding_interface}")
        else:  # bound
            interfaces.append(f"{binding_interface}!1")

        # Add required free interfaces if any
        for free_interface in required_free_interfaces:
            if free_interface != binding_interface:  # Don't duplicate the binding interface
                interfaces.append(free_interface)

        # Construct the molecule expression
        if interfaces:
            interface_str = "(" + ",".join(interfaces) + ")"
            return f"{mol_name}{interface_str}"
        else:
            return mol_name

    @property
    def reactant_molecule_types(self) -> Tuple[Optional['MoleculeType'], Optional['MoleculeType']]:
        """Get the molecule types for both reactant interfaces.

        Returns:
            Tuple of MoleculeType objects (or None if not set) for the
            two reactant interfaces.
        """
        return (
            self.reactant_interfaces[0].this_mol_type,
            self.reactant_interfaces[1].this_mol_type
        )

    def get_reactant_interface_names(self) -> Tuple[str, str]:
        """Get the names of both reactant interfaces.

        Returns:
            Tuple of interface names for the two reactants.
        """
        return (
            self.reactant_interfaces[0].get_name(),
            self.reactant_interfaces[1].get_name()
        )

    def to_dict(self) -> dict:
        """
        output a dictionary the content
        """
        result = {
            'expr': self.expr,
            'reactant_interfaces': [
                self.reactant_interfaces[0].get_name(),
                self.reactant_interfaces[1].get_name()
            ],
            'required_free': [
                list(self.required_free[0]),
                list(self.required_free[1])
            ],
            'ka': self.ka,
            'kb': self.kb
        }

        # Add geometry information if present
        if self.geometry is not None:
            result['geometry'] = {
                'theta1': self.geometry.theta1,
                'theta2': self.geometry.theta2,
                'phi1': self.geometry.phi1,
                'phi2': self.geometry.phi2,
                'omega': self.geometry.omega,
                'sigma_nm': self.geometry.sigma_nm,
                'norm1': self.geometry.norm1.tolist(),  # Convert numpy array to list
                'norm2': self.geometry.norm2.tolist()   # Convert numpy array to list
            }
        else:
            result['geometry'] = None

        return result
