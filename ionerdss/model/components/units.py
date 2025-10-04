"""
ionerdss.model.components.units

Unit definitions and conversions for molecular simulation parameters.

This module provides standardized unit definitions used throughout the ionerdss
simulation framework. It ensures consistency in physical units across different
components and provides conversion utilities when needed.

The standard units used in ionerdss differ from some common file formats:
- Coordinates: nanometers (nm) - Note: PDB files use Angstroms (Å), 1 nm = 10 Å
- Binding radius: nanometers (nm)
- Translational diffusion: nm²/μs
- Rotational diffusion: rad²/μs

Classes:
    Units: Dataclass that defines and tracks the units used for different
        physical quantities in the simulation.

Key Features:
    - Standardized unit definitions across the simulation framework
    - Serialization support for configuration files
    - Default values matching ionerdss conventions
    - Clear documentation of unit differences from common file formats

Unit Conversion Notes:
    When importing data from external sources, be aware of unit differences:
    - PDB files: coordinates in Angstroms (Å) → convert to nm by dividing by 10
    - Some diffusion data: may be in different time units (s vs μs)
    - Angles: ensure rotational quantities are in radians, not degrees

Example:
    ```python
    # Create default units
    units = Units()
    print(units.coords)  # "nm"
    
    # Create custom units
    custom_units = Units(
        coords="Å",  # For PDB compatibility
        D_trans="nm^2/s"  # Different time scale
    )
    
    # Serialize for configuration
    config = units.to_dict()
    restored_units = Units.from_dict(config)
    ```

See Also:
    ionerdss.utils.conversions: Unit conversion utilities (if available)
    ionerdss.io.pdb: PDB file handling with automatic unit conversion
"""

from dataclasses import dataclass

@dataclass
class Units:
    """Defines the units used for physical quantities in molecular simulations.
    
    This class standardizes the units used throughout the ionerdss framework
    and provides serialization capabilities for configuration management.
    
    Important: The default units in ionerdss use nanometers for length scales,
    which differs from PDB files that use Angstroms (1 nm = 10 Å).
    
    Attributes:
        coords: Unit for spatial coordinates. Default is "nm" (nanometers).
            Note: PDB files use Angstroms - conversion required.
        binding_radius: Unit for binding interaction distances. Default is "nm".
        D_trans: Unit for translational diffusion constants. Default is "nm^2/µs".
        D_rot: Unit for rotational diffusion constants. Default is "rad^2/µs".
    """

    coords: str = "nm"
    binding_radius: str = "nm"
    D_trans: str = "nm^2/µs"
    D_rot: str = "rad^2/µs"

    def to_dict(self) -> dict:
        """Convert the Units to a dictionary representation.
        
        Serializes the units data to a dictionary format suitable for
        JSON export, configuration files, or storage.
        
        Returns:
            Dictionary containing the units data with keys: 'coords',
            'binding_radius', 'diffusion_translation', and 'diffusion_rotation'.
            
        Note:
            The keys 'diffusion_translation' and 'diffusion_rotation' are used
            in the dictionary format for consistency with other components,
            while the class attributes use shorter names.
        """
        return {
            "coords": self.coords,
            "binding_radius": self.binding_radius,
            "diffusion_translation": self.D_trans,
            "diffusion_rotation": self.D_rot,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Units":
        """Create a Units instance from a dictionary.
        
        Deserializes units data from a dictionary format, typically
        loaded from JSON configuration files or other storage formats.
        Missing fields are filled with default values.
        
        Args:
            d: Dictionary containing units data. All keys are optional.
                Supported keys: 'coords', 'binding_radius', 
                'diffusion_translation', 'diffusion_rotation'.
                If d is None or empty, returns default Units instance.
        
        Returns:
            New Units instance populated with the dictionary data or defaults.
            
        Example:
            ```python
            config = {
                "coords": "Å",
                "diffusion_translation": "nm^2/s"
            }
            units = Units.from_dict(config)
            # units.coords == "Å"
            # units.binding_radius == "nm" (default)
            ```
        """
        if not d:
            return cls()
        return cls(
            coords=d.get("coords", "nm"),
            binding_radius=d.get("binding_radius", "nm"),
            D_trans=d.get("diffusion_translation", "nm^2/µs"),
            D_rot=d.get("diffusion_rotation", "rad^2/µs"),
        )
