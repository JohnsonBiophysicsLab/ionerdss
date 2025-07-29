import numpy as np

def compute_diffusion_constants_nm_us(radius_nm,
                                      temperature_kelvin=298.15,
                                      viscosity_pas=8.9e-4):
    """
    Compute translational and rotational diffusion constants based on
    einstein stokes equation.
    
    Args:
        radius_nm (float): radius in nanometers
        temperature_kelvin (float): temperature in Kelvin
        viscosity_pas (float): viscosity in Pa s (default: water)
    
    Returns:
        tuple: (translational diffusion constant in nm^2/μs,
        rotational diffusion constant in rad^2/μs)
    """
    if radius_nm <= 0:
        raise ValueError("Radius must be positive!")
    if temperature_kelvin <= 0:
        raise ValueError("Temperture (in Kelvin) must be positive!")
    if viscosity_pas <= 0:
        raise ValueError("Viscosity (in Pa s) must be positive!")

    boltzmann_constant = 1.380649e-23  # J/K
    radius_m = radius_nm * 1e-9  # convert nm to meters

    # Diffusion constants
    d_t_m2_per_s = boltzmann_constant * temperature_kelvin / (6 * np.pi * viscosity_pas * radius_m)
    d_r_rad2_per_s = boltzmann_constant * temperature_kelvin / (8 * np.pi * viscosity_pas * radius_m**3)

    # Convert units
    translational_diffusion_constant_nm2_per_us = d_t_m2_per_s * 1e12  # m²/s → nm²/μs
    rotational_diffusion_constant_rad2_per_us = d_r_rad2_per_s * 1e-6  # rad²/s → rad²/μs

    return translational_diffusion_constant_nm2_per_us, rotational_diffusion_constant_rad2_per_us
