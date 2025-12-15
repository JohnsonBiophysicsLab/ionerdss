"""
ODE solver module for reaction kinetics.

This module provides tools for parsing reaction strings and solving
ordinary differential equations for chemical reaction systems.
"""

from .reaction_string_parser import ReactionStringParser
from .reaction_ode_solver import (
    dydt,
    solve_reaction_ode,
    calculate_macroscopic_reaction_rates
)

__all__ = [
    'ReactionStringParser',
    'dydt',
    'solve_reaction_ode',
    'calculate_macroscopic_reaction_rates'
]
