"""
ionerdss: A user-friendly toolkit for setting up
NERDSS simulations and analyzing results.
================================================

Documentation is available in the docstrings and
online at https://ionerdss.readthedocs.io/en/

Subpackages
-----------
::

 Simulation                   --- Main class for running simulations.
 Analysis                     --- Main class for analyzing simulation data.
 Model                        --- The core model class for defining the system.
 MoleculeType                 --- Defines a type of molecule in the model.
 MoleculeInterface            --- Defines the binding interface for a molecule.
 ReactionType                 --- Defines a type of reaction in the model.
 Coords                       --- Represents 3D coordinates.
 PDBModel                     --- Creates a model from a PDB file.
 DesignModel                  --- A model for designing molecular structures.
 PlatonicSolid                --- Class for generating platonic solid geometries.
 generate_ode_model_from_pdb  --- Generates an ODE model from PDB complexes.
 ParseComplexes               --- Alias for generate_ode_model_from_pdb.
 ReactionStringParser         --- Parses reaction definitions from a string.
 solve_reaction_ode           --- Solves reaction kinetics using Ordinary Differential Equations (ODEs).
 reaction_dydt                --- The rate-of-change function (dy/dt) for the ODE solver.
 calculate_macroscopic_reaction_rates --- Calculates macroscopic reaction rates from microscopic parameters.
 SimpleGillespie              --- Implements the Gillespie stochastic simulation algorithm (SSA).
 AdaptiveRates                --- Implements adaptive rate constants for simulations.
 gui                          --- Launches the main graphical user interface.
 pdb_gui                      --- A specific GUI for PDB file manipulation and viewing.
 cube_face                    --- Component class for a cube face.
 cube_vert                    --- Component class for a cube vertex.
 dode_face                    --- Component class for a dodecahedron face.
 dode_vert                    --- Component class for a dodecahedron vertex.
 icos_face                    --- Component class for an icosahedron face.
 icos_vert                    --- Component class for an icosahedron vertex.
 octa_face                    --- Component class for an octahedron face.
 octa_vert                    --- Component class for an octahedron vertex.
 tetr_face                    --- Component class for a tetrahedron face.
 tetr_vert                    --- Component class for a tetrahedron vertex.
 convert_simularium           --- Converts simulation data to the Simularium format.
 DataIO                       --- Handles reading and writing of simulation data.

Public API in the main ionerdss namespace
-----------------------------------------
::
 
 __version__       --- SciPy version string

"""

# >>>>>>>>>>>>>>>> setup logging level >>>>>>>>>>>>>>>> 
# TODO: Take user input to setup logging level. 
# This might be done when creating Analysis instance.
import logging as _logging
_logging.basicConfig(level=_logging.WARNING)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>> Version information >>>>>>>>>>>>>>>>
try:
    import pkg_resources
    __version__ = pkg_resources.get_distribution("ioNERDSS").version
except:
    __version__ = "unknown (This might be a local copy)"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


import importlib as _importlib


# A mapping from public API names to their internal module locations.
# Structure:
#     'PublicAPIName': ['.internal.module.path', 'ClassName']
submodules = {
    'Model': ['.model.model', 'Model'],
    'MoleculeType': ['.model.model', 'MoleculeType'],
    'MoleculeInterface': ['.model.model', 'MoleculeInterface'],
    'ReactionType': ['.model.model', 'ReactionType'],
    'Coords': ['.model.coords', 'Coords'],
    'PDBModel': ['.model.pdb_model', 'PDBModel'],
    'DesignModel': ['.model.design_model', 'DesignModel'],
    'PlatonicSolid': ['.model.PlatonicSolids', 'PlatonicSolid'],
    'generate_ode_model_from_pdb': ['.model.complex', 'generate_ode_model_from_pdb'],
    'ParseComplexes': ['.model.complex', 'generate_ode_model_from_pdb'],
    'ReactionStringParser': ['.ode_solver.reaction_string_parser', 'ReactionStringParser'],
    'solve_reaction_ode': ['.ode_solver.reaction_ode_solver', 'solve_reaction_ode'],
    'reaction_dydt': ['.ode_solver.reaction_ode_solver', 'dydt'],
    'calculate_macroscopic_reaction_rates': ['.ode_solver.reaction_ode_solver', 'calculate_macroscopic_reaction_rates'],
    'SimpleGillespie': ['.gillespie_simulation.simple_gillespie', ''],
    'AdaptiveRates': ['.gillespie_simulation.adaptive_rates', ''],
    'gui': ['.nerdss_guis.gui', 'gui'],
    'pdb_gui': ['.nerdss_guis.nerdss', 'nerdss'],
    'cube_face': ['.model.platonic_solids.cube.cube_face', 'cube_face'],
    'cube_vert': ['.model.platonic_solids.cube.cube_vert', 'cube_vert'],
    'dode_face': ['.model.platonic_solids.dode.dode_face', 'dode_face'],
    'dode_vert': ['.model.platonic_solids.dode.dode_vert', 'dode_vert'],
    'icos_face': ['.model.platonic_solids.icos.icos_face', 'icos_face'],
    'icos_vert': ['.model.platonic_solids.icos.icos_vert', 'icos_vert'],
    'octa_face': ['.model.platonic_solids.octa.octa_face', 'octa_face'],
    'octa_vert': ['.model.platonic_solids.octa.octa_vert', 'octa_vert'],
    'tetr_face': ['.model.platonic_solids.tetr.tetr_face', 'tetr_face'],
    'tetr_vert': ['.model.platonic_solids.tetr.tetr_vert', 'tetr_vert'],
    'convert_simularium': ['.simularium_converter.simularium_converter', 'convert_simularium'],
    'Simulation': ['.nerdss_simulation', 'Simulation'],
    'Analyzer': ['.analysis', 'Analyzer'],
}

__all__ = list(submodules.keys()) + [
    '__version__',
]

def __dir__():
    return __all__

_module_loaded = {}

def _get_module(module_path):
    if module_path in _module_loaded:
        module = _module_loaded[module_path]
    else:
        module = _importlib.import_module(module_path, package='ionerdss')
        _module_loaded.update({module_path:module})
    return module

def __getattr__(name):
    if name in submodules:
        _module_path, _attribute = submodules[name]
        if bool(_attribute):
            _module = _get_module(_module_path)
            return getattr(_module, _attribute)
        else:
            return _get_module(_module_path)
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(
                f"Module 'ionerdss' has no attribute '{name}'"
            )

