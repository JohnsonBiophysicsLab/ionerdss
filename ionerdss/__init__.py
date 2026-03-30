"""
ionerdss: A user-friendly toolkit for setting up
NERDSS simulations and analyzing results.
================================================

Documentation is available in the docstrings and
online at https://johnsonbiophysicslab.github.io/ionerdss/
 
 __version__       --- SciPy version string

"""

# >>>>>>>>>>>>>>>> setup logging level >>>>>>>>>>>>>>>> 
# TODO: Take user input to setup logging level. 
# This might be done when creating Analysis instance.
from __future__ import annotations

import logging as _logging
_logging.basicConfig(level=_logging.WARNING)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>> Version information >>>>>>>>>>>>>>>>
try:
    from importlib.metadata import version
    __version__ = version("ioNERDSS")
except:
    __version__ = "unknown (This might be a local copy)"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


import importlib as _importlib

# A mapping from public API names to their internal module locations.
# Structure:
#     'PublicAPIName': ['.internal.module.path', 'ClassName']
submodules = {
    'System': ['.model.components.system', 'System'],
    'platonic_solid_generator': ['.model.PlatonicSolids', 'PlatonicSolidsModel'],
    'build_system_from_plat': ['.model.PlatonicSolids', 'build_system_from_plat'],
    'convert_simularium': ['.simularium_converter.simularium_converter', 'convert_simularium'],
    'Simulation': ['.nerdss_simulation', 'Simulation'],
    'Analyzer': ['.analysis', 'Analyzer'],
    'ODEPipelineConfig': ['.ode.ode_pipeline', 'ODEPipelineConfig'],
    'run_ode_pipeline': ['.ode.ode_pipeline', 'run_ode_pipeline'],
    'build_system_from_pdb': ['.api', 'build_system_from_pdb'],
    'prepare_structure_validation_for_system': ['.api', 'prepare_structure_validation_for_system'],
    'compare_structure_to_design': ['.api', 'compare_structure_to_design'],
    'StructureValidationConfig': ['.model.pdb.structure_validation', 'StructureValidationConfig'],
    'StructureValidationArtifacts': ['.model.pdb.structure_validation', 'StructureValidationArtifacts'],
    'StructureAlignmentResult': ['.model.pdb.structure_validation', 'StructureAlignmentResult'],
    'visualize_trajectory_ovito': ['.ovito_visualizer', 'visualize_trajectory_ovito']
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
