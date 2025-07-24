# ionerdss/__init__.py

# setup logging level
# TODO: Take user input to setup logging level. 
# This might be done when creating Analysis instance.
import logging

logging.basicConfig(level=logging.WARNING)

class LazyLoader:
    """Transparent lazy module loader"""
    def __init__(self, module_path, attribute=None):
        self._module_path = module_path
        self._attribute = attribute
        self._loaded_module = None
        self._loaded_attribute = None
    
    def __call__(self, *args, **kwargs):
        """Support calling the lazy-loaded object directly"""
        if self._loaded_attribute is None:
            if self._loaded_module is None:
                import importlib
                self._loaded_module = importlib.import_module(self._module_path, package="ionerdss")
            
            self._loaded_attribute = getattr(self._loaded_module, self._attribute) if self._attribute else self._loaded_module
        
        if callable(self._loaded_attribute):
            return self._loaded_attribute(*args, **kwargs)
        raise TypeError(f"{self._attribute} is not callable")
            
    def __getattr__(self, name):
        """Access attributes of the lazy-loaded module/object"""
        if self._loaded_module is None:
            import importlib
            self._loaded_module = importlib.import_module(self._module_path, package="ionerdss")
        
        if self._attribute:
            if self._loaded_attribute is None:
                self._loaded_attribute = getattr(self._loaded_module, self._attribute)
            return getattr(self._loaded_attribute, name)
        return getattr(self._loaded_module, name)
    
    # def __dir__(self):
    #     """
    #     Provide attributes for autocompletion by parsing __all__.
        
    #     This method is called by dir() and IDEs. It attempts to parse
    #     __all__ from the source file for speed. If that fails or is not
    #     available, it falls back to importing the module.
    #     """
    #     if self._attribute:

    #     else:
    #         if self._loaded_module is None:
    #             import importlib
    #             self._loaded_module = importlib.import_module(self._module_path, package="ionerdss")
        
    #     return 

# Lazily load core classes
Model = LazyLoader('.model.model', 'Model')
MoleculeType = LazyLoader('.model.model', 'MoleculeType')
MoleculeInterface = LazyLoader('.model.model', 'MoleculeInterface')
ReactionType = LazyLoader('.model.model', 'ReactionType')
Coords = LazyLoader('.model.coords', 'Coords')
PDBModel = LazyLoader('.model.pdb_model', 'PDBModel')
DesignModel = LazyLoader('.model.design_model', 'DesignModel')
PlatonicSolid = LazyLoader('.model.PlatonicSolids', 'PlatonicSolid')
# Standard, preferred name
generate_ode_model_from_pdb = LazyLoader('.model.complex', 'generate_ode_model_from_pdb')
# Backward-compatible alias
ParseComplexes = generate_ode_model_from_pdb
ReactionStringParser = LazyLoader('.ode_solver.reaction_string_parser', 'ReactionStringParser')
solve_reaction_ode = LazyLoader('.ode_solver.reaction_ode_solver', 'solve_reaction_ode')
reaction_dydt = LazyLoader('.ode_solver.reaction_ode_solver', 'dydt')
calculate_macroscopic_reaction_rates = LazyLoader('.ode_solver.reaction_ode_solver', 'calculate_macroscopic_reaction_rates')
# Gillespie simulation related
SimpleGillespie = LazyLoader('.gillespie_simulation.simple_gillespie')
AdaptiveRates = LazyLoader('.gillespie_simulation.adaptive_rates')
# GUI
gui = LazyLoader('.nerdss_guis.gui', 'gui')
pdb_gui = LazyLoader('.nerdss_guis.nerdss', 'nerdss')
cube_face = LazyLoader('.model.platonic_solids.cube.cube_face', 'cube_face')
cube_vert = LazyLoader('.model.platonic_solids.cube.cube_vert', 'cube_vert')
dode_face = LazyLoader('.model.platonic_solids.dode.dode_face', 'dode_face')
dode_vert = LazyLoader('.model.platonic_solids.dode.dode_vert', 'dode_vert')
icos_face = LazyLoader('.model.platonic_solids.icos.icos_face', 'icos_face')
icos_vert = LazyLoader('.model.platonic_solids.icos.icos_vert', 'icos_vert')
octa_face = LazyLoader('.model.platonic_solids.octa.octa_face', 'octa_face')
octa_vert = LazyLoader('.model.platonic_solids.octa.octa_vert', 'octa_vert')
tetr_face = LazyLoader('.model.platonic_solids.tetr.tetr_face', 'tetr_face')
tetr_vert = LazyLoader('.model.platonic_solids.tetr.tetr_vert', 'tetr_vert')
convert_simularium = LazyLoader('.simularium_converter.simularium_converter', 'convert_simularium')


# Lazily load simulation and analysis modules
Simulation = LazyLoader('.nerdss_simulation.simulation', 'Simulation')
Analysis = LazyLoader('.nerdss_analysis.analysis', 'Analysis')
DataIO = LazyLoader('.nerdss_analysis.data_readers', 'DataIO')

def configure_plotting():
    """Configure plotting styles - only call this when you're ready to plot."""
    import seaborn as sns
    fontsize = 12
    sns.set_style("ticks")
    sns.set_context("paper", rc={
        "font.size": fontsize,
        "axes.titlesize": fontsize,
        "axes.labelsize": fontsize,
        "xtick.labelsize": fontsize,
        "ytick.labelsize": fontsize,
        "legend.fontsize": fontsize,
        "font.family": "serif"
    })

# Version information - only computed when __version__ is accessed
def _get_version():
    try:
        import pkg_resources
        return pkg_resources.get_distribution("ioNERDSS").version
    except:
        return "unknown"

# Define __version__ as a property
class VersionProperty:
    def __get__(self, obj, objtype=None):
        return _get_version()

class _ModuleProperties:
    __version__ = VersionProperty()

# Apply the properties
import sys
sys.modules[__name__].__class__ = type("ionerdss", (sys.modules[__name__].__class__, _ModuleProperties), {})