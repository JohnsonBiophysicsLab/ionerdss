"""Components module - provides compatibility layer for pdb_model imports."""

# Import from the model.py file to provide backward compatibility
from ..model import MoleculeType, MoleculeInterface, ReactionType, Model

# These helper classes are defined in pdb_model.py itself
# We'll need to handle them separately
__all__ = ['MoleculeType', 'MoleculeInterface', 'ReactionType', 'Model']
