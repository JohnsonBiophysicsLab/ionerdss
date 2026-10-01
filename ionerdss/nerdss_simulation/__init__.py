# ionerdss/nerdss_simulation/__init__.py

from .simulation import Simulation
from .executable import resolve_nerdss_executable

# >>>>>> For autocompletion >>>>>>
__all__ = [
    'Simulation',
    'resolve_nerdss_executable',
]
def __dir__():
    return __all__
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<