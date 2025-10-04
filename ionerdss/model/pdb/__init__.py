"""
ionerdss.model.pdb.__init__.py

PDB to NERDSS parameter pipeline for molecular simulation model generation.

This package converts all-atom protein complexes from mmCIF/PDB files into 
coarse-grained representations suitable for NERDSS modeling. The pipeline
includes structure parsing, interface detection, chain grouping, template
generation, and system assembly.

Modules:
    hyperparameters: Simulation hyperparameters and configuration
    parser: PDB/mmCIF file parsing and structure processing
    coarse_graining: Interface detection and coarse-grained representation
    chain_grouping: Repeated chain detection and classification
    template_builder: Molecular and interface template generation
    system_builder: Final system assembly and validation

Key Features:
    - Automatic interface detection via KD-tree spatial queries
    - Multiple chain grouping strategies (header, sequence, structure)
    - Template-based regularization across symmetry mates
    - Steric clash detection for mutually exclusive interfaces
    - Complete system serialization support

Example Usage:
    ```python
    from ionerdss.model.pdb import PDBModelBuilder
    
    # Build system from PDB file
    builder = PDBModelBuilder("1ABC.pdb")
    system = builder.build_system(
        workspace_path="/path/to/workspace",
        distance_cutoff=0.6,  # nm
        matching_mode="default"
    )
    
    # Save complete system
    system.to_json("1ABC_system.json")
    ```
"""

from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser
from .coarse_graining import CoarseGrainer
from .chain_grouping import ChainGrouper
from .template_builder import TemplateBuilder
from .system_builder import SystemBuilder
from .main import PDBModelBuilder

__all__ = [
    'PDBModelHyperparameters',
    'PDBParser', 
    'CoarseGrainer',
    'ChainGrouper',
    'TemplateBuilder',
    'SystemBuilder',
    'PDBModelBuilder'
]
