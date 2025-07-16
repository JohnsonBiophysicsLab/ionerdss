"""
reaction.py

Constructs NERDSS-compatible ReactionTemplate objects that describe how protein interfaces
interact, bind, or transition based on coarse-grained structural data.

This module parses:
- Regularized model output (with chain/interface geometry)
- Interface pairing logic (based on proximity, energy, or repeated_chainy)
- Optionally energy thresholds to define reaction formation

Each reaction connects two InterfaceTemplate objects from different MoleculeTemplates
and defines a placeholder or computed reaction rate constant.

Functions
---------
- build_binding_reactions(model: dict, interface_templates: list, energy_cutoff: float)
    Generate ReactionTemplate objects between interface pairs meeting energy/geometry criteria.
"""

from ionerdss.model.components import ReactionTemplate

def build_binding_reactions(model, interface_templates, energy_cutoff=0.0):
    """
    Create binding reactions from the coarse-grain model.

    Parameters
    ----------
    model : dict
        Output of `regularize_model()` including 'interfaces', 'interface_energies', etc.
    interface_templates : list of InterfaceTemplate
        One-to-one list of interface objects per chain, structured as:
            interface_templates[i][j] == j-th interface of chain i
    energy_cutoff : float, optional
        Only form reactions between interfaces with energy below this threshold (e.g., -1.0).

    Returns
    -------
    reactions : list of ReactionTemplate
        Reactions formed between interface pairs.
    """
    reactions = []
    N = len(model["interfaces"])

    for i in range(N):
        for j, neighbor in enumerate(model["interfaces"][i]):
            # Get current interface and neighbor interface
            iface1 = interface_templates[i][j]

            # Identify back index in neighbor's interface list
            try:
                k = model["interfaces"][neighbor].index(i)
                iface2 = interface_templates[neighbor][k]
            except ValueError:
                continue  # Unidirectional or asymmetric entry; skip

            # Prevent duplication: always define (i < neighbor)
            if i < neighbor:
                energy = model["interface_energies"][i][j]
                if energy <= energy_cutoff:
                    reaction = ReactionTemplate(
                        reactants=[iface1, iface2],
                        products=[],  # Filled later in simulation setup
                        rate=1.0,     # Placeholder
                        metadata={"energy": energy}
                    )
                    reactions.append(reaction)

    return reactions
