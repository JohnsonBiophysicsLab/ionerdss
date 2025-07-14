# model/pdb/coarse_grain.py
"""
Here is a detailed documentation explaining the logic, steps, assumptions, and output structure of the `coarse_grain_structure` function. You can use this for both internal documentation and as a reference in project wikis or READMEs.

---

## Documentation: Coarse-Graining Protein Structure for NERDSS Modeling

### Overview

The `coarse_grain_structure()` function takes a Biopython `Structure` object (typically parsed from a `.cif` or `.pdb` file) and computes a **coarse-grained molecular representation** for use in NERDSS. This includes:

* Center-of-mass (COM) for each chain.
* Approximate molecular radii.
* Detection of **binding interfaces** between chains.
* Interface coordinates, residues, and pairwise contact energies.

This step is a foundational part of the NERDSS modeling pipeline that enables downstream creation of `Molecule`, `Interface`, and `Reaction` templates.

---

### Inputs

```python
coarse_grain_structure(structure, distance_cutoff=0.35, residue_cutoff=3, options=None)
```

* **`structure`**: A Biopython `Structure` object, parsed using `Bio.PDB.MMCIFParser` or `PDBParser`.
* **`distance_cutoff`** (`float`): Max atom-atom distance (in nanometers) for two residues to be considered interacting. Default is `0.35 nm` (i.e., 3.5 Å).
* **`residue_cutoff`** (`int`): Minimum number of contacting residues required to define a valid interface between two chains.
* **`options`** (`dict` or `None`): Reserved for future use (e.g., plotting flags).

---

### Steps Performed

#### 1. **Chain Identification and Filtering**

* All chains with at least one standard amino acid (`is_aa(res)`) are included.
* Chains are sorted by chain ID.

#### 2. **Center-of-Mass and Radius Calculation**

* For each chain:

  * Collect all atomic coordinates (across all amino acids).
  * Compute the average as the **COM**.
  * Compute RMS deviation from the COM as the **molecular radius**.

#### 3. **Bounding Box Pre-check**

* For each pair of chains, 3D bounding boxes are computed.
* Pairs whose bounding boxes are disjoint (with margin) are skipped early for efficiency.

#### 4. **Interface Detection**

* For each valid pair of chains:

  * All Cα atoms are extracted from each chain.
  * A KD-tree is built on chain B’s atoms to query neighbors within `distance_cutoff × 10` from each atom in chain A.
  * If enough residue pairs are found:

    * Define a **binding interface** for each chain (as the average of contact residue coordinates).
    * Record the list of contacting residues.
    * Look up interaction energies for each residue pair using the `energy_table` (e.g., ARG–GLU → −3.5 kcal/mol).
    * Store total contact energy.

---

### Output

Returns a dictionary:

```python
{
    'chains': List[Chain],                   # Biopython Chain objects
    'COMs': List[Coords],                    # Center of mass for each chain
    'radii': List[float],                    # RMS radius per chain
    'interfaces': List[List[str]],           # Neighbor chain IDs per chain
    'interface_coords': List[List[Coords]],  # Interface COMs per chain
    'interface_residues': List[List[int]],   # Residue IDs per interface
    'interface_energies': List[List[float]]  # Interaction energy per interface
}
```

Each list is ordered per-chain, with sublists per interface (i.e., indexed by interacting partner).

---

### Assumptions

* Cα atoms are representative of residue positions (i.e., no explicit side-chain modeling).
* Interface detection is symmetric: if A binds B, both A and B get interfaces.
* Residue-residue energies are additive and symmetric (based on a simplified energy table).
* Chain IDs are assumed to be unique within the structure.

---

### Example Use

```python
from Bio.PDB import MMCIFParser
from model.pdb.coarse_grain import coarse_grain_structure

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("TEST", "1abc.cif")
result = coarse_grain_structure(structure)

for i, com in enumerate(result['COMs']):
    print(f"Chain {result['chains'][i].id} COM: {com}")
```

---

### Possible Extensions

* Switch to residue centers instead of Cα atoms for more accuracy.
* Add electrostatic or hydrophobic scoring.
* Expose options for including DNA, RNA, or ligands.
* Visualize detected interfaces and residues using PyMOL.

---

Let me know if you'd like a graphical diagram or inline docstring version of this documentation.

"""
import numpy as np
from scipy.spatial import KDTree
from Bio.PDB.Polypeptide import is_aa
from ionerdss.math.coords import Coords  # Assuming your Coords class is already implemented
from .energy_table import get_default_energy_table

def coarse_grain_structure(structure, distance_cutoff=0.35, residue_cutoff=3, options=None):
    """
    Analyze a Biopython structure to identify chain COMs and binding interfaces.

    Args:
        structure : Bio.PDB.Structure.Structure
        distance_cutoff : float
            Max distance (nm) for atoms to be considered interacting.
        residue_cutoff : int
            Min number of contacting residues to define an interface.
        options : dict or None
            Optional config (e.g., plotting).

    Returns:
        dict : {
            'chains': list of chains,
            'COMs': list of Coords,
            'interfaces': list of interface data per chain,
        }
    """
    chains = sorted([chain for chain in structure.get_chains()
                     if any(is_aa(res) for res in chain.get_residues())],
                    key=lambda c: c.id)
    num_chains = len(chains)
    COMs = []
    chain_radii = []
    interfaces = [[] for _ in range(num_chains)]
    interface_coords = [[] for _ in range(num_chains)]
    interface_residues = [[] for _ in range(num_chains)]
    interface_energies = [[] for _ in range(num_chains)]

    energy_table = get_default_energy_table()

    # Precompute COMs and radii
    for chain in chains:
        atoms = [atom.coord for res in chain for atom in res if is_aa(res)]
        if not atoms:
            COMs.append(None)
            chain_radii.append(0.0)
            continue
        atoms = np.array(atoms)
        avg = atoms.mean(axis=0)
        COMs.append(Coords(*avg))
        radius = np.sqrt(np.mean(np.sum((atoms - avg) ** 2, axis=1)))
        chain_radii.append(radius)

    # Bounding boxes for fast exclusion
    def bounding_box(chain):
        coords = np.array([atom.coord for res in chain for atom in res if is_aa(res)])
        return coords.min(axis=0), coords.max(axis=0) if len(coords) else (None, None)

    boxes = [bounding_box(c) for c in chains]

    # Loop over chain pairs
    for i in range(num_chains - 1):
        for j in range(i + 1, num_chains):
            min1, max1 = boxes[i]
            min2, max2 = boxes[j]
            if min1 is None or min2 is None:
                continue
            if np.any(min2 > max1 + distance_cutoff * 10) or np.any(max2 < min1 - distance_cutoff * 10):
                continue

            chain_i = chains[i]
            chain_j = chains[j]
            residues_i = [(res.id[1], res.get_resname().upper(), res['CA'].coord)
                          for res in chain_i if is_aa(res) and 'CA' in res]
            residues_j = [(res.id[1], res.get_resname().upper(), res['CA'].coord)
                          for res in chain_j if is_aa(res) and 'CA' in res]

            if not residues_i or not residues_j:
                continue

            atoms_j = np.array([r[2] for r in residues_j])
            tree = KDTree(atoms_j)
            results = tree.query_ball_point([r[2] for r in residues_i], r=distance_cutoff * 10)

            iface_i_ids, iface_i_coords, iface_i_types = [], [], []
            iface_j_ids, iface_j_coords, iface_j_types = [], [], []
            residue_pairs = {}

            for idx_i, neighbors in enumerate(results):
                if neighbors:
                    res_i_id, res_i_type, ca_i = residues_i[idx_i]
                    if res_i_id not in iface_i_ids:
                        iface_i_ids.append(res_i_id)
                        iface_i_coords.append(ca_i)
                        iface_i_types.append(res_i_type)
                    for idx_j in neighbors:
                        res_j_id, res_j_type, ca_j = residues_j[idx_j]
                        if res_j_id not in iface_j_ids:
                            iface_j_ids.append(res_j_id)
                            iface_j_coords.append(ca_j)
                            iface_j_types.append(res_j_type)
                        key = (res_i_type, res_j_type)
                        residue_pairs[(res_i_id, res_j_id)] = energy_table.get(key, 0.0)

            total_energy = sum(residue_pairs.values())
            if len(iface_i_ids) >= residue_cutoff and len(iface_j_ids) >= residue_cutoff:
                COM_i = Coords(*np.mean(iface_i_coords, axis=0))
                COM_j = Coords(*np.mean(iface_j_coords, axis=0))

                interfaces[i].append(chain_j.id)
                interface_coords[i].append(COM_i)
                interface_residues[i].append(sorted(iface_i_ids))
                interface_energies[i].append(total_energy)

                interfaces[j].append(chain_i.id)
                interface_coords[j].append(COM_j)
                interface_residues[j].append(sorted(iface_j_ids))
                interface_energies[j].append(total_energy)

    return {
        'chains': chains,
        'COMs': COMs,
        'radii': chain_radii,
        'interfaces': interfaces,
        'interface_coords': interface_coords,
        'interface_residues': interface_residues,
        'interface_energies': interface_energies,
    }
