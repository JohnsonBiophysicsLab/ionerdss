# `ionerdss.model.pdb` Module Overview

## **How reaction / simulation systems are built up from a PDB or mmCIF file**

### Pipeline Step: `coarse_grain_structure(path_to_structure)`

This function (from `coarse_grain.py`) processes an input `.pdb` or `.cif` file to detect **interfaces** using spatial geometry and contact scoring:

### Detection Process:

* **Load structure** using `Bio.PDB.MMCIFParser` or `PDBParser`
* **Identify all chains** and compute each chain’s:

  * **Center of mass (COM)**
  * **Radius** (based on heavy atom spread or Cα)
* **Detect interfaces between chains** by:

  * Calculating **minimum inter-chain residue-residue distances**
  * Defining a **pair as interacting** if residues from two chains are within a cutoff (e.g., 6 Å)
  * For each interface:

    * Compute the interface **center**
    * Store list of **interacting residue indices**
    * Estimate a **binding energy** using a residue-pair energy lookup (`energy_table.py`)

### Output:

Each chain is assigned:

```python
{
    'COMs': [...],
    'radii': [...],
    'interfaces': [['B', 'C'], ...],
    'interface_coords': [[...], ...],
    'interface_residues': [[...], ...],
    'interface_energies': [[...], ...]
}
```

---

## **How interfaces are assigned types (repeated classification)**

### Module: `repeated_chain.py`

Chains are grouped into **repeated sets** using:

### Method A. **Header-based entity IDs** (from CIF header)

* Extract `_entity_poly.pdbx_strand_id` to map chains from same molecule entity
* Example: chains `A`, `B`, and `C` all belong to entity `1`, so are repeated

### Method B. **Structure-based fallback** (if header missing or .pdb file):

* Perform all-by-all **Cα alignment** using `Bio.PDB.Superimposer`
* Compute RMSD between chain pairs
* Assign chains to the same repeated_chain group if **RMSD < threshold (e.g., 1 Å)**

### Output:

* `chains_map`: `{ 'B': 'A', 'C': 'A', 'D': 'D' }` (map chain → representative)
* `chains_group`: `[ ['A', 'B', 'C'], ['D'] ]`

Used by `regularize.py` to merge redundant chains or build symmetry-aware models.

---

## **How to use the subpackage**

### Goal: From a `.cif` or `.pdb` file → NERDSS-ready model with interfaces, reactions, and molecule templates.

### Minimal Usage

```python
from model.pdb import coarse_grain, repeated_chain, regularize, templates, reaction

# 1. Step: Coarse-grain the structure to detect interfaces
raw = coarse_grain.coarse_grain_structure("example.cif")

# 2. Step: Find repeated chain groups (from CIF header or RMSD)
chains_map, chains_group = repeated_chain.identify_repeated_chains("example.cif", raw["structure"])

# 3. Step: Regularize chain and interface indices
model = regularize.regularize_model(raw, chains_map=chains_map)

# 4. Step: Convert to NERDSS Molecule / Interface templates
template_data = templates.build_templates(model)

# 5. Step: Generate Reaction templates
interface_list = template_data["interfaces"]  # flat list or nested per chain
reactions = reaction.build_binding_reactions(model, interface_list)

# Now: model + templates + reactions = ready for export
```

---

## Directory Structure Summary

```bash
model/
└── pdb/
    ├── coarse_grain.py      # Detect COMs, radii, and interfaces from structure
    ├── energy_table.py      # Heuristic residue-residue energy table
    ├── geometry.py          # Tools for alignment, clash checking, RMSD
    ├── repeated_chain.py           # Detect repeated chains via header or RMSD
    ├── molecule.py          # MoleculeTemplate and InterfaceTemplate classes
    ├── reaction.py          # ReactionTemplate logic between interfaces
    ├── regularize.py        # Reindex chains/interfaces, symmetry grouping
    ├── templates.py         # Build NERDSS-ready molecule/interface templates
```

---

# Coarse-Graining Protein Structures: Interface Detection Walkthrough

This document describes the algorithm implemented in `ionerdss.model.pdb.coarse_grain`. It detects chain-chain interfaces from a PDB structure using residue contact heuristics and computes coarse-grained interface representations suitable for downstream modeling.

## Overview

The coarse-graining process extracts the following from a Biopython `Structure` object:

* Center of mass (COM) of each chain
* Interacting residues between chain pairs (interfaces)
* Geometric location of interfaces (average of contacting residue CAs)
* Energetic estimate of interface strength based on residue pair potential

---

## Parameters

| Parameter         | Type                          | Description                                                |
| ----------------- | ----------------------------- | ---------------------------------------------------------- |
| `structure`       | `Bio.PDB.Structure.Structure` | Parsed PDB structure                                       |
| `distance_cutoff` | `float`                       | Max CA-CA distance (in nm) to count as contacting          |
| `residue_cutoff`  | `int`                         | Min number of contacting residues to consider an interface |

---

## Workflow Breakdown

### 1. Filter Valid Chains

```python
chains = sorted([chain for chain in structure.get_chains()
                 if any(is_aa(res) for res in chain.get_residues())],
                key=lambda c: c.id)
```

We keep only chains that contain amino acid residues. Sorting ensures consistent output.

---

### 2. Compute COM and Radius for Each Chain

```python
atoms = [atom.coord for res in chain for atom in res if is_aa(res)]
```

We collect coordinates of all atoms from amino acid residues to:

* Compute the center of mass (simple mean of all atom positions)
* Estimate chain spread (RMS distance from COM)

---

### 3. Precompute Bounding Boxes

```python
coords = np.array([atom.coord for res in chain for atom in res if is_aa(res)])
return coords.min(axis=0), coords.max(axis=0)
```

Bounding boxes are used to **exclude distant chain pairs** quickly.
If bounding boxes don’t overlap within a margin, skip the pair.

---

### 4. For Each Chain Pair (i, j)

We:

* Extract all residues that have alpha carbons
* Build a list of `(residue_id, type, CA_coord)`
* Call `compute_interface()` to:

  * Use KDTree to find CA-CA contacts within `distance_cutoff`
  * Track contacting residue IDs and types
  * Compute the mean CA coordinate as interface COM
  * Compute total interaction energy from a residue pair potential table

---

### 5. Store Detected Interfaces

```python
interfaces[i].append(chain_j.id)
interface_coords[i].append(com_i)
interface_residues[i].append(sorted(ids_i))
interface_energies[i].append(total_energy)
```

Interface metadata is saved in per-chain lists.
Both directions (i->j and j->i) are recorded symmetrically.

---

## Return Format

```python
{
  'chains': list of Bio.PDB.Chain objects,
  'COMs': list of Coords (per-chain center of mass),
  'radii': list of float (per-chain RMS radius),
  'interfaces': list of lists of chain IDs (who each chain interacts with),
  'interface_coords': list of Coords (per-chain list of interface locations),
  'interface_residues': list of lists of residue IDs,
  'interface_energies': list of floats
}
```

---

## Notes

* All distances are assumed in nanometers (consistent with simulation input units)
* Residue-residue energies are symmetric (i.e., E(A,B) = E(B,A))
* Function relies on `compute_interface()` and `get_default_energy_table()` to handle geometric and energetic logic

---

## Dependencies

* `Bio.PDB` for parsing
* `Coords` class (user-defined or from utility module)
* `KDTree` from `scipy.spatial`
* `numpy` for math operations


---

# Repeated Chain Detection Walkthrough


Here is a thorough and well-structured `.md` documentation file explaining the purpose, logic, and usage of the repeated_chain detection code in `.pdb.repeated_chain_detection`:

---

# Repeated Chain Detection in PDB Structures

`ionerdss.model.pdb.repeated_chain_detection` module provides utilities to **identify repeated chains** in a protein structure parsed with Biopython. It supports two strategies:

1. **Parsing metadata from mmCIF headers** (preferred and fast).
2. **Fallback to structural alignment of chain coordinates** (slower but robust).

## Purpose

Many PDB structures contain **repeated, symmetric subunits** such as dimers, trimers, or larger complexes. These subunits often consist of chains that are **structurally and/or sequence-wise repeated**.

To simplify downstream analysis (e.g., coarse-graining, interface detection, symmetry analysis), it’s useful to **group equivalent chains** and select **canonical representatives** for each group.

This module provides functions to:

* Identify repeated chain groups (e.g., `['A', 'B', 'C']`)
* Map each chain to its group representative (e.g., `{'A': 'A', 'B': 'A', 'C': 'A'}`)

---

## Functions Overview

### `identify_repeated_chains(pdb_id, structure)`

Primary function that attempts to identify repeated chains using two approaches:

1. **mmCIF Header Parsing** (fast, preferred):

   * Extracts chain grouping information from the `_entity_poly.pdbx_strand_id` and `_entity_poly.entity_id` fields in the mmCIF metadata.
   * Builds a `chains_map` and `chains_group` from this metadata.

2. **Fallback: Structural Alignment** (slow but general):

   * If CIF fields are missing or malformed, falls back to `_find_repeated_chains_by_alignment`.

**Parameters:**

* `pdb_id (str)` – PDB ID or path for error messages.
* `structure (Bio.PDB.Structure.Structure)` – Biopython structure object.

**Returns:**

* `chains_map (dict)` – Map each chain ID to its canonical representative.
* `chains_group (list[list[str]])` – List of repeated chain groups.

---

### `_find_repeated_chains_by_alignment(chains, rmsd_threshold=1.0)`

Fallback method that performs **pairwise structural alignment** using **C-alpha atoms only**. Chains are grouped based on RMSD similarity.

**Algorithm:**

* Loops through each chain.
* Extracts its Cα coordinates.
* Aligns it to each other unvisited chain.
* If RMSD < `rmsd_threshold`, considers the chains repeated.

**Parameters:**

* `chains (list)` – List of Biopython `Chain` objects.
* `rmsd_threshold (float)` – Max RMSD for chains to be grouped (default: 1.0 Å).

**Returns:**

* `chains_map (dict)` – Each chain mapped to a canonical ID.
* `chains_group (list[list[str]])` – Groups of equivalent chains.

---

### `_assign_original_chain_ids(chains_group, chains_map, structure)`

Sorts and reorders chain groups to follow the **original appearance order** in the PDB file.

**Purpose:** Ensure consistency and reproducibility by preserving order.

**Parameters:**

* `chains_group (list of lists)` – Groups of repeated chains.
* `chains_map (dict)` – Mapping from chain ID to representative.
* `structure (Structure)` – Biopython structure to extract chain order.

**Returns:**

* Updated `chains_map` and `chains_group`.

---

### `_validate_chain_mapping(chains_map, structure)`

Basic sanity check to ensure **every chain** in the structure has an assigned group.

**Raises:**

* `ValueError` if any chain is unmapped.

---

## Example Usage

```python
from Bio.PDB import MMCIFParser
from ionerdss.model.pdb.repeated_chain_detection import identify_repeated_chains

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("1ABC", "1abc.cif")

chains_map, chains_group = identify_repeated_chains("1abc", structure)

print("Chain Mapping:", chains_map)
print("Chain Groups:", chains_group)
```

**Example Output:**

```python
Chain Mapping: {'A': 'A', 'B': 'A', 'C': 'A', 'D': 'D'}
Chain Groups: [['A', 'B', 'C'], ['D']]
```

This means chains A, B, and C are repeated, and D is a singleton group.

---

## Notes and Considerations

* **CIF Header Parsing is faster and more reliable**, but not all mmCIF files include the necessary metadata.
* **Structural alignment** is slower, especially for large assemblies, but works on any structure containing valid Cα atoms.
* This module **does not depend on sequence alignment**, making it robust to chain ID or sequence name differences.
* The `_assign_original_chain_ids` utility ensures that analysis follows the **original chain ordering** in the PDB file.

