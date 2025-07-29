# 1. `ionerdss.model.pdb` Module Overview

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

# 2. Coarse-Graining Protein Structures: Interface Detection Walkthrough

This document describes the algorithm implemented in `ionerdss.model.pdb.coarse_grain`. It detects chain-chain interfaces from a PDB structure using residue contact heuristics and computes coarse-grained interface representations suitable for downstream modeling.

## Overview

The coarse-graining process extracts the following from a Biopython `Structure` object:

* Center of mass (COM) of each chain
* Interacting residues between chain pairs (interfaces)
* Geometric location of interfaces (average of contacting residue CAs)
* Energetic estimate of interface strength based on residue pair potential

Note that we are using **CA-CA** distance for now! (More options?)

## Parameters

| Parameter         | Type                          | Description                                                |
| ----------------- | ----------------------------- | ---------------------------------------------------------- |
| `distance_cutoff` | `float`                       | Max CA-CA distance (in nm) to count as contacting          |
| `residue_cutoff`  | `int`                         | Min number of contacting residues to consider an interface |

## Workflow Breakdown

### 1. Filter Valid Chains

```python
chains = sorted([chain for chain in structure.get_chains()
                 if any(is_aa(res) for res in chain.get_residues())],
                key=lambda c: c.id)
```

We keep only chains that contain amino acid residues. Sorting ensures consistent output.

### 2. Compute COM and Radius for Each Chain

```python
atoms = [atom.coord for res in chain for atom in res if is_aa(res)]
```

We collect coordinates of all atoms from amino acid residues to:

* Compute the center of mass (simple mean of all atom positions)
* Estimate chain spread (RMS distance from COM)

### 3. Precompute Bounding Boxes

```python
coords = np.array([atom.coord for res in chain for atom in res if is_aa(res)])
return coords.min(axis=0), coords.max(axis=0)
```

Bounding boxes are used to **exclude distant chain pairs** quickly.
If bounding boxes don’t overlap within a margin, skip the pair.

### 4. For Each Chain Pair (i, j)

We:

* Extract all residues that have alpha carbons
* Build a list of `(residue_id, type, CA_coord)`
* Call `compute_interface()` to:

  * Use KDTree to find CA-CA contacts within `distance_cutoff`
  * Track contacting residue IDs and types
  * Compute the mean CA coordinate as interface COM
  * Compute total interaction energy from a residue pair potential table

### 5. Store Detected Interfaces

```python
interfaces[i].append(chain_j.id)
interface_coords[i].append(com_i)
interface_residues[i].append(sorted(ids_i))
interface_energies[i].append(total_energy)
```

Interface metadata is saved in per-chain lists.
Both directions (i->j and j->i) are recorded symmetrically.

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

## Notes

* All distances are assumed in nanometers (consistent with simulation input units)
* Residue-residue energies are symmetric (i.e., E(A,B) = E(B,A))
* Function relies on `compute_interface()` and `get_default_energy_table()` to handle geometric and energetic logic

## Dependencies

* `Bio.PDB` for parsing
* `Coords` class (user-defined or from utility module)
* `KDTree` from `scipy.spatial`
* `numpy` for math operations


# 3. Repeated Chain Detection Walkthrough

This module provides utility functions to identify **repeated protein chains** in a biological structure (e.g., derived from PDB/mmCIF files). These chains may be identical or nearly identical due to symmetry, assembly, or template-based modeling. Identifying such equivalence is essential for **coarse-graining** and **downstream simulation** steps in tools like **NERDSS**.

## Purpose

The primary purpose of this module is to:

* Detect **structurally or sequence-wise equivalent chains** within a structure.
* Assign a **canonical representative** to each group of equivalent chains.
* Enable **template-based molecular modeling**, where symmetric copies are replaced by references to a single entity.

This is useful for:

* Molecular symmetry analysis
* Template construction for simulations
* Reducing model redundancy
* Enabling consistent labeling across assemblies

## Core Logic

### Main Entry Point

```python
identify_repeated_chains(pdb_id, structure, mode='default', ...)
```

This is the **main function** to detect repeated chains. It supports three modes of operation:

| Mode          | Description                                                                 |
| ------------- | --------------------------------------------------------------------------- |
| `'default'`   | Uses mmCIF header (if available), otherwise falls back to sequence identity |
| `'structure'` | Forces structural superimposition-based matching                            |
| `'sequence'`  | Forces sequence identity-based matching                                     |

The function returns:

* `chains_map`: a dictionary mapping each chain ID to a canonical representative
* `chains_group`: a list of lists representing chain equivalence groups

## Methods and Algorithms

### 1. Header-Based Detection (Default)

```python
if '_entity_poly.pdbx_strand_id' in mmcif_dict:
    # Use mmCIF annotations to group chains
```

* Extracts chain groups from `_entity_poly.entity_id` and `_entity_poly.pdbx_strand_id`
* Fast and direct, but depends on correct mmCIF annotations

### 2. Structural Superimposing (Fallback or Forced)

```python
_find_repeated_chains_by_structure_superimposing(chains, rmsd_threshold=1.0)
```

* Aligns each pair of chains using **C-alpha RMSD**
* Chains are considered equivalent if RMSD < threshold (default: 1.0 Å)
* Uses Biopython’s `Superimposer`

**Advantages:**

* Works without metadata
* Captures subtle differences in structure

**Limitations:**

* Assumes similar length and backbone topology
* Fails if chains differ significantly in number of residues

### 3. Sequence-Based Matching (Fallback or Forced)

```python
_find_repeated_chains_by_sequence(chains, seq_threshold=0.95)
```

* Uses `Bio.Align.PairwiseAligner` to align peptide sequences
* Groups chains if identity ≥ threshold (default: 95%)

You can also provide a **custom aligner** to modify scoring schemes (e.g., for gap penalties).

**Advantages:**

* Robust to small structural variations
* Works even when structure is missing or imprecise

**Limitations:**

* Requires valid sequences (must have N, CA, and C atoms to be parsed)
* Sensitive to missing backbone data

## Parameters

| Parameter        | Type                          | Description                                                    |
| ---------------- | ----------------------------- | -------------------------------------------------------------- |
| `mode`           | `str`                         | `'default'`, `'structure'`, or `'sequence'`                    |
| `rmsd_threshold` | `float`                       | Max C-alpha RMSD to consider chains as structurally equivalent |
| `seq_threshold`  | `float`                       | Min sequence identity to consider chains equivalent            |
| `custom_aligner` | `Bio.Align.PairwiseAligner`   | Optional user-defined aligner instance                         |

## Example Usage

```python
from Bio.PDB import MMCIFParser
from ionerdss.model.pdb.repeated_chain_detection import identify_repeated_chains

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("1ABC", "1abc.cif")

chains_map, chains_group = identify_repeated_chains("1ABC", structure, mode='default')

print(chains_map)
# {'A': 'A', 'B': 'A', 'C': 'C'}

print(chains_group)
# [['A', 'B'], ['C']]
```


## Internal Outputs

| Variable       | Type   | Description                           |
| -------------- | ------ | ------------------------------------- |
| `chains_map`   | `dict` | Each chain mapped to its canonical ID |
| `chains_group` | `list` | Groups of chain IDs                   |

## Error Handling

* If `mode='default'` and CIF header is missing, a warning is issued:

  ```python
  warnings.warn(f"{pdb_id}: Missing mmCIF entity_poly info. Falling back...")
  ```

* If an invalid mode is specified:

  ```python
  raise ValueError(...)
  ```

## Caveats

* `sequence` mode depends on proper peptide atom construction (`N`, `CA`, `C`). Synthetic structures must mimic real peptide geometry. i.e. currently only work for proteins or polypeptide-like structures. 
* `structure` mode assumes all chains are alignable — it skips chains with missing C-alpha atoms.
* This tool doesn't consider reflection/inversion symmetry — only RMSD/sequence identity.

Here's a detailed module-level documentation for the `regularize_repeats.py` file, explaining its **purpose, logic, and algorithmic pipeline** within the broader NERDSS framework:

---

# 4. `regularize_repeats.py` — molecule template regularization for repeated chains

## Purpose

This module processes coarse-grained molecular models that contain **repeated protein chains** (i.e. symmetric or repeated units in the same PDB structure). Its primary goal is to generate a **uniform, symmetry-aware representation** of:

* **Molecule templates**: shared geometry, diffusion constants, and structure.
* **Interface templates**: representative coordinates, signature geometry, and orientations.
* **Binding interactions**: consistent partner naming and directional symmetry (e.g., `A_2` means interface on molecule that binds another `A`).

---

## Problem Context

After **coarse-graining** a molecular structure:

* Each protein chain is abstracted into a center-of-mass (COM) and interaction sites.
* Chains may repeat (e.g., A1, A2, A3) but have different interface order or coordinate noise.
* Chain-chain interactions (e.g., dimers, trimers) may be structurally symmetric or asymmetric.

What’s needed:

* A **single consistent template** per molecular unit (e.g., A) across all its repeats.
* A **shared interface template** for each type of interaction (e.g., A+A vs A+B).
* **Canonical geometry** that averages or registers noisy replicas into one reference model.
* **Symmetry-aware interface registration** for use in NERDSS reaction network generation.

---

## Core Logic

The pipeline performs the following key steps:

### 1. **Group and Iterate Over Repeated Chains**

Each chain group (e.g., `[A, B, C]`) is assumed to consist of repeated or repeated chains.

```python
for group in chains_groups:
    process_chain_group(...)
```

---

### 2. **Assign Templates and Molecules**

* Assign a **MoleculeTemplate** object to each unique chain type.
* For each chain instance, create a **CoarseGrainedMolecule** with COM, radius, and diffusion constants.

```python
get_or_create_molecule_template(...)
get_or_create_molecule(...)
```

---

### 3. **Identify Interfaces Between Chains**

For each interaction:

* Extract interface coordinates (COM + site)
* Compute a **signature**: `{dA, dB, θA, θB}` — capturing geometry between two chains
* Invert the signature to check the reverse direction

Then:

* If this signature has been seen before, reuse the **existing interface templates**
* If not, generate **new interface templates** and register them with molecule templates

```python
build_signature(COM_A, I_A, COM_B, I_B)
invert_signature(...)
signature_hash(...)
```

---

### 4. **Determine Homodimer vs Heterodimer**

* If both chains map to the **same template**, and the geometry is symmetric:
  → It’s a **homodimer** (`A_2`)
* Otherwise, it’s a **heterodimer** (`B_1` on A; `A_1` on B)

```python
if chains_map[mol_a] == chains_map[mol_b] and geometry_is_symmetric(...):
    is_repeated = True
```

---

### 5. **Align All Chains to First Copy**

To ensure consistent geometry across all repeats:

* Rigidly align each repeated chain to the first in the group using `rigid_transform_chains`
* Apply the transform to COM and interface coordinates
* Store relative interface offsets for consistent modeling

```python
R, t = rigid_transform_chains(chain2, chain1)
apply_rigid_transform(...)
```

---

### 6. **Detect Required-Free Interface Conflicts**

To avoid steric clashes from simultaneous use of incompatible interfaces:

* Loop over all interfaces across repeated chains
* Superimpose potential partner positions
* If overlap is detected, mark those templates as **mutually exclusive** (required\_free)

```python
if check_clashes_between_two_sets(...):
    intf1.required_free_list.append(...)
```

---

### 7. **Build MoleculeTypes**

Generate `MoleculeType` objects with:

* Template name
* Interface coordinates (global, not relative)
* Diffusion constants

Used later for:

* Network generation
* Simulation setup

```python
MoleculeType(name, interfaces, ...)
```

---

## Interface Naming Convention

Interface templates are named by their **binding partner**:

* `E_2` on molecule `A` means “this site on A binds to E”
* `A_2` on `A` means “this site on A binds to another A” (i.e. dimer)

This convention ensures:

* Partner identity is clear
* Homodimeric interactions are recognizable
* Templates are grouped logically

---

## Customization & Parameters

* `dist_threshold_intra`: Tolerance for geometric symmetry (used in homo-dimer detection)
* `dist_threshold_inter`: Tolerance for different molecule types
* `angle_threshold`: Angular deviation for signature matching

---

## Returns

The final result is a structured model:

```python
{
    "molecule_templates": [...],
    "molecules": [...],
    "interface_templates": [...],
    "interfaces": [...],
    "binding_pairs": [...],
    "molecule_types": [...],
}
```

Used downstream for:

* Network export
* `.mol` and `.parm.inp` files
* PyMOL rendering and visual validation

---

## Summary

This module performs **symmetry-aware, structure-consistent regularization** of repeated molecular chains by:

* Recognizing **structural repetitions** across chains
* Computing **geometric interface signatures**
* **Aligning**, **merging**, and **registering** templates consistently
* Producing **uniform molecule/interface templates** for simulation and analysis
