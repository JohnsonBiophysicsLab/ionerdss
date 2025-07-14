Here is a **comprehensive walkthrough** of how the `model.pdb` subpackage works after your modular refactor. This includes:

1. **How interfaces are detected from a PDB (or mmCIF) file**
2. **How binding sites are grouped into types (homologs)**
3. **How a user uses the subpackage to generate a NERDSS model**

---

## **How binding sites / interfaces are detected from a PDB or mmCIF file**

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

## **How interfaces are assigned types (homologous classification)**

### Module: `homolog.py`

Chains are grouped into **homologous sets** using:

### Method A. **Header-based entity IDs** (from CIF header)

* Extract `_entity_poly.pdbx_strand_id` to map chains from same molecule entity
* Example: chains `A`, `B`, and `C` all belong to entity `1`, so are homologous

### Method B. **Structure-based fallback** (if header missing or .pdb file):

* Perform all-by-all **Cα alignment** using `Bio.PDB.Superimposer`
* Compute RMSD between chain pairs
* Assign chains to the same homolog group if **RMSD < threshold (e.g., 1 Å)**

### Output:

* `chains_map`: `{ 'B': 'A', 'C': 'A', 'D': 'D' }` (map chain → representative)
* `chains_group`: `[ ['A', 'B', 'C'], ['D'] ]`

Used by `regularize.py` to merge redundant chains or build symmetry-aware models.

---

## **How to use the subpackage**

### Goal: From a `.cif` or `.pdb` file → NERDSS-ready model with interfaces, reactions, and molecule templates.

### Minimal Usage

```python
from model.pdb import coarse_grain, homolog, regularize, templates, reaction

# 1. Step: Coarse-grain the structure to detect interfaces
raw = coarse_grain.coarse_grain_structure("example.cif")

# 2. Step: Find homologous chain groups (from CIF header or RMSD)
chains_map, chains_group = homolog.identify_homologous_chains("example.cif", raw["structure"])

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
    ├── homolog.py           # Detect homologous chains via header or RMSD
    ├── molecule.py          # MoleculeTemplate and InterfaceTemplate classes
    ├── reaction.py          # ReactionTemplate logic between interfaces
    ├── regularize.py        # Reindex chains/interfaces, symmetry grouping
    ├── templates.py         # Build NERDSS-ready molecule/interface templates
```
