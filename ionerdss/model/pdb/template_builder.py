"""
ionerdss.model.pdb.template_builder

Molecular and interface template generation with proper signature-based deduplication.

This module builds molecular templates and interface templates from grouped
chains, handles geometric signature calculation, template deduplication based
on signatures, and regularization across symmetry mates.


## Key Concepts

### Template-Based Modeling Rule

**Molecular Templates**: Reusable molecular components that capture the essential
geometric and interaction properties of protein chains, allowing multiple
instances to be created from a single template.

**Interface Templates**: Geometric and energetic descriptions of protein-protein
interaction sites that can be reused across multiple molecular instances.

**Signature-Based Deduplication**: Multiple interface types can exist between the
same molecule types if they have different geometric signatures, enabling complex
multi-site interactions.

### Template Hierarchy

```
System Level:
├── Molecule Types (Templates)
│   ├── Radius, diffusion constants
│   ├── Interface binding sites
│   └── Geometric properties
└── Interface Types (Templates)
    ├── Geometric signature
    ├── Binding energy
    ├── Local coordinates
    └── Partner relationships
```

### Coordinate System Management

**Input**: Coordinates in Angstroms (from structural data)
**Processing**: Automatic conversion to nanometers for NERDSS compatibility
**Local Coordinates**: Interface positions relative to molecule center of mass
**Absolute Coordinates**: Interface positions in global coordinate system

## Geometric Signatures

### Signature Components

The `GeometricSignature` captures the essential geometric relationship of an interface:

```python
@dataclass
class GeometricSignature:
    d_i: float      # COM-to-interface distance for chain i
    d_j: float      # COM-to-interface distance for chain j  
    theta_i: float  # Angle between COM-to-interface and COM-to-COM vectors (chain i)
    theta_j: float  # Angle between COM-to-interface and COM-to-COM vectors (chain j)
```

### Signature Calculation

**Distance Components**:
```python
d_i = ||interface_coord_i - COM_i||
d_j = ||interface_coord_j - COM_j||
```

**Angular Components**:
```python
offset_i = interface_coord_i - COM_i
com_vector_ij = COM_j - COM_i

theta_i = arccos(dot(offset_i, com_vector_ij) / (||offset_i|| * ||com_vector_ij||))
```

**Geometric Interpretation**:
- **d_i, d_j**: How far the interface is from each molecule's center
- **theta_i, theta_j**: The orientation of the interface relative to the intermolecular vector

### Signature Applications

**Similarity Detection**:
```python
signature1.is_similar_to(signature2, 
                        distance_threshold=5.0,  # Angstroms
                        angle_threshold=0.5)     # Radians (~30°)
```

**Homotypic Detection**:
```python
# Detect symmetric homodimer interfaces
is_homotypic = signature.is_homotypic(
    distance_threshold=1.0,  # Angstroms
    angle_threshold=0.2      # Radians (~11°)
)
```

## Template Generation Process

### 1. Molecular Template Creation

```python
def _build_molecule_template(self, group: ChainGroup) -> None:
```

**Process**:
1. **Representative Selection**: Use group representative as template basis
2. **Name Generation**: Create unique template names with conflict resolution
3. **Radius Calculation**: Convert from Angstroms to nanometers
4. **Diffusion Constants**: Calculate from molecular radius using Stokes-Einstein relation
5. **Metadata Storage**: Store grouping information and original chain names

**Template Properties**:
```python
molecule_template = MoleculeType(
    name="ProteinA",                    # Unique template name
    radius_nm=1.5,                      # Radius in nanometers
    diffusion_constants_calculated=True  # Auto-calculated from radius
)

# Metadata for traceability
molecule_template.signature = {
    'group_representative': 'A',
    'group_members': ['A', 'C', 'E'],
    'grouping_method': 'sequence_similarity',
    'original_chain_names': ['A', 'C', 'E']
}
```

### 2. Interface Template Creation

```python
def _build_all_interface_templates(self) -> None:
```

**Workflow**:
1. **Interface Processing**: Iterate through all detected interfaces
2. **Group Resolution**: Map chains to their representative groups
3. **Template Mapping**: Find corresponding molecule templates
4. **Signature Calculation**: Compute geometric signature for each interface
5. **Deduplication**: Check for existing similar interface types
6. **Template Creation**: Create new templates or reuse existing ones

### 3. Template Naming Strategy

**Molecule Templates**:
```python
# Preferred: Use representative chain name
"A" → "A"

# Conflict resolution: Add descriptive suffix
"A" (taken) → "A_group"

# Final fallback: Numeric suffix
"A_group" (taken) → "A_1"
```

**Interface Templates**:
```python
# Homodimer homotypic
"A" + "A" + index → "A_A_1"

# Heterotypic (bidirectional)
"A" + "B" + index → "A_B_1", "B_A_1"

# Multiple interfaces between same types
"A" + "B" + different_signature → "A_B_2", "B_A_2"
```

## Deduplication Strategy

### Signature-Based Matching

**Matching Criteria**:
```python
def _find_matching_interface_type(self, template_i: str, template_j: str,
                                  signature: GeometricSignature) -> Optional[str]:
```

**Process**:
1. **Template Pair Check**: Find interfaces between same molecule types
2. **Signature Comparison**: Use relaxed thresholds for matching
3. **Bidirectional Matching**: Handle both A→B and B→A orientations

**Matching Thresholds**:
```python
distance_threshold = 5.0   # 5 Angstroms tolerance
angle_threshold = 0.5      # ~30 degrees tolerance
```

### Multiple Interface Types

**Same Molecule Pair, Different Signatures**:
```python
# Interface 1: Close to molecule centers
signature_1 = GeometricSignature(d_i=3.0, d_j=3.0, theta_i=0.2, theta_j=0.2)
→ Creates: A_B_1, B_A_1

# Interface 2: Far from molecule centers  
signature_2 = GeometricSignature(d_i=8.0, d_j=8.0, theta_i=1.4, theta_j=1.4)
→ Creates: A_B_2, B_A_2
```

## Interface Types

### Homodimer Homotypic Interfaces

**Characteristics**:
- Same molecule type on both sides
- Symmetric geometric signature
- Single shared interface template

**Creation**:
```python
# Symmetric signature detection
is_homotypic = (abs(d_i - d_j) < threshold and 
                abs(theta_i - theta_j) < threshold)

if is_homotypic:
    # Create single shared template
    interface_template = InterfaceType(
        this_mol_type_name="A",
        partner_mol_type_name="A",
        interface_index=1,
        # ... coordinates and energy
    )
```

## Usage Examples

### Basic Template Building

```python
from ionerdss.model.pdb.template_builder import TemplateBuilder

# Build templates from processed components
builder = TemplateBuilder(
    parser=parser,
    coarse_grainer=coarse_grainer,
    chain_grouper=chain_grouper,
    hyperparams=hyperparams,
    workspace_manager=workspace_manager
)

# Get generated templates
molecule_templates = builder.get_molecule_templates()
interface_templates = builder.get_interface_templates()

print(f"Created {len(molecule_templates)} molecule templates:")
for name, template in molecule_templates.items():
    print(f"  {name}: radius={template.radius_nm:.2f} nm")

print(f"Created {len(interface_templates)} interface templates:")
for name, template in interface_templates.items():
    print(f"  {name}: {template.this_mol_type_name} ↔ {template.partner_mol_type_name}")
```

### Template Analysis

```python
# Get comprehensive summary
summary = builder.get_summary()

print("Template Building Summary:")
print(f"  Molecule Templates: {summary['num_molecule_templates']}")
print(f"  Interface Templates: {summary['num_interface_templates']}")

# Analyze interface type distribution
interface_counts = summary['interface_type_counts_by_molecule_pair']
for mol_pair, count in interface_counts.items():
    print(f"  {mol_pair[0]} ↔ {mol_pair[1]}: {count} interface types")

# Chain to template mapping
chain_mapping = summary['chain_name_mapping']
print("Chain → Template Mapping:")
for chain, template in chain_mapping.items():
    print(f"  Chain {chain} → Template {template}")
```


"""

from typing import Dict, List, Tuple, Set, Optional
from dataclasses import dataclass, field

import numpy as np
from Bio.PDB.Superimposer import Superimposer

from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.units import Units
from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser
from .coarse_graining import CoarseGrainer, InterfaceString, CoarseGrainedChain
from .chain_grouping import ChainGrouper, ChainGroup
from .file_manager import WorkspaceManager


import numpy as np
from typing import Dict, List, Tuple, Any

import numpy as np
from typing import Dict, List, Tuple, Any

def _kabsch_transform(P: np.ndarray, Q: np.ndarray):
    """
    Classic Kabsch: find R, t that minimize || R@P + t - Q ||_F
    P, Q: (N,3) point sets with N>=3 and in 1–1 correspondence.
    Returns:
        R (3x3), t (3,)
    """
    assert P.shape == Q.shape and P.shape[1] == 3, "P and Q must be (N,3)"

    # centroids
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)

    # centered
    P0 = P - Pc
    Q0 = Q - Qc

    # covariance
    H = P0.T @ Q0

    # SVD
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # proper rotation (avoid reflection)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    t = Qc - R @ Pc
    return R, t


def _collect_points_for_fit_with_normals(ref, mem, normal_lever: float = 5.0):
    """
    Build correspondence lists (P from ref, Q from mem) using:
      - COM (if present)
      - named sites (aa1.., or obj.sites dict)
      - interfaces matched by a stable 'key'
      - OPTIONAL: if both sides have a normal for the same interface key,
                  add an extra correspondence at coord + lever * normal
                  to give Kabsch roll sensitivity.

    Returns:
        P, Q  (each (N,3) np.ndarray; N may be >= 0)
    """
    P, Q = [], []

    # --- COM
    if hasattr(ref, "com") and hasattr(mem, "com"):
        P.append(np.asarray(ref.com, float))
        Q.append(np.asarray(mem.com, float))

    # --- named sites
    def _get_sites_dict_local(obj):
        d = {}
        if hasattr(obj, "sites") and isinstance(obj.sites, dict):
            for k, v in obj.sites.items():
                d[str(k)] = np.asarray(v, float)
        else:
            for k in ("aa1", "aa2", "aa3", "aa4"):
                if hasattr(obj, k):
                    d[k] = np.asarray(getattr(obj, k), float)
        return d

    s_ref = _get_sites_dict_local(ref)
    s_mem = _get_sites_dict_local(mem)
    for name in sorted(set(s_ref) & set(s_mem)):
        P.append(s_ref[name])
        Q.append(s_mem[name])

    # --- interfaces by stable key
    def _rich(obj):
        out = []
        if hasattr(obj, "interfaces") and obj.interfaces is not None:
            for i, raw in enumerate(obj.interfaces):
                # normalize dict-like access
                if isinstance(raw, dict):
                    coord = raw.get("coord")
                    normal = raw.get("normal", None)
                    key = raw.get("key", None)
                else:
                    coord = getattr(raw, "coord", None)
                    normal = getattr(raw, "normal", None)
                    key = getattr(raw, "key", None)
                if coord is None:
                    continue
                out.append({
                    "coord": np.asarray(coord, float),
                    "normal": None if normal is None else np.asarray(normal, float),
                    "key": str(key) if key is not None else f"index:{i}",
                })
        return out

    r_if = _rich(ref)
    m_if = _rich(mem)
    mem_by_key = {e["key"]: e for e in m_if}

    for e in r_if:
        k = e["key"]
        if k not in mem_by_key:
            continue
        # base coord pair
        P.append(e["coord"])
        Q.append(mem_by_key[k]["coord"])

        # add normal-lever pair if both have normals
        nr = e.get("normal")
        nm = mem_by_key[k].get("normal")
        if nr is not None and nm is not None:
            # use unit normals; if zero, skip
            nr = np.asarray(nr, float)
            nm = np.asarray(nm, float)
            nrm = np.linalg.norm(nr)
            nmm = np.linalg.norm(nm)
            if nrm > 1e-9 and nmm > 1e-9 and normal_lever > 0:
                Pr = e["coord"] + normal_lever * (nr / nrm)
                Qm = mem_by_key[k]["coord"] + normal_lever * (nm / nmm)
                P.append(Pr)
                Q.append(Qm)

    if len(P) == 0:
        return np.zeros((0, 3)), np.zeros((0, 3))

    return np.vstack(P), np.vstack(Q)


def _compute_fit_rmsd(P: np.ndarray, Q: np.ndarray, R: np.ndarray, t: np.ndarray) -> float:
    """RMSD between R@P + t and Q (Å)."""
    if P.size == 0:
        return 0.0
    Pfit = (R @ P.T).T + t
    diff = Pfit - Q
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _solve_best_com_shift_for_member(chain_obj, mean_dirs_by_key, reg_lambda: float = 1e-6) -> Optional[np.ndarray]:
    """
    Given current chain geometry (absolute coords) and group mean unit directions
    for each interface key, find a COM shift 'c' that minimizes angular discrepancy:
        min_c sum_k || (I - a_k a_k^T) * (v_k - c) ||^2
    where v_k = coord_k - COM (member), a_k = mean direction for key k.

    Returns the 3-vector shift (Å) or None if under-constrained.
    """
    if not hasattr(chain_obj, "com"):
        return None

    com = np.asarray(chain_obj.com, float)
    A = np.zeros((3, 3), float)
    b = np.zeros(3, float)
    have = 0

    for e in _get_interfaces_rich(chain_obj):
        key = e["key"]
        if key not in mean_dirs_by_key:
            continue
        a = np.asarray(mean_dirs_by_key[key], float)
        v = np.asarray(e["coord"], float) - com  # current member ray
        P = np.eye(3) - np.outer(a, a)           # projector to plane ⟂ a
        A += P
        b += P @ v
        have += 1

    if have == 0:
        return None

    # regularize in case the projector sum is ill-conditioned (e.g. all keys collinear)
    A += reg_lambda * np.eye(3)
    try:
        c = np.linalg.solve(A, b)  # Å
    except np.linalg.LinAlgError:
        return None
    return c

# --- Site access ---
def _get_sites_dict(obj: Any) -> Dict[str, np.ndarray]:
    sites: Dict[str, np.ndarray] = {}
    if hasattr(obj, "sites") and isinstance(obj.sites, dict):
        for k, v in obj.sites.items():
            sites[str(k)] = np.asarray(v, float)
    else:
        for k in ("aa1","aa2","aa3","aa4"):
            if hasattr(obj, k):
                sites[k] = np.asarray(getattr(obj, k), float)
    return sites

def _u(x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    n = float(np.linalg.norm(x))
    return x / n if n > eps else np.array([1.0, 0.0, 0.0], float)

# --- Interface key extraction ---
def _iface_key(raw: Any, idx: int) -> str:
    # Prefer semantic keys that are stable across instances
    for k in ("name","id","type","partner","binds","target","label","key"):
        if isinstance(raw, dict) and k in raw and raw[k]:
            return str(raw[k])
        if hasattr(raw, k) and getattr(raw, k) is not None:
            return str(getattr(raw, k))
    return f"index:{idx}"  # robust fallback

def _get_interfaces_rich(obj: Any) -> List[dict]:
    """Each entry: {'coord':(3,), 'normal':(3,)?,'key':str,'raw':raw,'index':i}"""
    out: List[dict] = []
    if hasattr(obj, "interfaces") and obj.interfaces is not None:
        for i, raw in enumerate(obj.interfaces):
            coord = raw.get("coord") if isinstance(raw, dict) else getattr(raw, "coord", None)
            normal = raw.get("normal") if isinstance(raw, dict) else getattr(raw, "normal", None)
            if coord is None:
                continue
            out.append({
                "coord": np.asarray(coord, float),
                "normal": None if normal is None else np.asarray(normal, float),
                "key": _iface_key(raw, i),
                "raw": raw,
                "index": i,
            })
    return out

# --- Fetch interface info to chain group ---

def _materialize_interfaces_on_chains_from_global(coarse_grainer, chains, group_members):
    for cid in group_members:
        chains[cid].interfaces = []

    per_pair_counters = {}
    cg_chains = coarse_grainer.get_coarse_grained_chains()

    for iface in coarse_grainer.get_interfaces():
        if iface.chain_i not in group_members and iface.chain_j not in group_members:
            continue

        base_key = str(getattr(iface, "interface_type", "")) or None
        if not base_key:
            pair = tuple(sorted((iface.chain_i, iface.chain_j)))
            per_pair_counters[pair] = per_pair_counters.get(pair, 0) + 1
            base_key = f"{pair[0]}-{pair[1]}-{per_pair_counters[pair]}"

        # σ is “toward partner COM”
        com_i = np.asarray(cg_chains[iface.chain_i].com, float)
        com_j = np.asarray(cg_chains[iface.chain_j].com, float)

        if iface.chain_i in group_members:
            chains[iface.chain_i].interfaces.append({
                "coord": np.asarray(iface.coord_i, float),
                "key": f"{base_key}:i",
                "sigma": com_j - com_i
            })
        if iface.chain_j in group_members:
            chains[iface.chain_j].interfaces.append({
                "coord": np.asarray(iface.coord_j, float),
                "key": f"{base_key}:j",
                "sigma": com_i - com_j
            })

# --- Means & radial normalize (keyed by interface identity) ---

def _enforce_identical_local_geometry_after_com(
    mem,
    R_rep_to_mem: np.ndarray,
    com_new: np.ndarray,
    site_canon_rep: dict,
    if_vec_canon_rep: dict,
    *,
    only_existing: bool = False,
) -> None:
    """
    Overwrite local geometry. If only_existing=True, update only the interfaces
    (and sites) that already exist on `mem` — never create new slots.
    """
    # Sites
    has_site_dict = hasattr(mem, "sites") and isinstance(mem.sites, dict)
    mem_sites = _get_sites_dict(mem)

    for name, v_rep in site_canon_rep.items():
        if only_existing and name not in mem_sites:
            continue
        v_world = R_rep_to_mem @ np.asarray(v_rep, float)
        p_world = com_new + v_world
        if has_site_dict:
            mem.sites[name] = p_world
        else:
            setattr(mem, name, p_world)

    # Interfaces
    if_list = getattr(mem, "interfaces", []) or []
    key_to_idx = {}
    for i, raw in enumerate(if_list):
        k = raw.get("key") if isinstance(raw, dict) else getattr(raw, "key", None)
        if k is not None:
            key_to_idx[str(k)] = i

    for key, v_rep in if_vec_canon_rep.items():
        if only_existing and key not in key_to_idx:
            continue
        v_world = R_rep_to_mem @ np.asarray(v_rep, float)
        p_world = com_new + v_world

        if key in key_to_idx:
            i = key_to_idx[key]
            raw = if_list[i]
            if isinstance(raw, dict):
                raw["coord"] = p_world
            else:
                setattr(raw, "coord", p_world)
        else:
            # only_existing=False path can create new
            if_list.append({"key": key, "coord": p_world})

    mem.interfaces = if_list

# -----------------------------------------------------------------------------
# Math helpers
# -----------------------------------------------------------------------------

def _unit(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < eps:
        return np.array([1.0, 0.0, 0.0])
    return v / n


# -----------------------------------------------------------------------------
# Weighted / RANSAC / IRLS Kabsch
# -----------------------------------------------------------------------------

def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    assert P.shape == Q.shape and P.shape[1] == 3
    Pc, Qc = P.mean(0), Q.mean(0)
    P0, Q0 = P - Pc, Q - Qc
    U, S, Vt = np.linalg.svd(P0.T @ Q0)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Qc - R @ Pc
    return R, t


def _kabsch_weighted(P: np.ndarray, Q: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Weighted Kabsch using weights w>=0 (shape (N,))."""
    w = np.clip(np.asarray(w, float), 0.0, np.inf)
    if P.shape[0] != Q.shape[0] or P.shape[1] != 3 or Q.shape[1] != 3:
        raise ValueError("P,Q shape mismatch")
    if w.shape != (P.shape[0],):
        raise ValueError("w shape mismatch")
    W = w[:, None]
    wsum = float(np.sum(w))
    if wsum <= 1e-12:
        return np.eye(3), np.zeros(3)
    Pc = (P * W).sum(0) / wsum
    Qc = (Q * W).sum(0) / wsum
    P0 = P - Pc
    Q0 = Q - Qc
    C = (W * P0).T @ Q0
    U, S, Vt = np.linalg.svd(C)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Qc - R @ Pc
    return R, t


def _ransac_kabsch(P: np.ndarray, Q: np.ndarray, thresh: float = 4.0,
                   max_trials: int = 200, min_inliers: Optional[int] = None,
                   rng: Optional[np.random.Generator] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """RANSAC to seed a robust rigid transform. Returns (R,t,inlier_mask)."""
    n = P.shape[0]
    if n < 3:
        return np.eye(3), np.zeros(3), np.ones(n, dtype=bool)
    if rng is None:
        rng = np.random.default_rng()
    if min_inliers is None:
        min_inliers = max(3, int(0.6*n))

    best_inliers = np.zeros(n, dtype=bool)
    best_R, best_t = np.eye(3), np.zeros(3)

    for _ in range(max_trials):
        idx = rng.choice(n, size=3, replace=False)
        R, t = _kabsch(P[idx], Q[idx])
        residuals = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
        inliers = residuals < thresh
        if inliers.sum() > best_inliers.sum():
            best_inliers = inliers
            best_R, best_t = R, t
            if best_inliers.sum() >= min_inliers:
                break

    # refine on inliers
    if best_inliers.sum() >= 3:
        R, t = _kabsch(P[best_inliers], Q[best_inliers])
        residuals = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
        best_inliers = residuals < thresh
        if best_inliers.sum() >= 3:
            R, t = _kabsch(P[best_inliers], Q[best_inliers])
            return R, t, best_inliers

    return best_R, best_t, best_inliers


def _irls_huber_kabsch(P: np.ndarray, Q: np.ndarray, R0: np.ndarray, t0: np.ndarray,
                       delta: float = 3.0, max_iters: int = 10) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Iteratively reweighted least squares with Huber loss. Returns R,t,weights."""
    R, t = R0.copy(), t0.copy()
    w = np.ones(P.shape[0], float)
    for _ in range(max_iters):
        R, t = _kabsch_weighted(P, Q, w)
        r = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
        # Huber weights
        w = np.where(r <= delta, 1.0, (delta / np.maximum(r, 1e-12)))
    return R, t, w


@dataclass
class RegularizationReport:
    group_id: str
    n_members: int
    r_mad: Dict[str, float] = field(default_factory=dict)
    theta_mad: Dict[str, float] = field(default_factory=dict)
    phi_mad: Dict[str, float] = field(default_factory=dict)
    r_median: Dict[str, float] = field(default_factory=dict)
    theta_median: Dict[str, float] = field(default_factory=dict)
    phi_median: Dict[str, float] = field(default_factory=dict)
    inlier_ratio: Optional[float] = None
    residual_rms: Optional[float] = None
    flags: List[str] = field(default_factory=list)  # e.g., ["LARGE_DEVIATION"]


# -----------------------------------------------------------------------------
# Public: robust Kabsch wrapper with diagnostics
# -----------------------------------------------------------------------------

def robust_kabsch(P: np.ndarray, Q: np.ndarray,
                  ransac_thresh: float = 4.0,
                  min_inlier_ratio: float = 0.55,
                  huber_delta: float = 3.0,
                  irls_max_iters: int = 10) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    n = P.shape[0]
    if n < 3:
        return np.eye(3), np.zeros(3), {"inlier_ratio": 1.0, "residual_rms": 0.0}

    R, t, inliers = _ransac_kabsch(P, Q, thresh=ransac_thresh,
                                   max_trials=200,
                                   min_inliers=max(3, int(min_inlier_ratio*n)))
    R, t, w = _irls_huber_kabsch(P[inliers], Q[inliers], R, t,
                                 delta=huber_delta, max_iters=irls_max_iters)
    # compute overall residuals under final model
    resid = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
    rms = float(np.sqrt(np.mean(resid**2))) if resid.size else 0.0
    report = {
        "inlier_ratio": float(inliers.mean()),
        "residual_rms": rms,
        "n": int(n),
        "n_inliers": int(inliers.sum()),
    }
    return R, t, report

def _stable_perp(u: np.ndarray) -> np.ndarray:
    """Return a deterministic unit vector ⟂ u."""
    u = _unit(u)
    cand = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = cand - u * np.dot(cand, u)
    n = np.linalg.norm(x)
    return x / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])


@dataclass
class GeometricSignature:
    """Geometric signature for interface template deduplication.

    Attributes:
        d_i: COM-to-interface distance for chain i.
        d_j: COM-to-interface distance for chain j.
        theta_i: Angle between COM-to-interface vector and COM-to-COM vector for chain i.
        theta_j: Angle between COM-to-interface vector and COM-to-COM vector for chain j.
    """
    d_i: float
    d_j: float
    theta_i: float
    theta_j: float

    def normalize(self, precision: int = 6) -> Tuple[float, float, float, float]:
        """Normalize signature to avoid floating-point errors.

        Args:
            precision: Number of decimal places for rounding.

        Returns:
            Normalized signature tuple.
        """
        return (
            round(self.d_i, precision),
            round(self.d_j, precision),
            round(self.theta_i, precision),
            round(self.theta_j, precision)
        )

    def is_similar_to(self, other: 'GeometricSignature',
                      distance_threshold: float, angle_threshold: float) -> bool:
        """Check if this signature is similar to another signature.

        Args:
            other: Other geometric signature to compare with.
            distance_threshold: Distance threshold for similarity.
            angle_threshold: Angle threshold for similarity.

        Returns:
            True if signatures are similar within thresholds.
        """
        return (abs(self.d_i - other.d_i) < distance_threshold and
                abs(self.d_j - other.d_j) < distance_threshold and
                abs(self.theta_i - other.theta_i) < angle_threshold and
                abs(self.theta_j - other.theta_j) < angle_threshold)

    def is_homotypic(self, distance_threshold: float, angle_threshold: float) -> bool:
        """Check if signature represents a homotypic interaction.

        Args:
            distance_threshold: Distance threshold for homodimer detection.
            angle_threshold: Angle threshold for homotypic detection.

        Returns:
            True if signature is nearly symmetric (homotypic).
        """
        return (abs(self.d_i - self.d_j) < distance_threshold and
                abs(self.theta_i - self.theta_j) < angle_threshold)
        
    def flipped(self):
        """Return the flipped geometric signature
        
        Returns: 
            the flipped geometric signature
        """
        return GeometricSignature(
            d_i=self.d_j,
            d_j=self.d_i,
            theta_i=self.theta_j,
            theta_j=self.theta_i,
        )


class TemplateBuilder:
    """Builder for molecular and interface templates.

    Generates molecular templates from chain groups and creates interface
    templates with geometric signatures for proper deduplication. Multiple
    interface types can exist between the same molecule types if they have
    different geometric signatures.

    Attributes:
        parser: PDB parser with structure data.
        coarse_grainer: Coarse-grainer with interface data.
        chain_grouper: Chain grouper with group information.
        hyperparams: Configuration parameters.
        units: Unit system for the model.
        workspace_manager: Workspace manager for logging (optional).
        molecule_templates: Dictionary of molecular templates by name.
        interface_templates: Dictionary of interface templates by name.
        interface_signatures: Dictionary mapping interface names to signatures.
        group_to_template: Mapping from group representative to template name.
        used_template_names: Set of already used template names.
        interface_type_counters: Counter for interface types between molecule pairs.
    """

    def __init__(self, parser: PDBParser, coarse_grainer: CoarseGrainer,
                 chain_grouper: ChainGrouper, hyperparams: PDBModelHyperparameters,
                 units: Optional[Units] = None, workspace_manager: Optional[WorkspaceManager] = None):
        """Initialize template builder.

        Args:
            parser: PDB parser with structure data.
            coarse_grainer: Coarse-grainer with processed data.
            chain_grouper: Chain grouper with group information.
            hyperparams: Configuration parameters.
            units: Unit system (defaults to standard units).
            workspace_manager: Workspace manager for logging (optional).
        """
        self.parser = parser
        self.coarse_grainer = coarse_grainer
        self.chain_grouper = chain_grouper
        self.hyperparams = hyperparams
        self.units = units or Units()
        self.interface_to_type_mapping = {}

        # Get workspace manager from parser if not provided
        if workspace_manager is None and hasattr(parser, 'workspace_manager'):
            self.workspace_manager = parser.workspace_manager
        else:
            self.workspace_manager = workspace_manager

        # Template storage
        self.molecule_templates: Dict[str, MoleculeType] = {}
        self.interface_templates: Dict[str, InterfaceType] = {}
        self.interface_signatures: Dict[str, GeometricSignature] = {}
        
        self._reference_geometry = {}  # Store reference geometry for each template

        # Mapping from group representative to template name
        self.group_to_template: Dict[str, str] = {}

        # Track used template names to avoid conflicts
        self.used_template_names: Set[str] = set()

        # Counter for interface types between molecule pairs
        self.interface_type_counters: Dict[Tuple[str, str], int] = {}
        
        # Track canonical orientation for homodimeric-heterotypic (HHT) signatures.
        # Keyed by (template_name, ordered_signature_tuple) -> {'canon_order': 'ij'|'ji', 'f': name_f, 'b': name_b, 'index': int}
        self.hht_catalog: Dict[
            Tuple[str, Tuple[float, float, float, float]],
            Dict[str, any]
        ] = {}



        # Build templates
        self._build_templates()

        # Debug interface detection
        self._debug_interface_detection()

        # Regularize across groups
        self._regularize_templates()
        
        # apply_prototype_method
        self.apply_prototype_method()

        # Detect steric clashes if enabled
        if self.hyperparams.steric_clash_mode == "auto":
            self._detect_steric_clashes()
        
        # debug print
        for intf_name, intf in self.interface_templates.items():
            intf_type = intf.get_name()
            self.workspace_manager.logger.info("Interface %s of type %s", intf_name, intf_type)
            
    def _sig_tuple(self, sig: GeometricSignature) -> Tuple[float, float, float, float]:
        """Ordered signature tuple (d_i, d_j, theta_i, theta_j) with rounding consistent with hyperparams."""
        return sig.normalize(self.hyperparams.signature_precision)

    def _sig_tuple_swapped(self, sig: GeometricSignature) -> Tuple[float, float, float, float]:
        """Swap i↔j in the signature (used to detect the same HHT under reversed ordering)."""
        di, dj, ti, tj = self._sig_tuple(sig)
        return (dj, di, tj, ti)
    
    def _match_existing_hht_signature(
        self,
        template_name: str,
        signature: "GeometricSignature",
        *,
        interface: "InterfaceString",        # current interface we’re trying to classify
        distance_threshold: float = 5.0,
        angle_threshold: float = 0.5,
    ):
        """
        Try to find an *existing* HHT entry for the same template whose signature
        this new one conforms to, under relaxed thresholds + residue validation.

        Returns:
            ("ij", catalog_entry)  -> matched stored orientation
            ("ji", catalog_entry)  -> matched reversed orientation
            (None, None)           -> no match
        """
        if not self.hht_catalog:
            return None, None

        for (tpl, stored_sig_tuple), cat in self.hht_catalog.items():
            if tpl != template_name:
                continue

            stored_iface = cat.get("exemplar_interface")
            stored_sig_obj = cat.get("exemplar_signature")

            # ------------------------------------------------------------------
            # CASE 1: old-style entry, only has the 4-tuple signature
            # ------------------------------------------------------------------
            if stored_iface is None or stored_sig_obj is None:
                d_i, d_j, th_i, th_j = stored_sig_tuple
                stored_sig_ij = GeometricSignature(d_i, d_j, th_i, th_j)

                # try direct
                if signature.is_similar_to(stored_sig_ij, distance_threshold, angle_threshold):
                    return "ij", cat

                # try flipped
                stored_sig_ji = GeometricSignature(d_j, d_i, th_j, th_i)
                if signature.is_similar_to(stored_sig_ji, distance_threshold, angle_threshold):
                    return "ji", cat

                # not this catalog entry → check next
                continue

            # ------------------------------------------------------------------
            # CASE 2: new-style entry, we have exemplar interface + signature
            # let the *robust* validator decide BOTH geometry + residue AND direction
            # ------------------------------------------------------------------
            ok, matched_order = self._is_homotypic_with_residue_validation(
                stored_iface,     # the one we stored when we first created the pair
                interface,        # the new one we’re trying to match
                stored_sig_obj,   # stored signature (canonical)
                signature,        # new signature
            )

            if ok:
                # matched_order is either "ij" or "ji"
                if matched_order == "ij":
                    return "ij", cat
                elif matched_order == "ji":
                    return "ji", cat
                else:
                    # very defensive: validator said ok but no order → treat as no match
                    if self.workspace_manager:
                        self.workspace_manager.logger.warning(
                            "HHT: homotypic validator returned True but no order; "
                            "catalog entry=%s, template=%s",
                            cat, template_name
                        )
                    return None, None

            # if not ok, just try next catalog entry

        # nothing matched
        return None, None




    
    def _ensure_hht_canonical_and_assign(
        self,
        interface: InterfaceString,
        template_name: str,
        signature: GeometricSignature,
    ) -> str:
        """
        HHT *with* conforming step.

        Flow:
        1. Try to match this signature to an existing HHT entry for `template_name`
           using relaxed thresholds (5 Å / 0.5 rad), in both (ij) and (ji) order.
           - If we match (ij): just return the existing 'f'
           - If we match (ji): flip the interface to canonical (ij), then return 'f'
        2. If no match: create a brand-new canonical pair (…f, …b), store it, return 'f'.

        This prevents the 1VYM case from spawning A_A_1*, A_A_2*, A_A_3* when
        all three copies are basically the same geometry.
        """
        # 1) try to reuse an existing HHT for this template
        match_dir, cat = self._match_existing_hht_signature(
            template_name,
            signature,
            interface=interface,
        )


        if match_dir is not None:
            # We already have a canonical pair for this geometry
            if match_dir == "ij":
                # interface is already in canonical order
                return cat["f"]
            else:  # match_dir == "ji"
                # we matched the reversed version, so normalize the actual interface
                interface.chain_i, interface.chain_j = interface.chain_j, interface.chain_i
                interface.coord_i, interface.coord_j = interface.coord_j, interface.coord_i
                return cat["f"]

        # 2) brand-new HHT pair (no existing one was close enough)
        sig_ord = self._sig_tuple(signature)   # canonical, rounded tuple

        # Get next interface index for this homodimer template pair
        template_pair = (template_name, template_name)
        next_index = self.interface_type_counters.get(tuple(sorted(template_pair)), 0) + 1
        self.interface_type_counters[tuple(sorted(template_pair))] = next_index

        # Construct names (A_A_#f / A_A_#b)
        name_f = f"{template_name}_{template_name}_{next_index}f"
        name_b = f"{template_name}_{template_name}_{next_index}b"

        # Build both sides (reuse nm conversion)
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[interface.chain_i]
        chain_j_data = self.coarse_grainer.get_coarse_grained_chains()[interface.chain_j]

        com_i_nm = self.parser.convert_coords_to_nm(chain_i_data.com)
        com_j_nm = self.parser.convert_coords_to_nm(chain_j_data.com)
        intf_i_nm = self.parser.convert_coords_to_nm(interface.coord_i)
        intf_j_nm = self.parser.convert_coords_to_nm(interface.coord_j)

        tmpl_f = InterfaceType(
            this_mol_type_name=template_name,
            partner_mol_type_name=template_name,
            interface_index=int(next_index),
            absolute_coord=intf_i_nm,
            local_coord=(intf_i_nm - com_i_nm),
            energy=interface.energy,
            tag="f",
        )
        tmpl_b = InterfaceType(
            this_mol_type_name=template_name,
            partner_mol_type_name=template_name,
            interface_index=int(next_index),
            absolute_coord=intf_j_nm,
            local_coord=(intf_j_nm - com_j_nm),
            energy=interface.energy,
            tag="b",
        )

        # Cross-link partners
        tmpl_f.partner_interface_type = tmpl_b
        tmpl_b.partner_interface_type = tmpl_f

        meta_base = {
            "original_chain_i": interface.chain_i,
            "original_chain_j": interface.chain_j,
            "interaction_type": "homodimeric_heterotypic",
            "geometric_signature": sig_ord,
            "this_side": "i",
            "complementary_interface": name_b,
        }
        tmpl_f.signature = {**meta_base, "interface_subtype": "barbed_end"}
        tmpl_b.signature = {
            **meta_base,
            "this_side": "j",
            "interface_subtype": "pointed_end",
            "complementary_interface": name_f,
        }

        # Store templates and signatures
        self.interface_templates[name_f] = tmpl_f
        self.interface_templates[name_b] = tmpl_b
        self.interface_signatures[name_f] = signature
        self.interface_signatures[name_b] = signature

        # Register in molecule template’s neighbor map
        self.molecule_templates[template_name].interfaces_neighbors_map[name_f] = template_name
        self.molecule_templates[template_name].interfaces_neighbors_map[name_b] = template_name

        # record canonical order as 'ij' for this ordered signature
        # right after you build tmpl_f / tmpl_b and before you write into catalog:
        feat_i = self._hht_side_features_from_interface(interface, "i")
        feat_j = self._hht_side_features_from_interface(interface, "j")

        # we must store the *canonicalized* interface as exemplar
        self.hht_catalog[(template_name, sig_ord)] = {
            "canon_order": "ij",
            "f": name_f,
            "b": name_b,
            "index": next_index,
            # NEW:
            "exemplar_interface": interface,     # this is now in canonical ij order
            "exemplar_signature": signature,     # GeometricSignature object
        }



        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "HHT: created canonical pair (%s, %s) for template %s; canon_order=ij; sig=%s",
                name_f,
                name_b,
                template_name,
                sig_ord,
            )

        return name_f

    def _hht_side_features_from_interface(self, interface: InterfaceString, side: str) -> dict:
        """
        Extracts the features we will use to tell 'this is the i-side' vs 'this is the j-side'
        for HHT re-use. Keep it small and hashable-ish.
        """
        if side == "i":
            residues = self._extract_interface_residues(interface, interface.chain_i, interface.coord_i)
            seq = interface.get_residue_sequence_i()
            comp = interface.get_residue_composition_i()
        else:
            residues = self._extract_interface_residues(interface, interface.chain_j, interface.coord_j)
            seq = interface.get_residue_sequence_j()
            comp = interface.get_residue_composition_j()

        # normalize to something comparable
        return {
            "res_types": sorted(list(self._extract_residue_types_enhanced(residues))),
            # sequences can be long; keep first N to avoid bloating catalog
            "seq_head": seq[:25] if seq else "",
            # composition keys, sorted, but we don't need counts for equality
            "comp_keys": sorted(list(comp.keys())) if comp else [],
        }


    def _hht_features_match(self, a: dict, b: dict) -> bool:
        """Very cheap equality/near-equality test for HHT sides."""
        if not a or not b:
            # if we cannot tell, do NOT claim it's the same orientation
            return False
        # must match residue types and the 'shape' of composition
        if a["res_types"] != b["res_types"]:
            return False
        if a["comp_keys"] != b["comp_keys"]:
            return False
        # sequence head is weaker signal; allow mismatch here if types+comp match
        return True

    def _get_reference_geometry(self, group: ChainGroup) -> Tuple[np.ndarray, float, Dict[str, np.ndarray]]:
        """Get reference geometry from group representative chain.
        
        Returns:
            Tuple of (reference_com, reference_radius, reference_interface_positions)
        """
        ref_chain_id = group.representative
        ref_chain = self.coarse_grainer.chains[ref_chain_id]
        
        # Get reference COM and radius
        ref_com = ref_chain.com.copy()
        ref_radius = ref_chain.radius
        
        # Get reference interface positions for this chain
        ref_interfaces = {}
        for interface in self.coarse_grainer.interfaces:
            if interface.chain_i == ref_chain_id:
                # This chain is the 'i' side of the interface
                key = f"{interface.chain_i}_{interface.chain_j}_i"
                ref_interfaces[key] = interface.coord_i.copy()
            elif interface.chain_j == ref_chain_id:
                # This chain is the 'j' side of the interface
                key = f"{interface.chain_j}_{interface.chain_i}_j"
                ref_interfaces[key] = interface.coord_j.copy()
        
        return ref_com, ref_radius, ref_interfaces

            
    def apply_prototype_method(self, groups: Optional[List[ChainGroup]] = None) -> None:
        """
        Prototype Method (public): treat each group's representative as the prototype,
        rigidly copy its geometry to all group members, then radial-regularize
        interfaces/sites across symmetry mates.

        This is a thin wrapper that exposes the behavior already implemented in
        `_regularize_group` without changing any of the underlying logic.

        Args:
            groups: Optional subset of ChainGroup objects to process. If None,
                    all groups from the ChainGrouper are used.
        """
        if groups is None:
            groups = self.chain_grouper.get_groups()

        # keep deterministic order (matches build path)
        groups_sorted = sorted(groups, key=lambda g: g.representative)

        for group in groups_sorted:
            try:
                self.regularize_group(group)
            except Exception as e:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Prototype method failed for group %s: %s",
                        getattr(group, "id", group.representative), str(e)
                    )


    def _generate_template_name(self, group: ChainGroup) -> str:
        """Generate a template name based on the chain group.

        Uses the representative chain name if unique, otherwise creates
        a descriptive name based on the group members.

        Args:
            group: Chain group to generate name for.

        Returns:
            Unique template name.
        """
        # Start with the representative chain name and force uppercase
        # This prevents filesystem collisions on case-insensitive systems (macOS)
        # where AA.mol and Aa.mol would be treated as the same file
        representative_name = group.representative.upper()

        # Check if the representative name is already used
        if representative_name not in self.used_template_names:
            self.used_template_names.add(representative_name)
            return representative_name

        # If representative name is taken, try variations
        base_name = representative_name

        # For groups with multiple members, try adding suffix
        if len(group.members) > 1:
            # Try adding "group" suffix (no underscore)
            candidate = f"{base_name}group"
            if candidate not in self.used_template_names:
                self.used_template_names.add(candidate)
                return candidate

        # If still conflicts, add numeric suffix (no underscore)
        counter = 1
        while True:
            candidate = f"{base_name}{counter}"
            if candidate not in self.used_template_names:
                self.used_template_names.add(candidate)
                return candidate
            counter += 1

    def _build_templates(self) -> None:
        """Build molecular and interface templates from chain groups."""
        groups = self.chain_grouper.get_groups()

        # Sort groups by representative name for deterministic ordering
        groups_sorted = sorted(groups, key=lambda g: g.representative)

        for group in groups_sorted:
            # Create molecular template from group representative
            self._build_molecule_template(group)

        # Build interface templates - process all interfaces
        self._build_all_interface_templates()

    def _build_molecule_template(self, group: ChainGroup) -> None:
        """Build molecular template from chain group using reference chain geometry."""
        
        representative_id = group.representative
        template_name = self._generate_template_name(group)

        # Get reference geometry from representative chain
        ref_com, ref_radius, ref_interfaces = self._get_reference_geometry(group)
        
        # Convert reference COM to nanometers
        ref_com_nm = self.parser.convert_coords_to_nm(ref_com)
        ref_radius_nm = ref_radius / 10.0  # Convert Å to nm directly
        
        # Create molecule template using reference geometry
        mol_template = MoleculeType(
            name=template_name,
            radius_nm=ref_radius_nm
        )
        
        # Set diffusion constants from radius (THIS WAS MISSING)
        mol_template.set_diffusion_constants_from_radius()
        
        # Store mapping from group representative to template name (THIS WAS MISSING)
        self.group_to_template[representative_id] = template_name
        
        # Store reference geometry for interface building
        self._reference_geometry[template_name] = {
            'com': ref_com,
            'com_nm': ref_com_nm,
            'radius': ref_radius,
            'interfaces': ref_interfaces
        }
        
        self.molecule_templates[template_name] = mol_template
        
        # Add group information to template for reference (THIS WAS MISSING)
        mol_template.signature = {
            'group_representative': representative_id,
            'group_members': group.members.copy(),
            'grouping_method': group.grouping_method,
            'original_chain_names': group.members.copy()
        }
        
        # Define ref coord (for safety)
        mol_template.ref1_local = np.array([1.0, 0.0, 0.0], float)
        mol_template.ref2_local = np.array([0.0, 0.0, 1.0], float)
        
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Created molecule template '%s' from reference chain '%s' with %d members",
                template_name, group.representative, len(group.members)
            )

    def _build_all_interface_templates(self) -> None:
        """Build interface templates for all detected interfaces."""
        interfaces = self.coarse_grainer.get_interfaces()

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Processing %d interfaces for template creation", len(
                    interfaces)
            )

        for interface in interfaces:
            # Get groups for both chains
            group_i = self.chain_grouper.get_group_for_chain(interface.chain_i)
            group_j = self.chain_grouper.get_group_for_chain(interface.chain_j)

            if not group_i or not group_j:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Missing group for interface %s <-> %s",
                        interface.chain_i, interface.chain_j
                    )
                continue

            # Get template names
            template_i = self.group_to_template.get(group_i.representative)
            template_j = self.group_to_template.get(group_j.representative)

            if not template_i or not template_j:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Missing template for interface %s <-> %s",
                        interface.chain_i, interface.chain_j
                    )
                continue

            # Calculate geometric signature
            signature = self._calculate_geometric_signature(interface)

            # Find or create interface type for this signature
            self._process_interface_with_signature(
                interface, template_i, template_j, signature)

    def _process_interface_with_signature(self, interface: InterfaceString,
                                      template_i: str, template_j: str,
                                      signature: GeometricSignature) -> None:
        """
        Process an interface and assign it to an interface type based on signature.
        HHT handling:
        - If template_i == template_j and NOT homotypic -> use canonical HHT registry:
            · Create canonical (…1f, …1b) for new signatures with order='ij'
            · If the same signature appears as (j,i), flip to canonical order and still assign 'f' to first
            · Persist the (possibly flipped) order in interface.chain_i / chain_j and coord_i / coord_j
        """
        # --- SPECIAL CASE: homodimeric heterotypic ---
        if template_i == template_j and not self._is_interface_homotypic(interface, template_i, template_j, signature):
            # HHT path, but now with conforming step inside the helper
            type_name_for_i = self._ensure_hht_canonical_and_assign(interface, template_i, signature)
            interface.interface_type = type_name_for_i
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "HHT: assigned %s to (first) %s after canonicalization; chain ordering is now (%s,%s)",
                    type_name_for_i,
                    interface.chain_i,
                    interface.chain_i,
                    interface.chain_j,
                )
            return


        # --- DEFAULT PATHS (heterodimeric or homodimeric-homotypic) remain unchanged ---
        matching_interface_name = self._find_matching_interface_type(template_i, template_j, signature)
        if matching_interface_name:
            interface.interface_type = matching_interface_name
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Interface %s <-> %s assigned to existing type: %s",
                    interface.chain_i, interface.chain_j, matching_interface_name
                )
            return

        # fall back to your existing creation logic for non-HHT
        new_interface_names = self._create_new_interface_type(interface, template_i, template_j, signature)
        if len(new_interface_names) == 1:
            interface.interface_type = new_interface_names[0]
        else:
            # heterodimeric path; preserve your selection logic
            interface_type_name = self._assign_interface_to_heterotypic_type(
                interface, new_interface_names, template_i, template_j
            )
            interface.interface_type = interface_type_name


    def _assign_interface_to_heterotypic_type(self, interface: InterfaceString,
                                              interface_type_names: List[str],
                                              template_i: str, template_j: str) -> str:
        """Assign an interface to one of the heterotypic interface types.

        Args:
            interface: Interface object to assign.
            interface_type_names: List of available interface type names.
            template_i: Template name for chain i.
            template_j: Template name for chain j.

        Returns:
            Selected interface type name.
        """
        if len(interface_type_names) != 2:
            return interface_type_names[0]

        # Get the two interface templates
        type_1_template = self.interface_templates[interface_type_names[0]]
        type_2_template = self.interface_templates[interface_type_names[1]]

        # Determine which type this interface should use based on the template's "this_side"
        # and the interface's chain information

        # Check which template corresponds to which side
        if (type_1_template.this_mol_type_name == template_i and
            hasattr(type_1_template, 'signature') and
                type_1_template.signature.get('this_side') == 'i'):
            return interface_type_names[0]
        elif (type_1_template.this_mol_type_name == template_j and
              hasattr(type_1_template, 'signature') and
              type_1_template.signature.get('this_side') == 'j'):
            return interface_type_names[0]
        else:
            return interface_type_names[1]

    def _store_interface_mapping(self, interface: InterfaceString, interface_type_name: str) -> None:
        """Store mapping from interface to interface type for later use in system building.

        Args:
            interface: Interface object.
            interface_type_name: Name of the interface type.
        """
        # Store mapping for system builder to use
        if not hasattr(self, 'interface_to_type_mapping'):
            self.interface_to_type_mapping = {}

        # For heterotypic interfaces, we need to store mappings for both sides
        # Create unique keys for both sides of the interface
        interface_key_i = f"{interface.chain_i}_{interface.chain_j}_{interface.coord_i[0]:.3f}_{interface.coord_i[1]:.3f}_{interface.coord_i[2]:.3f}"
        interface_key_j = f"{interface.chain_j}_{interface.chain_i}_{interface.coord_j[0]:.3f}_{interface.coord_j[1]:.3f}_{interface.coord_j[2]:.3f}"

        # Store the mapping
        self.interface_to_type_mapping[interface_key_i] = interface_type_name

        # For heterotypic interfaces, also store the reverse mapping if it's a different type
        if interface_type_name in self.interface_templates:
            interface_template = self.interface_templates[interface_type_name]
            if hasattr(interface_template, 'partner_interface_type') and interface_template.partner_interface_type:
                partner_type_name = interface_template.partner_interface_type.name if hasattr(
                    interface_template.partner_interface_type, 'name') else str(interface_template.partner_interface_type)
                self.interface_to_type_mapping[interface_key_j] = partner_type_name

    def _create_new_interface_type(self, interface: InterfaceString,
                                template_i: str, template_j: str,
                                signature: GeometricSignature) -> List[str]:
        """Create new interface type(s) for the given interface.

        Uses enhanced homotypic detection to determine whether to create
        a single shared template (homotypic) or separate templates (heterotypic).

        Args:
            interface: Interface object.
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature for this interface.

        Returns:
            List of created interface type names.
        """
        # Get next interface index for this template pair
        template_pair = tuple(sorted([template_i, template_j]))
        interface_index = self.interface_type_counters.get(template_pair, 0) + 1
        self.interface_type_counters[template_pair] = interface_index

        # Enhanced homotypic detection
        is_homodimeric_homotypic = self._is_interface_homotypic(
            interface, template_i, template_j, signature
        )

        created_names = []

        if is_homodimeric_homotypic:
            # Create single shared interface template
            interface_name = self._create_homotypic_interface_template(
                interface, template_i, signature, interface_index
            )
            created_names.append(interface_name)
            
            if self.workspace_manager:
                detection_method = self.hyperparams.homotypic_detection
                self.workspace_manager.logger.info(
                    "Created homotypic interface type %s using %s detection",
                    interface_name, detection_method
                )
        else:
            # Create separate interface templates for each side
            interface_names = self._create_heterotypic_interface_templates(
                interface, template_i, template_j, signature, interface_index
            )
            created_names.extend(interface_names)
            
            if self.workspace_manager:
                detection_method = self.hyperparams.homotypic_detection
                reason = "different templates" if template_i != template_j else "failed residue similarity"
                self.workspace_manager.logger.info(
                    "Created heterotypic interface types %s using %s detection (reason: %s)",
                    interface_names, detection_method, reason
                )

        return created_names

    def _is_interface_homotypic(self, interface: InterfaceString,
                            template_i: str, template_j: str,
                            signature: GeometricSignature) -> bool:
        """Determine if an interface should be treated as homotypic using enhanced detection.
        
        Args:
            interface: Interface object.
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature for this interface.
            
        Returns:
            True if interface should be treated as homotypic.
        """
        # Handle different detection modes
        if self.hyperparams.homotypic_detection == "off":
            return False  # Never treat as homotypic
        
        # Must be between same template types
        if template_i != template_j:
            return False
        
        # Check geometric signature symmetry
        is_geometrically_symmetric = signature.is_homotypic(
            self.hyperparams.homodimer_distance_threshold,
            self.hyperparams.homodimer_angle_threshold
        )
        
        if not is_geometrically_symmetric:
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Interface %s <-> %s failed geometric symmetry test: d_diff=%.3f, theta_diff=%.3f",
                    interface.chain_i, interface.chain_j,
                    abs(signature.d_i - signature.d_j),
                    abs(signature.theta_i - signature.theta_j)
                )
            return False
        
        # For "signature" mode, geometric symmetry is sufficient
        if self.hyperparams.homotypic_detection == "signature":
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Interface %s <-> %s passed signature-only homotypic test",
                    interface.chain_i, interface.chain_j
                )
            return True
        
        # For "auto" mode, also check residue similarity
        if self.hyperparams.homotypic_detection == "auto":
            return self._check_interface_residue_symmetry(interface)
        
        return False

    def _check_interface_residue_symmetry(self, interface: InterfaceString) -> bool:
        """Check if an interface has symmetric residue composition on both sides.
        
        For a truly homotypic interface, both sides should have similar residue
        compositions since they represent the same binding site type.
        
        Args:
            interface: Interface object to check.
            
        Returns:
            True if residue compositions are sufficiently similar.
        """
        try:
            # Extract interface residues for both sides
            residues_i = self._extract_interface_residues(
                interface, interface.chain_i, interface.coord_i
            )
            residues_j = self._extract_interface_residues(
                interface, interface.chain_j, interface.coord_j
            )
            
            # Calculate residue similarity between the two sides
            similarity = self._calculate_residue_similarity(residues_i, residues_j)
            
            # Check if similarity meets threshold
            is_symmetric = similarity >= self.hyperparams.homotypic_detection_residue_similarity_threshold
            
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Interface %s <-> %s residue symmetry check: similarity=%.3f, threshold=%.3f, result=%s",
                    interface.chain_i, interface.chain_j, similarity,
                    self.hyperparams.homotypic_detection_residue_similarity_threshold,
                    "PASS" if is_symmetric else "FAIL"
                )
                
                if not is_symmetric:
                    self.workspace_manager.logger.info(
                        "  Residues side i (%d): %s",
                        len(residues_i), sorted(list(residues_i))[:10]  # Show first 10
                    )
                    self.workspace_manager.logger.info(
                        "  Residues side j (%d): %s", 
                        len(residues_j), sorted(list(residues_j))[:10]  # Show first 10
                    )
            
            return is_symmetric
            
        except Exception as e:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Failed to check residue symmetry for interface %s <-> %s: %s",
                    interface.chain_i, interface.chain_j, str(e)
                )
            
            # Fallback to signature-only detection on error
            return True

    def _find_matching_interface_type(self, template_i: str, template_j: str,
                                      signature: GeometricSignature) -> Optional[str]:
        """Find existing interface type that matches the given signature.

        Args:
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature to match.

        Returns:
            Name of matching interface type, or None if no match found.
        """
        # Look for existing interface types between these templates
        for interface_name, existing_signature in self.interface_signatures.items():
            interface_template = self.interface_templates[interface_name]

            # Check if this interface is between the same template types (in either direction)
            templates_match = (
                (interface_template.this_mol_type_name == template_i and
                 interface_template.partner_mol_type_name == template_j) or
                (interface_template.this_mol_type_name == template_j and
                 interface_template.partner_mol_type_name == template_i)
            )

            if templates_match:
                # Use more relaxed thresholds for signature matching
                distance_threshold = 5.0  # 5 Angstroms tolerance
                angle_threshold = 0.5     # ~30 degrees tolerance

                # Check if signatures are similar
                if signature.is_similar_to(existing_signature, distance_threshold, angle_threshold):
                    if self.workspace_manager:
                        self.workspace_manager.logger.info(
                            "Found matching interface type %s for signature d_i=%.2f, d_j=%.2f, theta_i=%.3f, theta_j=%.3f",
                            interface_name, signature.d_i, signature.d_j, signature.theta_i, signature.theta_j
                        )
                    return interface_name

        return None

    def _create_homotypic_interface_template(self, interface: InterfaceString,
                                         template_name: str, signature: GeometricSignature,
                                         interface_index: int) -> str:
        """Create shared interface template for homodimeric homotypic interaction.
        
        Enhanced version that stores detailed residue information.

        Returns:
            Name of created interface template.
        """
        # Generate interface name using index
        interface_name = f"{template_name}_{template_name}_{interface_index}"

        # Convert coordinates to nanometers and calculate local coordinates
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[interface.chain_i]
        com_nm = self.parser.convert_coords_to_nm(chain_i_data.com)
        intf_coord_nm = self.parser.convert_coords_to_nm(interface.coord_i)
        local_coord_nm = intf_coord_nm - com_nm

        # Create interface template
        interface_template = InterfaceType(
            this_mol_type_name=template_name,
            partner_mol_type_name=template_name,
            interface_index=interface_index,
            absolute_coord=intf_coord_nm,
            local_coord=local_coord_nm,
            energy=interface.energy
        )

        # Enhanced metadata with detailed residue information
        interface_template.signature = {
            'original_chain_i': interface.chain_i,
            'original_chain_j': interface.chain_j,
            'interaction_type': 'homodimer',
            'contacting_residues_i': [str(r) for r in interface.residue_details_i],
            'contacting_residues_j': [str(r) for r in interface.residue_details_j],
            'residue_composition_i': interface.get_residue_composition_i(),
            'residue_composition_j': interface.get_residue_composition_j(),
            'residue_sequence_i': interface.get_residue_sequence_i(),
            'residue_sequence_j': interface.get_residue_sequence_j(),
            'geometric_signature': signature.normalize(self.hyperparams.signature_precision),
            'original_interface': interface,  # Store for future residue comparisons
            'homotypic_detection_method': self.hyperparams.homotypic_detection,
            'residue_similarity_used': self.hyperparams.homotypic_detection == "auto"
        }

        # Store template and signature
        self.interface_templates[interface_name] = interface_template
        self.interface_signatures[interface_name] = signature

        # Add to molecule template's interface map
        mol_template = self.molecule_templates[template_name]
        mol_template.interfaces_neighbors_map[interface_name] = template_name

        if self.workspace_manager:
            detection_info = f"method={self.hyperparams.homotypic_detection}"
            if self.hyperparams.homotypic_detection == "auto":
                detection_info += f", residue_threshold={self.hyperparams.homotypic_detection_residue_similarity_threshold}"
            
            # Show residue composition in log
            comp_i = interface.get_residue_composition_i()
            comp_j = interface.get_residue_composition_j()
            
            self.workspace_manager.logger.info(
                "Created homodimer interface type: %s (d_i=%.2f, d_j=%.2f, theta_i=%.3f, theta_j=%.3f, %s)",
                interface_name, signature.d_i, signature.d_j, signature.theta_i, signature.theta_j, detection_info
            )
            self.workspace_manager.logger.info(
                "  Residue composition - Chain %s: %s, Chain %s: %s",
                interface.chain_i, comp_i, interface.chain_j, comp_j
            )

        return interface_name


    def _create_heterotypic_interface_templates(self, interface: InterfaceString,
                                                template_i: str, template_j: str,
                                                signature: GeometricSignature,
                                                interface_index: int) -> List[str]:
        """Create separate interface templates for heterotypic interaction.
        
        For truly heterotypic interactions (different molecule types), creates:
        - A_B_1 and B_A_1 (bidirectional partners)
        
        For homodimeric heterotypic interactions (same molecule type but different residue composition), creates:
        - A_A_1f and A_A_1b (complementary partners for same molecule type)
        - forward and backward

        Returns:
            List of created interface template names.
        """
        created_names = []

        # Determine if this is a homodimeric heterotypic case (same molecule type but different interfaces)
        is_homodimeric_heterotypic = (template_i == template_j)

        if is_homodimeric_heterotypic:
            # For homodimeric heterotypic: create A_A_1 and A_A_2 (complementary interface types)
            interface_name_i = f"{template_i}_{template_j}_{interface_index}f"        # A_A_1f (e.g., barbed end)
            interface_name_j = f"{template_i}_{template_j}_{interface_index}b"        # A_A_1b (e.g., pointed end)
            
            # Update the counter to account for using two indices
            template_pair = tuple(sorted([template_i, template_j]))
            self.interface_type_counters[template_pair] = interface_index
            
        else:
            # For true heterotypic: create A_B_1 and B_A_1 (bidirectional)
            interface_name_i = f"{template_i}_{template_j}_{interface_index}"        # A_B_1
            interface_name_j = f"{template_j}_{template_i}_{interface_index}"        # B_A_1

        # Create interface template for side i
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[interface.chain_i]
        com_i_nm = self.parser.convert_coords_to_nm(chain_i_data.com)
        intf_i_nm = self.parser.convert_coords_to_nm(interface.coord_i)
        local_i_nm = intf_i_nm - com_i_nm

        # For homodimeric heterotypic, both interfaces have the same molecule types but different indices
        # For true heterotypic, they have different molecule types
        if is_homodimeric_heterotypic:
            # Both sides are the same molecule type, but represent different interface subtypes
            interface_template_i = InterfaceType(
                this_mol_type_name=template_i,
                partner_mol_type_name=template_i if is_homodimeric_heterotypic else template_j,
                interface_index=interface_index,
                absolute_coord=intf_i_nm,
                local_coord=local_i_nm,
                energy=interface.energy,
                tag=('f' if is_homodimeric_heterotypic else None),
            )
        else:
            interface_template_i = InterfaceType(
                this_mol_type_name=template_i,
                partner_mol_type_name=template_j,
                interface_index=interface_index,
                absolute_coord=intf_i_nm,
                local_coord=local_i_nm,
                energy=interface.energy
            )

        # Create interface template for side j
        chain_j_data = self.coarse_grainer.get_coarse_grained_chains()[interface.chain_j]
        com_j_nm = self.parser.convert_coords_to_nm(chain_j_data.com)
        intf_j_nm = self.parser.convert_coords_to_nm(interface.coord_j)
        local_j_nm = intf_j_nm - com_j_nm

        if is_homodimeric_heterotypic:
            # Both sides are the same molecule type, but represent different interface subtypes
            # build side j
            interface_template_j = InterfaceType(
                this_mol_type_name=template_j if not is_homodimeric_heterotypic else template_i,
                partner_mol_type_name=template_i if not is_homodimeric_heterotypic else template_i,
                interface_index=interface_index,         # ← keep SAME index
                absolute_coord=intf_j_nm,
                local_coord=local_j_nm,
                energy=interface.energy,
                tag=('b' if is_homodimeric_heterotypic else None),
            )
        else:
            interface_template_j = InterfaceType(
                this_mol_type_name=template_j,
                partner_mol_type_name=template_i,
                interface_index=interface_index,
                absolute_coord=intf_j_nm,
                local_coord=local_j_nm,
                energy=interface.energy
            )

        # Enhanced metadata including detection method information
        base_metadata = {
            'original_chain_i': interface.chain_i,
            'original_chain_j': interface.chain_j,
            'interaction_type': 'failed_homotypic' if is_homodimeric_heterotypic else 'heterotypic',
            'geometric_signature': signature.normalize(self.hyperparams.signature_precision),
            'homotypic_detection_method': self.hyperparams.homotypic_detection,
            'residue_composition_i': interface.get_residue_composition_i(),
            'residue_composition_j': interface.get_residue_composition_j(),
            'residue_sequence_i': interface.get_residue_sequence_i(),
            'residue_sequence_j': interface.get_residue_sequence_j(),
        }

        # Add side-specific metadata
        interface_template_i.signature = {
            **base_metadata,
            'this_side': 'i',
            'contacting_residues': [str(r) for r in interface.residue_details_i],
            'heterotypic_reason': self._get_heterotypic_reason(template_i, template_j, signature, interface),
            'interface_subtype': 'barbed_end' if is_homodimeric_heterotypic else 'side_i',
            'complementary_interface': interface_name_j,  # Store the partner interface name
        }

        interface_template_j.signature = {
            **base_metadata,
            'this_side': 'j',
            'contacting_residues': [str(r) for r in interface.residue_details_j],
            'heterotypic_reason': self._get_heterotypic_reason(template_i, template_j, signature, interface),
            'interface_subtype': 'pointed_end' if is_homodimeric_heterotypic else 'side_j',
            'complementary_interface': interface_name_i,  # Store the partner interface name
        }

        # Set up cross-references - this is crucial for homodimeric heterotypic cases
        if is_homodimeric_heterotypic:
            # For homodimeric heterotypic: A_A_1f partners with A_A_1b
            interface_template_i.partner_interface_type = interface_template_j
            interface_template_j.partner_interface_type = interface_template_i
            
            # Both reference the same molecule type
            interface_template_i.this_mol_type = self.molecule_templates[template_i]
            interface_template_i.partner_mol_type = self.molecule_templates[template_i]  # Same molecule type
            interface_template_j.this_mol_type = self.molecule_templates[template_i]
            interface_template_j.partner_mol_type = self.molecule_templates[template_i]  # Same molecule type
            
        else:
            # For true heterotypic: A_B_1 partners with B_A_1
            interface_template_i.partner_interface_type = interface_template_j
            interface_template_j.partner_interface_type = interface_template_i
            
            interface_template_i.this_mol_type = self.molecule_templates[template_i]
            interface_template_i.partner_mol_type = self.molecule_templates[template_j]
            interface_template_j.this_mol_type = self.molecule_templates[template_j]
            interface_template_j.partner_mol_type = self.molecule_templates[template_i]

        # Store templates and signatures
        self.interface_templates[interface_name_i] = interface_template_i
        self.interface_templates[interface_name_j] = interface_template_j
        self.interface_signatures[interface_name_i] = signature
        self.interface_signatures[interface_name_j] = signature

        # Add to molecule templates' interface maps
        self.molecule_templates[template_i].interfaces_neighbors_map[interface_name_i] = template_j
        self.molecule_templates[template_j].interfaces_neighbors_map[interface_name_j] = template_i

        created_names = [interface_name_i, interface_name_j]

        if self.workspace_manager:
            detection_info = f"method={self.hyperparams.homotypic_detection}"
            reason = self._get_heterotypic_reason(template_i, template_j, signature, interface)
            
            # Show residue composition in log
            comp_i = interface.get_residue_composition_i()
            comp_j = interface.get_residue_composition_j()
            
            interaction_type = "homodimeric heterotypic" if is_homodimeric_heterotypic else "heterotypic"
            
            if is_homodimeric_heterotypic:
                self.workspace_manager.logger.info(
                    "Created %s interface types: %s ↔ %s (complementary partners for %s, %s, reason=%s)",
                    interaction_type, interface_name_i, interface_name_j, template_i, detection_info, reason
                )
            else:
                self.workspace_manager.logger.info(
                    "Created %s interface types: %s, %s (d_i=%.2f, d_j=%.2f, theta_i=%.3f, theta_j=%.3f, %s, reason=%s)",
                    interaction_type, interface_name_i, interface_name_j, signature.d_i, signature.d_j, 
                    signature.theta_i, signature.theta_j, detection_info, reason
                )
            
            self.workspace_manager.logger.info(
                "  Residue composition - Chain %s: %s, Chain %s: %s",
                interface.chain_i, comp_i, interface.chain_j, comp_j
            )

        return created_names

    def get_interface_type_for_interface(self, interface: InterfaceString) -> Optional[str]:
        """Get the interface type name for a specific interface.

        Args:
            interface: Interface object.

        Returns:
            Interface type name if found, None otherwise.
        """
        if not hasattr(self, 'interface_to_type_mapping'):
            return None

        # Create the same key used during mapping
        interface_key = f"{interface.chain_i}_{interface.chain_j}_{interface.coord_i[0]:.3f}_{interface.coord_i[1]:.3f}_{interface.coord_i[2]:.3f}"
        return self.interface_to_type_mapping.get(interface_key)

    def _calculate_geometric_signature(self, interface: InterfaceString) -> GeometricSignature:
        """Calculate geometric signature for an interface.

        Args:
            interface: Interface to calculate signature for.

        Returns:
            GeometricSignature object.
        """
        # Get chain data
        chain_i_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_i]
        chain_j_data = self.coarse_grainer.get_coarse_grained_chains()[
            interface.chain_j]

        # COMs and interface coordinates
        com_i = chain_i_data.com
        com_j = chain_j_data.com
        intf_i = interface.coord_i
        intf_j = interface.coord_j

        # Calculate offset vectors
        offset_i = intf_i - com_i
        offset_j = intf_j - com_j
        com_vector_ij = com_j - com_i
        com_vector_ji = com_i - com_j

        # Calculate distances
        d_i = np.linalg.norm(offset_i)
        d_j = np.linalg.norm(offset_j)

        # Calculate angles (dot product normalized)
        if d_i > 0 and np.linalg.norm(com_vector_ij) > 0:
            theta_i = np.dot(offset_i, com_vector_ij) / \
                (d_i * np.linalg.norm(com_vector_ij))
            theta_i = np.arccos(np.clip(theta_i, -1.0, 1.0)
                                )  # Ensure valid range
        else:
            theta_i = 0.0

        if d_j > 0 and np.linalg.norm(com_vector_ji) > 0:
            theta_j = np.dot(offset_j, com_vector_ji) / \
                (d_j * np.linalg.norm(com_vector_ji))
            theta_j = np.arccos(np.clip(theta_j, -1.0, 1.0)
                                )  # Ensure valid range
        else:
            theta_j = 0.0

        return GeometricSignature(d_i, d_j, theta_i, theta_j)

    def _regularize_templates(self) -> None:
        for group in self.chain_grouper.get_groups():
            if len(group.members) > 0:
                self.regularize_group(group)
                
    # --- handle homodimeric heterotypic partners

    def _apply_iface_means_to_templates(self, iface_means_angstrom: Dict[str, float]) -> None:
        """Resize InterfaceType.local_coord to match enforced group means.
        Handles homodimeric heterotypic partners correctly.
        """
        for name, templ in self.interface_templates.items():
            # Get the interface key that matches how means were calculated
            interface_key = None
            
            # Check if this template has signature info to determine the key format
            if hasattr(templ, "signature") and isinstance(templ.signature, dict):
                side = templ.signature.get("this_side")  # 'i' or 'j'
                if side:
                    interface_key = f"{name}:{side}"
                else:
                    # Homodimer case - use template name directly
                    interface_key = name
            else:
                # Fallback - use template name
                interface_key = name
            
            # Apply the mean distance if we have it
            if interface_key in iface_means_angstrom:
                target_nm = iface_means_angstrom[interface_key] / 10.0  # Å → nm
                v = np.asarray(templ.local_coord, float)
                n = float(np.linalg.norm(v))
                if n > 1e-12:
                    templ.local_coord = (v / n) * target_nm
                    
                    # Update absolute_coord if it exists
                    if hasattr(templ, "absolute_coord") and templ.absolute_coord is not None:
                        templ.absolute_coord = templ.local_coord


    # -----------------------------------------------------------------------------
    # PATCH `_regularize_group` body (replace the alignment+consensus section)
    # -----------------------------------------------------------------------------

    def regularize_group(
        self,
        group,
        rmsd_warn_ang: float = 3.0,      # Å; warn if a member exceeds this Kabsch RMSD
        com_shift_warn_ang: float = 5.0,  # Å; warn if COM move exceeds this
        use_normals: bool = True,
        normal_lever: float = 5.0,        # Å; lever arm for normal pseudo-points
    ):
        """
        Enforce identical *local* geometry for all members of a species (group),
        robust to crystal discrepancies.

        Steps
        -----
        1) Populate per-chain interfaces (adds σ vectors to partners).
        2) For each member, compute Kabsch (R, t) that aligns its landmarks to the
        representative. Warn if RMSD is large.
        3) Bring every site/interface local vector into the *representative frame*,
        collect robust medians:
            - sites: median of local vectors
            - interfaces: median (r, θ, φ) relative to *per-key* mean σ
        4) Build canonical local vectors from medians.
        5) For every member:
            - world targets = COM + R @ canonical_local
            - choose COM_new with a capped, strength-weighted shift (blending an
            angular-consistent least-squares shift with the naive landmark mean).
        6) Warn if COM shift is large.

        Notes
        -----
        * "Identical local geometry" means the set {local vectors} is the same for
        all members, up to the member's rigid body orientation R from Kabsch.
        * We never diverge bindings by geometry classes; chemistry decides types.
        """

        # -------- master switch: allow disabling regularization entirely ----------
        lam_strength = float(getattr(self.hyperparams, "template_regularization_strength", 0.5))
        if lam_strength < 0:
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Template regularization disabled (template_regularization_strength < 0)."
                )
            return

        com_shift_cap = float(getattr(self.hyperparams, "com_shift_cap_ang", 6.0))  # Å
        com_shift_reg = float(getattr(self.hyperparams, "com_shift_reg", 1e-6))     # Å^2

        chains = self.coarse_grainer.get_coarse_grained_chains()
        members = list(group.members)
        if not members:
            return

        # 1) materialize σ (“toward partner COM”) and stable keys on each chain
        _materialize_interfaces_on_chains_from_global(
            self.coarse_grainer, chains, set(members)
        )

        rep_id = group.representative
        rep = chains[rep_id]

        # ----------------------------- local helpers -----------------------------
        def _sites_dict(obj):
            if hasattr(obj, "sites") and isinstance(obj.sites, dict):
                return {str(k): np.asarray(v, float) for k, v in obj.sites.items()}
            out = {}
            for k in ("aa1", "aa2", "aa3", "aa4"):
                if hasattr(obj, k):
                    out[k] = np.asarray(getattr(obj, k), float)
            return out

        def _rich_ifaces(obj):
            out = []
            for i, raw in enumerate(getattr(obj, "interfaces", []) or []):
                if isinstance(raw, dict):
                    coord = raw.get("coord")
                    normal = raw.get("normal", None)
                    key = raw.get("key", None)
                    sigma = raw.get("sigma", None)
                else:
                    coord = getattr(raw, "coord", None)
                    normal = getattr(raw, "normal", None)
                    key = getattr(raw, "key", None)
                    sigma = getattr(raw, "sigma", None)
                if coord is None:
                    continue
                out.append({
                    "coord": np.asarray(coord, float),
                    "normal": None if normal is None else np.asarray(normal, float),
                    "sigma": None if sigma is None else np.asarray(sigma, float),
                    "key": str(key) if key is not None else f"index:{i}",
                    "raw": raw,
                })
            return out

        def _u(x):
            n = float(np.linalg.norm(x))
            return x / n if n > 1e-12 else np.array([1.0, 0.0, 0.0], float)

        def _angles_from_sigma(v, s_hat):
            """Return (r, theta, phi) using a *per-key* sigma direction s_hat in the same frame as v."""
            vhat = _u(v); shat = _u(s_hat)
            r = float(np.linalg.norm(v))
            theta = float(np.arccos(np.clip(np.dot(vhat, shat), -1.0, 1.0)))
            # φ: use projected σ as the canonical normal (your original convention)
            s_perp = shat - vhat*np.dot(shat, vhat)
            if np.linalg.norm(s_perp) < 1e-8:
                return r, theta, 0.0
            x0 = _u(s_perp); y0 = _u(np.cross(vhat, x0))
            phi = float(np.arctan2(np.dot(_u(s_perp), y0), np.dot(_u(s_perp), x0)))
            return r, theta, phi

        def _ideal_vec_from_sigma(s_hat, r, theta, phi):
            """Build canonical vector with length r at angles (theta,phi) w.r.t. the given σ=s_hat."""
            z = _u(s_hat)
            x = np.array([1, 0, 0], float) if abs(z[0]) < 0.9 else np.array([0, 1, 0], float)
            x = _u(x - z * np.dot(x, z))
            # rotate z toward x by θ (deterministic azimuth=0)
            c, s = np.cos(theta), np.sin(theta)
            vdir = c * z + s * x
            # φ is around v (here φ only affects the canonical normal, but we keep interface vector = vdir)
            return r * vdir

        # ---------- 2) member→rep transforms with (optional) normal pseudo-points ----------
        member_RT = {rep_id: (np.eye(3), np.zeros(3))}
        member_RMSD = {rep_id: 0.0}
        underspecified = {}   # mid -> bool
        priority = {}         # mid -> "high" | "low"

        for mid in members:
            if mid == rep_id:
                priority[mid] = "high"
                underspecified[mid] = False
                continue

            P, Q = _collect_points_for_fit_with_normals(
                rep, chains[mid],
                normal_lever if use_normals else 0.0
            )

            if P.shape[0] >= 3:
                R, t = _kabsch_transform(P, Q)      # rep→member
                R_m2r = R.T
                t_m2r = -(R.T @ t)
                member_RT[mid] = (R, t)
                member_RMSD[mid] = _compute_fit_rmsd(P, Q, R, t)
                chains[mid].__RT_member_to_rep__ = (R_m2r, t_m2r)
                underspecified[mid] = False
                priority[mid] = "high"
            else:
                R = np.eye(3)
                t = np.asarray(chains[mid].com, float) - np.asarray(rep.com, float)
                member_RT[mid] = (R, t)
                member_RMSD[mid] = 0.0
                chains[mid].__RT_member_to_rep__ = (R.T, -(R.T @ t))
                underspecified[mid] = True
                priority[mid] = "low"

            if self.workspace_manager and member_RMSD[mid] > rmsd_warn_ang:
                self.workspace_manager.logger.warning(
                    "Large Kabsch RMSD (%.2f Å) for member %s in group %s; will use high-priority averaging.",
                    member_RMSD[mid], mid, getattr(group, "id", group.representative)
                )

        # ---------- 3) per-key mean σ in the representative frame ----------
        def _group_mean_sigmas_by_key():
            acc = {}  # key -> list of unit sigma vectors in rep frame
            for mid in members:
                c = chains[mid]
                Rinv, tinv = getattr(c, "__RT_member_to_rep__", (np.eye(3), np.zeros(3)))
                for e in _rich_ifaces(c):
                    s = e["sigma"]
                    if s is None:
                        continue
                    s = _u(np.asarray(s, float))
                    s_rep = Rinv @ s
                    acc.setdefault(e["key"], []).append(s_rep)
            means = {}
            for k, L in acc.items():
                m = np.sum(np.vstack(L), axis=0)
                nm = np.linalg.norm(m)
                if nm > 1e-12:
                    means[k] = m / nm
            return means

        mean_sigma_rep = _group_mean_sigmas_by_key()
        z_fallback = np.array([0.0, 0.0, 1.0])  # only used if a key has no σ data

        # ---------- collect local vectors in the rep frame & medians ----------
        site_acc_hi, site_acc_lo = {}, {}
        if_acc_hi,   if_acc_lo   = {}, {}

        # representative (HIGH)
        for name, p in _sites_dict(rep).items():
            site_acc_hi.setdefault(name, []).append(np.asarray(p, float) - np.asarray(rep.com, float))

        for e in _rich_ifaces(rep):
            v = e["coord"] - np.asarray(rep.com, float)
            s_hat = mean_sigma_rep.get(e["key"], z_fallback)
            r, th, ph = _angles_from_sigma(v, s_hat)
            if_acc_hi.setdefault(e["key"], []).append((r, th, ph))

        # other members → rep frame; bucket by priority
        for mid in members:
            if mid == rep_id:
                continue
            Rinv, tinv = chains[mid].__RT_member_to_rep__
            com_m = np.asarray(chains[mid].com, float)

            bucket_sites = site_acc_hi if priority[mid] == "high" else site_acc_lo
            bucket_if    = if_acc_hi   if priority[mid] == "high" else if_acc_lo

            for name, p in _sites_dict(chains[mid]).items():
                v_world = np.asarray(p, float) - com_m
                v_rep   = Rinv @ v_world
                bucket_sites.setdefault(name, []).append(v_rep)

            for e in _rich_ifaces(chains[mid]):
                v_world = e["coord"] - com_m
                v_rep   = Rinv @ v_world
                s_hat   = mean_sigma_rep.get(e["key"], z_fallback)
                r, th, ph = _angles_from_sigma(v_rep, s_hat)
                bucket_if.setdefault(e["key"], []).append((r, th, ph))

        def _median_vec(L):
            return np.median(np.vstack(L), axis=0)

        # site medians
        all_site_keys = set(site_acc_hi) | set(site_acc_lo)
        site_canon_rep = {}
        for k in sorted(all_site_keys):
            if k in site_acc_hi and len(site_acc_hi[k]) > 0:
                site_canon_rep[k] = _median_vec(site_acc_hi[k])
            elif k in site_acc_lo and len(site_acc_lo[k]) > 0:
                site_canon_rep[k] = _median_vec(site_acc_lo[k])

        # interface medians (r,θ,φ) - FIXED KEY GENERATION
        all_if_keys = set(if_acc_hi) | set(if_acc_lo)
        if_pose_canon = {}
        for k in sorted(all_if_keys):
            src = if_acc_hi.get(k) if (k in if_acc_hi and len(if_acc_hi[k]) > 0) else if_acc_lo.get(k, [])
            if not src:
                continue
            arr = np.asarray(src, float)
            r_med  = float(np.median(arr[:, 0]))
            th_med = float(np.median(arr[:, 1]))
            ph_med = float(np.median(arr[:, 2]))
            if_pose_canon[k] = (r_med, th_med, ph_med)
        
        # synthesize canonical interface vectors (rep frame) using per-key σ
        if_vec_canon_rep = {}
        for key, (r_med, th_med, ph_med) in if_pose_canon.items():
            s_hat = mean_sigma_rep.get(key, z_fallback)
            if_vec_canon_rep[key] = _ideal_vec_from_sigma(s_hat, r_med, th_med, ph_med)
            
        # FIXED: Apply interface means to templates with consistent key mapping
        iface_means_for_templates = {}
        for key, (r_med, th_med, ph_med) in if_pose_canon.items():
            # Convert interface key to template key format
            iface_means_for_templates[key] = r_med
        self._apply_iface_means_to_templates(iface_means_for_templates)

        # ---------- 4&5) write back to HIGH-priority members with capped COM shift ----------
        for mid in members:
            mem = chains[mid]
            R_rep_to_mem, t_rep_to_mem = member_RT.get(mid, (np.eye(3), np.zeros(3)))

            if priority[mid] == "low":
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "Group %s: member %s is LOW priority — skipping geometry overwrite.",
                        getattr(group, "id", group.representative), mid
                    )
                continue

            # targets and current landmarks for COM solve
            contributions = []

            for name, v_rep in site_canon_rep.items():
                cur = _sites_dict(mem).get(name, None)
                if cur is None:
                    continue
                v_world = R_rep_to_mem @ v_rep
                contributions.append((np.asarray(cur, float), v_world, ("site", name)))

            mem_if_list = _rich_ifaces(mem)
            mem_if_by_key = {e["key"]: e for e in mem_if_list}
            for key, v_rep in if_vec_canon_rep.items():
                e = mem_if_by_key.get(key)
                if e is None:
                    continue
                v_world = R_rep_to_mem @ v_rep
                contributions.append((e["coord"], v_world, ("iface", key)))

            # naive average shift (current - target)
            if contributions:
                naive = np.mean(np.vstack([cw - vw for (cw, vw, _tag) in contributions]), axis=0)
            else:
                naive = np.zeros(3)

            # angular-consistent LS shift using per-key σ (if available)
            # requires _solve_best_com_shift_for_member() defined elsewhere in this module
            mean_dirs_by_key = {k: _u(v) for k, v in mean_sigma_rep.items()}
            c_opt = _solve_best_com_shift_for_member(mem, mean_dirs_by_key, reg_lambda=com_shift_reg)

            if c_opt is None:
                proposal = naive
            else:
                proposal = 0.5 * naive + 0.5 * c_opt  # blend two estimates

            # cap magnitude
            d = float(np.linalg.norm(proposal))
            if d > com_shift_cap:
                proposal = proposal * (com_shift_cap / max(d, 1e-12))

            # blend with strength parameter
            com_target = np.asarray(mem.com, float) + proposal
            com_new = (1.0 - lam_strength) * np.asarray(mem.com, float) + lam_strength * com_target

            com_shift = float(np.linalg.norm(com_new - np.asarray(mem.com, float)))
            if self.workspace_manager and com_shift > com_shift_warn_ang:
                self.workspace_manager.logger.warning(
                    "Large COM shift (%.2f Å) for member %s in group %s during symmetry regularization.",
                    com_shift, mid, getattr(group, "id", group.representative)
                )

            # Rotate template-local reference directions into this member's world frame
            R = R_rep_to_mem  # rotation from representative/template frame to this member

            # @TODO
            #ref1_world = R @ tpl_ref1_local
            #ref2_world = R @ tpl_ref2_local

            # Defensive normalization & Gram–Schmidt to keep ref2 orthogonal to ref1
            #n1 = np.linalg.norm(ref1_world)
            #ref1_world = ref1_world / (n1 if n1 > 1e-12 else 1.0)

            # Make ref2 orthogonal to ref1 and normalize
            #ref2_world = ref2_world - ref1_world * np.dot(ref1_world, ref2_world)
            #n2 = np.linalg.norm(ref2_world)
            #if n2 > 1e-12:
            #    ref2_world = ref2_world / n2
            #else:
                # fallback if degenerate: choose any vector not parallel to ref1
            #    alt = np.array([1.0, 0.0, 0.0]) if abs(ref1_world[0]) < 0.9 else np.#array([0.0, 1.0, 0.0])
            #    ref2_world = alt - ref1_world * np.dot(ref1_world, alt)
            #    ref2_world = ref2_world / max(np.linalg.norm(ref2_world), 1e-12)

            # Store on the member (instance). Overwrite if present.
            #setattr(mem, "ref1", ref1_world.astype(float))
            #setattr(mem, "ref2", ref2_world.astype(float))

            # place COM and enforce identical local geometry
            mem.com = com_new
            _enforce_identical_local_geometry_after_com(
                mem,
                R_rep_to_mem,
                com_new,
                site_canon_rep,
                if_vec_canon_rep
            )

        # ---------- 6) summary log ----------
        if self.workspace_manager:
            used_high = any(len(v) for v in if_acc_hi.values()) or any(len(v) for v in site_acc_hi.values())
            n_hi = sum(1 for m in members if priority[m] == "high")
            n_lo = len(members) - n_hi
            self.workspace_manager.logger.info(
                "Regularized group %s: %d high-priority, %d low-priority members. Averaging used %s set (λ=%.2f, cap=%.1f Å).",
                getattr(group, "id", group.representative),
                n_hi, n_lo,
                "HIGH" if used_high else "LOW",
                lam_strength, com_shift_cap
            )

    def _debug_interface_detection(self) -> None:
        """Debug interface detection to ensure multiple types are found."""
        if self.workspace_manager:
            interface_types_by_template_pair = {}
            
            for name, template in self.interface_templates.items():
                pair = tuple(sorted([template.this_mol_type_name, template.partner_mol_type_name]))
                if pair not in interface_types_by_template_pair:
                    interface_types_by_template_pair[pair] = []
                interface_types_by_template_pair[pair].append(name)
            
            self.workspace_manager.logger.info("Interface types by template pair:")
            for pair, types in interface_types_by_template_pair.items():
                self.workspace_manager.logger.info(f"  {pair[0]} <-> {pair[1]}: {len(types)} types: {types}")
                
            # Check signatures
            self.workspace_manager.logger.info("Interface signatures:")
            for name, sig in self.interface_signatures.items():
                self.workspace_manager.logger.info(
                    f"  {name}: d_i={sig.d_i:.2f}, d_j={sig.d_j:.2f}, theta_i={sig.theta_i:.3f}, theta_j={sig.theta_j:.3f}"
                )

    def _compute_rigid_transform(self, reference_data: CoarseGrainedChain,
                                 member_data: CoarseGrainedChain) -> np.ndarray:
        """Compute rigid transform between two chains.

        Args:
            reference_data: Reference chain data.
            member_data: Member chain data to transform.

        Returns:
            4x4 transformation matrix.
        """
        # Get Cα coordinates for both chains
        ref_coords = self.parser.get_chain_data(
            reference_data.chain_id)['ca_coords']
        mem_coords = self.parser.get_chain_data(
            member_data.chain_id)['ca_coords']

        if len(ref_coords) == len(mem_coords) and len(ref_coords) > 0:
            # Use Superimposer to compute transformation
            sup = Superimposer()
            try:
                sup.set_atoms(ref_coords, mem_coords)
                # Extract rotation matrix and translation vector
                rotation = sup.rotran[0]
                translation = sup.rotran[1]

                # Build 4x4 transformation matrix
                transform = np.eye(4)
                transform[:3, :3] = rotation
                transform[:3, 3] = translation
                return transform
            except Exception:
                # Fallback to identity transform
                return np.eye(4)

        return np.eye(4)

    def _detect_steric_clashes(self) -> None:
        """Detect steric clashes between interface templates."""
        # This is a simplified implementation
        # Full implementation would require detailed Cα clash checking

        for template_name, interface_template in self.interface_templates.items():
            # For each interface, check for potential clashes with other interfaces
            # on the same molecule type
            mol_type = interface_template.this_mol_type
            if mol_type:
                other_interfaces = [
                    name for name, intf in self.interface_templates.items()
                    if (intf.this_mol_type == mol_type and name != template_name)
                ]

                # Simple distance-based clash detection
                # (Real implementation would use detailed atomic coordinates)
                for other_name in other_interfaces:
                    other_intf = self.interface_templates[other_name]
                    distance = np.linalg.norm(
                        interface_template.local_coord - other_intf.local_coord
                    )

                    # If interfaces are very close, mark as mutually exclusive
                    if distance < 0.5:  # 0.5 nm threshold (adjustable)
                        interface_template.required_free.append(other_name)
                        other_intf.required_free.append(template_name)

    def get_molecule_templates(self) -> Dict[str, MoleculeType]:
        """Get all molecular templates.

        Returns:
            Dictionary mapping template names to MoleculeType objects.
        """
        return self.molecule_templates.copy()

    def get_interface_templates(self) -> Dict[str, InterfaceType]:
        """Get all interface templates.

        Returns:
            Dictionary mapping interface names to InterfaceType objects.
        """
        return self.interface_templates.copy()
    
    def _get_interface_key_for_template(self, template_name: str, template: InterfaceType) -> str:
        """Get consistent interface key for template that matches regularization collection."""
        if hasattr(template, "signature") and isinstance(template.signature, dict):
            side = template.signature.get("this_side")
            if side:
                return f"{template_name}:{side}"
        
        # For homodimer or templates without side info
        return template_name

    def get_template_name_for_group(self, group_representative: str) -> Optional[str]:
        """Get the template name for a group representative.

        Args:
            group_representative: Chain ID of group representative.

        Returns:
            Template name if found, None otherwise.
        """
        return self.group_to_template.get(group_representative)

    def get_chain_name_mapping(self) -> Dict[str, str]:
        """Get mapping from original chain names to template names.

        Returns:
            Dictionary mapping chain IDs to their template names.
        """
        chain_to_template = {}
        for group in self.chain_grouper.get_groups():
            template_name = self.group_to_template.get(group.representative)
            if template_name:
                for chain_id in group.members:
                    chain_to_template[chain_id] = template_name
        return chain_to_template

    def get_summary(self) -> Dict[str, any]:
        """Get summary of template building results.

        Returns:
            Dictionary with template statistics and naming information.
        """
        # Count interface types by molecule pair
        interface_type_counts = {}
        for interface_name, interface_template in self.interface_templates.items():
            mol_pair = tuple(sorted([interface_template.this_mol_type_name,
                                     interface_template.partner_mol_type_name]))
            interface_type_counts[mol_pair] = interface_type_counts.get(
                mol_pair, 0) + 1

        return {
            "num_molecule_templates": len(self.molecule_templates),
            "num_interface_templates": len(self.interface_templates),
            "molecule_templates": list(self.molecule_templates.keys()),
            "interface_templates": list(self.interface_templates.keys()),
            "interface_type_counts_by_molecule_pair": interface_type_counts,
            "group_to_template_mapping": self.group_to_template.copy(),
            "chain_name_mapping": self.get_chain_name_mapping(),
            "template_naming_strategy": "signature_based"
        }


########################
# Key additions for custom homotypic interface detection
########################


    def _extract_interface_residues(self, interface: InterfaceString, chain_id: str, 
                                interface_coord: np.ndarray) -> Set[str]:
        """Extract residues with amino acid names from enhanced interface data.
        
        Args:
            interface: Interface object containing detailed residue information.
            chain_id: Chain ID to extract residues from.
            interface_coord: Interface coordinate center (not used with enhanced data).
            
        Returns:
            Set of residue identifiers in format "ALA123", "GLY45", etc.
        """
        interface_residues = set()
        
        try:
            if chain_id == interface.chain_i:
                # Use detailed residue information from enhanced InterfaceString
                for residue_info in interface.residue_details_i:
                    residue_id = f"{residue_info.residue_name}{residue_info.residue_id}"
                    interface_residues.add(residue_id)
                    
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "Extracted %d residues for chain %s from interface.residue_details_i", 
                        len(interface_residues), chain_id
                    )
                    
            elif chain_id == interface.chain_j:
                # Use detailed residue information from enhanced InterfaceString
                for residue_info in interface.residue_details_j:
                    residue_id = f"{residue_info.residue_name}{residue_info.residue_id}"
                    interface_residues.add(residue_id)
                    
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "Extracted %d residues for chain %s from interface.residue_details_j", 
                        len(interface_residues), chain_id
                    )
            else:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Chain %s not found in interface %s <-> %s", 
                        chain_id, interface.chain_i, interface.chain_j
                    )
                    
        except Exception as e:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Failed to extract residues from enhanced interface object for chain %s: %s", 
                    chain_id, str(e)
                )
            
            # Fallback to legacy residue extraction if enhanced data fails
            interface_residues = self._extract_residues_fallback(interface, chain_id)
        
        # Debug: show sample of residues and their types
        if self.workspace_manager and interface_residues:
            sample_residues = sorted(list(interface_residues))[:5]  # Show first 5
            self.workspace_manager.logger.info(
                "Sample residues for chain %s: %s", 
                chain_id, sample_residues
            )
        
        return interface_residues
    def _extract_residues_fallback(self, interface: InterfaceString, chain_id: str) -> Set[str]:
        """Fallback residue extraction using legacy residue IDs.
        
        Args:
            interface: Interface object.
            chain_id: Chain ID to extract residues from.
            
        Returns:
            Set of residue identifiers (may be less detailed than enhanced version).
        """
        interface_residues = set()
        
        try:
            if chain_id == interface.chain_i and hasattr(interface, 'residues_i'):
                # Convert legacy residue IDs to strings
                for residue_id in interface.residues_i:
                    interface_residues.add(str(residue_id))
                    
            elif chain_id == interface.chain_j and hasattr(interface, 'residues_j'):
                # Convert legacy residue IDs to strings
                for residue_id in interface.residues_j:
                    interface_residues.add(str(residue_id))
                    
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Fallback extraction: got %d residues for chain %s", 
                    len(interface_residues), chain_id
                )
                
        except Exception as e:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Fallback residue extraction also failed for chain %s: %s", 
                    chain_id, str(e)
                )
        
        return interface_residues

    def _calculate_residue_similarity(self, residues1: Set[str], residues2: Set[str]) -> float:
        """Calculate similarity between two sets of interface residues.
        
        Uses Jaccard similarity: |intersection| / |union|
        Now handles the enhanced format "ALA123", "GLY45", etc.
        
        Args:
            residues1: First set of residue identifiers.
            residues2: Second set of residue identifiers.
            
        Returns:
            Similarity score between 0.0 and 1.0.
        """
        if not residues1 and not residues2:
            return 1.0  # Both empty = perfect similarity
        
        if not residues1 or not residues2:
            return 0.0  # One empty, one not = no similarity
        
        # Extract residue types from enhanced format
        types1 = self._extract_residue_types_enhanced(residues1)
        types2 = self._extract_residue_types_enhanced(residues2)
        
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Residue similarity calculation: types1=%s, types2=%s", 
                sorted(list(types1))[:5], sorted(list(types2))[:5]  # Show first 5 for debugging
            )
        
        # Calculate Jaccard similarity
        intersection = len(types1.intersection(types2))
        union = len(types1.union(types2))
        
        if union == 0:
            return 1.0
        
        return intersection / union

    def _extract_residue_types_enhanced(self, residues: Set[str]) -> Set[str]:
        """Extract residue types from enhanced residue identifiers.
        
        Handles the new format "ALA123" -> "ALA", "GLY45" -> "GLY", etc.
        
        Args:
            residues: Set of residue identifiers in enhanced format.
            
        Returns:
            Set of three-letter amino acid codes.
        """
        types = set()
        
        for res in residues:
            try:
                if isinstance(res, str) and len(res) >= 3:
                    # Enhanced format: "ALA123" -> extract "ALA"
                    # Find where numbers start
                    for i, char in enumerate(res):
                        if char.isdigit():
                            amino_acid = res[:i].upper()
                            if amino_acid:
                                types.add(amino_acid)
                            break
                    else:
                        # No numbers found, assume entire string is amino acid
                        types.add(res.upper())
                else:
                    # Fallback for other formats
                    types.add(str(res).upper())
                    
            except Exception as e:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Failed to extract residue type from %s: %s", res, str(e)
                    )
                # Fallback - use string representation
                types.add(str(res))
        
        return types

    def _check_interface_residue_symmetry(self, interface: InterfaceString) -> bool:
        """Check if an interface has symmetric residue composition on both sides.
        
        Enhanced version that uses detailed residue information from InterfaceString.
        
        Args:
            interface: Interface object to check.
            
        Returns:
            True if residue compositions are sufficiently similar.
        """
        try:
            # Extract interface residues using enhanced data
            residues_i = self._extract_interface_residues(
                interface, interface.chain_i, interface.coord_i
            )
            residues_j = self._extract_interface_residues(
                interface, interface.chain_j, interface.coord_j
            )
            
            # Check if we have any residue data
            if not residues_i and not residues_j:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "No residue data available for interface %s <-> %s, falling back to signature-only",
                        interface.chain_i, interface.chain_j
                    )
                return True  # Fallback to signature-only behavior
            
            # Calculate residue similarity between the two sides
            similarity = self._calculate_residue_similarity(residues_i, residues_j)
            
            # Check if similarity meets threshold
            is_symmetric = similarity >= self.hyperparams.homotypic_detection_residue_similarity_threshold
            
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Interface %s <-> %s residue symmetry check: similarity=%.3f, threshold=%.3f, result=%s",
                    interface.chain_i, interface.chain_j, similarity,
                    self.hyperparams.homotypic_detection_residue_similarity_threshold,
                    "PASS" if is_symmetric else "FAIL"
                )
                
                if not is_symmetric and len(residues_i) > 0 and len(residues_j) > 0:
                    # Show detailed residue information for debugging
                    self._log_detailed_residue_comparison(interface, residues_i, residues_j)
            
            return is_symmetric
            
        except Exception as e:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Failed to check residue symmetry for interface %s <-> %s: %s. Falling back to signature-only.",
                    interface.chain_i, interface.chain_j, str(e)
                )
            
            # Fallback to signature-only detection on error
            return True

    def _log_detailed_residue_comparison(self, interface: InterfaceString, 
                                    residues_i: Set[str], residues_j: Set[str]) -> None:
        """Log detailed residue comparison for debugging.
        
        Args:
            interface: Interface object.
            residues_i: Residues from chain i.
            residues_j: Residues from chain j.
        """
        try:
            # Show sample residues for debugging
            sample_i = sorted(list(residues_i))[:10]  # Show first 10
            sample_j = sorted(list(residues_j))[:10]
            
            types_i = sorted(list(self._extract_residue_types_enhanced(residues_i)))
            types_j = sorted(list(self._extract_residue_types_enhanced(residues_j)))
            
            self.workspace_manager.logger.info(
                "  Chain %s residues (%d): %s",
                interface.chain_i, len(residues_i), sample_i
            )
            self.workspace_manager.logger.info(
                "  Chain %s residues (%d): %s", 
                interface.chain_j, len(residues_j), sample_j
            )
            self.workspace_manager.logger.info(
                "  Chain %s amino acids: %s", interface.chain_i, types_i
            )
            self.workspace_manager.logger.info(
                "  Chain %s amino acids: %s", interface.chain_j, types_j
            )
            
            # Show composition comparison
            comp_i = interface.get_residue_composition_i()
            comp_j = interface.get_residue_composition_j()
            
            self.workspace_manager.logger.info(
                "  Chain %s composition: %s", interface.chain_i, comp_i
            )
            self.workspace_manager.logger.info(
                "  Chain %s composition: %s", interface.chain_j, comp_j
            )
            
            # Show sequence comparison
            seq_i = interface.get_residue_sequence_i()
            seq_j = interface.get_residue_sequence_j()
            
            self.workspace_manager.logger.info(
                "  Chain %s sequence: %s", interface.chain_i, seq_i
            )
            self.workspace_manager.logger.info(
                "  Chain %s sequence: %s", interface.chain_j, seq_j
            )
            
        except Exception as e:
            self.workspace_manager.logger.warning(
                "Failed to log detailed residue comparison: %s", str(e)
            )

    def _is_homotypic_with_residue_validation(
        self,
        interface1: InterfaceString,
        interface2: InterfaceString,
        signature1: GeometricSignature,
        signature2: GeometricSignature,
    ):
        """
        Direction-aware homotypic check.
        
        Note: for canonical sequence, we:
           - Test both.
           - If only one passes → pick that one.
           - If both pass → pick the better one (higher residue similarity min), not automatically the canonical one.

        Returns:
            (ok: bool, order: str | None)
            order is "ij" if interface2 should follow interface1's (i->f, j->b)
            order is "ji" if interface2 must be reversed
        """
        # 1) geometric checks, both directions
        direct_geom_ok = signature2.is_similar_to(
            signature1,
            distance_threshold=self.hyperparams.homodimer_distance_threshold * 10,
            angle_threshold=self.hyperparams.homodimer_angle_threshold,
        )

        # build reversed signature from interface1
        rev_sig1 = GeometricSignature(
            signature1.d_j, signature1.d_i,
            signature1.theta_j, signature1.theta_i,
        )
        flipped_geom_ok = signature2.is_similar_to(
            rev_sig1,
            distance_threshold=self.hyperparams.homodimer_distance_threshold * 10,
            angle_threshold=self.hyperparams.homodimer_angle_threshold,
        )

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "[HOMO-VALID] geometric check: direct=%s, flipped=%s (dist_th=%.3f Å, ang_th=%.3f rad)",
                direct_geom_ok,
                flipped_geom_ok,
                self.hyperparams.homodimer_distance_threshold * 10,
                self.hyperparams.homodimer_angle_threshold,
            )

        # if detection mode is purely signature-based, just pick whichever passes
        if self.hyperparams.homotypic_detection != "auto":
            if direct_geom_ok:
                return True, "ij"
            if flipped_geom_ok:
                return True, "ji"
            return False, None

        # 2) residue-based (auto) → we need residue sets for BOTH orientations
        def extract_pair(iface_a, iface_b):
            ai = self._extract_interface_residues(iface_a, iface_a.chain_i, iface_a.coord_i)
            aj = self._extract_interface_residues(iface_a, iface_a.chain_j, iface_a.coord_j)
            bi = self._extract_interface_residues(iface_b, iface_b.chain_i, iface_b.coord_i)
            bj = self._extract_interface_residues(iface_b, iface_b.chain_j, iface_b.coord_j)
            return ai, aj, bi, bj

        # interface1 is the exemplar, interface2 is the new one
        i1_i, i1_j, i2_i, i2_j = extract_pair(interface1, interface2)

        # DIRECT: (i1,i2) and (j1,j2)
        direct_res_i = self._calculate_residue_similarity(i1_i, i2_i)
        direct_res_j = self._calculate_residue_similarity(i1_j, i2_j)
        direct_min = min(direct_res_i, direct_res_j)

        # FLIPPED: (i1, j2) and (j1, i2)
        flipped_res_i = self._calculate_residue_similarity(i1_i, i2_j)
        flipped_res_j = self._calculate_residue_similarity(i1_j, i2_i)
        flipped_min = min(flipped_res_i, flipped_res_j)

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "[HOMO-VALID] residue (direct): side_i=%.3f, side_j=%.3f, min=%.3f, thr=%.3f",
                direct_res_i, direct_res_j, direct_min, self.hyperparams.homotypic_detection_residue_similarity_threshold
            )
            self.workspace_manager.logger.info(
                "[HOMO-VALID] residue (flipped): side_i=%.3f, side_j=%.3f, min=%.3f, thr=%.3f",
                flipped_res_i, flipped_res_j, flipped_min, self.hyperparams.homotypic_detection_residue_similarity_threshold
            )

        thr = self.hyperparams.homotypic_detection_residue_similarity_threshold
        direct_pass = direct_geom_ok and (direct_min >= thr)
        flipped_pass = flipped_geom_ok and (flipped_min >= thr)

        # --- DECISION ---
        if direct_pass and not flipped_pass:
            if self.workspace_manager:
                self.workspace_manager.logger.info("[HOMO-VALID] → MATCH in canonical order (ij)")
            return True, "ij"

        if flipped_pass and not direct_pass:
            if self.workspace_manager:
                self.workspace_manager.logger.info("[HOMO-VALID] → MATCH in reversed order (ji)")
            return True, "ji"

        if direct_pass and flipped_pass:
            # PICK THE BETTER ONE
            # small epsilon to avoid float noise
            eps = 1e-6
            if flipped_min > direct_min + eps:
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "[HOMO-VALID] → BOTH pass, chose FLIPPED (ji) because flipped_min=%.3f > direct_min=%.3f",
                        flipped_min, direct_min
                    )
                return True, "ji"
            else:
                if self.workspace_manager:
                    self.workspace_manager.logger.info(
                        "[HOMO-VALID] → BOTH pass, chose CANONICAL (ij) (direct_min=%.3f, flipped_min=%.3f)",
                        direct_min, flipped_min
                    )
                return True, "ij"

        # neither orientation works
        return False, None



    
    def _get_heterotypic_reason(self, template_i: str, template_j: str, 
                           signature: GeometricSignature, interface: InterfaceString) -> str:
        """Get reason why interface was classified as heterotypic.
        
        Enhanced version that provides detailed reasoning including residue analysis.
        
        Args:
            template_i: Template name for chain i.
            template_j: Template name for chain j.
            signature: Geometric signature.
            interface: Interface object.
            
        Returns:
            Human-readable reason string.
        """
        if template_i != template_j:
            return "different_molecule_types"
        
        if not signature.is_homotypic(
            self.hyperparams.homodimer_distance_threshold,
            self.hyperparams.homodimer_angle_threshold
        ):
            return "asymmetric_geometry"
        
        if self.hyperparams.homotypic_detection == "auto":
            try:
                residues_i = self._extract_interface_residues(interface, interface.chain_i, interface.coord_i)
                residues_j = self._extract_interface_residues(interface, interface.chain_j, interface.coord_j)
                similarity = self._calculate_residue_similarity(residues_i, residues_j)
                
                if similarity < self.hyperparams.homotypic_detection_residue_similarity_threshold:
                    # Get composition details for more specific reason
                    comp_i = interface.get_residue_composition_i()
                    comp_j = interface.get_residue_composition_j()
                    
                    # Find the most different amino acids
                    all_aa = set(comp_i.keys()) | set(comp_j.keys())
                    max_diff = 0
                    diff_aa = None
                    
                    for aa in all_aa:
                        count_i = comp_i.get(aa, 0)
                        count_j = comp_j.get(aa, 0)
                        diff = abs(count_i - count_j)
                        if diff > max_diff:
                            max_diff = diff
                            diff_aa = aa
                    
                    return f"low_residue_similarity_{similarity:.2f}_diff_{diff_aa}_{max_diff}"
            except Exception:
                return "residue_analysis_failed"
        
        return "unknown"
    
