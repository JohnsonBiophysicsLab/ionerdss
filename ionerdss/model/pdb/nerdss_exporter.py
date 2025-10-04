"""
ionerdss.model.pdb.nerdss_exporter

Export ionerdss System to NERDSS simulation files.

This module converts a complete ionerdss System object into the file format
required by NERDSS simulations, including .mol files for each molecule type
and parms.inp with reaction parameters. It generates `.mol` files for
each molecule type and a `parms.inp` file with reaction parameters and simulation
settings.

## Key Concepts

### Maximum Binding Sites Selection

**Problem**: Structural data often contains only partial assemblies of large
protein complexes. For example, the 6BNO structure may show only a subset of all
possible binding interfaces that exist in the complete biological assembly.

**Solution**: The exporter selects the molecule instance with the maximum number
of binding sites for each molecule type. This ensures we capture the full binding
potential of each protein, even if most instances in the structure are partially
connected.

**Mathematical Rationale**:
```
For molecule type M with instances {I₁, I₂, ..., Iₙ}:
Representative = argmax(|interfaces(Iᵢ)|)
                  i∈{1,n}
```

**Limitations**: This approach may overestimate binding capacity for systems with
mutually exclusive binding states (e.g., GATOR protein complexes where competing
interfaces prevent simultaneous maximum connectivity).

### Homotypic Interface Mapping

**Purpose**: Ensure parameter consistency for interfaces of the same structural type.

**Implementation**: Create a mapping from each site label to a representative site label:
```python
homotypic_interface_map = {
    "a1": "a1",  # Representative for type A_A_1
    "a4": "a1",  # Also type A_A_1, maps to same representative
    "a2": "a2",  # Representative for type A_A_2  
    "a3": "a2"   # Also type A_A_2, maps to same representative
}
```

## Mathematical Foundations

### Geometric Signature Calculation

NERDSS requires precise geometric parameters for each reaction. The exporter
calculates five key angles:

#### Bond Length (σ)
```
σ = ||rᵢ - rⱼ||
```
Where `rᵢ` and `rⱼ` are interface coordinate vectors.

#### Association Angles

**θ₁ (COM-Interface-Interface angle for molecule 1)**:
```
θ₁ = arccos((rᶜᵒᵐ¹ - rᵢ) · (rⱼ - rᵢ) / (||rᶜᵒᵐ¹ - rᵢ|| × ||rⱼ - rᵢ||))
```

**θ₂ (COM-Interface-Interface angle for molecule 2)**:
```
θ₂ = arccos((rᶜᵒᵐ² - rⱼ) · (rᵢ - rⱼ) / (||rᶜᵒᵐ² - rⱼ|| × ||rᵢ - rⱼ||))
```

**φ₁, φ₂ (Dihedral angles for molecular orientation)**:
```
φ₁ = dihedral(rⱼ, rᵢ, rᶜᵒᵐ¹, rⁿᵒʳᵐ¹)
φ₂ = dihedral(rᵢ, rⱼ, rᶜᵒᵐ², rⁿᵒʳᵐ²)
```

**ω (Inter-molecular dihedral angle)**:
```
ω = dihedral(rᶜᵒᵐ², rⱼ, rᵢ, rᶜᵒᵐ¹)
```

Where:
- `rᶜᵒᵐ¹`, `rᶜᵒᵐ²`: Centers of mass
- `rᵢ`, `rⱼ`: Interface coordinates  
- `rⁿᵒʳᵐ¹`, `rⁿᵒʳᵐ²`: Normal vectors (default: [0,0,1])

### Parameter Caching Mathematics

**Cache Key Generation**:
```python
representative_site1 = homotypic_map[site1]
representative_site2 = homotypic_map[site2] 
cache_key = tuple(sorted([representative_site1, representative_site2]))
```

**Cache Hit Logic**:
```
If cache_key ∈ parameter_cache:
    return cached_parameters
Else:
    calculate_parameters()
    parameter_cache[cache_key] = calculated_parameters
    return calculated_parameters
```

## Binding Site Selection Strategy

### 1. Maximum Interface Selection

**Algorithm**:
```python
def select_representative_instance(molecule_type):
    max_interfaces = 0
    representative = None
    
    for instance in system.molecule_instances:
        if instance.molecule_type == molecule_type:
            interface_count = len(instance.interfaces_neighbors_map)
            if interface_count > max_interfaces:
                max_interfaces = interface_count
                representative = instance
    
    return representative
```

**Example (6BNO case)**:
- Instance A₁: 2 interfaces (partial assembly)
- Instance A₂: 4 interfaces (more complete)
- **Selected**: A₂ (captures maximum binding potential)

### 2. Interface Type Grouping

Interfaces are grouped by their structural type:
```python
interface_groups = {
    "A_A_1": [interface1, interface4],  # Same binding type
    "A_A_2": [interface2, interface3]   # Different binding type
}
```

### 3. Site Label Generation

**Base Label Algorithm**:
```python
def get_base_site_label(mol_name, interface_type_name):
    initial = mol_name[0].lower()  # "A" → "a"
    
    # Extract index from "A_A_1" → "1"
    if "_" in interface_type_name:
        parts = interface_type_name.split("_")
        if parts[-1].isdigit():
            return f"{initial}{parts[-1]}"  # "a1"
    
    return f"{initial}1"  # Default
```

**Uniqueness Enforcement**:
```python
def get_unique_site_label(base_label, used_labels):
    if base_label not in used_labels:
        return base_label
    
    # Extract base and increment: "a1" → "a2", "a3", ...
    base = ''.join(c for c in base_label if c.isalpha())
    num = int(''.join(c for c in base_label if c.isdigit()) or "1")
    
    counter = num + 1
    while f"{base}{counter}" in used_labels:
        counter += 1
    
    return f"{base}{counter}"
```

## Reaction Generation Logic

### Homotypic Reactions (A + A)

**Mathematical Basis**: For n sites of the same type, generate all unique
pairwise combinations without double-counting.

**Formula**: Total reactions = n(n+1)/2

**Algorithm**:
```python
for i, site1 in enumerate(all_sites):
    for j, site2 in enumerate(all_sites):
        if i <= j:  # Avoid duplicates: (i,j) yes, (j,i) no
            generate_reaction(site1, site2)
            
            # Cross-reactions (i≠j) get doubled rate
            is_cross_reaction = (i != j)
```

**Example**: Sites [a1, a2, a4] (where a1≡a4 structurally)
- Generated: (a1,a1), (a1,a2), (a1,a4), (a2,a2), (a2,a4), (a4,a4)
- Cross-reactions: (a1,a2), (a1,a4), (a2,a4) get 2× rate
- Self-reactions: (a1,a1), (a2,a2), (a4,a4) get 1× rate

**Rate Doubling Rationale**: Cross-reactions like A(a1)+A(a2) can occur in two
orientations in solution, effectively doubling the encounter rate.

### Heterotypic Reactions (A + B)

**Mathematical Basis**: For n sites on molecule A and m sites on molecule B,
generate all n×m combinations.

**Formula**: Total reactions = n × m

**Algorithm**:
```python
for site_a in molecule_a_sites:
    for site_b in molecule_b_sites:
        generate_reaction(site_a, site_b)
        # No rate doubling needed
```

**Example**: A has sites [a1, a2], B has sites [b1, b2]
- Generated: A(a1)+B(b1), A(a1)+B(b2), A(a2)+B(b1), A(a2)+B(b2)
- All reactions get same base rate (no doubling)

**No Rate Doubling**: Each A-B pair represents a unique molecular encounter.

## Parameter Calculation

### Coordinate System Setup

**Default Normal Vectors**: [0, 0, 1] (z-axis)
```python
norm1 = np.array([0.0, 0.0, 1.0])
norm2 = np.array([0.0, 0.0, 1.0])

# Convert to absolute coordinates
abs_norm1 = com1 + norm1
abs_norm2 = com2 + norm2
```

### Angle Calculation Pipeline

**1. Theta Angles (Bond angles)**:
```python
theta1 = angles_from_points(com1, intf1, intf2)
theta2 = angles_from_points(com2, intf2, intf1)
```

**2. Phi Angles (Orientation dihedrals)**:
```python
phi1 = dihedrals_from_points(intf2, intf1, com1, abs_norm1)
phi2 = dihedrals_from_points(intf1, intf2, com2, abs_norm2)
```

**3. Omega Angle (Inter-molecular dihedral)**:
```python
omega = dihedrals_from_points(com2, intf2, intf1, com1)
```

### Parameter Caching Strategy

**Cache Key Construction**:
```python
# Map sites to representatives
rep1 = homotypic_interface_map[site1]  # a4 → a1
rep2 = homotypic_interface_map[site2]  # a4 → a1

# Create canonical cache key
cache_key = tuple(sorted([rep1, rep2]))  # (a1, a1)
```

**Cache Benefits**:
- Ensures A(a1)+A(a1) and A(a4)+A(a4) have identical parameters
- Reduces computation for repeated interface types
- Maintains consistency across similar binding events

## File Format Specifications

### MOL File Format

```
Name = ProteinA

# translational diffusion constants
D = [1.000000, 1.000000, 1.000000]

# rotational diffusion constants  
Dr = [0.100000, 0.100000, 0.100000]

# Coordinates
COM   0.0000000   0.0000000   0.0000000
a1    1.0000000   0.0000000   0.0000000
a2    0.0000000   1.0000000   0.0000000

# bonds
bonds = 2
com a1
com a2
```

**Key Elements**:
- **Name**: Molecule type identifier
- **D**: Translational diffusion constants (nm²/μs)
- **Dr**: Rotational diffusion constants (rad²/μs)  
- **Coordinates**: COM at origin, interfaces at local coordinates
- **Bonds**: All interfaces bonded to COM (rigid body assumption)

### PARMS File Format

```
start parameters
    nItr = 100000.0
    timestep = 0.5
    timeWrite = 10000.0
    trajWrite = 100000.0
    restartWrite = 100000.0
    checkPoint = 100000.0
    pdbWrite = 100000.0
end parameters

start boundaries
    WaterBox = [100.0, 100.0, 100.0] #nm
end boundaries

start molecules
    ProteinA : 10
    ProteinB : 5
end molecules

start reactions
    
    # Binding reactions
    ProteinA(a1) + ProteinA(a1) <-> ProteinA(a1!1).ProteinA(a1!1)
    onRate3Dka = 100.0
    offRatekb = 1000.028
    norm1 = [0.0, 0.0, 1.0]
    norm2 = [0.0, 0.0, 1.0]
    sigma = 1.2345678
    assocAngles = [1.234,0.567,-0.891,0.234,-0.456]
    
    ProteinA(a1) + ProteinA(a2) <-> ProteinA(a1!1).ProteinA(a2!1)
    onRate3Dka = 200.0  # Doubled for cross-reaction
    offRatekb = 1000.028
    norm1 = [0.0, 0.0, 1.0]
    norm2 = [0.0, 0.0, 1.0]
    sigma = 1.5678901
    assocAngles = [0.789,1.234,-0.567,0.891,-0.234]
    
end reactions
```

**Rate Doubling Logic**:
- **Self-reactions** (a1+a1): Base rate
- **Cross-reactions** (a1+a2): 2× base rate
- **Heterotypic reactions** (A+B): Base rate (no doubling)

## Usage Examples

### Basic Export

```python
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter

# Create exporter
exporter = NERDSSExporter(system, workspace_manager)

# Export with defaults
output_files = exporter.export_all()

# Files created:
# - ProteinA.mol
# - ProteinB.mol  
# - parms.inp
```

### Custom Parameters

```python
# Custom molecule counts and simulation box
output_files = exporter.export_all(
    molecule_counts={"ProteinA": 50, "ProteinB": 25},
    box_nm=(200.0, 200.0, 200.0),
    parms_overrides={
        "nItr": 5e5,
        "timestep": 0.1,
        "onRate3Dka": 500.0,
        "offRatekb": 2000.0
    }
)
```

### Parameter Analysis

```python
# Access internal mappings for analysis
print("Interface to site mapping:")
for interface, site in exporter.interface_to_site_map.items():
    print(f"  {interface} → {site}")

print("Homotypic mapping:")
for site, representative in exporter.homotypic_interface_map.items():
    print(f"  {site} → {representative}")

print("Reaction metadata:")
for i, metadata in enumerate(exporter.reaction_metadata):
    print(f"  Reaction {i}: cross_reaction={metadata['is_cross_reaction']}")
```

"""

import re
from typing import Dict, Any, Iterable, Optional, List, Tuple
from pathlib import Path

import numpy as np

from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.utils.angles import angles_from_points, dihedrals_from_points
from .file_manager import WorkspaceManager


class NERDSSExporter:
    """Exporter for converting ionerdss System to NERDSS simulation files.

    Generates .mol files for each molecule type and parms.inp with simulation
    parameters and reaction definitions.

    Attributes:
        system: ionerdss System to export.
        workspace_manager: Workspace manager for file organization.
        output_dir: Directory for NERDSS output files.
        interface_to_site_map: Mapping from interface names to site labels.
        reaction_metadata: Metadata about reactions for rate calculation.
    """

    def __init__(self, system: System, workspace_manager: Optional[WorkspaceManager] = None):
        """Initialize NERDSS exporter."""
        self.system = system
        self.workspace_manager = workspace_manager

        # Store mapping from interface names to site labels
        self.interface_to_site_map: Dict[str, str] = {}

        # Store reaction metadata for rate calculation
        self.reaction_metadata: List[Dict[str, Any]] = []

        # Cache for homotypic interface mapping: site_label -> representative_site_label
        self.homotypic_interface_map: Dict[str, str] = {}

        # Cache for reaction parameters: (site1, site2) -> (sigma, angles)
        self.reaction_params_cache: Dict[Tuple[str, str],
                                         Tuple[float, Tuple[float, float, float, float, float]]] = {}

        # Create NERDSS output directory in workspace
        if workspace_manager:
            self.output_dir = workspace_manager.workspace_path / 'nerdss_files'
            self.output_dir.mkdir(exist_ok=True)
            workspace_manager.logger.info(
                "Created NERDSS export directory: %s", self.output_dir)
        else:
            self.output_dir = Path("nerdss_files")
            self.output_dir.mkdir(exist_ok=True)

    def export_all(self, molecule_counts: Optional[Dict[str, int]] = None,
                   box_nm: Tuple[float, float, float] = (100.0, 100.0, 100.0),
                   parms_overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Path]:
        """Export complete NERDSS simulation setup.

        Args:
            molecule_counts: Number of molecules per type. Defaults to 10 each.
            box_nm: Simulation box size in nm. Default (100, 100, 100).
            parms_overrides: Additional parameters for parms.inp.

        Returns:
            Dictionary mapping file types to output paths.
        """
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Exporting NERDSS simulation files...")

        output_files = {}

        # Set default molecule counts
        if molecule_counts is None:
            molecule_counts = {}
            for mol_type in self.system.molecule_types:
                molecule_counts[mol_type.name] = 10

        # Clear mappings for fresh export
        self.interface_to_site_map.clear()
        self.reaction_metadata.clear()

        # Export .mol files for each molecule type (this builds the mapping)
        for mol_type in self.system.molecule_types:
            mol_file_path = self._write_mol_file(mol_type)
            output_files[f"{mol_type.name}_mol"] = mol_file_path

        # Generate reactions using the stored mapping
        reactions = self._generate_reactions()

        # Calculate bond lengths and angles for reactions
        sigma_list, angles_list = self._calculate_reaction_parameters(
            reactions)

        # Export parms.inp
        parms_path = self._write_parms_file(
            reactions=reactions,
            molecule_counts=molecule_counts,
            box_nm=box_nm,
            sigma_list=sigma_list,
            angles_list=angles_list,
            parms_overrides=parms_overrides
        )
        output_files['parms'] = parms_path

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Exported %d NERDSS files", len(output_files))
            for file_type, file_path in output_files.items():
                self.workspace_manager.logger.info(
                    "  %s: %s", file_type, file_path)

        return output_files

    def _write_mol_file(self, mol_type: MoleculeType) -> Path:
        """Write .mol file for a molecule type and build interface-to-site mapping."""
        mol_file_path = self.output_dir / f"{mol_type.name}.mol"

        # Find a representative molecule instance of this type to get actual interface counts
        representative_instance = None
        max_interfaces = 0

        for mol_instance in self.system.molecule_instances:
            if mol_instance.molecule_type and mol_instance.molecule_type.name == mol_type.name:
                interface_count = len(mol_instance.interfaces_neighbors_map)
                if interface_count > max_interfaces:
                    max_interfaces = interface_count
                    representative_instance = mol_instance

        if not representative_instance:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "No instances found for molecule type %s", mol_type.name
                )
            return mol_file_path

        # Group interface instances by their interface type
        interface_groups = {}
        for interface_instance, partner_instance\
            in representative_instance.interfaces_neighbors_map.items():
            if interface_instance.interface_type:
                type_name = interface_instance.interface_type.get_name()

                if type_name not in interface_groups:
                    interface_groups[type_name] = []

                interface_groups[type_name].append({
                    'instance': interface_instance,
                    'coord': interface_instance.interface_type.local_coord,
                    'partner': partner_instance.molecule_type.name\
                        if partner_instance.molecule_type else "unknown",
                    'type_name': type_name
                })

        # Generate unique site labels and build homotypic mapping
        interfaces = []
        used_labels = set()
        # Maps interface type to its first (representative) site label
        type_to_representative = {}

        for type_name, type_interfaces in interface_groups.items():
            # Sort interfaces by coordinate to ensure deterministic ordering
            type_interfaces.sort(key=lambda x: tuple(x['coord']))

            for i, interface_data in enumerate(type_interfaces):
                # Generate unique site label
                base_label = self._get_base_site_label(
                    mol_type.name, type_name)
                site_label = self._get_unique_site_label_with_base(
                    base_label, used_labels)
                used_labels.add(site_label)

                # Store the mapping for later use
                interface_key = f"{type_name}_{i+1}"
                self.interface_to_site_map[interface_key] = site_label

                # Build homotypic interface mapping
                if type_name not in type_to_representative:
                    # First interface of this type becomes the representative
                    type_to_representative[type_name] = site_label
                    self.homotypic_interface_map[site_label] = site_label
                else:
                    # Map this interface to the representative
                    representative_site = type_to_representative[type_name]
                    self.homotypic_interface_map[site_label] = representative_site

                # Also map the type name for the first instance (for backward compatibility)
                if i == 0:
                    self.interface_to_site_map[type_name] = site_label

                interfaces.append({
                    'type_name': type_name,
                    'instance_key': interface_key,
                    'site_label': site_label,
                    'coord': interface_data['coord'],
                    'partner': interface_data['partner']
                })

                if self.workspace_manager:
                    self.workspace_manager.logger.debug(
                        "Mapped interface %s -> site %s -> representative %s",
                        interface_key, site_label, self.homotypic_interface_map[site_label])

        # Write .mol file
        with open(mol_file_path, 'w', encoding='utf-8') as f:
            f.write(f"Name = {mol_type.name}\n\n")

            # Diffusion constants
            D_t = mol_type.D_t_nm2_us
            D_r = mol_type.D_r_rad2_us

            f.write("# translational diffusion constants\n")
            f.write(f"D = [{D_t:.6g}, {D_t:.6g}, {D_t:.6g}]\n\n")
            f.write("# rotational diffusion constants\n")
            f.write(f"Dr = [{D_r:.6g}, {D_r:.6g}, {D_r:.6g}]\n\n")

            # Coordinates
            f.write("# Coordinates\n")
            f.write(f"COM   {self._format_vec([0.0, 0.0, 0.0])}\n")

            # Interface sites with unique labels
            for interface in interfaces:
                site_label = interface['site_label']
                coord = interface['coord']
                f.write(f"{site_label}  {self._format_vec(coord)}\n")

            # Bonds from COM to each site
            f.write("\n# bonds\n")
            f.write(f"bonds = {len(interfaces)}\n")
            for interface in interfaces:
                site_label = interface['site_label']
                f.write(f"com {site_label}\n")

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Wrote .mol file with %d interfaces: %s", len(interfaces), mol_file_path)
            self.workspace_manager.logger.info(
                "Homotypic interface mapping: %s", self.homotypic_interface_map)

        return mol_file_path

    def _get_base_site_label(self, mol_name: str, interface_type_name: str) -> str:
        """Get base site label from molecule name and interface type.

        Args:
            mol_name: Molecule name.
            interface_type_name: Interface type name.

        Returns:
            Base site label.
        """
        initial = mol_name[0].lower() if mol_name else "x"

        # Extract index from interface type name
        if "_" in interface_type_name:
            parts = interface_type_name.split("_")
            if len(parts) >= 3 and parts[-1].isdigit():
                return f"{initial}{parts[-1]}"

        return f"{initial}1"

    def _get_unique_site_label_with_base(self, base_label: str, used_labels: set) -> str:
        """Generate unique site label starting from base label.

        Args:
            base_label: Base label to start from.
            used_labels: Set of already used labels.

        Returns:
            Unique site label.
        """
        if base_label not in used_labels:
            return base_label

        # Extract base and number
        base = ''.join(c for c in base_label if c.isalpha())
        base_num = int(''.join(c for c in base_label if c.isdigit()) or "1")

        # Find unique label
        counter = base_num + 1
        while True:
            label = f"{base}{counter}"
            if label not in used_labels:
                return label
            counter += 1

    def _generate_reactions(self) -> List[str]:
        """Generate BNGL reaction strings with proper handling of multiple homotypic interfaces."""
        reactions = []
        processed_pairs = set()

        # Group interface types by molecule pair and index
        interface_pairs = {}
        for interface_type in self.system.interface_types:
            mol1 = interface_type.this_mol_type_name
            mol2 = interface_type.partner_mol_type_name
            index = interface_type.interface_index

            # For heterotypic, don't canonicalize - process both directions
            pair_key = (mol1, mol2, index)

            if pair_key not in interface_pairs:
                interface_pairs[pair_key] = []
            interface_pairs[pair_key].append(interface_type)

        for pair_key, interface_types in interface_pairs.items():
            mol1, mol2, index = pair_key

            # Create canonical key only for duplicate checking
            canonical_key = tuple(sorted([mol1, mol2]) + [index])

            if canonical_key in processed_pairs:
                continue
            processed_pairs.add(canonical_key)

            # Get sites for both molecule types
            type_name = interface_types[0].get_name()
            mol1_sites = []
            mol2_sites = []

            # Find sites for mol1
            for key, site_label in self.interface_to_site_map.items():
                if key.startswith(type_name + "_") or key == type_name:
                    if site_label not in mol1_sites:
                        mol1_sites.append(site_label)

            # Find sites for mol2
            partner_type_name = f"{mol2}_{mol1}_{index}"
            for key, site_label in self.interface_to_site_map.items():
                if key.startswith(partner_type_name + "_") or key == partner_type_name:
                    if site_label not in mol2_sites:
                        mol2_sites.append(site_label)

            if mol1 == mol2:  # Homotypic case
                all_sites = sorted(set(mol1_sites + mol2_sites))
                # Generate reactions: (i,j) where i <= j to avoid duplicates
                for i, site1 in enumerate(all_sites):
                    for j, site2 in enumerate(all_sites):
                        if i <= j:  # ✅ Correct: avoids (j,i) when (i,j) exists
                            reaction = f"{mol1}({site1}) + {mol2}({site2}) <-> {mol1}({site1}!1).{mol2}({site2}!1)"
                            reactions.append(reaction)

                            is_cross_reaction = (i != j)
                            self.reaction_metadata.append({
                                'reaction': reaction,
                                'is_cross_reaction': is_cross_reaction,
                                'mol1': mol1, 'mol2': mol2,
                                'site1': site1, 'site2': site2
                            })
            else:  # Heterotypic case
                # ✅ Generate ALL combinations: A(i) + B(j) for all i,j
                for site1 in mol1_sites:
                    for site2 in mol2_sites:
                        reaction = f"{mol1}({site1}) + {mol2}({site2}) <-> {mol1}({site1}!1).{mol2}({site2}!1)"
                        reactions.append(reaction)

                        self.reaction_metadata.append({
                            'reaction': reaction,
                            'is_cross_reaction': False,
                            'mol1': mol1, 'mol2': mol2,
                            'site1': site1, 'site2': site2
                        })

        return reactions

    def _calculate_reaction_parameters(self, reactions: List[str]) -> Tuple[List[float], List[Tuple[float, float, float, float, float]]]:
        """Calculate bond lengths and angles for reactions using homotypic mapping cache."""
        sigma_list = []
        angles_list = []

        # Regex to parse reaction format
        reaction_re = re.compile(
            r"^\s*([A-Za-z0-9_]+)\(([A-Za-z0-9_]+)\)\s*\+\s*([A-Za-z0-9_]+)\(([A-Za-z0-9_]+)\)"
        )

        for reaction in reactions:
            match = reaction_re.match(reaction)
            if not match:
                sigma_list.append(1.0)
                angles_list.append((0.0, 0.0, 0.0, 0.0, 0.0))
                continue

            mol1, site1, mol2, site2 = match.groups()

            # Map sites to their representatives using homotypic mapping
            representative_site1 = self.homotypic_interface_map.get(
                site1, site1)
            representative_site2 = self.homotypic_interface_map.get(
                site2, site2)

            # Create cache key using representative sites (canonical order)
            cache_key = tuple(
                sorted([representative_site1, representative_site2]))

            if self.workspace_manager:
                self.workspace_manager.logger.debug(
                    "Reaction %s(%s) + %s(%s) -> representatives %s + %s -> cache key %s",
                    mol1, site1, mol2, site2, representative_site1, representative_site2, cache_key
                )

            # Check cache first
            if cache_key in self.reaction_params_cache:
                sigma, angles = self.reaction_params_cache[cache_key]
                sigma_list.append(sigma)
                angles_list.append(angles)

                if self.workspace_manager:
                    self.workspace_manager.logger.debug(
                        "Using cached parameters for %s: sigma=%.3f", cache_key, sigma
                    )
                continue

            # Calculate parameters for the first time
            com1, com2, interface1_coord, interface2_coord = self._get_coms_interfaces(
                mol1, representative_site1, mol2, representative_site2)

            if interface1_coord is None or interface2_coord is None or com1 is None or com2 is None:
                sigma, angles = 1.0, (0.0, 0.0, 0.0, 0.0, 0.0)
            else:
                sigma, angles = self._generate_reaction_angles(
                    interface1_coord, interface2_coord, com1, com2
                )

            # Cache the result
            self.reaction_params_cache[cache_key] = (sigma, angles)

            sigma_list.append(sigma)
            angles_list.append(angles)

            if self.workspace_manager:
                self.workspace_manager.logger.debug(
                    "Calculated and cached parameters for %s: sigma=%.3f", cache_key, sigma
                )

        return sigma_list, angles_list

    def _get_coms_interfaces(self, mol1_name: str, site1: str, mol2_name: str, site2: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
        """Get center of mass and interface site coordinates for two molecules.

        Args:
            mol1_name: First molecule type name.
            site1: First site label.
            mol2_name: Second molecule type name.
            site2: Second site label.

        Returns:
            Tuple of (com1, com2, interface1_coord, interface2_coord).
        """
        # Convert site labels back to interface names using stored mapping
        interface1_name = None
        interface2_name = None

        # Find interface names from site labels using reverse mapping
        for intf_name, site_label in self.interface_to_site_map.items():
            if site_label == site1:
                # Check if this interface belongs to mol1
                if intf_name.startswith(mol1_name + "_"):
                    interface1_name = intf_name
            if site_label == site2:
                # Check if this interface belongs to mol2
                if intf_name.startswith(mol2_name + "_"):
                    interface2_name = intf_name

        if not interface1_name or not interface2_name:
            if self.workspace_manager:
                self.workspace_manager.logger.warning(
                    "Could not find interface names for sites %s:%s or %s:%s",
                    mol1_name, site1, mol2_name, site2
                )
            return None, None, None, None

        # Initialize
        coord1 = None
        coord2 = None
        com1 = None
        com2 = None

        # Look up molecule instances to get their COMs
        for mol_instance in self.system.molecule_instances:
            if mol_instance.molecule_type and mol_instance.molecule_type.name == mol1_name:
                com1 = mol_instance.com
                for this_interface, neighbor_instance in mol_instance.interfaces_neighbors_map.items():
                    if (neighbor_instance.molecule_type and
                        neighbor_instance.molecule_type.name == mol2_name and
                            this_interface.interface_type.get_name() == interface1_name):
                        com2 = neighbor_instance.com
                        coord1 = this_interface.absolute_coord
                        coord2 = this_interface.partner_interface.absolute_coord
                        break

            # Break early if we found both
            if com1 is not None and com2 is not None and coord1 is not None and coord2 is not None:
                break

        return com1, com2, coord1, coord2

    def _generate_reaction_angles(self, intf1: np.ndarray, intf2: np.ndarray,
                                  com1: np.ndarray, com2: np.ndarray) -> Tuple[float, Tuple[float, float, float, float, float]]:
        """Generate reaction angles and bond length.

        Args:
            intf1: Interface 1 coordinates.
            intf2: Interface 2 coordinates.
            com1: COM 1 coordinates.
            com2: COM 2 coordinates.

        Returns:
            Tuple of (bond_length, (theta1, theta2, phi1, phi2, omega)).
        """
        # Default normal vectors
        norm1 = np.array([0.0, 0.0, 1.0])
        norm2 = np.array([0.0, 0.0, 1.0])

        # Convert to absolute normal vectors
        abs_norm1 = com1 + norm1
        abs_norm2 = com2 + norm2

        # Calculate angles
        theta1 = angles_from_points(com1, intf1, intf2)
        theta2 = angles_from_points(com2, intf2, intf1)
        phi1 = dihedrals_from_points(intf2, intf1, com1, abs_norm1)
        phi2 = dihedrals_from_points(intf1, intf2, com2, abs_norm2)
        omega = dihedrals_from_points(com2, intf2, intf1, com1)

        # Bond length
        bond_length = np.linalg.norm(intf1 - intf2)

        return bond_length, (theta1, theta2, phi1, phi2, omega)

    def _write_parms_file(self, reactions: List[str], molecule_counts: Dict[str, int],
                          box_nm: Tuple[float, float, float], sigma_list: List[float],
                          angles_list: List[Tuple[float, float, float, float, float]],
                          parms_overrides: Optional[Dict[str, Any]] = None) -> Path:
        """Write parms.inp file with proper handling of homotypic cross-reactions.

        Args:
            reactions: List of reaction strings.
            molecule_counts: Number of molecules per type.
            box_nm: Simulation box size.
            sigma_list: Bond lengths for each reaction.
            angles_list: Angles for each reaction.
            parms_overrides: Additional parameters.

        Returns:
            Path to written parms.inp file.
        """
        parms_path = self.output_dir / "parms.inp"

        # Default parameters
        params = {
            'nItr': 1e5,
            'timestep': 0.5,
            'timeWrite': 1e4,
            'trajWrite': 1e5,
            'restartWrite': 1e5,
            'checkPoint': 1e5,
            'pdbWrite': 1e5,
            'onRate3Dka': 100.0,
            'offRatekb': 1000.028,
            'norm1': (0.0, 0.0, 1.0),
            'norm2': (0.0, 0.0, 1.0),
        }

        # Apply overrides
        if parms_overrides:
            params.update(parms_overrides)

        with open(parms_path, 'w', encoding='utf-8') as f:
            # Parameters section
            f.write("start parameters\n")
            f.write(f"    nItr = {params['nItr']}\n")
            f.write(f"    timestep = {params['timestep']}\n")
            f.write(f"    timeWrite = {params['timeWrite']}\n")
            f.write(f"    trajWrite = {params['trajWrite']}\n")
            f.write(f"    restartWrite = {params['restartWrite']}\n")
            f.write(f"    checkPoint = {params['checkPoint']}\n")
            f.write(f"    pdbWrite = {params['pdbWrite']}\n")
            f.write("end parameters\n\n")

            # Boundaries section
            f.write("start boundaries\n")
            f.write(
                f"    WaterBox = [{box_nm[0]}, {box_nm[1]}, {box_nm[2]}] #nm\n")
            f.write("end boundaries\n\n")

            # Molecules section
            f.write("start molecules\n")
            for mol_name, count in molecule_counts.items():
                f.write(f"    {mol_name} : {count}\n")
            f.write("end molecules\n\n")

            # Reactions section
            f.write("start reactions\n")
            f.write("    \n")
            f.write("    # Binding reactions\n")

            for i, reaction in enumerate(reactions):
                f.write(f"    {reaction}\n")

                # Determine onRate3Dka based on whether it's a cross-reaction
                base_on_rate = params['onRate3Dka']
                if i < len(self.reaction_metadata):
                    if self.reaction_metadata[i]['is_cross_reaction']:
                        on_rate = base_on_rate * 2.0  # Double for cross-reactions
                    else:
                        on_rate = base_on_rate
                else:
                    on_rate = base_on_rate

                f.write(f"    onRate3Dka = {on_rate}\n")
                f.write(f"    offRatekb = {params['offRatekb']}\n")

                norm1 = params['norm1']
                norm2 = params['norm2']
                f.write(f"    norm1 = [{norm1[0]}, {norm1[1]}, {norm1[2]}]\n")
                f.write(f"    norm2 = [{norm2[0]}, {norm2[1]}, {norm2[2]}]\n")

                # Use calculated sigma or default
                sigma = sigma_list[i] if i < len(sigma_list) else 1.0
                f.write(f"    sigma = {sigma}\n")

                # Use calculated angles or default
                angles = angles_list[i] if i < len(
                    angles_list) else (0.0, 0.0, 0.0, 0.0, 0.0)
                f.write(
                    f"    assocAngles = [{angles[0]},{angles[1]},{angles[2]},{angles[3]},{angles[4]}]\n")
                f.write("    \n")

            f.write("end reactions\n")

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Wrote parms.inp file: %s", parms_path)

        return parms_path

    def _format_vec(self, vec: Iterable[float], precision: int = 7) -> str:
        """Format vector for output files.

        Args:
            vec: Vector to format.
            precision: Number of decimal places.

        Returns:
            Formatted vector string.
        """
        return "   ".join(f"{float(v):.{precision}f}" for v in vec)
