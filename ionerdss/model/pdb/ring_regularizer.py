"""
ionerdss.model.pdb.ring_regularizer

Ring structure regularization for homo-n-mer assemblies.

This module detects and regularizes ring-like homo-n-mer structures to minimize
geometric deviations that can accumulate and cause artifacts in NERDSS simulations.
Small crystallographic deviations can be geometrically amplified when assembled
in rings, so this module aligns all subunits to ideal geometric configurations.

REGULARIZATION MODES:
- "off": No regularization (default)
- "separate": Each detected ring is fitted to its own ideal geometry
- "uniform": All rings are fitted to a single best-fit geometry

GEOMETRIC CONFIGURATIONS:
- "cylinder": Fit subunits to cylindrical surface with uniform z-plane
- "sphere": Fit subunits to spherical surface

DETECTION CRITERIA:
1. Homo-n-mer ring with n >= 3 subunits
2. Ring connectivity via same interface types (homotypic A_A_1 or heterotypic A_A_1<->A_A_2)
3. Closed ring topology: within the specific homo-n-mer ring, each subunit connects 
   to exactly 2 neighbors using the ring-forming interface types. Subunits may have 
   additional connections to other molecules (homomers or heteromers) outside the ring.

DETAILED METHODS AND USAGE

## Key Concepts

### Ring Structure Problem

**Issue**: Crystallographic structures often contain small positional errors
(±0.1-0.5 Å) that are negligible for individual proteins but become significant
when propagated around ring structures. In a 6-mer ring, a 0.2 Å deviation per
subunit can accumulate to 1.2 Å closure error, causing simulation instabilities.

**Solution**: Detect ring structures and regularize them to ideal geometric
configurations (cylinders or spheres) while preserving relative orientations
and interface relationships.

### Homo-n-mer Rings

**Definition**: Closed ring structures formed by n identical subunits
(homo-n-mers) where each subunit connects to exactly 2 neighbors using
the same interface types.

**Examples**:
- **Homo-trimer**: 3 identical subunits in triangular arrangement
- **Homo-hexamer**: 6 identical subunits in hexagonal arrangement  
- **Homo-dodecamer**: 12 identical subunits in large ring

### Ring Topology Validation

**Closed Ring Requirement**: Within the ring, each subunit must connect
to exactly 2 neighbors using ring-forming interface types. Additional
connections outside the ring are allowed.

**Interface Consistency**: All ring connections must use either:
- **Homotypic interfaces**: Same interface type (e.g., all A_A_1)
- **Complementary heterotypic**: Alternating complementary types
(e.g., A_B_1 ↔ B_A_1)

## Detection Criteria

### 1. Minimum Ring Size
```python
min_ring_size >= 3  # Default minimum (triangle)
```
Rings smaller than 3 subunits are not geometrically meaningful.

### 2. Homo-n-mer Requirement
All subunits in the ring must be of the same molecule type:
```python
all(subunit.molecule_type.name == reference_type for subunit in ring)
```

### 3. Closed Ring Topology
Each subunit connects to exactly 2 ring neighbors:
```python
for subunit in ring:
    ring_connections = count_ring_interface_connections(subunit)
    assert ring_connections == 2  # Previous and next in ring
```

### 4. Interface Consistency
Ring interfaces must be structurally compatible:
```python
# Valid: All same type
interface_types = {"A_A_1"}

# Valid: Complementary pair  
interface_types = {"A_B_1", "B_A_1"}

# Invalid: Multiple incompatible types
interface_types = {"A_A_1", "B_B_1", "C_C_1"}
```

## Regularization Modes

### "off" Mode (Default)
```python
regularizer = RingRegularizer(system, mode="off")
result = regularizer.regularize()  # Returns False, no processing
```

**Use Case**: When ring structures are already well-formed or regularization is not desired.

### "separate" Mode
```python
regularizer = RingRegularizer(system, mode="separate", geometry="cylinder")
```

**Behavior**: Each detected ring is fitted to its own ideal geometry independently.

**Use Case**: When rings have different sizes or conformations that should be preserved.

### "uniform" Mode
```python
regularizer = RingRegularizer(system, mode="uniform", geometry="sphere")
```

**Behavior**: All rings are fitted to a single best-fit geometry calculated from all ring positions.

**Use Case**: Viral capsids or other assemblies with multiple identical ring structures.

## Geometric Configurations

### Cylinder Geometry
```python
regularizer = RingRegularizer(system, geometry="cylinder")
```

**Mathematical Model**:
- **Center**: (cx, cy, cz)
- **Axis**: Unit vector (ax, ay, az)  
- **Radius**: Distance from axis to subunit centers

**Positioning Algorithm**:
```python
angles = np.linspace(0, 2π, n_subunits, endpoint=False)
for i, angle in enumerate(angles):
    # Position on unit cylinder
    local_pos = radius * [cos(angle), sin(angle), 0]
    
    # Rotate to align with cylinder axis
    rotated_pos = rotation_matrix @ local_pos
    
    # Translate to cylinder center
    new_position = center + rotated_pos
```

**Use Cases**:
- Membrane protein complexes
- Ring structures with preferred axis orientation
- Structures where subunits should lie in parallel planes

### Sphere Geometry
```python
regularizer = RingRegularizer(system, geometry="sphere")
```

**Mathematical Model**:
- **Center**: (cx, cy, cz)
- **Radius**: Distance from center to subunit centers
- **No preferred axis**: Spherically symmetric

**Positioning Algorithm**:
```python
angles = np.linspace(0, 2π, n_subunits, endpoint=False)
for i, angle in enumerate(angles):
    # Position on unit circle in ring plane
    local_pos = [cos(angle), sin(angle), 0]
    
    # Rotate to preserve ring plane orientation
    rotated_pos = ring_plane_rotation @ local_pos
    
    # Scale to sphere radius and translate
    new_position = center + radius * rotated_pos
```

**Use Cases**:
- Globular protein assemblies
- Structures without preferred orientation
- Ring structures on curved surfaces

## Usage Examples

### Basic Ring Detection and Regularization

```python
from ionerdss.model.pdb.ring_regularizer import RingRegularizer

# Create regularizer with default settings
regularizer = RingRegularizer(
    system=system,
    workspace_manager=workspace_manager,
    mode="separate",
    geometry="cylinder",
    min_ring_size=3
)

# Perform regularization
success = regularizer.regularize()

if success:
    print("Ring regularization completed")
    summary = regularizer.get_summary()
    print(f"Detected {summary['rings_detected']} rings")
    print(f"Average fit error: {summary['average_fit_error']:.3f} Å")
else:
    print("No rings detected or regularization disabled")
```

### Analyzing Detection Results

```python
# Get detailed summary
summary = regularizer.get_summary()

print(f"Mode: {summary['mode']}")
print(f"Geometry: {summary['geometry']}")
print(f"Rings detected: {summary['rings_detected']}")
print(f"Ring sizes: {summary['ring_sizes']}")
print(f"Fit errors: {summary['fit_errors']}")

# Access detected ring structures
for i, ring in enumerate(regularizer.detected_rings):
    print(f"Ring {i}:")
    print(f"  Subunits: {len(ring.molecules)}")
    print(f"  Interface types: {ring.interface_types}")
    print(f"  Center: {ring.ring_center}")
    print(f"  Radius: {ring.ring_radius:.2f} Å")
```

### Configuration for Different Systems

**Small Symmetric Assemblies**:
```python
# Detect small rings, use cylinder geometry
regularizer = RingRegularizer(
    system=system,
    mode="separate",
    geometry="cylinder", 
    min_ring_size=3
)
```

**Large Viral Capsids**:
```python
# Use uniform fitting for consistency across multiple rings
regularizer = RingRegularizer(
    system=system,
    mode="uniform",
    geometry="sphere",
    min_ring_size=5  # Larger minimum for capsid subunits
)
```

**Membrane Protein Complexes**:
```python
# Cylinder geometry for membrane-embedded rings
regularizer = RingRegularizer(
    system=system,
    mode="separate", 
    geometry="cylinder",
    min_ring_size=4
)
```

### Integration with Pipeline

```python
from ionerdss.model.pdb.parser import PDBParser
from ionerdss.model.pdb.coarse_graining import CoarseGrainer
from ionerdss.model.pdb.ring_regularizer import RingRegularizer

# Complete pipeline with ring regularization
parser = PDBParser("viral_capsid.pdb")
coarse_grainer = CoarseGrainer(parser, hyperparams)

# Create system
system = build_system(parser, coarse_grainer)

# Apply ring regularization before final processing
regularizer = RingRegularizer(
    system=system,
    workspace_manager=workspace_manager,
    mode="uniform",
    geometry="sphere"
)

if regularizer.regularize():
    print("Applied ring regularization")
    
# Continue with template building and export
template_builder = TemplateBuilder(system, ...)
```

## Algorithm Details

### Ring Detection Pipeline

**1. Connectivity Graph Construction**:
```python
graph = nx.DiGraph()
for molecule in system.molecule_instances:
    graph.add_node(molecule.name, molecule=molecule)
    
    for interface, partner in molecule.interfaces_neighbors_map.items():
        graph.add_edge(
            molecule.name, 
            partner.name,
            interface_type=interface.interface_type.get_name()
        )
```

**2. Cycle Detection**:
```python
cycles = list(nx.simple_cycles(graph))
valid_cycles = [cycle for cycle in cycles if len(cycle) >= min_ring_size]
```

**3. Ring Validation**:
```python
for cycle in valid_cycles:
    # Check homo-n-mer requirement
    molecule_types = {get_molecule_type(mol) for mol in cycle}
    if len(molecule_types) != 1:
        continue  # Not homo-n-mer
    
    # Check ring topology
    if not validate_ring_topology(cycle, graph):
        continue  # Invalid topology
        
    # Check interface consistency
    interface_types = get_ring_interface_types(cycle, graph)
    if not validate_interface_consistency(interface_types):
        continue  # Inconsistent interfaces
        
    # Create valid ring structure
    ring = create_ring_structure(cycle, graph)
    detected_rings.append(ring)
```

### Quality Assessment

**Fit Error Interpretation**:
```python
summary = regularizer.get_summary()
avg_error = summary['average_fit_error']

if avg_error < 0.5:
    print("Excellent fit - minimal regularization needed")
elif avg_error < 1.0:
    print("Good fit - small corrections applied")
elif avg_error < 2.0:
    print("Moderate fit - significant regularization applied")
else:
    print("Poor fit - check ring detection or geometry choice")
```

"""

from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from scipy.optimize import minimize
import numpy as np
import networkx as nx

from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.components.system import System
from .file_manager import WorkspaceManager


@dataclass
class RingStructure:
    """Represents a detected ring structure.

    Attributes:
        molecules: List of molecule instances in the ring (ordered)
        interface_types: Set of interface types used in ring connections
        ring_center: Calculated center of the ring
        ring_radius: Average radius of the ring
        ring_normal: Normal vector to the ring plane
    """
    molecules: List[MoleculeInstance]
    interface_types: Set[str]
    ring_center: np.ndarray
    ring_radius: float
    ring_normal: np.ndarray


@dataclass
class GeometricFit:
    """Represents a geometric fit (cylinder or sphere).

    Attributes:
        geometry_type: "cylinder" or "sphere"
        center: Center point of the geometry
        radius: Radius of the geometry
        axis: Axis vector (for cylinder only)
        fit_error: RMS error of the fit
    """
    geometry_type: str
    center: np.ndarray
    radius: float
    axis: Optional[np.ndarray] = None
    fit_error: float = 0.0


class RingRegularizer:
    """Regularizes ring-like homo-n-mer structures to ideal geometries.

    Attributes:
        system: ionerdss System to process
        workspace_manager: Workspace manager for logging
        mode: Regularization mode ("off", "separate", "uniform")
        geometry: Target geometry ("cylinder", "sphere")
        min_ring_size: Minimum ring size to consider (default 3)
    """

    def __init__(self, system: System, workspace_manager: Optional[WorkspaceManager] = None,
                 mode: str = "off", geometry: str = "cylinder", min_ring_size: int = 3):
        """Initialize ring regularizer.

        Args:
            system: ionerdss System to process
            workspace_manager: Workspace manager for logging
            mode: Regularization mode ("off", "separate", "uniform")
            geometry: Target geometry ("cylinder", "sphere")
            min_ring_size: Minimum ring size to consider
        """
        self.system = system
        self.workspace_manager = workspace_manager
        self.mode = mode.lower()
        self.geometry = geometry.lower()
        self.min_ring_size = max(3, min_ring_size)

        # Validate inputs
        if self.mode not in ["off", "separate", "uniform"]:
            raise ValueError(
                f"Invalid mode: {mode}. Must be 'off', 'separate', or 'uniform'")
        if self.geometry not in ["cylinder", "sphere"]:
            raise ValueError(
                f"Invalid geometry: {geometry}. Must be 'cylinder' or 'sphere'")

        # Storage for detected structures
        self.detected_rings: List[RingStructure] = []
        self.geometric_fits: Dict[int, GeometricFit] = {}  # ring_id -> fit

    def regularize(self) -> bool:
        """Perform ring regularization if mode is not 'off'.

        Returns:
            True if regularization was performed, False otherwise
        """
        if self.mode == "off":
            self.workspace_manager.logger.info(
                "Regularizer off"
            )
            return False

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Starting ring regularization: mode=%s, geometry=%s",
                self.mode, self.geometry
            )

        # Step 1: Detect ring structures
        self.detected_rings = self._detect_ring_structures()

        if not self.detected_rings:
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "No ring structures detected")
            return False

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Detected %d ring structures", len(self.detected_rings)
            )

        # Step 2: Calculate geometric fits
        if self.mode == "uniform":
            self._calculate_uniform_fit()
        else:  # separate mode
            self._calculate_separate_fits()

        # Step 3: Apply regularization
        self._apply_regularization()

        if self.workspace_manager:
            self.workspace_manager.logger.info("Ring regularization completed")

        return True

    def _detect_ring_structures(self) -> List[RingStructure]:
        """Detect homo-n-mer ring structures in the system.

        Returns:
            List of detected ring structures
        """
        rings = []

        # Build connectivity graph
        graph = self._build_connectivity_graph()

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Built connectivity graph with %d nodes and %d edges",
                graph.number_of_nodes(), graph.number_of_edges()
            )

        # Find all cycles in the graph
        try:
            # For directed graphs, find simple cycles
            cycles = list(nx.simple_cycles(graph))
        except Exception as e:
            # Fallback for undirected graph
            self.workspace_manager.logger.warning(
                "Fallback for undirected graph: %s", str(e)
            )
            cycles = []
            try:
                undirected_graph = graph.to_undirected()
                cycle_basis = nx.cycle_basis(undirected_graph)
                cycles = [cycle for cycle in cycle_basis if len(
                    cycle) >= self.min_ring_size]
            except Exception as e2:
                if self.workspace_manager:
                    self.workspace_manager.logger.warning(
                        "Failed to find cycles in graph: %s", str(e2)
                    )
                return rings

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Found %d potential cycles", len(cycles))

        # Filter and validate cycles
        for cycle in cycles:
            if len(cycle) >= self.min_ring_size:
                ring = self._validate_and_create_ring(cycle, graph)
                if ring:
                    rings.append(ring)
                    if self.workspace_manager:
                        self.workspace_manager.logger.info(
                            "Detected valid ring: %d molecules, interface types: %s",
                            len(ring.molecules), ring.interface_types
                        )

        return rings

    def _build_connectivity_graph(self) -> nx.DiGraph:
        """Build a connectivity graph from molecule instances and their interfaces.

        Returns:
            NetworkX directed graph representing molecular connectivity
        """
        graph = nx.DiGraph()

        # Add all molecule instances as nodes
        for mol_instance in self.system.molecule_instances:
            graph.add_node(mol_instance.name, molecule=mol_instance)

        # Add edges based on interface connections
        for mol_instance in self.system.molecule_instances:
            for interface_instance, partner_molecule in mol_instance.interfaces_neighbors_map.items():
                if interface_instance.interface_type and partner_molecule:
                    # Add edge with interface type information
                    interface_type_name = interface_instance.interface_type.get_name()
                    graph.add_edge(
                        mol_instance.name,
                        partner_molecule.name,
                        interface_type=interface_type_name,
                        interface_instance=interface_instance
                    )

        return graph

    def _validate_and_create_ring(self, cycle: List[str], graph: nx.DiGraph) -> Optional[RingStructure]:
        """Validate a cycle and create a RingStructure if valid.

        Args:
            cycle: List of molecule names forming a cycle
            graph: Connectivity graph

        Returns:
            RingStructure if valid, None otherwise
        """
        if len(cycle) < self.min_ring_size:
            return None

        # Check if all molecules are of the same type (homo-n-mer)
        molecule_instances = []
        molecule_types = set()

        for mol_name in cycle:
            if mol_name in graph.nodes:
                mol_instance = graph.nodes[mol_name]['molecule']
                molecule_instances.append(mol_instance)
                if mol_instance.molecule_type:
                    molecule_types.add(mol_instance.molecule_type.name)

        if len(molecule_types) != 1:
            return None  # Not homo-n-mer

        # Check interface consistency around the ring and validate ring topology
        interface_types_in_ring = set()
        ring_connections = {}  # mol_name -> list of (neighbor, interface_type)

        for i in range(len(cycle)):
            current = cycle[i]
            next_mol = cycle[(i + 1) % len(cycle)]

            if graph.has_edge(current, next_mol):
                edge_data = graph.get_edge_data(current, next_mol)
                interface_type = edge_data['interface_type']
                interface_types_in_ring.add(interface_type)

                # Track ring connections for topology validation
                if current not in ring_connections:
                    ring_connections[current] = []
                ring_connections[current].append((next_mol, interface_type))

        # Validate ring topology: each molecule should have exactly 2 connections
        # within the ring using the ring-forming interface types
        if not self._validate_ring_topology(cycle, ring_connections, interface_types_in_ring, graph):
            return None

        # Validate interface consistency
        if not self._validate_interface_consistency(interface_types_in_ring):
            return None

        # Calculate ring geometry
        positions = np.array([mol.com for mol in molecule_instances])
        center = np.mean(positions, axis=0)

        # Calculate ring normal using cross product of vectors
        if len(positions) >= 3:
            v1 = positions[1] - positions[0]
            v2 = positions[2] - positions[0]
            normal = np.cross(v1, v2)
            if np.linalg.norm(normal) > 1e-10:
                normal = normal / np.linalg.norm(normal)
            else:
                normal = np.array([0.0, 0.0, 1.0])  # Default
        else:
            normal = np.array([0.0, 0.0, 1.0])

        # Calculate average radius
        distances = np.linalg.norm(positions - center, axis=1)
        radius = np.mean(distances)

        return RingStructure(
            molecules=molecule_instances,
            interface_types=interface_types_in_ring,
            ring_center=center,
            ring_radius=radius,
            ring_normal=normal
        )

    def _validate_ring_topology(self, cycle: List[str], ring_connections: Dict[str, List[Tuple[str, str]]],
                                ring_interface_types: Set[str], graph: nx.DiGraph) -> bool:
        """Validate that the cycle forms a proper closed ring topology.

        Within the ring, each subunit should connect to exactly 2 neighbors using the 
        ring-forming interface types. Additional connections outside the ring are allowed.

        Args:
            cycle: List of molecule names in the cycle
            ring_connections: Mapping of molecule to its ring neighbors and interface types
            ring_interface_types: Set of interface types that form the ring
            graph: Full connectivity graph

        Returns:
            True if topology is valid, False otherwise
        """
        for mol_name in cycle:
            # Count connections within the ring using ring-forming interface types
            ring_neighbor_count = 0

            # Check all edges from this molecule
            if mol_name in graph:
                for neighbor in graph.neighbors(mol_name):
                    if neighbor in cycle:  # Neighbor is part of the ring
                        edge_data = graph.get_edge_data(mol_name, neighbor)
                        if edge_data and edge_data.get('interface_type') in ring_interface_types:
                            ring_neighbor_count += 1

            # Each molecule should have exactly 2 ring connections (previous and next in cycle)
            if ring_neighbor_count != 2:
                if self.workspace_manager:
                    self.workspace_manager.logger.debug(
                        "Invalid ring topology: molecule %s has %d ring connections (expected 2)",
                        mol_name, ring_neighbor_count
                    )
                return False

        # Additional validation: ensure the ring is actually closed
        # Check that each molecule connects to its expected neighbors in the cycle
        for i, mol_name in enumerate(cycle):
            prev_mol = cycle[(i - 1) % len(cycle)]
            next_mol = cycle[(i + 1) % len(cycle)]

            # Check if this molecule connects to both its ring neighbors
            has_prev_connection = False
            has_next_connection = False

            if graph.has_edge(mol_name, prev_mol):
                edge_data = graph.get_edge_data(mol_name, prev_mol)
                if edge_data and edge_data.get('interface_type') in ring_interface_types:
                    has_prev_connection = True

            if graph.has_edge(mol_name, next_mol):
                edge_data = graph.get_edge_data(mol_name, next_mol)
                if edge_data and edge_data.get('interface_type') in ring_interface_types:
                    has_next_connection = True

            # For undirected connections, also check reverse direction
            if not has_prev_connection and graph.has_edge(prev_mol, mol_name):
                edge_data = graph.get_edge_data(prev_mol, mol_name)
                if edge_data and edge_data.get('interface_type') in ring_interface_types:
                    has_prev_connection = True

            if not has_next_connection and graph.has_edge(next_mol, mol_name):
                edge_data = graph.get_edge_data(next_mol, mol_name)
                if edge_data and edge_data.get('interface_type') in ring_interface_types:
                    has_next_connection = True

            if not (has_prev_connection and has_next_connection):
                if self.workspace_manager:
                    self.workspace_manager.logger.debug(
                        "Incomplete ring: molecule %s missing connection to prev=%s or next=%s",
                        mol_name, prev_mol, next_mol
                    )
                return False

        return True

    def _validate_interface_consistency(self, interface_types: Set[str]) -> bool:
        """Validate that ring uses consistent interface types.

        Args:
            interface_types: Set of interface types used in the ring

        Returns:
            True if interfaces are consistent (all same or complementary pairs)
        """
        if len(interface_types) == 1:
            # All same interface type (homotypic)
            return True
        elif len(interface_types) == 2:
            # Check if they are complementary heterotypic interfaces
            types_list = list(interface_types)
            type1_parts = types_list[0].split('_')
            type2_parts = types_list[1].split('_')

            if len(type1_parts) >= 3 and len(type2_parts) >= 3:
                # Check if they are complementary (A_B_1 and B_A_1)
                if (type1_parts[0] == type2_parts[1] and
                    type1_parts[1] == type2_parts[0] and
                        type1_parts[2] == type2_parts[2]):
                    return True

        return False

    def _calculate_uniform_fit(self) -> None:
        """Calculate a single geometric fit for all rings."""
        if not self.detected_rings:
            return

        # Collect all ring positions
        all_positions = []
        for ring in self.detected_rings:
            positions = np.array([mol.com for mol in ring.molecules])
            all_positions.extend(positions)

        all_positions = np.array(all_positions)

        # Calculate best fit
        if self.geometry == "cylinder":
            fit = self._fit_cylinder(all_positions)
        else:  # sphere
            fit = self._fit_sphere(all_positions)

        # Apply the same fit to all rings
        for i in range(len(self.detected_rings)):
            self.geometric_fits[i] = fit

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                "Uniform %s fit: center=%s, radius=%.3f, error=%.3f",
                fit.geometry_type, fit.center, fit.radius, fit.fit_error
            )

    def _calculate_separate_fits(self) -> None:
        """Calculate separate geometric fits for each ring."""
        for i, ring in enumerate(self.detected_rings):
            positions = np.array([mol.com for mol in ring.molecules])

            if self.geometry == "cylinder":
                fit = self._fit_cylinder(positions)
            else:  # sphere
                fit = self._fit_sphere(positions)

            self.geometric_fits[i] = fit

            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    "Ring %d %s fit: center=%s, radius=%.3f, error=%.3f",
                    i, fit.geometry_type, fit.center, fit.radius, fit.fit_error
                )

    def _fit_cylinder(self, positions: np.ndarray) -> GeometricFit:
        """Fit a cylinder to the given positions.

        Args:
            positions: Array of 3D positions to fit

        Returns:
            GeometricFit object for cylinder
        """
        def cylinder_error(params):
            cx, cy, cz, ax, ay, az, radius = params
            center = np.array([cx, cy, cz])
            axis = np.array([ax, ay, az])
            axis = axis / (np.linalg.norm(axis) + 1e-10)  # Normalize

            # Calculate distances from points to cylinder axis
            errors = []
            for pos in positions:
                # Vector from center to point
                vec = pos - center
                # Project onto axis
                proj_length = np.dot(vec, axis)
                proj_point = center + proj_length * axis
                # Distance from point to axis
                dist = np.linalg.norm(pos - proj_point)
                errors.append((dist - radius) ** 2)

            return np.sum(errors)

        # Initial guess
        center_init = np.mean(positions, axis=0)
        axis_init = np.array([0.0, 0.0, 1.0])  # Default vertical axis
        radius_init = np.std(np.linalg.norm(positions - center_init, axis=1))

        initial_params = [
            center_init[0], center_init[1], center_init[2],
            axis_init[0], axis_init[1], axis_init[2],
            radius_init
        ]

        # Optimize
        try:
            result = minimize(cylinder_error, initial_params, method='BFGS')
            if result.success:
                cx, cy, cz, ax, ay, az, radius = result.x
                center = np.array([cx, cy, cz])
                axis = np.array([ax, ay, az])
                axis = axis / (np.linalg.norm(axis) + 1e-10)
                error = np.sqrt(result.fun / len(positions))
            else:
                raise RuntimeError("Optimization failed")
        except Exception:
            # Fallback to simple fit
            center = center_init
            axis = axis_init
            radius = radius_init
            error = 0.0

        return GeometricFit(
            geometry_type="cylinder",
            center=center,
            radius=abs(radius),
            axis=axis,
            fit_error=error
        )

    def _fit_sphere(self, positions: np.ndarray) -> GeometricFit:
        """Fit a sphere to the given positions.

        Args:
            positions: Array of 3D positions to fit

        Returns:
            GeometricFit object for sphere
        """
        def sphere_error(params):
            cx, cy, cz, radius = params
            center = np.array([cx, cy, cz])

            distances = np.linalg.norm(positions - center, axis=1)
            errors = (distances - radius) ** 2
            return np.sum(errors)

        # Initial guess
        center_init = np.mean(positions, axis=0)
        radius_init = np.mean(np.linalg.norm(positions - center_init, axis=1))

        initial_params = [center_init[0],
                          center_init[1], center_init[2], radius_init]

        # Optimize
        try:
            result = minimize(sphere_error, initial_params, method='BFGS')
            if result.success:
                cx, cy, cz, radius = result.x
                center = np.array([cx, cy, cz])
                error = np.sqrt(result.fun / len(positions))
            else:
                raise RuntimeError("Optimization failed")
        except Exception:
            # Fallback to simple fit
            center = center_init
            radius = radius_init
            error = 0.0

        return GeometricFit(
            geometry_type="sphere",
            center=center,
            radius=abs(radius),
            fit_error=error
        )

    def _apply_regularization(self) -> None:
        """Apply the calculated geometric fits to regularize ring structures."""
        for ring_id, ring in enumerate(self.detected_rings):
            if ring_id not in self.geometric_fits:
                continue

            fit = self.geometric_fits[ring_id]
            self._regularize_ring(ring, fit)

    def _regularize_ring(self, ring: RingStructure, fit: GeometricFit) -> None:
        """Regularize a single ring to the given geometric fit.

        Args:
            ring: Ring structure to regularize
            fit: Target geometric fit
        """
        n_molecules = len(ring.molecules)

        if fit.geometry_type == "cylinder":
            self._regularize_to_cylinder(ring, fit, n_molecules)
        else:  # sphere
            self._regularize_to_sphere(ring, fit, n_molecules)

    def _regularize_to_cylinder(self, ring: RingStructure, fit: GeometricFit, n_molecules: int) -> None:
        """Regularize ring molecules to cylinder surface.

        Args:
            ring: Ring structure to regularize
            fit: Cylinder fit parameters
            n_molecules: Number of molecules in ring
        """
        # Calculate ideal positions on cylinder
        angles = np.linspace(0, 2 * np.pi, n_molecules, endpoint=False)

        # Create rotation matrix to align with cylinder axis
        target_axis = fit.axis
        current_axis = np.array([0.0, 0.0, 1.0])

        if np.allclose(target_axis, current_axis):
            rotation_matrix = np.eye(3)
        else:
            # Calculate rotation matrix to align current_axis with target_axis
            v = np.cross(current_axis, target_axis)
            s = np.linalg.norm(v)
            c = np.dot(current_axis, target_axis)

            if s > 1e-10:
                vx = np.array(
                    [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                rotation_matrix = np.eye(
                    3) + vx + np.dot(vx, vx) * ((1 - c) / (s * s))
            else:
                rotation_matrix = np.eye(3)

        # Calculate new positions
        for i, molecule in enumerate(ring.molecules):
            angle = angles[i]

            # Position on unit cylinder
            local_pos = np.array(
                [np.cos(angle), np.sin(angle), 0.0]) * fit.radius

            # Rotate to align with cylinder axis
            rotated_pos = np.dot(rotation_matrix, local_pos)

            # Translate to cylinder center (all on same z-plane)
            new_com = fit.center + rotated_pos

            # Calculate translation vector
            translation = new_com - molecule.com

            # Apply translation to molecule
            self._translate_molecule(molecule, translation)

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Regularized ring to cylinder: center=%s, radius=%.3f, axis=%s",
                fit.center, fit.radius, fit.axis
            )

    def _regularize_to_sphere(self, ring: RingStructure, fit: GeometricFit, n_molecules: int) -> None:
        """Regularize ring molecules to sphere surface.

        Args:
            ring: Ring structure to regularize
            fit: Sphere fit parameters
            n_molecules: Number of molecules in ring
        """
        # Distribute molecules evenly on sphere surface
        # Use spherical coordinates with fixed latitude (ring on sphere)

        # Calculate current ring plane
        current_positions = np.array([mol.com for mol in ring.molecules])
        current_center = np.mean(current_positions, axis=0)

        # Project current positions to sphere
        angles = np.linspace(0, 2 * np.pi, n_molecules, endpoint=False)

        # Keep molecules in their current plane but project to sphere surface
        plane_normal = ring.ring_normal

        for i, molecule in enumerate(ring.molecules):
            angle = angles[i]

            # Create position on unit circle in xy-plane
            local_pos = np.array([np.cos(angle), np.sin(angle), 0.0])

            # Rotate to align with ring plane
            if not np.allclose(plane_normal, [0, 0, 1]):
                # Calculate rotation matrix
                target_normal = plane_normal
                current_normal = np.array([0.0, 0.0, 1.0])

                v = np.cross(current_normal, target_normal)
                s = np.linalg.norm(v)
                c = np.dot(current_normal, target_normal)

                if s > 1e-10:
                    vx = np.array(
                        [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                    rotation_matrix = np.eye(
                        3) + vx + np.dot(vx, vx) * ((1 - c) / (s * s))
                    local_pos = np.dot(rotation_matrix, local_pos)

            # Scale to sphere radius and translate to center
            new_com = fit.center + local_pos * fit.radius

            # Calculate translation vector
            translation = new_com - molecule.com

            # Apply translation to molecule
            self._translate_molecule(molecule, translation)

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                "Regularized ring to sphere: center=%s, radius=%.3f",
                fit.center, fit.radius
            )

    def _translate_molecule(self, molecule: MoleculeInstance, translation: np.ndarray) -> None:
        """Translate a molecule instance and all its interface coordinates.

        Args:
            molecule: Molecule instance to translate
            translation: Translation vector
        """
        # Translate COM
        molecule.com += translation

        # Translate all interface coordinates
        for interface_instance in molecule.interfaces_neighbors_map.keys():
            if hasattr(interface_instance, 'absolute_coord'):
                interface_instance.absolute_coord += translation

            # Also update interface type coordinates if they reference this molecule
            if interface_instance.interface_type:
                if hasattr(interface_instance.interface_type, 'absolute_coord'):
                    interface_instance.interface_type.absolute_coord += translation

    def get_summary(self) -> Dict[str, any]:
        """Get summary of ring regularization results.

        Returns:
            Dictionary with regularization statistics
        """
        return {
            "mode": self.mode,
            "geometry": self.geometry,
            "min_ring_size": self.min_ring_size,
            "rings_detected": len(self.detected_rings),
            "rings_regularized": len(self.geometric_fits),
            "ring_sizes": [len(ring.molecules) for ring in self.detected_rings],
            "fit_errors": [fit.fit_error for fit in self.geometric_fits.values()],
            "average_fit_error": np.mean([fit.fit_error\
                for fit in self.geometric_fits.values()])\
                    if self.geometric_fits else 0.0
        }
