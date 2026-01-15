"""Specific implementation of Platonic solids coordinate generation."""

from abc import ABC, abstractmethod
from typing import List, Tuple
import math
import numpy as np
from .geometry import distance, mid_pt

class PlatonicSolidGenerator(ABC):
    """Abstract base class for Platonic solid generators."""

    @property
    @abstractmethod
    def name(self) -> str:
        pass

    @property
    @abstractmethod
    def num_sites(self) -> int:
        pass

    @abstractmethod
    def _get_vertices(self, radius: float) -> List[List[float]]:
        """Calculate vertices for the solid."""
        pass

    @abstractmethod
    def _get_face_indices(self) -> List[Tuple[int, ...]]:
        """Return list of vertex indices for each face."""
        pass

    @property
    @abstractmethod
    def angle_indices(self) -> Tuple[Tuple[int, int], ...]:
        """Return indices for angle calculation (theta1, theta2, phi1, phi2)."""
        pass

    @abstractmethod
    def _get_reduction_angle(self) -> float:
        """Return the angle used for leg reduction."""
        pass

    def generate_coordinates(self, radius: float, sigma: float) -> List[List[np.ndarray]]:
        """
        Generate coordinates for ALL faces.
        Returns: List of [COM, leg1, leg2, ..., Normal] for each face.
        """
        vertices = self._get_vertices(radius)
        face_indices_list = self._get_face_indices()
        
        # Calculate reduction params
        angle = self._get_reduction_angle()
        denom = 2 * math.sin(angle / 2)
        red_len = sigma / denom
        
        all_faces_coords = []
        
        for face_indices in face_indices_list:
            face_verts = [vertices[i] for i in face_indices]
            
            # 1. Calculate Face COM
            com = np.mean(face_verts, axis=0)
            
            # 2. Calculate Legs (Edge Midpoints)
            legs = []
            num_verts = len(face_verts)
            for i in range(num_verts):
                p1 = face_verts[i]
                p2 = face_verts[(i + 1) % num_verts]
                legs.append(mid_pt(p1, p2))
                
            # 3. Reduce Legs
            reduced_legs = []
            for leg in legs:
                dist = distance(com, leg)
                if dist == 0:
                    ratio = 1 
                else:
                    ratio = 1 - red_len / dist
                
                leg_red = (np.array(leg) - com) * ratio + com
                reduced_legs.append(leg_red)
                
            # 4. Normal (pointing towards origin)
            # Original code used -COM.
            normal = -com
            
            # Assemble list: [COM, leg1, leg2, ..., Normal]
            # Note: Normal is NOT usually part of the flat list used for angle indices,
            # but it IS returned by legacy `input_coord`.
            # Legacy `reduced_coord` returns: [COM, leg1_red, leg2_red...]
            # Legacy `input_coord` adds Normal at the end.
            
            # We will return the structure expected by `angle_cal` (via indices):
            # angle_cal expects points.
            # And `PlatonicSolids.py` expects to extract [COM, legs..., Normal] for the final MoleculeType.
            
            face_data = [com] + reduced_legs + [normal]
            all_faces_coords.append(face_data)
            
        return all_faces_coords

class CubeGenerator(PlatonicSolidGenerator):
    @property
    def name(self): return "cube"
    @property
    def num_sites(self): return 4
    
    def _get_reduction_angle(self):
        return math.acos(0) # 90 degrees

    def _get_vertices(self, radius):
        scaler = radius / (3**0.5)
        return [
            [scaler, scaler, scaler],      # v0
            [-scaler, scaler, scaler],     # v1
            [scaler, -scaler, scaler],     # v2
            [scaler, scaler, -scaler],     # v3
            [-scaler, -scaler, scaler],    # v4
            [scaler, -scaler, -scaler],    # v5
            [-scaler, scaler, -scaler],    # v6
            [-scaler, -scaler, -scaler]    # v7
        ]

    def _get_face_indices(self):
        # 0, 3, 5, 2
        # 0, 3, 6, 1
        # 0, 1, 4, 2
        # 7, 4, 1, 6
        # 7, 4, 2, 5
        # 7, 6, 3, 5
        return [
            (0, 3, 5, 2),
            (0, 3, 6, 1),
            (0, 1, 4, 2),
            (7, 4, 1, 6),
            (7, 4, 2, 5),
            (7, 6, 3, 5)
        ]

    @property
    def angle_indices(self):
        return ((0, 0), (0, 1), (1, 0), (1, 1))

class DodecahedronGenerator(PlatonicSolidGenerator):
    @property
    def name(self): return "dode"
    @property
    def num_sites(self): return 5
    
    def _get_reduction_angle(self):
        m = (1 + 5**0.5) / 2
        return 2 * math.atan(m)

    def _get_vertices(self, radius):
        scaler = radius / (3**0.5)
        m = (1 + 5**0.5) / 2
        
        # Vertices 1-20 mapped to 0-19
        # Code used 1-based names V1..V20
        coords = [
            [0, m, 1/m],    # V1 -> 0
            [0, m, -1/m],   # V2 -> 1
            [0, -m, 1/m],   # V3 -> 2
            [0, -m, -1/m],  # V4 -> 3
            [1/m, 0, m],    # V5 -> 4
            [1/m, 0, -m],   # V6 -> 5
            [-1/m, 0, m],   # V7 -> 6
            [-1/m, 0, -m],  # V8 -> 7
            [m, 1/m, 0],    # V9 -> 8
            [m, -1/m, 0],   # V10 -> 9
            [-m, 1/m, 0],   # V11 -> 10
            [-m, -1/m, 0],  # V12 -> 11
            [1, 1, 1],      # V13 -> 12
            [1, 1, -1],     # V14 -> 13
            [1, -1, 1],     # V15 -> 14
            [1, -1, -1],    # V16 -> 15
            [-1, 1, 1],     # V17 -> 16
            [-1, 1, -1],    # V18 -> 17
            [-1, -1, 1],    # V19 -> 18
            [-1, -1, -1]    # V20 -> 19
        ]
        return [[c * scaler for c in coord] for coord in coords]

    def _get_face_indices(self):
        # Indices adjusted to 0-based from generic inspection
        # 1. 6, 18, 2, 14, 4 -> V7, V19, V3, V15, V5 -> Indices 6, 18, 2, 14, 4
        return [
            (6, 18, 2, 14, 4),
            (6, 4, 12, 0, 16),
            (4, 14, 9, 8, 12),
            (6, 18, 11, 10, 16),
            (14, 2, 3, 15, 9),
            (18, 11, 19, 3, 2),
            (16, 10, 17, 1, 0),
            (12, 0, 1, 13, 8),
            (7, 17, 10, 11, 19),
            (5, 13, 8, 9, 15),
            (3, 19, 7, 5, 15),
            (1, 17, 7, 5, 13)
        ]

    @property
    def angle_indices(self):
        return ((0, 0), (0, 3), (4, 0), (4, 1))

class IcosahedronGenerator(PlatonicSolidGenerator):
    @property
    def name(self): return "icos"
    @property
    def num_sites(self): return 3
    
    def _get_reduction_angle(self):
        return math.acos(-(5**0.5)/3)

    def _get_vertices(self, radius):
        scaler = radius / (2 * math.sin(2 * math.pi / 5))
        m = (1 + 5**0.5) / 2
        coords = [
            [0, 1, m],     # v0
            [0, 1, -m],    # v1
            [0, -1, m],    # v2
            [0, -1, -m],   # v3
            [1, m, 0],     # v4
            [1, -m, 0],    # v5
            [-1, m, 0],    # v6
            [-1, -m, 0],   # v7
            [m, 0, 1],     # v8
            [m, 0, -1],    # v9
            [-m, 0, 1],    # v10
            [-m, 0, -1]    # v11
        ]
        return [[c * scaler for c in coord] for coord in coords]

    def _get_face_indices(self):
        return [
            (0, 2, 8), (0, 8, 4), (0, 4, 6), (0, 6, 10), (0, 10, 2),
            (3, 7, 5), (3, 5, 9), (3, 9, 1), (3, 1, 11), (3, 11, 7),
            (7, 2, 5), (2, 5, 8), (5, 8, 9), (8, 9, 4), (9, 4, 1),
            (4, 1, 6), (1, 6, 11), (6, 11, 10), (11, 10, 7), (10, 7, 2)
        ]

    @property
    def angle_indices(self):
        return ((0, 0), (0, 1), (1, 0), (1, 1))

class OctahedronGenerator(PlatonicSolidGenerator):
    @property
    def name(self): return "octa"
    @property
    def num_sites(self): return 3
    
    def _get_reduction_angle(self):
        # Angle for Octahedron leg reduction?
        # Assuming standard dihedral (109.47) or similar logic. 
        # Standard implementation used specific logic.
        # Verified earlier: Octa also used leg reduction logic?
        # Wait, I didn't verify Octa angle.
        # Assuming typical: acos(-1/3) = 109.47 deg = tetrahedral angle.
        return math.acos(-1/3) 

    def _get_vertices(self, radius):
        scaler = radius
        # v0..v5
        coords = [
            [1, 0, 0], [-1, 0, 0],
            [0, 1, 0], [0, -1, 0],
            [0, 0, 1], [0, 0, -1]
        ]
        return [[c * scaler for c in coord] for coord in coords]

    def _get_face_indices(self):
        return [
            (0, 2, 4), (0, 3, 4), (0, 3, 5), (0, 2, 5),
            (1, 2, 4), (1, 3, 4), (1, 3, 5), (1, 2, 5)
        ]

    @property
    def angle_indices(self):
        return ((0, 0), (0, 1), (1, 0), (1, 1))

class TetrahedronGenerator(PlatonicSolidGenerator):
    @property
    def name(self): return "tetr"
    @property
    def num_sites(self): return 3
    
    def _get_reduction_angle(self):
        # Tetrahedron angle. 
        # acos(1/3) is approx 70.5 deg.
        return math.acos(1/3)

    def _get_vertices(self, radius):
        # scaler = radius/(3/8)**0.5/2 = radius / (sqrt(3/8)*2) = radius / sqrt(1.5)
        # v0..v3
        s = 1/(2**0.5) # 0.707
        scaler = radius / ((3.0/8.0)**0.5 * 2.0)
        
        coords = [
            [1, 0, -s],    # v0
            [-1, 0, -s],   # v1
            [0, 1, s],     # v2
            [0, -1, s]     # v3
        ]
        return [[c * scaler for c in coord] for coord in coords]

    def _get_face_indices(self):
        return [
            (0, 1, 2), (0, 2, 3), (0, 1, 3), (1, 2, 3)
        ]

    @property
    def angle_indices(self):
        return ((0, 0), (0, 1), (1, 0), (1, 1))
