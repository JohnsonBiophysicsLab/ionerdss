"""
ionerdss.model.pdb.parser

PDB/mmCIF file parsing and structure processing with workspace integration.

This module handles the parsing of PDB and mmCIF files using BioPython,
with proper file management and workspace organization. It provides a unified
interface for loading protein structures from local files or directly from the
Protein Data Bank, extracting essential geometric and sequence information for
downstream processing.


## Key Features

### Automatic Format Detection
- **PDB format**: `.pdb`, `.ent` files
- **mmCIF format**: `.cif`, `.mmcif` files
- **Auto-detection**: Based on file extension and content

### Multiple Input Sources
- **Local files**: Parse existing structure files
- **PDB database**: Automatic download by PDB ID
- **Workspace integration**: Organized file management

### Comprehensive Data Extraction
- **Geometric data**: Coordinates, center of mass, radius, bounding boxes
- **Sequence data**: Amino acid sequences from structure
- **Chain validation**: Automatic filtering of valid protein chains
- **Residue processing**: Standard amino acid extraction with Cα coordinates

## Basic Usage

### Parse Local File

```python
from ionerdss.model.pdb.parser import PDBParser

# Parse local PDB file
parser = PDBParser("structure.pdb")

# Parse local mmCIF file
parser = PDBParser("structure.cif")

# Get basic information
print(f"PDB ID: {parser.get_pdb_id()}")
print(f"Chains: {parser.get_chain_ids()}")
```

### Download from PDB Database

```python
# Download and parse by PDB ID
parser = PDBParser("1ABC", fetch_from_pdb=True)

# Specify format (default is mmCIF)
parser = PDBParser("1ABC", fetch_from_pdb=True, file_format="pdb")

# Auto-detection (4-character alphanumeric = PDB ID)
parser = PDBParser("1ABC")  # Automatically fetches if file doesn't exist locally
```

### With Workspace Management

```python
from ionerdss.model.pdb.file_manager import WorkspaceManager

# Create workspace
with WorkspaceManager("/path/to/workspace", "1ABC") as workspace:
    # Parser integrates with workspace
    parser = PDBParser("1ABC", fetch_from_pdb=True, workspace_manager=workspace)
    
    # Files are automatically organized in workspace
    # Downloaded structures go to: workspace/structures/downloaded/
    # Logs are written to: workspace/logs/pipeline.log
```

## File Format Support

### Supported Formats

| Format | Extensions | Description | Use Case |
|--------|------------|-------------|----------|
| **PDB** | `.pdb`, `.ent` | Legacy text format | Older structures, simple parsing |
| **mmCIF** | `.cif`, `.mmcif` | Modern structured format | New structures, complex assemblies |

### Format Selection for Downloads

```python
# Download mmCIF (recommended - more complete metadata)
parser = PDBParser("1ABC", fetch_from_pdb=True, file_format="mmcif")

# Download PDB format
parser = PDBParser("1ABC", fetch_from_pdb=True, file_format="pdb")
```

**Recommendation**: Use mmCIF format for new work as it contains more complete
structural metadata and handles large assemblies better.

## Data Extraction

### Chain Information

```python
# Get all valid chain IDs
chain_ids = parser.get_chain_ids()
print(f"Available chains: {chain_ids}")

# Get detailed chain data
for chain_id in chain_ids:
    chain_data = parser.get_chain_data(chain_id)
    
    print(f"Chain {chain_id}:")
    print(f"  Sequence: {chain_data['sequence']}")
    print(f"  Residues: {len(chain_data['residues'])}")
    print(f"  Center of mass: {chain_data['com']}")
    print(f"  Radius: {chain_data['radius']:.2f} Å")
```

### Chain Data Structure

Each chain contains the following information:

```python
chain_data = {
    'id': 'A',                           # Chain identifier
    'sequence': 'MKLAVQNCTGRLKDE...',     # Amino acid sequence
    'residues': [                        # List of residue information
        {
            'id': 1,                     # Residue number
            'name': 'MET',               # Three-letter code
            'ca_coord': [x, y, z]        # Cα coordinates (Å)
        },
        # ... more residues
    ],
    'ca_coords': np.array(...),          # All Cα coordinates (N×3)
    'all_coords': np.array(...),         # All atom coordinates (M×3)
    'com': np.array([x, y, z]),          # Center of mass (Å)
    'radius': 15.2,                      # RMS radius from COM (Å)
    'bbox_min': np.array([x, y, z]),     # Bounding box minimum (Å)
    'bbox_max': np.array([x, y, z])      # Bounding box maximum (Å)
}
```

### Geometric Calculations

**Center of Mass (COM)**:
```
COM = (1/N) × Σᵢ rᵢ
```
Unweighted geometric center of all atoms in the chain.

**Radius Calculation**:
```
radius = √[(1/N) × Σᵢ ||rᵢ - COM||²]
```
Root-mean-square distance from center of mass (coarse-grained radius).

**Bounding Box**:
```
bbox_min = [min(x), min(y), min(z)]
bbox_max = [max(x), max(y), max(z)]
```
Axis-aligned bounding box containing all atoms.

### Coordinate Access

```python
chain_data = parser.get_chain_data('A')

# Cα coordinates only (for backbone analysis)
ca_coords = chain_data['ca_coords']  # Shape: (n_residues, 3)

# All atom coordinates (for detailed calculations)
all_coords = chain_data['all_coords']  # Shape: (n_atoms, 3)

# Individual residue access
for residue in chain_data['residues']:
    residue_id = residue['id']
    residue_name = residue['name']
    ca_position = residue['ca_coord']
    print(f"Residue {residue_id} ({residue_name}): Cα at {ca_position}")
```

## Workspace Integration

### Automatic File Organization

When using a `WorkspaceManager`, files are automatically organized:

```
workspace/
├── structures/
│   ├── downloaded/
│   │   ├── 1abc.cif          # Downloaded from PDB
│   │   └── local_copy.pdb    # Copied from local file
│   └── processed/            # For processed structures
├── logs/
│   └── pipeline.log          # Parsing logs and progress
└── outputs/                  # For downstream results
```

### Logging Integration

```python
# Automatic logging with workspace
with WorkspaceManager("/workspace", "1ABC") as workspace:
    parser = PDBParser("1ABC", workspace_manager=workspace)
    
    # Log messages are automatically written:
    # "Downloading 1ABC in mmcif format..."
    # "Downloaded 1ABC to /workspace/structures/downloaded/1abc.cif"
    # "Structure parsed successfully"
    # "Processed 4 valid chains: ['A', 'B', 'C', 'D']"
```

### File Path Management

```python
# Get the actual file path used
file_path = parser.filepath
print(f"Structure file: {file_path}")

# Check if structure was downloaded
if parser.pdb_id:
    print(f"PDB ID: {parser.pdb_id}")
```

## Coordinate Systems

### Unit Handling

The parser maintains coordinates in **Angstroms** during parsing (BioPython
standard), with conversion utilities for downstream processing:

```python
# Coordinates are in Angstroms
chain_data = parser.get_chain_data('A')
com_angstrom = chain_data['com']  # [x, y, z] in Å

# Convert to nanometers for NERDSS
com_nm = parser.convert_coords_to_nm(com_angstrom)  # [x, y, z] in nm

# Convert distances
distance_nm = 2.5  # nm
distance_angstrom = parser.convert_distance_to_angstrom(distance_nm)  # 25.0 Å
```

### Coordinate Validation

```python
# Check coordinate completeness
chain_data = parser.get_chain_data('A')

print(f"Total atoms: {len(chain_data['all_coords'])}")
print(f"Cα atoms: {len(chain_data['ca_coords'])}")
print(f"Residues with Cα: {len(chain_data['residues'])}")

# Verify coordinate ranges
coords = chain_data['all_coords']
if len(coords) > 0:
    coord_range = coords.max() - coords.min()
    print(f"Coordinate range: {coord_range:.1f} Å")
    
    if coord_range > 1000:
        print("Warning: Unusually large coordinate range detected")
```


"""

from typing import Any, Dict, List, Tuple, Optional, Union
from pathlib import Path
import shutil
import tempfile

import numpy as np
from Bio.PDB import PDBParser as BioPDBParser, MMCIFParser, PDBList
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
# from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import PPBuilder, is_aa

from ionerdss.model.components.units import Units
from .file_manager import WorkspaceManager


class PDBParser:
    """Parser for PDB and mmCIF files with workspace integration.

    Handles file parsing, chain filtering, and coordinate extraction
    with proper file management through WorkspaceManager.

    Attributes:
        structure: BioPython Structure object from parsed file.
        units: Unit system (coordinates remain in Angstroms during parsing).
        chain_data: Processed chain information including coordinates and metadata.
        pdb_id: PDB identifier if fetched from database.
        filepath: Path to the structure file (in workspace).
        workspace_manager: Workspace manager for file organization.
    """

    def __init__(self, source: Union[str, Path], units: Optional[Units] = None,
                 fetch_from_pdb: bool = False, file_format: str = 'mmcif',
                 workspace_manager: Optional[WorkspaceManager] = None,
                 concat_all_frames=True, max_frames=None):
        """Initialize parser with PDB/mmCIF file or PDB ID for fetching.

        Args:
            source: Either a file path or PDB ID (4-character code).
            units: Unit system for the model. Defaults to standard units.
            fetch_from_pdb: If True, treat source as PDB ID and fetch from database.
            file_format: Format for fetching ('pdb' or 'mmcif'). Default 'mmcif'.
            workspace_manager: Workspace manager for file organization.
        """
        self.units = units or Units()
        self.structure: Optional[Structure] = None
        self.chain_data: Dict[str, dict] = {}
        self.pdb_id: Optional[str] = None
        self.filepath: Optional[Path] = None
        self.workspace_manager = workspace_manager

        # Dealing with PDB with multiple frames
        self.concat_all_frames = concat_all_frames
        self.max_frames = max_frames
        self.frames = {}  # Store multiple frames: {frame_num: structure_data}
        self.frame_count = 0

        # Handle source - either file path or PDB ID
        if fetch_from_pdb or self._looks_like_pdb_id(str(source)):
            self.pdb_id = str(source).upper()
            self.filepath = self._fetch_structure(self.pdb_id, file_format)
        else:
            self.filepath = Path(source)
            self.pdb_id = self._extract_pdb_id_from_filename()

            # Copy local file to workspace if workspace manager is provided
            if self.workspace_manager and self.filepath.exists():
                workspace_path = self.workspace_manager.copy_file_to_workspace(
                    self.filepath, 'downloaded'
                )
                self.filepath = workspace_path

        # Log file information
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                f"Using structure file: {self.filepath}")
            if self.filepath.exists():
                size = self.filepath.stat().st_size
                self.workspace_manager.logger.info(f"File size: {size} bytes")

        # Parse the structure
        self._parse_structure()

        # Extract chain data
        self._extract_chain_data()

    def _looks_like_pdb_id(self, source: str) -> bool:
        """Check if source looks like a PDB ID (4 alphanumeric characters).

        Args:
            source: Source string to check.

        Returns:
            True if source appears to be a PDB ID.
        """
        # Remove any file extensions
        source_clean = Path(source).stem.upper()

        # Check if it's exactly 4 alphanumeric characters
        return (len(source_clean) == 4 and
                source_clean.isalnum() and
                not Path(source).exists())  # And file doesn't exist locally

    def _fetch_structure(self, pdb_id: str, file_format: str = 'mmcif') -> Path:
        """Fetch structure from Protein Data Bank.

        Args:
            pdb_id: 4-character PDB identifier.
            file_format: Format to download ('pdb' or 'mmcif').

        Returns:
            Path to downloaded file in workspace.

        Raises:
            ValueError: If PDB ID is invalid or download fails.
        """
        if len(pdb_id) != 4 or not pdb_id.isalnum():
            raise ValueError(
                f"Invalid PDB ID: {pdb_id}. Must be 4 alphanumeric characters.")

        # Log download attempt
        if self.workspace_manager:
            self.workspace_manager.logger.info(
                f"Downloading {pdb_id} in {file_format} format...")

        # Get target path in workspace
        if self.workspace_manager:
            target_path = self.workspace_manager.get_structure_download_path(
                pdb_id, file_format)
        else:
            # Fallback to temp directory if no workspace manager
            temp_dir = Path(tempfile.mkdtemp(prefix=f"pdb_{pdb_id}_"))
            if file_format.lower() == 'mmcif':
                target_path = temp_dir / f"{pdb_id.lower()}.cif"
            else:
                target_path = temp_dir / f"{pdb_id.lower()}.pdb"

        # Create temporary directory for download
        temp_dir = Path(tempfile.mkdtemp(prefix=f"pdb_download_{pdb_id}_"))

        try:
            # Initialize PDB downloader with HTTPS server (more reliable than FTP)
            pdb_list = PDBList(server='https://files.rcsb.org')

            if file_format.lower() == 'mmcif':
                # Download mmCIF file
                downloaded_file = pdb_list.retrieve_pdb_file(
                    pdb_id,
                    pdir=str(temp_dir),
                    file_format='mmCif'
                )
            else:  # pdb format
                # Download PDB file
                downloaded_file = pdb_list.retrieve_pdb_file(
                    pdb_id,
                    pdir=str(temp_dir),
                    file_format='pdb'
                )

            # The downloaded file path returned by BioPython may not exist
            # Check what was actually downloaded in the temp directory
            downloaded_path = Path(downloaded_file)
            
            if not downloaded_path.exists():
                # BioPython may download with different naming (e.g., assembly files)
                # Search for any file in temp_dir that matches the pattern
                if file_format.lower() == 'mmcif':
                    pattern = f"{pdb_id.lower()}*.cif"
                else:
                    pattern = f"{pdb_id.lower()}*.pdb"
                
                matching_files = list(temp_dir.glob(pattern))
                
                if matching_files:
                    # Use the first matching file
                    downloaded_path = matching_files[0]
                    if self.workspace_manager:
                        self.workspace_manager.logger.info(
                            f"Found downloaded file: {downloaded_path.name}")
                else:
                    raise ValueError(
                        f"Failed to download PDB structure {pdb_id}. "
                        f"Expected file not found in {temp_dir}")

            # Move downloaded file to workspace
            target_path.parent.mkdir(parents=True, exist_ok=True)
            # Use shutil.move instead of rename to handle cross-filesystem moves
            shutil.move(str(downloaded_path), str(target_path))

            # Clean up temp directory
            shutil.rmtree(temp_dir, ignore_errors=True)

            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    f"Downloaded {pdb_id} to {target_path}")

            return target_path

        except Exception as e:
            # Clean up temp directory on failure
            shutil.rmtree(temp_dir, ignore_errors=True)

            if self.workspace_manager:
                self.workspace_manager.logger.error(
                    f"Failed to download {pdb_id}: {str(e)}")

            raise ValueError(f"Failed to fetch PDB {pdb_id}: {str(e)}") from e

    def _extract_pdb_id_from_filename(self) -> Optional[str]:
        """Extract PDB ID from filename if possible.

        Returns:
            PDB ID if extractable, None otherwise.
        """
        if not self.filepath:
            return None

        stem = self.filepath.stem.upper()

        # Try to extract 4-character PDB ID from filename
        if len(stem) >= 4:
            # Check if first 4 characters look like PDB ID
            potential_id = stem[:4]
            if potential_id.isalnum():
                return potential_id

        return None

    def _parse_structure(self) -> None:
        """Parse PDB or mmCIF file using appropriate BioPython parser."""
        if not self.filepath or not self.filepath.exists():
            raise FileNotFoundError(f"Structure file not found: {self.filepath}")

        if self.workspace_manager:
            self.workspace_manager.logger.info(f"Parsing structure file: {self.filepath.name}")

        try:
            if self.filepath.suffix.lower() in ['.pdb', '.ent']:
                parser = BioPDBParser(QUIET=True)
                self.structure = parser.get_structure('structure', self.filepath)
            elif self.filepath.suffix.lower() in ['.cif', '.mmcif']:
                parser = MMCIFParser(QUIET=True)
                self.structure = parser.get_structure('structure', self.filepath)
            else:
                # Try to guess format from content or use mmCIF as default for downloaded files
                if self.pdb_id:
                    # Assume mmCIF for downloaded files
                    parser = MMCIFParser(QUIET=True)
                    self.structure = parser.get_structure('structure', self.filepath)
                else:
                    raise ValueError(f"Unsupported file format: {self.filepath.suffix}")

            # Handle frame concatenation if requested
            if self.concat_all_frames:
                self._concatenate_all_frames()

            if self.workspace_manager:
                if self.concat_all_frames:
                    models = list(self.structure.get_models())
                    total_models = len(models)
                    self.workspace_manager.logger.info(
                        f"Structure parsed successfully - concatenated {total_models} frames into single frame")
                else:
                    self.workspace_manager.logger.info("Structure parsed successfully")

        except Exception as e:
            if self.workspace_manager:
                self.workspace_manager.logger.error(f"Failed to parse structure: {str(e)}")
            raise ValueError(f"Failed to parse structure file {self.filepath}: {str(e)}") from e
    
    def _extract_chain_data(self) -> None:
        """Extract and process chain data from parsed structure."""
        if not self.structure:
            raise RuntimeError("Structure not parsed")

        if self.workspace_manager:
            self.workspace_manager.logger.info("Extracting chain data...")

        # Get first model (most PDB files have only one)
        model = self.structure[0]

        # Process each chain
        valid_chains = 0
        for chain in model:
            if self._is_valid_chain(chain):
                chain_id = chain.get_id()
                self.chain_data[chain_id] = self._process_chain(chain)
                valid_chains += 1

        # Sort chain IDs for deterministic ordering
        self.chain_data = dict(sorted(self.chain_data.items()))

        if self.workspace_manager:
            self.workspace_manager.logger.info(
                f"Processed {valid_chains} valid chains: {list(self.chain_data.keys())}")

    def _is_valid_chain(self, chain: Chain) -> bool:
        """Check if chain contains at least one standard amino acid.

        Args:
            chain: BioPython Chain object.

        Returns:
            True if chain contains standard amino acids, False otherwise.
        """
        for residue in chain:
            if is_aa(residue, standard=True):
                return True
        return False

    def _process_chain(self, chain: Chain) -> dict:
        """Process a single chain to extract coordinates and metadata.

        Args:
            chain: BioPython Chain object.

        Returns:
            Dictionary containing chain data including coordinates,
            residues, COM, and radius.
        """
        # Extract amino acid residues with Cα atoms
        residues = []
        all_coords = []
        ca_coords = []

        for residue in chain:
            if is_aa(residue, standard=True):
                # Get Cα coordinate if available
                if 'CA' in residue:
                    ca_coord = residue['CA'].get_coord()
                    residue_data = {
                        'id': residue.get_id()[1],  # Residue number
                        'name': residue.get_resname(),
                        'ca_coord': ca_coord
                    }
                    residues.append(residue_data)
                    ca_coords.append(ca_coord)

                # Collect all atom coordinates for COM/radius calculation
                for atom in residue:
                    all_coords.append(atom.get_coord())

        # Convert to numpy arrays
        all_coords = np.array(all_coords) if all_coords else np.empty((0, 3))
        ca_coords = np.array(ca_coords) if ca_coords else np.empty((0, 3))

        # Calculate chain COM and radius (unweighted by mass)
        com = self._calculate_com(all_coords)
        radius = self._calculate_radius(all_coords, com)

        # Calculate axis-aligned bounding box
        bbox_min, bbox_max = self._calculate_bounding_box(all_coords)

        chain_data = {
            'id': chain.get_id(),
            'residues': residues,
            'ca_coords': ca_coords,
            'all_coords': all_coords,
            'com': com,
            'radius': radius,
            'bbox_min': bbox_min,
            'bbox_max': bbox_max,
            'sequence': self._extract_sequence(chain)
        }

        if self.workspace_manager:
            self.workspace_manager.logger.debug(
                f"Chain {chain.get_id()}: {len(residues)} residues, "
                f"radius={radius:.2f}Å, sequence_length={len(chain_data['sequence'])}"
            )

        return chain_data

    def _calculate_com(self, coords: np.ndarray) -> np.ndarray:
        """Calculate geometric center of mass (unweighted).

        Args:
            coords: Array of atomic coordinates (N, 3).

        Returns:
            Center of mass coordinates (3,).
        """
        if len(coords) == 0:
            return np.zeros(3)
        return coords.mean(axis=0)

    def _calculate_radius(self, coords: np.ndarray, com: np.ndarray) -> float:
        """Calculate coarse-grained radius as RMS distance from COM.

        Args:
            coords: Array of atomic coordinates (N, 3).
            com: Center of mass coordinates (3,).

        Returns:
            RMS radius in Angstroms.
        """
        if len(coords) == 0:
            return 0.0
        distances_sq = np.sum((coords - com) ** 2, axis=1)
        return np.sqrt(distances_sq.mean())

    def _calculate_bounding_box(self, coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate axis-aligned bounding box.

        Args:
            coords: Array of atomic coordinates (N, 3).

        Returns:
            Tuple of (min_coords, max_coords) for bounding box.
        """
        if len(coords) == 0:
            return np.zeros(3), np.zeros(3)
        return coords.min(axis=0), coords.max(axis=0)

    def _extract_sequence(self, chain: Chain) -> str:
        """Extract amino acid sequence from chain.

        Args:
            chain: BioPython Chain object.

        Returns:
            Single-letter amino acid sequence string.
        """
        ppb = PPBuilder()
        peptides = ppb.build_peptides(chain)
        if peptides:
            return str(peptides[0].get_sequence())
        return ""


    def get_chain_ids(self) -> List[str]:
        """Get sorted list of valid chain IDs.

        Returns:
            List of chain identifiers.
        """
        return list(self.chain_data.keys())

    def get_chain_data(self, chain_id: str) -> Dict[str, Any]:
        """Get data for a specific chain.

        Args:
            chain_id: Chain identifier (may include frame suffix like 'A_f2').

        Returns:
            Dictionary containing chain data.

        Raises:
            KeyError: If chain_id not found.
        """
        if chain_id not in self.chain_data:
            raise KeyError(f"Chain {chain_id} not found. Available chains: {list(self.chain_data.keys())}")
        
        return self.chain_data[chain_id]

    def get_pdb_id(self) -> Optional[str]:
        """Get the PDB identifier.

        Returns:
            PDB ID if available, None otherwise.
        """
        return self.pdb_id

    def get_strand_ids(self) -> Optional[Dict[str, List[str]]]:
        """Extract strand IDs from mmCIF header if available.

        Returns:
            Dictionary mapping entity IDs to lists of chain IDs,
            or None if not available or not mmCIF format.
        """
        if not isinstance(self.structure, Structure):
            return None

        # This would require accessing mmCIF-specific metadata
        # Implementation depends on BioPython's mmCIF header parsing
        # For now, return None (fallback to other grouping methods)
        return None

    def convert_coords_to_nm(self, coords_angstrom: np.ndarray) -> np.ndarray:
        """Convert coordinates from Angstroms to nanometers.

        Args:
            coords_angstrom: Coordinates in Angstroms.

        Returns:
            Coordinates in nanometers (divided by 10).
        """
        return coords_angstrom / 10.0

    def convert_distance_to_angstrom(self, distance_nm: float) -> float:
        """Convert distance from nanometers to Angstroms.

        Args:
            distance_nm: Distance in nanometers.

        Returns:
            Distance in Angstroms (multiplied by 10).
        """
        return distance_nm * 10.0

    # --------------------------------------------
    # Process Multi-frame PDBs
    # --------------------------------------------
    def _concatenate_all_frames(self) -> None:
        """Concatenate all MODEL frames into a single frame with unique chain IDs."""
        
        models = list(self.structure.get_models())
        
        if len(models) <= 1:
            # No concatenation needed for single model
            return
        
        # Apply max_frames limit if specified
        if self.max_frames:
            models = models[:self.max_frames]
            if self.workspace_manager:
                self.workspace_manager.logger.info(
                    f"Concatenating first {len(models)} models (max_frames={self.max_frames})")
        
        # Use the first model as the base
        base_model = models[0]
        
        # Collect all chains from all models with frame-aware naming
        frame_chain_count = {}  # Track chains per frame for naming
        
        for frame_idx, model in enumerate(models, 1):
            frame_chain_count[frame_idx] = 0
            
            for chain in model:
                chain_id = chain.get_id()
                
                # Create frame-aware chain ID for subsequent frames
                new_chain_id = f"{chain_id}F{frame_idx}"
                
                # Clone the chain with new ID and add to base model
                if frame_idx > 1:
                    new_chain = chain.copy()
                    new_chain.id = new_chain_id
                    base_model.add(new_chain)
                    frame_chain_count[frame_idx] += 1
        
        # Remove all other models, keeping only the concatenated base model
        models_to_remove = [model.get_id() for model in models[1:]]
        for model_id in models_to_remove:
            self.structure.detach_child(model_id)
        
        if self.workspace_manager:
            total_chains = sum(len(list(model.get_chains())) for model in models)
            self.workspace_manager.logger.info(
                f"Concatenated {len(models)} frames into single frame with {total_chains} chains")
            
            # Log chain naming details
            for frame_idx, count in frame_chain_count.items():
                if frame_idx == 1:
                    original_count = len(list(models[0].get_chains()))
                    self.workspace_manager.logger.debug(f"Frame 1: {original_count} chains (original IDs)")
                else:
                    self.workspace_manager.logger.debug(f"Frame {frame_idx}: {count} chains (with _f{frame_idx} suffix)")