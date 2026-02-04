"""
ionerdss.model.pdb.file_manager

File management utilities for PDB pipeline workspace organization.

This module provides centralized file management for the PDB to NERDSS pipeline,
ensuring all downloaded structures, generated files, and outputs are properly
organized within the workspace directory. It handles workspace creation, file 
operations, logging, and cleanup.

## Table of Contents

- [Workspace Structure](#workspace-structure)
- [Basic Usage](#basic-usage)
- [Understanding Output Files](#understanding-output-files)
- [File Organization](#file-organization)
- [Logging and Reports](#logging-and-reports)
- [Common Workflows](#common-workflows)

## Workspace Structure

When you run the ionerdss pipeline, a structured workspace is automatically created:

```
workspace_1ABC/
├── logs/
│   └── pipeline.log              # Complete pipeline execution log
├── structures/
│   ├── downloaded/
│   │   └── 1abc.cif             # Original structure from PDB
│   └── processed/
│       └── 1abc_cleaned.pdb     # Processed/cleaned structures
├── outputs/
│   ├── systems/
│   │   ├── 1ABC_system.json     # Final ionerdss system definition
│   │   └── 1ABC_final_system.json
│   └── reports/
│       ├── 1ABC_summary.txt     # Workspace summary report
│       └── 1ABC_validation.txt  # Validation reports
└── temp/                        # Temporary files (cleaned up automatically)
```

## Basic Usage

### Context Manager (Recommended)

```python
from ionerdss.model.pdb.file_manager import WorkspaceManager

# Automatic cleanup and report generation
with WorkspaceManager("/path/to/workspace", "1ABC") as manager:
    # Your pipeline code here
    structure_path = manager.get_structure_download_path("1ABC", "mmcif")
    # ... process structure ...
    output_path = manager.get_system_output_path("final")
    # ... save results ...
# Summary report automatically generated on exit
```

### Manual Management

```python
manager = WorkspaceManager("/path/to/workspace", "1ABC")

# Get paths for different file types
download_path = manager.get_structure_download_path("1ABC", "mmcif")
output_path = manager.get_system_output_path()
report_path = manager.get_report_path("validation")

# File operations
manager.copy_file_to_workspace("external_file.pdb", "downloaded")
manager.cleanup_temp_files()

# Generate final report
manager.generate_summary_report()
```

## Understanding Output Files

### System Files (`outputs/systems/`)

**Purpose:** Contains the final ionerdss system definitions in JSON format.

**File Naming:**
- `{PDB_ID}_system.json` - Basic system output
- `{PDB_ID}_{suffix}_system.json` - Processed variants (e.g., "final", "validated")

**Content:** Complete molecular system definition including:
- Molecule types and instances
- Interface definitions
- Binding parameters
- Geometric properties

**Example:**
```json
{
  "metadata": {
    "pdb_id": "1ABC",
    "workspace_path": "/path/to/workspace",
    "timestamp": "2023-01-01T12:00:00"
  },
  "molecule_types": [...],
  "molecule_instances": [...],
  "interface_types": [...],
  "interface_instances": [...]
}
```

### Structure Files (`structures/`)

**Downloaded (`structures/downloaded/`):**
- Original PDB/mmCIF files from the Protein Data Bank
- Unmodified reference structures
- Named as `{pdb_id}.{extension}` (e.g., `1abc.cif`, `1abc.pdb`)

**Processed (`structures/processed/`):**
- Cleaned or modified structure files
- Intermediate processing results
- May include chain renaming, residue filtering, etc.

### Reports (`outputs/reports/`)

**Summary Report (`{PDB_ID}_summary.txt`):**
```
PDB Pipeline Summary Report
==================================================

PDB ID: 1ABC
Workspace: /path/to/workspace
Generated: 2023-01-01T12:00:00

STRUCTURES_DOWNLOADED:
  - 1abc.cif (245760 bytes)

OUTPUTS_SYSTEMS:
  - 1ABC_system.json (15420 bytes)
  - 1ABC_final_system.json (15890 bytes)

OUTPUTS_REPORTS:
  - 1ABC_summary.txt (2048 bytes)
```

**Validation Reports:** Custom reports generated during pipeline execution containing:
- Structure quality assessments
- Interface detection statistics
- Processing warnings and errors

### Log Files (`logs/`)

**Pipeline Log (`pipeline.log`):**
```
2023-01-01 12:00:00 - ionerdss.pdb.1ABC - INFO - === PDB Pipeline Started for 1ABC ===
2023-01-01 12:00:01 - ionerdss.pdb.1ABC - INFO - Workspace: /path/to/workspace
2023-01-01 12:00:02 - ionerdss.pdb.1ABC - INFO - Downloaded structure: 1abc.cif
2023-01-01 12:00:05 - ionerdss.pdb.1ABC - INFO - Detected 12 interfaces
2023-01-01 12:00:08 - ionerdss.pdb.1ABC - INFO - Generated system: 1ABC_system.json
2023-01-01 12:00:10 - ionerdss.pdb.1ABC - INFO - === PDB Pipeline Completed for 1ABC ===
```

## File Organization

### File Type Categories

| Category | Purpose | Persistence | Examples |
|----------|---------|-------------|----------|
| **Downloaded** | Original PDB data | Permanent | `1abc.cif`, `2xyz.pdb` |
| **Processed** | Intermediate results | Permanent | `cleaned.pdb`, `filtered.cif` |
| **Systems** | Final outputs | Permanent | `1ABC_system.json` |
| **Reports** | Analysis results | Permanent | `summary.txt`, `validation.txt` |
| **Logs** | Execution traces | Permanent | `pipeline.log` |
| **Temporary** | Working files | Cleaned up (*) | `temp_*.tmp`, processing files |

(*) Temporary files are cleaned up after the pipeline completes.

### Path Generation

```python
# Structure files
mmcif_path = manager.get_structure_download_path("1ABC", "mmcif")
# → workspace/structures/downloaded/1abc.cif

pdb_path = manager.get_structure_download_path("1ABC", "pdb") 
# → workspace/structures/downloaded/1abc.pdb

# System outputs
system_path = manager.get_system_output_path()
# → workspace/outputs/systems/1ABC_system.json

system_path = manager.get_system_output_path("final")
# → workspace/outputs/systems/1ABC_final_system.json

# Reports
report_path = manager.get_report_path("validation")
# → workspace/outputs/reports/1ABC_validation.txt

# Temporary files
temp_path = manager.get_temp_path("processing")
# → workspace/temp/temp_processing/
```

## Logging and Reports

### Understanding Log Levels

**INFO:** Normal pipeline progress
```
INFO - Downloaded structure: 1abc.cif
INFO - Detected 12 interfaces between 4 chains
INFO - Generated system with 8 molecule instances
```

**WARNING:** Issues that don't stop processing
```
WARNING - Chain X has missing residues
WARNING - Interface detection threshold may be too strict
```

**ERROR:** Problems that halt processing
```
ERROR - Failed to download structure 1ABC
ERROR - Invalid PDB format detected
```

### Summary Report Interpretation

**File Counts:** Number of files in each category
- High counts in `temp/` suggest cleanup issues
- Missing files in `outputs/` indicate processing failures
- Multiple versions in `systems/` show processing iterations

**File Sizes:** Indicate content complexity
- Large structure files (>1MB) suggest complex assemblies
- Small system files (<10KB) may indicate simple structures
- Zero-byte files indicate processing failures

## Common Workflows

### Basic Pipeline Execution

```python
with WorkspaceManager("/my/workspace", "1ABC") as manager:
    # 1. Download structure
    structure_file = manager.get_structure_download_path("1ABC", "mmcif")
    download_structure("1ABC", structure_file)
    
    # 2. Process structure
    processed_file = manager.get_system_output_path("processed")
    system = process_structure(structure_file)
    save_system(system, processed_file)
    
    # 3. Generate final output
    final_file = manager.get_system_output_path("final")
    final_system = finalize_system(system)
    save_system(final_system, final_file)
    
    # 4. Create validation report
    report_file = manager.get_report_path("validation")
    generate_validation_report(final_system, report_file)

# Summary report automatically created at workspace/outputs/reports/1ABC_summary.txt
```

### File Import/Export

```python
manager = WorkspaceManager("/workspace", "1ABC")

# Import external files
external_structure = "/external/path/structure.pdb"
imported_path = manager.copy_file_to_workspace(
    external_structure, 
    "downloaded", 
    "1abc_external.pdb"
)


"""

import shutil
from pathlib import Path
from typing import Optional, Union
import logging
from datetime import datetime


class WorkspaceManager:
    """Manages file organization within a workspace directory.

    Creates and maintains a structured workspace for PDB processing with
    subdirectories for different file types and comprehensive logging.

    Workspace Structure:
        workspace/
        ├── logs/
        │   └── pipeline.log
        ├── structures/
        │   ├── downloaded/
        │   └── processed/
        ├── outputs/
        │   ├── systems/
        │   └── reports/
        └── temp/

    Attributes:
        workspace_path: Path to the main workspace directory.
        pdb_id: PDB identifier for this workspace.
        logger: Logger instance for this workspace.
        paths: Dictionary of important subdirectory paths.
    """

    def __init__(self, workspace_path: Union[str, Path], pdb_id: Optional[str] = None):
        """Initialize workspace manager.

        Args:
            workspace_path: Path to workspace directory.
            pdb_id: PDB identifier (used for naming and logging).
        """
        self.workspace_path = Path(workspace_path).resolve()
        self.pdb_id = pdb_id or "unknown"

        # Create workspace structure
        self.paths = self._create_workspace_structure()

        # Set up logging
        self.logger = self._setup_logging()

        self.logger.info("Initialized workspace for %s at %s",
                         self.pdb_id, self.workspace_path)

    def _create_workspace_structure(self) -> dict:
        """Create the workspace directory structure.

        Returns:
            Dictionary mapping directory names to Path objects.
        """
        # Main workspace directory
        self.workspace_path.mkdir(parents=True, exist_ok=True)

        # Subdirectories
        subdirs = {
            'logs': self.workspace_path / 'logs',
            'structures': self.workspace_path / 'structures',
            'structures_downloaded': self.workspace_path / 'structures' / 'downloaded',
            'structures_processed': self.workspace_path / 'structures' / 'processed',
            'outputs': self.workspace_path / 'outputs',
            'outputs_systems': self.workspace_path / 'outputs' / 'systems',
            'outputs_reports': self.workspace_path / 'outputs' / 'reports',
            'temp': self.workspace_path / 'temp'
        }

        # Create all subdirectories
        for path in subdirs.values():
            path.mkdir(parents=True, exist_ok=True)

        return subdirs

    def _setup_logging(self) -> logging.Logger:
        """Set up logging for the workspace.

        Returns:
            Configured logger instance.
        """
        # Create logger
        logger_name = f"ionerdss.pdb.{self.pdb_id}"
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.INFO)
        logger.propagate = False

        # Remove existing handlers to avoid duplicates
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        # File handler for workspace log
        log_file = self.paths['logs'] / 'pipeline.log'
        file_handler = logging.FileHandler(
            log_file, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.INFO)

        # Console handler - WARNING level to reduce output verbosity
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.WARNING)

        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add handlers
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        # Log workspace creation
        logger.info("=== PDB Pipeline Started for %s ===", self.pdb_id)
        logger.info("Workspace: %s", self.workspace_path)
        logger.info("Timestamp: %s", datetime.now().isoformat())

        return logger

    def get_structure_download_path(self, pdb_id: str, file_format: str = 'mmcif') -> Path:
        """Get path for downloading a structure file.

        Args:
            pdb_id: PDB identifier.
            file_format: File format ('pdb' or 'mmcif').

        Returns:
            Path where the structure should be downloaded.
        """
        if file_format.lower() == 'mmcif':
            filename = f"{pdb_id.lower()}.cif"
        else:
            filename = f"{pdb_id.lower()}.pdb"

        return self.paths['structures_downloaded'] / filename

    def get_system_output_path(self, suffix: str = "") -> Path:
        """Get path for system JSON output.

        Args:
            suffix: Optional suffix for filename.

        Returns:
            Path for system JSON file.
        """
        if suffix:
            filename = f"{self.pdb_id}_{suffix}_system.json"
        else:
            filename = f"{self.pdb_id}_system.json"

        return self.paths['outputs_systems'] / filename

    def get_report_path(self, report_type: str) -> Path:
        """Get path for a report file.

        Args:
            report_type: Type of report (e.g., 'summary', 'validation').

        Returns:
            Path for report file.
        """
        filename = f"{self.pdb_id}_{report_type}.txt"
        return self.paths['outputs_reports'] / filename

    def get_temp_path(self, suffix: str = "") -> Path:
        """Get path for temporary files.

        Args:
            suffix: Optional suffix for temp directory.

        Returns:
            Path for temporary files.
        """
        if suffix:
            temp_dir = self.paths['temp'] / f"temp_{suffix}"
        else:
            temp_dir = self.paths['temp'] / "temp"

        temp_dir.mkdir(exist_ok=True)
        return temp_dir

    def move_file_to_workspace(self, source_path: Union[str, Path],
                               destination_type: str, new_name: Optional[str] = None) -> Path:
        """Move a file into the workspace.

        Args:
            source_path: Path to source file.
            destination_type: Type of destination ('downloaded', 'processed', 'output').
            new_name: Optional new filename.

        Returns:
            Path to moved file.

        Raises:
            ValueError: If destination_type is invalid.
            FileNotFoundError: If source file doesn't exist.
        """
        source_path = Path(source_path)

        if not source_path.exists():
            raise FileNotFoundError(f"Source file not found: {source_path}")

        # Determine destination directory
        if destination_type == 'downloaded':
            dest_dir = self.paths['structures_downloaded']
        elif destination_type == 'processed':
            dest_dir = self.paths['structures_processed']
        elif destination_type == 'output':
            dest_dir = self.paths['outputs']
        else:
            raise ValueError(f"Invalid destination_type: {destination_type}")

        # Determine filename
        if new_name:
            dest_path = dest_dir / new_name
        else:
            dest_path = dest_dir / source_path.name

        # Move file
        shutil.move(str(source_path), str(dest_path))
        self.logger.info("Moved file: %s -> %s", source_path, dest_path)

        return dest_path

    def copy_file_to_workspace(self, source_path: Union[str, Path],
                               destination_type: str, new_name: Optional[str] = None) -> Path:
        """Copy a file into the workspace.

        Args:
            source_path: Path to source file.
            destination_type: Type of destination ('downloaded', 'processed', 'output').
            new_name: Optional new filename.

        Returns:
            Path to copied file.
        """
        source_path = Path(source_path)

        if not source_path.exists():
            raise FileNotFoundError(f"Source file not found: {source_path}")

        # Determine destination directory
        if destination_type == 'downloaded':
            dest_dir = self.paths['structures_downloaded']
        elif destination_type == 'processed':
            dest_dir = self.paths['structures_processed']
        elif destination_type == 'output':
            dest_dir = self.paths['outputs']
        else:
            raise ValueError(f"Invalid destination_type: {destination_type}")

        # Determine filename
        if new_name:
            dest_path = dest_dir / new_name
        else:
            dest_path = dest_dir / source_path.name

        # Copy file
        shutil.copy2(str(source_path), str(dest_path))
        self.logger.info("Copied file: %s -> %s", source_path, dest_path)

        return dest_path

    def cleanup_temp_files(self):
        """Clean up temporary files and directories."""
        temp_dir = self.paths['temp']
        if temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
                temp_dir.mkdir(exist_ok=True)  # Recreate empty temp dir
                self.logger.info("Cleaned up temporary files")
            except Exception as e:
                self.logger.warning("Failed to clean up temp files: %s", e)

    def generate_summary_report(self) -> Path:
        """Generate a summary report of the workspace contents.

        Returns:
            Path to the generated summary report.
        """
        report_path = self.get_report_path('summary')

        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("PDB Pipeline Summary Report\n")
            f.write(f"{'=' * 50}\n\n")
            f.write(f"PDB ID: {self.pdb_id}\n")
            f.write(f"Workspace: {self.workspace_path}\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n\n")

            # List files in each directory
            for dir_name, dir_path in self.paths.items():
                if dir_path.exists():
                    files = list(dir_path.glob('*'))
                    f.write(f"{dir_name.upper()}:\n")
                    if files:
                        for file_path in sorted(files):
                            if file_path.is_file():
                                size = file_path.stat().st_size
                                f.write(
                                    f"  - {file_path.name} ({size} bytes)\n")
                    else:
                        f.write("  (empty)\n")
                    f.write("\n")

        self.logger.info("Generated summary report: %s", report_path)
        return report_path

    def get_workspace_info(self) -> dict:
        """Get information about the workspace.

        Returns:
            Dictionary with workspace information.
        """
        info = {
            'workspace_path': str(self.workspace_path),
            'pdb_id': self.pdb_id,
            'paths': {name: str(path) for name, path in self.paths.items()},
            'files': {}
        }

        # Count files in each directory
        for dir_name, dir_path in self.paths.items():
            if dir_path.exists():
                files = [f for f in dir_path.glob('*') if f.is_file()]
                info['files'][dir_name] = {
                    'count': len(files),
                    'files': [f.name for f in files]
                }
            else:
                info['files'][dir_name] = {'count': 0, 'files': []}

        return info

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup and final logging."""
        if exc_type is not None:
            self.logger.error("Pipeline failed with %s: %s",
                              exc_type.__name__, exc_val)
        else:
            self.logger.info("Pipeline completed successfully")

        # Generate final summary
        self.generate_summary_report()

        # Log final statistics
        info = self.get_workspace_info()
        self.logger.info("Final workspace contents:")
        for dir_name, file_info in info['files'].items():
            self.logger.info("  %s: %d files", dir_name, file_info['count'])

        self.logger.info("=== PDB Pipeline Completed for %s ===", self.pdb_id)
