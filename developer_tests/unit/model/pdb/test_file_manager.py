"""
Unit tests for ionerdss.model.pdb.file_manager

Tests the WorkspaceManager class and its file management utilities.
"""

import unittest
import tempfile
import shutil
import logging
from pathlib import Path
from unittest.mock import patch, MagicMock

from ionerdss.model.pdb.file_manager import WorkspaceManager


class TestWorkspaceManager(unittest.TestCase):
    """Test cases for WorkspaceManager class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directory for testing
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        self.pdb_id = "1ABC"

        # Create mock datetime object
        self.mock_datetime = MagicMock()
        self.mock_datetime.isoformat.return_value = "2023-01-01T12:00:00"

    def tearDown(self):
        """Clean up test fixtures."""
        # Remove temporary directory
        if self.temp_path.exists():
            shutil.rmtree(self.temp_path)

    def test_initialization_with_pdb_id(self):
        """Test WorkspaceManager initialization with PDB ID."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Check basic attributes
        self.assertEqual(manager.workspace_path, workspace_path.resolve())
        self.assertEqual(manager.pdb_id, self.pdb_id)
        self.assertIsInstance(manager.logger, logging.Logger)

        # Check workspace directory was created
        self.assertTrue(workspace_path.exists())
        self.assertTrue(workspace_path.is_dir())

    def test_initialization_without_pdb_id(self):
        """Test WorkspaceManager initialization without PDB ID."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path)

        # Should use "unknown" as default PDB ID
        self.assertEqual(manager.pdb_id, "unknown")

    def test_workspace_structure_creation(self):
        """Test that all required directories are created."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Check all expected directories exist
        expected_dirs = [
            'logs',
            'structures',
            'structures_downloaded',
            'structures_processed',
            'outputs',
            'outputs_systems',
            'outputs_reports',
            'temp'
        ]

        for dir_name in expected_dirs:
            self.assertIn(dir_name, manager.paths)
            self.assertTrue(manager.paths[dir_name].exists())
            self.assertTrue(manager.paths[dir_name].is_dir())

    def test_workspace_structure_hierarchy(self):
        """Test that directory hierarchy is correct."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Check parent-child relationships
        self.assertEqual(
            manager.paths['structures_downloaded'].parent,
            manager.paths['structures']
        )
        self.assertEqual(
            manager.paths['structures_processed'].parent,
            manager.paths['structures']
        )
        self.assertEqual(
            manager.paths['outputs_systems'].parent,
            manager.paths['outputs']
        )
        self.assertEqual(
            manager.paths['outputs_reports'].parent,
            manager.paths['outputs']
        )

    def test_logging_setup(self):
        """Test logging configuration."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Check logger properties
        self.assertEqual(manager.logger.name, f"ionerdss.pdb.{self.pdb_id}")
        self.assertEqual(manager.logger.level, logging.INFO)

        # Check handlers
        self.assertEqual(len(manager.logger.handlers), 2)  # File + Console

        # Check log file was created
        log_file = manager.paths['logs'] / 'pipeline.log'
        self.assertTrue(log_file.exists())

    def test_get_structure_download_path_mmcif(self):
        """Test get_structure_download_path for mmCIF format."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_structure_download_path("1XYZ", "mmcif")
        expected_path = manager.paths['structures_downloaded'] / "1xyz.cif"

        self.assertEqual(path, expected_path)

    def test_get_structure_download_path_pdb(self):
        """Test get_structure_download_path for PDB format."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_structure_download_path("1XYZ", "pdb")
        expected_path = manager.paths['structures_downloaded'] / "1xyz.pdb"

        self.assertEqual(path, expected_path)

    def test_get_system_output_path_no_suffix(self):
        """Test get_system_output_path without suffix."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_system_output_path()
        expected_path = manager.paths['outputs_systems'] / \
            f"{self.pdb_id}_system.json"

        self.assertEqual(path, expected_path)

    def test_get_system_output_path_with_suffix(self):
        """Test get_system_output_path with suffix."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_system_output_path("processed")
        expected_path = manager.paths['outputs_systems'] / \
            f"{self.pdb_id}_processed_system.json"

        self.assertEqual(path, expected_path)

    def test_get_report_path(self):
        """Test get_report_path method."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_report_path("validation")
        expected_path = manager.paths['outputs_reports'] / \
            f"{self.pdb_id}_validation.txt"

        self.assertEqual(path, expected_path)

    def test_get_temp_path_no_suffix(self):
        """Test get_temp_path without suffix."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_temp_path()
        expected_path = manager.paths['temp'] / "temp"

        self.assertEqual(path, expected_path)
        self.assertTrue(path.exists())  # Should be created

    def test_get_temp_path_with_suffix(self):
        """Test get_temp_path with suffix."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        path = manager.get_temp_path("processing")
        expected_path = manager.paths['temp'] / "temp_processing"

        self.assertEqual(path, expected_path)
        self.assertTrue(path.exists())  # Should be created

    def test_move_file_to_workspace_downloaded(self):
        """Test move_file_to_workspace to downloaded directory."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create a test file
        test_file = self.temp_path / "test.pdb"
        test_file.write_text("test content")

        # Move file
        moved_path = manager.move_file_to_workspace(test_file, "downloaded")
        expected_path = manager.paths['structures_downloaded'] / "test.pdb"

        self.assertEqual(moved_path, expected_path)
        self.assertTrue(expected_path.exists())
        self.assertFalse(test_file.exists())  # Original should be gone
        self.assertEqual(expected_path.read_text(), "test content")

    def test_move_file_to_workspace_with_new_name(self):
        """Test move_file_to_workspace with new filename."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create a test file
        test_file = self.temp_path / "original.pdb"
        test_file.write_text("test content")

        # Move file with new name
        moved_path = manager.move_file_to_workspace(
            test_file, "processed", "renamed.pdb")
        expected_path = manager.paths['structures_processed'] / "renamed.pdb"

        self.assertEqual(moved_path, expected_path)
        self.assertTrue(expected_path.exists())
        self.assertFalse(test_file.exists())

    def test_move_file_to_workspace_invalid_destination(self):
        """Test move_file_to_workspace with invalid destination type."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create a test file
        test_file = self.temp_path / "test.pdb"
        test_file.write_text("test content")

        # Try to move with invalid destination
        with self.assertRaises(ValueError) as context:
            manager.move_file_to_workspace(test_file, "invalid_destination")

        self.assertIn("Invalid destination_type", str(context.exception))

    def test_move_file_to_workspace_nonexistent_file(self):
        """Test move_file_to_workspace with nonexistent source file."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Try to move nonexistent file
        nonexistent_file = self.temp_path / "nonexistent.pdb"

        with self.assertRaises(FileNotFoundError):
            manager.move_file_to_workspace(nonexistent_file, "downloaded")

    def test_copy_file_to_workspace(self):
        """Test copy_file_to_workspace method."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create a test file
        test_file = self.temp_path / "test.pdb"
        test_file.write_text("test content")

        # Copy file
        copied_path = manager.copy_file_to_workspace(test_file, "downloaded")
        expected_path = manager.paths['structures_downloaded'] / "test.pdb"

        self.assertEqual(copied_path, expected_path)
        self.assertTrue(expected_path.exists())
        self.assertTrue(test_file.exists())  # Original should still exist
        self.assertEqual(expected_path.read_text(), "test content")

    def test_cleanup_temp_files(self):
        """Test cleanup_temp_files method."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create some temp files
        temp_file1 = manager.paths['temp'] / "temp1.txt"
        temp_file2 = manager.paths['temp'] / "temp2.txt"
        temp_file1.write_text("temp content 1")
        temp_file2.write_text("temp content 2")

        # Create subdirectory with file
        temp_subdir = manager.paths['temp'] / "subdir"
        temp_subdir.mkdir()
        temp_file3 = temp_subdir / "temp3.txt"
        temp_file3.write_text("temp content 3")

        # Verify files exist
        self.assertTrue(temp_file1.exists())
        self.assertTrue(temp_file2.exists())
        self.assertTrue(temp_file3.exists())

        # Cleanup
        manager.cleanup_temp_files()

        # Verify temp directory is empty but still exists
        self.assertTrue(manager.paths['temp'].exists())
        self.assertEqual(len(list(manager.paths['temp'].iterdir())), 0)

    def test_cleanup_temp_files_error_handling(self):
        """Test cleanup_temp_files error handling."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Mock shutil.rmtree to raise an exception
        with patch('ionerdss.model.pdb.file_manager.shutil.rmtree') as mock_rmtree:
            mock_rmtree.side_effect = PermissionError("Access denied")

            # Should not raise exception, just log warning
            manager.cleanup_temp_files()

            # Verify rmtree was called
            mock_rmtree.assert_called_once()

    def test_generate_summary_report(self):
        """Test generate_summary_report method."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create some test files
        test_file1 = manager.paths['structures_downloaded'] / "test1.pdb"
        test_file2 = manager.paths['outputs_systems'] / "test2.json"
        test_file1.write_text("content1")
        test_file2.write_text("content2")

        # Generate report
        report_path = manager.generate_summary_report()
        expected_path = manager.paths['outputs_reports'] / \
            f"{self.pdb_id}_summary.txt"

        self.assertEqual(report_path, expected_path)
        self.assertTrue(report_path.exists())

        # Check report content
        report_content = report_path.read_text()
        self.assertIn("PDB Pipeline Summary Report", report_content)
        self.assertIn(f"PDB ID: {self.pdb_id}", report_content)
        self.assertIn("test1.pdb", report_content)
        self.assertIn("test2.json", report_content)

    def test_get_workspace_info(self):
        """Test get_workspace_info method."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Create some test files
        test_file1 = manager.paths['structures_downloaded'] / "test1.pdb"
        test_file2 = manager.paths['structures_downloaded'] / "test2.cif"
        test_file1.write_text("content1")
        test_file2.write_text("content2")

        info = manager.get_workspace_info()

        # Check info structure
        self.assertIn('workspace_path', info)
        self.assertIn('pdb_id', info)
        self.assertIn('paths', info)
        self.assertIn('files', info)

        # Check specific values
        self.assertEqual(info['pdb_id'], self.pdb_id)
        self.assertEqual(info['workspace_path'], str(workspace_path.resolve()))

        # Check file counts
        self.assertEqual(info['files']['structures_downloaded']['count'], 2)
        self.assertIn('test1.pdb', info['files']
                      ['structures_downloaded']['files'])
        self.assertIn('test2.cif', info['files']
                      ['structures_downloaded']['files'])

    def test_context_manager_success(self):
        """Test context manager with successful execution."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            with WorkspaceManager(workspace_path, self.pdb_id) as manager:
                # Do some work
                test_file = manager.paths['outputs'] / "test.txt"
                test_file.write_text("test")

        # Check that summary report was generated
        summary_path = manager.paths['outputs_reports'] / \
            f"{self.pdb_id}_summary.txt"
        self.assertTrue(summary_path.exists())

    def test_context_manager_with_exception(self):
        """Test context manager with exception."""
        workspace_path = self.temp_path / "test_workspace"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            try:
                with WorkspaceManager(workspace_path, self.pdb_id) as manager:
                    raise ValueError("Test exception")
            except ValueError:
                pass  # Expected

        # Check that summary report was still generated
        summary_path = manager.paths['outputs_reports'] / \
            f"{self.pdb_id}_summary.txt"
        self.assertTrue(summary_path.exists())

    def test_string_path_input(self):
        """Test that string paths are properly converted to Path objects."""
        workspace_path_str = str(self.temp_path / "test_workspace")

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path_str, self.pdb_id)

        # Should convert string to Path and resolve it
        self.assertIsInstance(manager.workspace_path, Path)
        self.assertEqual(manager.workspace_path,
                         Path(workspace_path_str).resolve())

    def test_existing_workspace_directory(self):
        """Test initialization with existing workspace directory."""
        workspace_path = self.temp_path / "existing_workspace"
        workspace_path.mkdir()

        # Create some existing content
        existing_file = workspace_path / "existing.txt"
        existing_file.write_text("existing content")

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager = WorkspaceManager(workspace_path, self.pdb_id)

        # Should not overwrite existing content
        self.assertTrue(existing_file.exists())
        self.assertEqual(existing_file.read_text(), "existing content")

        # Should still create required subdirectories
        self.assertTrue(manager.paths['logs'].exists())
        self.assertTrue(manager.paths['structures'].exists())

    def test_logger_name_uniqueness(self):
        """Test that different PDB IDs get different logger names."""
        workspace_path1 = self.temp_path / "workspace1"
        workspace_path2 = self.temp_path / "workspace2"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            manager1 = WorkspaceManager(workspace_path1, "1ABC")
            manager2 = WorkspaceManager(workspace_path2, "2XYZ")

        self.assertEqual(manager1.logger.name, "ionerdss.pdb.1ABC")
        self.assertEqual(manager2.logger.name, "ionerdss.pdb.2XYZ")
        self.assertNotEqual(manager1.logger.name, manager2.logger.name)


class TestWorkspaceManagerIntegration(unittest.TestCase):
    """Integration tests for WorkspaceManager."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

        # Create mock datetime object
        self.mock_datetime = MagicMock()
        self.mock_datetime.isoformat.return_value = "2023-01-01T12:00:00"

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_path.exists():
            shutil.rmtree(self.temp_path)

    def test_complete_workflow(self):
        """Test a complete file management workflow."""
        workspace_path = self.temp_path / "workflow_test"
        pdb_id = "1XYZ"

        with patch('ionerdss.model.pdb.file_manager.datetime') as mock_dt:
            mock_dt.now.return_value = self.mock_datetime

            with WorkspaceManager(workspace_path, pdb_id) as manager:
                # Simulate downloading a structure
                downloaded_file = manager.get_structure_download_path(
                    pdb_id, "mmcif")
                downloaded_file.write_text("structure content")

                # Simulate processing and creating output
                system_file = manager.get_system_output_path("final")
                system_file.write_text('{"system": "data"}')

                # Create temporary files
                temp_dir = manager.get_temp_path("processing")
                temp_file = temp_dir / "intermediate.tmp"
                temp_file.write_text("temporary data")

                # Clean up temp files
                manager.cleanup_temp_files()

                # Verify temp files are gone but others remain
                self.assertFalse(temp_file.exists())
                self.assertTrue(downloaded_file.exists())
                self.assertTrue(system_file.exists())

        # After context exit, summary should be generated
        summary_file = manager.paths['outputs_reports'] / \
            f"{pdb_id}_summary.txt"
        self.assertTrue(summary_file.exists())

        # Check summary content
        summary_content = summary_file.read_text()
        self.assertIn(f"{pdb_id.lower()}.cif", summary_content)
        self.assertIn(f"{pdb_id}_final_system.json", summary_content)


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)
