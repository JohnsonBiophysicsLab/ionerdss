#!/usr/bin/env python3
"""
Test script for biological assembly download functionality.
This script demonstrates the new bioassembly download feature.
"""

from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

def test_bioassembly1_default():
    """Test default bioassembly1 download."""
    print("="*60)
    print("Test 1: Default bioassembly1 download (5L93)")
    print("="*60)
    
    # Use default hyperparameters (pdb_file_format="bioassembly1")
    params = PDBModelHyperparameters()
    print(f"pdb_file_format: {params.pdb_file_format}")
    
    model = PDBModelBuilder(source="5L93", hyperparams=params)
    system = model.build_system(
        workspace_path="test_5L93_bioassembly1",
        hyperparams=params
    )
    
    print(f"✓ Successfully downloaded and parsed bioassembly1 for 5L93")
    print(f"  Chains: {system.get_summary()['num_chains']}")
    print()


def test_bioassembly2():
    """Test bioassembly2 download."""
    print("="*60)
    print("Test 2: bioassembly2 download (1A0R)")
    print("="*60)
    
    params = PDBModelHyperparameters(pdb_file_format="bioassembly2")
    print(f"pdb_file_format: {params.pdb_file_format}")
    
    model = PDBModelBuilder(source="1A0R", hyperparams=params)
    system = model.build_system(
        workspace_path="test_1A0R_bioassembly2",
        hyperparams=params
    )
    
    print(f"✓ Successfully downloaded and parsed bioassembly2 for 1A0R")
    print(f"  Chains: {system.get_summary()['num_chains']}")
    print()


def test_case_insensitive():
    """Test case-insensitive format matching."""
    print("="*60)
    print("Test 3: Case-insensitive format (BioAssembly1)")
    print("="*60)
    
    params = PDBModelHyperparameters(pdb_file_format="BioAssembly1")
    print(f"pdb_file_format: {params.pdb_file_format}")
    
    model = PDBModelBuilder(source="6BNO", hyperparams=params)
    system = model.build_system(
        workspace_path="test_6BNO_bioassembly1",
        hyperparams=params
    )
    
    print(f"✓ Case-insensitive format worked correctly")
    print(f"  Chains: {system.get_summary()['num_chains']}")
    print()


def test_standard_cif():
    """Test backward compatibility with standard cif format."""
    print("="*60)
    print("Test 4: Standard CIF format (backward compatibility)")
    print("="*60)
    
    params = PDBModelHyperparameters(pdb_file_format="cif")
    print(f"pdb_file_format: {params.pdb_file_format}")
    
    model = PDBModelBuilder(source="1ABC", hyperparams=params)
    system = model.build_system(
        workspace_path="test_1ABC_cif",
        hyperparams=params
    )
    
    print(f"✓ Standard CIF format still works")
    print(f"  Chains: {system.get_summary()['num_chains']}")
    print()


def test_nonexistent_assembly():
    """Test error handling for non-existent assembly."""
    print("="*60)
    print("Test 5: Non-existent assembly (error handling)")
    print("="*60)
    
    params = PDBModelHyperparameters(pdb_file_format="bioassembly99")
    print(f"pdb_file_format: {params.pdb_file_format}")
    
    try:
        model = PDBModelBuilder(source="5L93", hyperparams=params)
        system = model.build_system(
            workspace_path="test_5L93_bioassembly99",
            hyperparams=params
        )
        print("✗ Should have raised an error!")
    except ValueError as e:
        print(f"✓ Correctly raised error: {str(e)[:100]}...")
    print()


if __name__ == "__main__":
    print("\\nTesting Biological Assembly Download Functionality\\n")
    
    # Run tests
    try:
        test_bioassembly1_default()
    except Exception as e:
        print(f"✗ Test 1 failed: {e}\\n")
    
    try:
        test_bioassembly2()
    except Exception as e:
        print(f"✗ Test 2 failed: {e}\\n")
    
    try:
        test_case_insensitive()
    except Exception as e:
        print(f"✗ Test 3 failed: {e}\\n")
    
    try:
        test_standard_cif()
    except Exception as e:
        print(f"✗ Test 4 failed: {e}\\n")
    
    try:
        test_nonexistent_assembly()
    except Exception as e:
        print(f"✗ Test 5 failed: {e}\\n")
    
    print("\\n" + "="*60)
    print("Testing complete!")
    print("="*60)
