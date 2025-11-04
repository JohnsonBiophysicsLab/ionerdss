"""Test script for ProAffinity-GNN integration."""

import sys
import os

# Add ionerdss to path
sys.path.append('../../')

from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

def test_proaffinity_integration():
    """Test the complete ProAffinity-GNN integration."""
    
    print("=" * 60)
    print("Testing ProAffinity-GNN Integration")
    print("=" * 60)
    
    # Test parameters
    test_pdb = "1PPE"
    test_chains = "E,I"
    
    print(f"\nTest case: {test_pdb} with chains {test_chains}")
    print("-" * 60)
    
    # Test the complete pipeline
    print("\n✓ Testing complete pipeline (PDB download → PDBQT → Prediction)")
    
    try:
        result = predict_proaffinity_binding_energy(
            pdb_id=test_pdb,
            chains=test_chains,
            verbose=True
        )

        print("-" * 60)
        print(f"\n✓ Prediction successful!")
        print(f"Binding energy: {result:.4f} kJ/mol")
        print("-" * 60)
        
        # Validate result
        if result is not None and not isinstance(result, type(None)):
            print("\n✓ Result is valid (not None)")
            return True
        else:
            print("\n✗ Result is None or invalid")
            return False
            
    except Exception as e:
        print(f"\n✗ Error during test: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("ProAffinity-GNN Integration Test Suite")
    print("=" * 60 + "\n")
    
    success = test_proaffinity_integration()
    
    if success:
        print("\n" + "=" * 60)
        print("Integration test PASSED ✓")
        print("=" * 60)
    else:
        print("\n" + "=" * 60)
        print("Integration test FAILED ✗")
        print("=" * 60)
        sys.exit(1)
