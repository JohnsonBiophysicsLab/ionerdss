import numpy as np
import sys
sys.path.append("..")
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

for chainPair in [("A", "H"), ("A", "L"), ("H", "L")]:
    chain1, chain2 = chainPair
    # Predict binding energy from PDB file
    binding_energy = predict_proaffinity_binding_energy(
        pdb_id="8erq",
        chains=f"{chain1},{chain2}",
        verbose=False,
        adfr_path='/home/local/WIN/msang2/ADFRsuite-1.0/bin/prepare_receptor'
    )
    print()
    print(f"Predicted binding energy between chains {chain1} and {chain2}: {binding_energy} kJ/mol")
    R = 8.314 / 1000 # kJ/(mol*K)
    temperature = 298 # K
    K = np.exp(-binding_energy / (R * temperature))
    print(f"Predicted binding constant K: {K}")
    print()