import unittest
import json
import math
import tempfile
from pathlib import Path
# import sys
# sys.path.append('/home/local/WIN/msang2/mankun/GitHub/ionerdss')

from ionerdss.model.pdb.core import PDBModel

def is_number(val):
    try:
        float(val)
        return True
    except (TypeError, ValueError):
        return False


def angle_close(a, b, tol=0.01):
    """
    Compare two angles a and b in radians, accounting for periodicity.
    """
    diff = math.atan2(math.sin(a - b), math.cos(a - b))
    return abs(diff) < tol


def compare_values(val1, val2, tol=0.01, path="root"):
    def is_angle_path(path):
        # Identify fields likely to contain angles
        return any(keyword in path.lower() for keyword in ["binding_angles", "theta", "angle"])

    if isinstance(val1, dict) and isinstance(val2, dict):
        if set(val1.keys()) != set(val2.keys()):
            print(f"Key mismatch at {path}: {val1.keys()} != {val2.keys()}")
            return False
        return all(compare_values(val1[k], val2[k], tol, f"{path}.{k}") for k in val1)

    elif isinstance(val1, list) and isinstance(val2, list):
        if len(val1) != len(val2):
            print(f"List length mismatch at {path}: {len(val1)} != {len(val2)}")
            return False
        return all(
            compare_values(v1, v2, tol, f"{path}[{i}]")
            for i, (v1, v2) in enumerate(zip(val1, val2))
        )

    elif is_number(val1) and is_number(val2):
        f1, f2 = float(val1), float(val2)
        if is_angle_path(path):
            if not angle_close(f1, f2, tol):
                print(f"Angle mismatch at {path}: {f1} != {f2} (wrapped, tol={tol})")
                return False
        else:
            if not math.isclose(f1, f2, abs_tol=tol):
                print(f"Value mismatch at {path}: {f1} != {f2} (tol={tol})")
                return False
        return True

    else:
        if val1 != val2:
            print(f"Exact mismatch at {path}: {val1} != {val2}")
            return False
        return True


class TestPDBModelOutput(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.save_folder = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def build_pdb_model(self, pdb_id):
        pdb_model = PDBModel(pdb_id=pdb_id, save_dir=str(self.save_folder))
        return pdb_model

    def test_model_output_8y7s(self):
        self.build_pdb_model("8y7s")

    def test_model_output_8erq(self):
        self.build_pdb_model("8erq")

    def test_model_output_5va4(self):
        self.build_pdb_model("5va4")


if __name__ == "__main__":
    unittest.main()