import unittest
import json
import math
import tempfile
from pathlib import Path

from ionerdss import PDBModel, ParseComplexes


def is_number(val):
    try:
        float(val)
        return True
    except (TypeError, ValueError):
        return False


def compare_values(val1, val2, tol=0.01, path="root"):
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
        pdb_model.coarse_grain(
            distance_cutoff=0.35,
            residue_cutoff=3,
            show_coarse_grained_structure=False,
            save_pymol_script=False,
            standard_output=False,
        )
        pdb_model.regularize_homologous_chains(
            dist_threshold_intra=3.5,
            dist_threshold_inter=3.5,
            angle_threshold=25.0,
            show_coarse_grained_structure=False,
            save_pymol_script=False,
            standard_output=False,
        )
        return pdb_model

    def run_model_test(self, pdb_id, tol=0.01):
        # Get the path to the data directory relative to this test file
        test_dir = Path(__file__).parent.parent
        expected_path = test_dir / "data" / f"{pdb_id}_model.json"
        actual_path = self.save_folder / f"{pdb_id}_model.json"

        pdb_model = self.build_pdb_model(pdb_id)

        with open(expected_path, "r") as f_expected:
            expected_data = json.load(f_expected)

        with open(actual_path, "r") as f_actual:
            actual_data = json.load(f_actual)

        self.assertTrue(
            compare_values(expected_data, actual_data, tol=tol),
            f"The actual model output for {pdb_id} does not match the expected output within the tolerance.",
        )

    def test_model_output_8y7s(self):
        self.run_model_test("8y7s", tol=0.02)

    def test_model_output_8erq(self):
        self.run_model_test("8erq")

    def test_model_output_5va4(self):
        self.run_model_test("5va4")

    def test_parse_complexes_print_8erq(self):
        pdb_model = self.build_pdb_model("8erq")
        complex_list, complex_reaction_system = ParseComplexes(pdb_model)

        complex_list_length = len(complex_list)
        self.assertEqual(complex_list_length, 10, "Complex list length did not match expected.")

        complex_reaction_system_length = len(complex_reaction_system.reactions)
        self.assertEqual(complex_reaction_system_length, 24, "Complex reaction system length did not match expected.")

    def test_parse_complexes_print_5va4(self):
        pdb_model = self.build_pdb_model("5va4")
        complex_list, complex_reaction_system = ParseComplexes(pdb_model)

        complex_list_length = len(complex_list)
        self.assertEqual(complex_list_length, 4, "Complex list length did not match expected.")

        complex_reaction_system_length = len(complex_reaction_system.reactions)
        self.assertEqual(complex_reaction_system_length, 6, "Complex reaction system length did not match expected.")

    def test_parse_complexes_print_8y7s(self):
        pdb_model = self.build_pdb_model("8y7s")
        complex_list, complex_reaction_system = ParseComplexes(pdb_model)

        complex_list_length = len(complex_list)
        self.assertEqual(complex_list_length, 25, "Complex list length did not match expected.")

        complex_reaction_system_length = len(complex_reaction_system.reactions)
        self.assertEqual(complex_reaction_system_length, 114, "Complex reaction system length did not match expected.")

    # def test_model_output_7uhy(self):
    #     self.run_model_test("7uhy", tol=1)


if __name__ == "__main__":
    unittest.main()