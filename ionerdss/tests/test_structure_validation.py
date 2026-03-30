from pathlib import Path
import tempfile
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model import pdb
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.pdb.structure_validation import (
    align_structure_to_design,
    get_designed_structure_coordinates,
    get_structure_validation_counts,
)


def _make_instance(name: str, molecule_type: MoleculeType, com, interface_count: int) -> MoleculeInstance:
    return MoleculeInstance(
        name=name,
        molecule_type=molecule_type,
        com=np.asarray(com, dtype=float),
        norm=np.array([0.0, 0.0, 1.0]),
        ref1=np.array([1.0, 0.0, 0.0]),
        ref2=np.array([0.0, 1.0, 0.0]),
        interfaces_neighbors_map={f"iface_{i}": None for i in range(interface_count)},
    )


def test_validation_counts_use_one_copy_per_molecule_type():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    type_b = MoleculeType(name="B")
    system.molecule_types.add(type_a)
    system.molecule_types.add(type_b)

    system.molecule_instances.add(_make_instance("a_sparse", type_a, [0.0, 0.0, 0.0], interface_count=1))
    system.molecule_instances.add(_make_instance("a_dense", type_a, [1.0, 0.0, 0.0], interface_count=3))
    system.molecule_instances.add(_make_instance("b_only", type_b, [0.0, 2.0, 0.0], interface_count=2))

    counts = get_structure_validation_counts(system)
    coords = get_designed_structure_coordinates(system)

    assert counts == {"A": 1, "B": 1}
    assert coords == {"A": (1.0, 0.0, 0.0), "B": (0.0, 2.0, 0.0)}


def test_align_structure_to_design_recovers_rigid_transform():
    designed = {
        "A": [0.0, 0.0, 0.0],
        "B": [1.0, 0.0, 0.0],
        "C": [0.0, 1.0, 0.0],
    }

    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([5.0, -3.0, 2.0])

    designed_xyz = np.asarray([designed[key] for key in sorted(designed)], dtype=float)
    observed_xyz = (rotation @ designed_xyz.T).T + translation
    observed = {
        key: observed_xyz[idx].tolist()
        for idx, key in enumerate(sorted(designed))
    }

    result = align_structure_to_design(designed, observed)

    assert result.labels == ("A", "B", "C")
    assert result.rmsd < 1e-10
    assert np.allclose(result.aligned_observed_coordinates, result.designed_coordinates)


def test_validation_module_exposes_prepare_and_compare():
    assert hasattr(pdb, "validation")
    assert callable(pdb.validation.prepare)
    assert callable(pdb.validation.compare)
    assert callable(pdb.validation.setup_simulation)
    assert callable(pdb.validation.align_structure)


def test_setup_simulation_writes_titration_file():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    system.molecule_types.add(type_a)
    system.molecule_instances.add(_make_instance("a_only", type_a, [0.0, 0.0, 0.0], interface_count=0))

    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            import os
            os.chdir(tmpdir)
            artifacts = pdb.validation.setup_simulation(
                system,
                initial_molecule_count=1,
                titration_on_rate=2.0e-5,
            )
        finally:
            os.chdir(old_cwd)

    assert artifacts.molecule_counts == {"A": 1}
    assert "parms_titrate" in artifacts.nerdss_files
    assert artifacts.nerdss_files["parms_titrate"].name == "parms_titrate.inp"


def test_titration_sites_follow_mol_file_and_skip_com_ref():
    system = System(workspace_path=".")

    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            import os
            os.chdir(tmpdir)
            exporter = NERDSSExporter(system)
            mol_path = exporter.output_dir / "A.mol"
            mol_path.write_text(
                "\n".join(
                    [
                        "Name = A",
                        "",
                        "# Coordinates",
                        "COM 0.0 0.0 0.0",
                        "REF 1.0 0.0 0.0",
                        "al1 0.0 1.0 0.0",
                        "ah1 0.0 0.0 1.0",
                    ]
                ),
                encoding="utf-8",
            )
            sites = exporter._get_titration_site_labels("A")
        finally:
            os.chdir(old_cwd)

    assert sites == ["al1", "ah1"]
