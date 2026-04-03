from pathlib import Path
import tempfile
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ionerdss.model.components.instances import InterfaceInstance
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model import pdb
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.pdb.structure_validation import (
    _find_observed_com_coordinates_in_restart_snapshots,
    _extract_observed_com_coordinates,
    align_structure_to_design,
    get_designed_structure,
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


def test_validation_counts_use_full_designed_stoichiometry():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    type_b = MoleculeType(name="B")
    system.molecule_types.add(type_a)
    system.molecule_types.add(type_b)

    system.molecule_instances.add(_make_instance("a_sparse", type_a, [0.0, 0.0, 0.0], interface_count=1))
    system.molecule_instances.add(_make_instance("a_dense", type_a, [1.0, 0.0, 0.0], interface_count=3))
    system.molecule_instances.add(_make_instance("b_only", type_b, [0.0, 2.0, 0.0], interface_count=2))

    counts = get_structure_validation_counts(system)
    coords = get_designed_structure(system)

    assert counts == {"A": 2, "B": 1}
    assert coords == {
        "a_dense": {
            "instance": "a_dense",
            "type": "A",
            "global_com_coord": (1.0, 0.0, 0.0),
            "interfaces": [],
        },
        "a_sparse": {
            "instance": "a_sparse",
            "type": "A",
            "global_com_coord": (0.0, 0.0, 0.0),
            "interfaces": [],
        },
        "b_only": {
            "instance": "b_only",
            "type": "B",
            "global_com_coord": (0.0, 2.0, 0.0),
            "interfaces": [],
        },
    }


def test_designed_structure_coordinates_include_global_interface_metadata():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    type_b = MoleculeType(name="B")
    system.molecule_types.add(type_a)
    system.molecule_types.add(type_b)

    instance_a = _make_instance("a0", type_a, [0.0, 0.0, 0.0], interface_count=0)
    instance_b = _make_instance("b0", type_b, [2.0, 0.0, 0.0], interface_count=0)

    interface_ab = InterfaceInstance(
        absolute_coord=np.array([0.5, 0.0, 0.0]),
        this_mol=instance_a,
        this_mol_name="a0",
        partner_mol_name="b0",
        interface_index=0,
    )
    interface_ba = InterfaceInstance(
        absolute_coord=np.array([1.5, 0.0, 0.0]),
        this_mol=instance_b,
        this_mol_name="b0",
        partner_mol_name="a0",
        interface_index=0,
    )

    instance_a.interfaces_neighbors_map = {interface_ab: instance_b}
    instance_b.interfaces_neighbors_map = {interface_ba: instance_a}

    system.molecule_instances.add(instance_a)
    system.molecule_instances.add(instance_b)

    coords = get_designed_structure(system)

    assert coords["a0"] == {
        "instance": "a0",
        "type": "A",
        "global_com_coord": (0.0, 0.0, 0.0),
        "interfaces": [
            {
                "interface_instance": "a0_b0_0",
                "interface_type": None,
                "binding_partner_instance": "b0",
                "binding_partner_type": "B",
                "global_interface_coord": (0.5, 0.0, 0.0),
            }
        ],
    }
    assert coords["b0"] == {
        "instance": "b0",
        "type": "B",
        "global_com_coord": (2.0, 0.0, 0.0),
        "interfaces": [
            {
                "interface_instance": "b0_a0_0",
                "interface_type": None,
                "binding_partner_instance": "a0",
                "binding_partner_type": "A",
                "global_interface_coord": (1.5, 0.0, 0.0),
            }
        ],
    }


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


def test_setup_simulation_uses_supplied_designed_coordinates():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    system.molecule_types.add(type_a)
    system.molecule_instances.add(_make_instance("a_only", type_a, [0.0, 0.0, 0.0], interface_count=0))

    designed_coordinates = {"A": (9.0, 8.0, 7.0)}

    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            import os
            os.chdir(tmpdir)
            artifacts = pdb.validation.setup_simulation(
                system,
                initial_molecule_count=1,
                designed_coordinates=designed_coordinates,
            )
            import json
            payload = json.loads(Path(artifacts.target_file).read_text(encoding="utf-8"))
        finally:
            os.chdir(old_cwd)

    assert artifacts.designed_coordinates == designed_coordinates
    assert payload["designed_coordinates"] == {"A": [9.0, 8.0, 7.0]}


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


def test_observed_structure_extraction_uses_one_connected_component():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        psf_path = tmp_path / "system.psf"
        xyz_path = tmp_path / "final_coords.xyz"
        restart_path = tmp_path / "restart.dat"

        psf_path.write_text(
            "\n".join(
                [
                    "PSF",
                    "",
                    "    6 !NATOM",
                    "       1    A    0  COM    O    0          0    1.0         0",
                    "       2    H    4  COM    O    0          0    1.0         0",
                    "       3    L    2  COM    O    0          0    1.0         0",
                    "       4    A    3  COM    O    0          0    1.0         0",
                    "       5    H    1  COM    O    0          0    1.0         0",
                    "       6    L    5  COM    O    0          0    1.0         0",
                ]
            ),
            encoding="utf-8",
        )

        xyz_path.write_text(
            "\n".join(
                [
                    "6",
                    "mol output final",
                    "A 0.0 0.0 0.0",
                    "H 10.0 0.0 0.0",
                    "L 0.0 1.0 0.0",
                    "A 10.0 10.0 0.0",
                    "H 1.0 0.0 0.0",
                    "L 10.0 11.0 0.0",
                ]
            ),
            encoding="utf-8",
        )

        def _restart_block(mol_id: int, partners: list[int], coord) -> list[str]:
            bound_iface_line = " ".join([str(len(partners))] + [str(i) for i in range(len(partners))])
            return [
                f"{mol_id} 0 0 0 0",
                "1.0 0 0 0 0 0",
                f"{coord[0]} {coord[1]} {coord[2]}",
                "4 0 1 2 3",
                bound_iface_line,
                " ".join([str(len(partners))] + [str(pid) for pid in partners]),
                "4",
                "0 0 0 0 \0 0",
                "0.0 0.0 0.0",
                "1 1 0 0 \0 0",
                "0.0 0.0 0.0",
                "2 2 0 0 \0 0",
                "0.0 0.0 0.0",
                "3 3 0 0 \0 0",
                "0.0 0.0 0.0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
            ]

        restart_lines = [
            "#All Molecules and coordinates",
            "6 6",
        ]
        restart_lines.extend(_restart_block(0, [1, 2], (0.0, 0.0, 0.0)))
        restart_lines.extend(_restart_block(1, [0, 2], (1.0, 0.0, 0.0)))
        restart_lines.extend(_restart_block(2, [0, 1], (0.0, 1.0, 0.0)))
        restart_lines.extend(_restart_block(3, [4, 5], (10.0, 10.0, 0.0)))
        restart_lines.extend(_restart_block(4, [3, 5], (10.0, 0.0, 0.0)))
        restart_lines.extend(_restart_block(5, [3, 4], (10.0, 11.0, 0.0)))
        restart_path.write_text("\n".join(restart_lines), encoding="utf-8")

        observed = _extract_observed_com_coordinates(
            system_psf_file=psf_path,
            final_coords_file=xyz_path,
            restart_file=restart_path,
            target_counts={"A": 1, "H": 1, "L": 1},
        )

    assert observed == {
        "A": (0.0, 0.0, 0.0),
        "H": (1.0, 0.0, 0.0),
        "L": (0.0, 1.0, 0.0),
    }


def test_restart_snapshot_search_falls_back_to_restart_directory():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        data_dir = tmp_path / "DATA"
        restart_dir = tmp_path / "RESTART"
        data_dir.mkdir()
        restart_dir.mkdir()

        psf_path = data_dir / "system.psf"
        xyz_path = data_dir / "final_coords.xyz"
        primary_restart_path = data_dir / "restart.dat"
        older_restart_path = restart_dir / "restart_000100.dat"

        psf_path.write_text(
            "\n".join(
                [
                    "PSF",
                    "",
                    "    3 !NATOM",
                    "       1    A    0  COM    O    0          0    1.0         0",
                    "       2    H    1  COM    O    0          0    1.0         0",
                    "       3    L    2  COM    O    0          0    1.0         0",
                ]
            ),
            encoding="utf-8",
        )

        xyz_path.write_text(
            "\n".join(
                [
                    "3",
                    "mol output final",
                    "A 100.0 0.0 0.0",
                    "H 200.0 0.0 0.0",
                    "L 300.0 0.0 0.0",
                ]
            ),
            encoding="utf-8",
        )

        def _restart_block(mol_id: int, partners: list[int], coord) -> list[str]:
            bound_iface_line = " ".join([str(len(partners))] + [str(i) for i in range(len(partners))])
            return [
                f"{mol_id} 0 0 0 0",
                "1.0 0 0 0 0 0",
                f"{coord[0]} {coord[1]} {coord[2]}",
                "4 0 1 2 3",
                bound_iface_line,
                " ".join([str(len(partners))] + [str(pid) for pid in partners]),
                "4",
                "0 0 0 0 \0 0",
                "0.0 0.0 0.0",
                "1 1 0 0 \0 0",
                "0.0 0.0 0.0",
                "2 2 0 0 \0 0",
                "0.0 0.0 0.0",
                "3 3 0 0 \0 0",
                "0.0 0.0 0.0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
            ]

        primary_lines = [
            "#All Molecules and coordinates",
            "3 3",
        ]
        primary_lines.extend(_restart_block(0, [], (10.0, 0.0, 0.0)))
        primary_lines.extend(_restart_block(1, [], (20.0, 0.0, 0.0)))
        primary_lines.extend(_restart_block(2, [], (30.0, 0.0, 0.0)))
        primary_restart_path.write_text("\n".join(primary_lines), encoding="utf-8")

        older_lines = [
            "#All Molecules and coordinates",
            "3 3",
        ]
        older_lines.extend(_restart_block(0, [1, 2], (0.0, 0.0, 0.0)))
        older_lines.extend(_restart_block(1, [0, 2], (1.0, 0.0, 0.0)))
        older_lines.extend(_restart_block(2, [0, 1], (0.0, 1.0, 0.0)))
        older_restart_path.write_text("\n".join(older_lines), encoding="utf-8")

        observed, used_restart = _find_observed_com_coordinates_in_restart_snapshots(
            system_psf_file=psf_path,
            final_coords_file=xyz_path,
            restart_file=primary_restart_path,
            target_counts={"A": 1, "H": 1, "L": 1},
        )

    assert used_restart == older_restart_path
    assert observed == {
        "A": (0.0, 0.0, 0.0),
        "H": (1.0, 0.0, 0.0),
        "L": (0.0, 1.0, 0.0),
    }


def test_observed_structure_extraction_falls_back_to_restart_native_molecule_types():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        psf_path = tmp_path / "system.psf"
        restart_path = tmp_path / "restart.dat"

        # Deliberately mismap the COM entries so PSF-based composition matching fails.
        psf_path.write_text(
            "\n".join(
                [
                    "PSF",
                    "",
                    "    3 !NATOM",
                    "       1    A    0  COM    O    0          0    1.0         0",
                    "       2    A    1  COM    O    0          0    1.0         0",
                    "       3    A    2  COM    O    0          0    1.0         0",
                ]
            ),
            encoding="utf-8",
        )

        def _restart_block(mol_id: int, coord, iface_count: int, partners: list[int]) -> list[str]:
            lines = [
                f"{mol_id} 0 0 0 0",
                "1.0 0 0 0 0 0",
                f"{coord[0]} {coord[1]} {coord[2]}",
                " ".join([str(iface_count)] + [str(i) for i in range(iface_count)]),
                " ".join([str(len(partners))] + [str(i) for i in range(len(partners))]),
                " ".join([str(len(partners))] + [str(pid) for pid in partners]),
                str(iface_count),
            ]

            bound_iface_indexes = set(range(max(iface_count - len(partners), 0), iface_count))
            for iface_idx in range(iface_count):
                is_bound = 1 if iface_idx in bound_iface_indexes else 0
                lines.append(f"{iface_idx} {iface_idx} 0 0 0 {is_bound}")
                lines.append("0.0 0.0 0.0")
                if is_bound:
                    lines.append("0 0 0")

            lines.extend(["0", "0", "0", "0", "0", "0"])
            return lines

        restart_lines = [
            "#MolTemplates",
            "0 A 4",
            "1 B 2",
            "#All Molecules and coordinates",
            "3 3",
        ]
        restart_lines.extend(_restart_block(0, (0.0, 0.0, 0.0), 4, [2]))
        restart_lines.extend(_restart_block(1, (1.0, 0.0, 0.0), 4, [2]))
        restart_lines.extend(_restart_block(2, (0.0, 1.0, 0.0), 2, [0, 1]))
        restart_path.write_text("\n".join(restart_lines), encoding="utf-8")

        observed = _extract_observed_com_coordinates(
            system_psf_file=psf_path,
            final_coords_file=None,
            restart_file=restart_path,
            target_counts={"A": 2, "B": 1},
        )

    assert observed == {
        "A_0": (0.0, 0.0, 0.0),
        "A_1": (1.0, 0.0, 0.0),
        "B": (0.0, 1.0, 0.0),
    }
