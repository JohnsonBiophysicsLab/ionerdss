from pathlib import Path
import tempfile
import sys
import logging
import json

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ionerdss.model.components.instances import InterfaceInstance
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model import pdb
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.structure_validation import (
    _extract_observed_com_coordinates_from_complex_json,
    _find_observed_com_coordinates_in_complex_json_snapshots,
    _find_observed_com_coordinates_in_restart_snapshots,
    _extract_observed_com_coordinates,
    _format_disconnected_design_warning,
    _get_designed_connected_components,
    align_structure_to_design,
    get_disconnected_design_message,
    get_designed_structure,
    get_structure_validation_counts,
    collect_structure_validation_results,
    run_structure_validation_simulation,
    StructureValidationArtifacts,
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


def test_align_structure_to_design_keeps_exact_key_matching():
    designed = {
        "A_0": [0.0, 0.0, 0.0],
        "A_1": [1.0, 0.0, 0.0],
        "B_0": [0.0, 1.0, 0.0],
    }

    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([2.0, 3.0, -1.0])

    designed_xyz = np.asarray([designed[key] for key in sorted(designed)], dtype=float)
    observed_xyz = (rotation @ designed_xyz.T).T + translation
    observed = {
        key: observed_xyz[idx].tolist()
        for idx, key in enumerate(sorted(designed))
    }

    result = align_structure_to_design(designed, observed)

    assert result.labels == ("A_0", "A_1", "B_0")
    assert result.rmsd < 1e-10


def test_align_structure_to_design_matches_homomer_labels_by_type():
    designed = {
        "chainA_A": [0.0, 0.0, 0.0],
        "chainB_A": [2.0, 0.0, 0.0],
        "chainC_B": [0.0, 3.0, 0.0],
    }
    observed = {
        "A_0": [2.0, 0.0, 0.0],
        "A_1": [0.0, 0.0, 0.0],
        "B_0": [0.0, 3.0, 0.0],
    }

    result = align_structure_to_design(designed, observed)

    assert result.labels == ("A", "A", "B")
    assert result.rmsd < 1e-10
    assert np.allclose(result.aligned_observed_coordinates, result.designed_coordinates)


def test_align_structure_to_design_scales_to_large_homomer():
    # An 18-copy homomer like 5L93 has 18! relabelings, far past what can be enumerated.
    angles = np.arange(6) * (2.0 * np.pi / 6.0)
    asymmetric_unit = np.array([[6.0, 0.0, 0.0], [7.0, 1.0, 0.5], [6.5, -1.0, 1.5]])
    points = np.concatenate(
        [
            (
                np.array(
                    [
                        [np.cos(angle), -np.sin(angle), 0.0],
                        [np.sin(angle), np.cos(angle), 0.0],
                        [0.0, 0.0, 1.0],
                    ]
                )
                @ asymmetric_unit.T
            ).T
            for angle in angles
        ]
    )

    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([4.0, -2.0, 7.0])
    moved = (rotation @ points.T).T + translation

    designed = {f"chain{idx:02d}_A": points[idx] for idx in range(len(points))}
    # NERDSS numbers copies in assembly order, which need not follow the designed order.
    shuffled = np.roll(np.arange(len(points)), 7)
    observed = {f"A_{idx}": moved[shuffled[idx]] for idx in range(len(points))}

    result = align_structure_to_design(designed, observed)

    assert result.labels == ("A",) * len(points)
    assert result.rmsd < 1e-8
    assert np.allclose(result.aligned_observed_coordinates, result.designed_coordinates, atol=1e-8)


def test_validation_module_exposes_prepare_and_compare():
    assert hasattr(pdb, "validation")
    assert callable(pdb.validation.prepare)
    assert callable(pdb.validation.compare)
    assert callable(pdb.validation.setup_simulation)
    assert callable(pdb.validation.align_structure)


def test_designed_connected_components_detect_disconnected_subunits():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    type_b = MoleculeType(name="B")
    type_c = MoleculeType(name="C")
    type_d = MoleculeType(name="D")
    system.molecule_types.add(type_a)
    system.molecule_types.add(type_b)
    system.molecule_types.add(type_c)
    system.molecule_types.add(type_d)

    instance_a = _make_instance("A", type_a, [0.0, 0.0, 0.0], interface_count=0)
    instance_b = _make_instance("B", type_b, [1.0, 0.0, 0.0], interface_count=0)
    instance_c = _make_instance("C", type_c, [5.0, 0.0, 0.0], interface_count=0)
    instance_d = _make_instance("D", type_d, [6.0, 0.0, 0.0], interface_count=0)

    interface_ab = InterfaceInstance(
        absolute_coord=np.array([0.5, 0.0, 0.0]),
        this_mol=instance_a,
        this_mol_name="A",
        partner_mol_name="B",
        interface_index=0,
    )
    interface_ba = InterfaceInstance(
        absolute_coord=np.array([0.5, 0.0, 0.0]),
        this_mol=instance_b,
        this_mol_name="B",
        partner_mol_name="A",
        interface_index=0,
    )
    interface_cd = InterfaceInstance(
        absolute_coord=np.array([5.5, 0.0, 0.0]),
        this_mol=instance_c,
        this_mol_name="C",
        partner_mol_name="D",
        interface_index=0,
    )
    interface_dc = InterfaceInstance(
        absolute_coord=np.array([5.5, 0.0, 0.0]),
        this_mol=instance_d,
        this_mol_name="D",
        partner_mol_name="C",
        interface_index=0,
    )

    instance_a.interfaces_neighbors_map = {interface_ab: instance_b}
    instance_b.interfaces_neighbors_map = {interface_ba: instance_a}
    instance_c.interfaces_neighbors_map = {interface_cd: instance_d}
    instance_d.interfaces_neighbors_map = {interface_dc: instance_c}

    system.molecule_instances.add(instance_a)
    system.molecule_instances.add(instance_b)
    system.molecule_instances.add(instance_c)
    system.molecule_instances.add(instance_d)

    components = _get_designed_connected_components(system)

    assert components == [("A", "B"), ("C", "D")]
    assert _format_disconnected_design_warning(components) == (
        "Validation preflight warning: the designed assembly graph is disconnected, so it cannot form a "
        "single N-mer. Subunits A, B are disconnected from subunits C, D."
    )
    assert get_disconnected_design_message(system, prefix="Preflight error") == (
        "Preflight error: the designed assembly graph is disconnected, so it cannot form a single N-mer. "
        "Subunits A, B are disconnected from subunits C, D."
    )


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
    assert artifacts.preflight_warning_message is None


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


def test_setup_simulation_warns_for_disconnected_designed_assembly():
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    type_b = MoleculeType(name="B")
    system.molecule_types.add(type_a)
    system.molecule_types.add(type_b)
    system.molecule_instances.add(_make_instance("A", type_a, [0.0, 0.0, 0.0], interface_count=0))
    system.molecule_instances.add(_make_instance("B", type_b, [1.0, 0.0, 0.0], interface_count=0))

    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            import os
            os.chdir(tmpdir)
            with pytest.warns(RuntimeWarning, match="designed assembly graph is disconnected"):
                artifacts = pdb.validation.setup_simulation(system, initial_molecule_count=1)
        finally:
            os.chdir(old_cwd)

    assert artifacts.preflight_warning_message == (
        "Validation preflight warning: the designed assembly graph is disconnected, so it cannot form a "
        "single N-mer. Subunits A are disconnected from subunits B."
    )


def test_build_system_raises_early_for_disconnected_designed_assembly(monkeypatch, tmp_path):
    import ionerdss.model.pdb.main as pdb_main

    class FakeWorkspaceManager:
        def __init__(self, workspace_path, pdb_id):
            self.workspace_path = Path(workspace_path)
            self.workspace_path.mkdir(parents=True, exist_ok=True)
            self.logger = logging.getLogger(f"test-workspace-{pdb_id}")

        def get_system_output_path(self):
            return self.workspace_path / "system.json"

        def get_report_path(self, report_type):
            return self.workspace_path / f"{report_type}.txt"

    class FakeParser:
        def __init__(self, source, units, file_format, workspace_manager):
            self.source = source
            self.units = units
            self.file_format = file_format
            self.workspace_manager = workspace_manager
            self.chain_data = {
                "A": {"com": np.array([0.0, 0.0, 0.0])},
                "B": {"com": np.array([10.0, 0.0, 0.0])},
            }

        def get_pdb_id(self):
            return "TEST"

    class FakeCoarseGrainer:
        def __init__(self, parser, hyperparams):
            self.parser = parser
            self.hyperparams = hyperparams

        def get_summary(self):
            return {"num_interfaces": 0, "num_chains": 2}

    class FakeChainGrouper:
        def __init__(self, parser, coarse_grainer, hyperparams):
            self.parser = parser

        def get_summary(self):
            return {
                "groups": [{"representative": "A", "size": 1}, {"representative": "B", "size": 1}],
                "num_groups": 2,
                "grouping_method": "test",
            }

    class FakeTemplateBuilder:
        def __init__(self, parser, coarse_grainer, chain_grouper, hyperparams, units, workspace_manager):
            self.parser = parser

        def get_summary(self):
            return {"num_molecule_templates": 2, "num_interface_templates": 0}

    class FakeSystemBuilder:
        def __init__(self, parser, coarse_grainer, chain_grouper, template_builder, hyperparams, workspace_path, pdb_id, units, workspace_manager):
            self.system = System(workspace_path=workspace_path, pdb_id=pdb_id, units=units)
            type_a = MoleculeType(name="A")
            type_b = MoleculeType(name="B")
            self.system.molecule_types.add(type_a)
            self.system.molecule_types.add(type_b)
            self.system.molecule_instances.add(_make_instance("A", type_a, [0.0, 0.0, 0.0], interface_count=0))
            self.system.molecule_instances.add(_make_instance("B", type_b, [10.0, 0.0, 0.0], interface_count=0))

        def get_system(self):
            return self.system

    monkeypatch.setattr(pdb_main, "WorkspaceManager", FakeWorkspaceManager)
    monkeypatch.setattr(pdb_main, "PDBParser", FakeParser)
    monkeypatch.setattr(pdb_main, "CoarseGrainer", FakeCoarseGrainer)
    monkeypatch.setattr(pdb_main, "ChainGrouper", FakeChainGrouper)
    monkeypatch.setattr(pdb_main, "TemplateBuilder", FakeTemplateBuilder)
    monkeypatch.setattr(pdb_main, "SystemBuilder", FakeSystemBuilder)

    builder = PDBModelBuilder("test_input.pdb")

    with pytest.warns(RuntimeWarning, match="Preflight warning: the designed assembly graph is disconnected"):
        builder.build_system(
            workspace_path=str(tmp_path / "workspace"),
            generate_visualizations=False,
            generate_nerdss_files=False,
            ode_enabled=False,
        )


def test_build_system_allows_connected_designed_assembly(monkeypatch, tmp_path):
    import ionerdss.model.pdb.main as pdb_main

    class FakeWorkspaceManager:
        def __init__(self, workspace_path, pdb_id):
            self.workspace_path = Path(workspace_path)
            self.workspace_path.mkdir(parents=True, exist_ok=True)
            self.logger = logging.getLogger(f"test-workspace-{pdb_id}")

        def get_system_output_path(self):
            return self.workspace_path / "system.json"

        def get_report_path(self, report_type):
            return self.workspace_path / f"{report_type}.txt"

    class FakeParser:
        def __init__(self, source, units, file_format, workspace_manager):
            self.source = source
            self.units = units
            self.file_format = file_format
            self.workspace_manager = workspace_manager
            self.chain_data = {
                "A": {"com": np.array([0.0, 0.0, 0.0])},
                "B": {"com": np.array([1.0, 0.0, 0.0])},
            }

        def get_pdb_id(self):
            return "TEST"

    class FakeCoarseGrainer:
        def __init__(self, parser, hyperparams):
            self.parser = parser
            self.hyperparams = hyperparams

        def get_summary(self):
            return {"num_interfaces": 1, "num_chains": 2}

    class FakeChainGrouper:
        def __init__(self, parser, coarse_grainer, hyperparams):
            self.parser = parser

        def get_summary(self):
            return {
                "groups": [{"representative": "A", "size": 1}, {"representative": "B", "size": 1}],
                "num_groups": 2,
                "grouping_method": "test",
            }

    class FakeTemplateBuilder:
        def __init__(self, parser, coarse_grainer, chain_grouper, hyperparams, units, workspace_manager):
            self.parser = parser

        def get_summary(self):
            return {"num_molecule_templates": 2, "num_interface_templates": 1}

    class FakeSystemBuilder:
        def __init__(self, parser, coarse_grainer, chain_grouper, template_builder, hyperparams, workspace_path, pdb_id, units, workspace_manager):
            self.system = System(workspace_path=workspace_path, pdb_id=pdb_id, units=units)
            type_a = MoleculeType(name="A")
            type_b = MoleculeType(name="B")
            self.system.molecule_types.add(type_a)
            self.system.molecule_types.add(type_b)

            instance_a = _make_instance("A", type_a, [0.0, 0.0, 0.0], interface_count=0)
            instance_b = _make_instance("B", type_b, [1.0, 0.0, 0.0], interface_count=0)
            interface_ab = InterfaceInstance(
                absolute_coord=np.array([0.5, 0.0, 0.0]),
                this_mol=instance_a,
                this_mol_name="A",
                partner_mol_name="B",
                interface_index=0,
            )
            interface_ba = InterfaceInstance(
                absolute_coord=np.array([0.5, 0.0, 0.0]),
                this_mol=instance_b,
                this_mol_name="B",
                partner_mol_name="A",
                interface_index=0,
            )
            instance_a.interfaces_neighbors_map = {interface_ab: instance_b}
            instance_b.interfaces_neighbors_map = {interface_ba: instance_a}

            self.system.molecule_instances.add(instance_a)
            self.system.molecule_instances.add(instance_b)

        def get_system(self):
            return self.system

    monkeypatch.setattr(pdb_main, "WorkspaceManager", FakeWorkspaceManager)
    monkeypatch.setattr(pdb_main, "PDBParser", FakeParser)
    monkeypatch.setattr(pdb_main, "CoarseGrainer", FakeCoarseGrainer)
    monkeypatch.setattr(pdb_main, "ChainGrouper", FakeChainGrouper)
    monkeypatch.setattr(pdb_main, "TemplateBuilder", FakeTemplateBuilder)
    monkeypatch.setattr(pdb_main, "SystemBuilder", FakeSystemBuilder)

    builder = PDBModelBuilder("test_input.pdb")
    system = builder.build_system(
        workspace_path=str(tmp_path / "workspace"),
        generate_visualizations=False,
        generate_nerdss_files=False,
        ode_enabled=False,
    )

    assert system is builder.system
    assert len(system.molecule_instances) == 2


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

        def _restart_block(mol_id: int, iface_count: int, partners: list[int], coord) -> list[str]:
            bound_iface_line = " ".join([str(len(partners))] + [str(i) for i in range(len(partners))])
            return [
                f"{mol_id} 0 0 0 0",
                "1.0 0 0 0 0 0",
                f"{coord[0]} {coord[1]} {coord[2]}",
                " ".join([str(iface_count)] + [str(i) for i in range(iface_count)]),
                bound_iface_line,
                " ".join([str(len(partners))] + [str(pid) for pid in partners]),
                str(iface_count),
                *[
                    token
                    for iface_idx in range(iface_count)
                    for token in (
                        f"{iface_idx} {iface_idx} 0 0 \0 0",
                        "0.0 0.0 0.0",
                    )
                ],
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
            ]

        restart_lines = [
            "#MolTemplates",
            "0 A 4",
            "1 H 3",
            "2 L 2",
            "#All Molecules and coordinates",
            "6 6",
        ]
        restart_lines.extend(_restart_block(0, 4, [1, 2], (0.0, 0.0, 0.0)))
        restart_lines.extend(_restart_block(1, 3, [0, 2], (1.0, 0.0, 0.0)))
        restart_lines.extend(_restart_block(2, 2, [0, 1], (0.0, 1.0, 0.0)))
        restart_lines.extend(_restart_block(3, 4, [4, 5], (10.0, 10.0, 0.0)))
        restart_lines.extend(_restart_block(4, 3, [3, 5], (10.0, 0.0, 0.0)))
        restart_lines.extend(_restart_block(5, 2, [3, 4], (10.0, 11.0, 0.0)))
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

        def _restart_block(mol_id: int, iface_count: int, partners: list[int], coord) -> list[str]:
            bound_iface_line = " ".join([str(len(partners))] + [str(i) for i in range(len(partners))])
            return [
                f"{mol_id} 0 0 0 0",
                "1.0 0 0 0 0 0",
                f"{coord[0]} {coord[1]} {coord[2]}",
                " ".join([str(iface_count)] + [str(i) for i in range(iface_count)]),
                bound_iface_line,
                " ".join([str(len(partners))] + [str(pid) for pid in partners]),
                str(iface_count),
                *[
                    token
                    for iface_idx in range(iface_count)
                    for token in (
                        f"{iface_idx} {iface_idx} 0 0 \0 0",
                        "0.0 0.0 0.0",
                    )
                ],
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
            ]

        primary_lines = [
            "#MolTemplates",
            "0 A 4",
            "1 H 3",
            "2 L 2",
            "#All Molecules and coordinates",
            "3 3",
        ]
        primary_lines.extend(_restart_block(0, 4, [], (10.0, 0.0, 0.0)))
        primary_lines.extend(_restart_block(1, 3, [], (20.0, 0.0, 0.0)))
        primary_lines.extend(_restart_block(2, 2, [], (30.0, 0.0, 0.0)))
        primary_restart_path.write_text("\n".join(primary_lines), encoding="utf-8")

        older_lines = [
            "#MolTemplates",
            "0 A 4",
            "1 H 3",
            "2 L 2",
            "#All Molecules and coordinates",
            "3 3",
        ]
        older_lines.extend(_restart_block(0, 4, [1, 2], (0.0, 0.0, 0.0)))
        older_lines.extend(_restart_block(1, 3, [0, 2], (1.0, 0.0, 0.0)))
        older_lines.extend(_restart_block(2, 2, [0, 1], (0.0, 1.0, 0.0)))
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


def test_extract_observed_com_coordinates_from_complex_json():
    with tempfile.TemporaryDirectory() as tmpdir:
        json_path = Path(tmpdir) / "99999.json"
        json_path.write_text(
            json.dumps(
                [
                    {"names": ["A", "A"], "coords": [[9.0, 9.0, 9.0], [8.0, 8.0, 8.0]]},
                    {"names": ["A", "A", "B"], "coords": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]},
                ]
            ),
            encoding="utf-8",
        )

        observed = _extract_observed_com_coordinates_from_complex_json(json_path, {"A": 2, "B": 1})

    assert observed == {
        "A_0": (0.0, 0.0, 0.0),
        "A_1": (1.0, 0.0, 0.0),
        "B": (0.0, 1.0, 0.0),
    }


def test_find_observed_com_coordinates_in_complex_json_snapshots_prefers_latest():
    with tempfile.TemporaryDirectory() as tmpdir:
        complexes_dir = Path(tmpdir) / "COMPLEXES"
        complexes_dir.mkdir()
        (complexes_dir / "10000.json").write_text(
            json.dumps([{"names": ["A", "A", "B"], "coords": [[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]]}]),
            encoding="utf-8",
        )
        (complexes_dir / "99999.json").write_text(
            json.dumps([{"names": ["A", "A", "B"], "coords": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]}]),
            encoding="utf-8",
        )

        observed, used_json = _find_observed_com_coordinates_in_complex_json_snapshots(
            complexes_dir, {"A": 2, "B": 1}
        )

    assert used_json.name == "99999.json"
    assert observed == {
        "A_0": (0.0, 0.0, 0.0),
        "A_1": (1.0, 0.0, 0.0),
        "B": (0.0, 1.0, 0.0),
    }


def test_run_structure_validation_simulation_prefers_complex_json(monkeypatch, tmp_path):
    work_dir = tmp_path / "nerdss_files"
    data_dir = work_dir / "validation_output" / "1" / "DATA"
    complexes_dir = data_dir / "COMPLEXES"
    complexes_dir.mkdir(parents=True)

    parms_path = work_dir / "parms.inp"
    parms_path.parent.mkdir(parents=True, exist_ok=True)
    parms_path.write_text("parms", encoding="utf-8")
    (data_dir / "histogram_complexes_time.dat").write_text("Time (s): 0\n", encoding="utf-8")
    (data_dir / "final_coords.xyz").write_text("0\ncomment\n", encoding="utf-8")
    (data_dir / "system.psf").write_text("PSF\n", encoding="utf-8")
    (data_dir / "restart.dat").write_text("unused", encoding="utf-8")
    (complexes_dir / "99999.json").write_text(
        json.dumps([{"names": ["A", "A", "B"], "coords": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]}]),
        encoding="utf-8",
    )

    artifacts = StructureValidationArtifacts(
        molecule_counts={"A": 2, "B": 1},
        target_counts={"A": 2, "B": 1},
        designed_coordinates={"A_0": (0.0, 0.0, 0.0), "A_1": (1.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)},
        target_file=work_dir / "target.json",
        nerdss_files={"parms": parms_path},
    )

    class FakeSimulation:
        def __init__(self, path):
            self.path = path
            self.parmfile = None

        def run_new_simulations(self, **kwargs):
            return None

    import ionerdss.nerdss_simulation as nerdss_simulation
    import ionerdss.analysis.io.parser as parser_module

    monkeypatch.setattr(nerdss_simulation, "Simulation", FakeSimulation)
    monkeypatch.setattr(
        parser_module,
        "parse_complex_histogram",
        lambda path: (np.asarray([]), [], None),
    )

    result = run_structure_validation_simulation(artifacts, nerdss_dir=tmp_path / "bin")

    assert result.full_assembly_found is True
    assert result.warning_message is None
    assert result.restart_file.name == "99999.json"
    assert result.observed_coordinates == {
        "A_0": (0.0, 0.0, 0.0),
        "A_1": (1.0, 0.0, 0.0),
        "B": (0.0, 1.0, 0.0),
    }


def test_run_structure_validation_simulation_falls_back_to_restart_with_warning(monkeypatch, tmp_path):
    work_dir = tmp_path / "nerdss_files"
    data_dir = work_dir / "validation_output" / "1" / "DATA"
    data_dir.mkdir(parents=True)

    parms_path = work_dir / "parms.inp"
    parms_path.parent.mkdir(parents=True, exist_ok=True)
    parms_path.write_text("parms", encoding="utf-8")
    (data_dir / "histogram_complexes_time.dat").write_text("Time (s): 0\n", encoding="utf-8")
    (data_dir / "final_coords.xyz").write_text("0\ncomment\n", encoding="utf-8")
    (data_dir / "system.psf").write_text("PSF\n", encoding="utf-8")
    (data_dir / "restart.dat").write_text("unused", encoding="utf-8")

    artifacts = StructureValidationArtifacts(
        molecule_counts={"A": 2, "B": 1},
        target_counts={"A": 2, "B": 1},
        designed_coordinates={"A_0": (0.0, 0.0, 0.0), "A_1": (1.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)},
        target_file=work_dir / "target.json",
        nerdss_files={"parms": parms_path},
    )

    class FakeSimulation:
        def __init__(self, path):
            self.path = path
            self.parmfile = None

        def run_new_simulations(self, **kwargs):
            return None

    import ionerdss.nerdss_simulation as nerdss_simulation
    import ionerdss.analysis.io.parser as parser_module
    import ionerdss.model.pdb.structure_validation as validation_module

    monkeypatch.setattr(nerdss_simulation, "Simulation", FakeSimulation)
    monkeypatch.setattr(
        parser_module,
        "parse_complex_histogram",
        lambda path: (np.asarray([]), [], None),
    )
    monkeypatch.setattr(
        validation_module,
        "_find_observed_com_coordinates_in_restart_snapshots",
        lambda **kwargs: (
            {"A_0": (0.0, 0.0, 0.0), "A_1": (1.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)},
            Path(kwargs["restart_file"]),
        ),
    )

    result = run_structure_validation_simulation(artifacts, nerdss_dir=tmp_path / "bin")

    assert result.full_assembly_found is True
    assert "no COMPLEXES JSON snapshots were found" in (result.warning_message or "")
    assert result.restart_file.name == "restart.dat"


def test_run_structure_validation_simulation_passes_env_to_nerdss(monkeypatch, tmp_path):
    work_dir = tmp_path / "nerdss_files"
    data_dir = work_dir / "validation_output" / "1" / "DATA"
    complexes_dir = data_dir / "COMPLEXES"
    complexes_dir.mkdir(parents=True)

    parms_path = work_dir / "parms.inp"
    parms_path.write_text("parms", encoding="utf-8")
    (data_dir / "histogram_complexes_time.dat").write_text("Time (s): 0\n", encoding="utf-8")
    (data_dir / "final_coords.xyz").write_text("0\ncomment\n", encoding="utf-8")
    (data_dir / "system.psf").write_text("PSF\n", encoding="utf-8")
    (data_dir / "restart.dat").write_text("unused", encoding="utf-8")
    (complexes_dir / "99999.json").write_text(
        json.dumps([{"names": ["A", "B"], "coords": [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]}]),
        encoding="utf-8",
    )

    artifacts = StructureValidationArtifacts(
        molecule_counts={"A": 1, "B": 1},
        target_counts={"A": 1, "B": 1},
        designed_coordinates={"A": (0.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)},
        target_file=work_dir / "target.json",
        nerdss_files={"parms": parms_path},
    )

    captured = {}

    class FakeSimulation:
        def __init__(self, path):
            self.path = path
            self.parmfile = None

        def run_new_simulations(self, **kwargs):
            captured.update(kwargs)
            return None

    import ionerdss.nerdss_simulation as nerdss_simulation
    import ionerdss.analysis.io.parser as parser_module

    monkeypatch.setattr(nerdss_simulation, "Simulation", FakeSimulation)
    monkeypatch.setattr(
        parser_module,
        "parse_complex_histogram",
        lambda path: (np.asarray([]), [], None),
    )

    run_structure_validation_simulation(
        artifacts,
        nerdss_dir=tmp_path / "bin",
        env={"LD_LIBRARY_PATH": "/opt/gsl/lib"},
    )

    assert captured["env"] == {"LD_LIBRARY_PATH": "/opt/gsl/lib"}


def test_run_new_simulations_merges_env_onto_os_environ(monkeypatch, tmp_path):
    from ionerdss.nerdss_simulation.simulation import Simulation

    work_dir = tmp_path / "work"
    work_dir.mkdir()
    (work_dir / "parms.inp").write_text("parms", encoding="utf-8")

    nerdss_bin_dir = tmp_path / "NERDSS" / "bin"
    nerdss_bin_dir.mkdir(parents=True)
    (nerdss_bin_dir / "nerdss").write_text("#!/bin/sh\n", encoding="utf-8")

    captured = {}

    class FakeProcess:
        def wait(self):
            return 0

        def poll(self):
            return 0

    def fake_popen(cmd, **kwargs):
        captured["env"] = kwargs.get("env")
        return FakeProcess()

    monkeypatch.setattr("ionerdss.nerdss_simulation.simulation.subprocess.Popen", fake_popen)
    monkeypatch.setenv("IONERDSS_TEST_INHERITED", "inherited")

    simulation = Simulation(str(work_dir))
    simulation.run_new_simulations(
        sim_dir=str(tmp_path / "out"),
        nerdss_dir=str(tmp_path / "NERDSS"),
        progress=False,
        verbose=False,
        env={"LD_LIBRARY_PATH": "/opt/gsl/lib"},
    )

    assert captured["env"]["LD_LIBRARY_PATH"] == "/opt/gsl/lib"
    assert captured["env"]["IONERDSS_TEST_INHERITED"] == "inherited"

    simulation.run_new_simulations(
        sim_dir=str(tmp_path / "out"),
        nerdss_dir=str(tmp_path / "NERDSS"),
        progress=False,
        verbose=False,
    )

    assert captured["env"] is None


def _write_external_validation_run(run_dir):
    """Lay out the DATA tree a NERDSS run launched outside ioNERDSS would leave behind."""
    data_dir = run_dir / "DATA"
    complexes_dir = data_dir / "COMPLEXES"
    complexes_dir.mkdir(parents=True)
    (data_dir / "histogram_complexes_time.dat").write_text("Time (s): 0\n", encoding="utf-8")
    (data_dir / "final_coords.xyz").write_text("0\ncomment\n", encoding="utf-8")
    (data_dir / "system.psf").write_text("PSF\n", encoding="utf-8")
    (data_dir / "restart.dat").write_text("unused", encoding="utf-8")
    (complexes_dir / "99999.json").write_text(
        json.dumps([{"names": ["A", "A", "B"], "coords": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]}]),
        encoding="utf-8",
    )
    return data_dir


def _external_run_artifacts(work_dir):
    parms_path = work_dir / "parms.inp"
    parms_path.parent.mkdir(parents=True, exist_ok=True)
    parms_path.write_text("parms", encoding="utf-8")
    return StructureValidationArtifacts(
        molecule_counts={"A": 2, "B": 1},
        target_counts={"A": 2, "B": 1},
        designed_coordinates={"A_0": (0.0, 0.0, 0.0), "A_1": (1.0, 0.0, 0.0), "B": (0.0, 1.0, 0.0)},
        target_file=work_dir / "target.json",
        nerdss_files={"parms": parms_path},
    )


def test_collect_structure_validation_results_reads_external_run(monkeypatch, tmp_path):
    run_dir = tmp_path / "nerdss_files"
    _write_external_validation_run(run_dir)
    artifacts = _external_run_artifacts(run_dir)

    import ionerdss.analysis.io.parser as parser_module

    monkeypatch.setattr(
        parser_module,
        "parse_complex_histogram",
        lambda path: (np.asarray([]), [], None),
    )

    result = collect_structure_validation_results(artifacts, simulation_dir=run_dir)

    assert result.simulation_dir == run_dir
    assert result.full_assembly_found is True
    assert result.warning_message is None
    assert result.restart_file.name == "99999.json"
    assert result.observed_coordinates == {
        "A_0": (0.0, 0.0, 0.0),
        "A_1": (1.0, 0.0, 0.0),
        "B": (0.0, 1.0, 0.0),
    }


def test_collect_structure_validation_results_accepts_data_dir(monkeypatch, tmp_path):
    run_dir = tmp_path / "nerdss_files"
    data_dir = _write_external_validation_run(run_dir)
    artifacts = _external_run_artifacts(run_dir)

    import ionerdss.analysis.io.parser as parser_module

    monkeypatch.setattr(
        parser_module,
        "parse_complex_histogram",
        lambda path: (np.asarray([]), [], None),
    )

    result = collect_structure_validation_results(artifacts, simulation_dir=data_dir)

    assert result.simulation_dir == run_dir
    assert result.full_assembly_found is True
    assert result.histogram_file == data_dir / "histogram_complexes_time.dat"
