"""Tests for exporting the structure-validation deck.

The deck asks NERDSS for DATA/COMPLEXES snapshots, without which the result readers see
only the final restart snapshot and miss an assembly that formed and then grew. It is
written to its own directory, so it cannot replace a regular export in nerdss_files/.
When build_system exports it, it carries the builder's hyperparameters and targets every
designed molecule instance, not one per molecule type.
"""
from pathlib import Path
import json
import logging
import re

import numpy as np
import pytest

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import InterfaceType, MoleculeType
from ionerdss.model.pdb import validation
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.pdb.structure_validation import align_structure_to_design
from ionerdss.model.pdb.system_builder import SystemBuilder


DESIGNED_COMS = {"A_A": (0.0, 0.0, 0.0), "B_A": (4.0, 0.0, 0.0)}


def _homodimer_system() -> System:
    """Two copies of molecule type A bound through the homotypic interface AA1."""
    system = System(workspace_path=".")
    type_a = MoleculeType(
        name="A",
        radius_nm=2.0,
        D_t_nm2_us=1.0,
        D_r_rad2_us=0.1,
        ref1_local=np.array([1.0, 0.0, 0.0]),
        ref2_local=np.array([0.0, 0.0, 1.0]),
    )
    system.molecule_types.add(type_a)

    a_binds_a = InterfaceType(
        this_mol_type_name="A",
        partner_mol_type_name="A",
        interface_index=1,
        absolute_coord=np.array([2.0, 0.0, 0.0]),
        local_coord=np.array([2.0, 0.0, 0.0]),
        energy=-5.0,
        this_mol_type=type_a,
        partner_mol_type=type_a,
    )
    a_binds_a.partner_interface_type = a_binds_a
    system.interface_types.add(a_binds_a)

    first, second = (
        MoleculeInstance(
            name=name,
            molecule_type=type_a,
            com=np.asarray(com, dtype=float),
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([ref1_x, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0]),
        )
        for (name, com), ref1_x in zip(DESIGNED_COMS.items(), (1.0, -1.0))
    )
    system.molecule_instances.add(first)
    system.molecule_instances.add(second)

    first_site = InterfaceInstance(
        absolute_coord=np.array([1.9, 0.0, 0.0]),
        interface_type=a_binds_a,
        this_mol=first,
        this_mol_name="A",
        partner_mol_name="A",
        interface_index=1,
    )
    second_site = InterfaceInstance(
        absolute_coord=np.array([2.1, 0.0, 0.0]),
        interface_type=a_binds_a,
        this_mol=second,
        this_mol_name="A",
        partner_mol_name="A",
        interface_index=1,
    )
    first_site.partner_interface = second_site
    second_site.partner_interface = first_site
    first.interfaces_neighbors_map[first_site] = second
    second.interfaces_neighbors_map[second_site] = first
    return system


def _parameters(parms_path) -> dict:
    """The key = value pairs of a parms file's parameters block."""
    block = Path(parms_path).read_text(encoding="utf-8").split("end parameters")[0]
    return dict(re.findall(r"^\s*(\w+) = (\S+)", block, re.MULTILINE))


def _molecule_counts(parms_path) -> dict:
    text = Path(parms_path).read_text(encoding="utf-8")
    block = re.search(r"start molecules\n(.*?)end molecules", text, re.DOTALL).group(1)
    return {name: int(count) for name, count in re.findall(r"(\w+) : (\d+)", block)}


def _water_box(parms_path) -> str:
    return re.search(r"WaterBox = (\[.*?\])", Path(parms_path).read_text(encoding="utf-8")).group(1)


def _directory_contents(directory: Path) -> dict:
    return {path.name: path.read_text(encoding="utf-8") for path in sorted(directory.iterdir())}


class _FakeWorkspaceManager:
    def __init__(self, workspace_path, pdb_id=None):
        self.workspace_path = Path(workspace_path)
        self.workspace_path.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(f"test-validation-deck-{pdb_id}")

    def get_system_output_path(self):
        return self.workspace_path / "system.json"

    def get_report_path(self, report_type):
        return self.workspace_path / f"{report_type}.txt"


# ---------------------------------------------------------------------------
# bondedComplexWrite
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "parms_overrides, expected",
    [
        ({"nItr": 2_000_000}, 20_000),
        ({}, 1_000),
        ({"hyperparams": PDBModelHyperparameters(nerdss_n_itr=50_000)}, 500),
        ({"nItr": 2_000_000, "hyperparams": PDBModelHyperparameters(nerdss_n_itr=50_000)}, 20_000),
        ({"nItr": 50}, 1),
        ({"nItr": 2_000_000, "bondedComplexWrite": 7}, 7),
    ],
    ids=[
        "explicit-nItr",
        "exporter-default-nItr",
        "hyperparameter-nItr",
        "explicit-nItr-over-hyperparameters",
        "at-least-every-step",
        "explicit-bondedComplexWrite",
    ],
)
def test_validation_deck_asks_nerdss_for_complex_snapshots(tmp_path, monkeypatch, parms_overrides, expected):
    monkeypatch.chdir(tmp_path)
    caller_overrides = dict(parms_overrides)

    artifacts = validation.setup_simulation(_homodimer_system(), parms_overrides=parms_overrides)

    for deck in ("parms", "parms_titrate"):
        parameters = _parameters(artifacts.nerdss_files[deck])
        assert parameters["bondedComplexWrite"] == str(expected)
    assert parms_overrides == caller_overrides


# ---------------------------------------------------------------------------
# Deck directory
# ---------------------------------------------------------------------------

def test_validation_deck_leaves_a_regular_export_in_the_workspace_alone(tmp_path):
    system = _homodimer_system()
    workspace = _FakeWorkspaceManager(tmp_path / "workspace")
    NERDSSExporter(system, workspace).export_all(
        molecule_counts={"A": 10}, box_nm=(500.0, 500.0, 500.0)
    )
    regular_dir = workspace.workspace_path / "nerdss_files"
    regular_deck = _directory_contents(regular_dir)

    artifacts = validation.setup_simulation(system, workspace_manager=workspace)

    assert _directory_contents(regular_dir) == regular_deck
    deck_dir = workspace.workspace_path / "structure_validation"
    assert sorted(path.name for path in deck_dir.iterdir()) == [
        "A.mol",
        "parms.inp",
        "parms_titrate.inp",
        "structure_validation_target.json",
    ]
    assert {Path(path).parent for path in artifacts.nerdss_files.values()} == {deck_dir}
    assert artifacts.target_file == deck_dir / "structure_validation_target.json"
    assert _molecule_counts(artifacts.nerdss_files["parms"]) == {"A": 2}


def test_validation_deck_dir_is_relative_to_the_workspace_or_the_current_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    system = _homodimer_system()

    without_workspace = validation.setup_simulation(system)
    assert without_workspace.nerdss_files["parms"] == Path("structure_validation") / "parms.inp"
    assert not (tmp_path / "nerdss_files").exists()

    workspace = _FakeWorkspaceManager(tmp_path / "workspace")
    custom = validation.prepare(system, workspace_manager=workspace, deck_dir="two_copies")
    assert custom.nerdss_files["parms"] == workspace.workspace_path / "two_copies" / "parms.inp"
    assert custom.target_file == workspace.workspace_path / "two_copies" / "structure_validation_target.json"
    assert not (workspace.workspace_path / "nerdss_files").exists()


# ---------------------------------------------------------------------------
# Transition counting
# ---------------------------------------------------------------------------

def test_validation_deck_does_not_count_transitions(tmp_path, monkeypatch):
    # Titration keeps adding copies, so a complex can outgrow any transitionMatrixSize.
    # This one is too small even for the two starting copies, which the exporter
    # refuses whenever transitions are counted.
    monkeypatch.chdir(tmp_path)
    hyperparams = PDBModelHyperparameters(count_transition=True, transition_matrix_size=1)

    with pytest.warns(RuntimeWarning, match="count_transition off"):
        artifacts = validation.setup_simulation(
            _homodimer_system(), parms_overrides={"hyperparams": hyperparams}
        )

    mol_file = Path(artifacts.nerdss_files["A_mol"]).read_text(encoding="utf-8")
    assert "countTransition = false" in mol_file
    assert hyperparams.count_transition is True


# ---------------------------------------------------------------------------
# build_system(structure_validation=True)
# ---------------------------------------------------------------------------

class _FakeParser:
    # Chain COMs are kept in Angstrom, as the real parser keeps them.
    chain_data = {
        "A": {"com": np.array([0.0, 0.0, 0.0])},
        "B": {"com": np.array([40.0, 0.0, 0.0])},
    }

    def __init__(self, source, units, file_format, workspace_manager, **kwargs):
        pass

    def get_pdb_id(self):
        return "TEST"

    def get_chain_data(self, chain_id):
        return self.chain_data[chain_id]

    def convert_coords_to_nm(self, coords):
        return np.asarray(coords, dtype=float) / 10.0


class _FakeCoarseGrainer:
    def __init__(self, parser, hyperparams):
        pass

    def get_summary(self):
        return {"num_interfaces": 1, "num_chains": 2}


class _FakeChainGrouper:
    def __init__(self, parser, coarse_grainer, hyperparams):
        pass

    def get_summary(self):
        return {
            "groups": [{"representative": "A", "size": 2}],
            "num_groups": 1,
            "grouping_method": "test",
        }


class _FakeTemplateBuilder:
    def __init__(self, parser, coarse_grainer, chain_grouper, hyperparams, units, workspace_manager):
        pass

    def get_summary(self):
        return {"num_molecule_templates": 1, "num_interface_templates": 1}


class _PrebuiltSystemBuilder(SystemBuilder):
    """The real SystemBuilder exports, around the homodimer built by hand."""

    def __init__(self, parser, coarse_grainer, chain_grouper, template_builder, hyperparams,
                 workspace_path, pdb_id, units, workspace_manager):
        self.parser = parser
        self.coarse_grainer = coarse_grainer
        self.hyperparams = hyperparams
        self.workspace_manager = workspace_manager
        self.system = _homodimer_system()


@pytest.fixture
def fake_pipeline(monkeypatch):
    import ionerdss.model.pdb.main as pdb_main

    monkeypatch.setattr(pdb_main, "WorkspaceManager", _FakeWorkspaceManager)
    monkeypatch.setattr(pdb_main, "PDBParser", _FakeParser)
    monkeypatch.setattr(pdb_main, "CoarseGrainer", _FakeCoarseGrainer)
    monkeypatch.setattr(pdb_main, "ChainGrouper", _FakeChainGrouper)
    monkeypatch.setattr(pdb_main, "TemplateBuilder", _FakeTemplateBuilder)
    monkeypatch.setattr(pdb_main, "SystemBuilder", _PrebuiltSystemBuilder)


def _build(workspace: Path, **kwargs) -> PDBModelBuilder:
    builder = PDBModelBuilder("test_input.pdb")
    builder.build_system(
        workspace_path=str(workspace),
        generate_visualizations=False,
        ode_enabled=False,
        nerdss_total_molecule_count=12,
        nerdss_water_box=[300.0, 300.0, 300.0],
        nerdss_overlap_sep_limit=2.5,
        nerdss_n_itr=40_000,
        structure_validation=True,
        **kwargs,
    )
    return builder


def test_build_system_keeps_the_regular_deck_beside_the_validation_deck(fake_pipeline, tmp_path):
    workspace = tmp_path / "workspace"
    builder = _build(workspace)

    regular = workspace / "nerdss_files" / "parms.inp"
    assert _molecule_counts(regular) == {"A": 12}
    assert _water_box(regular) == "[300.0, 300.0, 300.0]"
    assert "0 -> " not in regular.read_text(encoding="utf-8")

    deck_dir = workspace / "structure_validation"
    artifacts = builder.structure_validation_artifacts
    assert Path(artifacts.nerdss_files["parms_titrate"]) == deck_dir / "parms_titrate.inp"
    assert _molecule_counts(deck_dir / "parms_titrate.inp") == {"A": 2}
    assert _water_box(deck_dir / "parms_titrate.inp") == "[100.0, 100.0, 100.0]"


def test_build_system_exports_the_validation_deck_with_the_builder_hyperparameters(fake_pipeline, tmp_path):
    builder = _build(tmp_path / "workspace")

    parameters = _parameters(builder.structure_validation_artifacts.nerdss_files["parms_titrate"])
    assert parameters["overlapSepLimit"] == "2.5"
    assert parameters["nItr"] == "40000"
    assert parameters["bondedComplexWrite"] == "400"


def test_build_system_validation_options_override_the_builder_hyperparameters(fake_pipeline, tmp_path):
    builder = _build(
        tmp_path / "workspace",
        structure_validation_options={
            "box_nm": (60.0, 60.0, 60.0),
            "parms_overrides": {"nItr": 8_000, "overlapSepLimit": 1.5},
        },
    )

    parms_titrate = builder.structure_validation_artifacts.nerdss_files["parms_titrate"]
    parameters = _parameters(parms_titrate)
    assert parameters["nItr"] == "8000"
    assert parameters["overlapSepLimit"] == "1.5"
    assert parameters["bondedComplexWrite"] == "80"
    assert _water_box(parms_titrate) == "[60.0, 60.0, 60.0]"


def test_build_system_validation_target_holds_every_designed_copy(fake_pipeline, tmp_path):
    builder = _build(tmp_path / "workspace")

    artifacts = builder.structure_validation_artifacts
    target = json.loads(Path(artifacts.target_file).read_text(encoding="utf-8"))
    assert target["molecule_counts"] == {"A": 2}
    assert target["designed_coordinates"] == {name: list(com) for name, com in DESIGNED_COMS.items()}
    assert artifacts.designed_coordinates == DESIGNED_COMS

    # NERDSS names the two copies of the assembled dimer A_0 and A_1, in its own order.
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = {name: rotation @ np.asarray(com) + [10.0, -5.0, 3.0] for name, com in DESIGNED_COMS.items()}
    observed = {"A_0": moved["B_A"], "A_1": moved["A_A"]}
    alignment = align_structure_to_design(artifacts.designed_coordinates, observed)
    assert alignment.rmsd < 1e-9
