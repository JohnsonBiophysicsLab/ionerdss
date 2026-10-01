"""The detailed summary must report the NERDSS box the regular deck was exported with.

``build_system`` exports with ``nerdss_water_box`` and falls back to its ``box_nm``
argument only when that hyperparameter is empty, but the summary used to print
``box_nm`` regardless, so a default build reported a 100 nm box for a 500 nm deck.
"""

import logging
from pathlib import Path

import numpy as np
import pytest

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model.pdb.main import PDBModelBuilder


def _instance(name, molecule_type, com):
    return MoleculeInstance(
        name=name,
        molecule_type=molecule_type,
        com=np.asarray(com, dtype=float),
        norm=np.array([0.0, 0.0, 1.0]),
        ref1=np.array([1.0, 0.0, 0.0]),
        ref2=np.array([0.0, 1.0, 0.0]),
    )


@pytest.fixture
def exported_boxes(monkeypatch):
    """Replace the PDB stages with fakes; return the boxes passed to the NERDSS export."""
    import ionerdss.model.pdb.main as pdb_main

    boxes = []

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
        def __init__(self, source, units, file_format, workspace_manager, **kwargs):
            self.chain_data = {
                "A": {"com": np.array([0.0, 0.0, 0.0])},
                "B": {"com": np.array([10.0, 0.0, 0.0])},
            }

        def get_pdb_id(self):
            return "TEST"

    class FakeCoarseGrainer:
        def __init__(self, parser, hyperparams):
            pass

        def get_summary(self):
            return {"num_interfaces": 1, "num_chains": 2}

    class FakeChainGrouper:
        def __init__(self, parser, coarse_grainer, hyperparams):
            pass

        def get_summary(self):
            return {
                "groups": [{"representative": "A", "size": 1}, {"representative": "B", "size": 1}],
                "num_groups": 2,
                "grouping_method": "test",
            }

    class FakeTemplateBuilder:
        def __init__(self, parser, coarse_grainer, chain_grouper, hyperparams, units, workspace_manager):
            pass

        def get_summary(self):
            return {"num_molecule_templates": 2, "num_interface_templates": 1}

    class FakeSystemBuilder:
        def __init__(self, parser, coarse_grainer, chain_grouper, template_builder, hyperparams,
                     workspace_path, pdb_id, units, workspace_manager):
            # A connected A-B dimer, so the build raises no preflight warnings
            self.system = System(workspace_path=workspace_path, pdb_id=pdb_id, units=units)
            type_a, type_b = MoleculeType(name="A"), MoleculeType(name="B")
            self.system.molecule_types.add(type_a)
            self.system.molecule_types.add(type_b)
            instance_a = _instance("A", type_a, [0.0, 0.0, 0.0])
            instance_b = _instance("B", type_b, [1.0, 0.0, 0.0])
            site_a = InterfaceInstance(absolute_coord=np.array([0.5, 0.0, 0.0]), this_mol=instance_a,
                                       this_mol_name="A", partner_mol_name="B", interface_index=0)
            site_b = InterfaceInstance(absolute_coord=np.array([0.5, 0.0, 0.0]), this_mol=instance_b,
                                       this_mol_name="B", partner_mol_name="A", interface_index=0)
            instance_a.interfaces_neighbors_map = {site_a: instance_b}
            instance_b.interfaces_neighbors_map = {site_b: instance_a}
            self.system.molecule_instances.add(instance_a)
            self.system.molecule_instances.add(instance_b)

        def get_system(self):
            return self.system

        def export_nerdss_files(self, molecule_counts, box_nm, parms_overrides):
            boxes.append(box_nm)
            return {}

    monkeypatch.setattr(pdb_main, "WorkspaceManager", FakeWorkspaceManager)
    monkeypatch.setattr(pdb_main, "PDBParser", FakeParser)
    monkeypatch.setattr(pdb_main, "CoarseGrainer", FakeCoarseGrainer)
    monkeypatch.setattr(pdb_main, "ChainGrouper", FakeChainGrouper)
    monkeypatch.setattr(pdb_main, "TemplateBuilder", FakeTemplateBuilder)
    monkeypatch.setattr(pdb_main, "SystemBuilder", FakeSystemBuilder)
    return boxes


@pytest.mark.parametrize(
    "build_kwargs, expected_box",
    [
        ({}, (500.0, 500.0, 500.0)),  # the nerdss_water_box default
        ({"nerdss_water_box": [300.0, 300.0, 300.0]}, (300.0, 300.0, 300.0)),
        ({"nerdss_water_box": [], "box_nm": (200.0, 200.0, 200.0)}, (200.0, 200.0, 200.0)),
    ],
    ids=["default-water-box", "custom-water-box", "box_nm-fallback"],
)
def test_detailed_summary_reports_the_exported_box(exported_boxes, tmp_path, build_kwargs, expected_box):
    workspace = tmp_path / "workspace"
    PDBModelBuilder("test_input.pdb").build_system(
        workspace_path=str(workspace),
        generate_visualizations=False,
        generate_nerdss_files=True,
        ode_enabled=False,
        **build_kwargs,
    )

    assert exported_boxes == [expected_box]
    summary = (workspace / "detailed_summary.txt").read_text(encoding="utf-8")
    assert f"Box size (nm): {expected_box}\n" in summary
