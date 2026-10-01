"""Tests for hyperparameters passed as keyword arguments.

``build_system``, ``build_system_from_pdb`` and ``set_hyperparameters`` take
``PDBModelHyperparameters`` fields as keywords. A name that is not a field, whether
misspelled or an older short form such as ``distance_cutoff``, must raise instead of
being dropped, which would leave the build running on the default value.
"""
from pathlib import Path
import logging
import sys

import pytest
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import ionerdss.model.pdb.main as pdb_main
from ionerdss.api import build_system_from_pdb
from ionerdss.model.components.units import Units
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.main import PDBModelBuilder


class _StopAfterHyperparameters(Exception):
    """Raised by the stand-in parser, the first step after the keywords are merged."""


class _StoppingParser:
    def __init__(self, **kwargs):
        raise _StopAfterHyperparameters


class _FakeWorkspaceManager:
    def __init__(self, workspace_path, pdb_id):
        self.logger = logging.getLogger(f"test-hyperparameter-keywords-{pdb_id}")


def test_build_system_rejects_an_unknown_keyword_before_creating_the_workspace(tmp_path):
    workspace = tmp_path / "workspace"

    with pytest.raises(TypeError, match=r"build_system\(\) got an unexpected keyword argument 'logger_level'"):
        PDBModelBuilder("1abc").build_system(
            workspace_path=str(workspace),
            generate_nerdss_files=False,
            logger_level=logging.WARNING,
        )

    assert not workspace.exists()


def test_build_system_names_the_current_field_for_an_old_short_name(tmp_path):
    with pytest.raises(
        TypeError,
        match=r"'distance_cutoff' \(did you mean 'interface_detect_distance_cutoff'\?\)",
    ):
        PDBModelBuilder("1abc").build_system(
            workspace_path=str(tmp_path / "workspace"), distance_cutoff=1.0
        )


def test_build_system_applies_keywords_to_a_copy_that_keeps_units_and_aligner(monkeypatch, tmp_path):
    monkeypatch.setattr(pdb_main, "WorkspaceManager", _FakeWorkspaceManager)
    monkeypatch.setattr(pdb_main, "PDBParser", _StoppingParser)
    units = Units(coords="angstrom")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    hyperparams = PDBModelHyperparameters(units=units, chain_grouping_custom_aligner=aligner)
    builder = PDBModelBuilder("1abc")

    with pytest.raises(_StopAfterHyperparameters):
        builder.build_system(
            workspace_path=str(tmp_path / "workspace"),
            hyperparams=hyperparams,
            interface_detect_distance_cutoff=1.0,
        )

    assert builder.hyperparams.interface_detect_distance_cutoff == 1.0
    assert builder.hyperparams.units is units
    assert builder.hyperparams.chain_grouping_custom_aligner is aligner
    assert hyperparams.interface_detect_distance_cutoff == 0.9


def test_build_system_from_pdb_rejects_an_unknown_keyword_the_same_way(tmp_path):
    workspace = tmp_path / "workspace"

    with pytest.raises(
        TypeError,
        match=r"build_system_from_pdb\(\) got an unexpected keyword argument "
              r"'matching_mode' \(did you mean 'chain_grouping_matching_mode'\?\)",
    ):
        build_system_from_pdb("1abc", workspace_path=str(workspace), matching_mode="sequence")

    assert not workspace.exists()


def test_set_hyperparameters_rejects_an_unknown_keyword_when_creating_and_updating():
    builder = PDBModelBuilder("1abc")

    with pytest.raises(TypeError, match=r"set_hyperparameters\(\) got an unexpected keyword argument 'logger_level'"):
        builder.set_hyperparameters(logger_level=logging.WARNING)
    assert builder.hyperparams is None

    builder.set_hyperparameters(interface_detect_distance_cutoff=1.0)
    aligner = builder.hyperparams.chain_grouping_custom_aligner

    with pytest.raises(TypeError, match=r"'rmsd_threshold' \(did you mean 'chain_grouping_rmsd_threshold'\?\)"):
        builder.set_hyperparameters(rmsd_threshold=1.0)

    builder.set_hyperparameters(chain_grouping_rmsd_threshold=1.0)
    assert builder.hyperparams.interface_detect_distance_cutoff == 1.0
    assert builder.hyperparams.chain_grouping_rmsd_threshold == 1.0
    assert builder.hyperparams.chain_grouping_custom_aligner is aligner
