"""Tests for exporting the structure-validation deck.

The deck asks NERDSS for DATA/COMPLEXES snapshots, without which the result readers see
only the final restart snapshot and miss an assembly that formed and then grew.
"""
from pathlib import Path
import re

import numpy as np
import pytest

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import InterfaceType, MoleculeType
from ionerdss.model.pdb import validation
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters


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
