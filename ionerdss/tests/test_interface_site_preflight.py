"""Tests for the interface-site geometry and box-fit preflight checks.

NERDSS defines a bond's theta and phi from the vector between a molecule's centre of
mass and the reacting interface site, and orients the molecule onto its template from
the same vectors. Sites that coincide or all sit on
it therefore leave the angles undefined and NERDSS exits at the first association.
"""
from pathlib import Path
import os
import sys
import tempfile
import warnings

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ionerdss.model import pdb
from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.structure_validation import (
    COINCIDENT_SITE_TOLERANCE_NM,
    DEFAULT_INTERFACE_COM_PROXIMITY_THRESHOLD_NM,
    DegenerateSiteLayout,
    NearComInterfaceSite,
    classify_site_layout,
    get_box_fit_message,
    get_degenerate_site_layouts,
    get_designed_assembly_extent,
    get_interface_com_proximity_message,
    get_near_com_interface_sites,
    prepare_structure_validation,
)


def _make_type(name: str, interface_names) -> MoleculeType:
    """A molecule type declaring the given interfaces (2+ makes the exporter write phi)."""
    return MoleculeType(name=name, interfaces_neighbors_map={iface: None for iface in interface_names})


def _make_instance(name: str, molecule_type: MoleculeType, com) -> MoleculeInstance:
    return MoleculeInstance(
        name=name,
        molecule_type=molecule_type,
        com=np.asarray(com, dtype=float),
        norm=np.array([0.0, 0.0, 1.0]),
        ref1=np.array([1.0, 0.0, 0.0]),
        ref2=np.array([0.0, 1.0, 0.0]),
    )


def _bind(instance_a: MoleculeInstance, instance_b: MoleculeInstance, site_a, site_b, index: int = 1) -> None:
    """Bind two instances through sites at the given absolute coordinates."""
    interface_ab = InterfaceInstance(
        absolute_coord=np.asarray(site_a, dtype=float),
        this_mol=instance_a,
        this_mol_name=instance_a.name,
        partner_mol_name=instance_b.name,
        interface_index=index,
    )
    interface_ba = InterfaceInstance(
        absolute_coord=np.asarray(site_b, dtype=float),
        this_mol=instance_b,
        this_mol_name=instance_b.name,
        partner_mol_name=instance_a.name,
        interface_index=index,
    )
    interface_ab.partner_interface = interface_ba
    interface_ba.partner_interface = interface_ab
    instance_a.interfaces_neighbors_map[interface_ab] = instance_b
    instance_b.interfaces_neighbors_map[interface_ba] = instance_a


def _system(*types_and_instances) -> System:
    system = System(workspace_path=".")
    for molecule_type, instances in types_and_instances:
        system.molecule_types.add(molecule_type)
        for instance in instances:
            system.molecule_instances.add(instance)
    return system


def _collagen_like_system() -> System:
    """Three chains of one type whose two sites each coincide on the chain's COM."""
    type_a = _make_type("A", ["AA1f", "AA1b"])
    a = _make_instance("A_A", type_a, [0.0, 0.0, 0.0])
    b = _make_instance("B_A", type_a, [0.55, 0.0, 0.0])
    c = _make_instance("C_A", type_a, [0.28, 0.2, 0.0])
    for first, second in ((a, b), (b, c), (c, a)):
        _bind(first, second, first.com + [0.01, 0.0, 0.0], second.com + [0.01, 0.0, 0.0])
    return _system((type_a, [a, b, c]))


def _healthy_dimer_system(site_offset=(0.8, 0.0, 0.0)) -> System:
    type_a = _make_type("A", ["AB1"])
    type_b = _make_type("B", ["BA1"])
    a = _make_instance("A_A", type_a, [0.0, 0.0, 0.0])
    b = _make_instance("B_B", type_b, [2.0, 0.0, 0.0])
    _bind(a, b, np.asarray(site_offset), [1.2, 0.0, 0.0])
    return _system((type_a, [a]), (type_b, [b]))


# ---------------------------------------------------------------------------
# classify_site_layout
# ---------------------------------------------------------------------------

def test_classify_site_layout_ignores_single_sites():
    assert classify_site_layout([[0.0, 0.0, 0.0]]) is None
    assert classify_site_layout([]) is None


def test_classify_site_layout_detects_coincident_sites_at_any_distance():
    far = [[0.39, 0.0, 0.0], [0.39, 0.0, 0.0]]
    assert classify_site_layout(far) == "coincident"
    near = [[0.015, 0.0, 0.01], [0.015, 0.0, 0.01]]
    assert classify_site_layout(near) == "coincident"


def test_classify_site_layout_accepts_sites_on_one_line_through_the_com():
    # Sites on one line through the COM are fine: NERDSS orients such a template from
    # its first site alone, and the amyloid segments of the PDB benchmark whose two
    # sites face opposite neighbours assemble.
    assert classify_site_layout([[5.0, 0.0, 0.0], [-5.0, 0.0, 0.0]]) is None
    assert classify_site_layout([[1.0, 1.0, 0.0], [2.0, 2.0, 0.0], [-3.0, -3.0, 0.0]]) is None


def test_classify_site_layout_detects_point_like_molecules():
    layout = classify_site_layout([[0.05, 0.0, 0.0], [0.0, 0.06, 0.0]], proximity_threshold_nm=0.15)
    assert layout == "point_like"
    assert classify_site_layout([[0.05, 0.0, 0.0], [0.0, 0.06, 0.0]], proximity_threshold_nm=0.05) is None


def test_classify_site_layout_accepts_well_separated_sites():
    assert classify_site_layout([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]) is None
    assert classify_site_layout([[0.1, 0.0, 0.0], [0.0, 1.0, 0.0]]) is None


# ---------------------------------------------------------------------------
# get_near_com_interface_sites / get_degenerate_site_layouts
# ---------------------------------------------------------------------------

def test_near_com_sites_lists_every_reacting_site_below_the_threshold_sorted_by_distance():
    system = _collagen_like_system()

    near = get_near_com_interface_sites(system, threshold_nm=0.15)

    assert len(near) == 6
    assert [site.distance_nm for site in near] == sorted(site.distance_nm for site in near)
    assert {site.molecule_type for site in near} == {"A"}
    assert all(isinstance(site, NearComInterfaceSite) for site in near)
    assert all(site.phi_defined for site in near)
    assert all(abs(site.distance_nm - 0.01) < 1e-9 for site in near)
    first = near[0]
    assert first.partner_instance in {"A_A", "B_A", "C_A"} and first.partner_type == "A"


def test_near_com_sites_skips_unbound_interfaces_and_far_sites():
    system = _healthy_dimer_system()
    instance_a = next(inst for inst in system.molecule_instances if inst.name == "A_A")
    dangling = InterfaceInstance(
        absolute_coord=np.array([0.001, 0.0, 0.0]),
        this_mol=instance_a,
        this_mol_name="A_A",
        partner_mol_name="nobody",
        interface_index=7,
    )
    instance_a.interfaces_neighbors_map[dangling] = None

    assert get_near_com_interface_sites(system, threshold_nm=0.15) == []
    assert [site.distance_nm for site in get_near_com_interface_sites(system, threshold_nm=float("inf"))] == [0.8, 0.8]


def test_near_com_sites_reports_single_interface_types_as_phi_undefined():
    system = _healthy_dimer_system(site_offset=(0.02, 0.0, 0.0))

    near = get_near_com_interface_sites(system, threshold_nm=0.15)

    assert [(site.molecule_instance, site.phi_defined) for site in near] == [("A_A", False)]


def test_degenerate_layouts_flag_coincident_sites_on_multi_interface_types():
    system = _collagen_like_system()

    layouts = get_degenerate_site_layouts(system)

    assert [layout.molecule_instance for layout in layouts] == ["A_A", "B_A", "C_A"]
    assert {layout.kind for layout in layouts} == {"coincident"}
    assert all(isinstance(layout, DegenerateSiteLayout) and layout.site_count == 2 for layout in layouts)
    assert layouts[0].interface_instances == ("A_A_B_A_1", "A_A_C_A_1")
    assert layouts[0].max_site_distance_nm == pytest.approx(0.01)


def test_degenerate_layouts_exempt_single_interface_types_and_healthy_layouts():
    # A peptide with one declared interface sitting on its COM: the exporter writes
    # phi = nan for it, so NERDSS never needs to orient it.
    assert get_degenerate_site_layouts(_healthy_dimer_system(site_offset=(0.0, 0.0, 0.0))) == []
    assert get_degenerate_site_layouts(_healthy_dimer_system()) == []


def test_degenerate_layouts_use_the_threshold_for_point_like_molecules():
    type_a = _make_type("A", ["AB1", "AC1"])
    type_b = _make_type("B", ["BA1"])
    type_c = _make_type("C", ["CA1"])
    a = _make_instance("A_A", type_a, [0.0, 0.0, 0.0])
    b = _make_instance("B_B", type_b, [2.0, 0.0, 0.0])
    c = _make_instance("C_C", type_c, [0.0, 2.0, 0.0])
    _bind(a, b, [0.08, 0.0, 0.0], [1.5, 0.0, 0.0])
    _bind(a, c, [0.0, 0.1, 0.0], [0.0, 1.5, 0.0])
    system = _system((type_a, [a]), (type_b, [b]), (type_c, [c]))

    assert [layout.kind for layout in get_degenerate_site_layouts(system, threshold_nm=0.15)] == ["point_like"]
    assert get_degenerate_site_layouts(system, threshold_nm=0.09) == []


# ---------------------------------------------------------------------------
# messages
# ---------------------------------------------------------------------------

def test_interface_com_proximity_message_is_none_for_a_healthy_design():
    assert get_interface_com_proximity_message(_healthy_dimer_system(), prefix="Preflight warning") is None


def test_interface_com_proximity_message_names_the_degenerate_type_and_lists_sites():
    message = get_interface_com_proximity_message(
        _collagen_like_system(), prefix="Validation preflight warning", threshold_nm=0.15
    )

    assert message.startswith(
        "Validation preflight warning: NERDSS cannot resolve the binding angles for molecule type A "
        "(its 2 interface sites coincide, within 0.010 nm of the COM; on A_A, B_A, C_A)."
    )
    assert "Cannot resolve phi angle" in message
    assert "6 reacting interface sites within 0.15 nm" in message
    assert "A_A A_A_B_A_1 at 0.010 nm (partner B_A)" in message
    assert "interface_site_placement='auto'" in message


def test_interface_com_proximity_message_reports_a_lone_near_site_without_the_fatal_part():
    message = get_interface_com_proximity_message(
        _healthy_dimer_system(site_offset=(0.02, 0.0, 0.0)), prefix="Preflight warning", threshold_nm=0.15
    )

    assert message.startswith("Preflight warning: 1 reacting interface site within 0.15 nm")
    assert "cannot resolve the binding angles" not in message
    assert "skips phi altogether" in message


# ---------------------------------------------------------------------------
# box fit
# ---------------------------------------------------------------------------

def test_designed_assembly_extent_measures_bounding_spheres_from_the_sites():
    system = _healthy_dimer_system()  # COMs at x=0 and x=2, sites 0.8 from each COM

    extent = get_designed_assembly_extent(system)

    assert extent["largest_molecule_diameter_nm"] == pytest.approx(1.6)
    assert extent["assembly_diameter_nm"] == pytest.approx(2.0 * (1.0 + 0.8))


def test_box_fit_message_warns_only_when_the_assembly_outgrows_the_shortest_edge():
    system = _healthy_dimer_system()

    assert get_box_fit_message(system, (50.0, 50.0, 50.0), prefix="Preflight warning") is None
    assert get_box_fit_message(system, (3.7, 3.7, 3.7), prefix="Preflight warning") is None

    message = get_box_fit_message(system, (10.0, 3.0, 10.0), prefix="Preflight warning")
    assert message.startswith("Preflight warning: the designed assembly spans 3.6 nm")
    assert "10 x 3 x 10 nm" in message
    assert "Molecule seems outside simulation volume" in message
    assert "at least 4 nm" in message

    message = get_box_fit_message(system, (1.0, 1.0, 1.0), prefix="Preflight warning")
    assert message.startswith("Preflight warning: the largest molecule alone spans 1.6 nm")


# ---------------------------------------------------------------------------
# prepare_structure_validation surfaces both warnings
# ---------------------------------------------------------------------------

def _setup_in_tmpdir(system, **kwargs):
    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            os.chdir(tmpdir)
            return pdb.validation.setup_simulation(system, initial_molecule_count=1, **kwargs)
        finally:
            os.chdir(old_cwd)


def test_setup_simulation_warns_and_records_interface_com_proximity():
    with pytest.warns(RuntimeWarning, match="cannot resolve the binding angles for molecule type A"):
        artifacts = _setup_in_tmpdir(_collagen_like_system())

    assert artifacts.interface_com_proximity_warning_message.startswith(
        "Validation preflight warning: NERDSS cannot resolve the binding angles"
    )
    assert artifacts.box_fit_warning_message is None
    assert artifacts.preflight_warning_message is None


def test_setup_simulation_threshold_argument_controls_the_near_site_report():
    system = _healthy_dimer_system(site_offset=(0.1, 0.0, 0.0))

    with pytest.warns(RuntimeWarning, match="1 reacting interface site within 0.15 nm"):
        flagged = _setup_in_tmpdir(system)
    assert flagged.interface_com_proximity_warning_message is not None

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        quiet = _setup_in_tmpdir(system, interface_com_proximity_threshold_nm=0.05)
    assert quiet.interface_com_proximity_warning_message is None
    assert not any("interface site" in str(warning.message) for warning in caught)


def test_setup_simulation_takes_the_threshold_from_hyperparameters_when_not_given():
    system = _healthy_dimer_system(site_offset=(0.1, 0.0, 0.0))
    hyperparams = PDBModelHyperparameters(interface_com_proximity_threshold=0.05)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        artifacts = _setup_in_tmpdir(system, parms_overrides={"hyperparams": hyperparams})

    assert artifacts.interface_com_proximity_warning_message is None
    assert not any("interface site" in str(warning.message) for warning in caught)


def test_setup_simulation_warns_when_the_assembly_does_not_fit_the_box():
    with pytest.warns(RuntimeWarning, match="the designed assembly spans 3.6 nm"):
        artifacts = _setup_in_tmpdir(_healthy_dimer_system(), box_nm=(3.0, 3.0, 3.0))

    assert "simulation box is 3 x 3 x 3 nm" in artifacts.box_fit_warning_message
    assert artifacts.interface_com_proximity_warning_message is None


def test_default_threshold_matches_the_hyperparameter_default():
    assert PDBModelHyperparameters().interface_com_proximity_threshold == DEFAULT_INTERFACE_COM_PROXIMITY_THRESHOLD_NM
    assert COINCIDENT_SITE_TOLERANCE_NM < DEFAULT_INTERFACE_COM_PROXIMITY_THRESHOLD_NM


def test_hyperparameters_validate_the_new_fields():
    assert PDBModelHyperparameters().validate() == []
    errors = PDBModelHyperparameters(interface_com_proximity_threshold=0.0, interface_site_placement="edge").validate()
    assert "interface_com_proximity_threshold must be positive" in errors
    assert "interface_site_placement must be 'centroid' or 'auto'" in errors

    round_trip = PDBModelHyperparameters.from_dict(
        PDBModelHyperparameters(interface_com_proximity_threshold=0.2, interface_site_placement="auto").to_dict()
    )
    assert round_trip.interface_com_proximity_threshold == 0.2
    assert round_trip.interface_site_placement == "auto"
