"""Tests for interface_site_placement='auto' and the exporter's collinear-phi rule."""
from pathlib import Path
import math
import os
import sys
import tempfile
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.pdb.system_builder import (
    SystemBuilder,
    principal_axis,
    project_site_to_chain_surface,
)


def _rod(length=80.0, radius=3.0, n_rings=50, per_ring=12, axis=(0.0, 0.0, 1.0)):
    """Atoms on rings around ``axis``, centred on the origin (an exactly symmetric rod)."""
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    seed = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(axis, seed); u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    t = np.linspace(-length / 2, length / 2, n_rings)
    angle = np.arange(per_ring) * 2 * np.pi / per_ring
    ring = radius * (np.outer(np.cos(angle), u) + np.outer(np.sin(angle), v))
    return np.concatenate([np.outer(np.full(per_ring, ti), axis) + ring for ti in t])


def _blob(radius=10.0, n=500):
    """Points filling a sphere of the given radius."""
    rng = np.random.default_rng(0)
    points = rng.normal(size=(n, 3))
    points /= np.linalg.norm(points, axis=1)[:, None]
    return points * radius * rng.random((n, 1)) ** (1 / 3)


def test_principal_axis_reports_direction_and_aspect_ratio():
    axis, aspect = principal_axis(_rod())
    assert abs(float(axis @ [0.0, 0.0, 1.0])) == pytest.approx(1.0)
    assert aspect > 3.0

    _, blob_aspect = principal_axis(_blob())
    assert blob_aspect < 1.5

    assert principal_axis(np.zeros((1, 3))) == (None, 1.0)


def test_projection_puts_the_site_on_the_surface_facing_the_partner():
    coords = _blob(radius=10.0)
    com = np.zeros(3)

    site = project_site_to_chain_surface(com, coords, partner_com=np.array([30.0, 0.0, 0.0]))

    expected = np.max(coords[:, 0])
    assert site == pytest.approx([expected, 0.0, 0.0])
    assert 8.0 < expected <= 10.0


def test_projection_keeps_an_alongside_partner_on_the_flank_of_a_rod():
    coords = _rod(length=80.0, radius=3.0)  # long axis z, thickness 3
    com = np.zeros(3)
    partner = np.array([0.8, 0.0, 5.0])  # staggered along the rod, slightly to +x

    site = project_site_to_chain_surface(com, coords, partner_com=partner)

    # On the flank (x = rod radius), not at the tip (z = 40).
    assert site[2] == pytest.approx(0.0, abs=1e-6)
    assert site[0] == pytest.approx(3.0, abs=1e-6)
    assert site[1] == pytest.approx(0.0, abs=1e-6)


def test_projection_follows_the_full_direction_when_the_partner_is_beside_a_compact_chain():
    coords = _rod(length=20.0, radius=8.0)  # aspect below the elongated threshold
    partner = np.array([0.0, 0.0, 40.0])

    site = project_site_to_chain_surface(np.zeros(3), coords, partner_com=partner)

    assert site[0] == pytest.approx(0.0, abs=1e-9) and site[1] == pytest.approx(0.0, abs=1e-9)
    assert site[2] == pytest.approx(10.0, abs=0.5)


def test_projection_falls_back_when_the_centres_of_mass_coincide():
    coords = _blob(radius=10.0)
    com = np.zeros(3)

    site = project_site_to_chain_surface(com, coords, partner_com=com, fallback_direction=np.array([0.0, 2.0, 0.0]))
    assert site[0] == pytest.approx(0.0, abs=1e-9) and site[2] == pytest.approx(0.0, abs=1e-9)
    assert site[1] == pytest.approx(np.max(coords[:, 1]))

    site = project_site_to_chain_surface(com, coords, partner_com=com)
    assert site[0] == pytest.approx(np.max(coords[:, 0]))


def test_find_degenerate_site_chains_mirrors_the_layout_classifier():
    chain_data = {
        "A": {"com": np.array([0.0, 0.0, 0.0])},
        "B": {"com": np.array([5.5, 0.0, 0.0])},
        "C": {"com": np.array([2.8, 20.0, 0.0])},
    }
    builder = SystemBuilder.__new__(SystemBuilder)
    builder.parser = SimpleNamespace(
        get_chain_data=lambda chain_id: chain_data[chain_id],
        convert_coords_to_nm=lambda coords: np.asarray(coords, dtype=float) / 10.0,
    )
    # A's two sites coincide on its COM; B's and C's are well separated and off-line.
    interfaces = [
        SimpleNamespace(chain_i="A", chain_j="B", coord_i=np.array([0.1, 0.0, 0.0]), coord_j=np.array([0.0, 0.0, 0.0])),
        SimpleNamespace(chain_i="A", chain_j="C", coord_i=np.array([0.1, 0.0, 0.0]), coord_j=np.array([2.8, 10.0, 0.0])),
        SimpleNamespace(chain_i="B", chain_j="C", coord_i=np.array([5.5, 8.0, 0.0]), coord_j=np.array([9.0, 12.0, 0.0])),
    ]

    degenerate = builder._find_degenerate_site_chains(interfaces, proximity_threshold_nm=0.15)

    assert degenerate == {"A": "coincident"}


def _two_type_system(site1, site2, com2=(3.0, 0.0, 0.0)):
    """A and B, each declaring two interfaces so the exporter computes phi, bound through one."""
    type_a = MoleculeType(name="A", interfaces_neighbors_map={"AB1": None, "AC1": None})
    type_b = MoleculeType(name="B", interfaces_neighbors_map={"BA1": None, "BC1": None})
    a = MoleculeInstance(name="A_A", molecule_type=type_a, com=np.zeros(3), norm=np.array([0.0, 0.0, 1.0]),
                         ref1=np.array([1.0, 0.0, 0.0]), ref2=np.array([0.0, 1.0, 0.0]))
    b = MoleculeInstance(name="B_B", molecule_type=type_b, com=np.asarray(com2, dtype=float), norm=np.array([0.0, 0.0, 1.0]),
                         ref1=np.array([1.0, 0.0, 0.0]), ref2=np.array([0.0, 1.0, 0.0]))
    iface_ab = InterfaceInstance(absolute_coord=np.asarray(site1, dtype=float), this_mol=a, this_mol_name="A_A", partner_mol_name="B_B", interface_index=1)
    iface_ba = InterfaceInstance(absolute_coord=np.asarray(site2, dtype=float), this_mol=b, this_mol_name="B_B", partner_mol_name="A_A", interface_index=1)
    a.interfaces_neighbors_map[iface_ab] = b
    b.interfaces_neighbors_map[iface_ba] = a
    system = System(workspace_path=".")
    system.molecule_types.add(type_a); system.molecule_types.add(type_b)
    system.molecule_instances.add(a); system.molecule_instances.add(b)
    return system, a, b, iface_ab, iface_ba


def _angles_for(system, a, b, iface_ab, iface_ba):
    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            os.chdir(tmpdir)
            exporter = NERDSSExporter(system)
            return exporter._generate_reaction_angles(
                iface_ab.absolute_coord, iface_ba.absolute_coord, a.com, b.com, "A", "B", "ab1", "ba1"
            )
        finally:
            os.chdir(old_cwd)


def test_exporter_writes_nan_phi_when_the_site_lies_on_the_line_joining_the_sites():
    # COM_A -- site_A -- site_B -- COM_B all on the x axis: theta is pi on both sides
    # and sigma has no component perpendicular to the COM-to-site axis to measure phi from.
    sigma, (theta1, theta2, phi1, phi2, omega) = _angles_for(*_two_type_system([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]))

    assert sigma == pytest.approx(1.0)
    assert theta1 == pytest.approx(math.pi) and theta2 == pytest.approx(math.pi)
    assert math.isnan(phi1) and math.isnan(phi2)
    assert not math.isnan(omega)


def test_exporter_keeps_phi_for_sites_off_the_line():
    sigma, (theta1, theta2, phi1, phi2, omega) = _angles_for(*_two_type_system([1.0, 0.0, 0.0], [2.0, 0.6, 0.0]))

    assert sigma == pytest.approx(math.hypot(1.0, 0.6))
    assert 0.0 < theta1 < math.pi and 0.0 < theta2 < math.pi
    assert not math.isnan(phi1) and not math.isnan(phi2)
