"""Regression tests for steric_clash_mode="auto".

Covers the two halves of the feature:

1. ``TemplateBuilder._detect_steric_clashes`` must actually populate
   ``required_free`` for a quasi-equivalent assembly (a T=3 capsid is the
   motivating case: one over-valent template that behaves as a 3-state
   molecule).
2. ``NERDSSExporter`` must round-trip the reactions that pass produces, which
   carry comma-separated ancillary (required-free) sites, e.g. ``A(aa1f,aa4b)``.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ionerdss.model.components.types import InterfaceType, MoleculeType
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.pdb.template_builder import TemplateBuilder


# --------------------------------------------------------------------------
# Fixtures: a synthetic quasi-equivalent assembly
# --------------------------------------------------------------------------

# 12 directional interfaces partitioned into three quasi-equivalent conformers
# of four, plus one shared contact (AA2) present on every chain. Each 'f' side
# and its 'b' partner deliberately land on different conformers.
CONFORMERS = {
    "A": ["AA1f", "AA3f", "AA5f", "AA7f"],
    "B": ["AA1b", "AA3b", "AA4f", "AA6f"],
    "C": ["AA4b", "AA5b", "AA6b", "AA7b"],
}
SHARED = "AA2"
DIRECTIONAL = [name for names in CONFORMERS.values() for name in names]


def _make_interface(name: str, local_coord) -> InterfaceType:
    intf = InterfaceType(
        this_mol_type_name="A",
        partner_mol_type_name="A",
        interface_index=1,
        absolute_coord=np.asarray(local_coord, dtype=float),
        local_coord=np.asarray(local_coord, dtype=float),
    )
    intf.set_name(f"A_A_{name[2:]}")
    return intf


def _make_builder(chain_to_types, interface_templates, radius_nm=0.0):
    """Build a TemplateBuilder wired up just enough to run clash detection.

    ``chain_to_types`` maps a chain id to the interface type names observed on
    it; it is replayed through the coarse-grainer records and the
    interface->type mapping that ``_detect_steric_clashes`` reconstructs.
    """
    builder = TemplateBuilder.__new__(TemplateBuilder)
    builder.workspace_manager = None
    builder.interface_templates = interface_templates
    builder.molecule_templates = {"A": MoleculeType(name="A", radius_nm=radius_nm)}

    records = []
    mapping = {}
    for tick, (chain, type_names) in enumerate(sorted(chain_to_types.items())):
        for offset, type_name in enumerate(sorted(type_names)):
            coord_i = np.array([float(tick), float(offset), 0.0])
            coord_j = np.array([float(tick), float(offset), 1.0])
            records.append(
                SimpleNamespace(
                    chain_i=chain, chain_j="_partner",
                    coord_i=coord_i, coord_j=coord_j,
                )
            )
            key_i = (
                f"{chain}__partner_{coord_i[0]:.3f}_{coord_i[1]:.3f}_{coord_i[2]:.3f}"
            )
            mapping[key_i] = type_name

    builder.coarse_grainer = SimpleNamespace(interfaces=records)
    builder.interface_to_type_mapping = mapping
    return builder


def _quasi_equivalent_builder():
    templates = {
        name: _make_interface(name, (float(i), 0.0, 0.0))
        for i, name in enumerate(DIRECTIONAL + [SHARED])
    }
    chain_to_types = {
        f"chain_{conf}": set(names) | {SHARED}
        for conf, names in CONFORMERS.items()
    }
    return _make_builder(chain_to_types, templates)


# --------------------------------------------------------------------------
# 1. The detector runs and finds the quasi-equivalence
# --------------------------------------------------------------------------

def test_auto_steric_detection_populates_required_free():
    """The pass must not silently no-op: it once keyed off ``this_mol_type``,
    which is still None at template-build time, so nothing was ever detected."""
    builder = _quasi_equivalent_builder()

    builder._detect_steric_clashes()

    assert any(t.required_free for t in builder.interface_templates.values())


def test_auto_steric_detection_recovers_quasi_equivalence_exactly():
    builder = _quasi_equivalent_builder()

    builder._detect_steric_clashes()

    templates = builder.interface_templates

    # The shared contact co-occurs with everything, so it excludes nothing.
    assert templates[SHARED].required_free == []
    for name in DIRECTIONAL:
        assert SHARED not in templates[name].required_free

    # Each directional interface excludes exactly the 8 belonging to the other
    # two conformers, and none of the 3 sharing its own conformer.
    for conformer, names in CONFORMERS.items():
        others = [n for c, ns in CONFORMERS.items() if c != conformer for n in ns]
        for name in names:
            required = set(templates[name].required_free)
            assert required == set(others), name

    pairs = {
        frozenset((name, other))
        for name, tmpl in templates.items()
        for other in tmpl.required_free
    }
    assert len(pairs) == 48


def test_auto_steric_detection_is_symmetric():
    builder = _quasi_equivalent_builder()

    builder._detect_steric_clashes()

    templates = builder.interface_templates
    for name, tmpl in templates.items():
        for other in tmpl.required_free:
            assert name in templates[other].required_free


def test_cooccurring_interfaces_are_never_excluded_despite_proximity():
    """Co-occurrence is the primary rule: two interfaces seen together on a
    real chain stay compatible even when they are geometrically adjacent."""
    templates = {
        "AA1f": _make_interface("AA1f", (0.0, 0.0, 0.0)),
        "AA1b": _make_interface("AA1b", (0.01, 0.0, 0.0)),
    }
    builder = _make_builder(
        {"chain_A": {"AA1f", "AA1b"}}, templates, radius_nm=5.0
    )

    builder._detect_steric_clashes()

    assert templates["AA1f"].required_free == []
    assert templates["AA1b"].required_free == []


# --------------------------------------------------------------------------
# 2. The geometric heuristic (secondary) no longer discards large assemblies
# --------------------------------------------------------------------------

def test_geometric_fallback_uses_radius_threshold_beyond_2_5_nm():
    """Without chain-level evidence, the radius-based threshold decides. A
    fixed 2.5 nm early-out used to discard candidates before the threshold was
    computed, which under-detects on assemblies whose interfaces spread wider
    than that from the COM."""
    templates = {
        "AA1f": _make_interface("AA1f", (0.0, 0.0, 0.0)),
        "AA3f": _make_interface("AA3f", (3.4, 0.0, 0.0)),
    }
    # No observed chain occupancy -> co-occurrence cannot decide.
    builder = _make_builder({}, templates, radius_nm=3.45)

    builder._detect_steric_clashes()

    # threshold = r1 + r2 = 6.9 nm > 3.4 nm separation
    assert templates["AA1f"].required_free == ["AA3f"]
    assert templates["AA3f"].required_free == ["AA1f"]


def test_geometric_fallback_leaves_distant_interfaces_compatible():
    templates = {
        "AA1f": _make_interface("AA1f", (0.0, 0.0, 0.0)),
        "AA3f": _make_interface("AA3f", (10.0, 0.0, 0.0)),
    }
    builder = _make_builder({}, templates, radius_nm=1.0)

    builder._detect_steric_clashes()

    assert templates["AA1f"].required_free == []
    assert templates["AA3f"].required_free == []


# --------------------------------------------------------------------------
# 3. The exporter round-trips reactions carrying ancillary sites
# --------------------------------------------------------------------------

ANCILLARY_REACTION = (
    "A(aa1f,aa4b) + A(aa1b,aa3f) <-> A(aa1f!1,aa4b).A(aa1b!1,aa3f)"
)


def _make_exporter(tmp_path: Path) -> NERDSSExporter:
    exporter = NERDSSExporter.__new__(NERDSSExporter)
    exporter.workspace_manager = None
    exporter.output_dir = tmp_path
    exporter.precalculated_rates = {}
    exporter.interface_to_site_map = {
        "AA1f": "aa1f", "AA1b": "aa1b", "AA3f": "aa3f", "AA4b": "aa4b",
    }
    interface_types = [
        InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="A",
            interface_index=1,
            absolute_coord=np.zeros(3),
            local_coord=np.zeros(3),
            energy=-20.0,
        )
    ]
    exporter.system = SimpleNamespace(
        molecule_types=[MoleculeType(name="A", radius_nm=2.0, D_t_nm2_us=1.0)],
        interface_types=interface_types,
    )
    exporter._get_representative_instance = lambda mol_name: object()
    exporter._local_x_with_degeneracy = lambda mol, site: np.array([1.0, 0.0, 0.0])
    return exporter


def test_exporter_round_trips_reaction_with_ancillary_sites(tmp_path):
    """The narrow ``[A-Za-z0-9_]+`` site regexes rejected ``A(aa1f,aa4b)``,
    dropping into a fallback that left ``mol1`` unbound and crashed."""
    exporter = _make_exporter(tmp_path)

    parms_path = exporter._write_parms_file(
        reactions=[ANCILLARY_REACTION],
        molecule_counts={"A": 100},
        box_nm=(100.0, 100.0, 100.0),
        sigma_list=[1.0],
        angles_list=[(0.0, 0.0, 0.0, 0.0, 0.0)],
    )

    text = parms_path.read_text(encoding="utf-8")
    assert ANCILLARY_REACTION in text
    # The parse succeeded, so the real per-site normals were written rather
    # than the [0, 0, 1] fallback.
    assert "norm1 = [1.0, 0.0, 0.0]" in text
    assert "norm2 = [1.0, 0.0, 0.0]" in text


def test_exporter_resolves_primary_site_ignoring_ancillary_sites(tmp_path):
    """Only the first site in the group is the binding site; the rest are
    steric-exclusion requirements and must not reach the site lookups."""
    exporter = _make_exporter(tmp_path)
    seen = []
    exporter._local_x_with_degeneracy = lambda mol, site: (
        seen.append((mol, site)) or np.array([1.0, 0.0, 0.0])
    )
    exporter.precalculated_rates = {("A", "aa1f", "A", "aa1b"): (5.0, 7.0)}

    parms_path = exporter._write_parms_file(
        reactions=[ANCILLARY_REACTION],
        molecule_counts={"A": 100},
        box_nm=(100.0, 100.0, 100.0),
        sigma_list=[1.0],
        angles_list=[(0.0, 0.0, 0.0, 0.0, 0.0)],
    )

    assert seen == [("A", "aa1f"), ("A", "aa1b")]
    text = parms_path.read_text(encoding="utf-8")
    # Rates were looked up under the primary sites, not the full site strings.
    assert "onRate3Dka = 5.0" in text
    assert "offRatekb = 7.0" in text


@pytest.mark.parametrize(
    "reaction",
    [
        ANCILLARY_REACTION,
        "A(aa1f) + A(aa1b) <-> A(aa1f!1).A(aa1b!1)",
    ],
)
def test_auto_time_step_parses_both_reaction_forms(tmp_path, reaction):
    """The time-step regex shares the same site pattern and must accept the
    ancillary form too, otherwise those reactions vanish from the estimate."""
    exporter = _make_exporter(tmp_path)

    dt = exporter._calculate_auto_time_step(
        reactions=[reaction],
        molecule_counts={"A": 100},
        box_nm=(100.0, 100.0, 100.0),
        sigma_list=[1.0],
    )

    assert dt is not None and dt > 0.0
