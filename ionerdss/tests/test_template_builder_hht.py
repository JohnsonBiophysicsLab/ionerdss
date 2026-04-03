from collections import defaultdict
from types import SimpleNamespace

from ionerdss.model.pdb.template_builder import GeometricSignature, TemplateBuilder


def _make_builder() -> TemplateBuilder:
    builder = TemplateBuilder.__new__(TemplateBuilder)
    builder.chain_assigned_types = defaultdict(set)
    builder.hht_catalog = {}
    builder.workspace_manager = None
    return builder


def test_hht_reversed_match_allows_complementary_partner_reuse():
    builder = _make_builder()
    interface = SimpleNamespace(chain_i="A", chain_j="E")
    signature = GeometricSignature(1.0, 2.0, 0.3, 0.7)
    cat = {
        "f": "AA1f",
        "b": "AA1b",
        "exemplar_interface": object(),
        "exemplar_signature": signature,
    }
    builder.hht_catalog[("A", (1.0, 2.0, 0.3, 0.7))] = cat
    builder.chain_assigned_types["A"].add("AA1f")
    builder._homotypic_orientation_candidates = lambda *args, **kwargs: ["ji"]

    match_dir, matched_cat = builder._match_existing_hht_signature("A", signature, interface=interface)

    assert match_dir == "ji"
    assert matched_cat is cat


def test_hht_reversed_match_blocks_exact_duplicate_assignment():
    builder = _make_builder()
    interface = SimpleNamespace(chain_i="A", chain_j="E")
    signature = GeometricSignature(1.0, 2.0, 0.3, 0.7)
    cat = {
        "f": "AA1f",
        "b": "AA1b",
        "exemplar_interface": object(),
        "exemplar_signature": signature,
    }
    builder.hht_catalog[("A", (1.0, 2.0, 0.3, 0.7))] = cat
    builder.chain_assigned_types["A"].add("AA1b")
    builder._homotypic_orientation_candidates = lambda *args, **kwargs: ["ji"]

    match_dir, matched_cat = builder._match_existing_hht_signature("A", signature, interface=interface)

    assert match_dir is None
    assert matched_cat is None


def test_hht_match_falls_back_to_alternate_valid_orientation_when_preferred_conflicts():
    builder = _make_builder()
    interface = SimpleNamespace(chain_i="A", chain_j="E")
    signature = GeometricSignature(1.0, 2.0, 0.3, 0.7)
    cat = {
        "f": "AA1f",
        "b": "AA1b",
        "exemplar_interface": object(),
        "exemplar_signature": signature,
    }
    builder.hht_catalog[("A", (1.0, 2.0, 0.3, 0.7))] = cat
    builder.chain_assigned_types["A"].add("AA1f")
    builder._homotypic_orientation_candidates = lambda *args, **kwargs: ["ij", "ji"]

    match_dir, matched_cat = builder._match_existing_hht_signature("A", signature, interface=interface)

    assert match_dir == "ji"
    assert matched_cat is cat


def test_hht_side_features_prevent_wrong_orientation_fallback():
    builder = _make_builder()
    interface = SimpleNamespace(chain_i="B", chain_j="C")
    signature = GeometricSignature(1.0, 2.0, 0.3, 0.7)
    cat = {
        "f": "AA1f",
        "b": "AA1b",
        "feat_i": {"id": "left"},
        "feat_j": {"id": "right"},
        "exemplar_interface": object(),
        "exemplar_signature": signature,
    }
    builder.hht_catalog[("A", (1.0, 2.0, 0.3, 0.7))] = cat
    builder.chain_assigned_types["B"].add("AA1f")

    def fake_side_features(iface, side):
        return {"id": "left" if side == "i" else "right"}

    builder._hht_side_features_from_interface = fake_side_features
    builder._hht_features_match = lambda a, b: a["id"] == b["id"]
    builder._homotypic_orientation_candidates = lambda *args, **kwargs: ["ij", "ji"]

    match_dir, matched_cat = builder._match_existing_hht_signature("A", signature, interface=interface)

    assert match_dir is None
    assert matched_cat is None
