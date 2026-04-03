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


def test_heterotypic_duplicate_check_uses_concrete_partner_type_per_chain():
    builder = _make_builder()
    builder.hyperparams = SimpleNamespace(
        interface_type_assignment_distance_threshold=0.1,
        interface_type_assignment_angle_threshold=0.1,
    )
    ab_sig = GeometricSignature(1.0, 2.0, 0.3, 0.7)
    ba_sig = ab_sig.flipped()
    builder.interface_signatures = {
        "AB1": ab_sig,
        "BA1": ba_sig,
    }
    builder.interface_templates = {
        "AB1": SimpleNamespace(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            partner_interface_type=SimpleNamespace(name="BA1"),
        ),
        "BA1": SimpleNamespace(
            this_mol_type_name="B",
            partner_mol_type_name="A",
            partner_interface_type=SimpleNamespace(name="AB1"),
        ),
    }

    # After assigning an A-B contact, chain B should be tracked as carrying BA1.
    builder.chain_assigned_types["A"].add("AB1")
    builder.chain_assigned_types["B"].add("BA1")
    interface = SimpleNamespace(chain_i="B", chain_j="E")

    match = builder._find_matching_interface_type("B", "A", ba_sig, interface)

    assert match is None
