"""
Unit tests for ionerdss.model.components.interface_naming.

Interface-type names are alphanumeric (NERDSS site names cannot contain
underscores) and must decode unambiguously even when molecule type names have
different lengths or contain digits.  Single-letter molecule names keep the
historical spelling (AB1, AA1f); any other name is written as a
length-prefixed token (AA -> 2AA).
"""

import numpy as np
import pytest

from ionerdss.model.components import interface_naming as naming
from ionerdss.model.components.types import InterfaceType


# (this_mol, partner_mol, index, tag, expected name)
ROUND_TRIP_CASES = [
    ("A", "B", 1, None, "AB1"),
    ("A", "A", 1, "f", "AA1f"),
    ("A", "A", 1, "b", "AA1b"),
    ("A", "A", 2, None, "AA2"),
    ("A", "B", 0, None, "AB0"),
    ("A", "A", 10, None, "AA10"),
    ("A", "B", 123, None, "AB123"),
    ("A", "AA", 1, None, "A2AA1"),
    ("AA", "A", 1, None, "2AAA1"),
    ("AA", "AA", 1, "f", "2AA2AA1f"),
    ("AA", "AA", 1, "b", "2AA2AA1b"),
    ("AAAA", "A", 1, None, "4AAAAA1"),
    ("AA0", "A", 3, None, "3AA0A3"),
    ("A", "AA0", 1, None, "A3AA01"),
    ("AA0", "AA1", 2, None, "3AA03AA12"),
    ("1", "A", 1, None, "11A1"),
    ("A", "1", 1, None, "A111"),
    ("f", "b", 1, "f", "fb1f"),
    ("Ab", "aB", 0, None, "2Ab2aB0"),
    ("cube", "cube", 3, None, "4cube4cube3"),
    ("ABCDEFGHI", "A", 1, None, "9ABCDEFGHIA1"),
    ("ABCDEFGHIJ", "A", 1, None, "010ABCDEFGHIJA1"),
    ("Dodecahedron", "Dodecahedron", 12, "b", "012Dodecahedron012Dodecahedron12b"),
]


@pytest.mark.parametrize("m1,m2,idx,tag,expected", ROUND_TRIP_CASES)
def test_make_and_parse_round_trip(m1, m2, idx, tag, expected):
    name = naming.make_interface_name(m1, m2, idx, tag)
    assert name == expected
    assert name.isalnum()

    parsed = naming.parse_interface_name(name)
    assert (parsed.this_mol, parsed.partner_mol, parsed.index, parsed.tag) == (m1, m2, idx, tag)
    assert parsed.name == name


def test_single_letter_names_keep_historical_spelling():
    assert naming.make_interface_name("A", "B", 1, None) == "AB1"
    assert naming.make_interface_name("A", "A", 1, "f") == "AA1f"
    assert naming.make_interface_name("A", "A", 1, None) == "AA1"


def test_names_of_different_length_molecules_do_not_collide():
    """A-AA and AA-A both used to be spelled "AAA1"."""
    pairs = [("A", "AA"), ("AA", "A"), ("A", "AAA"), ("AAA", "A"), ("AA", "AA"), ("A", "A")]
    names = [naming.make_interface_name(m1, m2, 1, None) for m1, m2 in pairs]
    assert len(set(names)) == len(pairs)
    for (m1, m2), name in zip(pairs, names):
        parsed = naming.parse_interface_name(name)
        assert (parsed.this_mol, parsed.partner_mol) == (m1, m2)


def test_index_digits_are_not_confused_with_molecule_names():
    assert naming.parse_interface_name("AA10") == naming.ParsedName("A", "A", 10, None)
    assert naming.parse_interface_name("A3AA01") == naming.ParsedName("A", "AA0", 1, None)
    assert naming.parse_interface_name("3AA0A3") == naming.ParsedName("AA0", "A", 3, None)


@pytest.mark.parametrize(
    "bad",
    [
        "",
        "A",           # no partner
        "AB",          # no index
        "A1",          # length prefix with nothing after it
        "AB1x",        # unknown tag
        "AB1ff",       # doubled tag
        "AB1f2",       # tag before digits
        "A_B_1",       # legacy underscore form
        "2A1",         # token "A1" then nothing
        "9AB1",        # length prefix longer than the name
        "0AB1",        # two-digit length prefix without digits
        "01AB1",       # two-digit length prefix with a letter in it
        "009ABCDEFGHI1",  # two-digit length below ten
        "2AAAB",       # no index after the tokens
        "AB-1",
        "ABé1",
    ],
)
def test_parse_rejects_malformed_names(bad):
    with pytest.raises(ValueError):
        naming.parse_interface_name(bad)


def test_parse_rejects_non_strings():
    with pytest.raises(ValueError):
        naming.parse_interface_name(None)


@pytest.mark.parametrize(
    "m1,m2,idx,tag",
    [
        ("A_1", "B", 1, None),   # underscore
        ("", "B", 1, None),      # empty
        ("A B", "B", 1, None),   # space
        ("A", "B", -1, None),    # negative index
        ("A", "B", 1, "x"),      # unknown tag
        ("A" * 100, "B", 1, None),  # too long for a two-digit length
    ],
)
def test_make_rejects_bad_inputs(m1, m2, idx, tag):
    with pytest.raises(ValueError):
        naming.make_interface_name(m1, m2, idx, tag)


def test_make_accepts_empty_tag_as_none():
    assert naming.make_interface_name("A", "B", 1, "") == "AB1"


@pytest.mark.parametrize(
    "name,expected",
    [
        ("AB1", "BA1"),
        ("2AAA1", "A2AA1"),
        ("A2AA1", "2AAA1"),
        ("AA1f", "AA1b"),
        ("AA1b", "AA1f"),
        ("2AA2AA3f", "2AA2AA3b"),
        ("AA1", "AA1"),
    ],
)
def test_complementary_interface_name(name, expected):
    assert naming.complementary_interface_name(name) == expected
    assert naming.complementary_interface_name(expected) == name


@pytest.mark.parametrize(
    "name,expected_type",
    [
        ("AB1", "het"),
        ("2AAA1", "het"),
        ("AA1f", "hom_het"),
        ("AA1b", "hom_het"),
        ("AA1", "hom_hom"),
        ("2AA2AA1", "hom_hom"),
    ],
)
def test_get_type(name, expected_type):
    assert naming.parse_interface_name(name).get_type() == expected_type


def test_are_complementary_homodimeric_heterotypic():
    assert naming.are_complementary_homodimeric_heterotypic("AA1f", "AA1b")
    assert naming.are_complementary_homodimeric_heterotypic("AA1b", "AA1f")
    assert naming.are_complementary_homodimeric_heterotypic("2AA2AA3b", "2AA2AA3f")
    assert not naming.are_complementary_homodimeric_heterotypic("AA1f", "AA2b")
    assert not naming.are_complementary_homodimeric_heterotypic("AA1f", "AA1f")
    assert not naming.are_complementary_homodimeric_heterotypic("AA1f", "AA1")
    assert not naming.are_complementary_homodimeric_heterotypic("AA1f", "BB1b")


def test_are_complementary_heterodimer():
    assert naming.are_complementary_heterodimer("AB1", "BA1")
    assert naming.are_complementary_heterodimer("A2AA1", "2AAA1")
    assert not naming.are_complementary_heterodimer("AB1", "BA2")
    assert not naming.are_complementary_heterodimer("AB1", "AB1")
    assert not naming.are_complementary_heterodimer("AA1f", "AA1b")


@pytest.mark.parametrize(
    "name,label",
    [
        ("AB1", "ab1"),
        ("AA1f", "aa1f"),
        ("AA1", "aa1"),
        ("A2AA1", "a2aa1"),
        ("2AAA1", "2aaa1"),
        ("2AA2AA1b", "2aa2aa1b"),
        ("4cube4cube1", "4cube4cube1"),
    ],
)
def test_make_site_label(name, label):
    site = naming.make_site_label(name)
    assert site == label
    assert site.isalnum()
    # the label decodes to the same interface, lower-cased
    parsed = naming.parse_interface_name(site)
    original = naming.parse_interface_name(name)
    assert parsed == naming.ParsedName(
        original.this_mol.lower(), original.partner_mol.lower(), original.index, original.tag
    )


def test_site_labels_on_one_molecule_stay_distinct():
    """AA/A index 1 and AA/AA index 1 used to be spelled the same way as A/AA."""
    labels = {
        naming.make_site_label(naming.make_interface_name(m1, m2, idx, tag))
        for m1, m2, idx, tag in [
            ("AA", "A", 1, None), ("AA", "AA", 1, None), ("AA", "AA", 1, "f"),
            ("AA", "AA", 1, "b"), ("AA", "A", 11, None), ("AA", "A1", 1, None),
        ]
    }
    assert len(labels) == 6


def _interface(this_mol, partner_mol, index, tag=None):
    return InterfaceType(
        this_mol_type_name=this_mol,
        partner_mol_type_name=partner_mol,
        interface_index=index,
        absolute_coord=np.zeros(3),
        local_coord=np.zeros(3),
        tag=tag,
    )


def test_interface_type_get_name_uses_the_scheme():
    assert _interface("A", "B", 1).get_name() == "AB1"
    assert _interface("A", "A", 1, "f").get_name() == "AA1f"
    assert _interface("A", "AA", 1).get_name() == "A2AA1"
    assert _interface("AA", "A", 1).get_name() == "2AAA1"
    assert _interface("AA", "AA", 2, "b").get_name() == "2AA2AA2b"


def test_interface_type_get_name_rejects_non_alphanumeric_molecule_names():
    with pytest.raises(ValueError):
        _interface("Protein_A", "B", 1).get_name()


@pytest.mark.parametrize(
    "name,expected",
    [
        ("2AAA1", ("AA", "A", 1, None)),
        ("A2AA1", ("A", "AA", 1, None)),
        ("AA1f", ("A", "A", 1, "f")),
        ("AB7", ("A", "B", 7, None)),
        ("A_A_2b", ("A", "A", 2, "b")),          # legacy underscore form
        ("NewMol_NewPartner_5", ("NewMol", "NewPartner", 5, None)),
    ],
)
def test_interface_type_set_name_accepts_current_and_legacy_forms(name, expected):
    intf = _interface("X", "Y", 9)
    intf.set_name(name)
    assert (intf.this_mol_type_name, intf.partner_mol_type_name, intf.interface_index, intf.tag) == expected


def test_interface_type_set_name_round_trips_through_get_name():
    intf = _interface("X", "Y", 9)
    for name in ("2AAA1", "AA1f", "AB7", "012DodecahedronA1"):
        intf.set_name(name)
        assert intf.get_name() == name


def test_pdb_module_re_exports_the_components_scheme():
    from ionerdss.model.pdb import interface_naming as pdb_naming

    assert pdb_naming.make_interface_name is naming.make_interface_name
    assert pdb_naming.parse_interface_name is naming.parse_interface_name
    assert pdb_naming.are_complementary_heterodimer is naming.are_complementary_heterodimer
    assert (
        pdb_naming.are_complementary_homodimeric_heterotypic
        is naming.are_complementary_homodimeric_heterotypic
    )
