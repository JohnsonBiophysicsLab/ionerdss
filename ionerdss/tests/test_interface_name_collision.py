"""
Regression tests for interface-type names when molecule template names have
different lengths (e.g. ``A`` next to ``AA``).

Interface names used to be the bare concatenation ``{this}{partner}{index}``,
so the A-side and AA-side of an A-AA interface were both called ``AAA1``: the
second template overwrote the first in the builder's dictionaries, the
exported deck declared ``AA(aaa1) + A(aaa1)`` and ``A.mol`` listed no such
site, and NERDSS exited with code 120 ("aaa1 is not a valid interface for
molecule template A").  Structures hitting this included 2btj, 2vvj, 5exc,
6wat, 6zt0, 7zas, 8bkz, 8bmt and 8q2m.
"""

from pathlib import Path
import re
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import InterfaceType, MoleculeType
from ionerdss.model.pdb import validation
from ionerdss.model.pdb.coarse_graining import InterfaceString
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.template_builder import GeometricSignature, TemplateBuilder
from ionerdss.model.titrate.parms_titrator import parse_mol_file


# --------------------------------------------------------------------------
# Template builder: names created for A/AA templates must stay distinct
# --------------------------------------------------------------------------

def _mock_interface(chain_i="A", chain_j="B"):
    interface = Mock(spec=InterfaceString)
    interface.chain_i = chain_i
    interface.chain_j = chain_j
    interface.coord_i = np.array([5.0, 0.0, 0.0])
    interface.coord_j = np.array([15.0, 0.0, 0.0])
    interface.residues_i = {1, 2}
    interface.residues_j = {3, 4}
    interface.energy = -3.0
    interface.residue_details_i = [Mock(id=1, name="ALA"), Mock(id=2, name="GLY")]
    interface.residue_details_j = [Mock(id=3, name="CYS"), Mock(id=4, name="VAL")]
    return interface


def _builder_with_templates(*template_names):
    builder = TemplateBuilder.__new__(TemplateBuilder)
    builder.workspace_manager = None
    builder.hyperparams = PDBModelHyperparameters()
    builder.parser = SimpleNamespace(convert_coords_to_nm=lambda coords: np.asarray(coords) / 10.0)
    chains = {
        "A": SimpleNamespace(com=np.array([0.0, 0.0, 0.0])),
        "B": SimpleNamespace(com=np.array([20.0, 0.0, 0.0])),
    }
    builder.coarse_grainer = SimpleNamespace(get_coarse_grained_chains=lambda: chains)
    builder.interface_templates = {}
    builder.interface_signatures = {}
    builder.interface_type_counters = {}
    builder.hht_catalog = {}
    builder.molecule_templates = {
        name: SimpleNamespace(interfaces_neighbors_map={}) for name in template_names
    }
    return builder


def _assert_templates_consistent(builder):
    for name, template in builder.interface_templates.items():
        assert template.get_name() == name
        assert template.partner_interface_type is not None
        partner_name = template.partner_interface_type.get_name()
        assert partner_name in builder.interface_templates
        assert builder.interface_templates[partner_name].partner_interface_type is template


def test_heterotypic_templates_for_a_and_aa_do_not_collide():
    builder = _builder_with_templates("A", "AA")
    signature = GeometricSignature(5.0, 5.0, 1.0, 1.2)

    names = builder._create_heterotypic_interface_templates(
        _mock_interface(), "A", "AA", signature, 1
    )

    assert sorted(names) == ["2AAA1", "A2AA1"]
    assert set(builder.interface_templates) == {"A2AA1", "2AAA1"}
    on_a = builder.interface_templates["A2AA1"]
    on_aa = builder.interface_templates["2AAA1"]
    assert (on_a.this_mol_type_name, on_a.partner_mol_type_name) == ("A", "AA")
    assert (on_aa.this_mol_type_name, on_aa.partner_mol_type_name) == ("AA", "A")
    _assert_templates_consistent(builder)
    assert builder.molecule_templates["A"].interfaces_neighbors_map == {"A2AA1": "AA"}
    assert builder.molecule_templates["AA"].interfaces_neighbors_map == {"2AAA1": "A"}


def test_heterotypic_templates_created_from_both_orders_stay_distinct():
    """Creating AA-A after A-AA must add two more templates, not overwrite two."""
    builder = _builder_with_templates("A", "AA")
    signature = GeometricSignature(5.0, 5.0, 1.0, 1.2)

    builder._create_heterotypic_interface_templates(_mock_interface(), "A", "AA", signature, 1)
    builder._create_heterotypic_interface_templates(_mock_interface("B", "A"), "AA", "A", signature, 2)

    assert set(builder.interface_templates) == {"A2AA1", "2AAA1", "2AAA2", "A2AA2"}
    _assert_templates_consistent(builder)


def test_homodimeric_heterotypic_pair_for_multi_letter_template():
    builder = _builder_with_templates("AA")
    signature = GeometricSignature(5.0, 6.0, 1.0, 1.2)

    names = builder._create_heterotypic_interface_templates(
        _mock_interface(), "AA", "AA", signature, 1
    )

    assert sorted(names) == ["2AA2AA1b", "2AA2AA1f"]
    _assert_templates_consistent(builder)
    assert builder.interface_templates["2AA2AA1f"].tag == "f"
    assert builder.interface_templates["2AA2AA1b"].tag == "b"


def test_canonical_hht_pair_for_multi_letter_template():
    """The hand-built ``{name}{name}{idx}f`` names in _ensure_hht_canonical_and_assign."""
    builder = _builder_with_templates("AA")
    builder._match_existing_hht_signature = lambda *args, **kwargs: (None, None)
    builder._sig_tuple = lambda signature: (5.0, 6.0, 1.0, 1.2)
    builder._hht_side_features_from_interface = lambda interface, side: None
    signature = GeometricSignature(5.0, 6.0, 1.0, 1.2)

    name_f = builder._ensure_hht_canonical_and_assign(_mock_interface(), "AA", signature)

    assert name_f == "2AA2AA1f"
    assert set(builder.interface_templates) == {"2AA2AA1f", "2AA2AA1b"}
    _assert_templates_consistent(builder)
    assert builder.hht_catalog[("AA", (5.0, 6.0, 1.0, 1.2))]["b"] == "2AA2AA1b"
    assert set(builder.molecule_templates["AA"].interfaces_neighbors_map) == {"2AA2AA1f", "2AA2AA1b"}


def test_homotypic_template_for_multi_letter_template():
    builder = _builder_with_templates("AA")
    signature = GeometricSignature(5.0, 5.0, 1.0, 1.0)

    name = builder._create_homotypic_interface_template(_mock_interface(), "AA", signature, 1)

    assert name == "2AA2AA1"
    assert builder.interface_templates[name].get_name() == name


# --------------------------------------------------------------------------
# Export: every site named in a reaction must exist in that molecule's .mol
# --------------------------------------------------------------------------

def _molecule_type(name):
    return MoleculeType(
        name=name,
        radius_nm=2.0,
        D_t_nm2_us=1.0,
        D_r_rad2_us=0.1,
        ref1_local=np.array([1.0, 0.0, 0.0]),
        ref2_local=np.array([0.0, 0.0, 1.0]),
    )


def _interface_type(this_type, partner_type, index, local_coord, tag=None):
    return InterfaceType(
        this_mol_type_name=this_type.name,
        partner_mol_type_name=partner_type.name,
        interface_index=index,
        absolute_coord=np.asarray(local_coord, dtype=float),
        local_coord=np.asarray(local_coord, dtype=float),
        energy=-5.0,
        this_mol_type=this_type,
        partner_mol_type=partner_type,
        tag=tag,
    )


def _link(type_a, type_b):
    type_a.partner_interface_type = type_b
    type_b.partner_interface_type = type_a


def _molecule(name, mol_type, com):
    return MoleculeInstance(
        name=name,
        molecule_type=mol_type,
        com=np.asarray(com, dtype=float),
        norm=np.array([0.0, 0.0, 1.0]),
        ref1=np.array([1.0, 0.0, 0.0]),
        ref2=np.array([0.0, 0.0, 1.0]),
    )


def _bind(mol_1, type_1, coord_1, mol_2, type_2, coord_2):
    """Create the two interface instances of one bond and link them to their molecules."""
    inst_1 = InterfaceInstance(
        absolute_coord=np.asarray(coord_1, dtype=float),
        interface_type=type_1,
        this_mol=mol_1,
        this_mol_name=mol_1.molecule_type.name,
        partner_mol_name=mol_2.molecule_type.name,
        interface_index=type_1.interface_index,
    )
    inst_2 = InterfaceInstance(
        absolute_coord=np.asarray(coord_2, dtype=float),
        interface_type=type_2,
        this_mol=mol_2,
        this_mol_name=mol_2.molecule_type.name,
        partner_mol_name=mol_1.molecule_type.name,
        interface_index=type_2.interface_index,
    )
    inst_1.partner_interface = inst_2
    inst_2.partner_interface = inst_1
    mol_1.interfaces_neighbors_map[inst_1] = mol_2
    mol_2.interfaces_neighbors_map[inst_2] = mol_1


def _build_a_aa_system():
    """Two A and two AA molecules: A-AA heterodimer, A-A homotypic, AA-AA f/b."""
    system = System(workspace_path=".")
    type_a = _molecule_type("A")
    type_aa = _molecule_type("AA")
    system.molecule_types.add(type_a)
    system.molecule_types.add(type_aa)

    a_binds_aa = _interface_type(type_a, type_aa, 1, [0.0, 2.0, 0.0])
    aa_binds_a = _interface_type(type_aa, type_a, 1, [0.0, -2.0, 0.0])
    _link(a_binds_aa, aa_binds_a)
    a_binds_a = _interface_type(type_a, type_a, 1, [2.0, 0.0, 0.0])
    a_binds_a.partner_interface_type = a_binds_a
    aa_f = _interface_type(type_aa, type_aa, 1, [0.0, 2.0, 0.0], tag="f")
    aa_b = _interface_type(type_aa, type_aa, 1, [0.0, -2.0, 0.0], tag="b")
    _link(aa_f, aa_b)
    for interface_type in (a_binds_aa, aa_binds_a, a_binds_a, aa_f, aa_b):
        system.interface_types.add(interface_type)

    a_1 = _molecule("A_1", type_a, [0.0, 0.0, 0.0])
    a_2 = _molecule("A_2", type_a, [4.0, 0.0, 0.0])
    aa_1 = _molecule("AA_1", type_aa, [0.0, 6.0, 0.0])
    aa_2 = _molecule("AA_2", type_aa, [0.0, 12.0, 0.0])
    for molecule in (a_1, a_2, aa_1, aa_2):
        system.molecule_instances.add(molecule)

    # The exporter walks each molecule's interfaces_neighbors_map, so the
    # interface instances need not be registered on the system.
    _bind(a_1, a_binds_aa, [0.0, 2.0, 0.0], aa_1, aa_binds_a, [0.0, 4.0, 0.0])
    _bind(a_1, a_binds_a, [2.0, 0.0, 0.0], a_2, a_binds_a, [2.0, 0.0, 0.0])
    _bind(aa_1, aa_f, [0.0, 8.0, 0.0], aa_2, aa_b, [0.0, 10.0, 0.0])
    return system


SPECIES_RE = re.compile(r"([A-Za-z0-9]+)\(([^)]*)\)")


def _reaction_sites(parms_text):
    """Every (molecule, site) named in the reactions block of a parms file."""
    sites = set()
    in_block = False
    for line in parms_text.splitlines():
        stripped = line.strip()
        if stripped == "start reactions":
            in_block = True
            continue
        if stripped == "end reactions":
            in_block = False
            continue
        if not in_block or "->" not in stripped:
            continue
        for molecule, site_group in SPECIES_RE.findall(stripped):
            for site in site_group.split(","):
                site = site.strip().split("!")[0].split("~")[0]
                if site:
                    sites.add((molecule, site))
    return sites


def _binding_reactions(parms_text):
    return [
        line.strip()
        for line in parms_text.splitlines()
        if "<->" in line
    ]


def test_a_and_aa_interface_types_get_distinct_names():
    system = _build_a_aa_system()
    names = sorted(interface_type.get_name() for interface_type in system.interface_types)
    assert names == ["2AA2AA1b", "2AA2AA1f", "2AAA1", "A2AA1", "AA1"]


def test_exported_reactions_only_reference_sites_listed_in_mol_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    system = _build_a_aa_system()

    artifacts = validation.setup_simulation(
        system,
        box_nm=(50.0, 50.0, 50.0),
        initial_molecule_count=1,
        titration_on_rate={"A": 1.0e-5, "AA": 1.0e-5},
    )

    parms_path = Path(artifacts.nerdss_files["parms_titrate"])
    nerdss_dir = parms_path.parent
    parms_text = parms_path.read_text(encoding="utf-8")

    mol_sites = {
        name: set(parse_mol_file(str(nerdss_dir / f"{name}.mol"))) for name in ("A", "AA")
    }
    assert mol_sites["A"] == {"aa1", "a2aa1"}
    assert mol_sites["AA"] == {"2aaa1", "2aa2aa1f", "2aa2aa1b"}

    referenced = _reaction_sites(parms_text)
    assert referenced, "no reactions were exported"
    missing = sorted(
        (molecule, site) for molecule, site in referenced if site not in mol_sites.get(molecule, set())
    )
    assert missing == []

    reactions = _binding_reactions(parms_text)
    assert "A(a2aa1) + AA(2aaa1) <-> A(a2aa1!1).AA(2aaa1!1)" in reactions
    assert "A(aa1) + A(aa1) <-> A(aa1!1).A(aa1!1)" in reactions
    assert "AA(2aa2aa1f) + AA(2aa2aa1b) <-> AA(2aa2aa1f!1).AA(2aa2aa1b!1)" in reactions
    assert len(reactions) == 3
    # the historical collision: both sides of the A-AA bond called "aaa1"
    assert ("A", "aaa1") not in referenced and ("AA", "aaa1") not in referenced
