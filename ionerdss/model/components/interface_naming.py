"""
ionerdss.model.components.interface_naming
==========================================

Build and parse the names of interface types.

An interface type is identified by the molecule type it sits on, the molecule
type it binds, an index that distinguishes several interface types between the
same pair of molecule types, and an optional ``f``/``b`` tag that marks the two
halves of a homodimeric heterotypic pair:

======================================  ==============================  ================
interaction                             name                            example
======================================  ==============================  ================
heterodimeric (``het``)                 ``{this}{partner}{index}``      ``AB1`` / ``BA1``
homodimeric heterotypic (``hom_het``)   ``{mol}{mol}{index}f`` / ``b``  ``AA1f`` / ``AA1b``
homodimeric homotypic (``hom_hom``)     ``{mol}{mol}{index}``           ``AA1``
======================================  ==============================  ================

NERDSS site names may only contain letters and digits, so the parts are
concatenated without a separator.  To keep that concatenation decodable when
molecule names have different lengths (``A`` next to ``AA``) or contain digits
(``AA0``), every molecule name is written as a self-delimiting *token*:

* a single letter is written as is: ``A`` -> ``A``;
* any other name is prefixed by its length: ``AA`` -> ``2AA``,
  ``AA0`` -> ``3AA0``, ``1`` -> ``11``; a leading ``0`` introduces a two-digit
  length for names of ten to ninety-nine characters:
  ``Dodecahedron`` -> ``012Dodecahedron``.

A token that starts with a letter is exactly one character long and a token
that starts with a digit announces its own length, so the parser always knows
where each part ends and what remains is ``{index}{tag}``.  ``A``/``AA`` and
``AA``/``A`` therefore get the distinct names ``A2AA1`` and ``2AAA1`` instead
of both collapsing to ``AAA1``, and an index of ten or more can no longer be
mistaken for part of a molecule name (``AA10`` is A/A index 10).  Names built
from single-letter molecule types keep their historical spelling.

The NERDSS site label of an interface type is the same spelling with the
molecule names in lower case (``AB1`` -> ``ab1``, ``2AAA1`` -> ``2aaa1``); see
:func:`make_site_label`.

Rule: every stage that touches an interface-type name must go through
:func:`make_interface_name` / :func:`parse_interface_name` and must preserve
the tag.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Tuple

__all__ = [
    "ParsedName",
    "TAGS",
    "are_complementary_heterodimer",
    "are_complementary_homodimeric_heterotypic",
    "complementary_interface_name",
    "make_interface_name",
    "make_site_label",
    "parse_interface_name",
]

TAGS = ("f", "b")

_MAX_NAME_LENGTH = 99
_REMAINDER = re.compile(r"(?P<index>\d+)(?P<tag>[fb])?")


def _check_molecule_name(mol_name: str) -> str:
    if not isinstance(mol_name, str) or not mol_name:
        raise ValueError("molecule type name must be a non-empty string")
    if not (mol_name.isascii() and mol_name.isalnum()):
        raise ValueError(
            f"molecule type name {mol_name!r} must contain only ASCII letters and "
            "digits (NERDSS names cannot contain '_' or other symbols)"
        )
    if len(mol_name) > _MAX_NAME_LENGTH:
        raise ValueError(
            f"molecule type name {mol_name!r} is longer than {_MAX_NAME_LENGTH} characters"
        )
    return mol_name


def _encode_token(mol_name: str) -> str:
    """Write a molecule type name as a self-delimiting token."""
    mol_name = _check_molecule_name(mol_name)
    if len(mol_name) == 1 and mol_name.isalpha():
        return mol_name
    if len(mol_name) < 10:
        return f"{len(mol_name)}{mol_name}"
    return f"0{len(mol_name)}{mol_name}"


def _decode_token(name: str, pos: int) -> Tuple[str, int]:
    """Read the token starting at ``name[pos]``; return (molecule name, next position)."""
    head = name[pos:pos + 1]
    if head.isalpha():
        return head, pos + 1
    if head == "0":
        length_digits = name[pos + 1:pos + 3]
        if len(length_digits) != 2 or not length_digits.isdigit() or length_digits[0] == "0":
            raise ValueError(f"Bad interface name: {name!r}")
        length, start = int(length_digits), pos + 3
    elif head.isdigit():
        length, start = int(head), pos + 1
    else:
        raise ValueError(f"Bad interface name: {name!r}")
    token = name[start:start + length]
    if len(token) != length or not token.isalnum():
        raise ValueError(f"Bad interface name: {name!r}")
    return token, start + length


@dataclass(frozen=True)
class ParsedName:
    """The parts of an interface-type name (see the module docstring for the scheme).

    - heterodimeric ("het"): ``{this_mol}{partner_mol}{index}`` on each side,
      e.g. ``AB1`` on A binds ``BA1`` on B.
    - homodimeric heterotypic ("hom_het"): ``{mol}{mol}{index}f`` binds
      ``{mol}{mol}{index}b``, e.g. ``AA1f`` and ``AA1b``.
    - homodimeric homotypic ("hom_hom"): one self-binding interface
      ``{mol}{mol}{index}``, e.g. ``AA1``.
    """

    this_mol: str
    partner_mol: str
    index: int
    tag: Optional[str]  # 'f' | 'b' | None

    def get_type(self) -> str:
        """Classify the interaction as ``het``, ``hom_het`` (tagged) or ``hom_hom``."""
        if self.this_mol != self.partner_mol:
            return "het"
        if self.tag in TAGS:
            return "hom_het"
        return "hom_hom"

    @property
    def name(self) -> str:
        """The interface-type name spelled by these parts."""
        return make_interface_name(self.this_mol, self.partner_mol, self.index, self.tag)

    def complementary(self) -> "ParsedName":
        """The parts of the interface type this one binds to."""
        if self.this_mol != self.partner_mol:
            return ParsedName(self.partner_mol, self.this_mol, self.index, self.tag)
        if self.tag in TAGS:
            return ParsedName(self.this_mol, self.partner_mol, self.index,
                              "b" if self.tag == "f" else "f")
        return self


def make_interface_name(m1: str, m2: str, idx: int, tag: Optional[str] = None) -> str:
    """Spell the name of the interface type on molecule ``m1`` that binds ``m2``.

    Args:
        m1: Name of the molecule type carrying the interface.
        m2: Name of the molecule type it binds.
        idx: Non-negative interface index.
        tag: ``"f"`` or ``"b"`` for the two halves of a homodimeric heterotypic
            pair, otherwise ``None``.

    Raises:
        ValueError: If a molecule name is not ASCII alphanumeric, the index is
            negative, or the tag is not ``f``/``b``/``None``.
    """
    if tag in (None, ""):
        tag = None
    elif tag not in TAGS:
        raise ValueError(f"interface tag must be 'f', 'b' or None, got {tag!r}")
    index = int(idx)
    if index < 0:
        raise ValueError(f"interface index must be non-negative, got {idx!r}")
    return f"{_encode_token(m1)}{_encode_token(m2)}{index}{tag or ''}"


def parse_interface_name(name: str) -> ParsedName:
    """Split an interface-type name back into its parts.

    Raises:
        ValueError: If ``name`` does not follow the scheme.
    """
    if not isinstance(name, str) or not name.isascii():
        raise ValueError(f"Bad interface name: {name!r}")
    this_mol, pos = _decode_token(name, 0)
    partner_mol, pos = _decode_token(name, pos)
    match = _REMAINDER.fullmatch(name, pos)
    if match is None:
        raise ValueError(f"Bad interface name: {name!r}")
    return ParsedName(this_mol, partner_mol, int(match.group("index")), match.group("tag"))


def complementary_interface_name(name: str) -> str:
    """Name of the interface type that ``name`` binds to.

    ``AB1`` -> ``BA1``, ``2AAA1`` -> ``A2AA1``, ``AA1f`` -> ``AA1b``, ``AA1`` -> ``AA1``.
    """
    return parse_interface_name(name).complementary().name


def make_site_label(name: str) -> str:
    """NERDSS site label for an interface type: its name with lower-case molecule names.

    ``AB1`` -> ``ab1``, ``AA1f`` -> ``aa1f``, ``A2AA1`` -> ``a2aa1``, ``2AAA1`` -> ``2aaa1``.
    The label keeps the self-delimiting tokens, so two interface types on the
    same molecule can never share a site label.
    """
    parsed = parse_interface_name(name)
    return make_interface_name(parsed.this_mol.lower(), parsed.partner_mol.lower(),
                               parsed.index, parsed.tag)


def are_complementary_homodimeric_heterotypic(a: str, b: str) -> bool:
    """True if ``a`` and ``b`` are the ``f``/``b`` halves of one homodimeric heterotypic pair."""
    pa, pb = parse_interface_name(a), parse_interface_name(b)
    return (
        pa.this_mol == pa.partner_mol == pb.this_mol == pb.partner_mol
        and pa.index == pb.index
        and {pa.tag, pb.tag} == {"f", "b"}
    )


def are_complementary_heterodimer(a: str, b: str) -> bool:
    """True if ``a`` and ``b`` are the two untagged sides of one heterodimeric interface."""
    pa, pb = parse_interface_name(a), parse_interface_name(b)
    return (
        pa.this_mol == pb.partner_mol
        and pa.partner_mol == pb.this_mol
        and pa.index == pb.index
        and pa.tag is None
        and pb.tag is None
    )
