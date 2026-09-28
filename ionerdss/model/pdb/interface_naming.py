"""
ionerdss/model/pdb/interface_naming.py

Compatibility re-export of :mod:`ionerdss.model.components.interface_naming`.

The naming scheme lives in the components package so that
``InterfaceType.get_name`` can use it without importing the PDB pipeline.
See that module for the scheme.

Rule: every stage that touches a type name must call
parse_interface_name and must preserve tag.
"""

from ionerdss.model.components.interface_naming import (  # noqa: F401
    ParsedName,
    TAGS,
    are_complementary_heterodimer,
    are_complementary_homodimeric_heterotypic,
    complementary_interface_name,
    make_interface_name,
    make_site_label,
    parse_interface_name,
)

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
