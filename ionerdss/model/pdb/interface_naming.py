"""
ionerdss/model/pdb/interface_naming.py

Rule: every stage that touches a type name must call
parse_interface_name and must preserve tag.
"""

import re
from dataclasses import dataclass
from typing import Optional

# NERDSS only supports alphanumeric characters (no underscores)
PATTERN = re.compile(r"^(?P<m1>[A-Za-z0-9]+)(?P<m2>[A-Za-z0-9]+)(?P<idx>\d+)(?P<tag>[fb])?$")

@dataclass(frozen=True)
class ParsedName:
    """
    The parsed naming of an interface type
    - heterodimeric interactions ("het"): create two interfaces
     {this_mol}{partner_mol}{index}
     e.g. AB1 and BA1, where AB1 is the interface on A that
     interacts with B.
    - homodimeric heterotypic interactions ("hom_het"): create two interfaces
     {mol}{mol}{index}f and {mol}{mol}{index}b.
      e.g. AA1f and AA1b
    - homodimeric homotypic interactions ("hom_hom"): create
     one interface {mol}{mol}{index}
    """
    this_mol: str
    partner_mol: str
    index: int
    tag: Optional[str]  # 'f' | 'b' | None
    
    def get_type(self):
        """
        return the type as string.
        "het": heterodimeric interactions
        "hom_het": homodimeric heterotypic interactions
        "hom_hom": homodimeric homotypic interactions
        """
        # heterodimeric case
        if self.this_mol != self.partner_mol:
            return "het"
        # homodimeric case
        elif self.tag == 'f' or self.tag == 'b':
            return "hom_hom"
        else:
            return "hom_het"

def parse_interface_name(name: str) -> ParsedName:
    m = PATTERN.match(name)
    if not m:
        raise ValueError(f"Bad interface name: {name}")
    return ParsedName(
        m.group('m1'), m.group('m2'), int(m.group('idx')), m.group('tag')
    )

def make_interface_name(m1: str, m2: str, idx: int, tag: Optional[str]) -> str:
    # No underscores - NERDSS only supports alphanumeric
    base = f"{m1}{m2}{idx}"
    return base if not tag else f"{base}{tag}"

def are_complementary_homodimeric_heterotypic(a: str, b: str) -> bool:
    pa, pb = parse_interface_name(a), parse_interface_name(b)
    return (
        pa.this_mol == pa.partner_mol == pb.this_mol == pb.partner_mol and
        pa.index == pb.index and
        {pa.tag, pb.tag} == {"f","b"}
    )

def are_complementary_heterodimer(a: str, b: str) -> bool:
    pa, pb = parse_interface_name(a), parse_interface_name(b)
    return (pa.this_mol == pb.partner_mol) and (pa.partner_mol == pb.this_mol) and (pa.index == pb.index) and (pa.tag is None and pb.tag is None)
