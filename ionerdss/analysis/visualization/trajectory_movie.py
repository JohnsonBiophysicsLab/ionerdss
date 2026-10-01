"""
Fast movies of NERDSS trajectories, rendered with numpy and Pillow.

Every molecule is drawn as a shaded sphere at its center of mass:

* the camera is fixed on the simulation box (``WaterBox`` in ``parms.inp``),
  which is drawn as a wireframe, so neither the box nor the framing moves
  when molecules do;
* the time label is ``iteration x timeStep`` in one unit chosen for the whole
  movie, printed at a fixed position in a fixed-width field;
* each molecule type has its own color, whether the molecule is free or in a
  complex;
* the sphere radius is the one ionerdss used to compute ``D`` and ``Dr``,
  recovered from the ``.mol`` files by inverting Stokes-Einstein, and drawn
  ``radius_scale`` times larger (3 by default) so molecules stay visible in a
  box hundreds of nanometres wide; ``radius_scale=1`` draws them to scale.

Spheres are rasterized with a vectorized z-buffer, so a frame of ~10^5
molecules takes a fraction of a second on one core, with no OpenGL context,
ray tracer or child process involved.
"""

from __future__ import annotations

import logging
import math
import os
import re
import shutil
import subprocess
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from ionerdss.utils.diffusion_constant import compute_diffusion_constants_nm_us

logger = logging.getLogger(__name__)

# Categorical slots, assigned in this fixed order. The first three stay apart
# under color-vision deficiency in any combination; past that the legend
# carries identity too. Types past the eighth get the extension colors, which
# are distinct but not CVD-checked; pass `colors=` to choose.
TYPE_COLORS: Tuple[str, ...] = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948",
)
_EXTENSION_COLORS: Tuple[str, ...] = (
    "#8c564b", "#17becf", "#7f7f7f", "#bcbd22", "#393b79", "#ad494a",
    "#637939", "#7b4173", "#3182bd", "#e6550d", "#31a354", "#756bb1",
)
_UNKNOWN_COLOR = "#9a9a9a"

_BACKGROUND = (255, 255, 255)
_INK = (30, 30, 30)
_BOX_BACK = (175, 175, 175)
_BOX_FRONT = (90, 90, 90, 150)

# Spheres are drawn this many times their molecular radius unless told otherwise:
# at true size, molecules in a 500 nm box are one or two pixels across.
DEFAULT_RADIUS_SCALE = 3.0

# D and Dr set by ionerdss agree to the 6 significant digits written to the
# .mol file; a hand-written pair is typically off by a large factor.
_STOKES_EINSTEIN_TOLERANCE = 0.05

_TIME_UNITS = (("s", 1e6), ("ms", 1e3), ("µs", 1.0))

_DITHER_NONE = Image.Dither.NONE if hasattr(Image, "Dither") else Image.NONE

_EMPTY_KEY = np.uint64(np.iinfo(np.uint64).max)
_ID_MASK = np.uint64(0xFFFFFFFF)
_DEPTH_LEVELS = float(2**31 - 1)

# GIF frames are held in memory until the file is written.
_GIF_MEMORY_WARNING_BYTES = 1 << 30

_FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


# --------------------------------------------------------------------------
# Simulation inputs
# --------------------------------------------------------------------------

@dataclass
class MoleculeTypeInfo:
    """What the movie needs to know about one molecule type."""

    name: str
    D: Optional[float] = None               # nm^2/us
    Dr: Optional[float] = None              # rad^2/us
    site_distances: List[float] = field(default_factory=list)  # nm

    @property
    def n_sites(self) -> int:
        return len(self.site_distances)


@dataclass
class SimulationInputs:
    """Box, time step and molecule types read from a NERDSS run directory."""

    box_nm: np.ndarray
    time_step_us: Optional[float]
    molecule_types: List[MoleculeTypeInfo]
    missing_mol_files: List[str] = field(default_factory=list)


def _strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def _numbers(text: str) -> List[float]:
    return [float(x) for x in re.findall(_FLOAT, text)]


def _key_value(line: str) -> Optional[Tuple[str, str]]:
    if "=" not in line:
        return None
    key, value = line.split("=", 1)
    # NERDSS lower-cases keys before matching them, so `timeStep`,
    # `timestep` and `TIMESTEP` are the same parameter.
    return key.strip().lower(), value.strip()


def find_parms_file(sim_dir: Union[str, Path]) -> Path:
    """Return ``parms.inp`` in `sim_dir`, or its only ``.inp`` file."""
    sim_dir = Path(sim_dir)
    preferred = sim_dir / "parms.inp"
    if preferred.is_file():
        return preferred
    candidates = sorted(sim_dir.glob("*.inp"))
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"No parms.inp (or other .inp file) in {sim_dir}.")
    raise RuntimeError(
        f"{sim_dir} has several .inp files and no parms.inp "
        f"({', '.join(p.name for p in candidates)}); pass parms_file=."
    )


def read_parms(parms_path: Union[str, Path]) -> Tuple[Optional[np.ndarray], Optional[float], List[str]]:
    """Read ``(WaterBox in nm, timeStep in us, molecule names)`` from a parms file."""
    box = None
    time_step = None
    names: List[str] = []
    section = None
    for raw in Path(parms_path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = _strip_comment(raw)
        if not line:
            continue
        lowered = line.lower()
        if lowered.startswith("start "):
            section = lowered.split(None, 1)[1].strip()
            continue
        if lowered.startswith("end "):
            section = None
            continue
        if section == "molecules":
            name = line.split(":", 1)[0].strip()
            if name:
                names.append(name)
            continue
        kv = _key_value(line)
        if kv is None:
            continue
        key, value = kv
        if key == "timestep":
            time_step = float(_numbers(value)[0])
        elif key == "waterbox":
            box = np.array(_numbers(value)[:3], dtype=float)
    return box, time_step, names


def read_mol_file(mol_path: Union[str, Path]) -> MoleculeTypeInfo:
    """Read name, D, Dr and the COM-to-site distances from a NERDSS ``.mol`` file."""
    mol_path = Path(mol_path)
    info = MoleculeTypeInfo(name=mol_path.stem)
    com = None
    in_coords = False
    for raw in mol_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = _strip_comment(raw)
        if not line:
            continue
        tokens = line.split()
        if in_coords:
            values = _numbers(" ".join(tokens[1:])) if len(tokens) == 4 else []
            if len(values) == 3:
                info.site_distances.append(float(np.linalg.norm(np.array(values) - com)))
                continue
            in_coords = False
        if tokens[0].lower() == "com" and "=" not in line:
            com = np.array(_numbers(" ".join(tokens[1:]))[:3], dtype=float)
            in_coords = True
            continue
        kv = _key_value(line)
        if kv is None:
            continue
        key, value = kv
        values = _numbers(value)
        if key == "name":
            info.name = value
        elif key == "d" and values:
            info.D = float(np.mean(values))
        elif key == "dr" and values:
            info.Dr = float(np.mean(values))
    return info


def read_simulation_inputs(sim_dir: Union[str, Path], parms_file: Optional[Union[str, Path]] = None) -> SimulationInputs:
    """Read the box, time step and molecule types of a NERDSS run directory."""
    sim_dir = Path(sim_dir)
    parms_path = Path(parms_file) if parms_file else find_parms_file(sim_dir)
    box, time_step, names = read_parms(parms_path)
    if box is None:
        raise ValueError(f"No WaterBox in {parms_path}.")

    mol_dir = parms_path.parent
    mol_files = {p.stem: p for p in mol_dir.glob("*.mol")}
    ordered = [mol_files[n] for n in names if n in mol_files]
    if not ordered:
        ordered = [mol_files[k] for k in sorted(mol_files)]
    missing = [n for n in names if n not in mol_files]
    if missing:
        logger.warning("No .mol file for %s in %s; their radius falls back to 1 nm.",
                       ", ".join(missing), mol_dir)
    types = [read_mol_file(p) for p in ordered]
    types += [MoleculeTypeInfo(name=n) for n in missing]
    return SimulationInputs(box_nm=box, time_step_us=time_step, molecule_types=types,
                            missing_mol_files=missing)


def stokes_einstein_radii(D: Optional[float], Dr: Optional[float]) -> Tuple[Optional[float], Optional[float]]:
    """Radii (nm) implied by D and by Dr, inverting ionerdss's Stokes-Einstein relation."""
    D_1nm, Dr_1nm = compute_diffusion_constants_nm_us(radius_nm=1.0)
    r_translational = D_1nm / D if D and D > 0 else None
    r_rotational = (Dr_1nm / Dr) ** (1.0 / 3.0) if Dr and Dr > 0 else None
    return r_translational, r_rotational


def molecule_radius_nm(info: MoleculeTypeInfo) -> Tuple[float, str]:
    """Return ``(radius_nm, source)`` for one molecule type.

    The radius ionerdss used to compute D and Dr when the two agree under
    Stokes-Einstein; otherwise (D and Dr were set by hand) the NERDSS radius,
    i.e. the largest COM-to-interface distance.
    """
    r_t, r_r = stokes_einstein_radii(info.D, info.Dr)
    if r_t and r_r and abs(r_t - r_r) <= _STOKES_EINSTEIN_TOLERANCE * max(r_t, r_r):
        return r_t, "D, Dr (Stokes-Einstein)"
    if info.site_distances and max(info.site_distances) > 0:
        radius = max(info.site_distances)
        if r_t and r_r:
            logger.warning(
                "%s: D and Dr are not a Stokes-Einstein pair (they imply %.3g nm and "
                "%.3g nm), so they were set by hand; drawing the NERDSS radius (largest "
                "COM-to-interface distance, %.3g nm) instead. Pass radii= to override.",
                info.name, r_t, r_r, radius,
            )
        return radius, "COM-interface distance"
    if r_t:
        return r_t, "D (Stokes-Einstein)"
    if r_r:
        return r_r, "Dr (Stokes-Einstein)"
    return 1.0, "default"


# --------------------------------------------------------------------------
# Trajectory readers
# --------------------------------------------------------------------------

@dataclass
class Frame:
    """Molecule centers of one frame, in nm with the box centered on the origin."""

    iteration: int
    labels: np.ndarray       # molecule type label per molecule (as written by NERDSS)
    positions: np.ndarray    # (N, 3)


def list_pdb_frames(pdb_dir: Union[str, Path]) -> List[Tuple[int, Path]]:
    """``(iteration, path)`` of every ``<iteration>.pdb`` in `pdb_dir`, in order."""
    frames = [(int(p.stem), p) for p in Path(pdb_dir).glob("*.pdb") if p.stem.isdigit()]
    return sorted(frames)


# Column layout of the ATOM records written by NERDSS's write_pdb.cpp: name
# (COM) at 12-15, residue name (molecule type) at 17-19, then %8.1f x, y, z.
_PDB_COORD_COLUMNS = np.array([29, 37, 45])
_PDB_DECIMAL_POINTS = _PDB_COORD_COLUMNS + 6


def _read_pdb_com_fixed_width(data: bytes) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Vectorized COM extraction for NERDSS's fixed-width ATOM records.

    Returns None when any COM record deviates from that layout (e.g. files
    from older NERDSS versions), so the caller can parse line by line.
    """
    if not data.endswith(b"\n"):
        data += b"\n"
    buf = np.frombuffer(data, dtype=np.uint8)
    ends = np.flatnonzero(buf == ord("\n"))
    starts = np.concatenate(([0], ends[:-1] + 1))
    starts = starts[ends - starts >= 53]
    com = starts[(buf[starts] == ord("A")) & (buf[starts + 3] == ord("M"))
                 & (buf[starts + 12] == ord("C")) & (buf[starts + 13] == ord("O"))
                 & (buf[starts + 14] == ord("M")) & (buf[starts + 15] == ord(" "))]
    if com.size == 0:
        return np.zeros(0, dtype=object), np.zeros((0, 3))
    layout_ok = ((buf[com + 16] == ord(" ")).all() and (buf[com + 20] == ord(" ")).all()
                 and all((buf[com + d] == ord(".")).all() for d in _PDB_DECIMAL_POINTS))
    if not layout_ok:
        return None
    fields = buf[com[:, None, None] + _PDB_COORD_COLUMNS[None, :, None] + np.arange(8)[None, None, :]]
    positions = np.ascontiguousarray(fields).view("S8").reshape(-1, 3).astype(float)
    names = np.ascontiguousarray(buf[com[:, None] + np.arange(17, 20)]).view("S3").ravel()
    unique, inverse = np.unique(names, return_inverse=True)
    labels = np.array([u.decode(errors="replace").strip() for u in unique], dtype=object)[inverse]
    return labels, positions


def read_pdb_frame(path: Union[str, Path], iteration: int, box_nm: np.ndarray) -> Frame:
    """Read the COM atoms of one NERDSS PDB snapshot.

    NERDSS writes the molecule type (first 3 characters) as the residue name
    and shifts coordinates by half the box; the shift is undone here.
    """
    with open(path, "rb") as handle:
        parsed = _read_pdb_com_fixed_width(handle.read())
    if parsed is not None:
        labels, positions = parsed
        return Frame(iteration, labels, positions - box_nm / 2.0)

    labels = []
    coords = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            # The atom name sits in columns 13-16 in every NERDSS version.
            if line[12:15] != "COM" or not line.startswith("ATOM"):
                continue
            parts = line.split()
            labels.append(parts[3])
            coords.append(parts[5:8])
    positions = np.array(coords, dtype=float).reshape(-1, 3) - box_nm / 2.0
    return Frame(iteration, np.array(labels, dtype=object), positions)


def scan_xyz_iterations(path: Union[str, Path]) -> List[int]:
    """Iteration numbers of every frame in a NERDSS ``trajectory.xyz``."""
    pattern = re.compile(rb"^iteration:\s*(\d+)", re.M)
    iterations = []
    with open(path, "rb") as handle:
        tail = b""
        while True:
            chunk = handle.read(1 << 24)
            if not chunk:
                break
            block = tail + chunk
            cut = block.rfind(b"\n") + 1
            iterations += [int(m) for m in pattern.findall(block[:cut])]
            tail = block[cut:]
        iterations += [int(m) for m in pattern.findall(tail)]
    return iterations


def iter_xyz_frames(
    path: Union[str, Path],
    sites_per_label: Mapping[str, int],
    wanted: Optional[set] = None,
) -> Iterator[Frame]:
    """Yield the molecule centers of each frame of a NERDSS ``trajectory.xyz``.

    NERDSS writes each molecule's COM followed by one line per interface, all
    labelled with the molecule type, so COM lines are recovered from the
    per-type interface counts of the ``.mol`` files. The atom count in each
    frame header is not trusted: frames are delimited by ``iteration:`` lines.
    """
    def parse(iteration: int, body: List[str]) -> Frame:
        labels = []
        coords = []
        i = 0
        n = len(body)
        while i < n:
            parts = body[i].split()
            if len(parts) < 4:
                i += 1
                continue
            label = parts[0]
            if label not in sites_per_label:
                raise ValueError(
                    f"{path}: molecule type '{label}' at iteration {iteration} has no "
                    "matching .mol file, so its COM lines cannot be told apart from its "
                    "interface lines. Use source='pdb' or add the .mol file."
                )
            labels.append(label)
            coords.append(parts[1:4])
            i += 1 + sites_per_label[label]
        positions = np.array(coords, dtype=float).reshape(-1, 3)
        return Frame(iteration, np.array(labels, dtype=object), positions)

    iteration = None
    body: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("iteration:"):
                if iteration is not None and (wanted is None or iteration in wanted):
                    # The last body line is the next frame's atom count.
                    yield parse(iteration, body[:-1] if body and body[-1].strip().isdigit() else body)
                iteration = int(line.split(":", 1)[1])
                body = []
            elif iteration is not None:
                body.append(line)
    if iteration is not None and (wanted is None or iteration in wanted):
        yield parse(iteration, body)


# --------------------------------------------------------------------------
# Labels, colors, camera
# --------------------------------------------------------------------------

def _hex_to_rgb(color: Union[str, Sequence[float]]) -> np.ndarray:
    if isinstance(color, str):
        text = color.lstrip("#")
        if len(text) != 6:
            raise ValueError(f"Colors must be '#rrggbb' strings or RGB tuples, got {color!r}.")
        return np.array([int(text[i:i + 2], 16) for i in (0, 2, 4)], dtype=float) / 255.0
    rgb = np.array(color, dtype=float)[:3]
    return rgb / 255.0 if rgb.max() > 1.0 else rgb


def time_label_format(times_us: Sequence[float], unit: Optional[str] = None) -> Tuple[str, float, int, int]:
    """Choose ``(unit, us_per_unit, decimals, width)`` for the whole movie.

    One unit for every frame (the largest in which the final time is at least
    1) and the fewest decimals, up to 3, that show the usual spacing between
    frames exactly, or, when no such number exists (an irregular timeStep),
    that still tell neighbouring frames apart. `width` fits the longest value,
    so the label never changes length.
    """
    times = np.asarray(times_us, dtype=float)
    t_max = float(times.max()) if times.size else 0.0
    units = dict(_TIME_UNITS)
    if unit is not None:
        key = "µs" if unit in ("us", "μs") else unit
        if key not in units:
            raise ValueError(f"time_unit must be one of 's', 'ms', 'us', got {unit!r}.")
        unit_name, scale = key, units[key]
    else:
        unit_name, scale = next(((n, s) for n, s in _TIME_UNITS if t_max >= s), _TIME_UNITS[-1])
    values = times / scale
    # NERDSS writes its last frame off the regular grid (e.g. at nItr - 1), so
    # take the fewest decimals that show at least half of the frame spacings
    # exactly, rather than every value.
    steps = np.diff(np.unique(values))
    if steps.size == 0:
        steps = values[:1]
    decimals = None
    for d in range(4):
        exact = np.abs(np.round(steps, d) - steps) <= 1e-6 * np.maximum(1.0, np.abs(steps))
        if steps.size == 0 or exact.mean() >= 0.5:
            decimals = d
            break
    if decimals is None:
        # An irregular spacing, e.g. from an automatically chosen timeStep:
        # enough decimals to tell neighbouring frames apart.
        spacing = float(np.median(steps))
        decimals = min(3, max(0, math.ceil(-math.log10(spacing)))) if spacing > 0 else 3
    width = max(len(f"{v:.{decimals}f}") for v in values) if values.size else 1
    return unit_name, scale, decimals, width


def _rotation(azimuth_deg: float, elevation_deg: float) -> np.ndarray:
    """Rows are the screen right, screen up, and toward-the-viewer axes."""
    az = math.radians(azimuth_deg)
    el = math.radians(elevation_deg)
    toward_viewer = np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])
    right = np.array([-math.sin(az), math.cos(az), 0.0])
    up = np.cross(toward_viewer, right)
    return np.vstack([right, up, toward_viewer])


def _load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    """DejaVu Sans Mono ships with matplotlib, so digits keep a fixed width."""
    try:
        import matplotlib

        name = "DejaVuSansMono-Bold.ttf" if bold else "DejaVuSansMono.ttf"
        path = Path(matplotlib.get_data_path()) / "fonts" / "ttf" / name
        return ImageFont.truetype(str(path), size)
    except Exception:
        try:
            return ImageFont.load_default(size=size)
        except TypeError:  # Pillow < 10.1
            return ImageFont.load_default()


# --------------------------------------------------------------------------
# Renderer
# --------------------------------------------------------------------------

# Key light from the upper left, in front of the scene (screen right, up,
# toward the viewer), and the Blinn half-vector for its highlight.
_LIGHT = np.array([-0.45, 0.55, 0.70]) / np.linalg.norm([-0.45, 0.55, 0.70])
_HALF = (_LIGHT + [0.0, 0.0, 1.0]) / np.linalg.norm(_LIGHT + [0.0, 0.0, 1.0])


class TrajectoryRenderer:
    """Renders frames of one movie with a fixed camera, box, legend and label.

    Parameters:
        box_nm: Simulation box edge lengths (``WaterBox``), in nm.
        type_names: Molecule type names in legend/color order.
        radii_nm: Molecular radius per type name, in nm.
        colors: Color per type name; defaults to the categorical slots in order.
        labels: Map from the label NERDSS writes in a trajectory to a type name.
        size: Frame size in pixels (rounded up to even numbers for video codecs).
        view: Camera ``(azimuth, elevation)`` in degrees.
        time_step_us: ``timeStep`` in us; without it the label shows iterations.
        iterations: Every iteration that will be rendered, so the label can pick
            one unit and one width for the whole movie.
        show_legend: Also list each type's color (and molecular radius) next
            to the time.
        time_font_size: Height of the time label in pixels; about 6.5% of the
            smaller frame side by default.
        radius_scale: Every sphere is drawn this many times its molecular
            radius (3 by default); 1 draws molecules to scale.
        supersample: Spheres are rasterized at this multiple of the frame size
            and averaged down, which smooths their edges.
    """

    def __init__(
        self,
        box_nm: Sequence[float],
        type_names: Sequence[str],
        radii_nm: Mapping[str, float],
        *,
        colors: Optional[Mapping[str, Union[str, Sequence[float]]]] = None,
        labels: Optional[Mapping[str, str]] = None,
        size: Tuple[int, int] = (800, 800),
        view: Tuple[float, float] = (35.0, 25.0),
        time_step_us: Optional[float] = None,
        iterations: Optional[Sequence[int]] = None,
        time_unit: Optional[str] = None,
        show_box: bool = True,
        show_time: bool = True,
        show_legend: bool = False,
        show_radius_in_legend: bool = True,
        time_font_size: Optional[int] = None,
        radius_scale: float = DEFAULT_RADIUS_SCALE,
        supersample: int = 2,
    ) -> None:
        if radius_scale <= 0:
            raise ValueError(f"radius_scale must be positive, got {radius_scale}.")
        self.box = np.asarray(box_nm, dtype=float)
        self.radius_scale = float(radius_scale)
        self.width, self.height = (int(size[0]) + int(size[0]) % 2, int(size[1]) + int(size[1]) % 2)
        self.type_names = list(type_names)
        self.show_box = show_box
        self.show_time = show_time
        self.time_step_us = time_step_us
        self.ss = max(1, int(supersample))

        n_types = len(self.type_names)
        colors = dict(colors or {})
        palette = list(TYPE_COLORS) + list(_EXTENSION_COLORS)
        if n_types > len(TYPE_COLORS) and any(n not in colors for n in self.type_names[len(TYPE_COLORS):]):
            logger.warning(
                "%d molecule types: colors past the first %d are distinct but not "
                "validated for color-vision deficiency. Pass colors= to choose them.",
                n_types, len(TYPE_COLORS),
            )
        # The last row is for labels that match no known type.
        rgb = [
            _hex_to_rgb(colors.get(name, palette[i] if i < len(palette) else _UNKNOWN_COLOR))
            for i, name in enumerate(self.type_names)
        ]
        self.type_rgb = np.vstack(rgb + [_hex_to_rgb(_UNKNOWN_COLOR)])
        self._type_rgb32 = self.type_rgb.astype(np.float32)
        known_radii = [float(radii_nm[name]) for name in self.type_names]
        self.radii_nm = np.array(known_radii + [float(np.median(known_radii)) if known_radii else 1.0])
        # Molecular radii go in the legend; spheres are drawn at the scaled ones.
        self.drawn_radii_nm = self.radii_nm * self.radius_scale
        self.label_to_type = {name: i for i, name in enumerate(self.type_names)}
        for label, name in (labels or {}).items():
            if name in self.label_to_type:
                self.label_to_type[label] = self.label_to_type[name]
        self._unknown_labels: set = set()
        self._dropped_warned = False

        self.rotation = _rotation(*view)
        self.font_size = max(12, round(0.028 * min(self.width, self.height)))
        self.font = _load_font(self.font_size)
        self.time_font_size = int(time_font_size) if time_font_size else max(16, round(0.065 * min(self.width, self.height)))
        self.time_font = _load_font(self.time_font_size, bold=True)
        self._label_format = None
        if show_time:
            last = max(iterations) if iterations else 0
            if time_step_us is not None and iterations:
                self._label_format = time_label_format(
                    [it * time_step_us for it in iterations], time_unit)
            self._time_example = self.time_text(last * (time_step_us or 0.0), iteration=last)

        self._layout(show_legend, show_radius_in_legend)
        self._build_sprites()
        self._build_static_layers()
        self.palette_image = self._build_gif_palette()

    # -- layout ------------------------------------------------------------

    def time_text(self, time_us: Optional[float], iteration: int) -> str:
        if self._label_format is None:
            return f"step {iteration:,}"
        unit, scale, decimals, width = self._label_format
        return f"t = {time_us / scale:>{width}.{decimals}f} {unit}"

    def _layout(self, show_legend: bool, show_radius: bool) -> None:
        W, H = self.width, self.height
        fs = self.font_size
        pad = round(0.3 * self.time_font_size) if self.show_time else round(0.6 * fs)
        line_h = round(1.45 * fs)
        time_h = round(1.25 * self.time_font_size) if self.show_time else 0
        draw = ImageDraw.Draw(Image.new("RGB", (8, 8)))

        # The legend flows to the right of the time label and wraps under itself.
        self.legend_entries: List[Tuple[int, int, int, str]] = []
        x_start = pad
        if self.show_time:
            x_start += int(draw.textlength(self._time_example, font=self.time_font)) + 2 * fs
        rows = 0
        if show_legend and self.type_names:
            x, row = x_start, 0
            disc = round(0.8 * fs)
            for i, name in enumerate(self.type_names):
                text = f"{name} (r = {self.radii_nm[i]:.3g} nm)" if show_radius else name
                item_w = disc + round(0.4 * fs) + int(draw.textlength(text, font=self.font)) + fs
                if x + item_w > W - pad and x > x_start:
                    x, row = x_start, row + 1
                self.legend_entries.append((i, x, pad + row * line_h, text))
                x += item_w
            rows = row + 1
        self.line_h = line_h
        self.time_h = time_h
        self.pad = pad
        content = max(time_h, rows * line_h)
        band = 2 * pad + content if content else 0

        # Fit the projected box into the area under the band, leaving room for
        # spheres that sit on its faces.
        half = self.box / 2.0
        corners = np.array([[sx, sy, sz] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]) * half
        self._corners = corners
        proj = corners @ self.rotation[:2].T
        extent = proj.max(axis=0) - proj.min(axis=0)
        max_r = float(self.drawn_radii_nm.max())
        margin = 0.03 * min(W, H)
        for _ in range(3):
            avail_w = W - 2 * margin
            avail_h = H - band - 2 * margin
            if avail_w <= 0 or avail_h <= 0:
                raise ValueError(f"size={W, H} leaves no room for the box; use a larger size.")
            scale = min(avail_w / extent[0], avail_h / extent[1])
            margin = max(0.03 * min(W, H), max_r * scale + 2)
        self.scale = scale
        mid = (proj.max(axis=0) + proj.min(axis=0)) / 2.0
        self.cx = W / 2.0 - scale * mid[0]
        self.cy = band + (H - band) / 2.0 + scale * mid[1]

        self.radii_px = np.maximum(self.drawn_radii_nm * scale, 1.0)
        small = [n for n, r in zip(self.type_names, self.radii_px) if r < 2.0]
        if small:
            logger.warning(
                "Spheres of %s are under 2 px in radius at this size; increase size= or "
                "use radius_scale= to make them easier to see.", ", ".join(small),
            )

        # Depth runs from the nearest to the farthest point of any sphere in the box.
        depth_half = float(np.abs(self.rotation[2]) @ half) + max_r
        self.depth_min = -depth_half
        self.depth_span = 2 * depth_half

    # -- sprites -----------------------------------------------------------

    def _build_sprites(self) -> None:
        """Pixel offsets and depth bumps of each distinct sphere radius, supersampled."""
        ss = self.ss
        # Supersampled pixel j covers frame coordinates [(j - (ss-1)/2 - 0.5)/ss, ...),
        # so the ss x ss block of a frame pixel averages back onto it.
        self._scale_ss = self.scale * ss
        self._cx_ss = self.cx * ss + (ss - 1) / 2.0
        self._cy_ss = self.cy * ss + (ss - 1) / 2.0
        self.radii_ss = self.radii_px * ss
        self.pad_px = int(math.ceil(self.radii_ss.max())) + 2
        self.buf_w = self.width * ss + 2 * self.pad_px
        self.buf_h = self.height * ss + 2 * self.pad_px
        self.sprites: Dict[float, Tuple[np.ndarray, np.ndarray]] = {}
        # Lighting depends only on a pixel's offset from its sphere's center, so
        # it is tabulated once per radius: table[start + (dy + n) * (2n + 1) + dx + n].
        shade_tables, spec_tables = [], []
        table_start = {}
        start = 0
        for r in np.unique(self.radii_ss):
            n = int(math.ceil(r))
            dy, dx = np.mgrid[-n:n + 1, -n:n + 1]
            inside = dx * dx + dy * dy <= r * r
            offsets = (dy[inside] * self.buf_w + dx[inside]).astype(np.int64)
            bump_nm = np.sqrt(np.maximum(r * r - dx[inside] ** 2 - dy[inside] ** 2, 0.0)) / self._scale_ss
            self.sprites[float(r)] = (offsets, bump_nm)

            ux, uy = (dx / r).ravel(), (dy / r).ravel()
            rho2 = np.minimum(ux * ux + uy * uy, 1.0)
            nz = np.sqrt(1.0 - rho2)
            diffuse = np.clip(ux * _LIGHT[0] - uy * _LIGHT[1] + nz * _LIGHT[2], 0.0, 1.0)
            shade = np.minimum(0.35 + 0.75 * diffuse, 1.0)
            # A one-pixel darker rim separates touching spheres of the same color.
            if r >= 3.0:
                shade = np.where(rho2 > ((r - 1.0) / r) ** 2, shade * 0.6, shade)
            spec = np.clip(ux * _HALF[0] - uy * _HALF[1] + nz * _HALF[2], 0.0, 1.0) ** 40
            shade_tables.append(shade.astype(np.float32))
            spec_tables.append((0.55 * spec).astype(np.float32))
            table_start[float(r)] = (start, n)
            start += (2 * n + 1) ** 2
        self._shade_table = np.concatenate(shade_tables)
        self._spec_table = np.concatenate(spec_tables)
        starts_n = np.array([table_start[float(r)] for r in self.radii_ss], dtype=np.int64)
        self._table_start = starts_n[:, 0]
        self._table_n = starts_n[:, 1]

    # -- static layers -----------------------------------------------------

    def _project(self, points: np.ndarray) -> np.ndarray:
        cam = points @ self.rotation.T
        return np.column_stack([self.cx + self.scale * cam[:, 0], self.cy - self.scale * cam[:, 1]])

    def _build_static_layers(self) -> None:
        """Background (with the hidden box edges) and overlay (visible edges, legend)."""
        W, H = self.width, self.height
        background = Image.new("RGB", (W, H), _BACKGROUND)
        self.overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        if self.show_box:
            corners = self._corners
            xy = self._project(corners)
            far = int(np.argmin(corners @ self.rotation[2]))
            back_draw = ImageDraw.Draw(background)
            front_draw = ImageDraw.Draw(self.overlay)
            line_w = max(1, round(self.font_size / 12))
            for a in range(8):
                for b in range(a + 1, 8):
                    # Box edges join corners that differ in exactly one coordinate.
                    if np.count_nonzero(corners[a] != corners[b]) != 1:
                        continue
                    segment = [tuple(xy[a]), tuple(xy[b])]
                    # The three edges at the far corner are behind the contents.
                    if far in (a, b):
                        back_draw.line(segment, fill=_BOX_BACK, width=line_w)
                    else:
                        front_draw.line(segment, fill=_BOX_FRONT, width=line_w)
        draw = ImageDraw.Draw(self.overlay)
        disc = round(0.8 * self.font_size)
        for i, x, y, text in self.legend_entries:
            top = y + (self.line_h - disc) // 2
            color = tuple(int(round(c * 255)) for c in self.type_rgb[i])
            draw.ellipse([x, top, x + disc, top + disc], fill=color + (255,),
                         outline=tuple(int(c * 0.6) for c in color) + (255,))
            draw.text((x + disc + round(0.4 * self.font_size), y + self.line_h // 2), text,
                      font=self.font, fill=_INK + (255,), anchor="lm")

        self._background_u8 = np.asarray(background, dtype=np.uint8).reshape(-1, 3).copy()
        self._background_f = self._background_u8.astype(np.float64) / 255.0
        overlay = np.asarray(self.overlay, dtype=np.float64).reshape(-1, 4)
        self._overlay_index = np.flatnonzero(overlay[:, 3] > 0)
        alpha = overlay[self._overlay_index, 3:] / 255.0
        self._overlay_alpha = alpha
        self._overlay_premultiplied = overlay[self._overlay_index, :3] * alpha

    @staticmethod
    def _shade_rgb(base: np.ndarray, shade: np.ndarray, highlight: np.ndarray) -> np.ndarray:
        """Dark rim -> base color at full light -> white highlight."""
        shade = np.asarray(shade, dtype=np.float32)
        highlight = np.asarray(highlight, dtype=np.float32)
        lit = base * (0.25 + 0.75 * np.clip(shade, 0.0, 1.0))[..., None]
        return lit + (1.0 - lit) * np.clip(highlight, 0.0, 1.0)[..., None]

    def _build_gif_palette(self) -> Image.Image:
        """One palette for every frame, so GIF colors never flicker.

        Per type: the shading ramp, the highlight ramp, and blends with the
        background for antialiased edges; plus a grey ramp for box and text.
        """
        bg = np.array(_BACKGROUND, dtype=np.float32) / 255.0
        ink = np.array(_INK, dtype=np.float32) / 255.0
        entries = [bg + (ink - bg) * a for a in np.linspace(0.0, 1.0, 16)]
        per_type = max(6, (256 - 16) // len(self.type_rgb))
        n_fade = per_type * 3 // 10
        n_high = per_type // 5
        n_shade = per_type - n_fade - n_high
        for base in self.type_rgb.astype(np.float32):
            entries += list(self._shade_rgb(base[None, :], np.linspace(0.0, 1.0, n_shade), np.zeros(n_shade)))
            entries += list(self._shade_rgb(base[None, :], np.ones(n_high), np.linspace(0.1, 0.6, n_high)))
            mid = self._shade_rgb(base[None, :], np.array([0.6]), np.array([0.0]))[0]
            entries += [mid + (bg - mid) * a for a in np.linspace(0.2, 0.85, n_fade)]
        flat = np.clip(np.rint(np.array(entries[:256]) * 255), 0, 255).astype(np.uint8).ravel().tolist()
        palette_image = Image.new("P", (1, 1))
        palette_image.putpalette(flat + [0] * (768 - len(flat)))
        return palette_image

    # -- per frame ---------------------------------------------------------

    def type_indices(self, labels: Sequence[str]) -> np.ndarray:
        unknown = len(self.type_names)
        lookup = self.label_to_type.get
        out = np.fromiter((lookup(label, unknown) for label in labels), dtype=np.int64, count=len(labels))
        if (out == unknown).any():
            for label in set(np.asarray(labels, dtype=object)[out == unknown]) - self._unknown_labels:
                self._unknown_labels.add(label)
                logger.warning("Molecule type '%s' is not in the .mol files; drawing it grey.", label)
        return out

    def rasterize(self, positions: np.ndarray, types: np.ndarray) -> np.ndarray:
        """Z-buffer of packed ``(depth << 32) | molecule index`` over the padded frame."""
        zbuf = np.full(self.buf_w * self.buf_h, _EMPTY_KEY, dtype=np.uint64)
        self._px = self._py = self._far = np.zeros(0)
        if len(positions) == 0:
            return zbuf
        cam = positions @ self.rotation.T
        px = np.rint(self._cx_ss + self._scale_ss * cam[:, 0]).astype(np.int64) + self.pad_px
        py = np.rint(self._cy_ss - self._scale_ss * cam[:, 1]).astype(np.int64) + self.pad_px
        depth = -cam[:, 2]
        radius = self.radii_ss[types]
        reach = np.ceil(radius).astype(np.int64)
        visible = ((px - reach >= 0) & (px + reach < self.buf_w)
                   & (py - reach >= 0) & (py + reach < self.buf_h))
        if not visible.all() and not self._dropped_warned:
            self._dropped_warned = True
            logger.warning("%d molecules lie too far outside the box to draw; they are skipped.",
                           int((~visible).sum()))
        self._px, self._py = px, py
        self._far = np.clip((depth - self.depth_min) / self.depth_span, 0.0, 1.0).astype(np.float32)
        for r, (offsets, bump_nm) in self.sprites.items():
            members = np.flatnonzero(visible & (radius == r))
            if members.size == 0:
                continue
            base = py[members] * self.buf_w + px[members]
            ids = members.astype(np.uint64)
            d0 = depth[members]
            for off, bump in zip(offsets, bump_nm):
                q = np.clip((d0 - bump - self.depth_min) / self.depth_span, 0.0, 1.0)
                keys = ((q * _DEPTH_LEVELS).astype(np.uint64) << np.uint64(32)) | ids
                np.minimum.at(zbuf, base + off, keys)
        return zbuf

    def shade(self, zbuf: np.ndarray, types: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Frame pixels touched by spheres, their premultiplied color and their coverage.

        Only covered supersampled pixels are shaded and averaged down, so the
        cost follows how much of the frame the spheres fill.
        """
        ss, W, H = self.ss, self.width, self.height
        hit = np.flatnonzero(zbuf != _EMPTY_KEY)
        y, x = np.divmod(hit, self.buf_w)
        inside = ((x >= self.pad_px) & (x < self.pad_px + W * ss)
                  & (y >= self.pad_px) & (y < self.pad_px + H * ss))
        if not inside.all():
            hit, x, y = hit[inside], x[inside], y[inside]
        if hit.size == 0:
            return np.zeros(0, dtype=np.int64), np.zeros((0, 3)), np.zeros(0)
        mol = (zbuf[hit] & _ID_MASK).astype(np.int64)
        t = types[mol]
        n = self._table_n[t]
        entry = self._table_start[t] + (y - self._py[mol] + n) * (2 * n + 1) + (x - self._px[mol] + n)
        # Farther molecules are darker, which reads as depth.
        far = self._far[mol]
        shade = self._shade_table[entry] * (1.0 - 0.3 * far)
        highlight = self._spec_table[entry] * (1.0 - 0.5 * far)
        rgb = self._shade_rgb(self._type_rgb32[t], shade, highlight)

        # Average each ss x ss block onto its frame pixel.
        x -= self.pad_px
        y -= self.pad_px
        pixel = (y // ss) * W + (x // ss)
        n_pixels = W * H
        weight = 1.0 / (ss * ss)
        cover = np.bincount(pixel, minlength=n_pixels) * weight
        touched = np.flatnonzero(cover)
        color = np.column_stack([
            np.bincount(pixel, weights=rgb[:, c], minlength=n_pixels)[touched] * weight for c in range(3)
        ])
        return touched, color, cover[touched]

    def render(self, positions: np.ndarray, labels: Sequence[str], iteration: int = 0) -> Image.Image:
        """Render one frame. `positions` are molecule centers in nm, box centered on 0."""
        types = self.type_indices(labels)
        zbuf = self.rasterize(np.asarray(positions, dtype=float).reshape(-1, 3), types)
        touched, color, cover = self.shade(zbuf, types)
        out = self._background_u8.copy()
        out[touched] = np.rint(255.0 * (color + (1.0 - cover)[:, None] * self._background_f[touched]))
        # Visible box edges and the legend go over the spheres.
        o = self._overlay_index
        out[o] = np.rint(out[o] * (1.0 - self._overlay_alpha) + self._overlay_premultiplied)
        frame = Image.fromarray(out.reshape(self.height, self.width, 3), "RGB")
        if self.show_time:
            time_us = iteration * self.time_step_us if self.time_step_us is not None else None
            ImageDraw.Draw(frame).text(
                (self.pad, self.pad + self.time_h // 2), self.time_text(time_us, iteration),
                font=self.time_font, fill=_INK, anchor="lm",
            )
        return frame


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------

def _find_ffmpeg() -> Optional[str]:
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg

        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        return None


class _GifWriter:
    def __init__(self, path: Path, fps: float, loop: int, palette_image: Image.Image) -> None:
        self.path, self.fps, self.loop, self.palette_image = path, fps, loop, palette_image
        self.frames: List[Image.Image] = []
        self._bytes = 0
        self._warned = False

    def append(self, frame: Image.Image) -> None:
        self.frames.append(frame.quantize(palette=self.palette_image, dither=_DITHER_NONE))
        self._bytes += frame.width * frame.height
        if self._bytes > _GIF_MEMORY_WARNING_BYTES and not self._warned:
            self._warned = True
            logger.warning("GIF frames already take over 1 GB of memory; an .mp4 output or "
                           "frame_stride= would keep it down.")

    def close(self) -> None:
        if not self.frames:
            raise ValueError("No frames were rendered.")
        self.frames[0].save(
            self.path, save_all=True, append_images=self.frames[1:],
            duration=max(20, int(round(1000.0 / self.fps))), loop=self.loop, optimize=False,
        )


class _Mp4Writer:
    def __init__(self, path: Path, fps: float, size: Tuple[int, int], crf: int = 18) -> None:
        exe = _find_ffmpeg()
        if exe is None:
            raise RuntimeError(
                "Writing .mp4 needs ffmpeg: install it (e.g. `conda install ffmpeg` or "
                "`pip install imageio-ffmpeg`), or write a .gif instead."
            )
        self.path = path
        self.proc = subprocess.Popen(
            [exe, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
             "-s", f"{size[0]}x{size[1]}", "-r", str(fps), "-i", "-",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", str(crf), str(path)],
            stdin=subprocess.PIPE, stderr=subprocess.PIPE,
        )

    def append(self, frame: Image.Image) -> None:
        try:
            self.proc.stdin.write(frame.tobytes())
        except BrokenPipeError:
            self.close()

    def close(self) -> None:
        if self.proc.stdin and not self.proc.stdin.closed:
            try:
                self.proc.stdin.close()
            except BrokenPipeError:
                pass
        stderr = self.proc.stderr.read().decode(errors="replace") if self.proc.stderr else ""
        if self.proc.wait() != 0:
            raise RuntimeError(f"ffmpeg failed writing {self.path}:\n{stderr.strip()}")


# --------------------------------------------------------------------------
# Frame loop (serial or across worker processes)
# --------------------------------------------------------------------------

# A task is ("pdb", iteration, path) or ("arrays", iteration, (positions, labels)).
Task = Tuple[str, int, object]

_WORKER_STATE: Dict[str, object] = {}


def _render_task(renderer: TrajectoryRenderer, task: Task, box_nm: np.ndarray,
                 png_path: Optional[Path]) -> Image.Image:
    kind, iteration, payload = task
    if kind == "pdb":
        frame = read_pdb_frame(payload, iteration, box_nm)
        positions, labels = frame.positions, frame.labels
    else:
        positions, labels = payload
    image = renderer.render(positions, labels, iteration)
    if png_path is not None:
        image.save(png_path, compress_level=1)
    return image


def _init_worker(renderer_kwargs: dict, box_nm: np.ndarray) -> None:
    # The parent already reported what building the renderer has to say.
    logger.setLevel(logging.ERROR)
    _WORKER_STATE["renderer"] = TrajectoryRenderer(**renderer_kwargs)
    _WORKER_STATE["box"] = box_nm
    logger.setLevel(logging.NOTSET)


def _worker_render(task: Task, png_path: Optional[Path], want_image: bool) -> Optional[bytes]:
    image = _render_task(_WORKER_STATE["renderer"], task, _WORKER_STATE["box"], png_path)
    return image.tobytes() if want_image else None


def _render_frames(
    renderer: TrajectoryRenderer,
    renderer_kwargs: dict,
    tasks: Iterator[Task],
    box_nm: np.ndarray,
    frames_dir: Optional[Path],
    want_images: bool,
    n_jobs: int,
) -> Iterator[Optional[Image.Image]]:
    """Yield rendered frames in order, rendering up to `n_jobs` at once."""
    def png_path(index: int) -> Optional[Path]:
        return frames_dir / f"frame_{index:05d}.png" if frames_dir is not None else None

    if n_jobs == 1:
        for index, task in enumerate(tasks):
            yield _render_task(renderer, task, box_nm, png_path(index))
        return

    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures.process import BrokenProcessPool

    size = (renderer.width, renderer.height)
    # spawn: forking a process that holds threads (Jupyter, BLAS) is unsafe.
    pool = ProcessPoolExecutor(n_jobs, mp_context=multiprocessing.get_context("spawn"),
                               initializer=_init_worker, initargs=(renderer_kwargs, box_nm))
    pending: deque = deque()
    yielded = 0

    def next_result() -> Optional[Image.Image]:
        # Read before dequeuing: if the pool broke, this frame's task must stay
        # queued for the serial fallback below.
        data = pending[0][0].result()
        pending.popleft()
        return Image.frombytes("RGB", size, data) if data is not None else None

    try:
        # Keep only a few frames in flight so memory stays flat for long runs.
        for index, task in enumerate(tasks):
            pending.append((pool.submit(_worker_render, task, png_path(index), want_images), index, task))
            if len(pending) >= 2 * n_jobs:
                image = next_result()
                yielded += 1
                yield image
        while pending:
            image = next_result()
            yielded += 1
            yield image
    except BrokenProcessPool:
        pool.shutdown(wait=False, cancel_futures=True)
        if yielded:
            raise
        # Typically a script without an `if __name__ == "__main__":` guard, or
        # one read from stdin, which spawned workers cannot import.
        logger.warning("Worker processes could not start; rendering in this process instead.")
        remaining = [(i, task) for _, i, task in pending] if pending else []
        pending.clear()
        next_index = remaining[-1][0] + 1 if remaining else 0
        for index, task in remaining + list(enumerate(tasks, start=next_index)):
            yield _render_task(renderer, task, box_nm, png_path(index))
        return
    except BaseException:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown()


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def _match_labels(type_names: Sequence[str], width: int) -> Dict[str, str]:
    """Map the truncated labels NERDSS writes back to full type names."""
    by_label: Dict[str, List[str]] = {}
    for name in type_names:
        by_label.setdefault(name[:width], []).append(name)
    clashes = {label: names for label, names in by_label.items() if len(names) > 1}
    if clashes:
        detail = "; ".join(f"{label}: {', '.join(names)}" for label, names in clashes.items())
        logger.warning(
            "NERDSS truncates molecule names to %d characters in this output, so these "
            "types share a label and share a color: %s.", width, detail,
        )
    return {label: names[0] for label, names in by_label.items()}


def render_trajectory_movie(
    sim_dir: Union[str, Path] = ".",
    output: Optional[Union[str, Path]] = "trajectory.gif",
    *,
    source: str = "auto",
    parms_file: Optional[Union[str, Path]] = None,
    size: Tuple[int, int] = (800, 800),
    fps: float = 10,
    loop: int = 0,
    frame_stride: int = 1,
    view: Tuple[float, float] = (35.0, 25.0),
    radii: Optional[Union[float, Mapping[str, float]]] = None,
    radius_scale: float = DEFAULT_RADIUS_SCALE,
    colors: Optional[Mapping[str, Union[str, Sequence[float]]]] = None,
    time_unit: Optional[str] = None,
    show_box: bool = True,
    show_time: bool = True,
    show_legend: bool = False,
    time_font_size: Optional[int] = None,
    supersample: int = 2,
    frames_dir: Optional[Union[str, Path]] = None,
    n_jobs: int = 1,
    progress: bool = True,
) -> Optional[Path]:
    """Render a NERDSS run as a GIF or MP4 movie of molecule-center spheres.

    Parameters:
        sim_dir: NERDSS run directory (holds ``parms.inp``, the ``.mol`` files and
            ``PDB/`` or ``DATA/trajectory.xyz``).
        output: ``.gif`` or ``.mp4`` path (MP4 needs ffmpeg), or None to only
            write PNG frames to `frames_dir`.
        source: ``"pdb"`` (``PDB/<iteration>.pdb``), ``"xyz"``
            (``DATA/trajectory.xyz``), or ``"auto"``: PDB when present.
        parms_file: Parameter file, when it is not ``sim_dir/parms.inp``.
        size: Frame size in pixels (rounded up to even numbers).
        fps: Frames per second of the movie.
        loop: GIF loop count; 0 loops forever.
        frame_stride: Render every n-th saved frame.
        view: Camera ``(azimuth, elevation)`` in degrees.
        radii: Molecular radius in nm, for all types or per type name. By
            default the radius ionerdss used for D and Dr (see
            `molecule_radius_nm`).
        radius_scale: Every sphere is drawn this many times its molecular
            radius; 3 by default so molecules stay visible in a large box, 1
            draws them to scale. The legend always lists molecular radii.
        colors: Color per type name (``"#rrggbb"`` or RGB).
        time_unit: ``"us"``, ``"ms"`` or ``"s"``; chosen from the final time by default.
        show_box, show_time: Toggle the box wireframe and the time label.
        show_legend: Also list each type's color and radius next to the time.
        time_font_size: Height of the time label in pixels; about 6.5% of the
            smaller frame side by default.
        supersample: Rasterize spheres at this multiple of `size` for smooth
            edges; 1 is fastest.
        frames_dir: Also write each frame as a PNG here.
        n_jobs: Frames rendered at once in worker processes; -1 uses every
            core. Workers take a second or two to start, so this pays off for
            large systems or long runs. In a script, call this function under
            ``if __name__ == "__main__":`` when n_jobs is not 1.
        progress: Show a progress bar.

    Returns:
        The movie path, or None when only PNG frames were written.
    """
    sim_dir = Path(sim_dir)
    if output is None and frames_dir is None:
        raise ValueError("Give an output movie path, frames_dir, or both.")
    if frame_stride < 1:
        raise ValueError("frame_stride must be at least 1.")
    output_path = Path(output) if output is not None else None
    if output_path is not None and output_path.suffix.lower() not in (".gif", ".mp4"):
        raise ValueError(f"output must end in .gif or .mp4, got {output_path.name}.")

    inputs = read_simulation_inputs(sim_dir, parms_file)
    type_names = [t.name for t in inputs.molecule_types]
    if isinstance(radii, (int, float)):
        radius_by_name = {name: float(radii) for name in type_names}
    else:
        radius_by_name = {}
        for info in inputs.molecule_types:
            if radii is not None and info.name in radii:
                radius_by_name[info.name] = float(radii[info.name])
                logger.info("%s: radius %.3g nm (given)", info.name, radius_by_name[info.name])
            else:
                r, how = molecule_radius_nm(info)
                radius_by_name[info.name] = r
                logger.info("%s: radius %.3g nm (%s)", info.name, r, how)

    pdb_dir = sim_dir / "PDB"
    xyz_path = sim_dir / "DATA" / "trajectory.xyz"
    if source == "auto":
        source = "pdb" if pdb_dir.is_dir() and list_pdb_frames(pdb_dir) else "xyz"
    if source == "pdb":
        pdb_frames = list_pdb_frames(pdb_dir)
        if not pdb_frames:
            raise FileNotFoundError(f"No <iteration>.pdb files in {pdb_dir}.")
        pdb_frames = pdb_frames[::frame_stride]
        iterations = [it for it, _ in pdb_frames]
        labels = _match_labels(type_names, 3)
        tasks: Iterator[Task] = iter([("pdb", it, str(p)) for it, p in pdb_frames])
    elif source == "xyz":
        if not xyz_path.is_file():
            raise FileNotFoundError(f"No PDB frames or trajectory.xyz under {sim_dir}.")
        iterations = scan_xyz_iterations(xyz_path)[::frame_stride]
        labels = _match_labels(type_names, 4)
        no_mol = [t.name for t in inputs.molecule_types if t.name in inputs.missing_mol_files]
        if no_mol:
            raise FileNotFoundError(
                f"trajectory.xyz cannot be read without the .mol file of {', '.join(no_mol)}: "
                "it holds the interface count that separates COM lines from interface "
                "lines. Use source='pdb' or restore the .mol files."
            )
        sites = {label: next(t.n_sites for t in inputs.molecule_types if t.name == name)
                 for label, name in labels.items()}
        tasks = (("arrays", f.iteration, (f.positions, f.labels))
                 for f in iter_xyz_frames(xyz_path, sites, wanted=set(iterations)))
    else:
        raise ValueError(f"source must be 'auto', 'pdb' or 'xyz', got {source!r}.")
    if not iterations:
        raise ValueError(f"No frames found for source '{source}' in {sim_dir}.")

    if inputs.time_step_us is None:
        logger.warning("No timeStep in the parameter file; the label shows iterations.")
    renderer_kwargs = dict(
        box_nm=inputs.box_nm, type_names=type_names, radii_nm=radius_by_name,
        colors=colors, labels=labels, size=size, view=view,
        time_step_us=inputs.time_step_us, iterations=iterations, time_unit=time_unit,
        show_box=show_box, show_time=show_time, show_legend=show_legend,
        time_font_size=time_font_size, radius_scale=radius_scale, supersample=supersample,
    )
    renderer = TrajectoryRenderer(**renderer_kwargs)
    if n_jobs < 1:
        n_jobs = os.cpu_count() or 1
    n_jobs = max(1, min(n_jobs, len(iterations)))

    writer = None
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix.lower() == ".gif":
            writer = _GifWriter(output_path, fps, loop, renderer.palette_image)
        else:
            writer = _Mp4Writer(output_path, fps, (renderer.width, renderer.height))
    if frames_dir is not None:
        frames_dir = Path(frames_dir)
        frames_dir.mkdir(parents=True, exist_ok=True)

    bar = None
    if progress:
        try:
            from tqdm.auto import tqdm

            bar = tqdm(total=len(iterations), desc="Rendering frames", unit="frame")
        except Exception:
            bar = None
    images = _render_frames(renderer, renderer_kwargs, tasks, inputs.box_nm, frames_dir,
                            writer is not None, n_jobs)
    try:
        for image in images:
            if writer is not None:
                writer.append(image)
            if bar is not None:
                bar.update(1)
    except BaseException:
        if isinstance(writer, _Mp4Writer):
            writer.proc.kill()
        raise
    finally:
        if bar is not None:
            bar.close()
    if writer is not None:
        writer.close()
        logger.info("Saved %d frames to %s", len(iterations), output_path.resolve())
        return output_path
    return None
