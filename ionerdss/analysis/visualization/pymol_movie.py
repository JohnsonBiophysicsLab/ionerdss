"""
PyMOL movie export helpers for numeric-frame NERDSS PDB trajectories.
"""

from __future__ import annotations

import ast
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

from PIL import Image, ImageDraw, ImageFont


DEFAULT_TYPE_COLORS: Sequence[Tuple[str, Tuple[float, float, float]]] = (
    ("ionerdss_type_01", (0.894, 0.102, 0.110)),
    ("ionerdss_type_02", (0.216, 0.494, 0.722)),
    ("ionerdss_type_03", (0.302, 0.686, 0.290)),
    ("ionerdss_type_04", (0.596, 0.306, 0.639)),
    ("ionerdss_type_05", (1.000, 0.498, 0.000)),
    ("ionerdss_type_06", (1.000, 1.000, 0.200)),
    ("ionerdss_type_07", (0.651, 0.337, 0.157)),
    ("ionerdss_type_08", (0.969, 0.506, 0.749)),
    ("ionerdss_type_09", (0.600, 0.600, 0.600)),
    ("ionerdss_type_10", (0.090, 0.745, 0.812)),
    ("ionerdss_type_11", (0.121, 0.470, 0.705)),
    ("ionerdss_type_12", (0.682, 0.780, 0.910)),
)

TIMESTEP_PATTERN = re.compile(
    r"^\s*timestep\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)


def get_font(size: int) -> ImageFont.ImageFont:
    """Try common bold sans fonts before falling back to Pillow's default."""
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size)
        except Exception:
            continue
    return ImageFont.load_default()


def find_numeric_pdb_files(pdb_dir: Path) -> list[Path]:
    """Return numeric PDB files sorted by their integer stem."""
    pdb_files = sorted(
        [path for path in Path(pdb_dir).glob("*.pdb") if path.stem.isdigit()],
        key=lambda path: int(path.stem),
    )
    if not pdb_files:
        raise FileNotFoundError(
            f"No PDB files with numeric filenames were found in {Path(pdb_dir)}."
        )
    return pdb_files


def parse_timestep_width_us(base_dir: Path) -> tuple[float, Path]:
    """Parse the unique `timestep = ...` value from a `.inp` file in `base_dir`."""
    inp_files = list(Path(base_dir).glob("*.inp"))
    if not inp_files:
        raise FileNotFoundError(f"No .inp file found in {Path(base_dir)}.")
    if len(inp_files) > 1:
        raise RuntimeError(
            "Expected exactly one .inp file, found "
            f"{len(inp_files)}: {', '.join(path.name for path in inp_files)}"
        )

    inp_path = inp_files[0]
    with inp_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line_no_comment = line.split("#", 1)[0]
            match = TIMESTEP_PATTERN.search(line_no_comment)
            if match:
                return float(match.group(1)), inp_path

    raise ValueError(f"Could not find a 'timestep = ...' entry in {inp_path.name}.")


def format_timestamp(ntimestep: int, timestep_width_us: float) -> str:
    """Format a frame timestamp in seconds from a numeric frame stem."""
    time_s = ntimestep * timestep_width_us * 1e-6
    return f"t = {time_s:.2f} s"


def parse_chain_type_mapping_from_system_json(system_json_path: Path) -> Dict[str, str]:
    """Read `chain instance -> chain type` from an exported ionerdss system JSON."""
    with Path(system_json_path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    instances = payload.get("registries", {}).get("molecule_instances", [])
    chain_type_map: Dict[str, str] = {}
    for instance in instances:
        chain_type = instance.get("type")
        instance_name = instance.get("name", "")
        if not chain_type or not instance_name:
            continue

        suffix = f"_{chain_type}"
        if instance_name.endswith(suffix):
            chain_id = instance_name[: -len(suffix)]
        else:
            chain_id = instance_name
        chain_type_map[chain_id] = chain_type

    return chain_type_map


def parse_chain_type_mapping_from_summary(summary_path: Path) -> Dict[str, str]:
    """Read `chain_name_mapping` from a detailed summary text report."""
    pattern = re.compile(r"^\s*chain_name_mapping:\s*(\{.*\})\s*$")
    with Path(summary_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            match = pattern.match(line)
            if match:
                parsed = ast.literal_eval(match.group(1))
                if isinstance(parsed, dict):
                    return {str(key): str(value) for key, value in parsed.items()}
    return {}


def extract_chain_ids_from_pdb(pdb_path: Path) -> list[str]:
    """Extract chain identifiers from COM atoms in a NERDSS coarse-grained PDB."""
    chain_ids: list[str] = []
    seen: set[str] = set()

    with Path(pdb_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            parts = line.split()
            if len(parts) < 6 or parts[2] != "COM":
                continue

            chain_id = parts[4]
            if chain_id not in seen:
                seen.add(chain_id)
                chain_ids.append(chain_id)

    return chain_ids


def infer_chain_type_mapping_from_pdb(chain_ids: Iterable[str]) -> Dict[str, str]:
    """
    Conservative fallback when exported metadata is unavailable.

    Each chain is treated as its own type rather than guessing from naming.
    """
    return {chain_id: chain_id for chain_id in chain_ids}


def _candidate_metadata_files(base_dir: Path, pattern: str) -> list[Path]:
    direct = sorted(base_dir.glob(pattern))
    outputs = sorted((base_dir / "outputs").rglob(pattern)) if (base_dir / "outputs").exists() else []
    return [*direct, *outputs]


def resolve_chain_type_mapping(base_dir: Path, first_pdb_path: Path) -> Dict[str, str]:
    """
    Resolve `chain instance -> chain type` using exported metadata when available.

    Resolution order:
    1. `*_system.json`
    2. `*detailed_summary*.txt`
    3. first-frame PDB fallback
    """
    chain_ids = extract_chain_ids_from_pdb(first_pdb_path)

    for system_json in _candidate_metadata_files(Path(base_dir), "*_system.json"):
        mapping = parse_chain_type_mapping_from_system_json(system_json)
        if mapping and all(chain_id in mapping for chain_id in chain_ids):
            return {chain_id: mapping[chain_id] for chain_id in chain_ids}

    for summary_path in _candidate_metadata_files(Path(base_dir), "*detailed_summary*.txt"):
        mapping = parse_chain_type_mapping_from_summary(summary_path)
        if mapping and all(chain_id in mapping for chain_id in chain_ids):
            return {chain_id: mapping[chain_id] for chain_id in chain_ids}

    return infer_chain_type_mapping_from_pdb(chain_ids)


def build_chain_type_groups(chain_type_map: Dict[str, str]) -> Dict[str, list[str]]:
    """Group chain instances by chain type."""
    grouped: Dict[str, list[str]] = defaultdict(list)
    for chain_id, chain_type in chain_type_map.items():
        grouped[chain_type].append(chain_id)
    return {chain_type: sorted(chain_ids) for chain_type, chain_ids in sorted(grouped.items())}


def build_chain_type_color_map(chain_types: Iterable[str]) -> Dict[str, str]:
    """Assign a stable PyMOL color name to each chain type."""
    unique_types = sorted(set(chain_types))
    color_map: Dict[str, str] = {}
    for color_index, chain_type in enumerate(unique_types):
        color_name, _ = DEFAULT_TYPE_COLORS[color_index % len(DEFAULT_TYPE_COLORS)]
        color_map[chain_type] = color_name
    return color_map


def _estimate_content_anchor(
    image_path: Path,
    *,
    margin_x: int,
    margin_y: int,
    box_width: int,
    box_height: int,
    threshold: int = 245,
) -> tuple[int, int]:
    """Estimate a stable lower-right overlay anchor from a single frame."""
    image = Image.open(image_path).convert("L")
    width, height = image.size
    pixels = image.load()

    xs: list[int] = []
    ys: list[int] = []

    for y in range(height):
        for x in range(width):
            if pixels[x, y] < threshold:
                xs.append(x)
                ys.append(y)

    if xs and ys:
        content_right = max(xs)
        content_bottom = max(ys)
        return max(0, content_right - box_width), max(0, content_bottom - box_height)

    return max(0, width - margin_x - box_width), max(0, height - margin_y - box_height)


def add_timestamp_overlay_to_frame(
    image_path: Path | str,
    label_text: str,
    *,
    margin_x: int = 40,
    margin_y: int = 40,
    font_size: int = 42,
    pad_x: int = 18,
    pad_y: int = 10,
    anchor: Optional[tuple[int, int]] = None,
) -> tuple[int, int]:
    """
    Draw a timestamp box onto a rendered frame.

    If `anchor` is provided, that exact box origin is reused. Otherwise the box
    position is estimated from this frame and returned to the caller.
    """
    image_path = Path(image_path)
    image = Image.open(image_path).convert("RGBA")
    draw = ImageDraw.Draw(image)
    font = get_font(font_size)

    bbox = draw.textbbox((0, 0), label_text, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]

    box_width = text_width + 2 * pad_x
    box_height = text_height + 2 * pad_y

    if anchor is None:
        anchor = _estimate_content_anchor(
            image_path,
            margin_x=margin_x,
            margin_y=margin_y,
            box_width=box_width,
            box_height=box_height,
        )

    x0, y0 = anchor
    x1 = x0 + box_width
    y1 = y0 + box_height

    draw.rounded_rectangle(
        [x0, y0, x1, y1],
        radius=8,
        fill=(255, 255, 255, 235),
        outline=(0, 0, 0, 255),
        width=2,
    )

    text_x = x0 + pad_x
    text_y = y0 + pad_y
    for dx, dy in ((0, 0), (1, 0), (0, 1), (1, 1)):
        draw.text((text_x + dx, text_y + dy), label_text, font=font, fill=(0, 0, 0, 255))

    image.save(image_path)
    return anchor


def export_pymol_pdb_movie(
    *,
    base_dir: Path | str = ".",
    pdb_dir: Path | str = ".",
    add_timestamp_overlay: bool = True,
    timestamp_margin_x: int = 20,
    timestamp_margin_y: int = 20,
    timestamp_font_size: int = 96,
    timestamp_pad_x: int = 28,
    timestamp_pad_y: int = 28,
    viewport: tuple[int, int] = (2400, 1800),
    movie_fps: int = 5,
    out_dir_name: str = "frames",
) -> Path:
    """
    Export a one-frame-per-PDB movie from numeric NERDSS PDB snapshots.

    This function is intended to run inside a PyMOL Python session.
    """
    from pymol import cmd

    base_dir = Path(base_dir)
    pdb_dir = Path(pdb_dir)

    # 1. Load all PDBs
    print(f"Loading PDB files from: {pdb_dir.resolve()}")
    pdb_files = find_numeric_pdb_files(pdb_dir)
    timestep_width_us, _ = parse_timestep_width_us(base_dir)
    frame_chain_type_maps = {
        pdb_path.stem: resolve_chain_type_mapping(base_dir, pdb_path) for pdb_path in pdb_files
    }
    chain_type_color_map = build_chain_type_color_map(
        chain_type
        for chain_type_map in frame_chain_type_maps.values()
        for chain_type in chain_type_map.values()
    )

    obj_names: list[str] = []
    for pdb_path in pdb_files:
        obj_name = pdb_path.stem
        obj_names.append(obj_name)
        cmd.load(str(pdb_path), obj_name)

    cmd.hide("everything", "all")
    cmd.bg_color("white")
    cmd.set("antialias", 4)
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_frames", 1)
    cmd.set("orthoscopic", 1)
    cmd.viewport(*viewport)
    cmd.set("specular", 0.15)
    cmd.set("shininess", 10)
    cmd.set("hash_max", 300)

    # 2. Define Colors
    print("Defining Colors...")
    for color_index, chain_type in enumerate(sorted(chain_type_color_map)):
        color_name = chain_type_color_map[chain_type]
        _, rgb = DEFAULT_TYPE_COLORS[color_index % len(DEFAULT_TYPE_COLORS)]
        cmd.set_color(color_name, list(rgb))

    for obj_name in obj_names:
        chain_type_groups = build_chain_type_groups(frame_chain_type_maps[obj_name])
        cmd.show("spheres", f"{obj_name} and name COM")
        cmd.alter(f"{obj_name} and name COM", "vdw=1.5")
        cmd.rebuild(f"{obj_name} and name COM")
        cmd.set("sphere_scale", 1.0, f"{obj_name} and name COM")

        for chain_type, chain_ids in chain_type_groups.items():
            color_name = chain_type_color_map[chain_type]
            for chain_id in chain_ids:
                cmd.color(color_name, f"/{obj_name}///{chain_id}/COM")

    cmd.reset()
    cmd.zoom("all", buffer=1.0, complete=0)

    # 3. Create Animation
    print("Creating Animation...")
    n_frames = len(obj_names)
    cmd.disable("all")
    cmd.mset(f"1 x{n_frames}")
    for frame_index, obj_name in enumerate(obj_names, start=1):
        cmd.mdo(frame_index, f"disable all; enable {obj_name}; frame {frame_index}")
    cmd.set("movie_fps", movie_fps)

    cmd.disable("all")
    cmd.enable(obj_names[0])
    cmd.frame(1)

    out_dir = base_dir / out_dir_name
    out_dir.mkdir(exist_ok=True)
    cmd.mpng(str(out_dir / "frame"), 1, n_frames)

    if add_timestamp_overlay:
        first_label = format_timestamp(int(obj_names[0]), timestep_width_us)
        first_frame_path = out_dir / "frame0001.png"
        anchor = None
        if first_frame_path.exists():
            anchor = add_timestamp_overlay_to_frame(
                first_frame_path,
                first_label,
                margin_x=timestamp_margin_x,
                margin_y=timestamp_margin_y,
                font_size=timestamp_font_size,
                pad_x=timestamp_pad_x,
                pad_y=timestamp_pad_y,
            )

        for frame_index, obj_name in enumerate(obj_names[1:], start=2):
            frame_path = out_dir / f"frame{frame_index:04d}.png"
            if not frame_path.exists():
                continue
            add_timestamp_overlay_to_frame(
                frame_path,
                format_timestamp(int(obj_name), timestep_width_us),
                margin_x=timestamp_margin_x,
                margin_y=timestamp_margin_y,
                font_size=timestamp_font_size,
                pad_x=timestamp_pad_x,
                pad_y=timestamp_pad_y,
                anchor=anchor,
            )

    # Print ffmpeg commands to create a movie
    print(f"Exported frames to: {out_dir.resolve()}")
    print(f"To create a movie with ffmpeg:")
    print("Low-res MP4:")
    print("  ffmpeg -framerate 10 -i " + str(out_dir.resolve()) + "/frame%04d.png -vf scale=1280:-2 -c:v libx264 -crf 28 -pix_fmt yuv420p movie_lowres.mp4")
    print("\nHigh-res MP4:")
    print("  ffmpeg -framerate 10 -i " + str(out_dir.resolve()) + "/frame%04d.png -c:v libx264 -crf 16 -preset slow -pix_fmt yuv420p movie_highres.mp4")

    return out_dir
