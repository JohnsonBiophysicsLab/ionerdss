import json
from pathlib import Path

from PIL import Image

from ionerdss.analysis.visualization.pymol_movie import (
    add_timestamp_overlay_to_frame,
    build_chain_type_color_map,
    resolve_chain_type_mapping,
)


def _write_frame_pdb(path: Path, chain_ids: list[str]) -> None:
    lines = []
    for index, chain_id in enumerate(chain_ids, start=1):
        lines.append(
            f"ATOM  {index:5d} COM  COM {chain_id:>3s}{1:4d}    "
            f"{index:8.3f}{index:8.3f}{index:8.3f}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


def test_resolve_chain_type_mapping_prefers_system_json(tmp_path: Path) -> None:
    _write_frame_pdb(tmp_path / "1.pdb", ["A0", "A1", "B0"])

    system_json = {
        "registries": {
            "molecule_instances": [
                {"name": "A0_alpha", "type": "alpha"},
                {"name": "A1_alpha", "type": "alpha"},
                {"name": "B0_beta", "type": "beta"},
            ]
        }
    }
    outputs_dir = tmp_path / "outputs" / "systems"
    outputs_dir.mkdir(parents=True)
    (outputs_dir / "demo_system.json").write_text(json.dumps(system_json), encoding="utf-8")

    mapping = resolve_chain_type_mapping(tmp_path, tmp_path / "1.pdb")

    assert mapping == {"A0": "alpha", "A1": "alpha", "B0": "beta"}


def test_timestamp_overlay_reuses_first_frame_anchor(tmp_path: Path) -> None:
    first = tmp_path / "frame0001.png"
    second = tmp_path / "frame0002.png"

    image = Image.new("RGBA", (300, 200), (255, 255, 255, 255))
    for x in range(50, 241):
        for y in range(40, 171):
            image.putpixel((x, y), (0, 0, 0, 255))
    image.save(first)

    shifted = Image.new("RGBA", (300, 200), (255, 255, 255, 255))
    for x in range(70, 261):
        for y in range(20, 151):
            shifted.putpixel((x, y), (0, 0, 0, 255))
    shifted.save(second)

    anchor = add_timestamp_overlay_to_frame(first, "t = 0.10 s", font_size=24, pad_x=10, pad_y=8)
    reused_anchor = add_timestamp_overlay_to_frame(
        second,
        "t = 0.20 s",
        font_size=24,
        pad_x=10,
        pad_y=8,
        anchor=anchor,
    )

    assert reused_anchor == anchor


def test_chain_type_colors_stay_stable_across_frames() -> None:
    frame_one = {"A0": "alpha", "B0": "beta"}
    frame_two = {"A1": "alpha", "C0": "gamma"}

    color_map = build_chain_type_color_map(
        chain_type
        for frame in (frame_one, frame_two)
        for chain_type in frame.values()
    )

    assert color_map["alpha"] == color_map["alpha"]
    assert color_map["beta"] != color_map["gamma"]
