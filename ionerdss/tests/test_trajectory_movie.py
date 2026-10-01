import logging
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from ionerdss.analysis.visualization import trajectory_movie as tm
from ionerdss.analysis.visualization.trajectory_movie import (
    MoleculeTypeInfo,
    TrajectoryRenderer,
    iter_xyz_frames,
    molecule_radius_nm,
    read_parms,
    read_pdb_frame,
    render_trajectory_movie,
    time_label_format,
)
from ionerdss.utils.diffusion_constant import compute_diffusion_constants_nm_us


# --------------------------------------------------------------------------
# A small NERDSS run written the way NERDSS writes it
# --------------------------------------------------------------------------

TYPES = {"A": (3.0, [(2.5, 0.0, 0.0)]), "Bcd": (4.0, [(-3.0, 1.0, 0.0), (0.0, 3.5, 0.5)])}


def _pdb_atom(serial, name, mol_type, mol_index, xyz, complex_id):
    # write_pdb.cpp: name %-4s, resname = type[:3] %3s, resSeq %4d, x y z %8.1f.
    return (f"ATOM  {serial % 100000:>5} {name:<4} {mol_type[:3]:>3} {mol_index % 10000:>4}    "
            f"{xyz[0]:>8.1f}{xyz[1]:>8.1f}{xyz[2]:>8.1f}{1.0:>6.1f}{complex_id:>6}{mol_type[:4]:>12}\n")


def _write_run(root: Path, frames, box=100.0, time_step="timeStep = 0.5", xyz_count_offset=0):
    """frames: list of (iteration, [(type, (x, y, z)), ...]) with box-centered coordinates."""
    (root / "PDB").mkdir(parents=True)
    (root / "DATA").mkdir()
    (root / "parms.inp").write_text(
        "start parameters\n"
        f"\tnItr = 100000\n\t{time_step}   # us\n\tpdbWrite = 1000\n\ttrajWrite = 1000\n"
        "end parameters\n\n"
        f"start boundaries\n\tWaterBox = [{box}, {box}, {box}]\nend boundaries\n\n"
        "start molecules\n\tBcd : 1\n\tA : 2\nend molecules\n"
    )
    for name, (radius, sites) in TYPES.items():
        D, Dr = compute_diffusion_constants_nm_us(radius)
        lines = [f"Name = {name}", "checkOverlap = true", "",
                 f"D = [{D:.6g}, {D:.6g}, {D:.6g}]", "", f"Dr = [{Dr:.6g}, {Dr:.6g}, {Dr:.6g}]", "",
                 "COM\t0.0000\t0.0000\t0.0000"]
        lines += [f"s{i}\t{x:.6f}\t{y:.6f}\t{z:.6f}" for i, (x, y, z) in enumerate(sites)]
        lines += ["", f"bonds = {len(sites)}"] + [f"com s{i}" for i in range(len(sites))]
        (root / f"{name}.mol").write_text("\n".join(lines) + "\n")

    with open(root / "DATA" / "trajectory.xyz", "w") as traj:
        for iteration, molecules in frames:
            n_units = sum(1 + len(TYPES[t][1]) for t, _ in molecules)
            traj.write(f"{n_units + xyz_count_offset}\niteration: {iteration}\n")
            pdb = [f"TITLE  PDB TIMESTEP {iteration} CREATED Mon Sep 28 12:00:00 2026\n",
                   f"CRYST1{box:<9g}{box:<9g}{box:<9g}{90:<7}{90:<7}{90:<7} P 1\n"]
            serial = 0
            for index, (mol_type, com) in enumerate(molecules):
                com = np.array(com, dtype=float)
                label = mol_type[:4].rjust(4)
                traj.write(f"{label} {com[0]:12.6f}{com[1]:12.6f}{com[2]:12.6f}\n")
                pdb.append(_pdb_atom(serial, "COM", mol_type, index, com + box / 2, index))
                serial += 1
                for k, site in enumerate(TYPES[mol_type][1]):
                    q = com + np.array(site)
                    traj.write(f"{label} {q[0]:12.6f}{q[1]:12.6f}{q[2]:12.6f}\n")
                    pdb.append(_pdb_atom(serial, f"s{k}", mol_type, index, q + box / 2, index))
                    serial += 1
            pdb.append("END\n")
            (root / "PDB" / f"{iteration}.pdb").write_text("".join(pdb))
    return root


FRAMES = [
    (0, [("Bcd", (-30.0, 10.0, 5.0)), ("A", (20.0, -15.0, 0.0)), ("A", (0.0, 30.0, -20.0))]),
    (1000, [("Bcd", (-25.0, 12.0, 5.0)), ("A", (-18.0, 12.0, 5.0)), ("A", (5.0, 28.0, -18.0))]),
    (2000, [("Bcd", (-20.0, 14.0, 6.0)), ("A", (-13.0, 14.0, 6.0)), ("A", (8.0, 25.0, -15.0))]),
]


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------

@pytest.mark.parametrize("line", ["timeStep = 0.5", "timestep = 0.5", "TIMESTEP=0.5 # us"])
def test_read_parms_matches_keys_case_insensitively_like_nerdss(tmp_path, line):
    _write_run(tmp_path, FRAMES[:1], time_step=line)
    box, time_step, names = read_parms(tmp_path / "parms.inp")
    assert time_step == 0.5
    np.testing.assert_allclose(box, [100.0, 100.0, 100.0])
    assert names == ["Bcd", "A"]


def test_radius_is_the_one_behind_D_and_Dr():
    # Values ionerdss wrote for 1ca4 (radius 1.592 nm), at the 6 digits of the .mol file.
    info = MoleculeTypeInfo("A", D=154.078, Dr=45.5648, site_distances=[0.9])
    radius, source = molecule_radius_nm(info)
    assert radius == pytest.approx(1.592, abs=1e-3)
    assert source == "D, Dr (Stokes-Einstein)"


def test_hand_set_D_and_Dr_fall_back_to_the_nerdss_radius(caplog):
    # D = 10, Dr = 0.1 imply 24.5 nm and 12.3 nm: not a Stokes-Einstein pair.
    info = MoleculeTypeInfo("C", D=10.0, Dr=0.1, site_distances=[1.2, 2.21])
    with caplog.at_level(logging.WARNING):
        radius, source = molecule_radius_nm(info)
    assert radius == pytest.approx(2.21)
    assert source == "COM-interface distance"
    assert "not a Stokes-Einstein pair" in caplog.text


def test_radius_from_D_alone_when_there_is_nothing_else():
    D, _ = compute_diffusion_constants_nm_us(2.5)
    radius, _ = molecule_radius_nm(MoleculeTypeInfo("X", D=D))
    assert radius == pytest.approx(2.5)


# --------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------

def test_pdb_reader_returns_box_centered_coms(tmp_path):
    _write_run(tmp_path, FRAMES)
    frame = read_pdb_frame(tmp_path / "PDB" / "1000.pdb", 1000, np.array([100.0] * 3))
    assert list(frame.labels) == ["Bcd", "A", "A"]
    np.testing.assert_allclose(frame.positions, [c for _, c in FRAMES[1][1]], atol=0.051)


def test_pdb_reader_falls_back_for_older_nerdss_layouts(tmp_path):
    # Layout of NERDSS PDB files from early 2025 (not the current fixed columns).
    path = tmp_path / "10000.pdb"
    path.write_text(
        "TITLE  PDB TIMESTEP                                    0 CREATED Tue Mar 18 01:43:18 2025\n"
        "CRYST1  200      200      200      90     90     90      P1   \n"
        "ATOM      0 COM    C     0       128.7    99.9    58.6     0     0CL\n"
        "ATOM      1 A1     C     0       126.7    99.7    57.6     0     0CL\n"
        "ATOM      3 COM    A     1       137.9    94.4    77.2     0     0CL\n"
    )
    assert tm._read_pdb_com_fixed_width(path.read_bytes()) is None
    frame = read_pdb_frame(path, 10000, np.array([200.0] * 3))
    assert list(frame.labels) == ["C", "A"]
    np.testing.assert_allclose(frame.positions, [[28.7, -0.1, -41.4], [37.9, -5.6, -22.8]])


def test_xyz_reader_skips_interface_lines_and_ignores_header_counts(tmp_path):
    # NERDSS's frame header can overstate the lines written; frames are split
    # on the `iteration:` lines instead.
    _write_run(tmp_path, FRAMES, xyz_count_offset=7)
    frames = list(iter_xyz_frames(tmp_path / "DATA" / "trajectory.xyz", {"A": 1, "Bcd": 2}))
    assert [f.iteration for f in frames] == [0, 1000, 2000]
    for frame, (_, molecules) in zip(frames, FRAMES):
        assert list(frame.labels) == [t for t, _ in molecules]
        np.testing.assert_allclose(frame.positions, [c for _, c in molecules], atol=1e-6)
    only = list(iter_xyz_frames(tmp_path / "DATA" / "trajectory.xyz", {"A": 1, "Bcd": 2}, wanted={1000}))
    assert [f.iteration for f in only] == [1000]


# --------------------------------------------------------------------------
# The three movie requirements
# --------------------------------------------------------------------------

def test_time_label_uses_one_unit_and_enough_decimals_for_the_frame_spacing():
    # 0.2 us steps written every 10000 iterations, plus NERDSS's final frame at nItr - 1.
    times = [i * 10000 * 0.2 for i in range(100)] + [999999 * 0.2]
    assert time_label_format(times)[:3] == ("ms", 1e3, 0)
    # Still true when frame_stride leaves only two spacings, one of them off-grid.
    assert time_label_format(times[::50])[:3] == ("ms", 1e3, 0)
    # An automatically chosen timeStep gives irregular spacings (here 0.2837 ms):
    # one decimal already tells neighbouring frames apart.
    auto = [i * 500 * 0.5674508470329174 for i in range(201)]
    assert time_label_format(auto)[:3] == ("ms", 1e3, 1)
    assert time_label_format([0.0, 2474.173903465271])[:3] == ("ms", 1e3, 1)
    # 1.4998 ms apart: one decimal, so the label steps 1.5 ms each frame
    # instead of alternating between 1 and 2 ms.
    slow = [i * 6729 * 0.22288708680117642 for i in range(301)]
    assert time_label_format(slow)[:3] == ("ms", 1e3, 1)
    assert time_label_format([0, 500, 1000, 1500])[:3] == ("ms", 1e3, 1)
    assert time_label_format([0.0, 10.0, 20.0])[:3] == ("µs", 1.0, 0)
    assert time_label_format([0, 5e5, 1e6, 1.5e6])[:3] == ("s", 1e6, 1)
    assert time_label_format([0.0, 10.0, 20.0], unit="ms")[:3] == ("ms", 1e3, 2)
    with pytest.raises(ValueError):
        time_label_format([0.0], unit="minutes")


def test_time_label_keeps_its_length_and_position():
    iterations = [0, 1000, 20000, 999999]
    renderer = TrajectoryRenderer([100.0] * 3, ["A"], {"A": 2.0}, time_step_us=0.5,
                                  iterations=iterations, show_box=False)
    texts = [renderer.time_text(it * 0.5, it) for it in iterations]
    assert len({len(t) for t in texts}) == 1
    assert len({renderer.time_font.getbbox(t)[2] for t in texts}) == 1
    assert texts[0].startswith("t = ") and texts[0].endswith(" ms")
    # Time only by default, in a label much larger than legend text.
    assert renderer.legend_entries == []
    assert renderer.time_font_size >= 2 * renderer.font_size

    # Without a time step, the label counts iterations in a fixed-width field too.
    no_dt = TrajectoryRenderer([100.0] * 3, ["A"], {"A": 2.0}, iterations=iterations,
                               show_box=False, show_legend=False)
    assert no_dt.time_text(None, 20000) == "step 20,000"


def _sphere_footprint(renderer, positions, margin=2.0):
    """Pixels within a sphere's projected radius (plus a margin) of any center."""
    centers = renderer._project(np.asarray(positions, dtype=float))
    yy, xx = np.mgrid[0:renderer.height, 0:renderer.width]
    near = np.zeros((renderer.height, renderer.width), dtype=bool)
    for cx, cy in centers:
        near |= (xx + 0.5 - cx) ** 2 + (yy + 0.5 - cy) ** 2 <= (renderer.radii_px.max() + margin) ** 2
    return near


def test_box_and_framing_do_not_move_with_the_molecules():
    renderer = TrajectoryRenderer([100.0] * 3, ["A"], {"A": 3.0}, size=(300, 300),
                                  show_time=False, show_legend=False)
    empty = np.asarray(renderer.render(np.zeros((0, 3)), []))
    for positions in ([[-45.0, -45.0, -45.0]], [[40.0, 40.0, 40.0], [0.0, 0.0, 0.0], [-10.0, 30.0, 0.0]]):
        frame = np.asarray(renderer.render(np.array(positions), ["A"] * len(positions)))
        changed = (frame != empty).any(axis=2)
        # Only pixels under the spheres change; the box is drawn identically.
        assert changed.any()
        assert not (changed & ~_sphere_footprint(renderer, positions)).any()


def test_each_type_has_one_color_free_or_in_a_complex():
    # Camera on the +x axis: screen right is +y, screen up is +z, all at equal depth.
    renderer = TrajectoryRenderer([100.0] * 3, ["A", "B"], {"A": 3.0, "B": 3.0}, size=(300, 300),
                                  view=(0.0, 0.0), show_time=False, show_legend=False, show_box=False,
                                  radius_scale=1.0)
    positions = np.array([[0.0, -30.0, 0.0],   # A, free
                          [0.0, 10.0, 0.0],    # A, bound to the B next to it
                          [0.0, 16.0, 0.0]])   # B
    image = np.asarray(renderer.render(positions, ["A", "A", "B"])).astype(int)
    centers = np.floor(renderer._project(positions)).astype(int)
    free_a, bound_a, b = (image[y, x] for x, y in centers)
    # Same type, same depth: only sub-pixel sampling of the center differs.
    assert np.abs(free_a - bound_a).max() <= 12
    assert np.abs(free_a - b).max() > 40
    # Each center is closest to its own type's color.
    type_colors = renderer.type_rgb[:2] * 255
    for pixel, expected in ((free_a, 0), (bound_a, 0), (b, 1)):
        hue = pixel / max(pixel.max(), 1)
        refs = type_colors / type_colors.max(axis=1, keepdims=True)
        assert np.argmin(np.abs(refs - hue).sum(axis=1)) == expected


@pytest.mark.parametrize("radius_scale", [1.0, 3.0])
def test_sphere_radius_is_the_molecule_radius_times_radius_scale(radius_scale):
    renderer = TrajectoryRenderer([100.0] * 3, ["A"], {"A": 5.0}, size=(400, 400), view=(0.0, 0.0),
                                  show_time=False, show_legend=False, show_box=False,
                                  radius_scale=radius_scale)
    image = np.asarray(renderer.render(np.zeros((1, 3)), ["A"]))
    covered = (image < 200).any(axis=2)
    width_px = covered.any(axis=0).sum()
    assert width_px == pytest.approx(2 * 5.0 * radius_scale * renderer.scale, abs=2)


def test_spheres_are_drawn_three_times_their_radius_by_default():
    renderer = TrajectoryRenderer([100.0] * 3, ["A"], {"A": 2.5}, size=(300, 300), show_legend=True)
    assert renderer.radius_scale == 3.0
    np.testing.assert_allclose(renderer.radii_px[0], 3 * 2.5 * renderer.scale)
    # The legend lists the molecular radius, not the drawn one.
    assert renderer.legend_entries[0][3] == "A (r = 2.5 nm)"
    with pytest.raises(ValueError):
        TrajectoryRenderer([100.0] * 3, ["A"], {"A": 2.5}, radius_scale=0)


def test_unknown_types_are_drawn_grey_with_one_warning(caplog):
    renderer = TrajectoryRenderer([100.0] * 3, ["A"], {"A": 3.0}, size=(200, 200),
                                  show_time=False, show_legend=False, show_box=False)
    with caplog.at_level(logging.WARNING):
        renderer.render(np.zeros((2, 3)), ["Z", "Z"])
        renderer.render(np.zeros((1, 3)), ["Z"])
    assert caplog.text.count("'Z'") == 1


# --------------------------------------------------------------------------
# End to end
# --------------------------------------------------------------------------

@pytest.mark.parametrize("source", ["pdb", "xyz"])
def test_render_trajectory_movie_writes_a_gif(tmp_path, source):
    run = _write_run(tmp_path / "run", FRAMES)
    out = render_trajectory_movie(run, tmp_path / "movie.gif", source=source, size=(240, 200),
                                  progress=False)
    with Image.open(out) as gif:
        assert gif.n_frames == len(FRAMES)
        assert gif.size == (240, 200)


def test_frames_only_output_and_frame_stride(tmp_path):
    run = _write_run(tmp_path / "run", FRAMES)
    frames_dir = tmp_path / "frames"
    result = render_trajectory_movie(run, None, frames_dir=frames_dir, frame_stride=2,
                                     size=(160, 160), progress=False)
    assert result is None
    assert sorted(p.name for p in frames_dir.iterdir()) == ["frame_00000.png", "frame_00001.png"]


def test_mp4_without_ffmpeg_says_how_to_get_it(tmp_path, monkeypatch):
    run = _write_run(tmp_path / "run", FRAMES)
    monkeypatch.setattr(tm, "_find_ffmpeg", lambda: None)
    with pytest.raises(RuntimeError, match="ffmpeg"):
        render_trajectory_movie(run, tmp_path / "movie.mp4", size=(160, 160), progress=False)


def test_mp4_frames_are_streamed_to_ffmpeg(tmp_path, monkeypatch):
    # A stand-in ffmpeg that records the raw frames it is sent and its arguments.
    fake = tmp_path / "ffmpeg"
    fake.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "data = sys.stdin.buffer.read()\n"
        "open(sys.argv[-1], 'wb').write(data)\n"
        "open(sys.argv[-1] + '.args', 'w').write(' '.join(sys.argv[1:]))\n"
    )
    fake.chmod(0o755)
    monkeypatch.setattr(tm, "_find_ffmpeg", lambda: str(fake))
    run = _write_run(tmp_path / "run", FRAMES)
    out = render_trajectory_movie(run, tmp_path / "movie.mp4", size=(161, 120), progress=False)
    # Sizes are rounded up to even numbers for yuv420p.
    assert out.stat().st_size == len(FRAMES) * 162 * 120 * 3
    assert "-s 162x120" in (tmp_path / "movie.mp4.args").read_text()


@pytest.mark.skipif(tm._find_ffmpeg() is None, reason="ffmpeg not installed")
def test_render_trajectory_movie_writes_an_mp4(tmp_path):
    run = _write_run(tmp_path / "run", FRAMES)
    out = render_trajectory_movie(run, tmp_path / "movie.mp4", size=(161, 161), progress=False)
    assert out.stat().st_size > 0


def _break_process_pools(monkeypatch, accepted=None):
    """Replace process pools with one whose workers die on start-up.

    The first `accepted` frames (all by default) are queued and then fail;
    later ones are refused, as a real pool refuses them once it knows it is
    broken.
    """
    import concurrent.futures
    from concurrent.futures.process import BrokenProcessPool

    class BrokenPool:
        def __init__(self, *args, **kwargs):
            self.queued = 0

        def submit(self, fn, *args, **kwargs):
            if self.queued == accepted:
                raise BrokenProcessPool("the pool is not usable anymore")
            self.queued += 1
            future = concurrent.futures.Future()
            future.set_exception(BrokenProcessPool("workers died"))
            return future

        def shutdown(self, *args, **kwargs):
            pass

    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", BrokenPool)


@pytest.mark.parametrize("accepted", [None, 1], ids=["queued-frames-fail", "later-frames-refused"])
def test_all_frames_render_when_worker_processes_cannot_start(tmp_path, monkeypatch, caplog, accepted):
    _break_process_pools(monkeypatch, accepted)
    run = _write_run(tmp_path / "run", FRAMES)
    with caplog.at_level(logging.WARNING):
        render_trajectory_movie(run, None, frames_dir=tmp_path / "frames", size=(160, 160),
                                n_jobs=2, progress=False)
    assert "rendering in this process instead" in caplog.text
    # Every frame is rendered, including the one whose failure was seen first.
    assert sorted(p.name for p in (tmp_path / "frames").iterdir()) == [
        "frame_00000.png", "frame_00001.png", "frame_00002.png"]


def test_serial_fallback_reads_frames_only_as_it_renders_them(monkeypatch):
    _break_process_pools(monkeypatch)
    kwargs = dict(box_nm=[100.0] * 3, type_names=["A"], radii_nm={"A": 3.0}, size=(80, 80),
                  show_time=False)
    read = []

    def tasks():
        for iteration in range(10):
            read.append(iteration)
            yield ("arrays", iteration, (np.zeros((1, 3)), ["A"]))

    frames = tm._render_frames(TrajectoryRenderer(**kwargs), kwargs, tasks(), kwargs["box_nm"],
                               None, True, n_jobs=2)
    next(frames)
    # The first frame did not wait for the whole trajectory to be read.
    assert len(read) < 10
    assert sum(1 for _ in frames) == 9


def test_parallel_rendering_matches_serial(tmp_path):
    run = _write_run(tmp_path / "run", FRAMES)
    for jobs in (1, 2):
        render_trajectory_movie(run, None, frames_dir=tmp_path / f"jobs{jobs}", size=(160, 160),
                                n_jobs=jobs, progress=False)
    for name in ("frame_00000.png", "frame_00001.png", "frame_00002.png"):
        serial = np.asarray(Image.open(tmp_path / "jobs1" / name))
        parallel = np.asarray(Image.open(tmp_path / "jobs2" / name))
        np.testing.assert_array_equal(serial, parallel)


def test_truncated_trajectory_labels_map_back_to_full_type_names(caplog):
    # NERDSS writes type[:3] in PDB files and type[:4] in trajectory.xyz.
    assert tm._match_labels(["Clathrin", "AP2"], 3) == {"Cla": "Clathrin", "AP2": "AP2"}
    with caplog.at_level(logging.WARNING):
        labels = tm._match_labels(["GagPol", "GagX"], 3)
    assert labels == {"Gag": "GagPol"}
    assert "share a label" in caplog.text
    renderer = TrajectoryRenderer([100.0] * 3, ["Clathrin", "AP2"], {"Clathrin": 3.0, "AP2": 2.0},
                                  labels={"Cla": "Clathrin"}, size=(200, 200), show_legend=True)
    np.testing.assert_array_equal(renderer.type_indices(["Cla", "AP2", "Clathrin"]), [0, 1, 0])
    assert [text.split(" (")[0] for _, _, _, text in renderer.legend_entries] == ["Clathrin", "AP2"]


# --------------------------------------------------------------------------
# Deprecated pipelines
# --------------------------------------------------------------------------

def test_visualize_trajectory_ovito_is_deprecated(tmp_path):
    from ionerdss.ovito_visualizer import visualize_trajectory_ovito

    with pytest.warns(DeprecationWarning, match="render_trajectory_movie"):
        with pytest.raises(FileNotFoundError):
            visualize_trajectory_ovito(str(tmp_path / "missing.xyz"))


def test_export_pymol_pdb_movie_is_deprecated(tmp_path):
    from ionerdss.analysis.visualization import export_pymol_pdb_movie

    with pytest.warns(DeprecationWarning, match="render_trajectory_movie"):
        # Without PyMOL the import fails; with it, there are no PDB files here.
        with pytest.raises((ImportError, FileNotFoundError)):
            export_pymol_pdb_movie(base_dir=tmp_path, pdb_dir=tmp_path)
