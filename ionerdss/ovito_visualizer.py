"""
ionerdss.ovito_visualizer

Module for visualizing NERDSS xyz trajectories using OVITO.
This module requires the optional `ovito_rendering` dependencies.

Rendering runs in a child process (see `_ovito_render_worker.py`). OVITO aborts
the interpreter rather than raising when it cannot render -- typically on a
headless machine, where the OpenGL path has no offscreen context -- and a native
abort cannot be caught by `except Exception`. Run in-process, that silently
kills a Jupyter kernel; isolated, it becomes a Python error naming the renderer
that failed, and the next renderer is tried automatically.
"""

import importlib.util
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile

logger = logging.getLogger(__name__)

_WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "_ovito_render_worker.py")

# Qt bindings that ship their own copy of the Qt libraries. OVITO bundles its
# own PySide6/Qt6 build, and loading two different Qt6 copies into one process
# aborts the interpreter. Only relevant when this module imports ovito itself,
# i.e. on the `save_gif=False` path.
_CONFLICTING_QT_MODULES = ("PyQt6", "PyQt5", "PySide2")

# Exit codes reported by the worker.
_EXIT_NO_FRAMES = 3
_EXIT_NO_RENDERER = 4

# Rendered frames land on disk one PNG per frame, so warn when a trajectory is
# long enough for that to be slow and bulky.
_LONG_TRAJECTORY_WARNING = 500

_RENDER_SIZE = (800, 600)


def _check_qt_conflict():
    """Fail loudly *before* importing ovito if another Qt build is loaded."""
    loaded = [m for m in _CONFLICTING_QT_MODULES if m in sys.modules]
    if loaded:
        raise RuntimeError(
            f"{', '.join(loaded)} is already imported in this process. OVITO ships its "
            "own PySide6/Qt6 build, and loading two Qt copies into one interpreter "
            "crashes it (in Jupyter this shows up as a dead kernel).\n"
            "Restart the kernel and run this function before importing anything that "
            "pulls in PyQt, or install ioNERDSS's ovito extra into an environment "
            "without PyQt."
        )


def _require(*packages):
    """Check optional dependencies without importing them."""
    missing = [p for p in packages if importlib.util.find_spec(p) is None]
    if missing:
        raise ImportError(
            f"OVITO visualization requires the missing package(s): {', '.join(missing)}.\n"
            "Please install the required dependencies using: "
            "pip install \"ioNERDSS[ovito_rendering]\""
        )


def _run_worker(trajectory_path, out_dir, show_simulation_box, exclude):
    """Render every frame in a child process.

    Returns (returncode, renderer_name_or_None, num_frames_or_None, stderr).
    A negative returncode means the child was killed by a signal, i.e. OVITO
    took the process down -- exactly what this isolation exists to contain.
    """
    config_path = os.path.join(out_dir, "render_config.json")
    with open(config_path, "w") as f:
        json.dump({
            "trajectory_path": os.path.abspath(trajectory_path),
            "out_dir": out_dir,
            "size": list(_RENDER_SIZE),
            "show_simulation_box": bool(show_simulation_box),
            "exclude": sorted(exclude),
        }, f)

    proc = subprocess.Popen(
        [sys.executable, _WORKER, config_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )

    renderer_name = None
    num_frames = None
    for line in proc.stdout:
        line = line.strip()
        if line.startswith("VERSION "):
            logger.info("OVITO version: %s", line.split(None, 1)[1])
        elif line.startswith("FRAMES "):
            num_frames = int(line.split()[1])
        elif line.startswith("RENDERER "):
            renderer_name = line.split(None, 1)[1]
            logger.info("Using OVITO renderer: %s", renderer_name)
        elif line.startswith("SKIP "):
            logger.debug("Renderer unavailable: %s", line.split(None, 1)[1])
        elif line.startswith("FRAME "):
            logger.info("Rendering frame %d/%s", int(line.split()[1]) + 1, num_frames)
        elif line:
            logger.debug("worker: %s", line)

    stderr = proc.stderr.read()
    proc.stdout.close()
    proc.stderr.close()
    return proc.wait(), renderer_name, num_frames, stderr


def _render_frames(trajectory_path, out_dir, show_simulation_box):
    """Render all frames, retrying with the next renderer if one takes the child down."""
    excluded = set()
    last_error = ""
    while True:
        attempt_dir = os.path.join(out_dir, f"attempt_{len(excluded)}")
        os.makedirs(attempt_dir)
        code, renderer, num_frames, stderr = _run_worker(
            trajectory_path, attempt_dir, show_simulation_box, excluded
        )

        if code == 0:
            if num_frames and num_frames > _LONG_TRAJECTORY_WARNING:
                logger.warning(
                    "Rendered %d frames; increase `trajWrite` in parms.inp to keep "
                    "future trajectories shorter.", num_frames
                )
            return attempt_dir, num_frames

        shutil.rmtree(attempt_dir, ignore_errors=True)
        last_error = stderr.strip()

        if code == _EXIT_NO_FRAMES:
            raise ValueError(
                f"OVITO found no frames in '{trajectory_path}'. The file may be "
                "empty or truncated."
            )
        if code == _EXIT_NO_RENDERER or renderer is None:
            # Either nothing is left to try, or the child never got as far as
            # picking a renderer -- retrying cannot help.
            raise RuntimeError(_failure_message(trajectory_path, excluded, last_error))

        how = f"was killed by signal {-code}" if code < 0 else f"exited with code {code}"
        logger.warning(
            "The %s renderer %s; retrying with the next available renderer.",
            renderer, how,
        )
        excluded.add(renderer)


def _failure_message(trajectory_path, excluded, stderr):
    if excluded:
        message = (
            f"OVITO could not render '{trajectory_path}' with any available renderer "
            f"(tried: {', '.join(sorted(excluded))}).\n"
            "This usually means OVITO cannot render offscreen on this machine. Headless "
            "rendering works out of the box from OVITO 3.16 onwards, so try: "
            "pip install -U ovito\n"
            "(If pip refuses to upgrade, a numpy<2 pin is holding it back -- current "
            "OVITO requires numpy>=2.)"
        )
    else:
        # The child never got as far as picking a renderer, so this is a problem
        # with the trajectory or the OVITO installation, not with rendering.
        message = (
            f"The OVITO render process failed before rendering started, while loading "
            f"'{trajectory_path}'."
        )
    if stderr:
        message += f"\n\nOutput from the render process:\n{stderr}"
    return message


def _open_gif_writer(imageio, gif_name, fps, loop):
    """Open a streaming GIF writer, tolerating imageio's fps/duration churn."""
    try:
        return imageio.get_writer(gif_name, mode="I", fps=fps, loop=loop)
    except TypeError:
        return imageio.get_writer(gif_name, mode="I", duration=1.0 / fps, loop=loop)


def _write_gif(frame_dir, gif_name, fps, loop):
    import imageio

    # Frames are appended as they are read back, so peak memory stays at one
    # frame regardless of trajectory length.
    imread = getattr(getattr(imageio, "v2", imageio), "imread")
    frames = sorted(f for f in os.listdir(frame_dir) if f.endswith(".png"))
    with _open_gif_writer(imageio, gif_name, fps, loop) as writer:
        for filename in frames:
            path = os.path.join(frame_dir, filename)
            writer.append_data(imread(path))
            os.remove(path)


def _load_into_scene(trajectory_path, show_simulation_box):
    """Load the trajectory into OVITO's scene in this process (no rendering)."""
    _check_qt_conflict()
    import warnings
    warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
    from ovito.io import import_file
    from ovito.vis import Viewport

    pipeline = import_file(trajectory_path)
    pipeline.add_to_scene()
    if not show_simulation_box:
        logger.info("Disabling simulation box visualization by default.")
        data = pipeline.compute()
        if data and data.cell:
            data.cell.vis.enabled = False
    Viewport(type=Viewport.Type.PERSPECTIVE).zoom_all()
    return pipeline


def visualize_trajectory_ovito(
    trajectory_path: str,
    save_gif: bool = True,
    gif_name: str = "trajectory.gif",
    fps: int = 10,
    loop: int = 0,
    show_simulation_box: bool = False
):
    """
    Visualizes a trajectory from an XYZ file and optionally saves it as a GIF.
    This functionality requires the `ovito_rendering` optional extra.

    Parameters:
        trajectory_path (str): Path to the XYZ trajectory file.
        save_gif (bool): If True, saves the trajectory animation as a GIF.
        gif_name (str): Name of the output GIF file (if save_gif is True).
        fps (int): Frames per second for the GIF animation.
        show_simulation_box (bool): If True, shows the simulation cell bounding box.
    """
    if not os.path.exists(trajectory_path):
        raise FileNotFoundError(f"Trajectory file '{trajectory_path}' not found.")

    _require("ovito", "imageio", "PIL")

    if not save_gif:
        logger.info("save_gif is False. Only loading into scene. Note that GUI rendering requires ovito Pro.")
        _load_into_scene(trajectory_path, show_simulation_box)
        return

    logger.info("Rendering %s with OVITO in a separate process...", trajectory_path)
    work_dir = tempfile.mkdtemp(prefix="ionerdss_ovito_")
    try:
        frame_dir, _ = _render_frames(trajectory_path, work_dir, show_simulation_box)
        logger.info("Compiling GIF at %s fps...", fps)
        _write_gif(frame_dir, gif_name, fps, loop)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    logger.info(f"Successfully saved animation to {os.path.abspath(gif_name)}")
