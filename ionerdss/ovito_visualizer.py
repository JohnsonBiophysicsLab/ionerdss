"""
ionerdss.ovito_visualizer

Module for visualizing NERDSS xyz trajectories using OVITO.
This module requires the optional `ovito_rendering` dependencies.
"""

import os
import sys
import shutil
import tempfile
import logging

logger = logging.getLogger(__name__)

# Qt bindings that ship their own copy of the Qt libraries. OVITO bundles its
# own PySide6/Qt6 build, and loading two different Qt6 copies into one process
# aborts the interpreter (this is what kills a Jupyter kernel outright, since a
# native abort cannot be caught by `except Exception`).
_CONFLICTING_QT_MODULES = ("PyQt6", "PyQt5", "PySide2")

# OVITO only stopped creating the global Qt application object eagerly at import
# time in 3.11.0; before that it clashes with any other Qt-based package.
_MIN_SAFE_OVITO_VERSION = (3, 11, 0)

# Headless `render_image()` via the OpenGL/Standard renderer only became
# reliable in 3.16.0. Below that we prefer a software renderer.
_HEADLESS_OPENGL_OK_VERSION = (3, 16, 0)

# Renderers are tried in this order. The software ray tracers need no OpenGL
# context at all, so they are attempted first: a renderer that cannot be
# constructed raises a catchable Python error, whereas an OpenGL renderer with
# no usable context can take the whole process down on some platforms.
_SOFTWARE_RENDERERS = ("TachyonRenderer", "OSPRayRenderer")
_GPU_RENDERERS = ("StandardRenderer", "OpenGLRenderer")

# Rendered frames are streamed to the GIF one at a time, but warn when a
# trajectory is long enough that the render itself will take a while.
_LONG_TRAJECTORY_WARNING = 500


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


def _check_ovito_version(ovito):
    """Warn about OVITO versions with known crash-in-Jupyter behaviour."""
    version = getattr(ovito, "version", None)
    if not isinstance(version, tuple):
        return None
    if version < _MIN_SAFE_OVITO_VERSION:
        logger.warning(
            "OVITO %s initializes its Qt application object at import time and is "
            "known to crash the interpreter when other Qt-based packages are present. "
            "Please upgrade: pip install -U ovito",
            ".".join(map(str, version)),
        )
    return version


def _select_renderer(vis, ovito_version):
    """Return (renderer, name), preferring backends that need no GL context."""
    candidates = list(_SOFTWARE_RENDERERS) + list(_GPU_RENDERERS)
    if ovito_version is not None and ovito_version >= _HEADLESS_OPENGL_OK_VERSION:
        # 3.16+ renders offscreen without a windowing system, and it is much
        # faster than the ray tracers.
        candidates = list(_GPU_RENDERERS) + list(_SOFTWARE_RENDERERS)

    for name in candidates:
        cls = getattr(vis, name, None)
        if cls is None:
            continue
        try:
            renderer = cls()
        except Exception as e:  # e.g. renderer only available in OVITO Pro
            logger.debug("Renderer %s unavailable: %s", name, e)
            continue
        logger.info("Using OVITO renderer: %s", name)
        return renderer, name

    logger.warning(
        "No renderer could be instantiated; falling back to OVITO's default. "
        "If the process dies here, the OpenGL renderer has no usable offscreen "
        "context -- upgrade OVITO to 3.16 or newer."
    )
    return None, "default"


def _render_frame(vp, output_path, size, frame, renderer):
    """Render one frame, retrying once without an explicit renderer."""
    if renderer is None:
        vp.render_image(filename=output_path, size=size, frame=frame)
        return None
    try:
        vp.render_image(filename=output_path, size=size, frame=frame,
                        renderer=renderer)
        return renderer
    except Exception as e:
        logger.warning(
            "Rendering with the selected renderer failed (%s); retrying with "
            "OVITO's default renderer.", e
        )
        vp.render_image(filename=output_path, size=size, frame=frame)
        return None


def _open_gif_writer(imageio, gif_name, fps, loop):
    """Open a streaming GIF writer, tolerating imageio's fps/duration churn."""
    try:
        return imageio.get_writer(gif_name, mode="I", fps=fps, loop=loop)
    except TypeError:
        return imageio.get_writer(gif_name, mode="I", duration=1.0 / fps, loop=loop)


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

    _check_qt_conflict()

    try:
        import warnings
        warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
        import imageio
        import ovito
        from ovito.io import import_file
        from ovito.vis import Viewport
        import ovito.vis as vis
    except ImportError as e:
        logger.error(
            "OVITO visualization requires the `ovito`, `imageio`, and `Pillow` packages.\n"
            "Please install the required dependencies using: pip install \"ioNERDSS[ovito_rendering]\""
        )
        raise e

    ovito_version = _check_ovito_version(ovito)

    logger.info(f"Loading trajectory from {trajectory_path} into OVITO...")
    pipeline = import_file(trajectory_path)
    pipeline.add_to_scene()

    # `save_gif=False` deliberately leaves the pipeline in the scene for
    # interactive use; every other exit path takes it back out so that repeated
    # calls do not render previously loaded trajectories on top of each other.
    keep_in_scene = not save_gif
    try:
        if not show_simulation_box:
            logger.info("Disabling simulation box visualization by default.")
            data = pipeline.compute()
            if data and data.cell:
                data.cell.vis.enabled = False

        vp = Viewport(type=Viewport.Type.PERSPECTIVE)
        vp.zoom_all()

        if not save_gif:
            logger.info("save_gif is False. Only loading into scene. Note that GUI rendering requires ovito Pro.")
            return

        num_frames = pipeline.source.num_frames
        if num_frames == 0:
            raise ValueError(
                f"OVITO found no frames in '{trajectory_path}'. The file may be "
                "empty or truncated."
            )
        if num_frames > _LONG_TRAJECTORY_WARNING:
            logger.warning(
                "Trajectory has %d frames; rendering will take a while. Reduce it by "
                "increasing `trajWrite` in parms.inp before rerunning NERDSS.",
                num_frames,
            )

        renderer, _ = _select_renderer(vis, ovito_version)

        logger.info("Rendering frames and compiling GIF at %s fps...", fps)
        temp_dir = tempfile.mkdtemp()
        try:
            # Frames are appended to the GIF as they are rendered, so peak memory
            # stays at one frame regardless of trajectory length.
            imread = getattr(getattr(imageio, "v2", imageio), "imread")
            with _open_gif_writer(imageio, gif_name, fps, loop) as writer:
                for frame in range(num_frames):
                    logger.info(f"Rendering frame {frame+1}/{num_frames}")
                    output_path = os.path.join(temp_dir, f"frame_{frame:04d}.png")
                    renderer = _render_frame(
                        vp, output_path, (800, 600), frame, renderer
                    )
                    writer.append_data(imread(output_path))
                    os.remove(output_path)
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

        logger.info(f"Successfully saved animation to {os.path.abspath(gif_name)}")
    finally:
        if not keep_in_scene:
            pipeline.remove_from_scene()
