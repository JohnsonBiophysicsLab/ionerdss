"""
Out-of-process OVITO frame renderer.

Run as a standalone script by ionerdss.ovito_visualizer:

    python _ovito_render_worker.py <config.json>

This is deliberately NOT imported by the parent process. OVITO can abort the
interpreter outright (no usable offscreen context, Qt library conflicts), and a
native abort cannot be caught with `except`. Inside a Jupyter kernel that means
a dead kernel with no traceback, so all rendering happens in a child process
whose death the parent can survive and report.

Protocol: progress is written to stdout as single-token lines the parent parses
(`VERSION`, `FRAMES`, `RENDERER`, `SKIP`, `FRAME`); everything else, including
OVITO's own fatal messages, goes to stderr and is surfaced by the parent.

Exit codes: 0 success, 1 unhandled exception, 3 empty trajectory,
4 every renderer excluded or unavailable.
"""

import json
import os
import sys

# Software ray tracers need no OpenGL context, so they are tried first on the
# OVITO versions whose OpenGL path cannot render offscreen.
SOFTWARE_RENDERERS = ("OSPRayRenderer", "TachyonRenderer")
GPU_RENDERERS = ("StandardRenderer", "OpenGLRenderer")

# 3.16 renders offscreen without a windowing system and is much faster than the
# ray tracers, so it goes first there.
HEADLESS_OPENGL_OK_VERSION = (3, 16, 0)

# Last resort: let OVITO choose. Kept as an explicit candidate so the parent can
# exclude it after it fails like any other.
DEFAULT_RENDERER = "default"

EXIT_NO_FRAMES = 3
EXIT_NO_RENDERER = 4


def candidate_renderers(ovito_version):
    if ovito_version is not None and ovito_version >= HEADLESS_OPENGL_OK_VERSION:
        ordered = list(GPU_RENDERERS) + list(SOFTWARE_RENDERERS)
    else:
        ordered = list(SOFTWARE_RENDERERS) + list(GPU_RENDERERS)
    return ordered + [DEFAULT_RENDERER]


def select_renderer(vis, ovito_version, exclude):
    """Return (renderer_object_or_None, name), skipping excluded/unavailable ones."""
    for name in candidate_renderers(ovito_version):
        if name in exclude:
            continue
        if name == DEFAULT_RENDERER:
            return None, DEFAULT_RENDERER
        cls = getattr(vis, name, None)
        if cls is None:
            print(f"SKIP {name} not available in this OVITO build", flush=True)
            continue
        try:
            return cls(), name
        except Exception as e:  # e.g. renderer requires OVITO Pro
            print(f"SKIP {name} {e}", flush=True)
    return None, None


def main(argv):
    with open(argv[1]) as f:
        cfg = json.load(f)

    import ovito
    import ovito.vis as vis
    from ovito.io import import_file
    from ovito.vis import Viewport

    ovito_version = ovito.version if isinstance(ovito.version, tuple) else None
    if ovito_version:
        print("VERSION " + ".".join(map(str, ovito_version)), flush=True)

    pipeline = import_file(cfg["trajectory_path"])
    pipeline.add_to_scene()

    if not cfg["show_simulation_box"]:
        data = pipeline.compute()
        if data and data.cell:
            data.cell.vis.enabled = False

    vp = Viewport(type=Viewport.Type.PERSPECTIVE)
    vp.zoom_all()

    num_frames = pipeline.source.num_frames
    print(f"FRAMES {num_frames}", flush=True)
    if num_frames == 0:
        return EXIT_NO_FRAMES

    renderer, name = select_renderer(vis, ovito_version, set(cfg["exclude"]))
    if name is None:
        return EXIT_NO_RENDERER
    # Printed before the first render so the parent knows what to blame, and to
    # exclude on retry, if this process dies mid-frame.
    print(f"RENDERER {name}", flush=True)

    size = tuple(cfg["size"])
    for frame in range(num_frames):
        print(f"FRAME {frame}", flush=True)
        output_path = os.path.join(cfg["out_dir"], f"frame_{frame:05d}.png")
        if renderer is None:
            vp.render_image(filename=output_path, size=size, frame=frame)
        else:
            vp.render_image(filename=output_path, size=size, frame=frame,
                            renderer=renderer)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
