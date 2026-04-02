"""
ionerdss.ovito_visualizer

Module for visualizing NERDSS xyz trajectories using OVITO.
This module requires the optional `ovito_rendering` dependencies.
"""

import os
import tempfile
import logging

logger = logging.getLogger(__name__)

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

    try:
        import warnings
        warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
        import imageio
        from ovito.io import import_file
        from ovito.vis import Viewport
    except ImportError as e:
        logger.error(
            "OVITO visualization requires the `ovito`, `imageio`, and `Pillow` packages.\n"
            "Please install the required dependencies using: pip install \"ioNERDSS[ovito_rendering]\""
        )
        raise e

    logger.info(f"Loading trajectory from {trajectory_path} into OVITO...")
    pipeline = import_file(trajectory_path)
    pipeline.add_to_scene()
    
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

    logger.info("Rendering frames...")
    temp_dir = tempfile.mkdtemp()
    frame_paths = []

    for frame in range(pipeline.source.num_frames):
        logger.info(f"Rendering frame {frame+1}/{pipeline.source.num_frames}")
        output_path = os.path.join(temp_dir, f"frame_{frame:04d}.png")
        vp.render_image(filename=output_path, size=(800, 600), frame=frame)
        frame_paths.append(output_path)

    logger.info(f"Compiling GIF at {fps} fps...")
    images = []
    for filename in frame_paths:
        images.append(imageio.imread(filename))
    
    imageio.mimsave(gif_name, images, fps=fps, loop=loop)
    logger.info(f"Successfully saved animation to {os.path.abspath(gif_name)}")

    # Clean up temporary frames
    for filename in frame_paths:
        try:
            os.remove(filename)
        except OSError:
            pass
    try:
        os.rmdir(temp_dir)
    except OSError:
        pass
