"""
Legacy compatibility interface for ionerdss analysis.
Bridges between new modular API and existing plot_figure interface.
"""

import os
import warnings
from typing import Dict, Any, List, Optional, Tuple


class LegacyPlotInterface:
    """
    Legacy compatibility bridge for existing plot_figure interface.
    
    Translates old-style plot_figure calls to new modular API while
    maintaining complete backward compatibility.
    """
    
    def __init__(self, analysis):
        """
        Initialize legacy interface with reference to main Analysis object.
        
        Parameters:
            analysis: Main Analysis object with new API
        """
        self.analysis = analysis
        
    def plot_figure(self, figure_type: str, **kwargs):
        """
        Legacy plot_figure method with full backward compatibility.
        
        This method maintains the exact same interface as the original
        plot_figure method while internally using the new modular API.
        
        Parameters:
            figure_type (str): Type of figure to plot
            **kwargs: All original plot_figure parameters
            
        Returns:
            Result from the plotting operation
        """
        # Extract common parameters
        simulations = kwargs.get('simulations', None)
        x = kwargs.get('x', 'time')
        y = kwargs.get('y', 'species')
        z = kwargs.get('z', None)
        legend = kwargs.get('legend', None)
        
        # Extract plotting parameters
        bins = kwargs.get('bins', 10)
        time_bins = kwargs.get('time_bins', 10)
        time_frame = kwargs.get('time_frame', None)
        frequency = kwargs.get('frequency', False)
        normalize = kwargs.get('normalize', False)
        show_type = kwargs.get('show_type', 'both')
        font_size = kwargs.get('font_size', 12)
        figure_size = kwargs.get('figure_size', (10, 6))
        seaborn_style = kwargs.get('seaborn_style', 'ticks')
        seaborn_context = kwargs.get('seaborn_context', 'paper')
        user_file_name = kwargs.get('user_file_name', None)
        
        # Validate inputs
        if not legend:
            raise ValueError("Legend must be provided.")
        
        # Get configured data and plot objects
        data = self.analysis.load_simulations(
            simulations=simulations,
            time_frame=time_frame
        )
        
        plot_config = self.analysis.set_plot(
            style='seaborn' if seaborn_style else 'matplotlib',
            context=seaborn_context,
            font_size=font_size,
            figure_size=figure_size
        )
        
        # Map legacy plot configurations to new methods
        plot_mapping = {
            # Line plots
            ("line", "time", "species", None): self._plot_line_species_time,
            ("line", "time", "count", None): self._plot_line_count_time,
            ("line", "time", "maximum_assembly", None): self._plot_line_max_assembly_time,
            ("line", "time", "average_assembly", None): self._plot_line_avg_assembly_time,
            ("line", "time", "fraction_of_monomers_assembled", None): self._plot_line_fraction_time,
            ("line", "size", "free_energy", None): self._plot_line_free_energy_size,
            ("line", "size", "symmetric_association_probability", None): self._plot_line_sym_assoc_prob,
            ("line", "size", "asymmetric_association_probability", None): self._plot_line_asym_assoc_prob,
            ("line", "size", "symmetric_dissociation_probability", None): self._plot_line_sym_dissoc_prob,
            ("line", "size", "asymmetric_dissociation_probability", None): self._plot_line_asym_dissoc_prob,
            ("line", "size", "growth_probability", None): self._plot_line_growth_prob,
            ("line", "size", "lifetime", None): self._plot_line_lifetime,
            
            # Histogram plots
            ("hist", "size", "complex_count", None): self._plot_hist_complex_size,
            ("hist", "size", "monomer_count", None): self._plot_hist_monomer_count,
            
            # 3D histogram plots
            ("3dhist", "size", "time", "complex_count"): self._plot_3dhist_complex_size,
            ("3dhist", "size", "time", "monomer_count"): self._plot_3dhist_monomer_count,
            
            # Heatmap plots
            ("heatmap", "size", "time", "complex_count"): self._plot_heatmap_complex_size_time,
            ("heatmap", "size", "time", "monomer_count"): self._plot_heatmap_monomer_count_time,
            ("heatmap", "size", "size", "complex_count"): self._plot_heatmap_species_correlation,
            
            # Stacked histogram plots
            ("stacked", "size", "complex_count", None): self._plot_stacked_complex_size,
        }
        
        plot_config_key = (figure_type, x, y, z)
        
        if plot_config_key not in plot_mapping:
            raise ValueError(f"Unsupported plot configuration: {plot_config_key}")
        
        # Call the appropriate plotting method
        plot_method = plot_mapping[plot_config_key]
        return plot_method(data, plot_config, legend, kwargs)
    
    # Line plot implementations
    def _plot_line_species_time(self, data, plot_config, legend, kwargs):
        """Plot species copy number vs time."""
        return plot_config.line_speciescopy_vs_time(
            data=data, 
            legend=legend,
            user_file_name=kwargs.get('user_file_name'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_count_time(self, data, plot_config, legend, kwargs):
        """Plot complex count vs time."""
        return plot_config.line_complex_count_vs_time(
            data=data,
            target_complexes=legend,
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_max_assembly_time(self, data, plot_config, legend, kwargs):
        """Plot maximum assembly size vs time."""
        return plot_config.line_assembly_size_vs_time(
            data=data,
            assembly_type="maximum",
            legend=legend,
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_avg_assembly_time(self, data, plot_config, legend, kwargs):
        """Plot average assembly size vs time."""
        return plot_config.line_assembly_size_vs_time(
            data=data,
            assembly_type="average",
            legend=legend,
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_fraction_time(self, data, plot_config, legend, kwargs):
        """Plot fraction of monomers assembled vs time."""
        return plot_config.line_assembly_size_vs_time(
            data=data,
            assembly_type="fraction",
            legend=legend,
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_free_energy_size(self, data, plot_config, legend, kwargs):
        """Plot free energy vs cluster size."""
        return plot_config.line_free_energy_vs_size(
            data=data,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_sym_assoc_prob(self, data, plot_config, legend, kwargs):
        """Plot symmetric association probability vs size."""
        return plot_config.line_probability_vs_size(
            data=data,
            probability_type="symmetric_association",
            legend=legend,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_asym_assoc_prob(self, data, plot_config, legend, kwargs):
        """Plot asymmetric association probability vs size."""
        return plot_config.line_probability_vs_size(
            data=data,
            probability_type="asymmetric_association",
            legend=legend,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_sym_dissoc_prob(self, data, plot_config, legend, kwargs):
        """Plot symmetric dissociation probability vs size."""
        return plot_config.line_probability_vs_size(
            data=data,
            probability_type="symmetric_dissociation",
            legend=legend,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_asym_dissoc_prob(self, data, plot_config, legend, kwargs):
        """Plot asymmetric dissociation probability vs size."""
        return plot_config.line_probability_vs_size(
            data=data,
            probability_type="asymmetric_dissociation",
            legend=legend,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_growth_prob(self, data, plot_config, legend, kwargs):
        """Plot growth probability vs size."""
        return plot_config.line_probability_vs_size(
            data=data,
            probability_type="growth",
            legend=legend,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_line_lifetime(self, data, plot_config, legend, kwargs):
        """Plot lifetime vs size."""
        return plot_config.line_probability_vs_size(
            data=data,
            probability_type="lifetime",
            legend=legend,
            time_frame=kwargs.get('time_frame'),
            show_type=kwargs.get('show_type', 'both')
        )
    
    # Histogram plot implementations
    def _plot_hist_complex_size(self, data, plot_config, legend, kwargs):
        """Plot complex size histogram."""
        return plot_config.histogram_complex_size(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_frame=kwargs.get('time_frame'),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def _plot_hist_monomer_count(self, data, plot_config, legend, kwargs):
        """Plot monomer count histogram."""
        return plot_config.histogram_monomer_count(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_frame=kwargs.get('time_frame'),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False),
            show_type=kwargs.get('show_type', 'both')
        )
    
    # 3D histogram implementations
    def _plot_3dhist_complex_size(self, data, plot_config, legend, kwargs):
        """Plot 3D complex size histogram."""
        return plot_config.histogram_3d_complex_size(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_bins=kwargs.get('time_bins', 10),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False)
        )
    
    def _plot_3dhist_monomer_count(self, data, plot_config, legend, kwargs):
        """Plot 3D monomer count histogram."""
        return plot_config.histogram_3d_monomer_count(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_bins=kwargs.get('time_bins', 10),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False)
        )
    
    # Heatmap implementations
    def _plot_heatmap_complex_size_time(self, data, plot_config, legend, kwargs):
        """Plot complex size vs time heatmap."""
        return plot_config.heatmap_complex_size_time(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_bins=kwargs.get('time_bins', 10),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False)
        )
    
    def _plot_heatmap_monomer_count_time(self, data, plot_config, legend, kwargs):
        """Plot monomer count vs time heatmap."""
        return plot_config.heatmap_monomer_count_time(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_bins=kwargs.get('time_bins', 10),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False)
        )
    
    def _plot_heatmap_species_correlation(self, data, plot_config, legend, kwargs):
        """Plot species A vs species B heatmap."""
        return plot_config.heatmap_species_correlation(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_bins=kwargs.get('time_bins', 10),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False)
        )
    
    # Stacked histogram implementation
    def _plot_stacked_complex_size(self, data, plot_config, legend, kwargs):
        """Plot stacked complex size histogram."""
        return plot_config.histogram_stacked(
            data=data,
            legend=legend,
            bins=kwargs.get('bins', 10),
            time_frame=kwargs.get('time_frame'),
            frequency=kwargs.get('frequency', False),
            normalize=kwargs.get('normalize', False),
            show_type=kwargs.get('show_type', 'both')
        )
    
    def visualize_trajectory(self, **kwargs):
        """
        Legacy trajectory visualization with enhanced error handling.
        
        Parameters:
            trajectory_path (str, optional): Path to XYZ trajectory file
            save_gif (bool): Whether to save as GIF
            gif_name (str): Output GIF filename
            fps (int): Frames per second
            **kwargs: Additional visualization parameters
        """
        import warnings
        import tempfile
        import os
        import sys

        # Ignore OVITO warning
        warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
        
        try:
            from ovito.io import import_file
            from ovito.vis import Viewport
            import imageio 
            from PIL import Image
        except ImportError:
            msg = (
                "OVITO, imageio, and/or Pillow are required for trajectory visualization but not found."
                "These are optional dependencies. Please install them to enable this feature."
                "If using pip, you can install them with: pip install ionerdss[ovito_rendering]"
                "If using Conda, ensure ovito, imageio, and pillow are installed, for example from conda-forge and conda.ovito.org:"
                "  conda install -c conda.ovito.org -c conda-forge ovito imageio pillow"
            )
            raise ImportError(msg)
        
        # Extract parameters
        trajectory_path = kwargs.get('trajectory_path', None)
        save_gif = kwargs.get('save_gif', False)
        gif_name = kwargs.get('gif_name', 'trajectory.gif')
        fps = kwargs.get('fps', 10)
        
        # Find trajectory file if not specified
        if trajectory_path is None:
            if self.analysis.simulation_dirs:
                trajectory_path = os.path.join(
                    self.analysis.simulation_dirs[0], 
                    "DATA", 
                    "trajectory.xyz"
                )
            else:
                raise ValueError("No simulation directories found")
        
        if not os.path.exists(trajectory_path):
            raise FileNotFoundError(f"Trajectory file '{trajectory_path}' not found.")
        
        try:
            # Import trajectory
            pipeline = import_file(trajectory_path)
            pipeline.add_to_scene()
            vp = Viewport(type=Viewport.Type.PERSPECTIVE)
            vp.zoom_all()

            # Create temporary directory for frames
            with tempfile.TemporaryDirectory() as temp_dir_path:
                frame_paths = []

                # Render frames
                for frame in range(pipeline.source.num_frames):
                    output_path = os.path.join(temp_dir_path, f"frame_{frame:04d}.png")
                    vp.render_image(size=(800, 600), filename=output_path, frame=frame)
                    frame_paths.append(output_path)

                # Create GIF
                gif_path = os.path.join(temp_dir_path, "trajectory.gif")
                imageio.mimsave(gif_path, [imageio.imread(frame) for frame in frame_paths], fps=fps)

                try:
                    from IPython.display import display, Image as IPImage
                    # Display GIF
                    display(IPImage(filename=gif_path))
                except ImportError:
                    print("IPython is not available. GIF will not be displayed in this environment.")
                    print(f"GIF is available at: {gif_path}")

                # Save GIF if requested
                if save_gif:
                    if not hasattr(self, 'save_dir') or not self.save_dir:
                        gif_save_destination_dir = os.getcwd() 
                        print(f"Warning: self.save_dir not set. Saving GIF to current directory: {gif_save_destination_dir}")
                    else:
                        gif_save_destination_dir = self.save_dir
                    
                    if os.path.exists(gif_path):
                        import shutil
                        shutil.move(gif_path, gif_save_destination_dir)
                        print(f"Trajectory GIF saved at: {gif_save_destination_dir}")
                    else:
                        print(f"Error: Temporary GIF file {gif_path} not found. Cannot save.")
                elif not save_gif:
                    if 'IPython.display' not in sys.modules:
                        print(f"Note: GIF was created at {gif_path} but not saved and IPython is not available for display. It will be deleted.")
            
            return gif_path if not save_gif else gif_save_destination_dir
            
        except Exception as e:
            print(f"Error during trajectory visualization: {e}")
            # Cleanup on error
            if 'temp_dir' in locals():
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)
            raise
    
    def _warn_deprecated(self, old_method: str, new_method: str = None):
        """Issue deprecation warning for legacy methods."""
        if new_method:
            message = f"{old_method} is deprecated. Use {new_method} instead."
        else:
            message = f"{old_method} is deprecated."
        
        warnings.warn(message, DeprecationWarning, stacklevel=3)
    
    def get_legacy_help(self) -> str:
        """
        Get help text for legacy interface usage.
        
        Returns:
            str: Formatted help text with examples
        """
        help_text = """
Legacy Plot Interface Help
==========================

The legacy plot_figure interface is maintained for backward compatibility.
For new code, consider using the modular API:

# Legacy style (still supported):
analysis.plot_figure(figure_type="line", x="time", y="species", legend=[["A"]])

# New modular style (recommended):
data = analysis.load_simulations(simulations=[0,1,2])
plot = analysis.set_plot(figure_size=(12,8))
plot.line_speciescopy_vs_time(data=data, legend=[["A"]])

Supported figure_type combinations:
- line + time + species: Species copy numbers over time
- line + time + count: Complex counts over time  
- line + size + free_energy: Free energy landscape
- hist + size + complex_count: Complex size histogram
- heatmap + size + time + complex_count: Size vs time heatmap
- 3dhist + size + time + complex_count: 3D size vs time histogram

For full documentation, see the new API methods.
        """
        return help_text.strip()
    
    def validate_legacy_parameters(self, figure_type: str, **kwargs) -> bool:
        """
        Validate legacy plot_figure parameters.
        
        Parameters:
            figure_type (str): Type of figure
            **kwargs: Plot parameters to validate
            
        Returns:
            bool: True if parameters are valid
            
        Raises:
            ValueError: If parameters are invalid
        """
        required_params = {
            'line': ['legend'],
            'hist': ['legend'],
            'heatmap': ['legend'],
            '3dhist': ['legend'],
            'stacked': ['legend']
        }
        
        if figure_type not in required_params:
            raise ValueError(f"Unsupported figure_type: {figure_type}")
        
        # Check required parameters
        for param in required_params[figure_type]:
            if param not in kwargs or kwargs[param] is None:
                raise ValueError(f"Parameter '{param}' is required for figure_type '{figure_type}'")
        
        # Validate legend format
        legend = kwargs.get('legend')
        if legend is not None:
            if not isinstance(legend, list):
                raise ValueError("Legend must be a list")
            
            # Check for nested list structure for certain plot types
            if figure_type == 'line' and kwargs.get('y') == 'species':
                if not all(isinstance(item, list) for item in legend):
                    raise ValueError("For species plots, legend must be a list of lists")
        
        return True
