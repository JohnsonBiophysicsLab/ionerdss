"""
ionerdss.model.pdb.main

Main PDB model builder interface with workspace management and NERDSS export.

This module provides the high-level PDBModelBuilder class that orchestrates
the complete pipeline with proper file organization, logging, and NERDSS export.
"""

from typing import Optional, Union, Dict, Any, Tuple
from pathlib import Path

from ionerdss.model.components.system import System
from ionerdss.model.components.units import Units
from .hyperparameters import PDBModelHyperparameters
from .parser import PDBParser
from .coarse_graining import CoarseGrainer
from .chain_grouping import ChainGrouper
from .template_builder import TemplateBuilder
from .system_builder import SystemBuilder
from .file_manager import WorkspaceManager


class PDBModelBuilder:
    """High-level builder for creating ionerdss Systems from PDB files or PDB IDs.

    Orchestrates the complete pipeline with proper workspace management,
    file organization, comprehensive logging, and optional NERDSS export.

    Attributes:
        source: PDB ID or file path.
        fetch_format: Format for downloading structures.
        workspace_manager: Workspace manager for file organization.
        pdb_id: Extracted or provided PDB ID.
        parser: PDB parser instance (created during build).
    """

    def __init__(self, source: Union[str, Path], fetch_format: str = None,
                 hyperparams: PDBModelHyperparameters = None,):
        """Initialize PDB model builder.

        Args:
            source: PDB ID (4 characters) or path to PDB/mmCIF file.
            fetch_format: Format for downloading. If None, uses hyperparams.pdb_file_format. 
                          Kept for backwards compatibility.
        """
        self.source = source # set source from either PDB ID or local path
        self.fetch_format = fetch_format # format, can be None to use hyperparameter
        self.workspace_manager: Optional[WorkspaceManager] = None
        self.pdb_id: Optional[str] = None # pdb_id to be set from source
        self.parser: Optional[PDBParser] = None
        self.hyperparams = hyperparams

    def build_system(self, workspace_path: str,
                     hyperparams: PDBModelHyperparameters = None,
                     molecule_counts: Optional[Dict[str, int]] = None,
                     box_nm: Tuple[float, float, float] = (
                         100.0, 100.0, 100.0),
                     nerdss_params: Optional[Dict[str, Any]] = None,
                     **kwargs) -> System:
        """Build complete ionerdss System from PDB source.

        Args:
            workspace_path: Path for workspace directory.
            distance_cutoff: Contact search radius in nm. Default 0.6.
            residue_cutoff: Minimum contacting residues per chain. Default 3.
            rmsd_threshold: RMSD threshold for structure grouping in Å. Default 2.0.
            seq_threshold: Sequence similarity threshold. Default 0.5.
            matching_mode: Chain grouping mode. Default "default".
            steric_clash_mode: Steric clash detection mode. Default "off".
            units: Unit system. Defaults to standard units.
            generate_visualizations: Whether to generate visualization outputs. Default True.
            generate_nerdss_files: Whether to generate NERDSS simulation files. Default False.
            molecule_counts: Number of molecules per type for NERDSS. Default 10 each.
            box_nm: Simulation box size in nm for NERDSS. Default (100, 100, 100).
            nerdss_params: Additional NERDSS parameters. Default None.
            **kwargs: Additional hyperparameters.

        Returns:
            Complete System object ready for simulation.
        """
        # Extract PDB ID for workspace naming
        if self._looks_like_pdb_id(str(self.source)):
            pdb_id = str(self.source).upper()
        else:
            # Try to extract from filename
            pdb_id = self._extract_pdb_id_from_path(Path(self.source))

        # Create workspace manager
        self.workspace_manager = WorkspaceManager(workspace_path, pdb_id)
        self.pdb_id = pdb_id

        try:
            # Get hyperparameters: use provided, then builder's, then default
            if hyperparams is None:
                hyperparams = self.hyperparams
            if hyperparams is None:
                hyperparams = PDBModelHyperparameters()

            self.workspace_manager.logger.info(
                "Hyperparameters: %s", hyperparams)

            # Use provided units or create default
            units = hyperparams.units
            if units is None:
                units = Units()

            # Step 1: Parse PDB file or fetch from database
            self.workspace_manager.logger.info(
                "Step 1: Processing structure source: %s", self.source)
            
            # Determine file format: use fetch_format if provided, otherwise use hyperparameter
            file_format = self.fetch_format if self.fetch_format is not None else hyperparams.pdb_file_format
            
            self.parser = PDBParser(
                source=self.source,
                units=units,
                file_format=file_format,
                workspace_manager=self.workspace_manager
            )
            self.pdb_id = self.parser.get_pdb_id() or pdb_id

            # Step 2: Coarse-grain structure
            self.workspace_manager.logger.info(
                "Step 2: Detecting interfaces...")
            coarse_grainer = CoarseGrainer(self.parser, hyperparams)

            coarse_summary = coarse_grainer.get_summary()
            self.workspace_manager.logger.info("Found %d interfaces between %d chains",
                                               coarse_summary['num_interfaces'],
                                               coarse_summary['num_chains'])

            # Step 3: Group repeated chains
            self.workspace_manager.logger.info(
                "Step 3: Grouping repeated chains...")
            chain_grouper = ChainGrouper(
                self.parser, coarse_grainer, hyperparams)

            group_summary = chain_grouper.get_summary()
            self.workspace_manager.logger.info("Created %d chain groups using %s method",
                                               group_summary['num_groups'],
                                               group_summary['grouping_method'])

            # Step 4: Build templates
            self.workspace_manager.logger.info(
                "Step 4: Building molecular templates...")
            template_builder = TemplateBuilder(
                parser=self.parser,
                coarse_grainer=coarse_grainer,
                chain_grouper=chain_grouper,
                hyperparams=hyperparams,
                units=units,
                workspace_manager=self.workspace_manager
            )

            template_summary = template_builder.get_summary()
            self.workspace_manager.logger.info("Created %d molecule templates and %d interface templates",
                                               template_summary['num_molecule_templates'],
                                               template_summary['num_interface_templates'])

            # Step 5: Assemble final system
            self.workspace_manager.logger.info(
                "Step 5: Assembling final system...")
            system_builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=coarse_grainer,
                chain_grouper=chain_grouper,
                template_builder=template_builder,
                hyperparams=hyperparams,
                workspace_path=str(self.workspace_manager.workspace_path),
                pdb_id=self.pdb_id,
                units=units,
                workspace_manager=self.workspace_manager
            )

            system = system_builder.get_system()

            # Step 6: Generate visualizations (if requested)
            if hyperparams.generate_visualizations:
                self.workspace_manager.logger.info(
                    "Step 6: Generating visualizations...")
                viz_outputs = system_builder.generate_visualizations()

                for viz_type, viz_path in viz_outputs.items():
                    self.workspace_manager.logger.info(
                        "Generated %s: %s", viz_type, viz_path)

            # Step 7: Export NERDSS files (if requested)
            if hyperparams.generate_visualizations:
                self.workspace_manager.logger.info(
                    "Step 7: Exporting NERDSS simulation files...")

                # Set default molecule counts if not provided
                if molecule_counts is None:
                    molecule_counts = {}
                    for mol_type in system.molecule_types:
                        molecule_counts[mol_type.name] = 10

                # Export NERDSS files
                # Add hyperparameters to parms_overrides for transition matrix config
                if nerdss_params is None:
                    nerdss_params = {}
                nerdss_params['hyperparams'] = hyperparams
                
                # Use water box size from hyperparameters
                box_size = tuple(hyperparams.nerdss_water_box) if hyperparams.nerdss_water_box else box_nm
                
                nerdss_files = system_builder.export_nerdss_files(
                    molecule_counts=molecule_counts,
                    box_nm=box_size,
                    parms_overrides=nerdss_params
                )

                for file_type, file_path in nerdss_files.items():
                    self.workspace_manager.logger.info(
                        "Generated NERDSS file %s: %s", file_type, file_path)

            # Step 7.5: Run ODE pipeline (if enabled)
            if hyperparams.ode_enabled:
                step_num_ode = 8 if hyperparams.generate_nerdss_files else 7
                self.workspace_manager.logger.info(
                    "Step %d: Running ODE pipeline...", step_num_ode)
                
                try:
                    # Import ODE pipeline module and new System-compatible generator
                    from ionerdss.ode_pipeline import run_ode_pipeline, ODEPipelineConfig
                    from ionerdss.system_ode_generator import generate_ode_model_from_system
                    
                    # Generate complex reaction system using new System-compatible function
                    complex_list, complex_reaction_system = generate_ode_model_from_system(system)
                    
                    self.workspace_manager.logger.info(
                        "Generated %d complexes and %d reactions for ODE",
                        len(complex_list), len(complex_reaction_system.reactions))
                    
                    if len(complex_list) == 0:
                        raise ValueError("No complexes generated from system")
                    if len(complex_reaction_system.reactions) == 0:
                        self.workspace_manager.logger.warning("No reactions generated, ODE will have trivial dynamics")
                    
                    # Create ODE configuration from hyperparameters
                    ode_config = ODEPipelineConfig(
                        t_span=hyperparams.ode_time_span,
                        solver_method=hyperparams.ode_solver_method,
                        atol=hyperparams.ode_atol,
                        plot=hyperparams.ode_plot,
                        save_csv=hyperparams.ode_save_csv,
                        initial_concentrations=hyperparams.ode_initial_concentrations
                    )
                    
                    # Create ODE output directory
                    ode_output_dir = self.workspace_manager.workspace_path / "ode_results"
                    
                    # Run ODE pipeline
                    time, concentrations, species_names, saved_files = run_ode_pipeline(
                        complex_reaction_system,
                        ode_output_dir,
                        config=ode_config,
                        filename_prefix="ode_solution"
                    )
                    
                    self.workspace_manager.logger.info(
                        "ODE pipeline completed. Found %d species, solved for %d time points",
                        len(species_names), len(time))
                    
                    for file_type, file_path in saved_files.items():
                        self.workspace_manager.logger.info(
                            "Generated ODE %s: %s", file_type, file_path)
                            
                except Exception as ode_error:
                    self.workspace_manager.logger.warning(
                        "ODE pipeline failed (continuing with normal workflow): %s", str(ode_error))
                    import traceback
                    self.workspace_manager.logger.debug(traceback.format_exc())


            # Step 8: Save system and generate reports
            step_num = 8 if hyperparams.generate_nerdss_files else 7
            self.workspace_manager.logger.info(
                "Step %d: Saving outputs...", step_num)

            # Save system JSON
            system_path = self.workspace_manager.get_system_output_path()
            system.to_json(str(system_path))
            self.workspace_manager.logger.info(
                "Saved system to: %s", system_path)

            # Validate system and save validation report
            validation = system.validate_system()
            validation_path = self.workspace_manager.get_report_path(
                'validation')
            with open(validation_path, 'w') as f:
                f.write("System Validation Report\n")
                f.write("=" * 30 + "\n\n")
                f.write(f"Errors: {len(validation['errors'])}\n")
                for error in validation['errors']:
                    f.write(f"  - {error}\n")
                f.write(f"\nWarnings: {len(validation['warnings'])}\n")
                for warning in validation['warnings']:
                    f.write(f"  - {warning}\n")

            # Save detailed summary
            summary_path = self.workspace_manager.get_report_path(
                'detailed_summary')
            with open(summary_path, 'w') as f:
                f.write("Detailed Pipeline Summary\n")
                f.write("=" * 30 + "\n\n")
                f.write(f"PDB ID: {self.pdb_id}\n")
                f.write(f"Source: {self.source}\n")
                f.write(
                    f"Workspace: {self.workspace_manager.workspace_path}\n\n")

                f.write("Coarse-graining Results:\n")
                for key, value in coarse_summary.items():
                    f.write(f"  {key}: {value}\n")

                f.write("\nChain Grouping Results:\n")
                for key, value in group_summary.items():
                    f.write(f"  {key}: {value}\n")

                f.write("\nTemplate Building Results:\n")
                for key, value in template_summary.items():
                    f.write(f"  {key}: {value}\n")

                f.write("\nFinal System Summary:\n")
                system_summary = system.get_summary()
                for key, value in system_summary.items():
                    f.write(f"  {key}: {value}\n")

                # Add NERDSS export info if generated
                if hyperparams.generate_nerdss_files:
                    f.write(f"\nNERDSS Export:\n")
                    f.write(f"  Molecule counts: {molecule_counts}\n")
                    f.write(f"  Box size (nm): {box_nm}\n")
                    if nerdss_params:
                        f.write(f"  Custom parameters: {nerdss_params}\n")

            self.workspace_manager.logger.info(
                "Pipeline completed successfully!")
            return system

        except Exception as e:
            if self.workspace_manager:
                self.workspace_manager.logger.error(
                    "Pipeline failed: %s", str(e))
            raise e

    def _looks_like_pdb_id(self, source: str) -> bool:
        """Check if source looks like a PDB ID."""
        source_clean = Path(source).stem.upper()
        return (len(source_clean) == 4 and
                source_clean.isalnum() and
                not Path(source).exists())

    def _extract_pdb_id_from_path(self, path: Path) -> Optional[str]:
        """Extract PDB ID from file path."""
        stem = path.stem.upper()
        if len(stem) >= 4:
            potential_id = stem[:4]
            if potential_id.isalnum():
                return potential_id
        return "unknown"

    def get_pdb_id(self) -> Optional[str]:
        """Get the PDB ID."""
        return self.pdb_id

    def get_workspace_path(self) -> Optional[Path]:
        """Get the workspace path."""
        return self.workspace_manager.workspace_path if self.workspace_manager else None

    def cleanup(self) -> None:
        """Clean up workspace temporary files."""
        if self.workspace_manager:
            self.workspace_manager.cleanup_temp_files()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with workspace cleanup."""
        if self.workspace_manager:
            self.workspace_manager.__exit__(exc_type, exc_val, exc_tb)

    def set_hyperparameters(self, **kwargs) -> PDBModelHyperparameters:
        """Set or update hyperparameters for this builder instance.
        
        Convenience method that wraps the API function. See the API function
        documentation for complete parameter descriptions.
        
        Args:
            **kwargs: Hyperparameter field names and values to set or update.
        
        Returns:
            The updated PDBModelHyperparameters instance.
        
        Examples:
            >>> builder = PDBModelBuilder("1ABC")
            >>> builder.set_hyperparameters(
            ...     interface_detect_distance_cutoff=0.8,
            ...     ode_enabled=True
            ... )
        """
        from .api import set_hyperparameters as _set_hyperparameters
        return _set_hyperparameters(self, **kwargs)
    
    def export_hyperparameters(self, filepath: str):
        """Export builder's hyperparameters to JSON file.
        
        Convenience method that wraps the API function.
        
        Args:
            filepath: Path to save JSON file.
        
        Examples:
            >>> builder = PDBModelBuilder("1ABC")
            >>> builder.set_hyperparameters(interface_detect_distance_cutoff=0.8)
            >>> builder.export_hyperparameters("config.json")
        """
        from .api import export_hyperparameters as _export_hyperparameters
        return _export_hyperparameters(self, filepath)
    
    def import_hyperparameters(self, filepath: str) -> PDBModelHyperparameters:
        """Import hyperparameters from JSON file and set on this builder.
        
        Convenience method that wraps the API function.
        
        Args:
            filepath: Path to JSON file containing hyperparameters.
        
        Returns:
            The loaded PDBModelHyperparameters instance.
        
        Examples:
            >>> builder = PDBModelBuilder("1ABC")
            >>> builder.import_hyperparameters("config.json")
        """
        from .api import import_hyperparameters as _import_hyperparameters
        return _import_hyperparameters(self, filepath)
    
    def print_hyperparameters(self) -> str:
        """Print builder's hyperparameters in a human-readable format.
        
        Convenience method that wraps the API function.
        
        Returns:
            String representation of hyperparameters.
        
        Examples:
            >>> builder = PDBModelBuilder("1ABC")
            >>> builder.set_hyperparameters()
            >>> print(builder.print_hyperparameters())
        """
        from .api import print_hyperparameters as _print_hyperparameters
        return _print_hyperparameters(self)
