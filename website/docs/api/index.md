# API Reference

The current public API is concentrated around a small set of entry points:

- `ionerdss.build_system_from_pdb` for the high-level structure-to-system workflow
- `ionerdss.prepare_structure_validation_for_system` for exporting the one-copy validation setup
- `ionerdss.compare_structure_to_design` for rigid-alignment RMSD checks against the design target
- `ionerdss.System` for the serialized molecular system container
- `ionerdss.Simulation` for editing and running NERDSS input workspaces
- `ionerdss.Analyzer` for post-processing simulation outputs
- `ionerdss.ODEPipelineConfig` and `ionerdss.run_ode_pipeline` for kinetic precomputation
- `ionerdss.visualize_trajectory_ovito` for OVITO-based trajectory rendering
- `ionerdss.convert_simularium` for Simularium export

Additional top-level dataclasses exposed for structure validation include:

- `ionerdss.StructureValidationConfig`
- `ionerdss.StructureValidationArtifacts`
- `ionerdss.StructureAlignmentResult`

Use the pages in this section for the currently supported import surface rather than the removed Read the Docs layout.
