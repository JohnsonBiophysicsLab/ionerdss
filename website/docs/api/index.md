# API Reference

The current public API is concentrated around a small set of entry points:

- `ionerdss.build_system_from_pdb` for the high-level structure-to-system workflow
- `ionerdss.System` for the serialized molecular system container
- `ionerdss.Simulation` for editing and running NERDSS input workspaces
- `ionerdss.Analyzer` for post-processing simulation outputs
- `ionerdss.ODEPipelineConfig` and `ionerdss.run_ode_pipeline` for kinetic precomputation

The pages in this section render directly from the package docstrings where possible, so they stay aligned with the current code instead of the removed Read the Docs content.
