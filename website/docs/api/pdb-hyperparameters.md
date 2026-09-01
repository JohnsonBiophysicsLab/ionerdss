# PDBModelHyperparameters

`PDBModelHyperparameters` is the main configuration object for the PDB-to-NERDSS pipeline. It controls structure parsing thresholds, template generation behavior, exported NERDSS settings, optional ProAffinity prediction, ODE generation, and transition-matrix output.

## Typical usage

```python
from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

hyperparams = PDBModelHyperparameters(
    interface_detect_distance_cutoff=0.8,
    chain_grouping_matching_mode="sequence",
    generate_visualizations=True,
    generate_nerdss_files=True,
    ode_enabled=True,
    count_transition=True,
)

builder = PDBModelBuilder("6bno", hyperparams=hyperparams)
system = builder.build_system(workspace_path="6bno_dir")
```

## Parameter groups

### Core detection parameters

#### `interface_detect_distance_cutoff`

- type: `float`
- default: `0.9`
- units: `nm`
- purpose: contact search radius used when detecting candidate interfaces from atomic positions

Increase this when interface detection is too conservative. Decrease it when too many weak contacts are being merged into interfaces.

#### `interface_detect_n_residue_cutoff`

- type: `int`
- default: `2`
- units: `residues`
- purpose: minimum number of contacting residues required on each chain before an interface is accepted

Higher values make interface calls stricter.

#### `min_chain_length`

- type: `int`
- default: `4`
- units: `residues`
- purpose: filters out very small chains or ligands before the model-building pipeline proceeds

### Interface type assignment parameters

#### `interface_type_assignment_distance_threshold`

- type: `float`
- default: `2.0`
- units: `Å`
- purpose: merge interfaces into the same interface type when their distances are sufficiently similar

#### `interface_type_assignment_angle_threshold`

- type: `float`
- default: `0.2`
- units: `radians`
- purpose: angular tolerance used when grouping interfaces into the same type

#### `split_spatial_interface_patches`

- type: `bool`
- default: `False`
- purpose: split interfaces into spatial patches to speed up template building in some systems

### Chain grouping parameters

#### `chain_grouping_rmsd_threshold`

- type: `float`
- default: `2.0`
- units: `Å`
- purpose: RMSD threshold for deciding whether chains are structurally repeated

#### `chain_grouping_seq_threshold`

- type: `float`
- default: `0.5`
- purpose: sequence identity threshold used in sequence-based chain grouping

#### `chain_grouping_custom_aligner`

- type: `Optional[PairwiseAligner]`
- default: `None`
- purpose: provide a custom Biopython aligner instead of using the built-in default global aligner

If omitted, the class creates a default `PairwiseAligner` in `__post_init__`.

#### `chain_grouping_matching_mode`

- type: `"default" | "sequence" | "structure"`
- default: `"default"`
- purpose: choose the repeated-chain grouping strategy

Modes:

- `default`: use mmCIF/header information with sequence fallback
- `sequence`: use sequence identity comparisons
- `structure`: use structural superposition

### Steric clash detection

#### `steric_clash_mode`

- type: `"off" | "auto" | "custom"`
- default: `"off"`
- purpose: control whether steric clash detection is disabled, automatic, or user-specified

### Template building parameters

#### `signature_precision`

- type: `int`
- default: `6`
- units: decimal places
- purpose: normalize geometric signatures with consistent precision to avoid floating-point instability

#### `homodimer_distance_threshold`

- type: `float`
- default: `0.5`
- units: `nm`
- purpose: distance cutoff for homodimer detection logic

#### `homodimer_angle_threshold`

- type: `float`
- default: `0.5`
- units: `radians`
- purpose: angular cutoff for homodimer detection logic

### Homotypic detection parameters

#### `homotypic_detection`

- type: `"auto" | "signature" | "off"`
- default: `"auto"`
- purpose: choose how homotypic binding is identified

#### `homotypic_detection_residue_similarity_threshold`

- type: `float`
- default: `0.7`
- purpose: minimum residue-level similarity used when checking homotypic candidates

#### `homotypic_detection_interface_radius`

- type: `float`
- default: `8.0`
- units: `Å`
- purpose: search radius for homotypic interface detection

### Geometric regularization

#### `is_on_sphere`

- type: `bool`
- default: `False`
- purpose: project molecules onto a best-fit sphere, useful for spherical assemblies

#### `template_regularization_strength`

- type: `float`
- default: `0.0`
- purpose: regularization strength for template fitting

### Output generation

#### `generate_visualizations`

- type: `bool`
- default: `True`
- purpose: create visualization outputs in the workspace

#### `generate_nerdss_files`

- type: `bool`
- default: `True`
- purpose: export `.mol` files and `parms.inp` for NERDSS

### NERDSS export parameters

#### `nerdss_water_box`

- type: `list`
- default: `[500.0, 500.0, 500.0]`
- units: `nm`
- purpose: simulation box dimensions used when exporting NERDSS input files

#### `nerdss_total_molecule_count`

- type: `int`
- default: `75`
- units: molecules
- purpose: target total molecule count distributed across molecule types according to stoichiometry

#### `nerdss_time_step`

- type: `Optional[float]`
- default: `None`
- units: `us`
- purpose: explicit NERDSS time step; if `None`, the pipeline computes one automatically

#### `nerdss_n_itr`

- type: `int`
- default: `100000`
- units: steps
- purpose: number of NERDSS iterations written into the exported parameters

#### `nerdss_overlap_sep_limit`

- type: `float`
- default: `2.0`
- units: `nm`
- purpose: minimum allowed separation distance between molecule centers to avoid overlap-related artifacts

### ProAffinity binding energy prediction

#### `predict_affinity`

- type: `bool`
- default: `False`
- purpose: enable ProAffinity-GNN inference during model generation

#### `adfr_path`

- type: `Optional[str]`
- default: `None`
- purpose: path to the ADFR `prepare_receptor` tool if auto-detection is not sufficient

#### `proaffinity_backend`

- type: `str`
- default: `"auto"`
- purpose: where ProAffinity runs. `"auto"` uses the separate environment when one is configured and this interpreter otherwise, `"sidecar"` requires one, `"in_process"` never spawns one. ProAffinity pins numpy 1.x, so it usually lives in its own environment

#### `proaffinity_python`

- type: `Optional[str]`
- default: `None`
- purpose: interpreter of the environment ProAffinity is installed in. Defaults to `$IONERDSS_PROAFFINITY_PYTHON`

### Structure source options

#### `pdb_file_format`

- type: `str`
- default: `"bioassembly1"`
- purpose: choose how remote structures are fetched from the PDB

Examples include:

- `pdb`
- `cif`
- `mmcif`
- `bioassembly1`
- `bioassembly2`

### ODE pipeline parameters

#### `ode_enabled`

- type: `bool`
- default: `False`
- purpose: enable generation and solving of the ODE reaction system

#### `max_complex_size_ode`

- type: `int`
- default: `12`
- purpose: skip ODE generation when the assembly size becomes too large

#### `ode_time_span`

- type: `Optional[tuple]`
- default: `None`
- units: seconds
- purpose: explicit ODE integration interval `(start, end)`

If omitted, the code derives a time span from the NERDSS time step and iteration count.

#### `ode_solver_method`

- type: `str`
- default: `"BDF"`
- purpose: solver used for ODE integration

#### `ode_atol`

- type: `float`
- default: `1e-4`
- purpose: absolute tolerance for the ODE solver

#### `ode_plot`

- type: `bool`
- default: `True`
- purpose: save plots for ODE results

#### `ode_save_csv`

- type: `bool`
- default: `True`
- purpose: save CSV output for ODE results

#### `ode_initial_concentrations`

- type: `Optional[dict]`
- default: `None`
- purpose: explicit initial concentrations by species name

If omitted, the pipeline uses concentrations derived from the NERDSS export assumptions.

#### `ode_plot_species_indices`

- type: `Optional[list]`
- default: `None`
- purpose: restrict plotting to selected species indices

#### `ode_plot_sample_points`

- type: `int`
- default: `1000`
- purpose: number of sampled time points used in plotting

#### `ode_species_labels`

- type: `Optional[dict]`
- default: `None`
- purpose: custom labels for plotted species indices

### Kinetic parameter defaults

#### `default_on_rate_3d_ka`

- type: `float`
- default: `120.0`
- units: `nm^3/us`
- purpose: default 3D association rate used for diffusion-limited reactions

### Transition matrix output

#### `count_transition`

- type: `bool`
- default: `False`
- purpose: enable transition matrix collection during simulation

#### `transition_matrix_size`

- type: `Optional[int]`
- default: `None`
- purpose: set the dimensions of the tracked transition matrix, i.e. the largest number of copies of a single molecule type NERDSS can track inside one complex

If left as `None`, each molecule type gets a matrix sized to its own molecule count. NERDSS does not bounds check this value, so a complex holding more copies of a molecule type than `transition_matrix_size` corrupts memory and crashes the simulation. When `count_transition=True`, the export raises a `ValueError` if an explicit `transition_matrix_size` is smaller than the number of copies of a molecule type in the simulation.

#### `transition_write`

- type: `Optional[int]`
- default: `None`
- purpose: interval for writing transition matrix output files

If omitted, the export code defaults this interval to roughly `nItr / 10`.

### Units

#### `units`

- type: `Units`
- default: `Units()`
- purpose: internal unit system object used by the pipeline

This field is intentionally skipped in JSON serialization.

## Helper methods on the dataclass

### `to_dict()`

Serialize the hyperparameters into a regular dictionary. The custom aligner is converted into a simple parameter dictionary, and `units` is omitted.

### `from_dict(data)`

Reconstruct a `PDBModelHyperparameters` instance from serialized data. This also restores tuple handling for `ode_time_span` and rebuilds a `PairwiseAligner` when aligner settings are provided.

### `validate()`

Validate the current hyperparameter values and report invalid combinations or ranges.
