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
- purpose: contact search radius between Cα atoms used when detecting candidate interfaces

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

#### `include_modified_residues`

- type: `bool`
- default: `True`
- purpose: keep non-standard residues that belong to a polymer chain, such as modified
  amino acids (selenomethionine, D-amino acids) and terminal caps. They are stored as
  HETATM records and often form the binding interface. They are recognised from the
  mmCIF `_entity_poly_seq` header, or else by BioPython's non-standard amino-acid check
  (the only test available for PDB-format files). A kept residue without a Cα atom,
  typically a cap, counts towards the chain's centre of mass but not towards interface
  detection. Set `False` for the older standard-amino-acids-only behaviour

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
- purpose: split the contacts between two chains into spatially separate patches, each its own interface

When `False`, all contacts between two chains form a single interface.

### Interface site placement and preflight

#### `interface_com_proximity_threshold`

- type: `float`
- default: `0.15`
- units: `nm`
- purpose: distance below which a reacting interface site counts as sitting on its
  molecule's centre of mass. NERDSS defines the binding angles theta and phi from the
  COM-to-site vector and orients a molecule onto its template from the same vectors, so
  a multi-interface molecule type whose sites all lie within this distance (or coincide)
  makes NERDSS exit at the first association with `Cannot resolve phi angle`. The
  preflight check in `build_system` and the validation export reports every such site
  as a `RuntimeWarning` (the validation export reads this field only from
  `parms_overrides['hyperparams']`; see
  [Top-Level API](top-level.md#interface-at-centre-of-mass-preflight-warning)).
  Typical of chains that contact a partner along their whole length: collagen-like
  triple helices, peptides in a groove, amyloid segments

#### `interface_site_placement`

- type: `Literal["centroid", "auto"]`
- default: `"centroid"`
- purpose: where an interface site is placed on its chain. `"centroid"` is the mean
  position of the contacting Cα atoms. `"auto"` keeps the centroid except for the sites
  the proximity preflight flags, which are moved onto the chain surface facing the
  partner's centre of mass (the point where the chain's atoms end along the COM-to-COM
  direction). For an elongated chain (longest extent at least three times the next)
  whose partner lies within 45° of its long axis, as the strands of a collagen triple
  helix do, the axial component is dropped so the site sits on the flank facing the
  partner rather than at the rod's tip. The COM-to-site vector is then as long as the
  chain is thick and points at the partner, so NERDSS can orient the molecule and define
  theta; sigma runs along the same line, so phi is undefined for that bond and the
  exporter writes it as `nan`, and the thetas and omega fix the relative orientation.
  Rebuilding the x5-benchmark entries that NERDSS aborted on with this setting made
  1CGD, 4DMT, 1D5M and 1FIP assemble (COM RMSD 0.07, 0.09, 0.00 and 0.2-0.5 nm)

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

- type: `"default" | "sequence" | "structure" | "sequence_structure"`
- default: `"default"`
- purpose: choose the repeated-chain grouping strategy

Modes:

- `default`: use mmCIF/header information with sequence fallback
- `sequence`: use sequence identity comparisons
- `structure`: use structural superposition
- `sequence_structure`: both must pass; separates quasi-equivalent conformers of one sequence

### Steric clash detection

#### `steric_clash_mode`

- type: `"off" | "auto" | "custom"`
- default: `"off"`
- purpose: control whether steric clash detection is disabled, automatic, or user-specified

Only `"auto"` runs clash detection. No hyperparameter carries user-specified clash lists yet, so `"custom"` currently behaves like `"off"`.

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

This value is range-checked by `validate()` but not read by the current pipeline, so changing it has no effect.

### Geometric regularization

#### `geometric_regularization`

- type: `str`
- default: `"off"`
- values: `"off"`, `"auto"`
- purpose: detect the assembly's point group and snap its geometry onto it

  A cyclic homomer built from head-to-tail self-binding interfaces (`AA1f` binds
  `AA1b`) is a polymer subunit that happens to close into a ring. It only closes in
  simulation if the generator transform `T` taking subunit *k* to subunit *k+1*
  satisfies `T**n == I`. When the deposited geometry is only approximately n-fold,
  the closure error accumulates over *n* bonds until the last one falls outside
  NERDSS's binding tolerance, the chain elongates instead of closing, and the run
  reports over-assembly.

  `"auto"` places the subunits at exact `2*pi/n` spacing and — the important part —
  **synthesises each subunit's orientation from the group element** rather than
  recovering it by structural alignment.

  Only a single cyclic ring (`Cn`, n >= 3) is regularized. Filaments are detected and
  deliberately left alone: an actin filament genuinely extends past the deposited
  asymmetric unit, so forcing closure would be wrong. A dihedral point group qualifies
  only when all of its subunits form one n-fold ring, and then only the `Cn` rotation
  about the principal axis is applied; dihedral assemblies built from stacked rings,
  and cubic groups, are detected and reported but not regularized. The detected group
  is recorded on the returned system as `system.symmetry_detection`.

#### `symmetry_fold_tolerance`

- type: `float`
- default: `0.15`
- purpose: accept an n-fold symmetry only if rotating the subunit centres of mass by
  `2*pi/n` maps the set onto itself to within this fraction of the assembly radius

  This is the guard that stops a D2 tetramer from being forced into a C4 ring. Raise
  it to regularize looser assemblies, lower it to be stricter.

#### `com_shift_cap_ang`

- type: `float`
- default: `6.0`
- units: `Å`
- purpose: refuse geometric regularization if it would move any subunit centre of
  mass further than this; also caps template regularization

  Both read the value in Å. Geometric regularization is all or nothing: one subunit
  over the cap leaves the whole assembly as deposited and logs a warning. Template
  regularization instead shortens each member's shift to the cap.

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
- default: `0.1`
- units: `nm`
- purpose: minimum allowed separation distance between molecule centers to avoid overlap-related artifacts

Written to `parms.inp` as `overlapSepLimit`. Unless `disable_overlap_sep_limit_check` is set, `build_system` lowers a value above 0.9 × the smallest distance between chain centres of mass in the structure to that bound and logs a warning.

#### `disable_overlap_sep_limit_check`

- type: `bool`
- default: `False`
- purpose: skip the safety check that caps `nerdss_overlap_sep_limit` at 0.9 × the smallest chain centre-of-mass distance

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

A `fetch_format` passed to `PDBModelBuilder` or `build_system_from_pdb` takes precedence over this field.

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
- purpose: explicit initial concentrations by species name, in µM

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

If omitted, `transitionWrite` is left out of `parms.inp` and NERDSS falls back to `nItr / 10`.

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

Reconstruct a `PDBModelHyperparameters` instance from serialized data. This also restores tuple handling for `ode_time_span` and rebuilds a `PairwiseAligner` when aligner settings are provided. Keys that are not fields are ignored, which is also why a misspelled keyword passed to `PDBModelBuilder.build_system(**kwargs)` has no effect.

### `validate()`

Check the current values and return a list of error messages, empty when every value is valid. It does not raise.
