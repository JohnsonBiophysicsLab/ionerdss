# Top-Level API

These are the main objects and helpers exposed from `import ionerdss as ion`.

## Public entry points

- `build_system_from_pdb`: convenience wrapper for the structure-to-system pipeline.
- `System`: top-level molecular system container with registries and JSON serialization.
- `Simulation`: helper for editing and running a NERDSS workspace.
- `Analyzer`: entry point for post-processing simulation outputs.
- `ODEPipelineConfig`: configuration object for ODE solves.
- `run_ode_pipeline`: execute the ODE workflow.
- `platonic_solid_generator`: build a platonic-solid-based model.
- `build_system_from_plat`: convenience helper for platonic solid workflows.
- `prepare_structure_validation_for_system`: export the irreversible one-copy validation setup.
- `compare_structure_to_design`: align observed coordinates against the design target and report RMSD.
- `StructureValidationConfig`: configuration dataclass for the validation workflow.
- `StructureValidationArtifacts`: generated files and metadata from validation export.
- `StructureAlignmentResult`: output object returned by rigid alignment.
- `visualize_trajectory_ovito`: render XYZ trajectories to GIFs with OVITO.
- `convert_simularium`: export simulation outputs to Simularium.

## `build_system_from_pdb`

```python
import ionerdss as ion

system = ion.build_system_from_pdb(
    source="4v6x",
    workspace_path="4v6x_dir",
    ode_enabled=True,
    count_transition=True,
)
```

### Parameters

- `source`: PDB identifier such as `"4v6x"` or a local PDB/mmCIF path.
- `workspace_path`: output directory for generated files. Defaults to `<source>_dir`.
- `fetch_format`: optional remote structure format.
- `molecule_counts`: optional explicit counts for NERDSS export.
- `structure_validation`: if `True`, export a one-copy-per-type validation setup during the build.
- `structure_validation_options`: optional settings forwarded to the validation export workflow.
- `**hyperparams_kwargs`: any field accepted by `PDBModelHyperparameters`.

### Returns

A populated `System` object ready for export, simulation setup, or analysis.

## Structure validation helpers

### `prepare_structure_validation_for_system`

Prepare the special validation deck that exports one representative copy per molecule type, forces off-rates to zero, adds titration behavior, and writes the designed target coordinates for later comparison.

Typical return values are packaged in `StructureValidationArtifacts`, including:

- chosen molecule counts
- designed coarse-grained coordinates
- target JSON file path
- generated NERDSS input files
- `preflight_warning_message`: set when the designed assembly graph is disconnected
- `free_interface_warning_message`: set when the design leaves binding capacity unused
- `interface_com_proximity_warning_message`: set when a reacting interface site sits on
  its molecule's centre of mass, or a molecule type's sites coincide, so NERDSS cannot
  define the binding angles
- `box_fit_warning_message`: set when the designed assembly is larger than the
  simulation box

### Over-assembly preflight warning

A coarse-grained design can leave interfaces unbound. A molecule *type* declares every
interface it can bind through, but a given copy in the deposited assembly only realizes
the contacts observed in the PDB; anything declared but unrealized is a **free slot**.
Because a fresh copy arrives with all of its own interfaces free, that slot can bind
another subunit and grow the assembly past the deposited stoichiometry — over-assembly.

`get_free_interface_capacity(system)` returns
`{molecule type: {interface type: number of copies leaving it free}}`, and
`get_free_interface_message(system, prefix=...)` formats it. Both live in
`ionerdss.model.pdb.structure_validation`. The message is also raised as a
`RuntimeWarning` during validation export and carried on the artifacts.

Free interfaces are not necessarily a defect. A filament such as actin genuinely
nucleates beyond the deposited asymmetric unit, so the warning says a larger structure
*may* form, not that the model is wrong.

Two caveats worth knowing:

- The check is **static**, over the finished design. Over-assembly can also arise
  kinetically: during assembly, partial complexes transiently expose interfaces, and
  with irreversible binding they can fuse before completing. A saturated design is not
  a guarantee.
- Over-assembly is **unobservable at one copy** of the deposited stoichiometry, since
  the largest possible assembly is then the target itself. Supply more copies to test
  for it.

### Interface-at-centre-of-mass preflight warning

NERDSS defines the angles of a bond from the vector between a molecule's centre of
mass (COM) and the reacting interface site: theta is the angle between that vector and
sigma, and phi is the rotation about it. Before enforcing them it orients the molecule
onto its `.mol` template from the same site vectors, which needs a site that is
not on the COM. Coarse-graining places a site at the
centroid of the contacting residues, so when a chain contacts its partners along its
whole length — collagen-like triple helices, a peptide lying in a groove, amyloid
segments — every site falls on the chain's COM, all of the chain's sites coincide, and
NERDSS exits at the first association with `Cannot resolve phi angle ... Exiting`.

`get_near_com_interface_sites(system, threshold_nm=...)` lists every reacting interface
site closer than the threshold to its molecule's COM, and
`get_degenerate_site_layouts(system, threshold_nm=...)` lists the molecule instances
whose sites coincide or all lie within the threshold.
`get_interface_com_proximity_message(system, prefix=..., threshold_nm=...)` formats
both. The threshold is the `interface_com_proximity_threshold` hyperparameter
(default 0.15 nm); the validation helpers also accept it as
`interface_com_proximity_threshold_nm`. The message is raised as a `RuntimeWarning`
by `build_system` and during validation export, and carried on the artifacts.

Two things are worth knowing:

- Proximity alone is a symptom, not the cause. In the x5 benchmark the closest reacting
  site of a crashing entry lies a median 0.07 nm from its COM against 0.72 nm in
  size-matched controls, but 8% of healthy multi-site models also have a site within
  0.15 nm and simulate fine. What NERDSS cannot handle is a multi-site molecule whose
  sites coincide: 92% of the crashing entries have one, no control does. The benchmark
  driver therefore skips a model (status `IC`) only for a degenerate layout, and just
  prints the message for a lone site near the COM.
- A single-interface molecule type is exempt, because the exporter writes `phi = nan`
  for it and NERDSS then skips the phi rotation; such a subunit binds, but with an
  arbitrary orientation.

Setting the `interface_site_placement` hyperparameter to `'auto'` moves the flagged
sites onto the chain surface facing the partner's COM, which gives NERDSS a usable
geometry (see the hyperparameter page).

### Box-fit preflight warning

NERDSS keeps every molecule inside the box: it reflects a complex that reaches a wall,
cancels an association whose product would span the box, and exits with
`Molecule seems outside simulation volume` when a molecule is pushed out anyway.
`get_designed_assembly_extent(system)` returns the bounding-sphere diameters of the
designed assembly and of its largest molecule, measured from the interface sites the
way NERDSS measures its template radius, and
`get_box_fit_message(system, box_nm, prefix=...)` warns when the assembly is larger
than the shortest box edge. The check is static: a filament that keeps growing past
the deposited stoichiometry (see the over-assembly warning) can still leave the box.

### `compare_structure_to_design`

Rigidly align observed coarse-grained coordinates onto the design target and compute RMSD.

This function accepts either:

- dictionaries keyed by molecule type
- ordered coordinate arrays

It returns a `StructureAlignmentResult` containing labels, RMSD, the rotation and translation, and the aligned coordinates.

## Rendering and export helpers

### `visualize_trajectory_ovito`

Render an XYZ trajectory using OVITO and optionally save it as a GIF. This requires the `ovito_rendering` optional extra.

### `convert_simularium`

Convert supported simulation outputs into Simularium files for interactive 3D viewing.
