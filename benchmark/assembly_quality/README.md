# Assembly quality of a generated model

Does a model the PDB pipeline generates actually rebuild, in NERDSS, the assembly it was
generated from? The validation benchmark (`benchmark/run_validation_benchmark.py`) answers
that by asking whether the target N-mer appears. These scripts answer the other half: when
it does not, *what* did the run build instead, and which hyperparameter is responsible.

- `assembly_metrics.py` — scores one NERDSS run. Importable for a single run.
- `sweep.py` — builds each named hyperparameter setting, simulates it over several seeds,
  writes one CSV row per run.
- `configs_6bno.json`, `results_6bno.csv` — the 56 settings and the 340 runs behind the
  6BNO result below. The runs span four `--iterations` values (20,000 / 100,000 / 400,000
  / 1,000,000), so reproducing the whole table means one sweep per value, not one sweep.
  Three settings (`comcap*`) override `com_shift_cap_ang`, which only exists once the
  geometric-regularization work lands; against an older build the sweep records those
  three as build failures and carries on with the other 53. None of the conclusions below
  rests on them.

## What is measured

**Mis-assembly.** When the design is a linear polymer — the deposited subunits lie along
one axis and every copy of an interface type joins subunits a fixed number of positions
apart — each simulated complex is embedded back into that lattice. A bond contradicting
the embedding, or two subunits landing in one lattice position, is a branch or a doubled
strand. Steric overlap, two subunit centres closer than `--clash_nm`, is checked for every
design, lattice or not.

**Stretched bonds.** Bonds whose sites sit far outside the loop-closure window. NERDSS's
rigid-body placement cannot produce one, so these are counted separately rather than mixed
into mis-assembly; see the caveat at the end.

**Geometric closure.** NERDSS places a bond between two complexes exactly, but a bond
*inside* one complex only forms when its two sites already lie within `bindRadSameCom`
(1.1) times the binding radius, and nothing moves to help. So a model whose bond
geometries disagree can never close its rings: subunits keep free sites, and other
subunits bind them in the wrong place. The report composes short steps against the long
one and gives the leftover site gap against that window.

**Topology.** Each setting's site count, reaction count and self-binding reaction count.

## Running it

```bash
python benchmark/assembly_quality/sweep.py \
    --source 6bno \
    --configs benchmark/assembly_quality/configs_6bno.json \
    --nerdss_dir ~/Workspace/NERDSS \
    --workspace_root /tmp/6bno_sweep \
    --output results_6bno.csv \
    --seeds 1 2 3 --iterations 100000
```

The config file maps a name to hyperparameter overrides; everything else stays at the
pipeline defaults. Builds and simulations are kept under `--workspace_root` for inspection.

To score a run you already have:

```python
from assembly_metrics import analyse_run, closure_report
analyse_run("path/to/nerdss_files", "workspace/outputs/systems/6BNO_system.json")
```

## Limits

- The lattice check needs a linear-polymer design. It is derived from the deposited
  structure and verified by re-embedding the design itself, so rings, cages and branched
  designs are rejected rather than mismeasured — those runs still get the overlap and
  stretched-bond checks, and `lattice_available` says which applied.
- The loop-closure window assumes NERDSS's default `bindRadSameCom` of 1.1, which ioNERDSS
  does not export.
- Bond transforms come from isolated dimers in the run itself, so closure is only reported
  for bond types that formed at least one dimer.
- `--clash_nm` defaults to 3.0 nm, below the 4.2 nm nearest-neighbour distance in 6BNO. A
  much larger or smaller assembly needs a different value.

## 6BNO: two hyperparameters decide it

6BNO is an eight-subunit actin filament. The model should have four heterotypic sites —
`aa1f`/`aa1b` across strands, `aa2f`/`aa2b` along the long-pitch helix — and its runs
should grow filaments, not clumps. At the pipeline defaults it gets neither.

| Settings (rest default) | Sites | Runs | Mis-assembled |
|---|---|---|---|
| defaults (0.9 nm cutoff, 0.1 nm limit) | 5 | 3 | 3 |
| `interface_detect_distance_cutoff=1.0` only | 4 | 6 | 6 |
| `nerdss_overlap_sep_limit=3.5` only | 5 | 6 | 6 |
| `homotypic_detection="off"` + limit 3.5 | 4 | 6 | 2 |
| cutoff 1.0 + limit 1.5 | 4 | 13 | 3 |
| cutoff 1.0 + limit 2.0 | 4 | 13 | 1 |
| cutoff 1.0 + limit 2.5 | 4 | 13 | 0 |
| **cutoff 1.0 + limit 3.0** | **4** | **54** | **0** |
| cutoff 1.0 + limit 3.5 / 3.75 | 4 | 13 / 13 | 0 / 0 |
| cutoff 1.05 / 1.1 / 1.15 / 1.2 + limit 3.0 | 4 | 13 each | 0 |

75 copies in a 500 nm box (pipeline defaults). Run counts in the table pool every
iteration length used: 100,000 (125 runs), 1,000,000 (150), 400,000 (40) and 20,000 (25).
NERDSS was built from the local repository on 2026-09-01.

**Why the cutoff.** It decides both the sites and the geometry. At 0.9 nm two of the six
long-pitch contacts land 2.2 Å from the first contact's geometric signature, just past the
2.0 Å merge tolerance; they also pass the homodimer symmetry test (their `d_i`/`d_j`
differ by less than the 5 Å `homodimer_distance_threshold`) and the 0.7 residue-similarity
test, so they become a fifth, self-binding site `aa3` — a contact the deposited filament
does not contain. That geometry also cannot close: two cross-strand steps leave the
long-pitch sites 0.90 nm apart against a 0.80 nm closure window. At 1.0 nm all six merge
into one `aa2f`/`aa2b` pair and the gap (0.84 nm) fits the window (0.96 nm). Fixing only
the site count at 0.9 nm — `homotypic_detection="off"`, a lower
`homodimer_distance_threshold`, or a higher `interface_type_assignment_distance_threshold`,
which all produce the identical model — leaves the closure failure, and those runs branch.

**Why the overlap limit.** NERDSS cancels an association only when it puts two subunit
centres closer than `overlapSepLimit`. Filaments joining out of register are offset by the
leftover closure error, 0.1–2 nm, which the 0.1 nm default lets through. The limit has to
stay under the 4.2 nm cross-strand neighbour distance; `main.py` already caps it at 0.9 ×
the minimum chain COM distance, 3.79 nm here.

Settings that do not matter for 6BNO: `interface_detect_n_residue_cutoff` and
`chain_grouping_seq_threshold` leave the model unchanged; `template_regularization_strength`
makes closure worse (rotation error 8–16° against 2.1°); `geometric_regularization` skips
filaments by design.

About 81–98% of the lattice bonds form at the recommended settings. The rest sit just
outside the loop-closure window, so filaments are clean but not fully zipped.

## Caveat: stretched bonds come from NERDSS, not from the model

At the recommended settings, 5 of 43 runs at 100,000 iterations and 3 of 11 at 1,000,000
still report two clean filaments as a single complex, held together by one bond whose
sites are 5–23σ apart. Two associations between the same pair of complexes fire in the
same timestep; the first moves and merges them, and the second is then carried out by
`associate_box`'s loop-closure branch, which does not move anything and does not re-check
the distance. The serial reaction loop in `EXEs/nerdss.cpp` does not skip molecules whose
complex already moved that step, while the MPI path in `perform_bimolecular_reactions.cpp`
does.

A smaller `nerdss_time_step` makes it rarer without fixing it: over 40 runs each at the
same simulated duration, stretched bonds appeared in 5 runs at the automatic step, 3 at a
quarter of it and 1 at a tenth, for 2.6× and 5.9× the runtime. Mis-assembly stayed at zero
in all three. The fix belongs in NERDSS.
