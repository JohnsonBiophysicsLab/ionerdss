> **SUPERSEDED — do not cite these numbers.**
>
> The full paired rerun (`pdb_benchmark/BENCHMARK_v224_x5_auto.md`) shows neither change
> moves the benchmark. The improvement below is an artifact of two errors:
>
> 1. `poc_run.py` ran only the 100,000-iteration fast probe, while the "before" labels
>    came from the benchmark's full protocol, which reruns at 1,000,000 iterations when
>    the target does not appear — ten times longer for complexes to fuse. The two arms
>    were never comparable.
> 2. The alignment fix corrects `MoleculeInstance.ref1`/`ref2`, which the NERDSS export
>    and validation paths never read, so it cannot change a simulation outcome.
>
> Kept for the record.

# POC: over-assembly after the orientation and regularization fixes

120 PDB IDs sampled at random (seed 20260908) from the 53,377-entry x5 benchmark set,
rerun locally at **5 copies** of the deposited stoichiometry in an 85.5 nm box — the
same conditions as the full x5 run. 119 produced a result in every arm.

Three arms:

| arm | code |
|---|---|
| `original x5` | the labels from the full x5 run: before either fix |
| `auto OFF` | residue-intersection alignment fix only |
| `auto ON` | alignment fix + `geometric_regularization="auto"` |

## Result

| arm | OA | Success | UA |
|---|---|---|---|
| original x5 (old code) | **47.1%** | 48.7% | 4.2% |
| new code, `auto` OFF | **16.8%** | 78.2% | 5.0% |
| new code, `auto` ON | **14.3%** | 78.2% | 7.6% |

Raw per-run outcomes are in `poc_geometric_regularization_results.csv`.

## Reading it

**The alignment fix is what moved the number.** Over-assembly falls 47.1% -> 16.8%,
a 30.3 pp drop, from pairing Cα atoms by residue number instead of requiring equal
counts. Copies that previously kept the representative's frame now get a real
orientation, so a cyclic design is an actual ring rather than a straight polymer.

**`auto` added little here, and what it did add was not success.** OA fell a further
2.5 pp (20 -> 17 assemblies), but Success did not move at all: 93 in both arms. The
three assemblies that stopped over-assembling became **under-assembly**, not success.
Snapping geometry onto exact `Cn` can push an interface far enough that its partner no
longer binds, trading one failure mode for another. Three of 119 is well inside noise,
so treat this as "no measured benefit", not "measured harm".

On the subset where `auto` actually fired (17 assemblies detected as a `Cn` ring), the
original run was 58.8% OA and both new arms are **0%** — the alignment fix alone had
already fixed every one, leaving the regularizer nothing to contribute.

## Consequence

Ship the alignment fix; it is a clear, large win and it also benefits the dihedral and
cubic assemblies the regularizer declines to touch.

Keep `geometric_regularization` defaulting to `"off"`. Its value is unproven on this
sample, and the OA -> UA conversions should be understood before it is turned on by
default. A larger run stratified on `Cn`-detected assemblies would settle it; this
sample contains too few (17) to conclude anything.
