"""
Sidecar entry point for ProAffinity-GNN binding energy prediction.

Run as a standalone script by ionerdss.model.proaffinity_runner:

    <sidecar-python> _proaffinity_worker.py <request.json>

ProAffinity pins numpy 1.x / torch 2.2, which cannot share an environment with
packages that need numpy 2 (OVITO among them). Running it through this script
lets it live in an interpreter of its own: only JSON crosses the boundary, so
the two dependency stacks never meet.

The request/response contract is deliberately tiny, because the underlying call
already is -- a PDB file plus chain pairs in, one energy per pair out:

    request  {"predictions": [{"pdb_file": str, "chains": "A,B"}, ...],
              "adfr_path": str|null, "model_weights_path": str|null,
              "verbose": bool, "response_path": str}
    response {"energies": [float|null, ...]}          # null == prediction failed

A null energy means that one pair failed; the caller substitutes its default
binding energy. Failures that prevent any prediction exit non-zero instead, and
the parent reports stderr.
"""

import json
import os
import sys

EXIT_BAD_REQUEST = 2


def main(argv):
    if len(argv) != 2:
        print(f"usage: {os.path.basename(argv[0])} <request.json>", file=sys.stderr)
        return EXIT_BAD_REQUEST

    with open(argv[1]) as f:
        request = json.load(f)

    # Import ionerdss from the source tree this worker ships with, so the
    # sidecar interpreter only has to provide ProAffinity's dependencies rather
    # than a second ionerdss installation.
    package_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)

    from ionerdss.model.proaffinity_predictor import (
        predict_proaffinity_binding_energy_batch,
    )

    energies = predict_proaffinity_binding_energy_batch(
        predictions_list=[{"pdb_file": p["pdb_file"], "chains": p["chains"]}
                          for p in request["predictions"]],
        model_weights_path=request.get("model_weights_path"),
        adfr_path=request.get("adfr_path"),
        verbose=request.get("verbose", False),
    )

    # NaN has no JSON representation; null carries the same "this pair failed"
    # meaning back to the parent.
    energies = [None if e is None or e != e else float(e) for e in energies]
    with open(request["response_path"], "w") as f:
        json.dump({"energies": energies}, f)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
