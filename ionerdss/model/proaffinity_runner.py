"""
Backend selection for ProAffinity-GNN binding energy prediction.

ProAffinity pins numpy 1.x and torch 2.2. Anything needing numpy 2 -- OVITO
rendering, most notably -- cannot live in the same environment, and pip
resolves the clash by silently installing broken versions of one side or the
other. Rather than pick a winner, run ProAffinity in a separate interpreter and
pass JSON across (see `_proaffinity_worker.py`).

Backends:
    "in_process"  import ProAffinity here; requires numpy 1.x in this env.
    "sidecar"     run it in another interpreter, named by `python_executable`
                  or $IONERDSS_PROAFFINITY_PYTHON.
    "auto"        sidecar if one is configured, otherwise in-process.

`predict_binding_energies` returns one energy per requested pair, with NaN where
that pair failed -- the same contract the in-process function has always had, so
callers keep their existing fallback to a default binding energy. Failures that
prevent any prediction raise, which callers already treat the same way.
"""

import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile

logger = logging.getLogger(__name__)

_WORKER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "_proaffinity_worker.py")

#: Interpreter to run ProAffinity in when no backend is configured explicitly.
SIDECAR_PYTHON_ENV_VAR = "IONERDSS_PROAFFINITY_PYTHON"

BACKENDS = ("auto", "in_process", "sidecar")

_SIDECAR_HELP = (
    "Create one with:\n"
    "    python -m venv ~/.ionerdss-proaffinity\n"
    "    ~/.ionerdss-proaffinity/bin/pip install \"ioNERDSS[proaffinity]\"\n"
    f"then set {SIDECAR_PYTHON_ENV_VAR}=~/.ionerdss-proaffinity/bin/python "
    "(or pass proaffinity_python)."
)


def resolve_sidecar_python(python_executable=None):
    """Return the configured sidecar interpreter, or None if there isn't one."""
    candidate = python_executable or os.environ.get(SIDECAR_PYTHON_ENV_VAR, "")
    if not candidate:
        return None
    candidate = os.path.expanduser(candidate)
    if not os.path.isfile(candidate):
        raise FileNotFoundError(
            f"ProAffinity sidecar interpreter '{candidate}' does not exist.\n"
            + _SIDECAR_HELP
        )
    return candidate


def predict_binding_energies(predictions_list, adfr_path=None,
                             model_weights_path=None, backend="auto",
                             python_executable=None, verbose=False):
    """Predict binding energies for chain pairs, in-process or via a sidecar.

    Args:
        predictions_list: dicts with 'pdb_file' and 'chains' (e.g. 'A,B'). Extra
            keys are ignored, so callers may keep their own bookkeeping in them.
        adfr_path: Path to ADFR's prepare_receptor.
        model_weights_path: Path to model.pkl (defaults to the packaged one).
        backend: One of BACKENDS.
        python_executable: Sidecar interpreter; overrides the environment.
        verbose: Passed through to the predictor.

    Returns:
        list[float]: energies in kJ/mol, NaN where that pair could not be predicted.
    """
    if backend not in BACKENDS:
        raise ValueError(f"Unknown ProAffinity backend '{backend}'; "
                         f"expected one of {', '.join(BACKENDS)}.")
    if not predictions_list:
        return []

    sidecar = None
    if backend in ("auto", "sidecar"):
        sidecar = resolve_sidecar_python(python_executable)
        if sidecar is None and backend == "sidecar":
            raise RuntimeError(
                "backend='sidecar' but no ProAffinity interpreter is configured.\n"
                + _SIDECAR_HELP
            )
        if sidecar is not None and os.path.samefile(sidecar, sys.executable):
            # Pointing the sidecar at this interpreter defeats the purpose and
            # just pays for a subprocess.
            logger.warning(
                "ProAffinity sidecar interpreter is this interpreter; running in-process."
            )
            sidecar = None

    if sidecar is None:
        return _predict_in_process(predictions_list, adfr_path,
                                   model_weights_path, verbose)
    return _predict_in_sidecar(sidecar, predictions_list, adfr_path,
                               model_weights_path, verbose)


def _predict_in_process(predictions_list, adfr_path, model_weights_path, verbose):
    logger.info("Running ProAffinity in this interpreter.")
    try:
        from .proaffinity_predictor import predict_proaffinity_binding_energy_batch
    except ImportError as e:
        # The usual cause is that ProAffinity's stack was never installed here,
        # because its numpy 1.x pin conflicts with something else in this
        # environment. That is exactly what the sidecar backend is for.
        raise ImportError(
            f"ProAffinity is not importable in this environment ({e}).\n"
            + _SIDECAR_HELP
        ) from e

    return predict_proaffinity_binding_energy_batch(
        predictions_list=predictions_list,
        model_weights_path=model_weights_path,
        adfr_path=adfr_path,
        verbose=verbose,
    )


def _predict_in_sidecar(sidecar, predictions_list, adfr_path,
                        model_weights_path, verbose):
    import numpy as np

    logger.info("Running ProAffinity in a separate interpreter: %s", sidecar)
    work_dir = tempfile.mkdtemp(prefix="ionerdss_proaffinity_")
    request_path = os.path.join(work_dir, "request.json")
    response_path = os.path.join(work_dir, "response.json")
    try:
        with open(request_path, "w") as f:
            json.dump({
                "predictions": [
                    {"pdb_file": os.path.abspath(p["pdb_file"]),
                     "chains": p["chains"]}
                    for p in predictions_list
                ],
                "adfr_path": adfr_path,
                "model_weights_path": model_weights_path,
                "verbose": bool(verbose),
                "response_path": response_path,
            }, f)

        proc = subprocess.run([sidecar, _WORKER, request_path],
                              capture_output=True, text=True)
        if proc.stdout.strip():
            logger.info("ProAffinity output:\n%s", proc.stdout.strip())
        if proc.returncode != 0:
            raise RuntimeError(
                f"The ProAffinity sidecar ({sidecar}) exited with code "
                f"{proc.returncode}.\n{proc.stderr.strip()}"
            )

        with open(response_path) as f:
            energies = json.load(f)["energies"]
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    if len(energies) != len(predictions_list):
        raise RuntimeError(
            f"The ProAffinity sidecar returned {len(energies)} energies for "
            f"{len(predictions_list)} chain pairs."
        )
    return [np.nan if e is None else float(e) for e in energies]
