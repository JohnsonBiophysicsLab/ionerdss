from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace


BENCHMARK_SCRIPT = Path(__file__).resolve().parents[2] / "benchmark" / "run_validation_benchmark.py"
_SPEC = spec_from_file_location("run_validation_benchmark", BENCHMARK_SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
benchmark = module_from_spec(_SPEC)
_SPEC.loader.exec_module(benchmark)


def test_status_for_failed_validation_returns_ua_when_largest_assembly_is_too_small():
    sim_result = SimpleNamespace(
        largest_observed_assembly_size=3,
        warning_message="Validation warning: no full assembly matching target.",
        simulation_dir=Path("/tmp/does_not_exist"),
    )

    assert benchmark._status_for_failed_validation(sim_result, target_assembly_size=4) == "UA"


def test_status_for_failed_validation_returns_oa_when_largest_assembly_meets_or_exceeds_target():
    sim_result = SimpleNamespace(
        largest_observed_assembly_size=4,
        warning_message="Validation warning: no full assembly matching target.",
        simulation_dir=Path("/tmp/does_not_exist"),
    )

    assert benchmark._status_for_failed_validation(sim_result, target_assembly_size=4) == "OA"


def test_status_for_failed_validation_returns_nc_when_warning_contains_crash_signature():
    sim_result = SimpleNamespace(
        largest_observed_assembly_size=2,
        warning_message="NERDSS terminated with Segmentation fault (core dumped).",
        simulation_dir=Path("/tmp/does_not_exist"),
    )

    assert benchmark._status_for_failed_validation(sim_result, target_assembly_size=4) == "NC"


def test_status_for_failed_validation_returns_nc_when_output_log_contains_crash_signature(tmp_path):
    simulation_dir = tmp_path / "validation_output" / "1"
    simulation_dir.mkdir(parents=True)
    (simulation_dir / "output.log").write_text(
        "fatal error: abort trap\nSegmentation fault\n",
        encoding="utf-8",
    )
    sim_result = SimpleNamespace(
        largest_observed_assembly_size=2,
        warning_message="Validation warning: simulation did not produce a full assembly.",
        simulation_dir=simulation_dir,
    )

    assert benchmark._status_for_failed_validation(sim_result, target_assembly_size=4) == "NC"


def test_status_for_failed_validation_falls_back_to_old_failed_assembly_label_without_size_metadata():
    sim_result = SimpleNamespace(
        warning_message="Validation warning: simulation did not produce a full assembly.",
        simulation_dir=Path("/tmp/does_not_exist"),
    )

    assert benchmark._status_for_failed_validation(sim_result, target_assembly_size=4) == "Failed_Assembly"


def test_contains_nerdss_crash_signature_matches_common_runtime_crash_messages():
    assert benchmark._contains_nerdss_crash_signature("Segmentation fault (core dumped)")
    assert benchmark._contains_nerdss_crash_signature("Received signal 11 while running NERDSS")
    assert not benchmark._contains_nerdss_crash_signature("Validation warning: no full assembly was found")
