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


def test_status_from_partial_builder_returns_fp_when_partial_coarse_summary_has_too_few_chains():
    builder = SimpleNamespace(
        coarse_summary={"num_chains": 1},
        group_summary={"num_groups": 0},
        system=None,
    )

    assert benchmark._status_from_partial_builder(builder) == "FP"


def test_status_from_partial_builder_returns_dc_when_partial_system_is_disconnected():
    builder = SimpleNamespace(
        coarse_summary={"num_chains": 3},
        group_summary={"num_groups": 3},
        system=SimpleNamespace(),
    )

    original = benchmark.get_disconnected_design_message
    benchmark.get_disconnected_design_message = lambda system, prefix: "Validation preflight warning: disconnected"
    try:
        assert benchmark._status_from_partial_builder(builder) == "DC"
    finally:
        benchmark.get_disconnected_design_message = original


def _site_geometry_system(site_offsets):
    """One molecule type with two declared interfaces; instances bound through the given site offsets."""
    import numpy as np

    from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
    from ionerdss.model.components.system import System
    from ionerdss.model.components.types import MoleculeType

    molecule_type = MoleculeType(name="A", interfaces_neighbors_map={"AA1f": None, "AA1b": None})
    instances = [
        MoleculeInstance(
            name=name,
            molecule_type=molecule_type,
            com=np.asarray(com, dtype=float),
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 1.0, 0.0]),
        )
        for name, com in (("A_A", [0.0, 0.0, 0.0]), ("B_A", [1.0, 0.0, 0.0]))
    ]
    for instance, offsets in zip(instances, site_offsets):
        partner = instances[1] if instance is instances[0] else instances[0]
        for index, offset in enumerate(offsets, start=1):
            interface = InterfaceInstance(
                absolute_coord=instance.com + np.asarray(offset, dtype=float),
                this_mol=instance,
                this_mol_name=instance.name,
                partner_mol_name=partner.name,
                interface_index=index,
            )
            instance.interfaces_neighbors_map[interface] = partner

    system = System(workspace_path=".")
    system.molecule_types.add(molecule_type)
    for instance in instances:
        system.molecule_instances.add(instance)
    return system


def test_site_geometry_status_returns_ic_for_coincident_sites():
    system = _site_geometry_system([[(0.01, 0.0, 0.0), (0.01, 0.0, 0.0)], [(0.0, 0.5, 0.0), (0.5, 0.0, 0.0)]])

    status, message = benchmark._site_geometry_status(system, 0.15)

    assert status == "IC"
    assert "cannot resolve the binding angles for molecule type A" in message


def test_site_geometry_status_only_reports_a_lone_site_near_the_com():
    system = _site_geometry_system([[(0.05, 0.0, 0.0), (0.0, 0.8, 0.0)], [(0.0, 0.5, 0.0), (0.5, 0.0, 0.0)]])

    status, message = benchmark._site_geometry_status(system, 0.15)

    assert status is None
    assert message.startswith("Validation preflight warning: 1 reacting interface site within 0.15 nm")


def test_site_geometry_status_is_silent_for_a_healthy_design():
    system = _site_geometry_system([[(0.6, 0.0, 0.0), (0.0, 0.8, 0.0)], [(0.0, 0.5, 0.0), (0.5, 0.0, 0.0)]])

    assert benchmark._site_geometry_status(system, 0.15) == (None, None)
    assert benchmark._closest_site_com_distance(system) == 0.5


def test_status_from_partial_builder_returns_ic_when_the_partial_system_has_coincident_sites():
    builder = SimpleNamespace(
        coarse_summary={"num_chains": 2},
        group_summary={"num_groups": 1},
        system=_site_geometry_system([[(0.01, 0.0, 0.0), (0.01, 0.0, 0.0)], [(0.0, 0.5, 0.0), (0.5, 0.0, 0.0)]]),
        hyperparams=None,
    )

    assert benchmark._status_from_partial_builder(builder) == "IC"


def test_interface_com_proximity_threshold_prefers_the_builder_hyperparameters():
    assert benchmark._interface_com_proximity_threshold(SimpleNamespace(hyperparams=None)) == 0.15
    assert benchmark._interface_com_proximity_threshold(
        SimpleNamespace(hyperparams=SimpleNamespace(interface_com_proximity_threshold=0.3))
    ) == 0.3
