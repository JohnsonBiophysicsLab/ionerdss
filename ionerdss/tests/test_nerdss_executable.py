"""Tests for how ``nerdss_path`` is resolved to the NERDSS executable that gets run."""

import os

import pytest

from ionerdss.nerdss_simulation import Simulation, resolve_nerdss_executable


def _make_executable(path, mode=0o755):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n", encoding="utf-8")
    path.chmod(mode)
    return path


def test_executable_file_is_used_as_is(tmp_path):
    exe = _make_executable(tmp_path / "custom" / "my_nerdss_build")
    assert resolve_nerdss_executable(exe) == exe
    assert resolve_nerdss_executable(str(exe)) == exe


@pytest.mark.parametrize(
    "layout, pass_dir, expected",
    [
        (["bin/nerdss"], ".", "bin/nerdss"),  # NERDSS checkout
        (["bin/nerdss"], "bin", "bin/nerdss"),  # its bin/ folder
        (["nerdss", "bin/nerdss"], ".", "nerdss"),  # ./ before ./bin/
        (["bin/nerdss", "bin/nerdss_mpi"], ".", "bin/nerdss"),  # serial before MPI
        (["nerdss_mpi", "bin/nerdss"], ".", "bin/nerdss"),  # serial anywhere before MPI
        (["bin/nerdss_mpi"], ".", "bin/nerdss_mpi"),  # MPI-only build
    ],
)
def test_directory_is_searched(tmp_path, layout, pass_dir, expected):
    for relative in layout:
        _make_executable(tmp_path / relative)
    assert resolve_nerdss_executable(tmp_path / pass_dir) == tmp_path / expected


def test_tilde_is_expanded(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    exe = _make_executable(tmp_path / "NERDSS" / "bin" / "nerdss")
    assert resolve_nerdss_executable("~/NERDSS") == exe


def test_bare_command_name_is_looked_up_on_path(tmp_path, monkeypatch):
    exe = _make_executable(tmp_path / "bin" / "nerdss_mpi")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    monkeypatch.chdir(tmp_path)
    assert resolve_nerdss_executable("nerdss_mpi") == exe


def test_none_prefers_default_root_then_path(tmp_path, monkeypatch):
    on_path = _make_executable(tmp_path / "path_bin" / "nerdss")
    monkeypatch.setenv("PATH", str(tmp_path / "path_bin"))

    assert resolve_nerdss_executable() == on_path
    assert resolve_nerdss_executable(default_root=tmp_path / "missing") == on_path

    installed = _make_executable(tmp_path / "NERDSS" / "bin" / "nerdss")
    assert resolve_nerdss_executable(default_root=tmp_path / "NERDSS") == installed


def test_non_executable_file_raises_permission_error(tmp_path):
    exe = _make_executable(tmp_path / "nerdss", mode=0o644)
    with pytest.raises(PermissionError, match="not executable"):
        resolve_nerdss_executable(exe)


def test_missing_executable_lists_every_location_tried(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    _make_executable(tmp_path / "bin" / "nerdss", mode=0o644)

    with pytest.raises(FileNotFoundError) as excinfo:
        resolve_nerdss_executable(tmp_path)
    message = str(excinfo.value)
    for relative in ("nerdss", "bin/nerdss", "nerdss_mpi", "bin/nerdss_mpi"):
        assert str(tmp_path / relative) in message
    assert "exists but is not executable" in message

    with pytest.raises(FileNotFoundError, match="nerdss on PATH"):
        resolve_nerdss_executable()

    with pytest.raises(FileNotFoundError, match="does_not_exist"):
        resolve_nerdss_executable(tmp_path / "does_not_exist")


class _FakeProcess:
    def wait(self):
        return 0

    def poll(self):
        return 0


def _capture_popen(monkeypatch):
    calls = []

    def fake_popen(cmd, **kwargs):
        calls.append((cmd, kwargs["cwd"]))
        return _FakeProcess()

    monkeypatch.setattr("ionerdss.nerdss_simulation.simulation.subprocess.Popen", fake_popen)
    return calls


def _work_dir(tmp_path):
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    (work_dir / "parms.inp").write_text("parms", encoding="utf-8")
    return work_dir


def test_run_new_simulations_runs_executable_in_place(tmp_path, monkeypatch):
    exe = _make_executable(tmp_path / "NERDSS" / "bin" / "nerdss_mpi")
    calls = _capture_popen(monkeypatch)

    Simulation(str(_work_dir(tmp_path))).run_new_simulations(
        sim_dir=str(tmp_path / "out"), nerdss_path=str(tmp_path / "NERDSS"), progress=False, verbose=False
    )

    (cmd, cwd), = calls
    assert cmd == [str(exe), "-f", "parms.inp"]
    assert cwd == str(tmp_path / "out" / "1")
    assert sorted(os.listdir(cwd)) == ["output.log", "parms.inp"]


def test_run_new_simulations_defaults_to_work_dir_install(tmp_path, monkeypatch):
    work_dir = _work_dir(tmp_path)
    exe = _make_executable(work_dir / "NERDSS" / "bin" / "nerdss")
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    calls = _capture_popen(monkeypatch)

    Simulation(str(work_dir)).run_new_simulations(sim_dir=str(tmp_path / "out"), progress=False, verbose=False)

    assert calls[0][0][0] == str(exe)


def test_run_restart_simulations_runs_executable_in_place(tmp_path, monkeypatch):
    exe = _make_executable(tmp_path / "nerdss")
    restart = tmp_path / "out" / "1" / "DATA" / "restart.dat"
    restart.parent.mkdir(parents=True)
    restart.write_text("restart", encoding="utf-8")
    calls = _capture_popen(monkeypatch)

    Simulation(str(_work_dir(tmp_path))).run_restart_simulations(
        sim_dir=str(tmp_path / "out"), nerdss_path=str(exe)
    )

    (cmd, cwd), = calls
    assert cmd == [str(exe), "-r", "restart.dat"]
    assert sorted(os.listdir(cwd)) == ["output.log", "restart.dat"]


def test_nerdss_dir_is_a_deprecated_alias(tmp_path, monkeypatch):
    exe = _make_executable(tmp_path / "NERDSS" / "bin" / "nerdss")
    calls = _capture_popen(monkeypatch)
    simulation = Simulation(str(_work_dir(tmp_path)))

    with pytest.warns(DeprecationWarning, match="nerdss_path"):
        simulation.run_new_simulations(
            sim_dir=str(tmp_path / "out"), nerdss_dir=str(tmp_path / "NERDSS"), progress=False, verbose=False
        )
    assert calls[0][0][0] == str(exe)

    with pytest.raises(TypeError, match="both nerdss_path and nerdss_dir"):
        simulation.run_new_simulations(
            sim_dir=str(tmp_path / "out"), nerdss_path=str(exe), nerdss_dir=str(tmp_path / "NERDSS"),
            progress=False, verbose=False,
        )


def test_validation_run_simulation_forwards_nerdss_path(tmp_path, monkeypatch):
    from ionerdss.model.pdb import validation

    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(validation, "run_structure_validation_simulation", fake_run)

    validation.run_simulation(None, "~/NERDSS")
    assert captured["nerdss_path"] == "~/NERDSS"

    with pytest.warns(DeprecationWarning, match=r"run_simulation\(nerdss_dir=\.\.\.\)"):
        validation.run_simulation(None, nerdss_dir="~/old")
    assert captured["nerdss_path"] == "~/old"
