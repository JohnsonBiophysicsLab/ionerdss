"""Tests for bounding and retrying structure downloads.

Neither urllib.request.urlretrieve nor BioPython's PDBList accepts a timeout, and
neither has one by default, so a connection to RCSB that opens and then stalls
blocks its thread for ever. In a worker pool that is a wedged run rather than one
slow download: every worker parks on a socket read and the job sits at 0% CPU.

Nothing here touches the network; the transport is mocked throughout.
"""

import io
import socket
import urllib.error
from pathlib import Path

import pytest

from ionerdss.model.pdb.parser import (
    DOWNLOAD_ATTEMPTS,
    DOWNLOAD_TIMEOUT_SECONDS,
    _socket_timeout,
    download_with_retry,
)


class _Response(io.BytesIO):
    """Minimal stand-in for the object urlopen returns."""

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
        return False


def _urlopen(monkeypatch, behaviour):
    """Install a fake urlopen driven by `behaviour`, and record the timeouts seen."""
    seen = []

    def fake_urlopen(url, timeout=None):
        seen.append(timeout)
        result = behaviour(len(seen))
        if isinstance(result, BaseException):
            raise result
        return _Response(result)

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    return seen


def test_download_passes_a_timeout_to_every_attempt(monkeypatch, tmp_path):
    seen = _urlopen(monkeypatch, lambda n: b"CIF")
    target = tmp_path / "s.cif.gz"

    download_with_retry("https://example/s.cif.gz", target,
                        timeout=12.5, attempts=DOWNLOAD_ATTEMPTS)

    assert target.read_bytes() == b"CIF"
    assert seen == [12.5]


def test_a_stalled_connection_is_retried_then_reported(monkeypatch, tmp_path):
    # TimeoutError is what a socket read raises once the timeout bites.
    seen = _urlopen(monkeypatch, lambda n: TimeoutError("timed out"))
    monkeypatch.setattr("time.sleep", lambda _s: None)

    with pytest.raises(OSError):
        download_with_retry("https://example/s.cif.gz", tmp_path / "s.cif.gz",
                            timeout=1.0, attempts=3)

    assert len(seen) == 3, "every attempt should be spent before giving up"


def test_a_transient_failure_recovers_without_reaching_the_caller(monkeypatch, tmp_path):
    seen = _urlopen(monkeypatch,
                    lambda n: ConnectionResetError("reset") if n == 1 else b"CIF")
    monkeypatch.setattr("time.sleep", lambda _s: None)
    target = tmp_path / "s.cif.gz"

    download_with_retry("https://example/s.cif.gz", target, timeout=1.0, attempts=3)

    assert target.read_bytes() == b"CIF"
    assert len(seen) == 2


def test_a_missing_file_is_not_retried(monkeypatch, tmp_path):
    # 404 is the server saying the assembly does not exist. Retrying cannot change
    # that, and on a 52k-entry batch the wasted backoff is the whole cost.
    err = urllib.error.HTTPError("u", 404, "Not Found", {}, None)
    seen = _urlopen(monkeypatch, lambda n: err)

    with pytest.raises(urllib.error.HTTPError):
        download_with_retry("https://example/s.cif.gz", tmp_path / "s.cif.gz",
                            timeout=1.0, attempts=3)

    assert len(seen) == 1


def test_a_server_error_is_retried(monkeypatch, tmp_path):
    err = urllib.error.HTTPError("u", 503, "Unavailable", {}, None)
    seen = _urlopen(monkeypatch, lambda n: err)
    monkeypatch.setattr("time.sleep", lambda _s: None)

    with pytest.raises(urllib.error.HTTPError):
        download_with_retry("https://example/s.cif.gz", tmp_path / "s.cif.gz",
                            timeout=1.0, attempts=2)

    assert len(seen) == 2


def test_a_failed_attempt_leaves_no_partial_file(monkeypatch, tmp_path):
    target = tmp_path / "s.cif.gz"

    def behaviour(n):
        target.write_bytes(b"half")      # a transfer that died partway
        return TimeoutError("timed out")

    _urlopen(monkeypatch, behaviour)
    monkeypatch.setattr("time.sleep", lambda _s: None)

    with pytest.raises(OSError):
        download_with_retry("https://example/s.cif.gz", target, timeout=1.0, attempts=2)

    assert not target.exists()


def test_socket_timeout_scope_is_restored():
    # PDBList offers no timeout parameter, so the socket default is the only lever;
    # it is process-global, so leaving it changed would affect unrelated callers.
    before = socket.getdefaulttimeout()
    try:
        socket.setdefaulttimeout(None)
        with _socket_timeout(7.0):
            assert socket.getdefaulttimeout() == 7.0
        assert socket.getdefaulttimeout() is None

        socket.setdefaulttimeout(3.0)
        with pytest.raises(RuntimeError):
            with _socket_timeout(7.0):
                raise RuntimeError("boom")
        assert socket.getdefaulttimeout() == 3.0
    finally:
        socket.setdefaulttimeout(before)


def test_bioassembly_download_uses_the_parsers_settings(monkeypatch, tmp_path):
    """The wiring, not just the helper: a stall in the real fetch path is bounded."""
    from ionerdss.model.pdb.parser import PDBParser

    parser = object.__new__(PDBParser)
    parser.workspace_manager = None
    parser.download_timeout = 9.0
    parser.download_attempts = 4

    calls = {}

    def fake_download(url, destination, timeout, attempts, backoff=None, logger=None):
        calls.update(url=url, timeout=timeout, attempts=attempts)
        import gzip
        with gzip.open(destination, "wb") as handle:
            handle.write(b"data_test\n")

    monkeypatch.setattr("ionerdss.model.pdb.parser.download_with_retry", fake_download)

    result = parser._download_bioassembly("1ABC", 1)

    assert calls["timeout"] == 9.0
    assert calls["attempts"] == 4
    assert calls["url"].endswith("1abc-assembly1.cif.gz")
    assert Path(result).read_bytes() == b"data_test\n"


def test_a_stalled_bioassembly_download_raises_rather_than_hanging(monkeypatch):
    from ionerdss.model.pdb.parser import PDBParser

    parser = object.__new__(PDBParser)
    parser.workspace_manager = None
    parser.download_timeout = 0.5
    parser.download_attempts = 1

    def stall(url, destination, timeout, attempts, backoff=None, logger=None):
        raise TimeoutError("timed out")

    monkeypatch.setattr("ionerdss.model.pdb.parser.download_with_retry", stall)

    with pytest.raises(ValueError, match="Failed to download"):
        parser._download_bioassembly("1ABC", 1)


def test_parser_exposes_the_download_settings():
    from ionerdss.model.pdb.parser import PDBParser
    import inspect

    params = inspect.signature(PDBParser.__init__).parameters
    assert params["download_timeout"].default == DOWNLOAD_TIMEOUT_SECONDS
    assert params["download_attempts"].default == DOWNLOAD_ATTEMPTS
