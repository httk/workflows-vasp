"""Shared fixtures for the workflows-vasp test suite.

Ported from httk-workflow's ``tests/conftest.py`` (``git show
15b0b36:tests/conftest.py``), keeping only what these package-directory tests
need: workspace/config isolation, the machine-owned workspace registry
helper, and the normal/extended test-depth knob. The remote-adapter stand-ins
(fake ``ssh``/``sbatch``) that file also carried are unrelated to VASP
workflows and were dropped.
"""

import logging
import os
from collections.abc import Iterator
from dataclasses import dataclass

import pytest
from httk.workflow.registry import register_workspace

# Several tests import httk.atomistic (NumPy) while short-lived runner
# processes are spawned; keeping each BLAS/OMP runtime to one thread avoids
# multiplying that across many concurrent runners.
for _thread_limit in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_thread_limit, "1")


@pytest.fixture(autouse=True)
def _isolated_httk_config(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """Give every test its own httk config and data home.

    This keeps the global workspace registry (``$XDG_CONFIG_HOME/httk/workspaces.json``)
    from leaking between tests or into the developer's real configuration.
    """

    monkeypatch.setenv("HTTK_CONFIG_HOME", str(tmp_path_factory.mktemp("httk-config")))
    monkeypatch.setenv("HTTK_DATA_HOME", str(tmp_path_factory.mktemp("httk-store")))


@pytest.fixture(autouse=True)
def _isolated_workflow_logging() -> Iterator[None]:
    """Restore the ``httk.workflow`` logger after every test."""

    logger = logging.getLogger("httk.workflow")
    propagate = logger.propagate
    handlers = list(logger.handlers)
    level = logger.level
    yield
    logger.propagate = propagate
    logger.handlers[:] = handlers
    logger.setLevel(level)


def register_ws(path: object, name: str = "ws") -> str:
    """Register *path* under *name* and return the name."""

    register_workspace(name, str(path))
    return name


@dataclass(frozen=True)
class TestProfile:
    """Select normal or full-depth values without duplicating a test body."""

    name: str

    @property
    def extended(self) -> bool:
        return self.name == "extended"


@pytest.fixture(scope="session", autouse=True)
def test_profile() -> TestProfile:
    """The test-depth knob shared by every profiled test.

    Normal is deliberately the default for direct pytest invocations. Extended
    runs set ``HTTK_TEST_PROFILE=extended``.
    """

    name = os.environ.get("HTTK_TEST_PROFILE", "normal")
    if name not in {"normal", "extended"}:
        raise pytest.UsageError("HTTK_TEST_PROFILE must be 'normal' or 'extended'")
    return TestProfile(name)
