"""Each packaged workflow directory matches its matching built-in provider."""

import os
from pathlib import Path

import pytest
from httk.core.plugins.manifest import parse_plugin_manifest
from httk.workflow.packages import load_workflow_package
from httk.workflow.vasp.workflows import PROVIDERS

REPO_ROOT = Path(__file__).resolve().parent.parent

_DIRECTORY_TO_BUILTIN_ID = {
    "vasp-relax": "httk.vasp.relax",
    "vasp-relax-bash": "httk.vasp.relax-bash",
    "vasp-static": "httk.vasp.static",
    "vasp-relax-static": "httk.vasp.relax-static",
}
_BUILTIN_BY_ID = {provider.workflow_id: provider for provider in PROVIDERS}


@pytest.mark.parametrize("directory", sorted(_DIRECTORY_TO_BUILTIN_ID))
def test_package_matches_its_builtin_provider(directory: str) -> None:
    builtin = _BUILTIN_BY_ID[_DIRECTORY_TO_BUILTIN_ID[directory]]
    provider = load_workflow_package(REPO_ROOT / directory, register=False)

    assert provider.steps == builtin.steps
    assert provider.initial_step == builtin.initial_step
    assert provider.data_mode == builtin.data_mode
    assert provider.inputs == builtin.inputs
    assert provider.outputs == builtin.outputs
    assert provider.postprocess_scripts == builtin.postprocess_scripts
    assert provider.declarations["workflow"] == builtin.declarations["workflow"]


@pytest.mark.parametrize("directory", sorted(_DIRECTORY_TO_BUILTIN_ID))
def test_package_run_file_is_executable(directory: str) -> None:
    assert os.access(REPO_ROOT / directory / "run", os.X_OK)


def test_plugin_manifest_lists_exactly_the_four_directories() -> None:
    manifest = parse_plugin_manifest(REPO_ROOT)
    assert set(manifest.workflows) == set(_DIRECTORY_TO_BUILTIN_ID)
