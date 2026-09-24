"""Each packaged workflow directory: loads, matches its own runner, and is complete.

Rewritten because the workflow used to compare each directory against a
built-in ``httk.workflow.vasp.workflows.PROVIDERS`` entry; those built-in
providers are gone now that the four workflows live only here. This checks
each package against itself instead: it loads, its manifest steps equal the
set of steps its own ``run`` file reports through ``--describe``, its
``declarations["workflow"]`` equals the committed ``declaration.json``, its
declared inputs/outputs/postprocess entries are well-formed, and the plugin
manifest lists exactly the four directories.
"""

import json
import os
from pathlib import Path

import pytest
from httk.core.plugins.manifest import parse_plugin_manifest
from httk.workflow.packages import load_workflow_package
from httk.workflow.scaffold import describe_runner

REPO_ROOT = Path(__file__).resolve().parent.parent
_DIRECTORIES = ("vasp-relax", "vasp-relax-bash", "vasp-static", "vasp-relax-static")


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_package_loads(directory: str) -> None:
    provider = load_workflow_package(REPO_ROOT / directory, register=False)
    assert provider.directory == REPO_ROOT / directory
    assert provider.entry == "run"


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_package_steps_equal_the_runners_own_description(directory: str) -> None:
    provider = load_workflow_package(REPO_ROOT / directory, register=False)
    described = describe_runner(REPO_ROOT / directory / "run")
    assert set(described["steps"]) == set(provider.steps)
    assert described["workflow"] == provider.workflow_id


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_package_declaration_matches_the_committed_file(directory: str) -> None:
    provider = load_workflow_package(REPO_ROOT / directory, register=False)
    declared = json.loads(
        (REPO_ROOT / directory / "declaration.json").read_text(encoding="utf-8")
    )
    assert provider.declarations["workflow"] == declared


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_package_inputs_outputs_and_postprocess_entries_are_well_formed(
    directory: str,
) -> None:
    provider = load_workflow_package(REPO_ROOT / directory, register=False)

    assert provider.inputs, "every packaged VASP workflow declares at least one input"
    for name, destination in provider.inputs.items():
        assert isinstance(name, str) and name
        assert destination is None or (isinstance(destination, str) and destination)

    assert provider.outputs, "every packaged VASP workflow declares at least one output"
    for name, metadata in provider.outputs.items():
        assert isinstance(name, str) and name
        assert isinstance(metadata.get("entry_type"), str) and metadata["entry_type"]
        assert isinstance(metadata.get("ref"), str) and metadata["ref"]
        assert isinstance(metadata.get("description"), str) and metadata["description"]

    for name, script in provider.postprocess_scripts.items():
        assert isinstance(name, str) and name
        member = script.get("file")
        assert isinstance(member, str) and member
        path = REPO_ROOT / directory / member
        assert path.is_file()
        assert os.access(path, os.X_OK)
        assert isinstance(script.get("description"), str) and script["description"]


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_package_run_file_is_executable(directory: str) -> None:
    assert os.access(REPO_ROOT / directory / "run", os.X_OK)


def test_plugin_manifest_lists_exactly_the_four_directories() -> None:
    manifest = parse_plugin_manifest(REPO_ROOT)
    assert set(manifest.workflows) == set(_DIRECTORIES)
