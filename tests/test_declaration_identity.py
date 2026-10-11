"""Each package's declaration is identified by a GitHub raw URL at a declaration-only tag.

The tag namespace is ``wf-decl-v<version>``, where ``<version>`` is the declaration's own
``x-httk-definition.version``. Declaration versions form one increasing sequence shared with the tags: a changed
declaration gets a version above every existing ``wf-decl-*`` tag, and that tag is created, pushed and never moved.
Unchanged declarations keep their version and ``$id``. The URL is an identity that httk never fetches; it names the
declaration document, not the workflow code.
"""

import json
import re
import subprocess
import tomllib
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
# vasp-relax-bash ships vasp-relax's declaration, so its identity names vasp-relax's directory.
_OWNER = {"vasp-relax-bash": "vasp-relax"}
_DIRECTORIES = tuple(
    sorted(tomllib.loads((REPO_ROOT / "httk_plugin.toml").read_text(encoding="utf-8"))["plugin"]["workflows"])
)
_ID_PATTERN = re.compile(
    r"^https://raw\.githubusercontent\.com/httk/workflows-vasp/refs/tags/wf-decl-v(\d+\.\d+\.\d+)"
    r"/([a-z0-9-]+)/declaration\.json$"
)


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_declaration_identity(directory: str) -> None:
    declaration = json.loads((REPO_ROOT / directory / "declaration.json").read_text(encoding="utf-8"))
    manifest = tomllib.loads((REPO_ROOT / directory / "httk_workflow.toml").read_text(encoding="utf-8"))
    declaration_id = declaration["$id"]
    assert declaration_id == manifest["workflow"]["declaration_uri"]
    match = _ID_PATTERN.fullmatch(declaration_id)
    assert match, declaration_id
    assert match.group(1) == declaration["x-httk-definition"]["version"]
    assert match.group(2) == _OWNER.get(directory, directory)
    if directory == "vasp-relax-bash":
        assert (REPO_ROOT / directory / "declaration.json").read_bytes() == (
            REPO_ROOT / "vasp-relax" / "declaration.json"
        ).read_bytes()


@pytest.mark.parametrize("directory", _DIRECTORIES)
def test_declaration_matches_its_tag(directory: str) -> None:
    declaration = json.loads((REPO_ROOT / directory / "declaration.json").read_text(encoding="utf-8"))
    tag = f"wf-decl-v{declaration['x-httk-definition']['version']}"
    try:
        found = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "-q", "--verify", f"refs/tags/{tag}"],
            capture_output=True,
            check=False,
        )
    except OSError:
        pytest.skip("git is unavailable")
    if found.returncode != 0:
        pytest.skip(f"tag {tag} does not exist yet")
    owner = _OWNER.get(directory, directory)
    shown = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "show", f"{tag}:{owner}/declaration.json"], capture_output=True, check=False
    )
    assert shown.returncode == 0, shown.stderr
    assert shown.stdout == (REPO_ROOT / directory / "declaration.json").read_bytes()
