"""Each package's committed declaration matches the schema source that defines it.

Ported from httk-workflow's ``tests/test_z_vasp_declaration_source.py`` (``git
show 15b0b36:tests/test_z_vasp_declaration_source.py``). The original compared
a built-in provider's ``declarations["workflow"]`` mapping; the built-in
providers are gone, so this compares each package directory's committed
``declaration.json`` file directly. The schema source lives in the sibling
``httk-schemas-source`` checkout (AGENTS.md: authored definitions, generated
into vendored ``httk-schemas`` copies); when that checkout or its generated
``output/`` is absent the test skips, same as the original did for its
(differently named) sibling checkout.
"""

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
_SCHEMA_ROOT = REPO_ROOT.parent / "httk-schemas-source" / "output" / "defs" / "v0.1" / "workflows"

# vasp-relax-bash is the same relaxation authored in Bash: its manifest points
# at vasp-relax's declaration_uri and ships the identical declaration.json.
_DIRECTORY_TO_SOURCE_NAME = {
    "vasp-relax": "vasp-relax",
    "vasp-relax-bash": "vasp-relax",
    "vasp-static": "vasp-static",
    "vasp-relax-static": "vasp-relax-static",
}


@pytest.mark.skipif(
    not _SCHEMA_ROOT.is_dir(),
    reason="workspace-only: sibling httk-schemas-source checkout required",
)
@pytest.mark.parametrize("directory", sorted(_DIRECTORY_TO_SOURCE_NAME))
def test_package_declaration_equals_published_schema_source(directory: str) -> None:
    declared = json.loads((REPO_ROOT / directory / "declaration.json").read_text(encoding="utf-8"))
    expected = json.loads((_SCHEMA_ROOT / f"{_DIRECTORY_TO_SOURCE_NAME[directory]}.json").read_text(encoding="utf-8"))
    assert declared == expected
