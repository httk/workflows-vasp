"""Collect hook for the packaged ``vasp.eos`` workflow.

``fit`` writes ``eos.json`` and ``relaxed/POSCAR`` into this job's workdir (and,
with ``publish_data``, under ``data/<data_prefix>/``); this hook reads only those.
The fit is rebuilt from the stored fitted values, never refitted, so the typed
records do not change with BLAS bits. The energy RMSE stays in ``eos.json``. The
volume children collect as child runs of their own.
"""

import json
from pathlib import Path

from httk.codes.vasp.collect import read_structure
from httk.core import FileRecord
from httk.core.digests import sha256_file


def _file(record, path: Path) -> FileRecord:
    """Describe one result file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        # Only the exact name: a compressed copy is not JSON.
        media_type="application/json" if path.name == "eos.json" else "application/octet-stream",
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the relaxed structure, the four fitted equation-of-state values, and ``eos.json``.

    :param record: The collected job record.
    :return: Extracted output roles.
    :raises ValueError: If httk-analyse is not installed, or ``eos.json`` holds no fit.
    """
    try:
        from httk.analyse.matsci.eos import BirchMurnaghanFit
        from httk.analyse.records import bound_values, records
    except ImportError as exc:
        raise ValueError(
            "collecting vasp.eos needs httk-analyse for its typed records; install it where you collect"
        ) from exc
    prefix = record.parameter("data_prefix", "") or ""
    path = record.result_file("eos.json", data_prefix=prefix)
    stored = json.loads(path.read_text(encoding="utf-8")).get("fit")
    if not isinstance(stored, dict):
        raise ValueError(f"{path}: no fitted equation of state")  # noqa: TRY004 (a malformed document is a collect ValueError)
    fit = BirchMurnaghanFit(**stored)
    outputs = {
        "relaxed_structure": read_structure(record.result_file("relaxed/POSCAR", data_prefix=prefix)),
        "eos_result": _file(record, path),
    }
    for bound, value in zip(bound_values(fit), records(fit), strict=True):
        # The RMSE statistic (a derivation of total_energy) stays in eos.json.
        if bound.binding.derivation is None:
            outputs[bound.binding.definition.rsplit("/", 1)[1]] = value
    return outputs
