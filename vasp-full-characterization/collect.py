"""Collect hook for the packaged ``vasp.full-characterization`` workflow.

The job publishes ``characterization.json``, ``relaxed/POSCAR`` and
``standardized/POSCAR`` to its data. Every property is collected on the child
runs, each under its own package and declaration.
"""

from httk.codes.vasp.collect import read_structure
from httk.core import FileRecord
from httk.core.digests import sha256_file


def collect(record):
    """Return the relaxed and standardized structures and the characterization document.

    :param record: The collected job record.
    :return: The ``relaxed_structure``, ``standardized_structure`` and ``characterization`` roles.
    """
    path = record.result_file("characterization.json")
    return {
        "relaxed_structure": read_structure(record.result_file("relaxed/POSCAR")),
        "standardized_structure": read_structure(record.result_file("standardized/POSCAR")),
        "characterization": FileRecord(
            path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
            path.name,
            size=path.stat().st_size,
            # Only the exact name: a compressed copy is not JSON.
            media_type="application/json" if path.name == "characterization.json" else None,
            sha256=sha256_file(path),
        ),
    }
