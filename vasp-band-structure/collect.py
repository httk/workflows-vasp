"""Collect hook for the packaged ``vasp.band-structure`` workflow.

``run_scf`` archives the self-consistent run under ``scf/`` in the workdir, and
the band stage then runs in the workdir itself, next to ``band_path.json``. With
``publish_data``, ``publish`` puts the two stages under
``data/<data_prefix>/scf/`` and ``data/<data_prefix>/bands/``.
"""

from pathlib import Path

from httk.codes.vasp.collect import read_structure, read_total_energy
from httk.core import FileRecord
from httk.core.digests import sha256_file


def _file(record, path: Path) -> FileRecord:
    """Describe one result file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        # Only the exact name: a compressed copy is not JSON.
        media_type="application/json" if path.name == "band_path.json" else "application/octet-stream",
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the HPKOT cell, the scf total energy, and the band eigenvalue and path files.

    :param record: The collected job record.
    :return: Extracted output roles from the scf and band stages.
    """
    prefix = record.parameter("data_prefix", "") or ""
    return {
        "band_structure_cell": read_structure(
            record.result_file("POSCAR", data_prefix=prefix, published="bands/POSCAR")
        ),
        "total_energy": read_total_energy(record.result_file("scf/OUTCAR", data_prefix=prefix)),
        "band_eigenvalues": _file(
            record, record.result_file("EIGENVAL", data_prefix=prefix, published="bands/EIGENVAL")
        ),
        "band_path": _file(
            record, record.result_file("band_path.json", data_prefix=prefix, published="bands/band_path.json")
        ),
    }
