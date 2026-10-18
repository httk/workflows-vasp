"""Collect hook for the packaged ``vasp.dos`` workflow.

``run_scf`` archives the self-consistent run under ``scf/`` in the workdir, and
the nscf stage then runs in the workdir itself, where ``analyse`` writes
``band_edges.json``. With ``publish_data``, ``publish`` puts the two stages
under ``data/<data_prefix>/scf/`` and ``data/<data_prefix>/nscf/``. The gaps are
rebuilt from ``band_edges.json``, never re-analysed.
"""

import importlib.util
import json
from pathlib import Path

from httk.codes.vasp.collect import read_total_energy
from httk.core import FileRecord
from httk.core.digests import sha256_file


def _file(record, path: Path) -> FileRecord:
    """Describe one result file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        # Only the exact name: a compressed copy is not JSON.
        media_type="application/json" if path.name == "band_edges.json" else None,
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the scf total energy, the typed DOS, Fermi energy and gaps, and the DOS, eigenvalue and edge files.

    :param record: The collected job record.
    :return: Extracted output roles from the scf and nscf stages.
    :raises ValueError: If httk-analyse, which registers the typed records, is not installed.
    """
    if importlib.util.find_spec("httk.analyse") is None:
        # Without httk-analyse the values would be stored as generic rows and never served.
        raise ValueError("collecting vasp.dos needs httk-analyse for its typed records; install it where you collect")
    from httk.analyse.integrations.vasp import dos_from_vasp
    from httk.analyse.matsci.electronic import BandEdges
    from httk.analyse.records import bound_values, records
    from httk.atomistic.integrations.vasp.io.doscar import read_doscar

    prefix = record.parameter("data_prefix", "") or ""
    doscar = record.result_file("DOSCAR", data_prefix=prefix, published="nscf/DOSCAR")
    eigenval = record.result_file("EIGENVAL", data_prefix=prefix, published="nscf/EIGENVAL")
    edges_path = record.result_file("band_edges.json", data_prefix=prefix, published="nscf/band_edges.json")
    document = json.loads(edges_path.read_text(encoding="utf-8"))
    outputs = {
        "dos": _file(record, doscar),
        "eigenvalues": _file(record, eigenval),
        "band_edges": _file(record, edges_path),
    }
    # A given density's energy belongs to the job that computed it.
    if not document["given_density"]:
        outputs["total_energy"] = read_total_energy(record.result_file("scf/OUTCAR", data_prefix=prefix))
    spin_mode = "collinear" if document["ispin"] == 2 else "nonmagnetic"
    edges = BandEdges(**{name: document[name] for name in BandEdges.__dataclass_fields__})
    for result in (dos_from_vasp(read_doscar(doscar, spin_mode=spin_mode)), edges):
        for bound, value in zip(bound_values(result), records(result), strict=True):
            outputs[bound.binding.definition.rsplit("/", 1)[1]] = value
    return outputs
