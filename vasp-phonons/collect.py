"""Collect hook for the packaged ``vasp.phonons`` workflow.

``analyse`` writes ``phonons.json`` (every input of the analysis and its
results), ``phonopy_disp.yaml`` and ``FORCE_SETS``, and ``displace`` the phonon
unit cell ``unitcell/POSCAR`` and the relaxed structure ``relaxed/POSCAR`` into
the workdir; with ``publish_data`` they are also under ``data/<data_prefix>/``.
The records are rebuilt from the stored values, never recomputed, so collecting
twice gives the same content ids. The displacement children collect as child
runs of their own.
"""

import json
from pathlib import Path

from httk.codes.vasp.collect import read_structure
from httk.core import FileRecord
from httk.core.digests import sha256_file

_MEDIA_TYPES = {"phonons.json": "application/json", "phonopy_disp.yaml": "application/yaml"}


def _file(record, path: Path) -> FileRecord:
    """Describe one result file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        # Only the exact names: a compressed copy is not JSON or YAML.
        media_type=_MEDIA_TYPES.get(path.name),
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the phonon cell, its band structure, DOS and, when computed, thermodynamics, and the result files.

    :param record: The collected job record.
    :return: Extracted output roles.
    :raises ValueError: If httk-analyse is not installed, or ``phonons.json`` holds no phonon results.
    """
    try:
        from httk.analyse import definitions, property_records
        from httk.analyse.matsci.phonons import HarmonicThermodynamics
        from httk.analyse.records import bound_values, records
    except ImportError as exc:
        # Without httk-analyse the values would be stored as untyped rows a store must later rebuild.
        raise ValueError(
            "collecting vasp.phonons needs httk-analyse; pip install httk-analyse where you collect"
        ) from exc
    prefix = record.parameter("data_prefix", "") or ""
    path = record.result_file("phonons.json", data_prefix=prefix)
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document.get("band_structure"), dict) or not isinstance(document.get("dos"), dict):
        raise ValueError(f"{path}: no phonon results")  # noqa: TRY004 (a malformed document is a collect ValueError)
    kinds = property_records.RECORD_KINDS
    outputs = {
        "phonon_cell": read_structure(record.result_file("unitcell/POSCAR", data_prefix=prefix)),
        "phonon_band_structure": kinds[definitions.PHONON_BAND_STRUCTURE].from_value(document["band_structure"]),
        "phonon_density_of_states": kinds[definitions.PHONON_DENSITY_OF_STATES].from_value(document["dos"]),
        "phonons_result": _file(record, path),
        "phonopy_disp": _file(record, record.result_file("phonopy_disp.yaml", data_prefix=prefix)),
        "force_sets": _file(record, record.result_file("FORCE_SETS", data_prefix=prefix)),
    }
    if document["parameters"]["relax"]:
        outputs["relaxed_structure"] = read_structure(record.result_file("relaxed/POSCAR", data_prefix=prefix))
    stored = document.get("thermodynamics")
    # Omitted (imaginary modes, or none left): the roles stay unfulfilled, never null records.
    if stored is not None:
        thermodynamics = HarmonicThermodynamics(**stored)
        for bound, value in zip(bound_values(thermodynamics), records(thermodynamics), strict=True):
            outputs[bound.binding.definition.rsplit("/", 1)[1]] = value
    return outputs
