"""Collect hook for the packaged ``vasp.dielectric`` workflow.

Reads ``dielectric.json``, ``POSCAR`` and ``OUTCAR`` of the workdir (or of the published
``data/<data_prefix>/``). The typed records are those of ``atomate2.dielectric``; a missing
*httk-analyse* is an error here, not a skipped output.
"""

import importlib.util
import json
from pathlib import Path

from httk.codes.vasp.collect import read_structure, read_total_energy
from httk.core import FileRecord
from httk.core.digests import sha256_file

#: Voigt columns xx yy zz yz xz xy, taken from VASP's printed order XX YY ZZ XY YZ ZX.
_VASP_TO_VOIGT = [0, 1, 2, 4, 5, 3]


def _file(record, path: Path) -> FileRecord:
    """Describe the JSON document as a file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        media_type="application/json",
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the cell, energy, response file, permittivities, Born charges and piezoelectric tensor.

    :param record: The collected job record.
    :return: Extracted output roles; those whose data the run did not produce are omitted.
    :raises ValueError: If *httk-analyse* is not importable.
    """
    if importlib.util.find_spec("httk.analyse") is None:
        raise ValueError("collecting vasp.dielectric needs httk-analyse for its typed records; install it")
    import numpy as np
    from httk.analyse import definitions, property_records
    from httk.analyse.matsci.electronic import summarize_dielectric
    from httk.analyse.records import bound_values, records

    prefix = record.parameter("data_prefix", "vasp") or ""
    path = record.result_file("dielectric.json", data_prefix=prefix)
    document = json.loads(path.read_text(encoding="utf-8"))
    outputs = {
        "dielectric_cell": read_structure(record.result_file("POSCAR", data_prefix=prefix)),
        "total_energy": read_total_energy(record.result_file("OUTCAR", data_prefix=prefix)),
        "dielectric_response": _file(record, path),
    }
    kinds = property_records.RECORD_KINDS
    tensors = {"high_frequency": document["electronic"], "static": document.get("static")}
    for kind, tensor in tensors.items():
        if tensor is None:
            continue
        # The mean (the stored value) is the trace over three; symmetrizing only absorbs print rounding.
        summary = summarize_dielectric((np.asarray(tensor) + np.asarray(tensor).T) / 2)
        for bound, value in zip(bound_values(summary, kind=kind), records(summary, kind=kind), strict=True):
            outputs[bound.binding.definition.rsplit("/", 1)[1]] = value
    if document.get("born_charges") is not None:
        outputs["born_effective_charges"] = kinds[definitions.BORN_EFFECTIVE_CHARGES].from_value(
            document["born_charges"]
        )
    electronic, ionic = document.get("piezoelectric_electronic"), document.get("piezoelectric_ionic")
    if electronic is not None and ionic is not None:
        total = np.add(electronic, ionic)
        outputs["piezoelectric_tensor"] = kinds[definitions.PIEZOELECTRIC_TENSOR].from_value(
            total[:, _VASP_TO_VOIGT].tolist()
        )
    return outputs
