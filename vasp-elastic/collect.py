"""Collect hook for the packaged ``vasp.elastic`` workflow.

``fit`` writes ``elastic.json`` (every fit input and the fitted values) and
``deform`` the strained structure ``relaxed/POSCAR`` into the workdir; with
``publish_data`` both are also under ``data/<data_prefix>/``. The records are
rebuilt from the stored stiffness, never re-fitted, so collecting twice gives
the same content ids. The deformation children collect as child runs of their own.
"""

import json
from pathlib import Path

from httk.codes.vasp.collect import read_structure
from httk.core import FileRecord
from httk.core.digests import sha256_file


def _file(record, path: Path) -> FileRecord:
    """Describe the result document as a file, by its path relative to the workspace root."""

    return FileRecord(
        path.resolve().relative_to(record.workspace_root.resolve()).as_posix(),
        path.name,
        size=path.stat().st_size,
        # Only the exact name: a compressed copy is not JSON.
        media_type="application/json" if path.name == "elastic.json" else None,
        sha256=sha256_file(path),
    )


def collect(record):
    """Return the strained structure, the result document, the elastic tensor and, when stable, the moduli.

    :param record: The collected job record.
    :return: Extracted output roles.
    :raises ValueError: If httk-analyse is not installed, or ``elastic.json`` holds no fit.
    """
    try:
        from httk.analyse import definitions, property_records
        from httk.analyse.matsci.elasticity import ElasticTensor
        from httk.analyse.records import bound_values, records
    except ImportError as exc:
        # Without httk-analyse the values would be stored as untyped rows a store must later rebuild.
        raise ValueError(
            "collecting vasp.elastic needs httk-analyse; pip install httk-analyse where you collect"
        ) from exc
    prefix = record.parameter("data_prefix", "") or ""
    path = record.result_file("elastic.json", data_prefix=prefix)
    fit = json.loads(path.read_text(encoding="utf-8")).get("fit")
    if not isinstance(fit, dict):
        raise ValueError(f"{path}: no fitted elastic tensor")  # noqa: TRY004 (a malformed document is a collect ValueError)
    outputs = {
        "relaxed_structure": read_structure(record.result_file("relaxed/POSCAR", data_prefix=prefix)),
        "elastic_result": _file(record, path),
    }
    tensor = ElasticTensor(fit["stiffness"])
    if fit["stable"]:
        for bound, value in zip(bound_values(tensor), records(tensor), strict=True):
            outputs[bound.binding.definition.rsplit("/", 1)[1]] = value
        return outputs
    # An unstable stiffness has no Reuss moduli: only the tensors are collected.
    kinds = property_records.RECORD_KINDS
    outputs["elastic_tensor"] = kinds[definitions.ELASTIC_TENSOR].from_value(fit["stiffness"])
    if fit.get("compliance") is not None:
        outputs["compliance_tensor"] = kinds[definitions.COMPLIANCE_TENSOR].from_value(fit["compliance"])
    return outputs
