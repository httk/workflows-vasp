"""Collect hook for the packaged ``vasp.general`` workflow.

The calculation leaves its outputs in the persistent workdir; with
``publish_data``, ``publish`` also puts them under ``data/<data_prefix>/``.
Since the job decides what kind of calculation it is, an output whose file the
calculation did not write is left out rather than failing the collection: a
calculation that writes no CONTCAR has no ``final_structure``.
"""

from httk.codes.vasp.collect import read_structure, read_total_energy


def _result(record, name, prefix):
    try:
        return record.result_file(name, data_prefix=prefix)
    except ValueError:
        return None


def collect(record):
    """Return the final structure and the total energy, for the files that exist.

    :param record: The collected job record.
    :return: Extracted output roles of the general calculation.
    """
    prefix = record.parameter("data_prefix", "vasp") or ""
    outputs = {}
    contcar = _result(record, "CONTCAR", prefix)
    if contcar is not None and contcar.stat().st_size:
        outputs["final_structure"] = read_structure(contcar)
    outcar = _result(record, "OUTCAR", prefix)
    if outcar is not None:
        outputs["total_energy"] = read_total_energy(outcar)
    return outputs
