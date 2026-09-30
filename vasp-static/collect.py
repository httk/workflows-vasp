"""Collect hook for the packaged ``vasp.static`` workflow.

The single point leaves its OUTCAR in the persistent workdir; with
transactional data, ``publish`` also puts it under ``data/<data_prefix>/``.
"""

from httk.codes.vasp.collect import read_total_energy


def collect(record):
    """Return the total energy of the single point.

    :param record: The collected job record.
    :return: Extracted output roles for the static workflow.
    """
    prefix = record.parameter("data_prefix", "vasp") or ""
    return {"total_energy": read_total_energy(record.result_file("OUTCAR", data_prefix=prefix))}
