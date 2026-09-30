"""Collect hook for the packaged ``vasp.static`` workflow.

The single point leaves its OUTCAR in the persistent workdir; with
transactional data, ``publish`` also puts it under ``data/<data_prefix>/``.
"""

from httk.codes.vasp.collect import job_parameter, read_total_energy, result_file


def collect(record):
    """Return the total energy of the single point.

    :param record: The collected job record.
    :return: Extracted output roles for the static workflow.
    """
    prefix = job_parameter(record, "data_prefix", "vasp")
    return {"total_energy": read_total_energy(result_file(record, "OUTCAR", data_prefix=prefix))}
