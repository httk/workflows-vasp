"""Collect hook for the packaged ``vasp.relax-bash`` workflow.

The relaxation leaves CONTCAR and OUTCAR in the persistent workdir; with
transactional data, ``publish`` in ``run.sh`` also puts them under
``data/<data_prefix>/``.
"""

from httk.codes.vasp.collect import job_parameter, read_structure, read_total_energy, result_file


def collect(record):
    """Return the relaxed structure and the final energy of the relaxation.

    :param record: The collected job record.
    :return: Extracted output roles for the relaxation workflow.
    """
    prefix = job_parameter(record, "data_prefix", "vasp")
    return {
        "relaxed_structure": read_structure(result_file(record, "CONTCAR", data_prefix=prefix)),
        "total_energy": read_total_energy(result_file(record, "OUTCAR", data_prefix=prefix)),
    }
