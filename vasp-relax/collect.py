"""Collect hook for the packaged ``vasp.relax`` workflow.

The relaxation leaves CONTCAR and OUTCAR in the persistent workdir; with
``publish_data``, ``publish`` also puts them under ``data/<data_prefix>/``.
"""

from httk.codes.vasp.collect import read_structure, read_total_energy


def collect(record):
    """Return the relaxed structure and the final energy of the relaxation.

    :param record: The collected job record.
    :return: Extracted output roles for the relaxation workflow.
    """
    prefix = record.parameter("data_prefix", "vasp") or ""
    return {
        "relaxed_structure": read_structure(record.result_file("CONTCAR", data_prefix=prefix)),
        "total_energy": read_total_energy(record.result_file("OUTCAR", data_prefix=prefix)),
    }
