"""Collect hook for the packaged ``vasp.relax-static-hse`` workflow.

``after_hse`` copies the HSE stage's results into this job's workdir in the
layout of ``vasp-relax-static-coarse``: the HSE coarse evaluation at the final
geometry under ``relax/``, the HSE final-mesh static at the workdir top. With
``publish_data``, the same files go under ``data/<data_prefix>/relax/`` and
``data/<data_prefix>/static/``. The PBE and HSE stages are child jobs and
collect as child runs of their own.
"""

from httk.codes.vasp.collect import read_structure, read_total_energy


def collect(record):
    """Return the HSE-relaxed structure and the total energy of the HSE final-mesh static.

    :param record: The collected job record.
    :return: Extracted output roles from the HSE stage.
    """
    prefix = record.parameter("data_prefix", "") or ""
    return {
        # The HSE coarse evaluation at the final geometry: the same structure the HSE final static ran.
        "relaxed_structure": read_structure(record.result_file("relax/CONTCAR", data_prefix=prefix)),
        # Not relax/OUTCAR: the energy is the HSE final-mesh single point's.
        "total_energy": read_total_energy(record.result_file("OUTCAR", data_prefix=prefix, published="static/OUTCAR")),
    }
