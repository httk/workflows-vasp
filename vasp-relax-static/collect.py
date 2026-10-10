"""Collect hook for the packaged ``vasp.relax-static`` workflow.

``promote`` archives the relaxation under ``relax/`` in the workdir, and the
static stage then runs in the workdir itself. With ``publish_data``,
``publish`` puts the two stages under ``data/<data_prefix>/relax/`` and
``data/<data_prefix>/static/``.
"""

from httk.codes.vasp.collect import read_structure, read_total_energy


def collect(record):
    """Return the relaxed structure and the total energy of the static stage.

    :param record: The collected job record.
    :return: Extracted output roles from the relaxation and static stages.
    """
    prefix = record.parameter("data_prefix", "") or ""
    return {
        "relaxed_structure": read_structure(record.result_file("relax/CONTCAR", data_prefix=prefix)),
        # Not relax/OUTCAR: the energy is the single point's, of the relaxed cell.
        "total_energy": read_total_energy(record.result_file("OUTCAR", data_prefix=prefix, published="static/OUTCAR")),
    }
