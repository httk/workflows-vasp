"""Collect hook for the packaged ``vasp.relax-static-coarse`` workflow.

``finish`` archives the coarse evaluation at the final geometry under
``relax/`` in the workdir, and the final-mesh static is the last evaluation,
left at the workdir top. With ``publish_data``, ``publish`` puts the two under
``data/<data_prefix>/relax/`` and ``data/<data_prefix>/static/``, and the
driver's report (every evaluation and offset, with its provenance) under
``data/<data_prefix>/relax-driver/report.json``.

The job succeeds when its last final static has run, also when that static
misses the force and stress limits. Only a relaxation whose report says
``within_limits`` is harvested: an out-of-limits result stays in the workdir
and its report, and has no outputs, so a consumer of the outputs never takes a
provisional structure or energy for an accepted one.
"""

import json

from httk.codes.vasp.collect import read_structure, read_total_energy


def collect(record):
    """Return the relaxed structure and the total energy of the final-mesh static, if within the limits.

    :param record: The collected job record.
    :return: Extracted output roles from the coarse relaxation and the final static; none when the final static
        missed the limits.
    """
    prefix = record.parameter("data_prefix", "") or ""
    report = json.loads(record.result_file("relax-driver/report.json", data_prefix=prefix).read_text(encoding="utf-8"))
    if report.get("within_limits") is not True:
        return {}
    return {
        # The coarse evaluation at the final geometry: the same structure the final static ran.
        "relaxed_structure": read_structure(record.result_file("relax/CONTCAR", data_prefix=prefix)),
        # Not relax/OUTCAR: the energy is the final-mesh single point's.
        "total_energy": read_total_energy(record.result_file("OUTCAR", data_prefix=prefix, published="static/OUTCAR")),
    }
