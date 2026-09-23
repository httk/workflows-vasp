"""Collect hook for the packaged ``vasp.static`` workflow."""

from httk.workflow.vasp.collect import collect_vasp_static


def collect(record):
    """Extract the final energy from the job record.

    :param record: The collected job record.
    :return: Extracted output roles for the static workflow.
    """
    return collect_vasp_static(record)
