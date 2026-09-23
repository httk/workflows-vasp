"""Collect hook for the packaged ``vasp.relax-bash`` workflow."""

from httk.workflow.vasp.collect import collect_vasp_relax


def collect(record):
    """Extract the relaxed structure and final energy from the job record.

    :param record: The collected job record.
    :return: Extracted output roles for the relaxation workflow.
    """
    return collect_vasp_relax(record)
