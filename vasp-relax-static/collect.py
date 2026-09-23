"""Collect hook for the packaged ``vasp.relax-static`` workflow."""

from httk.workflow.vasp.collect import collect_vasp_relax_static


def collect(record):
    """Extract the relaxed structure and final static energy from the job record.

    :param record: The collected job record.
    :return: Extracted output roles from the relaxation and static stages.
    """
    return collect_vasp_relax_static(record)
