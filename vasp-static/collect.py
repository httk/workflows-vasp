"""Collect hook for the packaged ``vasp.static`` workflow.

The single point leaves its OUTCAR in the persistent workdir; with
``publish_data``, ``publish`` also puts it under ``data/<data_prefix>/``.

The forces are a product of the staged input structure, so they are collected
only when the calculation kept that geometry: one ionic configuration, no
rattle after a remedy, and the OUTCAR cell and positions equal to the job's
input POSCAR (``files/POSCAR``, or the ``poscar`` parameter) within 1e-4 Å.
The input POSCAR, not the workdir copy, is compared: a remedy (a rattle,
``scale_lattice``, ``contcar_to_poscar``) rewrites the workdir copy.
"""

import json
import logging
import math
from itertools import islice

from httk.codes.vasp import job_remedy_history_path
from httk.codes.vasp.collect import read_structure, read_total_energy
from httk.core.property_records import CartesianSiteForcesRecord

TOLERANCE = 1e-4  # Å

log = logging.getLogger(__name__)


def _dot(left, right):
    return sum(x * y for x, y in zip(left, right, strict=True))


def _cross(left, right):
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]


def _flat(rows):
    return [x for row in rows for x in row]


def kept_geometry_forces(record, outcar_path):
    """Return the OUTCAR forces when the run kept the staged geometry, else ``None``.

    :param record: The collected job record.
    :param outcar_path: The OUTCAR of the single point.
    :return: The forces as N×3 nested lists in POSCAR site order, or ``None`` (logged with the reason).
    """
    from httk.atomistic.integrations.vasp.io import OutcarFile

    def omit(reason):
        log.info("vasp.static %s: cartesian_site_forces not collected: %s", record.job_id, reason)

    try:
        with OutcarFile(outcar_path) as outcar:
            frames = list(islice(outcar.frames(), 2))
    except ValueError as exc:
        return omit(f"the OUTCAR frames cannot be read: {exc}")
    if len(frames) != 1 or frames[0].forces is None or frames[0].positions is None:
        return omit("the OUTCAR does not hold exactly one ionic configuration with forces")
    frame = frames[0]
    history = job_remedy_history_path(record.payload)
    remedied = history.is_file() and json.loads(history.read_text(encoding="utf-8")).get("events")
    # run.py rattles the POSCAR after every applied remedy when rattle_amplitude > 0.
    if remedied and float(record.parameter("rattle_amplitude", 0.0) or 0.0) > 0:
        return omit("a remedy rattled the structure")
    poscar = record.payload / str(record.parameter("poscar", "files/POSCAR") or "files/POSCAR")
    if not poscar.is_file():
        return omit(f"the staged structure {poscar} is missing")
    try:
        staged = read_structure(poscar, precision=TOLERANCE)
    except ValueError as exc:
        return omit(str(exc))
    cell = staged.lattice_vectors
    fractional = staged.fractional_site_positions
    forces = frame.forces_floats()
    positions = frame.positions_floats()
    if not len(forces) == len(positions) == len(fractional):
        return omit(f"the OUTCAR has {len(forces)} sites, the staged structure {len(fractional)}")
    if frame.cell is not None and math.dist(_flat(frame.cell_floats()), _flat(cell)) > TOLERANCE:
        return omit("the OUTCAR cell differs from the staged structure's")
    # The reciprocal rows (without 2π) turn a Cartesian difference into a fractional one.
    a, b, c = cell
    volume = _dot(a, _cross(b, c))
    reciprocal = [
        [x / volume for x in _cross(b, c)],
        [x / volume for x in _cross(c, a)],
        [x / volume for x in _cross(a, b)],
    ]
    for site, (position, staged_fractional) in enumerate(zip(positions, fractional, strict=True)):
        delta = [_dot(position, row) - f for row, f in zip(reciprocal, staged_fractional, strict=True)]
        # Wrapped to the nearest periodic image, measured in Å.
        wrapped = [sum((d - round(d)) * vector[k] for d, vector in zip(delta, cell, strict=True)) for k in range(3)]
        if math.hypot(*wrapped) > TOLERANCE:
            return omit(f"the OUTCAR position of site {site + 1} differs from the staged structure's")
    return [list(row) for row in forces]


def collect(record):
    """Return the total energy and, when the staged geometry was kept, the site forces of the single point.

    :param record: The collected job record.
    :return: Extracted output roles for the static workflow.
    """
    prefix = record.parameter("data_prefix", "vasp") or ""
    outcar = record.result_file("OUTCAR", data_prefix=prefix)
    outputs = {"total_energy": read_total_energy(outcar)}
    forces = kept_geometry_forces(record, outcar)
    if forces is not None:
        outputs["cartesian_site_forces"] = CartesianSiteForcesRecord.from_value(forces)
    return outputs
