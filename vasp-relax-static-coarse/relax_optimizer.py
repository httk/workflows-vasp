"""Trust-radius BFGS in reduced (symmetry-adapted) coordinates; standard library only.

Used by relax_driver_replay.py offline and, copied verbatim, by the in-job driver of the
vasp.relax-static-coarse workflow (owen, 2026-10-17). The optimizer knows nothing about
crystals: it works on a coordinate vector y, its energy E and gradient g = dE/dy, and asks a
`step_scale(p, trust_factor)` callback how far a proposed step p may go (the caller maps p
back to atom displacements and strains and enforces the trust radius in those units).

Rules:
- Newton step p = -B^-1 g on a positive definite Hessian model B (damped BFGS update, Powell
  0.2 rule, so B stays positive definite through soft or negative curvature).
- An evaluation whose energy rises more than `rise_tol` above the current point is rejected:
  the optimizer stays at the current point, keeps the curvature it learned from the rejected
  step, and halves the trust factor. An accepted step doubles it again, up to 1.
- Which energy rise: with acceptance "energy" the recorded one; with "gradient" the rise the
  gradients predict along the step, 0.5 (g0 + g1) . s (trapezoid; exact on a quadratic), so a
  step is judged on the surface the gradients belong to even where the recorded energy is
  inconsistent with them (owen, 2026-10-18: J-0210, an isotropic 2.15 kbar mismatch between
  energy and stress rejected 6 of 11 steps that lowered forces and stress). With "gradient" a
  step whose recorded energy rises more than `rise_ceiling` is rejected all the same. Every
  tell() records both rises.
- The update is skipped when the step is shorter than 1e-12, so a repeated geometry does not
  corrupt B.
"""


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def matvec(m, v):
    return [dot(row, v) for row in m]


def solve(m, b):
    """Gaussian elimination with partial pivoting, for small dense systems."""
    n = len(b)
    a = [list(row) + [b[i]] for i, row in enumerate(m)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(a[r][c]))
        if abs(a[piv][c]) < 1e-300:
            raise ZeroDivisionError('singular Hessian model')
        a[c], a[piv] = a[piv], a[c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            if f:
                for k in range(c, n + 1):
                    a[r][k] -= f * a[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (a[r][n] - sum(a[r][k] * x[k] for k in range(r + 1, n))) / a[r][r]
    return x


class TrustBFGS:
    def __init__(self, hessian0, step_scale, rise_tol, rise_ceiling=None):
        self.B = [list(map(float, row)) for row in hessian0]
        self.step_scale = step_scale
        self.rise_tol = float(rise_tol)
        self.rise_ceiling = None if rise_ceiling is None else float(rise_ceiling)
        self.trust_factor = 1.0
        self.current = None  # (y, E, g) of the accepted point
        self.history = []  # one dict per tell(): accepted flag, energy, trust factor

    def _update(self, s, d):
        if dot(s, s) < 1e-24:
            return False
        Bs = matvec(self.B, s)
        sBs = dot(s, Bs)
        sd = dot(s, d)
        if sBs <= 0:
            return False
        theta = 1.0 if sd >= 0.2 * sBs else 0.8 * sBs / (sBs - sd)
        r = [theta * di + (1 - theta) * bi for di, bi in zip(d, Bs)]
        sr = dot(s, r)
        n = len(s)
        for i in range(n):
            for j in range(n):
                self.B[i][j] += r[i] * r[j] / sr - Bs[i] * Bs[j] / sBs
        return theta < 1.0

    def tell(self, y, energy, grad, acceptance="energy"):
        """Record one evaluation. Returns True if the point was accepted as the new current point."""
        if acceptance not in ("energy", "gradient"):
            raise ValueError(f"unknown acceptance {acceptance!r}")
        y, grad = [float(v) for v in y], [float(v) for v in grad]
        if self.current is None:
            self.current = (y, float(energy), grad)
            self.history.append(dict(accepted=True, energy=float(energy), trust_factor=self.trust_factor, damped=False))
            return True
        y0, e0, g0 = self.current
        s = [a - b for a, b in zip(y, y0)]
        damped = self._update(s, [a - b for a, b in zip(grad, g0)])
        rise = float(energy) - e0
        predicted = 0.5 * dot([a + b for a, b in zip(g0, grad)], s)
        if acceptance == "energy":
            accepted = rise <= self.rise_tol
        else:
            accepted = predicted <= self.rise_tol and (self.rise_ceiling is None or rise <= self.rise_ceiling)
        if accepted:
            self.current = (y, float(energy), grad)
            self.trust_factor = min(1.0, 2 * self.trust_factor)
        else:
            self.trust_factor *= 0.5
        self.history.append(
            dict(
                accepted=accepted,
                energy=float(energy),
                trust_factor=self.trust_factor,
                damped=damped,
                acceptance=acceptance,
                rise=rise,
                predicted_rise=predicted,
            )
        )
        return accepted

    def ask(self):
        """The next coordinates to evaluate."""
        y0, _, g0 = self.current
        p = [-v for v in solve(self.B, g0)]
        if dot(p, g0) >= 0:  # not a descent direction (should not happen with B positive definite)
            p = [-v for v in g0]
        scale = min(1.0, self.step_scale(p, self.trust_factor))
        return [a + scale * b for a, b in zip(y0, p)]
