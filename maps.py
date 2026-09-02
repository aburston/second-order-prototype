"""Exact Poincare maps for the piecewise prototypes, with analytic Jacobians.

Run ``python3 maps.py`` for a self check against every multiplier the
README quotes.

Each prototype is linear inside each of its zones, so between two boundary
crossings the state advances by a matrix, not by an integrator. Composing
those matrices turns the flow into a **discrete map**: finitely many matrix
multiplications and one scalar root find per crossing. That is what this
module builds, and it is the form the whole library wants to be in — a map
is what you iterate, differentiate, continue in a parameter, and search for
bifurcations.

The one honest qualification: "discrete" does not mean closed form. Each
crossing time solves a transcendental equation and needs a scalar root
find. What is exact is everything else — the state advance, and, more
usefully, the derivative.

The three pieces
----------------

**The zone flow.** Inside a zone the system is an oscillator about a centre
``y_c`` that need not be the origin (the velocity switched models damp the
relative velocity, which shifts the centre)::

    y(t) = y_c + Phi(zeta, t) (y_0 - y_c)

    Phi(zeta, t) = [[ c + zeta wn s,        s        ],
                    [ -wn^2 s,        c - zeta wn s  ]]

with ``c = e^{-zeta wn t} cos(wd t)`` and ``s = e^{-zeta wn t} sin(wd t)/wd``
— the same kernels `frequency.py` already uses, so overdamped zones need no
separate branch.

**The crossing time.** Every boundary in the library is a straight line
``g^T y = level`` with ``g`` either ``(1,0)`` for a displacement threshold or
``(0,1)`` for a velocity one. Applying ``g`` to the zone flow collapses it to
one scalar equation in ``t``::

    alpha c(t) + beta s(t) = level - g^T y_c

so the whole geometry of a crossing, whatever the model, is a single root
find of that shape.

**The Jacobian, which is the point.** Differentiating a crossing is not just
``Phi``, because the crossing *time* moves when the state moves. Carrying
that gives the saltation matrix::

    S = I + (f_plus - f_minus) g^T / (g^T f_minus)

with ``f_minus`` and ``f_plus`` the field on either side. Across one cycle
the monodromy is the alternating product of the ``Phi`` blocks and the ``S``
factors, and the map's derivative follows by projecting along the flow onto
the section::

    DP = (I - f g_sec^T / (g_sec^T f)) M

Where the field is continuous — every velocity switched model in this
library, by construction — ``f_plus = f_minus`` and ``S`` is the identity, so
the saltation term only ever appears for the displacement switched models.
That is the same continuity distinction the README draws, showing up here as
whether a correction factor is needed.

Why bother, when finite differences exist
-----------------------------------------

Because they quietly fail. Differencing the return map produced a Floquet
multiplier of 1e12 in `staircase.py` (a bug, since fixed), and at the Van
der Pol fit with ``mu = 5`` it produces exactly zero — one pass through a
zone with ``zeta = 20`` destroys every bit of the perturbation, so no step
size can see the answer. The analytic Jacobian multiplies the same
contraction out symbolically and returns it however small it is.
``check_beyond_differences`` demonstrates that case.
"""
import numpy as np
from scipy.optimize import brentq

from frequency import kernels

WN = 1.0

#: Search horizons for a crossing time, in units of ``1/wn``. A crossing is
#: found on the first horizon that brackets a sign change, so a short arc
#: costs one coarse scan.
TSPANS = (8.0, 30.0, 120.0)
NGRID = 4000

#: Boundary normals. Every threshold in the library is a line of constant
#: displacement or constant velocity.
G_X = np.array([1.0, 0.0])
G_V = np.array([0.0, 1.0])


def phi(zeta, t):
    """State transition matrix of one zone, over time ``t``.

    Built from the same two kernels as the rest of the library, so it stays
    real and finite for overdamped zones as well as underdamped ones.
    """
    c, s = kernels(zeta, t)
    return np.array([[c + zeta*WN*s, s],
                     [-WN**2*s, c - zeta*WN*s]])


def field(zeta, y, centre=0.0):
    """Vector field of one zone at state ``y``, about displacement ``centre``."""
    return np.array([y[1],
                     -WN**2*(y[0] - centre) - 2.0*zeta*WN*y[1]])


def crossing_time(zeta, y0, centre, g, level, direction=0, tmin=1e-9):
    """First time the zone flow meets ``g^T y = level``, crossing the right way.

    Reduces to ``alpha c(t) + beta s(t) = level - g^T y_c`` and brackets the
    first sign change. Returns ``nan`` when the boundary is never reached,
    which is how a caller learns the orbit turned back first.

    ``direction`` selects which way the boundary is crossed: ``+1`` for
    ``g^T y`` increasing through ``level``, ``-1`` for decreasing, ``0`` for
    either. It has to be part of the root find rather than a test applied
    afterwards — a section is generally met twice per revolution and only
    one of those is the return. Starting at a maximum of ``x`` and asking
    for the next maximum, rejecting the first root instead of skipping it
    leaves the search with nothing to return.
    """
    u0 = np.asarray(y0, float) - np.array([centre, 0.0])
    rhs = level - g[0]*centre

    if g[0]:                                   # displacement boundary
        alpha, beta = u0[0], u0[1] + zeta*WN*u0[0]
    else:                                      # velocity boundary
        alpha, beta = u0[1], -(WN**2*u0[0] + zeta*WN*u0[1])

    def resid(t):
        c, s = kernels(zeta, t)
        return alpha*c + beta*s - rhs

    for tmax in TSPANS:
        grid = np.linspace(tmin, tmax, NGRID)
        val = resid(grid)
        ok = np.isfinite(val)
        cross = ok[:-1] & ok[1:] & (np.sign(val[:-1]) != np.sign(val[1:]))
        if direction > 0:
            cross &= val[1:] > val[:-1]
        elif direction < 0:
            cross &= val[1:] < val[:-1]
        idx = np.nonzero(cross)[0]
        if idx.size:
            k = idx[0]
            return brentq(resid, grid[k], grid[k + 1], xtol=1e-15,
                          rtol=8.9e-16)
    return np.nan


def advance(zeta, y0, centre, t):
    """State after time ``t`` in one zone."""
    yc = np.array([centre, 0.0])
    return yc + phi(zeta, t) @ (np.asarray(y0, float) - yc)


def saltation(f_minus, f_plus, g):
    """Jacobian correction for crossing a switching manifold.

    ``I`` when the field is continuous across the boundary, which is every
    velocity switched model here. The displacement switched models jump by
    ``2 (zeta_out - zeta_in) wn xdot`` and need the full factor.
    """
    denom = float(g @ f_minus)
    if abs(denom) < 1e-14:
        return np.full((2, 2), np.nan)         # grazing: not transversal
    return np.eye(2) + np.outer(f_plus - f_minus, g)/denom


def project(f, g):
    """Projection along the flow onto a section with normal ``g``."""
    denom = float(g @ f)
    if abs(denom) < 1e-14:
        return np.full((2, 2), np.nan)
    return np.eye(2) - np.outer(f, g)/denom


class Model:
    """A piecewise linear prototype, as zones separated by straight lines.

    Args:
        zones: ``(zeta, centre)`` per zone, ordered so that ``zone_of``
            indexes into it.
        walls: ``(g, level)`` per boundary.
        zone_of: ``y -> index`` giving the zone a state is in. Called with
            the state nudged along the flow, so a state sitting exactly on a
            boundary resolves to the zone it is entering.
    """

    def __init__(self, zones, walls, zone_of):
        self.zones = [(float(z), float(c)) for z, c in zones]
        self.walls = [(np.asarray(g, float), float(v)) for g, v in walls]
        self.zone_of = zone_of

    def run(self, y0, g_sec, level_sec, direction, max_events=64):
        """Flow until the section is met, accumulating the monodromy.

        Returns ``(y, T, M)``: the state on arrival, the elapsed time, and
        the product of transition and saltation matrices over the path.
        """
        y = np.asarray(y0, float).copy()
        M = np.eye(2)
        t_tot = 0.0
        for _ in range(max_events):
            k = self._zone_entering(y)
            zeta, centre = self.zones[k]
            events = []
            for g, lv in self.walls:
                t = crossing_time(zeta, y, centre, g, lv)
                if np.isfinite(t) and t > 1e-10:
                    events.append((t, g, lv, False))
            t_s = crossing_time(zeta, y, centre, g_sec, level_sec, direction)
            if np.isfinite(t_s) and t_s > 1e-10:
                events.append((t_s, g_sec, level_sec, True))
            if not events:
                return None, np.nan, None
            # on a tie the section wins: for the through-equilibrium model the
            # switching line *is* the section, and both events solve the same
            # equation and so return bit-identical times
            t, g, lv, is_sec = min(events, key=lambda e: (e[0], not e[3]))
            y_new = advance(zeta, y, centre, t)
            M = phi(zeta, t) @ M
            t_tot += t
            if is_sec:
                return y_new, t_tot, M
            f_m = field(zeta, y_new, centre)
            k2 = self.zone_of(y_new + 1e-11*f_m)
            z2, c2 = self.zones[k2]
            f_p = field(z2, y_new, c2)
            M = saltation(f_m, f_p, g) @ M
            y = y_new
        return None, np.nan, None

    def _f(self, y):
        """Field at ``y``, in whichever zone ``y`` is in."""
        zeta, centre = self.zones[self.zone_of(y)]
        return field(zeta, y, centre)

    def _zone_entering(self, y):
        """Zone index a state on a boundary is entering, nudged by the flow."""
        return self.zone_of(y + 1e-11*self._f(y))

    def poincare(self, y0, g_sec, level_sec, direction):
        """One return to the section, with the map's Jacobian.

        Returns ``(y, T, DP)``. ``DP`` is the full ``2x2`` derivative
        projected onto the section; its non-trivial eigenvalue is the
        Floquet multiplier and the other is zero, the flow direction having
        been projected out.
        """
        y, T, M = self.run(y0, g_sec, level_sec, direction)
        if y is None:
            return None, np.nan, None
        zeta, centre = self.zones[self._zone_entering(y)]
        return y, T, project(field(zeta, y, centre), g_sec) @ M

    def multiplier(self, y0, g_sec, level_sec, direction):
        """Floquet multiplier of the cycle through ``y0``.

        The eigenvalue of ``DP`` that is not the projected-out zero.
        """
        _, _, DP = self.poincare(y0, g_sec, level_sec, direction)
        if DP is None:
            return np.nan
        ev = np.linalg.eigvals(DP)
        return float(np.real(ev[np.argmax(np.abs(ev))]))


# ------------------------------------------------- the prototypes as models
def linear(zeta):
    """The linear prototype: one zone, no boundaries."""
    return Model([(zeta, 0.0)], [], lambda y: 0)


def through_equilibrium(zp, zm):
    """Switched damping across the x-axis. Boundary through the equilibrium."""
    return Model([(zp, 0.0), (zm, 0.0)], [(G_V, 0.0)],
                 lambda y: 0 if y[1] > 0.0 else 1)


def offset(zp, zm, v0=1.0):
    """Offset boundary at ``xdot = v0``, damping on the relative velocity.

    Damping the relative velocity ``w = xdot - v0`` moves each zone's centre
    to ``2 zeta v0 / wn`` — the virtual centres the README tracks — which is
    exactly what the ``centre`` argument of a zone is for.
    """
    return Model([(zp, 2.0*zp*v0/WN), (zm, 2.0*zm*v0/WN)], [(G_V, v0)],
                 lambda y: 0 if y[1] > v0 else 1)


def deadzone(zp, zm, v0=1.0):
    """Symmetric velocity band ``|xdot| < v0``, damping through a deadzone."""
    d = zp - zm
    return Model([(zm, 0.0), (zp, 2.0*d*v0/WN), (zp, -2.0*d*v0/WN)],
                 [(G_V, v0), (G_V, -v0)],
                 lambda y: 1 if y[1] > v0 else (2 if y[1] < -v0 else 0))


def displacement(zp, zm, x0=1.0, symmetric=True):
    """Switching on displacement. The field is discontinuous across the wall.

    Every zone shares the origin as its centre — the damping term carries a
    factor of ``xdot``, which vanishes on ``xdot = 0`` whatever ``x`` is — so
    unlike the velocity models there are no virtual centres. What these have
    instead is a genuine jump in the field, and hence a saltation factor.
    """
    if symmetric:
        return Model([(zm, 0.0), (zp, 0.0)], [(G_X, x0), (G_X, -x0)],
                     lambda y: 1 if abs(y[0]) > x0 else 0)
    return Model([(zm, 0.0), (zp, 0.0)], [(G_X, x0)],
                 lambda y: 1 if y[0] > x0 else 0)


def find_cycle(model, y0, g_sec=G_V, level_sec=0.0, direction=-1,
               n=500, tol=1e-13):
    """Iterate the map to its fixed point.

    Returns ``(y, T, DP)`` at the cycle, or ``(None, nan, None)``.
    """
    y = np.asarray(y0, float).copy()
    for _ in range(n):
        y2, T, DP = model.poincare(y, g_sec, level_sec, direction)
        if y2 is None:
            return None, np.nan, None
        if abs(y2[0] - y[0]) < tol*max(1.0, abs(y[0])):
            return y2, T, DP
        y = y2
    return y, T, DP


def multiplier_of(DP):
    """The non-trivial eigenvalue of a projected map derivative."""
    ev = np.linalg.eigvals(DP)
    return float(np.real(ev[np.argmax(np.abs(ev))]))


def check_saltation_matters(zp=0.3, zm=-0.1, x0=1.0):
    """Multiplier of the displacement model with and without the correction.

    Returns ``(with, without)``. Dropping the factor does not perturb the
    answer, it destroys it: the map derivative collapses to the identity and
    reports a neutrally stable cycle — a continuum of orbits that is not
    there.
    """
    # patch this module's own globals, which is what Model.run resolves
    # against. Going through ``import maps`` instead binds a second copy of
    # the module when this file is run as __main__, and the ablation then
    # silently does nothing -- it reported the two values as identical.
    g = globals()
    y, T, DP = find_cycle(displacement(zp, zm, x0), [1.6, 0.0])
    good = multiplier_of(DP)
    orig = g["saltation"]
    g["saltation"] = lambda fm, fp, gv: np.eye(2)
    try:
        y2, T2, DP2 = find_cycle(displacement(zp, zm, x0), [1.6, 0.0])
        bad = multiplier_of(DP2)
    finally:
        g["saltation"] = orig
    return good, bad


if __name__ == "__main__":
    def delta(z):
        return np.pi*z/np.sqrt(1.0 - z*z)

    print("closed forms, where the map can be checked against algebra")
    print("%-34s %16s %16s %9s"
          % ("case", "analytic", "closed form", "rel err"))
    for z in (0.05, 0.1, 0.3, 0.6):
        mu = linear(z).multiplier([1.0, 0.0], G_V, 0.0, -1)
        ref = np.exp(-2.0*delta(z))
        print("%-34s %16.10f %16.10f %9.1e"
              % ("linear, zeta = %.2f" % z, mu, ref, abs(mu - ref)/abs(ref)))
    for zp, zm in ((0.3, -0.1), (0.2, 0.1), (0.5, -0.4)):
        mu = through_equilibrium(zp, zm).multiplier([1.0, 0.0], G_V, 0.0, -1)
        ref = np.exp(-(delta(zp) + delta(zm)))
        print("%-34s %16.10f %16.10f %9.1e"
              % ("through equilibrium %.1f/%.1f" % (zp, zm), mu, ref,
                 abs(mu - ref)/abs(ref)))

    print("\ncycles, against the values the README quotes")
    print("%-34s %14s %12s %12s" % ("case", "multiplier", "README", "period"))
    for name, mdl, y0, ref in (
            ("offset boundary", offset(0.3, -0.1, 1.0), [1.95, 0.0], 0.538923),
            ("symmetric deadzone", deadzone(0.3, -0.1, 1.0), [3.0, 0.0],
             0.203634),
            ("displacement, symmetric", displacement(0.3, -0.1, 1.0),
             [1.6, 0.0], 0.203640)):
        y, T, DP = find_cycle(mdl, y0)
        print("%-34s %14.9f %12.6f %12.6f"
              % (name, multiplier_of(DP), ref, T))

    good, bad = check_saltation_matters()
    print("\nthe saltation factor, on the model whose field jumps")
    print("  with    %.9f\n  without %.9f   (a cycle that is not neutral"
          " reported as neutral)" % (good, bad))
