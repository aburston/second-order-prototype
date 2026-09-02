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


def arc_matrix(zeta, centre, t):
    """One arc as a ``3x3`` matrix on the augmented state ``[x, xdot, 1]``.

    A zone's flow is affine, ``y -> y_c + Phi (y - y_c)``, and an affine map
    becomes a linear one by carrying a constant alongside the state::

        [ y ]        [ Phi   (I - Phi) y_c ] [ y ]
        [   ]     =  [                     ] [   ]
        [ 1 ]_{k+1}  [  0            1     ] [ 1 ]_k

    So each arc is a genuine linear difference equation whose output is the
    next one's input, and a whole cycle is a matrix *product* rather than a
    composition of functions. The virtual centres that the velocity switched
    models introduce live entirely in the third column.

    The two scalars that are not linear stay outside: which zone the state
    is in, and the dwell ``t``. Everything else is exact matrix arithmetic.
    """
    p = phi(zeta, t)
    yc = np.array([centre, 0.0])
    a = np.eye(3)
    a[:2, :2] = p
    a[:2, 2] = yc - p @ yc
    return a


def saltation_matrix(f_minus, f_plus, g):
    """The saltation correction as a ``3x3`` factor, to sit in the product.

    A boundary crossing is linear in the state too, so it multiplies into
    the same chain as the arcs rather than needing separate treatment. The
    identity for a continuous field.
    """
    a = np.eye(3)
    a[:2, :2] = saltation(f_minus, f_plus, g)
    return a


def section_matrix(f, g):
    """Projection onto the section as a ``3x3`` factor."""
    a = np.eye(3)
    a[:2, :2] = project(f, g)
    return a


def drive_response(zeta, om, amp=1.0):
    """Steady state response matrix ``P`` to a unit drive of frequency ``om``.

    The particular solution of ``xddot + 2 zeta wn xdot + wn^2 x =
    amp cos(om t)`` is a fixed linear image of the drive's own state
    ``(cos(om t), sin(om t))``::

        y_p(phase) = P [cos, sin]^T

    which is what lets a forced arc stay linear.
    """
    a2 = WN**2 - om**2
    c2 = 2.0*zeta*WN*om
    den = a2**2 + c2**2
    return amp/den*np.array([[a2, c2],
                             [c2*om, -a2*om]])


def forced_arc_matrix(zeta, centre, t, amp, om):
    """One forced arc as a ``5x5`` matrix on ``[x, xdot, cos, sin, 1]``.

    A sinusoidal drive is itself a linear system — its state
    ``(cos(om t), sin(om t))`` rotates — so carrying it alongside keeps the
    arc linear rather than affine-plus-forcing::

        [ y   ]     [ Phi   P R - Phi P   (I - Phi) y_c ] [ y   ]
        [ d   ]  =  [  0         R              0       ] [ d   ]
        [ 1   ]_k+1 [  0         0              1       ] [ 1   ]_k

    with ``d = (cos, sin)`` the drive phase, ``R`` its rotation over the
    arc, and ``P`` from :func:`drive_response`. So the forced prototypes are
    linear difference equations too, on a state of five components instead
    of three, and every nonlinearity is still only the choice of zone and
    the dwell.

    This is the form that covers **chaos**. The autonomous planar
    prototypes cannot be chaotic — Poincare-Bendixson caps them at cycles —
    so a chaotic case needs the drive, and with the drive the natural
    section is stroboscopic: sample once per drive period and the recurrence
    steps from one sample to the next.
    """
    p = phi(zeta, t)
    yc = np.array([centre, 0.0])
    pr = drive_response(zeta, om, amp)
    r = np.array([[np.cos(om*t), -np.sin(om*t)],
                  [np.sin(om*t), np.cos(om*t)]])
    m = np.eye(5)
    m[:2, :2] = p
    m[:2, 2:4] = pr @ r - p @ pr
    m[:2, 4] = yc - p @ yc
    m[2:4, 2:4] = r
    return m


def forced_state(y, phase):
    """Pack a state and drive phase into the augmented vector."""
    return np.array([y[0], y[1], np.cos(phase), np.sin(phase), 1.0])


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

    def run(self, y0, g_sec, level_sec, direction, max_events=64,
            record=None):
        """Flow until the section is met, accumulating the monodromy.

        Returns ``(y, T, M)``: the state on arrival, the elapsed time, and
        the product of transition and saltation matrices over the path.

        Pass a list as ``record`` to collect the factor chain as it is
        built: ``("arc", zone, zeta, centre, t, A)`` and
        ``("salt", g, level, S)`` entries, each ``3x3`` on the augmented
        state, in the order they multiply. That is the difference equation
        the flow actually is, written out.
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
            if record is not None:
                record.append(("arc", k, zeta, centre, t,
                               arc_matrix(zeta, centre, t)))
            M = phi(zeta, t) @ M
            t_tot += t
            if is_sec:
                return y_new, t_tot, M
            f_m = field(zeta, y_new, centre)
            k2 = self.zone_of(y_new + 1e-11*f_m)
            z2, c2 = self.zones[k2]
            f_p = field(z2, y_new, c2)
            if record is not None:
                record.append(("salt", g, lv, saltation_matrix(f_m, f_p, g)))
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


def cycle_matrix(model, r, g_sec=G_V, level_sec=0.0, direction=-1):
    """One complete cycle as a single ``3x3`` matrix, from a point on the section.

    Multiplies out every arc of one full return, so the recurrence steps
    cycle to cycle rather than arc to arc::

        [y, 1]_{k+1} = C(r_k) [y, 1]_k

    On the section ``xdot = 0`` the state is ``(r, 0)``, so this collapses
    further to a scalar affine recurrence in the amplitude alone::

        r_{k+1} = a(r_k) r_k + b(r_k),   a = C[0,0],  b = C[0,2]

    Scope: this assumes the cycle closes, which is what an underdamped
    prototype or one carrying a limit cycle does. It is not defined for a
    trajectory that leaves without returning.

    Returns ``(r_next, C, a, b)``.
    """
    y0 = np.array([float(r), 0.0])
    rec = []
    y, T, M = model.run(y0, g_sec, level_sec, direction, record=rec)
    if y is None:
        return np.nan, None, np.nan, np.nan
    c = np.eye(3)
    for item in rec:
        if item[0] == "arc":                   # state chain: arcs only
            c = item[5] @ c
    return float(y[0]), c, float(c[0, 0]), float(c[0, 2])


def cycle_coefficients(model, rs, **kw):
    """``(a, b)`` of the once-per-cycle recurrence, at several amplitudes.

    Constant coefficients mean the prototype is a genuine linear difference
    equation and, by the argument in ``MAPS.md``, cannot hold an isolated
    limit cycle. Coefficients that move with ``r`` are what an isolated
    cycle is made of.
    """
    return [(r,) + cycle_matrix(model, r, **kw)[2:] for r in rs]


def forced_state_at(zeta, y0, centre, phase, amp, om, t):
    """State after time ``t`` in one forced zone, vectorised over ``t``.

    ``y(t) = y_c + Phi(t) u + P (cos(phase + om t), sin(phase + om t))``
    with ``u = y_0 - y_c - P (cos phase, sin phase)`` constant. Writing it
    this way keeps the crossing search vectorisable — the alternative,
    building a matrix per sample, is far slower.
    """
    c, s = kernels(zeta, t)
    p = drive_response(zeta, om, amp)
    d0 = np.array([np.cos(phase), np.sin(phase)])
    yc = np.array([centre, 0.0])
    u = np.asarray(y0, float) - yc - p @ d0
    ph = phase + om*np.asarray(t)
    w0 = p[0, 0]*np.cos(ph) + p[0, 1]*np.sin(ph)
    w1 = p[1, 0]*np.cos(ph) + p[1, 1]*np.sin(ph)
    x = centre + (c + zeta*WN*s)*u[0] + s*u[1] + w0
    v = -WN**2*s*u[0] + (c - zeta*WN*s)*u[1] + w1
    return x, v


def forced_crossing(zeta, y0, centre, phase, amp, om, g, level, tmax,
                    ngrid=400, tmin=1e-10):
    """First crossing of ``g^T y = level`` within ``tmax``, under forcing.

    Same idea as the unforced version, but the residual now carries the
    particular solution as well, so it is no longer two terms. It is still
    one scalar equation in ``t``, which is all the search needs.
    """
    def resid(t):
        x, v = forced_state_at(zeta, y0, centre, phase, amp, om, t)
        return (x if g[0] else v) - level

    grid = np.linspace(tmin, tmax, ngrid)
    val = resid(grid)
    ok = np.isfinite(val)
    cross = ok[:-1] & ok[1:] & (np.sign(val[:-1]) != np.sign(val[1:]))
    idx = np.nonzero(cross)[0]
    if not idx.size:
        return np.nan
    k = idx[0]
    return brentq(resid, grid[k], grid[k + 1], xtol=1e-14, rtol=8.9e-16)


def strobe_step(model, y, phase, amp, om, max_events=48):
    """Advance exactly one drive period, and return the tangent Jacobian.

    The stroboscopic map samples once per drive period, so unlike the
    autonomous case the step ends at a *fixed time* rather than on a
    section. There is therefore no projection at the end and no saltation
    for the final partial arc — only for the boundaries actually crossed on
    the way.

    The drive phase is exogenous: it is not perturbed, so the tangent space
    stays two dimensional and ``J`` is ``2x2``.

    Returns ``(y_next, phase_next, J)``.
    """
    td = 2.0*np.pi/om
    left = td
    y = np.asarray(y, float).copy()
    ph = float(phase)
    j = np.eye(2)
    for _ in range(max_events):
        k = model.zone_of(y + 1e-11*np.array([y[1], 0.0]))
        zeta, centre = model.zones[k]
        best_t, best_g, best_l = np.nan, None, None
        for g, lv in model.walls:
            t = forced_crossing(zeta, y, centre, ph, amp, om, g, lv, left)
            if np.isfinite(t) and (not np.isfinite(best_t) or t < best_t):
                best_t, best_g, best_l = t, g, lv
        if not np.isfinite(best_t):
            x, v = forced_state_at(zeta, y, centre, ph, amp, om, left)
            j = phi(zeta, left) @ j
            return np.array([float(x), float(v)]), ph + om*left, j
        x, v = forced_state_at(zeta, y, centre, ph, amp, om, best_t)
        y_new = np.array([float(x), float(v)])
        j = phi(zeta, best_t) @ j
        f_m = field(zeta, y_new, centre) + np.array([0.0, amp*np.cos(
            ph + om*best_t)])
        k2 = model.zone_of(y_new + 1e-11*f_m)
        z2, c2 = model.zones[k2]
        f_p = field(z2, y_new, c2) + np.array([0.0, amp*np.cos(
            ph + om*best_t)])
        j = saltation(f_m, f_p, best_g) @ j
        y, ph, left = y_new, ph + om*best_t, left - best_t
    return y, ph, j


def forced_lyapunov(model, y0, amp, om, n_skip=200, n=800, phase0=0.0):
    """Largest Lyapunov exponent from the exact Jacobian product.

    The stroboscopic map's Jacobian is a product of known matrices, so the
    exponent is that product's growth rate — accumulated by renormalising a
    tangent vector each step. No twin trajectory, no separation to choose,
    and no noise floor: the estimator this replaces has a floor near 0.008,
    which was large enough to flip four chaotic verdicts out of seven.

    Returns the exponent per unit time.
    """
    td = 2.0*np.pi/om
    y = np.asarray(y0, float).copy()
    ph = float(phase0)
    for _ in range(n_skip):
        y, ph, _ = strobe_step(model, y, ph, amp, om)
    v = np.array([1.0, 0.0])
    total = 0.0
    for _ in range(n):
        y, ph, j = strobe_step(model, y, ph, amp, om)
        v = j @ v
        nrm = float(np.linalg.norm(v))
        if nrm == 0.0 or not np.isfinite(nrm):
            return np.nan
        total += np.log(nrm)
        v /= nrm
    return total/(n*td)


def staircase_model(levels, edges):
    """The multi-threshold displacement prototype, as zones and walls.

    Same recipe as :func:`displacement` with more of each — which is the
    point of building the models this way. Zones are ordered innermost
    first, matching ``staircase.py``.
    """
    zones = [(float(z), 0.0) for z in levels]
    walls = [(G_X, float(e)) for e in edges] + [(G_X, -float(e))
                                                for e in edges]
    ed = np.asarray(edges, float)

    def zone_of(y):
        return int(np.searchsorted(ed, abs(y[0]), "right"))

    return Model(zones, walls, zone_of)


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


def cycles_bracketed(model, rmin, rmax, n=240, g_sec=G_V, level_sec=0.0,
                     direction=-1):
    """Every cycle in a range, by bracketing ``P(r) - r``, stable or not.

    Iterating the map finds only attractors: an unstable cycle repels, so
    forward iteration slides off it — started near the inner cycle of the
    bistable staircase, iteration converges to the origin instead. Root
    finding on the residual sees both.

    Returns a list of ``(r, multiplier, stable)``.
    """
    def residual(r):
        y, _, _ = model.poincare(np.array([r, 0.0]), g_sec, level_sec,
                                 direction)
        return np.nan if y is None else y[0] - r

    grid = np.linspace(rmin, rmax, n)
    vals = np.array([residual(r) for r in grid])
    out = []
    for k in range(len(grid) - 1):
        a, b = vals[k], vals[k + 1]
        if np.isfinite(a) and np.isfinite(b) and a*b < 0:
            r = brentq(residual, grid[k], grid[k + 1], xtol=1e-13)
            _, _, dp = model.poincare(np.array([r, 0.0]), g_sec, level_sec,
                                      direction)
            mu = multiplier_of(dp)
            out.append((r, mu, abs(mu) < 1.0))
    return out


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
