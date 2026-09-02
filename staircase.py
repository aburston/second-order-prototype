"""Two thresholds instead of one: a staircase of damping levels.

Run ``python3 staircase.py`` for a self check.

Every prototype so far switches the damping ratio **once**, between two
values, across a single boundary or a symmetric pair of them. This module
switches it more than once. Keeping the symmetric displacement case — the
one `displacement.py` calls a piecewise constant Van der Pol — and adding a
second threshold gives five zones in the phase plane::

    |x| > b        outer      zeta_2
    a < |x| < b    middle     zeta_1
    |x| < a        core       zeta_0
    a < |x| < b    middle     zeta_1     (mirror)
    |x| > b        outer      zeta_2     (mirror)

so three damping levels and two thresholds ``0 < a < b``. The code takes an
arbitrary number of levels, because the two interesting questions both need
that generality: whether a second threshold buys new behaviour, and how
closely a staircase of ``n`` levels can imitate a genuinely smooth
nonlinearity.

    xddot + 2 zeta(x) wn xdot + wn^2 x = A cos(Om t)

    zeta(x) = levels[k],  k = how many thresholds |x| exceeds

What the second threshold buys
------------------------------

**A second limit cycle, and with it bistability.** The single threshold
model has an effective damping that runs monotonically from ``zeta_0`` at
small amplitude to ``zeta_1`` at large amplitude, so it can cross zero only
once and has at most one limit cycle. Three levels let the effective damping
turn twice. With ``zeta_0 > 0``, ``zeta_1 < 0``, ``zeta_2 > 0`` the origin
attracts, a band of negative damping sits above it, and heavy damping caps
it: an unstable cycle nested inside a stable one, with the origin and the
outer cycle both attracting and the inner cycle the boundary between their
basins.

That is hard excitation — a system that sits quietly until something knocks
it past a threshold and then runs away to a large oscillation, and stays
there. It cannot be represented by any of the four earlier prototypes, and
it is a common failure mode in the field.

The effective damping, which is where the analysis lives
-------------------------------------------------------

For a near circular orbit of radius ``R`` the damping force does work at a
rate proportional to ``zeta(x) xdot^2``, so the cycle averaged damping
weights each zone by its share of ``xdot^2``. Writing ``x = R cos(theta)``,
the zone ``|x| > e`` occupies the angles where ``|cos theta| > e/R``, and
its share of the weight is

    w(e/R) = (2 phi - sin 2 phi) / pi,    phi = arccos(e/R)

which is the same ``2 phi - sin 2 phi`` the single threshold models already
use. Stacking the levels,

    <zeta>(R) = zeta_0 + sum_k (zeta_k - zeta_{k-1}) w(e_k / R)

and a limit cycle sits wherever this crosses zero. For two levels it
reduces to ``2 phi - sin 2 phi = pi rho`` exactly as before, so the earlier
amplitude equation is the ``n = 2`` case of this one. ``check_reduces``
verifies that against `displacement.py`.

This is an averaging result, exact only as the damping goes to zero. The
exact answer comes from ``half_map``, which composes the arcs analytically,
and the two are compared in the self check rather than assumed to agree.
"""
import numpy as np
from scipy.optimize import brentq

from frequency import kernels, TMAX_STEPS

WN = 1.0
NGRID = 24000

#: The worked three level example: quiet at the origin, self exciting in a
#: band, heavily damped beyond it. Bistable.
BISTABLE_LEVELS, BISTABLE_EDGES = (0.15, -0.25, 0.40), (0.6, 1.6)


def zeta_at(x, levels, edges):
    """Damping ratio at displacement ``x``."""
    return levels[int(np.searchsorted(np.asarray(edges), abs(x), "right"))]


def field(levels, edges, amp=0.0, om=1.0):
    """Right hand side of the staircase oscillator, forced or not."""
    lv, ed = np.asarray(levels, float), np.asarray(edges, float)

    def f(t, y):
        z = lv[int(np.searchsorted(ed, abs(y[0]), "right"))]
        return [y[1], -WN**2*y[0] - 2.0*z*WN*y[1] + amp*np.cos(om*t)]
    return f


def _arc_to_x(zeta, xi, vi, xtarget):
    """Advance within one zone until the displacement reaches ``xtarget``.

    Every zone is an oscillation about the origin — the damping term
    vanishes at ``xdot = 0`` whatever ``x`` is — so no centre offset is
    needed and the kernels apply directly.

    Returns ``(t, v)`` on arrival, or ``(nan, nan)`` if the target is never
    reached, which is how the caller learns the orbit turned around first.
    """
    def xt(t):
        c, s = kernels(zeta, t)
        return xi*c + (vi + zeta*WN*xi)*s - xtarget

    for tmax in TMAX_STEPS:
        g = np.linspace(1e-9, tmax, NGRID)
        val = xt(g)
        ok = np.isfinite(val)
        idx = np.nonzero(ok & (np.sign(val) == -np.sign(val[0])))[0]
        if len(idx):
            k = idx[0]
            t = brentq(xt, g[k - 1], g[k], xtol=1e-15, rtol=8.9e-16)
            c, s = kernels(zeta, t)
            return t, vi*c - (WN**2*xi + zeta*WN*vi)*s
    return np.nan, np.nan


def _arc_to_turn(zeta, xi, vi):
    """Advance within one zone until the velocity returns to zero.

    Returns ``(t, x)`` at the turning point. Used for the last arc of a half
    cycle, and for any arc where the orbit turns before reaching the next
    threshold.
    """
    def vt(t):
        c, s = kernels(zeta, t)
        return vi*c - (WN**2*xi + zeta*WN*vi)*s

    for tmax in TMAX_STEPS:
        g = np.linspace(1e-9, tmax, NGRID)
        val = vt(g)
        ok = np.isfinite(val)
        idx = np.nonzero(ok & (np.sign(val) == -np.sign(val[0])))[0]
        if len(idx):
            k = idx[0]
            t = brentq(vt, g[k - 1], g[k], xtol=1e-15, rtol=8.9e-16)
            c, s = kernels(zeta, t)
            return t, xi*c + (vi + zeta*WN*xi)*s
    return np.nan, np.nan


def half_map(amp0, levels, edges):
    """Exact half cycle map: peak ``amp0`` to the next peak, and its duration.

    Starts at ``(amp0, 0)`` and runs to the next zero of the velocity, which
    by oddness is at ``-amp1``. Displacement decreases monotonically over
    the whole half cycle — within a zone the flow is a linear oscillation
    about the origin, whose only velocity zeros are its own extremes — so
    the thresholds are met in a fixed order and each arc is solved exactly.

    Returns ``(amp1, t_half)``, or ``(nan, nan)`` if an arc fails.
    """
    walls = sorted([-e for e in edges] + list(edges), reverse=True)
    x, v, t_tot = float(amp0), 0.0, 0.0
    for _ in range(2*len(edges) + 2):
        eps = 1e-12*max(1.0, abs(x))
        z = zeta_at(x - eps, levels, edges)
        nxt = next((w for w in walls if w < x - eps), None)
        # Both events have to be timed and the earlier one taken. Solving
        # only for the wall is wrong: the arc formula is the linear solution
        # for this zone, which goes on oscillating for ever, so its first
        # arrival at a distant wall can be several oscillations after the
        # physical trajectory turned round and left the zone. That produced
        # half cycles of 9.06 against a half period of pi, an amplitude
        # jumping from 0.59 to 2.33 across a step of 0.01 in the initial
        # peak, and a Floquet multiplier of 1e12.
        t_turn, xturn = _arc_to_turn(z, x, v)
        t_wall, vwall = (_arc_to_x(z, x, v, nxt) if nxt is not None
                         else (np.nan, np.nan))
        if np.isfinite(t_wall) and (not np.isfinite(t_turn)
                                    or t_wall < t_turn):
            x, v, t_tot = nxt, vwall, t_tot + t_wall
            continue
        if not np.isfinite(t_turn):
            return np.nan, np.nan
        return -xturn, t_tot + t_turn
    return np.nan, np.nan


def w_share(c):
    """Share of the cycle's ``xdot^2`` weight lying beyond ``|x| = c R``."""
    if c >= 1.0:
        return 0.0
    if c <= 0.0:
        return 1.0
    phi = np.arccos(c)
    return (2.0*phi - np.sin(2.0*phi))/np.pi


def mean_damping(r, levels, edges):
    """Cycle averaged damping ratio of a near circular orbit of radius ``r``.

    Zero where a limit cycle sits, positive where the amplitude decays and
    negative where it grows, so its sign is the local stability of that
    amplitude and its zero crossings are the cycles.
    """
    z = levels[0]
    for k, e in enumerate(edges):
        z += (levels[k + 1] - levels[k])*w_share(e/r)
    return z


def cycles_predicted(levels, edges, rmax=40.0, n=4000):
    """Amplitudes where the averaged damping changes sign.

    Returns a list of ``(radius, stable)``. Stability is the sign of the
    slope: the averaged damping rising through zero means larger orbits are
    damped and smaller ones are driven, which is a stable cycle.
    """
    grid = np.linspace(min(edges)*0.2, rmax, n)
    val = np.array([mean_damping(r, levels, edges) for r in grid])
    out = []
    for k in range(len(grid) - 1):
        if np.isfinite(val[k]) and np.isfinite(val[k + 1]) \
                and val[k]*val[k + 1] < 0:
            r = brentq(mean_damping, grid[k], grid[k + 1],
                       args=(levels, edges), xtol=1e-13)
            out.append((r, val[k + 1] > val[k]))
    return out


def cycles_exact(levels, edges, rmax=40.0, n=600):
    """Amplitudes fixed by the exact half cycle map, with their multipliers.

    A symmetric cycle is a fixed point of the half cycle map, because the
    field is odd. Returns a list of ``(radius, multiplier, stable)``, the
    multiplier being that of the **full** cycle, ``(dH/dA)^2``.
    """
    grid = np.linspace(min(edges)*0.2, rmax, n)
    res = np.array([half_map(a, levels, edges)[0] - a for a in grid])
    out = []
    for k in range(len(grid) - 1):
        if np.isfinite(res[k]) and np.isfinite(res[k + 1]) \
                and res[k]*res[k + 1] < 0:
            r = brentq(lambda a: half_map(a, levels, edges)[0] - a,
                       grid[k], grid[k + 1], xtol=1e-13)
            h = 1e-6*r
            d = (half_map(r + h, levels, edges)[0]
                 - half_map(r - h, levels, edges)[0])/(2.0*h)
            out.append((r, d*d, abs(d) < 1.0))
    return out


def period(r, levels, edges):
    """Exact period of the cycle at radius ``r`` — twice the half cycle."""
    return 2.0*half_map(r, levels, edges)[1]


# ----------------------------------------------------- imitating Van der Pol
def vdp_zeta(x, mu):
    """Van der Pol's damping ratio, ``-mu (1 - x^2) / 2``.

    Van der Pol is written ``xddot - mu(1 - x^2) xdot + x = 0``; matching it
    to ``xddot + 2 zeta wn xdot + wn^2 x = 0`` with ``wn = 1`` gives this.
    """
    return -0.5*mu*(1.0 - x**2)


def vdp_staircase(mu, n, xmax=3.0):
    """Staircase of ``n`` levels approximating Van der Pol out to ``xmax``.

    Thresholds are placed evenly over ``(0, xmax)`` and each level is the
    mean of Van der Pol's damping ratio across its zone, which is the choice
    that makes the staircase's averaged damping match at large radius rather
    than at a point.

    The outermost zone is unbounded and its level is therefore constant: a
    staircase saturates however many levels it has, and Van der Pol does
    not. Raising ``n`` and ``xmax`` pushes the saturation outwards without
    ever removing it, and that is the structural difference the comparison
    is about.

    Returns ``(levels, edges)``.
    """
    edges = tuple(np.linspace(0.0, xmax, n + 1)[1:-1])
    bounds = (0.0,) + edges + (xmax,)
    levels = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        xs = np.linspace(lo, hi, 400)
        levels.append(float(np.mean(vdp_zeta(xs, mu))))
    return tuple(levels), edges


def check_reduces(zp=0.4, zm=-0.15, x0=1.0):
    """Confirm the averaging here is the two level equation already in use.

    `displacement.py` solves ``2 phi - sin 2 phi = pi rho`` for the
    symmetric displacement model. That is this module's ``mean_damping``
    with two levels, so both should give the same cycle radius.
    """
    import displacement
    here = cycles_predicted((zm, zp), (x0,))
    return (here[0][0] if here else np.nan,
            displacement.amplitude(zp, zm, x0, symmetric=True))


if __name__ == "__main__":
    print("the two level case must reproduce displacement.py")
    a, b = check_reduces()
    print("  staircase %.10f   displacement.py %.10f   diff %.2e\n"
          % (a, b, abs(a - b)))

    print("three levels, bistable: levels %s edges %s"
          % (BISTABLE_LEVELS, BISTABLE_EDGES))
    print("  averaged prediction:")
    for r, st in cycles_predicted(BISTABLE_LEVELS, BISTABLE_EDGES):
        print("    r = %8.5f   %s" % (r, "stable" if st else "unstable"))
    print("  exact half cycle map:")
    for r, m, st in cycles_exact(BISTABLE_LEVELS, BISTABLE_EDGES):
        print("    r = %8.5f   multiplier %9.6f   %s   T = %.6f"
              % (r, m, "stable" if st else "unstable",
                 period(r, BISTABLE_LEVELS, BISTABLE_EDGES)))

    print("\nstaircase approximations to Van der Pol, mu = 1")
    print("  %3s %10s %10s" % ("n", "levels", "cycle r"))
    for n in (2, 3, 5, 9, 17):
        lv, ed = vdp_staircase(1.0, n)
        ex = cycles_exact(lv, ed)
        r = ex[-1][0] if ex else float("nan")
        print("  %3d %10d %10.5f" % (n, len(lv), r))
