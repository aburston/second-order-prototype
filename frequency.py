"""Period and frequency of the offset-boundary limit cycle.

Run ``python3 frequency.py`` for a self check: the exact reduction against
direct integration, the period's independence of the offset, and the error
of the closed form.

Three results, in decreasing order of exactness:

1. ``period_exact`` reduces the period to two transcendental transit-time
   equations and one matching condition. It agrees with direct integration
   to about 1e-12 and costs a few root finds rather than a long
   integration.
2. ``period_series`` is an explicit weak-damping formula, accurate to
   better than 1% while both damping ratios stay inside about 0.35.
3. The limits in the module docstring below, which are exact.

The structural result behind all of them: with ``u`` constant the system is
invariant under ``(x, v0) -> (lambda x, lambda v0)``, and that scaling does
not touch time. The period therefore cannot depend on the offset ``v0``, on
``u``, or on amplitude. It is fixed by ``wn`` and the two damping ratios
alone, and ``wn`` only sets the timescale, so

    T = F(zeta_plus, zeta_minus) / wn

with ``F`` dimensionless. Two corners of ``F`` are exact:

``mean damping -> 0+``
    The cycle grows without bound relative to ``v0``, the boundary becomes
    negligible, and each half-plane takes half a revolution:
    ``T -> pi/wd_plus + pi/wd_minus``.
``zeta_minus -> 0-``
    The cycle shrinks onto the boundary and spends the whole revolution in
    the lower region: ``T -> 2 pi / wd_minus``.

In both, ``wd = wn sqrt(1 - zeta^2)`` is the damped natural frequency of
that half-plane. Because the correction to the undamped period is
proportional to ``zeta^2``, it has the same sign whichever way ``zeta``
points: the limit cycle is always *slower* than ``wn``, however the damping
is split.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

WN = 1.0

# Arc search windows, tried shortest first. An overdamped arc has no period
# to scale by, and near the offset variant's existence boundary a single arc
# runs past t = 100, so the horizon has to be able to grow; starting short
# keeps the resolution fine for ordinary arcs.
TMAX_STEPS = (16.0/WN, 60.0/WN, 240.0/WN)
NGRID = 24000


def wd(zeta):
    """Damped natural frequency of a half-plane with this damping ratio.

    Imaginary once ``|zeta| > 1``; use ``kernels`` rather than this directly
    if the damping may be overdamped.
    """
    return WN*np.sqrt(1 - zeta**2)


def kernels(zeta, t):
    """``exp(-zeta wn t)`` times ``cos(wd t)`` and ``sin(wd t)/wd``.

    These two products are what the transit equations actually use, and
    forming them directly is what keeps the overdamped branch usable. Both
    factors are entire functions of ``wd^2``, so they stay real when ``wd``
    turns imaginary — ``cos`` becomes ``cosh``, ``sin(wd t)/wd`` becomes
    ``sinh(mu t)/mu`` — which means one expression covers underdamped and
    overdamped alike and the transit equations need no separate form.

    Computing ``cosh`` first and multiplying by the decay afterwards
    overflows for strongly overdamped arcs: ``cosh(mu t)`` reaches infinity
    while the product is perfectly finite, and ``0 * inf`` then poisons the
    search grid with NaN. Writing the product as a sum of exponentials of
    the two characteristic roots avoids forming the large intermediate at
    all.

    Args:
        zeta: damping ratio, any magnitude.
        t: time, scalar or array.

    Returns:
        Tuple of the two damped kernels, real, shaped like ``t``.
    """
    d = np.sqrt(complex(zeta**2 - 1))*WN          # roots at -zeta*wn +/- d
    lp, lm = -zeta*WN + d, -zeta*WN - d
    with np.errstate(over="ignore", invalid="ignore"):
        ep, em = np.exp(lp*t), np.exp(lm*t)
        c = np.real((ep + em)/2)
        s = np.real((ep - em)/(2*d)) if abs(d) > 1e-12 else np.real(t*ep)
    return c, s


def transit(zeta, xi0, v0=1.0):
    """Time to cross one half-plane, and where the arc leaves it.

    Inside a half-plane the system is an ordinary linear oscillator about
    that region's own centre, so ``xi = x - centre`` obeys
    ``xi'' + 2 zeta wn xi' + wn^2 xi = 0``. The arc starts on the boundary
    with ``xdot = v0`` and ends at the next time ``xdot`` returns to ``v0``.

    ``t = 0`` always solves that condition, since the arc starts on the
    boundary; the root taken is the next one after it.

    Args:
        zeta: damping ratio of this half-plane.
        xi0: entry displacement, measured from this region's centre.
        v0: boundary velocity, the value ``xdot`` takes on entry and exit.

    Returns:
        Tuple ``(t, xi)`` of transit time and exit displacement from the
        same centre, or ``(nan, nan)`` if the arc never returns.
    """
    def vel(t):
        c, s = kernels(zeta, t)
        return v0*c - (xi0 + zeta*WN*v0)*s - v0
    for tmax in TMAX_STEPS:
        grid = np.linspace(1e-7, tmax, NGRID)
        vals = vel(grid)
        ok = np.isfinite(vals)
        idx = np.nonzero(ok & (np.sign(vals) == -np.sign(vals[0])))[0]
        if len(idx):
            k = idx[0]
            t = brentq(vel, grid[k-1], grid[k], xtol=1e-15, rtol=8.9e-16)
            c, s = kernels(zeta, t)
            return t, xi0*c + (v0 + zeta*WN*xi0)*s
    return np.nan, np.nan


def period_exact(zp, zm, v0=1.0, n=600, tol=1e-13):
    """Return the limit cycle period, and the time spent in each half-plane.

    The cycle is two arcs joined on the boundary, one per half-plane, each
    an arc of a linear oscillator about that region's own centre at
    ``2 zeta v0 / wn``. Requiring each arc to start and end at
    ``xdot = v0`` fixes its transit time; requiring the two arcs to join
    closes the cycle. That leaves a one dimensional fixed point in the
    entry displacement, iterated here.

    Entry to the upper region lies at ``x < 0``: on the boundary the
    damping term vanishes, so ``xddot = -wn^2 x``, and the velocity can
    only be increasing through ``v0`` where ``x`` is negative.

    Args:
        zp: damping ratio where ``xdot > v0``.
        zm: damping ratio where ``xdot < v0``. Must be negative, with
            ``zp + zm > 0``, for a limit cycle to exist.
        v0: boundary offset. The period does not depend on it; the
            argument exists so that independence can be tested.
        n: iteration cap on the fixed point.
        tol: relative convergence tolerance.

    Returns:
        Tuple ``(T, (t_plus, t_minus))``, or ``(nan, nan)`` if no cycle.
    """
    cp, cm = 2*zp*v0/WN, 2*zm*v0/WN
    a = -2.0*v0
    tp = tm = np.nan
    for _ in range(n):
        tp, xip = transit(zp, a - cp, v0)
        if not np.isfinite(tp):
            return np.nan, np.nan
        b = cp + xip
        tm, xim = transit(zm, b - cm, v0)
        if not np.isfinite(tm):
            return np.nan, np.nan
        a_new = cm + xim
        if abs(a_new - a) < tol*max(1.0, abs(a)):
            break
        a = a_new
    return tp + tm, (tp, tm)


def alpha(zp, zm):
    """Half-angle of the boundary chord, from the cycle-averaged energy balance.

    Treating the cycle as near circular of radius ``R`` in ``(x, xdot/wn)``,
    the boundary ``xdot = v0`` cuts it at ``sin(alpha) = v0 / (wn R)``, so
    the upper region occupies ``pi - 2 alpha`` of the revolution.

    The cycle is where the energy put in over a revolution cancels the
    energy taken out. The damping force is proportional to ``w``, not to
    ``xdot``, so the power is proportional to ``w xdot`` and the balance is
    *not* simply the time spent in each region weighted by its damping
    ratio. Doing the integral properly gives

        pi - 2 alpha - sin(2 alpha) = 2 pi rho,    rho = -zm / (zp - zm)

    and ``rho`` runs over exactly ``(0, 1/2)`` on the existence region
    ``zm < 0 < (zp + zm)/2``, mapping to ``alpha`` in ``(0, pi/2)``.

    Args:
        zp: damping ratio where ``xdot > v0``.
        zm: damping ratio where ``xdot < v0``.

    Returns:
        The angle in radians.
    """
    rho = -zm/(zp - zm)
    return brentq(lambda a: np.pi - 2*a - np.sin(2*a) - 2*np.pi*rho,
                  1e-13, np.pi/2 - 1e-13, xtol=1e-15)


def amplitude(zp, zm, v0=1.0):
    """Near-circular radius of the cycle in ``(x, xdot/wn)``.

    ``R = v0 / (wn sin(alpha))``, which makes the exact proportionality of
    amplitude to offset explicit. Good to a couple of percent while the
    damping ratios stay small; it degrades as they grow and the orbit stops
    being circular.
    """
    return v0/(WN*np.sin(alpha(zp, zm)))


def period_series(zp, zm):
    """Explicit weak-damping formula for the period.

    Each half-plane contributes its own damped period, weighted by the
    share of the revolution spent in it. Expanding
    ``1/sqrt(1 - zeta^2)`` to second order gives

        T = (2 pi / wn) [ 1 + (theta_plus zp^2 + theta_minus zm^2) / (4 pi) ]

    with ``theta_plus = pi - 2 alpha`` and ``theta_minus = pi + 2 alpha``.
    The bracket says the fractional slowing is half the phase-weighted mean
    square damping ratio. Since that is a mean of squares it is positive
    whatever the signs of the damping ratios, which is why the cycle is
    always slower than ``wn``.

    Accurate to better than 1% while both damping ratios stay within about
    0.35, degrading to a few percent by 0.5. Use ``period_exact`` when the
    number matters.
    """
    a = alpha(zp, zm)
    return 2*np.pi/WN*(1 + ((np.pi - 2*a)*zp**2 + (np.pi + 2*a)*zm**2)/(4*np.pi))


def frequency(zp, zm):
    """Cyclic frequency of the limit cycle, from the exact reduction."""
    T, _ = period_exact(zp, zm)
    return np.nan if not np.isfinite(T) else 1.0/T


def _period_integrated(zp, zm, v0=1.0, n=400, tol=1e-12):
    """Period by direct integration, used only to check ``period_exact``."""
    def field(t, y):
        w = y[1] - v0
        return [y[1], -WN**2*y[0] - 2*(zp if w > 0 else zm)*WN*w]
    xe = 2*zm*v0/WN
    def ev(t, y):
        return y[1]
    ev.direction = -1
    r, T = 2.2*v0, np.nan
    for _ in range(n):
        s = solve_ivp(field, (0, 60), [xe + r, 0.0], events=ev,
                      rtol=1e-12, atol=1e-14)
        i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
        if not i:
            return np.nan
        rn, T = s.y_events[0][i[0]][0] - xe, s.t_events[0][i[0]]
        if abs(rn - r) < tol*max(1.0, abs(r)):
            return T
        r = rn
    return T


if __name__ == "__main__":
    print("exact reduction vs direct integration:")
    for zp, zm in [(0.3, -0.1), (0.2, -0.15), (0.6, -0.3), (0.15, -0.02)]:
        Te, _ = period_exact(zp, zm)
        Ti = _period_integrated(zp, zm)
        print(f"  z+={zp:+.2f} z-={zm:+.2f}:  reduction {Te:.9f}   "
              f"integrated {Ti:.9f}   diff {abs(Te-Ti):.1e}")

    print("\nperiod does not depend on the offset:")
    for v0 in [0.25, 1.0, 4.0, 16.0]:
        T, _ = period_exact(0.3, -0.1, v0)
        print(f"  v0={v0:>6}:  T = {T:.9f}")

    print("\nclosed form against the exact reduction:")
    for zp, zm in [(0.1, -0.05), (0.2, -0.1), (0.35, -0.1), (0.6, -0.3)]:
        Te, _ = period_exact(zp, zm)
        Ts = period_series(zp, zm)
        print(f"  z+={zp:+.2f} z-={zm:+.2f}:  exact {Te:.6f}   series {Ts:.6f}"
              f"   {100*(Ts-Te)/Te:+.2f}%")
