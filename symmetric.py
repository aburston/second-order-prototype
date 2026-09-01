"""Symmetric version of the prototype: a deadzone rather than one boundary.

Run ``python3 symmetric.py`` for a self check.

The offset prototype puts the destabilising damping on one side of a single
boundary at ``xdot = v0``. This version makes the transition symmetric about
the axis: the damping ratio is ``zm`` inside a band ``|xdot| < v0`` and
``zp`` outside it, on both sides.

Keeping the field continuous means the switched term has to vanish on both
boundaries at once, so it acts through a **deadzone** rather than through a
relative velocity:

    xddot + 2 wn [ zm xdot + (zp - zm) dz(xdot) ] + wn^2 x = 0

    dz(v) = v - v0   for v >  v0
            0        for |v| <= v0
            v + v0   for v < -v0

Inside the band the system is the linear oscillator with damping ``zm``
about the origin. Outside it, it is the linear oscillator with damping
``zp`` about a virtual centre at ``+/- 2 (zp - zm) v0 / wn``.

What changes against the single-boundary version:

*The field is odd*, ``f(-z) = -f(z)``, so the equilibrium sits at the origin
and the cycle is symmetric under ``(x, xdot) -> (-x, -xdot)``. A half cycle
determines the whole thing, which is what closes the exact reduction here.

*The existence condition loses the mean.* A large orbit now spends almost
all of its time outside the band rather than half of it, so the effective
damping runs from ``zm`` at small amplitude to ``zp`` at large amplitude,
not to the mean. A limit cycle therefore exists exactly when

    zm < 0 < zp

with no condition on ``(zp + zm)/2``. Pairs like ``zp = 0.1, zm = -0.3``,
which escape in the single-boundary version because their mean damping is
negative, have a limit cycle here.

*The amplitude equation keeps its shape and halves its right hand side.*
Single boundary: ``pi - 2a - sin 2a = 2 pi rho``. Deadzone:
``pi - 2b - sin 2b = pi rho``, with the same
``rho = -zm / (zp - zm)``. Because the left side runs from ``pi`` down to
zero, the first needs ``rho < 1/2`` and the second only ``rho < 1`` — which
are precisely the two existence conditions.

What does not change: the period is still fixed by ``wn`` and the two
damping ratios alone, the amplitude is still exactly proportional to
``v0``, and the Floquet multiplier is still ``exp(2 Lambda)`` with
``Lambda`` the dwell weighted sum of the pole real parts.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

WN = 1.0


def deadzone(v, v0):
    """The deadzone nonlinearity: zero inside the band, offset outside."""
    if v > v0:
        return v - v0
    if v < -v0:
        return v + v0
    return 0.0


def field(zp, zm, v0):
    """Right hand side of the symmetric prototype, as ``f(t, y)``."""
    d = zp - zm
    def f(t, y):
        return [y[1], -WN**2*y[0] - 2*WN*(zm*y[1] + d*deadzone(y[1], v0))]
    return f


def beta(zp, zm):
    """Half angle of the band chord, from the cycle averaged energy balance.

    For a near circular orbit of radius ``R`` the band ``|xdot| < v0`` is
    crossed at ``sin(beta) = v0 / (wn R)``, so the orbit spends ``4 beta``
    of each revolution inside it. Balancing the energy over a revolution,
    with the damping force acting through the deadzone, gives

        pi - 2 beta - sin(2 beta) = pi rho,   rho = -zm / (zp - zm)

    ``rho`` covers ``(0, 1)`` exactly on the existence region ``zm < 0 < zp``,
    and the equation has no root outside it — which is the existence
    condition falling out of the formula rather than being imposed on it.

    Returns:
        The angle in radians, or NaN when no limit cycle exists.
    """
    if not (zm < 0 < zp):
        return np.nan
    rho = -zm/(zp - zm)
    return brentq(lambda b: np.pi - 2*b - np.sin(2*b) - np.pi*rho,
                  1e-13, np.pi/2 - 1e-13, xtol=1e-15)


def amplitude(zp, zm, v0=1.0):
    """Predicted cycle radius ``v0 / (wn sin beta)``. Good to about 1%.

    NaN where no cycle exists.
    """
    b = beta(zp, zm)
    return np.nan if not np.isfinite(b) else v0/(WN*np.sin(b))


def _arc(zeta, xi0, u0, utarget, ngrid=4000):
    """One arc of a linear oscillator, from ``xdot = u0`` to ``xdot = utarget``.

    Returns the transit time and the displacement at the end, both measured
    from that region's own centre. ``t = 0`` solves the condition whenever
    ``u0 == utarget``, so the root taken is always the next one.
    """
    w = WN*np.sqrt(1 - zeta**2)
    def vel(t):
        return np.exp(-zeta*WN*t)*(u0*np.cos(w*t)
                                   - ((xi0 + zeta*WN*u0)/w)*np.sin(w*t)) - utarget
    g = np.linspace(1e-9, 2.5*2*np.pi/w, ngrid)
    val = vel(g)
    idx = np.nonzero(np.sign(val) == -np.sign(val[0]))[0]
    if len(idx) == 0:
        return np.nan, np.nan
    k = idx[0]
    t = brentq(vel, g[k-1], g[k], xtol=1e-15, rtol=8.9e-16)
    xi = np.exp(-zeta*WN*t)*(xi0*np.cos(w*t)
                             + ((u0 + zeta*WN*xi0)/w)*np.sin(w*t))
    return t, xi


def period_exact(zp, zm, v0=1.0):
    """Exact period of the limit cycle, and the dwell times per half cycle.

    A half cycle is one arc outside the band followed by one arc inside it.
    Oddness closes it: the inner arc must end at minus where the outer arc
    began. That is a single scalar condition in the entry displacement, so
    no fixed point iteration is needed — one root find does it.

    Args:
        zp: damping ratio outside the band.
        zm: damping ratio inside the band; negative for a cycle to exist.
        v0: half width of the band in velocity. The period does not depend
            on it; the argument exists so that independence can be tested.

    Returns:
        Tuple ``(T, (t_out, t_in))`` with the dwell times for one half
        cycle, or ``(nan, None)`` if no cycle closes.
    """
    c = 2*(zp - zm)*v0/WN
    def closure(a):
        t1, xi1 = _arc(zp, a - c, v0, v0)
        if not np.isfinite(t1):
            return np.nan, None
        t2, xi2 = _arc(zm, c + xi1, v0, -v0)
        if not np.isfinite(t2):
            return np.nan, None
        return xi2 + a, (t1, t2)
    R = amplitude(zp, zm, v0)
    if not np.isfinite(R):
        return np.nan, None
    grid = np.linspace(-1.8*R - 2*v0, -0.05*R, 120)
    prev, pt = closure(grid[0])[0], grid[0]
    for t in grid[1:]:
        cur = closure(t)[0]
        if np.isfinite(prev) and np.isfinite(cur) and prev*cur < 0:
            a = brentq(lambda z: closure(z)[0], pt, t, xtol=1e-14)
            _, dw = closure(a)
            return 2*(dw[0] + dw[1]), dw
        prev, pt = cur, t
    return np.nan, None


def period_series(zp, zm):
    """Weak damping closed form for the period.

    Same shape as the single-boundary version — the fractional slowing is
    half the phase weighted mean square damping ratio — but with the band
    taking ``4 beta`` of the revolution and the outside taking the rest:

        T = (2 pi / wn) [ 1 + (4 b zm^2 + (2 pi - 4 b) zp^2) / (4 pi) ]

    Less accurate here than for the single boundary, because the two
    regions' circularising frames differ more: the measured phase angles
    sum to noticeably less than ``2 pi``, which this expression assumes.
    Around 1% while both ratios stay within 0.2, a few percent by 0.5.
    """
    b = beta(zp, zm)
    return 2*np.pi/WN*(1 + (4*b*zm**2 + (2*np.pi - 4*b)*zp**2)/(4*np.pi))


def _preturn(zp, zm, v0, r):
    """One return to the section ``{xdot = 0, x > 0}`` by integration."""
    def ev(t, y):
        return y[1]
    ev.direction = -1
    s = solve_ivp(field(zp, zm, v0), (0, 15), [r, 0.0], events=ev,
                  rtol=1e-12, atol=1e-14)
    i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
    if not i:
        return np.nan, np.nan
    return s.y_events[0][i[0]][0], s.t_events[0][i[0]]


def cycle_integrated(zp, zm, v0=1.0, n=300, tol=1e-13):
    """Radius and period of the limit cycle by iterating the return map."""
    r = amplitude(zp, zm, v0)
    if not np.isfinite(r):
        r = 2.5*v0                      # no cycle predicted; seed and watch
    T = np.nan
    for _ in range(n):
        rn, T = _preturn(zp, zm, v0, r)
        if not np.isfinite(rn) or rn > 1e9 or rn < 1e-9:
            return rn, np.nan
        if abs(rn - r) < tol*max(1.0, abs(r)):
            return rn, T
        r = rn
    return r, T


if __name__ == "__main__":
    print("existence: predicted iff  zeta- < 0 < zeta+   (mean plays no part)")
    for zp, zm in [(0.3, -0.1), (0.1, -0.3), (0.05, -0.5),
                   (0.3, 0.1), (-0.1, -0.3)]:
        r, T = cycle_integrated(zp, zm)
        got = ("escapes" if (not np.isfinite(r) or r > 1e6)
               else "decays" if r < 1e-6 else f"cycle r* = {r:.4f}")
        print(f"  z+={zp:+.2f} z-={zm:+.2f}  mean={(zp+zm)/2:+.3f}  ->  {got}")

    print("\nexact reduction, amplitude and series against integration:")
    print(f"{'z+':>6}{'z-':>7}{'T integrated':>14}{'T exact':>12}"
          f"{'R meas':>10}{'R energy':>10}{'T series':>11}{'err %':>8}")
    for zp, zm in [(0.3, -0.1), (0.1, -0.3), (0.5, -0.05), (0.2, -0.2)]:
        r, T = cycle_integrated(zp, zm)
        Te, _ = period_exact(zp, zm)
        print(f"{zp:>6.2f}{zm:>7.2f}{T:>14.7f}{Te:>12.7f}{r:>10.4f}"
              f"{amplitude(zp, zm):>10.4f}{period_series(zp, zm):>11.5f}"
              f"{100*(period_series(zp,zm)-T)/T:>8.2f}")

    print("\nperiod does not depend on the band half width:")
    for v0 in [0.25, 1.0, 4.0, 16.0]:
        r, T = cycle_integrated(0.3, -0.1, v0)
        print(f"  v0={v0:>6}:  r*/v0 = {r/v0:.9f}   T = {T:.9f}")

    print("\nFloquet multiplier vs exp(2 Lambda):")
    for zp, zm in [(0.3, -0.1), (0.1, -0.3), (0.2, -0.2)]:
        r, _ = cycle_integrated(zp, zm)
        h = 1e-6
        m = (_preturn(zp, zm, 1.0, r+h)[0] - _preturn(zp, zm, 1.0, r-h)[0])/(2*h)
        _, (t_out, t_in) = period_exact(zp, zm)
        lam = -WN*(zp*2*t_out + zm*2*t_in)
        print(f"  z+={zp:+.2f} z-={zm:+.2f}:  measured {m:.6f}   "
              f"exp(2L) {np.exp(2*lam):.6f}   diff {abs(m-np.exp(2*lam)):.1e}")
