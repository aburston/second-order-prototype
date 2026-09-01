"""Switching on displacement rather than velocity.

Run ``python3 displacement.py`` for a self check.

`frequency.py` and `symmetric.py` switch the damping on **velocity**: a
single boundary at ``xdot = v0``, or a band ``|xdot| < v0``. This module is
the same pair with the boundary on **displacement**: a single boundary at
``x = x0``, or a band ``|x| < x0``. Together the four are the complete set.

    xddot + 2 zeta(x) wn xdot + wn^2 x = 0

    asymmetric   zeta = zp for x > x0,    zm otherwise
    symmetric    zeta = zp for |x| > x0,  zm inside

The symmetric one is a piecewise constant Van der Pol: negative damping
near the origin, positive damping outside, switched on displacement.

Three structural differences from the velocity pair:

*The field is discontinuous.* The velocity models could keep it continuous
because the switched term carries a factor that vanishes on the boundary.
Here the term is ``2 zeta(x) wn xdot`` and the boundary is a line of
constant ``x``, so the jump across it is ``2 (zp - zm) wn xdot``, zero only
where the boundary meets ``xdot = 0``.

*No sliding, despite that.* The boundary is a vertical line in the phase
plane, so the component of the field normal to it is ``xdot`` — which is
the same on both sides, being continuous. The two sides therefore never
point at each other, and Filippov sliding cannot occur. Every crossing is
transversal except at the two tangency points ``(+/- x0, 0)``.

*Both regions share a centre.* The damping term vanishes at ``xdot = 0``
whatever ``x`` is, so the equilibrium is the origin in every region and
there are no virtual centres to track. The arcs are all oscillations about
the same point, which makes the exact reduction simpler than the velocity
case.

What carries over unchanged: the period depends only on ``wn`` and the two
damping ratios, amplitude is exactly proportional to ``x0``, and the
existence conditions keep their shape — the asymmetric model needs positive
mean damping, the symmetric one only needs ``zp > 0``.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from frequency import kernels, TMAX_STEPS

WN = 1.0
NGRID = 24000


def field(zp, zm, x0, symmetric):
    """Right hand side of either displacement-switched model."""
    def f(t, y):
        outside = (abs(y[0]) > x0) if symmetric else (y[0] > x0)
        z = zp if outside else zm
        return [y[1], -WN**2*y[0] - 2*z*WN*y[1]]
    return f


def _arc(zeta, xi, vi, xtarget):
    """Advance one arc within a region until the displacement reaches a target.

    Every region is an oscillation about the origin, so no centre offset is
    needed. Returns the transit time and the velocity on arrival, or NaN if
    the target is never reached.
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
            t = brentq(xt, g[k-1], g[k], xtol=1e-15, rtol=8.9e-16)
            c, s = kernels(zeta, t)
            return t, vi*c - (WN**2*xi + zeta*WN*vi)*s
    return np.nan, np.nan


def phi(zp, zm, symmetric):
    """Half angle of the boundary chord, from the cycle averaged energy balance.

    For a near circular orbit of radius ``R`` the boundary ``x = x0`` is met
    where ``cos(phi) = x0 / R``, so the orbit spends ``4 phi`` of each
    revolution outside a symmetric band, or ``2 phi - sin 2 phi`` worth of
    ``xdot^2`` weight beyond a single boundary.

    The damping force is ``2 zeta(x) wn xdot``, so the power is proportional
    to ``xdot^2`` and the balance weights each region by that. It comes out
    as

        2 phi - sin 2 phi = 2 pi rho   (asymmetric)
        2 phi - sin 2 phi =   pi rho   (symmetric)

    with the same ``rho = -zm / (zp - zm)`` as the velocity models, and the
    same halving between asymmetric and symmetric. The left side runs from
    zero to ``pi`` over ``phi`` in ``(0, pi/2)``, so the asymmetric form
    needs ``rho < 1/2`` and the symmetric one only ``rho < 1`` — exactly the
    two existence conditions.

    This is the same function of angle as the velocity models' equation
    under ``alpha = pi/2 - phi``: switching on displacement instead of
    velocity measures the chord from the other axis.

    Returns:
        The angle in radians, or NaN where no limit cycle exists.
    """
    if not (zm < 0 < zp):
        return np.nan
    rho = -zm/(zp - zm)
    rhs = (1.0 if symmetric else 2.0)*np.pi*rho
    if rhs >= np.pi:
        return np.nan
    return brentq(lambda p: 2*p - np.sin(2*p) - rhs, 1e-13, np.pi/2 - 1e-13,
                  xtol=1e-15)


def amplitude(zp, zm, x0=1.0, symmetric=True):
    """Predicted cycle radius ``x0 / cos(phi)``, NaN where none exists."""
    p = phi(zp, zm, symmetric)
    return np.nan if not np.isfinite(p) else x0/np.cos(p)


def period_series(zp, zm, symmetric=True):
    """Weak damping closed form for the period.

    Same shape as the velocity models — the fractional slowing is half the
    phase weighted mean square damping — with the outside taking ``4 phi``
    of the revolution in the symmetric model, or ``2 phi`` in the
    asymmetric one.
    """
    p = phi(zp, zm, symmetric)
    if not np.isfinite(p):
        return np.nan
    out = (4 if symmetric else 2)*p
    return 2*np.pi/WN*(1 + (out*zp**2 + (2*np.pi - out)*zm**2)/(4*np.pi))


def period_exact(zp, zm, x0=1.0, symmetric=True):
    """Exact period, by matching the arcs across the boundary crossings.

    Asymmetric: one arc outside ``x = x0`` and one inside, both returning to
    the boundary, closed by a fixed point in the entry velocity.

    Symmetric: the field is odd, so a half cycle is one outer arc from
    ``(x0, v)`` back to ``x0``, then an inner arc across to ``-x0``, and the
    cycle closes when the arrival velocity is ``-v``. One root find.

    Returns:
        ``(T, (t_out, t_in))`` with dwell times for one half cycle in the
        symmetric case and one full cycle in the asymmetric case, or
        ``(nan, None)``.
    """
    R = amplitude(zp, zm, x0, symmetric)
    if not np.isfinite(R):
        return np.nan, None
    v_guess = WN*np.sqrt(max(R**2 - x0**2, 1e-12))

    if symmetric:
        def closure(v):
            t1, v1 = _arc(zp, x0, v, x0)          # outside, back to +x0
            if not np.isfinite(t1):
                return np.nan, None
            t2, v2 = _arc(zm, x0, v1, -x0)        # inside, across to -x0
            if not np.isfinite(t2):
                return np.nan, None
            return v2 + v, (t1, t2)
        for span in (1.6, 4.0, 16.0, 64.0):
            grid = np.linspace(1e-6, span*v_guess + 4*WN*x0, 260)
            prev, pv = closure(grid[0])[0], grid[0]
            for v in grid[1:]:
                cur = closure(v)[0]
                if np.isfinite(prev) and np.isfinite(cur) and prev*cur < 0:
                    vs = brentq(lambda z: closure(z)[0], pv, v, xtol=1e-14)
                    _, dw = closure(vs)
                    return 2*(dw[0] + dw[1]), dw
                prev, pv = cur, v
        return np.nan, None

    v = v_guess
    for _ in range(400):
        t1, v1 = _arc(zp, x0, v, x0)              # outside
        if not np.isfinite(t1):
            return np.nan, None
        t2, v2 = _arc(zm, x0, v1, x0)             # inside, all the way round
        if not np.isfinite(t2):
            return np.nan, None
        if abs(v2 - v) < 1e-13*max(1.0, abs(v)):
            return t1 + t2, (t1, t2)
        v = v2
    return t1 + t2, (t1, t2)


def cycle_integrated(zp, zm, x0=1.0, symmetric=True, T=800.0):
    """Radius and period of the attracting cycle by direct integration."""
    R = amplitude(zp, zm, x0, symmetric)
    r0 = 2.5*x0 if not np.isfinite(R) else R
    def ev(t, y):
        return y[1]
    ev.direction = -1
    def big(t, y):
        return np.hypot(*y) - 1e8
    big.terminal = True
    def tiny(t, y):
        return np.hypot(*y) - 1e-7
    tiny.terminal = True
    s = solve_ivp(field(zp, zm, x0, symmetric), (0, T), [r0, 0.0],
                  events=[ev, big, tiny], rtol=1e-11, atol=1e-13)
    if len(s.t_events[1]) or len(s.t_events[2]) or len(s.t_events[0]) < 4:
        return np.nan, np.nan
    return s.y_events[0][-1][0], s.t_events[0][-1] - s.t_events[0][-2]


if __name__ == "__main__":
    print("existence: asymmetric needs mean damping > 0, symmetric needs zp > 0")
    print(f"{'zp':>6}{'zm':>7}{'mean':>7}   {'asymmetric':>22}{'symmetric':>22}")
    for zp, zm in [(0.3, -0.1), (0.1, -0.3), (0.05, -0.5), (0.3, 0.1)]:
        ra, _ = cycle_integrated(zp, zm, symmetric=False)
        rs, _ = cycle_integrated(zp, zm, symmetric=True)
        fa = f"cycle r*={ra:.4f}" if np.isfinite(ra) else "none"
        fs = f"cycle r*={rs:.4f}" if np.isfinite(rs) else "none"
        print(f"{zp:>6.2f}{zm:>7.2f}{(zp+zm)/2:>7.2f}   {fa:>22}{fs:>22}")

    for sym in (False, True):
        name = "symmetric" if sym else "asymmetric"
        print(f"\n{name}: exact reduction, amplitude and series vs integration")
        print(f"{'zp':>6}{'zm':>7}{'T integ':>12}{'T exact':>12}{'R meas':>10}"
              f"{'R energy':>10}{'T series':>11}{'err %':>8}")
        cases = ([(0.3, -0.1), (0.2, -0.15), (0.5, -0.2)] if not sym else
                 [(0.3, -0.1), (0.1, -0.3), (0.5, -0.05), (0.2, -0.2)])
        for zp, zm in cases:
            r, T = cycle_integrated(zp, zm, symmetric=sym)
            Te, _ = period_exact(zp, zm, symmetric=sym)
            Rp = amplitude(zp, zm, symmetric=sym)
            Ts = period_series(zp, zm, symmetric=sym)
            print(f"{zp:>6.2f}{zm:>7.2f}{T:>12.7f}{Te:>12.7f}{r:>10.4f}"
                  f"{Rp:>10.4f}{Ts:>11.5f}{100*(Ts-T)/T:>8.2f}")

    print("\nperiod does not depend on the boundary position:")
    for x0 in [0.25, 1.0, 4.0, 16.0]:
        r, T = cycle_integrated(0.3, -0.1, x0, symmetric=True)
        print(f"  x0={x0:>6}:  r*/x0 = {r/x0:.9f}   T = {T:.9f}")
