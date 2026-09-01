"""Classification of the switched-damping prototype over the damping plane.

Run ``python3 stability.py`` to check the rule below against direct
integration.

For the boundary through the equilibrium the field is positively
homogeneous, so the phase portrait is scale invariant and behaviour depends
only on the *direction* of the initial state. That makes a complete
classification possible by sampling directions.

Each half plane is an ordinary second order system with poles at
``wn(-zeta +/- sqrt(zeta^2 - 1))``. Once ``|zeta| >= 1`` those poles are
real, and real poles bring invariant rays ``x2 = lambda x1``. The sign
always puts at least one inside its own half plane:

``zeta <= -1``
    Both eigenvalues are positive, so an *escaping ray* lies in that half
    plane. A trajectory that reaches it leaves along it and never returns.
``zeta >= 1``
    Both are negative, so a *decaying sector* lies in that half plane. For
    ``zeta > 1`` the whole open wedge between the two rays is invariant and
    runs to the equilibrium without ever crossing the boundary.

Whether each exists gives the four cases in ``classify``. Where neither
exists — both half planes underdamped — every trajectory rotates, the
return map applies, and the sign of the mean damping decides globally.
"""
import numpy as np
from scipy.integrate import solve_ivp

WN = 1.0


def classify(zp, zm):
    """Classify the behaviour of the prototype at one pair of damping ratios.

    Args:
        zp: damping ratio where ``xdot > 0``.
        zm: damping ratio where ``xdot < 0``.

    Returns:
        ``"decays"``, ``"escapes"``, ``"neutral"`` (a continuum of closed
        orbits, on the line of zero mean damping), or ``"mixed"`` — the one
        case where the outcome depends on the initial condition, a
        separatrix dividing trajectories that decay from those that escape.
    """
    escaping = (zp <= -1) or (zm <= -1)
    decaying = (zp >= 1) or (zm >= 1)
    if escaping and decaying:
        return "mixed"
    if escaping:
        return "escapes"
    if decaying:
        return "decays"
    mean = zp + zm
    return "decays" if mean > 1e-9 else ("escapes" if mean < -1e-9 else "neutral")


def growth_rate(zp, zm, angle, T=150.0):
    """Exponential growth rate of the state norm from one initial direction.

    Integrates from the unit vector at ``angle`` and returns the slope of
    ``log|state|``, or ``+/-1`` if the trajectory escapes or decays past a
    threshold first.

    The decay threshold is kept far above the integrator's absolute
    tolerance on purpose. Setting the two equal lets a decayed trajectory
    dissolve into integrator noise and get thrown across the boundary into
    an unstable half plane, which reports a spurious escape — that bug
    produced a wrong classification of the mixed region before it was
    caught.
    """
    def f(t, y):
        return [y[1], -WN**2*y[0] - 2*(zp if y[1] > 0 else zm)*WN*y[1]]

    def big(t, y):
        return np.hypot(*y) - 1e6
    big.terminal = True

    def tiny(t, y):
        return np.hypot(*y) - 1e-6
    tiny.terminal = True

    s = solve_ivp(f, (0, T), [np.cos(angle), np.sin(angle)], events=[big, tiny],
                  rtol=1e-11, atol=1e-14, dense_output=True)
    if len(s.t_events[0]):
        return +1.0
    if len(s.t_events[1]):
        return -1.0
    ts = np.linspace(s.t[-1]/2, s.t[-1], 400)
    return np.polyfit(ts, np.log(np.hypot(*s.sol(ts)) + 1e-300), 1)[0]


def observe(zp, zm, n=32, tol=1e-4):
    """Classify by integration, sampling ``n`` initial directions."""
    r = np.array([growth_rate(zp, zm, a)
                  for a in np.linspace(0, 2*np.pi, n, endpoint=False)])
    if (r < -tol).all():
        return "decays"
    if (r > tol).all():
        return "escapes"
    if (np.abs(r) <= tol).all():
        return "neutral"
    return "mixed"


if __name__ == "__main__":
    short = {"decays": "S", "escapes": "U", "neutral": "N", "mixed": "M"}
    zs = [-2.0, -1.2, -1.0, -0.5, 0.0, 0.5, 1.0, 1.2, 2.0]
    print("S decays   U escapes   N neutral   M mixed")
    print("rule, or rule/observed where they disagree\n")
    print("rows: zeta-   cols: zeta+")
    print("      " + "".join(f"{z:>8.1f}" for z in zs))
    bad = []
    for zm in zs:
        row = ""
        for zp in zs:
            r, o = classify(zp, zm), observe(zp, zm)
            if r == o:
                row += f"{short[r]:>8}"
            else:
                bad.append((zp, zm, r, o))
                row += f"{short[r] + '/' + short[o]:>8}"
        print(f"{zm:>6.1f}{row}", flush=True)
    print(f"\n{len(zs)**2} cells, {len(bad)} disagreements")
    for zp, zm, r, o in bad:
        rates = [growth_rate(zp, zm, a)
                 for a in np.linspace(0, 2*np.pi, 32, endpoint=False)]
        print(f"   z+={zp} z-={zm}: rule {r}, observed {o}; "
              f"max|rate| = {max(abs(x) for x in rates):.1e}")
