"""Van der Pol in polar coordinates: integrating one revolution of the cycle.

Run ``python3 polar.py`` to print every number ``VANDERPOL.md`` quotes and
to write the figure it embeds into ``figures/`` as ``polar-*.png``, one light
and one dark rendering, in the style of ``figures.py``. ``python3 polar.py
checks`` prints the numbers only and ``python3 polar.py figures`` writes the
figure only. The series algebra needs ``sympy``; nothing else here does.

The question this module answers is whether one revolution of the Van der
Pol oscillator,

    xddot - mu (1 - x^2) xdot + x = 0,

can be integrated in closed form in the phase plane, from theta = 0 round
to theta = 2 pi, the way each arc of the switched prototypes can. In polar
coordinates ``x = r cos(theta)``, ``xdot = -r sin(theta)`` the equation of
the orbit is one first order equation,

    dr/dtheta = mu r (1 - r^2 cos^2 theta) sin^2 theta
                / (1 + mu (1 - r^2 cos^2 theta) sin theta cos theta),

and the answer is no: in ``s = r^2`` this is an Abel equation of the second
kind, one class beyond the last one (Riccati) that linearises, and its
solution is not an elementary or Liouvillian function. What *can* be done
exactly is the expansion in ``mu``. Every order is a polynomial in theta
times a trigonometric polynomial, so every order integrates in closed form,
and :func:`series` carries that out to any order with exact rational
arithmetic. The result is the one-revolution map as a polynomial in
``r0`` at each order in ``mu``.

Three formulations of the same revolution are kept side by side so that
each checks the others:

``revolution``
    the polar orbit equation integrated numerically in theta, carrying the
    variational equation (for the map's derivative) and the time (for the
    period) alongside;
``cartesian_return``
    the plain ``(x, xdot)`` integration to the next crossing of the
    positive ``x`` axis, which is the same map by a different route;
``series``
    the exact expansion, which must agree with both as ``mu -> 0``.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RTOL, ATOL = 1e-12, 1e-14
METHOD = "DOP853"


# ---------------------------------------------------------------- the flow
def cartesian_field(mu):
    """Right hand side of the unforced Van der Pol oscillator, ``wn = 1``."""
    def f(t, y):
        return [y[1], mu*(1.0 - y[0]**2)*y[1] - y[0]]
    return f


def polar_field(mu):
    """The orbit equation in theta, with its variational equation and time.

    State ``y = [r, dr/dr0, t]``. ``dr/dtheta`` is the ratio of the polar
    rates ``rdot = mu r (1 - r^2 cos^2) sin^2`` and
    ``thetadot = 1 + mu (1 - r^2 cos^2) sin cos``; the second component
    integrates the derivative of that ratio with respect to ``r`` along the
    orbit, which is what makes the map's slope available without a finite
    difference; the third accumulates ``dt = dtheta / thetadot``.
    """
    def f(th, y):
        r, dr, _ = y
        c, s = np.cos(th), np.sin(th)
        q = 1.0 - r*r*c*c
        num = mu*r*q*s*s
        den = 1.0 + mu*q*s*c
        dnum = mu*(q - 2.0*r*r*c*c)*s*s
        dden = -2.0*mu*r*c*c*s*c
        return [num/den, (dnum*den - num*dden)/den**2*dr, 1.0/den]
    return f


def revolution(mu, r0, rtol=RTOL, atol=ATOL):
    """Integrate one revolution, theta from 0 to 2 pi, starting at ``(r0, 0)``.

    Returns:
        ``(r, dr_dr0, T)``: the radius on return to the positive ``x`` axis,
        the derivative of that radius with respect to ``r0``, and the time
        the revolution took.
    """
    sol = solve_ivp(polar_field(mu), (0.0, 2.0*np.pi), [r0, 1.0, 0.0],
                    method=METHOD, rtol=rtol, atol=atol)
    return tuple(float(v) for v in sol.y[:, -1])


def cartesian_return(mu, r0, tmax=None):
    """The same map from the ``(x, xdot)`` integration, as a check.

    Starts at ``(r0, 0)`` and stops at the next crossing of ``xdot = 0``
    with ``xdot`` falling, which is the next visit to the positive ``x``
    axis. Returns ``(r, T)``.
    """
    if tmax is None:
        tmax = 200.0 + 4.0*mu*r0**2       # the crawl in from large r0 is slow
    f = cartesian_field(mu)

    def top(t, y):
        return y[1]
    top.direction = -1.0

    sol = solve_ivp(f, (0.0, tmax), [r0, 0.0], method=METHOD, rtol=RTOL,
                    atol=ATOL, events=top, dense_output=True)
    te = sol.t_events[0]
    te = te[te > 1e-8]
    return float(sol.sol(te[0])[0]), float(te[0])


def fixed_point(mu, lo=1.5, hi=2.5):
    """The limit cycle as the fixed point of the one-revolution map.

    Returns ``(r_star, multiplier, T)``: the amplitude (the radius on the
    positive ``x`` axis is the peak of ``x``), the Floquet multiplier
    ``P'(r_star)``, and the period.
    """
    rs = brentq(lambda r: revolution(mu, r)[0] - r, lo, hi, xtol=1e-14)
    _, dr, T = revolution(mu, rs)
    return rs, dr, T


# ------------------------------------------------- the cycle seen in theta
def cycle_polar(mu, n=20001, t_settle=None):
    """The limit cycle as ``(theta, r, thetadot)`` over one period.

    Integrated in Cartesian coordinates, so that it does not assume the
    polar chart is valid: ``theta`` is unwrapped and whether it is monotone
    is exactly what the caller wants to know.
    """
    if t_settle is None:
        t_settle = 60.0*max(1.0, mu) + 200.0
    f = cartesian_field(mu)
    warm = solve_ivp(f, (0.0, t_settle), [2.0, 0.0], method=METHOD,
                     rtol=RTOL, atol=ATOL)

    def top(t, y):
        return y[1]
    top.direction = -1.0

    sol = solve_ivp(f, (0.0, 30.0*max(1.0, mu) + 40.0), warm.y[:, -1],
                    method=METHOD, rtol=RTOL, atol=ATOL, events=top,
                    dense_output=True)
    te = sol.t_events[0]
    te = te[te > 1e-8]
    ts = np.linspace(te[0], te[1], n)
    x, v = sol.sol(ts)
    th = np.unwrap(np.arctan2(-v, x))
    th -= th[0]
    thd = 1.0 + mu*(1.0 - x**2)*np.sin(th)*np.cos(th)
    return th, np.hypot(x, v), thd


def min_thetadot(mu, r0):
    """Smallest ``thetadot`` met on one revolution started at ``(r0, 0)``.

    Positive means the polar chart describes that revolution as a graph
    ``r(theta)``; zero or negative means the angle stalls or reverses and
    the orbit equation in theta breaks down there.
    """
    f = cartesian_field(mu)

    def top(t, y):
        return y[1]
    top.direction = -1.0

    sol = solve_ivp(f, (0.0, 200.0 + 4.0*mu*r0**2), [r0, 0.0],
                    method=METHOD, rtol=1e-9, atol=1e-11, events=top,
                    dense_output=True)
    te = sol.t_events[0]
    te = te[te > 1e-8]
    x, v = sol.sol(np.linspace(0.0, te[0], 4001))
    th = np.arctan2(-v, x)
    return float((1.0 + mu*(1.0 - x**2)*np.sin(th)*np.cos(th)).min())


def chart_range(mu, r_top=8.0):
    """Interval of ``r0`` around the cycle on which the chart holds.

    Steps outward from ``r0 = 2`` until a revolution stalls, then bisects;
    inward, checks a grid and bisects if the small radii fail (they do
    once the origin is a node rather than a focus, ``mu >= 2``). Returns
    ``(r_min, r_max)`` with ``r_max = inf`` when nothing up to ``r_top``
    fails.
    """
    g = lambda r: min_thetadot(mu, r)
    hi, step = 2.0, 0.5
    while hi < r_top and g(hi + step) > 0.0:
        hi += step
    r_max = brentq(g, hi, hi + step, xtol=1e-4) if hi < r_top else np.inf
    if g(0.05) > 0.0:
        bad = [r for r in np.linspace(0.05, 2.0, 40) if g(r) <= 0.0]
        r_min = max(bad) if bad else 0.0
    else:
        r_min = brentq(g, 0.05, 2.0, xtol=1e-4)
    return r_min, r_max


# ------------------------------------------------------------- the series
def series(N=5, verbose=False):
    """The one-revolution map and period as exact series in ``mu``.

    Writes ``r(theta) = r0 + mu r_1(theta) + ... + mu^N r_N(theta)`` and
    integrates each order in closed form. Every ``r_k`` is a sum of terms
    ``theta^m exp(i n theta)`` with coefficients polynomial in ``r0`` over
    the Gaussian rationals, and that class is closed under the three
    operations needed — multiplication, integration from 0 (by parts,
    exactly), and evaluation at ``2 pi`` — so no approximation enters.

    Returns a dict with sympy expressions in ``r0`` and ``mu``:

    ``P``
        list ``[P_1, ..., P_N]``, the map ``r(2 pi) = r0 + sum mu^k P_k``;
    ``T``
        list ``[T_0, ..., T_N]``, the revolution time
        ``sum mu^k T_k(r0)``;
    ``rstar``
        the fixed point ``2 + ...`` to order ``mu^(N-1)``;
    ``mult``
        the multiplier ``P'(rstar)`` to order ``mu^(N-1)``;
    ``lnmult``
        its logarithm to the same order, which is the form that converges;
    ``Tstar``
        the period at the fixed point to order ``mu^N``.
    """
    import sympy as sp
    r0, mu = sp.symbols("r0 mu")      # no assumptions: they only slow sympy down
    I = sp.I

    def add(a, b):
        out = dict(a)
        for k, v in b.items():
            out[k] = sp.expand(out.get(k, 0) + v)
            if out[k] == 0:
                del out[k]
        return out

    def scale(a, c):
        return {k: sp.expand(v*c) for k, v in a.items()}

    def mul(a, b):
        out = {}
        for (m1, n1), v1 in a.items():
            for (m2, n2), v2 in b.items():
                k = (m1 + m2, n1 + n2)
                out[k] = out.get(k, 0) + v1*v2
        out = {k: sp.expand(v) for k, v in out.items()}
        return {k: v for k, v in out.items() if v != 0}

    def integ(a):
        """Integral from 0 to theta of ``sum c theta^m e^{i n theta}``."""
        out = {}
        for (m, n), v in a.items():
            if n == 0:
                k = (m + 1, 0)
                out[k] = out.get(k, 0) + v/sp.Integer(m + 1)
                continue
            for j in range(m + 1):
                c = (-1)**j*sp.factorial(m)/sp.factorial(m - j)/(I*n)**(j + 1)
                k = (m - j, n)
                out[k] = out.get(k, 0) + v*c
            c0 = (-1)**m*sp.factorial(m)/(I*n)**(m + 1)
            out[(0, 0)] = out.get((0, 0), 0) - v*c0
        out = {k: sp.expand(v) for k, v in out.items()}
        return {k: v for k, v in out.items() if v != 0}

    def at2pi(a):
        return sp.expand(sum(v*(2*sp.pi)**m for (m, n), v in a.items()))

    def ser_add(A, B):
        n = max(len(A), len(B))
        return [add(A[i] if i < len(A) else {}, B[i] if i < len(B) else {})
                for i in range(n)]

    def ser_mul(A, B, K):
        out = [{} for _ in range(K + 1)]
        for i, a in enumerate(A):
            for j, b in enumerate(B):
                if i + j <= K:
                    out[i + j] = add(out[i + j], mul(a, b))
        return out

    def ser_fn(A, fn):
        return [mul(a, fn) for a in A]

    ONE = {(0, 0): sp.Integer(1)}
    C = {(0, 1): sp.Rational(1, 2), (0, -1): sp.Rational(1, 2)}
    S = {(0, 1): -I/2, (0, -1): I/2}
    S2, C2 = mul(S, S), mul(C, C)
    C2S2, SC = mul(C2, S2), mul(S, C)
    C3S = mul(mul(C2, C), S)

    def f_of(r, K):
        r3 = ser_mul(ser_mul(r, r, K), r, K)
        return ser_add(ser_fn(r, S2), ser_fn(r3, scale(C2S2, -1)))

    def g_of(r, K):
        return ser_add([SC], ser_fn(ser_mul(r, r, K), scale(C3S, -1)))

    r = [{(0, 0): r0}] + [{} for _ in range(N)]
    for k in range(1, N + 1):
        f, g = f_of(r, k), g_of(r, k)
        neg_g = [scale(x, -1) for x in g]
        term, total = [ONE], {}
        for j in range(k):
            prod = ser_mul(f, term, k - 1)
            if k - 1 - j < len(prod):
                total = add(total, prod[k - 1 - j])
            term = ser_mul(term, neg_g, k - 1)
        r[k] = integ(total)
        if verbose:
            print("  order %d: %d terms" % (k, len(r[k])))
    P = [sp.factor(at2pi(r[k])) for k in range(1, N + 1)]

    g = g_of(r, N)
    gj = [[ONE]]
    for j in range(1, N + 1):
        gj.append(ser_mul(gj[-1], g, N))
    T = []
    for k in range(N + 1):
        tot = 0
        for j in range(k + 1):
            if k - j < len(gj[j]):
                tot += (-1)**j*at2pi(integ(gj[j][k - j]))
        T.append(sp.factor(sp.expand(tot)))

    Pmap = sp.expand(r0 + sum(mu**k*P[k - 1] for k in range(1, N + 1)))
    Tser = sp.expand(sum(mu**k*T[k] for k in range(N + 1)))

    def truncate(expr, order):
        expr = sp.Poly(sp.expand(expr), mu)
        return sum(expr.coeff_monomial(mu**k)*mu**k for k in range(order + 1))

    # fixed point 2 + a_1 mu + a_2 mu^2 + ..., one linear solve per order
    a = sp.symbols("a1:%d" % (N + 1))
    rstar = 2 + sum(a[i]*mu**(i + 1) for i in range(N - 1))
    eq = sp.Poly(sp.expand((Pmap - r0).subs(r0, rstar)), mu)
    solv = {}
    for i in range(N - 1):
        c = sp.expand(eq.coeff_monomial(mu**(i + 2)).subs(solv))
        solv[a[i]] = sp.expand(-c.subs(a[i], 0)/c.coeff(a[i]))
    rstar = sp.expand(rstar.subs(solv))
    dP = sp.Poly(Pmap, r0).diff(r0).as_expr()
    mult = truncate(dP.subs(r0, rstar), N - 1)
    u = sp.expand(mult - 1)                     # ln P' = ln(1 + u), u = O(mu)
    lnmult = truncate(sum((-1)**(k + 1)*u**k/k for k in range(1, N)), N - 1)
    Tstar = truncate(Tser.subs(r0, rstar), N)
    return dict(P=P, T=T, rstar=rstar, mult=mult, lnmult=lnmult, Tstar=Tstar,
                r0=r0, mu=mu, N=N)


def series_map(ser, mu, r0, order=None):
    """Evaluate the series map ``r0 -> r(2 pi)`` truncated at ``order``."""
    import sympy as sp
    order = ser["N"] if order is None else order
    expr = ser["r0"] + sum(ser["mu"]**k*ser["P"][k - 1]
                           for k in range(1, order + 1))
    return float(expr.subs({ser["r0"]: r0, ser["mu"]: mu}))


def series_period(ser, mu, r0, order=None):
    """Evaluate the series revolution time truncated at ``order``."""
    order = ser["N"] if order is None else order
    expr = sum(ser["mu"]**k*ser["T"][k] for k in range(order + 1))
    return float(expr.subs({ser["r0"]: r0, ser["mu"]: mu}))


# ------------------------------------------- an integrable model instead
def stuart_landau_map(mu, r0, R=2.0, T=2.0*np.pi):
    """The one-revolution map of the Hopf normal form, in closed form.

    ``rdot = (mu/2) r (1 - r^2/R^2)``, ``thetadot = 2 pi / T``: the radial
    law separates, so the map is ``R / sqrt(1 + (R^2/r0^2 - 1) exp(-mu T))``
    with multiplier ``exp(-mu T)`` at the fixed point ``R``. With ``R = 2``
    and ``T = 2 pi`` its map is the first order term of Van der Pol's
    series exactly.
    """
    return R/np.sqrt(1.0 + (R**2/r0**2 - 1.0)*np.exp(-mu*T))


def tuned_stuart_landau(mu):
    """Pin the normal form's three parameters to Van der Pol's cycle.

    ``R`` from the amplitude, ``T`` from the period, and ``mu_eff`` from the
    multiplier, so that the tuned map has Van der Pol's fixed point, slope
    and revolution time by construction. What it predicts elsewhere is the
    test. Returns ``(mu_eff, R, T)``.
    """
    R, m, T = fixed_point(mu)
    return -np.log(m)/T, R, T


# ------------------------------------------------------------- the checks
def checks():
    import sympy as sp
    print("theta along the limit cycle (x = r cos, xdot = -r sin)")
    print("%6s %12s %12s %12s %12s %10s" % ("mu", "T", "amplitude",
                                           "min thdot", "min r", "monotone"))
    for mu in (0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0):
        th, r, thd = cycle_polar(mu)
        T = fixed_point(mu)[2] if mu <= 5.0 else float("nan")
        print("%6.2f %12.5f %12.5f %12.5f %12.5f %10s"
              % (mu, T, r[0], thd.min(), r.min(),
                 bool(np.all(np.diff(th) > 0.0))))

    print("\nr0 range on which one revolution from (r0, 0) keeps thdot > 0")
    print("%6s %10s %10s" % ("mu", "r_min", "r_max"))
    for mu in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
        lo, hi = chart_range(mu)
        print("%6.2f %10.4f %10s"
              % (mu, lo, "> 8" if not np.isfinite(hi) else "%.4f" % hi))

    print("\npolar revolution against the cartesian return, r and T")
    for mu in (0.1, 1.0, 2.0):
        for r0 in (0.5, 2.0, 3.0):
            rp, _, Tp = revolution(mu, r0)
            rc, Tc = cartesian_return(mu, r0)
            print("  mu=%.1f r0=%.1f  r %.12f  T %.12f  |dr| %.1e |dT| %.1e"
                  % (mu, r0, rp, Tp, abs(rp - rc), abs(Tp - Tc)))

    print("\nfixed point of the revolution: amplitude, multiplier, period")
    print("%6s %16s %14s %14s %14s" % ("mu", "r*", "P'(r*)", "T",
                                       "exp(-2 pi mu)"))
    fp = {}
    for mu in (0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0):
        fp[mu] = fixed_point(mu)
        print("%6.2f %16.10f %14.6e %14.8f %14.6e"
              % (mu, fp[mu][0], fp[mu][1], fp[mu][2], np.exp(-2*np.pi*mu)))
    for mu in (3.0, 5.0):
        print("  log10 P'(r*) at mu=%.0f: %.4f, checked to six figures across"
              " DOP853 at 1e-12 and 1e-9, Radau at 1e-10" % (mu, np.log10(fp[mu][1])))

    print("\nthe exact series, order 5 (a minute or two)")
    ser = series(5, verbose=True)
    r0, mu = ser["r0"], ser["mu"]
    for k, Pk in enumerate(ser["P"], 1):
        print("  P_%d(r0) = %s" % (k, Pk))
    for k, Tk in enumerate(ser["T"]):
        print("  T_%d(r0) = %s" % (k, Tk))
    print("  r*      =", ser["rstar"])
    print("  P'(r*)  =", ser["mult"])
    print("  exp(-2 pi mu) =", sp.series(sp.exp(-2*sp.pi*mu), mu, 0, 5).removeO())
    print("  ln P'(r*) =", ser["lnmult"])
    print("  T(r*)   =", ser["Tstar"])
    print("  T(r*)/2pi =", sp.expand(ser["Tstar"]/(2*sp.pi)))
    print("  P_1 at r0 = 2:", ser["P"][0].subs(r0, 2),
          "  T_2 at r0 = 2:", ser["T"][2].subs(r0, 2),
          "  T_3 at r0 = 2:", ser["T"][3].subs(r0, 2))

    print("\nseries map against the integrated revolution, |error| in r(2 pi)")
    print("%6s %6s %12s %12s %12s %12s" % ("mu", "r0", "order 1", "order 3",
                                          "order 5", "numeric r"))
    for m in (0.1, 0.3, 0.5, 1.0):
        for rr in (1.0, 2.0, 3.0):
            rn = revolution(m, rr)[0]
            errs = [abs(series_map(ser, m, rr, k) - rn) for k in (1, 3, 5)]
            print("%6.2f %6.1f %12.2e %12.2e %12.2e %12.8f"
                  % (m, rr, errs[0], errs[1], errs[2], rn))

    print("\nfixed point, multiplier and period: series against integration")
    print("%6s %12s %12s %12s %12s %12s %12s %12s"
          % ("mu", "r* num", "r* series", "P' num", "P' series",
             "exp(ln ser)", "T num", "T series"))
    for m in (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        rs, dr, T = fp[m] if m in fp else fixed_point(m)
        print("%6.2f %12.7f %12.7f %12.4e %12.4e %12.4e %12.6f %12.6f"
              % (m, rs, float(ser["rstar"].subs(mu, m)), dr,
                 float(ser["mult"].subs(mu, m)),
                 float(sp.exp(ser["lnmult"].subs(mu, m))), T,
                 float(ser["Tstar"].subs(mu, m))))

    print("\nStuart-Landau against Van der Pol, |error| in r(2 pi), untuned")
    print("%6s %12s %12s %12s %12s" % ("mu", "r0 = 1", "r0 = 2", "r0 = 3",
                                      "r0 = 5"))
    for m in (0.1, 0.3, 0.5, 1.0):
        errs = [abs(stuart_landau_map(m, rr) - revolution(m, rr)[0])
                for rr in (1.0, 2.0, 3.0, 5.0)]
        print("%6.2f " % m + " ".join("%12.1e" % e for e in errs))

    print("\ntuned to the amplitude, period and multiplier: relative error"
          " elsewhere")
    print("%6s %8s %8s %12s %12s %12s %12s"
          % ("mu", "mu_eff", "R", "r0 = 0.5", "r0 = 1", "r0 = 3", "r0 = 5"))
    for m in (0.3, 0.5, 1.0, 2.0):
        me, R, T = tuned_stuart_landau(m)
        errs = []
        for rr in (0.5, 1.0, 3.0, 5.0):
            pv = revolution(m, rr)[0]
            errs.append(abs(stuart_landau_map(me, rr, R, T) - pv)/pv)
        print("%6.2f %8.4f %8.4f " % (m, me, R)
              + " ".join("%12.1e" % e for e in errs))
    return ser, fp


# ------------------------------------------------------------- the figure
def fig_polar(th, name, ser=None):
    """One revolution in theta, the map, and how far the series reaches."""
    import matplotlib
    matplotlib.use("Agg")
    from figures import style, legend, save, newfig
    if ser is None:
        ser = series(5)

    fig, axes = newfig(th, 1, 3, figsize=(14.5, 4.6))

    ax = axes[0]
    for i, mu in enumerate((0.1, 1.0, 5.0)):
        t, r, _ = cycle_polar(mu)
        ax.plot(np.degrees(t), r, color=th["series"][i], linewidth=1.8,
                label="$\\mu = %g$" % mu, zorder=3)
    ax.axhline(2.0, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)),
               zorder=2)
    ax.text(0.985, 2.0, "$r = 2$", transform=ax.get_yaxis_transform(),
            fontsize=8, color=th["ink2"], va="bottom", ha="right", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_ylim(1.0, 2.6)
    style(ax, th, "$\\theta$ (degrees)", "$r$",
          "The limit cycle as $r(\\theta)$ over one revolution")
    legend(ax, th, loc="lower right")

    ax = axes[1]
    mu = 0.3
    rr = np.linspace(0.2, 3.0, 141)
    num = np.array([revolution(mu, r)[0] - r for r in rr])
    ax.plot(rr, num, color=th["series"][0], linewidth=2.2,
            label="integrated", zorder=4)
    for k, ls, col in ((1, (0, (2, 2)), th["series"][1]),
                       (3, (0, (5, 3)), th["series"][2]),
                       (5, (0, (8, 3, 2, 3)), th["ink2"])):
        s = np.array([series_map(ser, mu, r, k) - r for r in rr])
        ax.plot(rr, s, color=col, linewidth=1.4, linestyle=ls,
                label="series, order %d" % k, zorder=3)
    ax.axhline(0.0, color=th["axis"], linewidth=0.8, zorder=1)
    ax.set_ylim(-1.5, 0.6)
    style(ax, th, "$r_0$", "$P(r_0) - r_0$",
          "One revolution at $\\mu = %g$" % mu)
    legend(ax, th, loc="lower left")

    ax = axes[2]
    mus = np.linspace(0.05, 2.0, 40)
    mnum = [fixed_point(m)[1] for m in mus]
    ax.semilogy(mus, mnum, color=th["series"][0], linewidth=2.2,
                label="integrated $P'(r^*)$", zorder=4)
    ax.semilogy(mus, np.exp(-2*np.pi*mus), color=th["series"][1],
                linewidth=1.4, linestyle=(0, (2, 2)),
                label="$e^{-2\\pi\\mu}$, order 1", zorder=3)
    ms = [float(np.exp(float(ser["lnmult"].subs(ser["mu"], m)))) for m in mus]
    ax.semilogy(mus, ms, color=th["series"][2], linewidth=1.4,
                linestyle=(0, (5, 3)),
                label="$e^{-2\\pi\\mu - \\pi\\mu^3/4}$, series to order 4",
                zorder=3)
    ax.set_ylim(1e-9, 2.0)
    style(ax, th, "$\\mu$", "Floquet multiplier",
          "How far the series reaches")
    legend(ax, th, loc="lower left")

    fig.suptitle("Van der Pol in polar coordinates: one revolution, "
                 "integrated and expanded", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "polar")


if __name__ == "__main__":
    import sys
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    ser = None
    if what in ("all", "checks"):
        ser, _ = checks()
        print()
    if what not in ("all", "figures"):
        sys.exit(0)
    from figures import THEMES
    if ser is None:
        ser = series(5)
    for name, th in THEMES.items():
        print(f"{name}:")
        fig_polar(th, name, ser)
