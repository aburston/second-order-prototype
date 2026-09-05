"""The Lotka-Volterra prototype: predator and prey with piecewise linear rates.

Run ``python3 lotka.py`` to print every number ``LOTKA.md`` quotes and to
write its figures into ``figures/`` as ``lotka-*.png`` in both themes.
``python3 lotka.py checks`` prints the numbers only and ``python3 lotka.py
figures`` writes the figures only. Nothing is cached; a full run takes
about a minute.

The model
---------

Lotka-Volterra, ``u' = alpha u - beta u v``, ``v' = -gamma v + delta u v``,
has one interior equilibrium ``(u*, v*) = (gamma/delta, alpha/beta)``. In
log deviations from it, ``xi = ln(u/u*)`` and ``eta = ln(v/v*)``, it reads

    xi'  = -alpha (e^eta - 1)
    eta' =  gamma (e^xi  - 1)

and the prototype replaces the exponential by a straight line with a floor::

    xi'  = -alpha phi(eta)
    eta' =  gamma phi(xi)          phi(s) = max(s, -s0)

Three parameters carry the dynamics and two more set the population scales:

``alpha``   slope of the prey's per-capita growth rate in ln(predators) at
            the equilibrium, in 1/time
``gamma``   the same for the predator in ln(prey)
``s0``      the floor: a population below e^-s0 of its equilibrium value has
            no further effect on the other. The prey's greatest growth rate
            is alpha s0 and the predator's greatest death rate is gamma s0.
            Lotka-Volterra has those at alpha and gamma, so s0 = 1 is the
            rate matched fit.
``u*, v*``  the equilibrium populations; they scale the populations and
            nothing else

The frequency of small oscillations is ``w0 = sqrt(alpha gamma)`` in both
systems. The floors cut the plane into four regions, each with a linear
field, and every orbit is a chain of elementary arcs: an ellipse where both
populations are within e^-s0 of equilibrium, a parabola where one of them is
scarce, a straight line where both are. ``circuit`` walks those arcs with
closed form transit times and is the exact period.

A third piece, a steeper slope ``k`` above a knee ``s1``, refines the upper
half of the law where the exponential is convex; the walk below handles any
such piecewise linear law, and ``S1, K`` carry the fit the document leads
with (period matched), ``S1_PEAK, K_PEAK`` the alternative (peak matched).

Two extensions use the same pieces. ``c > 0`` adds the prey's own density
dependence ``-c phi(xi)``, the log form of logistic growth, which is a linear
damper on the predator oscillator and makes the equilibrium a stable focus.
``hump=(zeta+, zeta-, xi1)`` makes that density dependence a tent, negative
below a best density ``u* e^xi1`` and positive above it, and inside the
inner region that is the README's offset boundary prototype written in
predator coordinates, limit cycle and all.

Everything below is in units ``alpha = gamma = 1`` unless a line says
otherwise; ``alpha`` and ``gamma`` scale the two axes and the clock.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from figures import THEMES, style, legend, save, newfig

ALPHA, GAMMA, S0 = 1.0, 1.0, 1.0
RTOL, ATOL = 1e-11, 1e-13
TWO_PI = 2.0*np.pi

#: The README's offset boundary cycle at zeta+ = 0.3, zeta- = -0.1, v0 = 1,
#: wn = 1: radius on the section, period, multiplier (README.md, limit_cycle.py).
README_R, README_T, README_M = 2.150651224, 6.367077, 0.5390
ZP, ZM = 0.3, -0.1
XI1_SMALL, XI1_LARGE = 0.3, 1.5
C_DAMP = 0.05


# ------------------------------------------------------------------ model
#: The third piece: above the knee ``S1`` the law rises with slope ``K``.
#: Fitted in ``fit_third_piece`` so the prey peak matches Lotka-Volterra's at
#: the two troughs ``FIT_TROUGHS``; ``K = 1`` is the two piece law.
FIT_TROUGHS = (-2.0, -8.0)      # the peak matched fit equates the prey peaks here
PERIOD_TROUGHS = (-3.0, -8.0)   # the period matched fit equates the periods here
S1_PEAK, K_PEAK = 0.446758, 3.675246     # fit_third_piece(); checks() re-derives both pairs
S1, K = 1.606541, 8.370444               # fit_period_piece(); the fit the document leads with


def pieces(s0=S0, s1=np.inf, k=1.0):
    """The law as intervals ``(lo, hi, a, b)`` on which ``phi(s) = a s + b``."""
    return ((-np.inf, -s0, 0.0, -s0), (-s0, s1, 1.0, 0.0), (s1, np.inf, k, s1*(1.0 - k)))


def phi(s, s0=S0, s1=np.inf, k=1.0):
    """The rate law: a line through the origin, a floor at ``-s0``, slope ``k`` above ``s1``."""
    s = np.asarray(s, float)
    return np.where(s > s1, s1 + k*(s - s1), np.maximum(s, -s0))


def Phi(s, s0=S0, s1=np.inf, k=1.0):
    """Integral of ``phi`` from zero: piecewise quadratic and continuous."""
    s = np.asarray(s, float)
    low = -s0*s - 0.5*s0*s0
    mid = 0.5*s*s
    high = 0.5*s1*s1 + s1*(s - s1) + 0.5*k*(s - s1)**2
    return np.where(s > s1, high, np.where(s >= -s0, mid, low))


def Phi_inv(h, sign, s0=S0, s1=np.inf, k=1.0):
    """The ``s`` of the given sign with ``Phi(s) = h``."""
    if sign < 0:
        return -np.sqrt(2.0*h) if h <= 0.5*s0*s0 else -(h + 0.5*s0*s0)/s0
    if h <= 0.5*s1*s1:
        return np.sqrt(2.0*h)
    # k/2 d^2 + s1 d + s1^2/2 - h = 0 for d = s - s1 > 0
    return s1 + (-s1 + np.sqrt(s1*s1 - 2.0*k*(0.5*s1*s1 - h)))/k


def energy(xi, eta, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """The conserved quantity ``H = gamma Phi(xi) + alpha Phi(eta)``."""
    return gamma*Phi(xi, s0, s1, k) + alpha*Phi(eta, s0, s1, k)


def lv_Phi(s):
    """Lotka-Volterra's counterpart of ``Phi``: ``e^s - 1 - s``."""
    return np.expm1(s) - s


def lv_energy(xi, eta, alpha=ALPHA, gamma=GAMMA):
    return gamma*lv_Phi(xi) + alpha*lv_Phi(eta)


def density(xi, alpha, gamma, s0, c=0.0, hump=None):
    """The prey's own density dependence ``D(xi)``, subtracted from ``xi'``.

    ``c phi(xi)`` is logistic growth in log form, floored like every other
    rate. ``hump = (zp, zm, xi1)`` adds the tent ``2 w0 zeta (phi(xi) - xi1)``
    with ``zeta = zm`` below the best density ``xi1`` and ``zp`` above it.
    Both use the two piece law: the third piece is a refinement of the
    conservative model's upper half and is not carried into them.
    """
    p = float(phi(xi, s0))
    D = c*p
    if hump is not None:
        zp, zm, xi1 = hump
        D = D + 2.0*np.sqrt(alpha*gamma)*(zp if p > xi1 else zm)*(p - xi1)
    return D


def field(alpha=ALPHA, gamma=GAMMA, s0=S0, c=0.0, hump=None, s1=np.inf, k=1.0):
    """Right hand side of the prototype in ``y = [xi, eta]`` for ``solve_ivp``."""
    def f(t, y):
        xi, eta = y
        return [-alpha*float(phi(eta, s0, s1, k)) - density(xi, alpha, gamma, s0, c, hump),
                gamma*float(phi(xi, s0, s1, k))]
    return f


def lv_field(alpha=ALPHA, gamma=GAMMA, c=0.0):
    """Lotka-Volterra in log coordinates, with logistic prey when ``c > 0``."""
    def f(t, y):
        xi, eta = y
        return [-alpha*np.expm1(eta) - c*np.expm1(xi), gamma*np.expm1(xi)]
    return f


def jacobian(f, y, h=1e-6):
    """Central difference Jacobian of ``f(t, y)`` at ``y``."""
    y = np.asarray(y, float)
    J = np.zeros((2, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = h
        J[:, j] = (np.asarray(f(0.0, y + e)) - np.asarray(f(0.0, y - e)))/(2*h)
    return J


# ------------------------------------------------------- exact pieces
# The breakpoints of the law on each axis cut the plane into rectangles,
# and in each the field is affine: ``xi' = p + q eta``, ``eta' = r + s xi``.
# With both slopes nonzero the arc is an ellipse about a centre, with one
# zero it is a parabola, with both zero a straight line. The orbit runs
# counterclockwise. Two piece regions are named 0 inner, 1 predators
# scarce, 2 prey scarce, 3 both scarce, and ``region_name`` keeps those.
def _piece_index(s, ds, s0, s1, k, tol=1e-9):
    """Which interval of the law ``s`` is in, resolved by its motion on a breakpoint."""
    for j, b in enumerate((-s0, s1)):
        if abs(s - b) < tol:
            return j + 1 if ds > 0 else j
    return 0 if s < -s0 else (1 if s < s1 else 2)


def region(xi, eta, alpha, gamma, s0, s1=np.inf, k=1.0):
    """The region as ``(i_xi, i_eta)``, piece indices 0 floor, 1 line, 2 steep."""
    dxi = -alpha*float(phi(eta, s0, s1, k))
    deta = gamma*float(phi(xi, s0, s1, k))
    return _piece_index(xi, dxi, s0, s1, k), _piece_index(eta, deta, s0, s1, k)


def region_name(reg):
    """0 inner, 1 predators scarce, 2 prey scarce, 3 both scarce, as in the two piece law."""
    i, j = reg
    return (2 if i == 0 else 0) + (1 if j == 0 else 0)


def _roots(a, b, c, eps):
    """Real roots of ``a t^2 + b t + c = 0`` greater than ``eps``, ascending."""
    if a == 0.0:
        return [] if b == 0.0 else [t for t in (-c/b,) if t > eps]
    disc = b*b - 4.0*a*c
    if disc < 0.0:
        return []
    q = -0.5*(b + np.copysign(np.sqrt(disc), b if b != 0.0 else 1.0))
    roots = (q/a, c/q) if q != 0.0 else (0.0, 0.0)
    return sorted(t for t in roots if t > eps)


def _affine(reg, alpha, gamma, s0, s1, k):
    """``(p, q, r, s, bounds)``: the region's field and its rectangle."""
    px = pieces(s0, s1, k)
    lo_x, hi_x, ax, bx = px[reg[0]]
    lo_y, hi_y, ay, by = px[reg[1]]
    return -alpha*by, -alpha*ay, gamma*bx, gamma*ax, (lo_x, hi_x, lo_y, hi_y)


def step(xi, eta, alpha, gamma, s0, s1=np.inf, k=1.0, stop=None, eps=1e-12):
    """Cross one region exactly. Returns ``(region, dt, xi1, eta1, exit)``.

    ``exit`` names the boundary crossed (``'xi'`` or ``'eta'``), or
    ``'start'`` when the orbit returns to the section point given by
    ``stop = (region, angle)``, or ``'loop'`` after a full turn.
    """
    reg = region(xi, eta, alpha, gamma, s0, s1, k)
    p, q, r, s, (lo_x, hi_x, lo_y, hi_y) = _affine(reg, alpha, gamma, s0, s1, k)
    if q != 0.0 and s != 0.0:
        xc, yc = -r/s, -p/q
        W = np.sqrt(-q*s)
        X, Y = np.sqrt(s)*(xi - xc), np.sqrt(-q)*(eta - yc)
        R, th0 = np.hypot(X, Y), np.arctan2(Y, X)
        cands = []
        for b, kind in ((lo_x, "xi"), (hi_x, "xi")):
            c = np.sqrt(s)*(b - xc)
            if np.isfinite(b) and abs(c) < R:      # a tangency is not a crossing
                ac = np.arccos(c/R)
                cands += [(ac, kind, b), (-ac, kind, b)]
        for b, kind in ((lo_y, "eta"), (hi_y, "eta")):
            c = np.sqrt(-q)*(b - yc)
            if np.isfinite(b) and abs(c) < R:
                asn = np.arcsin(c/R)
                cands += [(asn, kind, b), (np.pi - asn, kind, b)]
        if stop is not None and stop[0] == reg:
            cands.append((stop[1], "start", None))
        d, kind, b = TWO_PI, "loop", None
        for th, kd, bb in cands:
            dd = (th - th0) % TWO_PI
            if dd < 1e-7 or dd > TWO_PI - 1e-7:
                continue          # the crossing we are standing on
            if dd < d:
                d, kind, b = dd, kd, bb
        th1 = th0 + d
        xi1, eta1 = xc + R*np.cos(th1)/np.sqrt(s), yc + R*np.sin(th1)/np.sqrt(-q)
        if kind == "xi":
            xi1 = b
        elif kind == "eta":
            eta1 = b
        elif kind == "start":
            eta1 = 0.0
        return reg, d/W, xi1, eta1, kind
    # polynomial motion: xi(t) = xi + (p + q eta) t + q r t^2 / 2 when s == 0, etc.
    cx = (xi, p + q*eta, 0.5*q*r if s == 0.0 else 0.0)
    cy = (eta, r + s*xi, 0.5*s*p if q == 0.0 else 0.0)
    if s != 0.0:          # eta linear in xi which is quadratic: xi is linear here
        cx = (xi, p + q*eta, 0.0)
        cy = (eta, r + s*xi, 0.5*s*(p + q*eta))
    if q != 0.0:
        cy = (eta, r + s*xi, 0.0)
        cx = (xi, p + q*eta, 0.5*q*(r + s*xi))
    best = (np.inf, None, None)
    for coeffs, bounds, kind in ((cx, (lo_x, hi_x), "xi"), (cy, (lo_y, hi_y), "eta")):
        for b in bounds:
            if not np.isfinite(b):
                continue
            ts = _roots(coeffs[2], coeffs[1], coeffs[0] - b, eps)
            if ts and ts[0] < best[0]:
                best = (ts[0], kind, b)
    dt, kind, b = best
    xi1 = cx[0] + cx[1]*dt + cx[2]*dt*dt
    eta1 = cy[0] + cy[1]*dt + cy[2]*dt*dt
    if kind == "xi":
        xi1 = b
    else:
        eta1 = b
    return reg, dt, xi1, eta1, kind


def sample(reg, xi, eta, dt, alpha, gamma, s0, s1=np.inf, k=1.0, n=200):
    """The closed form arc through one region, sampled at ``n`` times."""
    t = np.linspace(0.0, dt, n)
    p, q, r, s, _ = _affine(reg, alpha, gamma, s0, s1, k)
    if q != 0.0 and s != 0.0:
        xc, yc = -r/s, -p/q
        W = np.sqrt(-q*s)
        X, Y = np.sqrt(s)*(xi - xc), np.sqrt(-q)*(eta - yc)
        R, th = np.hypot(X, Y), np.arctan2(Y, X) + W*t
        return t, xc + R*np.cos(th)/np.sqrt(s), yc + R*np.sin(th)/np.sqrt(-q)
    if s != 0.0:
        return t, xi + (p + q*eta)*t, eta + (r + s*xi)*t + 0.5*s*(p + q*eta)*t*t
    if q != 0.0:
        return t, xi + (p + q*eta)*t + 0.5*q*(r + s*xi)*t*t, eta + (r + s*xi)*t
    return t, xi + p*t, eta + r*t


def circuit(H, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """One period from the prey peak: ``(T, [(region, dt, xi, eta), ...])``."""
    xi, eta = Phi_inv(H/gamma, +1, s0, s1, k), 0.0
    reg0 = region(xi, eta, alpha, gamma, s0, s1, k)
    p, q, r, s, _ = _affine(reg0, alpha, gamma, s0, s1, k)
    th0 = np.arctan2(np.sqrt(-q)*(eta + p/q), np.sqrt(s)*(xi + r/s))
    pieces_, total = [], 0.0
    for _ in range(24):
        reg, dt, xi1, eta1, kind = step(xi, eta, alpha, gamma, s0, s1, k, stop=(reg0, th0))
        pieces_.append((reg, dt, xi, eta))
        total += dt
        xi, eta = xi1, eta1
        if kind in ("start", "loop"):
            break
    return total, pieces_


def orbit(H, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0, n=300):
    """The closed orbit of energy ``H`` as arrays ``(t, xi, eta)``."""
    T, pieces_ = circuit(H, alpha, gamma, s0, s1, k)
    ts, xs, ys, t0 = [], [], [], 0.0
    for reg, dt, xi, eta in pieces_:
        t, x, y = sample(reg, xi, eta, dt, alpha, gamma, s0, s1, k, n)
        ts.append(t + t0)
        xs.append(x)
        ys.append(y)
        t0 += dt
    return np.concatenate(ts), np.concatenate(xs), np.concatenate(ys)


def period(H, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """Exact period, the sum of the circuit's transit times."""
    return circuit(H, alpha, gamma, s0, s1, k)[0]


def period_formula(H, alpha=ALPHA, gamma=GAMMA, s0=S0):
    """The two piece law's period assembled by hand from the four transit times.

    Below both floors it is ``2 pi / w0``; touching one floor swaps an arc
    of the ellipse for a parabola; beyond the corner energy ``Hc`` all four
    regions are visited and the corner adds a straight line.
    """
    w0 = np.sqrt(alpha*gamma)
    r2 = 2.0*H
    Hc = 0.5*(alpha + gamma)*s0*s0
    wall, floor = r2 > gamma*s0*s0, r2 > alpha*s0*s0
    if not (wall or floor):
        return TWO_PI/w0
    sx, sy = np.sqrt(gamma)*s0/np.sqrt(r2), np.sqrt(alpha)*s0/np.sqrt(r2)
    eta_w = np.sqrt((r2 - gamma*s0*s0)/alpha) if wall else 0.0     # eta where the ellipse meets the wall
    xi_f = np.sqrt((r2 - alpha*s0*s0)/gamma) if floor else 0.0     # |xi| where it meets the floor
    if H <= Hc:
        T = TWO_PI/w0
        if wall:
            T += 2.0*eta_w/(gamma*s0) - 2.0*np.arccos(sx)/w0
        if floor:
            T += 2.0*xi_f/(alpha*s0) - 2.0*np.arccos(sy)/w0
        return T
    return ((0.5*np.pi + np.arcsin(sx) + np.arcsin(sy))/w0
            + (s0 + xi_f)/(alpha*s0) + (s0 + eta_w)/(gamma*s0)
            + (H - Hc)/(alpha*gamma*s0*s0))


def period_integrated(H, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """The period by direct integration from the prey peak, for checking."""
    f = field(alpha, gamma, s0, s1=s1, k=k)
    xi0 = Phi_inv(H/gamma, +1, s0, s1, k)

    def ev(t, y):
        return y[1]
    ev.direction = 1
    Tg = period(H, alpha, gamma, s0, s1, k)
    s = solve_ivp(f, (0.0, 1.5*Tg + 1.0), [xi0, 0.0], events=ev, rtol=RTOL, atol=ATOL,
                  method="DOP853", max_step=0.02)
    t = [t for t in s.t_events[0] if t > 1e-6]
    return t[0] if t else np.nan


# ------------------------------------------------------- amplitudes
def H_from_trough(xi_min, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """Energy of the prototype orbit whose prey trough is ``xi_min``."""
    return gamma*float(Phi(xi_min, s0, s1, k))


def extremes(H, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """``(xi_min, xi_max, eta_min, eta_max)`` of the prototype orbit."""
    return (Phi_inv(H/gamma, -1, s0, s1, k), Phi_inv(H/gamma, +1, s0, s1, k),
            Phi_inv(H/alpha, -1, s0, s1, k), Phi_inv(H/alpha, +1, s0, s1, k))


def fit_third_piece(troughs=FIT_TROUGHS, alpha=ALPHA, gamma=GAMMA, s0=S0):
    """``(s1, k)`` such that the prey peak equals Lotka-Volterra's at both troughs.

    At a trough ``xi_min`` the prototype's energy is ``gamma Phi(xi_min)``,
    unchanged by the third piece, and its peak solves ``Phi(xi_max) = H/gamma``;
    Lotka-Volterra's peak solves ``e^x - x - 1 = e^xi_min - xi_min - 1``.
    """
    from scipy.optimize import fsolve
    targets = []
    for xm in troughs:
        h = float(Phi(xm, s0))
        xl = brentq(lambda x, c: lv_Phi(x) - c, 0.0, 60.0, args=(lv_Phi(xm),), xtol=1e-14)
        targets.append((h, xl))

    def resid(v):
        s1, k = v
        return [float(Phi(xl, s0, s1, k)) - h for h, xl in targets]
    s1, k = fsolve(resid, [0.5, 3.0], xtol=1e-13)
    return float(s1), float(k)


def fit_period_piece(troughs=PERIOD_TROUGHS, alpha=ALPHA, gamma=GAMMA, s0=S0, lv_T=None):
    """``(s1, k)`` such that the period equals Lotka-Volterra's at both troughs.

    ``lv_T`` may carry the Lotka-Volterra periods already integrated, keyed
    by trough; otherwise they are integrated here.
    """
    from scipy.optimize import fsolve
    if lv_T is None:
        lv_T = {xm: lv_run(xm, alpha, gamma)[0] for xm in troughs}

    def resid(v):
        s1, k = v
        return [period(H_from_trough(xm, alpha, gamma, s0, s1, k), alpha, gamma, s0, s1, k) - lv_T[xm]
                for xm in troughs]
    s1, k = fsolve(resid, [1.5, 7.0], xtol=1e-12)
    return float(s1), float(k)


def lv_extremes(xi_min, alpha=ALPHA, gamma=GAMMA):
    """``(H, xi_max, eta_min, eta_max)`` of the Lotka-Volterra orbit with that prey trough."""
    H = gamma*lv_Phi(xi_min)

    def g(s, c):
        return lv_Phi(s) - c
    xi_max = brentq(g, 0.0, 60.0, args=(H/gamma,), xtol=1e-14)
    eta_max = brentq(g, 0.0, 60.0, args=(H/alpha,), xtol=1e-14)
    eta_min = brentq(g, -1e4, 0.0, args=(H/alpha,), xtol=1e-14)
    return H, xi_max, eta_min, eta_max


def lv_run(xi_min, alpha=ALPHA, gamma=GAMMA, n=3000):
    """One Lotka-Volterra period from the prey peak: ``(T, lag, t, xi, eta)``.

    ``lag`` is the time from the prey peak to the predator peak.
    """
    H, xi_max, eta_min, _ = lv_extremes(xi_min, alpha, gamma)
    f = lv_field(alpha, gamma)

    def back(t, y):
        return y[1]
    back.direction = 1

    def peak(t, y):
        return y[0]
    peak.direction = -1
    tmax = 4.0*(abs(xi_min)/alpha + abs(eta_min)/gamma) + 40.0/np.sqrt(alpha*gamma)
    s = solve_ivp(f, (0.0, tmax), [xi_max, 0.0], events=(back, peak), rtol=RTOL, atol=ATOL,
                  method="DOP853", dense_output=True)
    T = [t for t in s.t_events[0] if t > 1e-6][0]
    lag = s.t_events[1][0]
    t = np.linspace(0.0, T, n)
    y = s.sol(t)
    return T, lag, t, y[0], y[1]


def lag_integrated(H, alpha=ALPHA, gamma=GAMMA, s0=S0, s1=np.inf, k=1.0):
    """Prey peak to predator peak in the prototype, by integration."""
    f = field(alpha, gamma, s0, s1=s1, k=k)

    def peak(t, y):
        return y[0]
    peak.direction = -1
    s = solve_ivp(f, (0.0, 10.0/np.sqrt(alpha*gamma)), [Phi_inv(H/gamma, +1, s0, s1, k), 0.0],
                  events=peak, rtol=RTOL, atol=ATOL, method="DOP853", max_step=0.02)
    return s.t_events[0][0]


# ------------------------------------------------------- damped, and the cycle
def peaks(f, y0, n, tmax, direction=-1):
    """Successive crossings of ``xi = 0`` with ``xi`` falling: the predator peaks.

    Returns ``(times, states)`` of the first ``n`` after the start.
    """
    def ev(t, y):
        return y[0]
    ev.direction = direction
    s = solve_ivp(f, (0.0, tmax), y0, events=ev, rtol=RTOL, atol=ATOL, method="DOP853",
                  max_step=0.05)
    keep = [i for i, t in enumerate(s.t_events[0]) if t > 1e-6][:n]
    return s.t_events[0][keep], s.y_events[0][keep]


def troughs(f, y0, n, tmax):
    """Successive prey troughs (``eta = 0`` with ``eta`` falling): ``(times, xi)``."""
    def ev(t, y):
        return y[1]
    ev.direction = -1
    s = solve_ivp(f, (0.0, tmax), y0, events=ev, rtol=RTOL, atol=ATOL, method="DOP853",
                  max_step=0.05)
    keep = [i for i, t in enumerate(s.t_events[0]) if t > 1e-6][:n]
    return s.t_events[0][keep], s.y_events[0][keep][:, 0]


def energy_loss(H, c, alpha=ALPHA, gamma=GAMMA, s0=S0):
    """Averaging estimate of the energy lost per cycle, ``c gamma oint phi(xi)^2 dt``.

    Evaluated on the conservative orbit of energy ``H``, arc by arc.
    """
    _, pieces = circuit(H, alpha, gamma, s0)
    loss = 0.0
    for reg, dt, xi, eta in pieces:
        t, x, _ = sample(reg, xi, eta, dt, alpha, gamma, s0, n=4001)
        loss += np.trapezoid(phi(x, s0)**2, t)
    return c*gamma*loss


def trough_map(xi_min, c, n, alpha=ALPHA, gamma=GAMMA, s0=S0):
    """Successive prey troughs predicted by the energy map, from ``xi_min``."""
    out, H = [], H_from_trough(xi_min, alpha, gamma, s0)
    for _ in range(n):
        H = H - energy_loss(H, c, alpha, gamma, s0)
        if H <= 0:
            out.append(0.0)
            continue
        out.append(extremes(H, alpha, gamma, s0)[0])
    return out


def eta_eq(hump, alpha=ALPHA, gamma=GAMMA):
    """Predator equilibrium with the tent: ``2 zeta- w0 xi1 / alpha``."""
    zp, zm, xi1 = hump
    return 2.0*zm*np.sqrt(alpha*gamma)*xi1/alpha


def return_map(r, hump, alpha=ALPHA, gamma=GAMMA, s0=S0):
    """One return to the section ``{xi = 0, eta > eta_eq}``: ``(r', T, xi_min, eta_min)``."""
    f = field(alpha, gamma, s0, hump=hump)
    ye = eta_eq(hump, alpha, gamma)

    def ev(t, y):
        return y[0]
    ev.direction = -1
    s = solve_ivp(f, (0.0, 60.0/np.sqrt(alpha*gamma)), [0.0, ye + r], events=ev, rtol=RTOL,
                  atol=ATOL, method="DOP853", max_step=0.02, dense_output=True)
    i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
    if not i:
        return np.nan, np.nan, np.nan, np.nan
    T = s.t_events[0][i[0]]
    y = s.sol(np.linspace(0.0, T, 4000))
    return s.y_events[0][i[0]][1] - ye, T, y[0].min(), y[1].min()


def limit_cycle(hump, r0, alpha=ALPHA, gamma=GAMMA, s0=S0, n=400, tol=1e-11):
    """Iterate the return map to its fixed point: ``(r*, T, multiplier, xi_min, eta_min)``."""
    r = r0
    for _ in range(n):
        rn, T, xmin, ymin = return_map(r, hump, alpha, gamma, s0)
        if not np.isfinite(rn) or rn > 1e6 or rn < 1e-9:
            return rn, T, np.nan, xmin, ymin
        if abs(rn - r) < tol*max(1.0, abs(r)):
            r = rn
            break
        r = rn
    h = 1e-6
    m = (return_map(r + h, hump, alpha, gamma, s0)[0]
         - return_map(r - h, hump, alpha, gamma, s0)[0])/(2*h)
    return r, T, m, xmin, ymin


# ------------------------------------------------- the same equation elsewhere
# Each entry builds the system as its own field writes it, in its own
# variables, and returns ``(f, (u*, v*), mapping)`` with the dictionary onto
# the prototype's ``alpha``, ``gamma`` and density dependence ``c``.
def native_lv(alpha, beta, gamma, delta):
    """Ecology: ``u' = alpha u - beta u v``, ``v' = -gamma v + delta u v``."""
    def f(t, y):
        u, v = y
        return [alpha*u - beta*u*v, -gamma*v + delta*u*v]
    return f, (gamma/delta, alpha/beta), dict(alpha=alpha, gamma=gamma, c=0.0)


def native_lotka(k1A, k2, k3):
    """Lotka's mechanism ``A + X -> 2X``, ``X + Y -> 2Y``, ``Y -> B`` with ``A`` held fixed."""
    def f(t, y):
        X, Y = y
        return [k1A*X - k2*X*Y, k2*X*Y - k3*Y]
    return f, (k3/k2, k1A/k2), dict(alpha=k1A, gamma=k3, c=0.0)


def native_sir(mu, gamma_r, R0, N=1.0):
    """SIR with births and deaths: ``S' = mu N - mu S - beta S I``, ``I' = beta S I - (gamma + mu) I``."""
    beta = R0*(gamma_r + mu)
    Sstar, Istar = N/R0, mu*(R0 - 1.0)/beta

    def f(t, y):
        S, I = y
        return [mu*N - mu*S - beta*S*I, beta*S*I - (gamma_r + mu)*I]
    return f, (Sstar, Istar), dict(alpha=beta*Istar, gamma=gamma_r + mu, c=mu*R0)


def native_laser(gpar, kappa, p, G=1.0):
    """Class B laser rate equations: ``N' = P - g N - G N n``, ``n' = (G N - kappa) n``."""
    Nstar, nstar = kappa/G, gpar*(p - 1.0)/G
    P = p*gpar*Nstar

    def f(t, y):
        N, n = y
        return [P - gpar*N - G*N*n, (G*N - kappa)*n]
    return f, (Nstar, nstar), dict(alpha=G*nstar, gamma=kappa, c=gpar*p)


def native_goodwin(sigma, a, n, rho, gw):
    """Goodwin's growth cycle: employment rate ``e`` and wage share ``w``.

    ``e'/e = (1 - w)/sigma - (a + n)``, ``w'/w = rho e - (a + gw)``: capital
    to output ratio ``sigma``, productivity growth ``a``, labour force growth
    ``n``, Phillips curve ``rho e - gw``.
    """
    estar, wstar = (a + gw)/rho, 1.0 - sigma*(a + n)

    def f(t, y):
        e, w = y
        return [e*((1.0 - w)/sigma - (a + n)), w*(rho*e - (a + gw))]
    return f, (estar, wstar), dict(alpha=wstar/sigma, gamma=rho*estar, c=0.0)


#: The cases the tables and the figure use. Times: SIR in years (life
#: expectancy 70 y, infectious period two weeks, R0 = 15); the laser in
#: seconds (upper state lifetime 230 us, cavity lifetime 20 ns, pumped at
#: twice threshold: an Nd:YAG rod); Goodwin in years.
EXAMPLES = (
    ("ecology, Lotka-Volterra", native_lv(1.0, 1.0, 1.0, 1.0)),
    ("chemistry, Lotka's mechanism", native_lotka(1.0, 1.0, 1.0)),
    ("epidemics, SIR with births", native_sir(1.0/70.0, 365.0/14.0, 15.0)),
    ("laser, class B rate equations", native_laser(1.0/230e-6, 1.0/20e-9, 2.0)),
    ("economics, Goodwin's cycle", native_goodwin(3.0, 0.02, 0.01, 0.5, 0.44)),
)


def measure(f, eq, mapping, amp=0.01, n=4):
    """Period and decrement of small oscillations about ``eq`` in the native system.

    Kicks the predator variable by the fraction ``amp`` and records the next
    ``n`` predator peaks, the crossings of ``u = u*`` with ``u`` falling.
    Returns ``(mean period, mean ratio of successive peak deviations)``.
    """
    ustar, vstar = eq
    w0 = np.sqrt(mapping["alpha"]*mapping["gamma"])
    T0 = TWO_PI/w0

    def ev(t, y):
        return y[0] - ustar
    ev.direction = -1
    s = solve_ivp(f, (0.0, 1.3*(n + 1)*T0), [ustar, vstar*(1.0 + amp)], events=ev,
                  rtol=1e-11, atol=[1e-13*ustar, 1e-13*vstar], method="DOP853", max_step=T0/50)
    t = [tt for tt in s.t_events[0] if tt > 1e-6*T0]
    v = np.array([yy[1] for tt, yy in zip(s.t_events[0], s.y_events[0]) if tt > 1e-6*T0])
    t = np.array(t)
    return np.diff(t).mean(), ((v[1:] - vstar)/(v[:-1] - vstar)).mean()


def predicted(mapping):
    """``(w0, zeta, damped period, decrement per cycle)`` from the dictionary."""
    w0 = np.sqrt(mapping["alpha"]*mapping["gamma"])
    z = mapping["c"]/(2.0*w0)
    return w0, z, TWO_PI/(w0*np.sqrt(1.0 - z*z)), np.exp(-TWO_PI*z/np.sqrt(1.0 - z*z))


def native_versus_prototype(f, eq, mapping, y0_log, T, n=4000, s0=S0, s1=np.inf, k=1.0):
    """Integrate the native system and the prototype from the same log start.

    Returns ``(t, native eta, prototype eta)``: the predator's log deviation
    in each, the native one converted from its own variables.
    """
    ustar, vstar = eq
    y0 = [ustar*np.exp(y0_log[0]), vstar*np.exp(y0_log[1])]
    t = np.linspace(0.0, T, n)
    s = solve_ivp(f, (0.0, T), y0, t_eval=t, rtol=1e-10, atol=[1e-14*ustar, 1e-14*vstar],
                  method="DOP853", max_step=T/2000)
    p = solve_ivp(field(mapping["alpha"], mapping["gamma"], s0, c=mapping["c"], s1=s1, k=k), (0.0, T), y0_log,
                  t_eval=t, rtol=1e-10, atol=1e-12, method="DOP853", max_step=T/2000)
    return t, np.log(s.y[1]/vstar), p.y[1]


# ---------------------------------------------------------------- figures
def _floors(ax, th, s0=S0):
    for fn in (ax.axvline, ax.axhline):
        fn(-s0, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)


def fig_law(th, name):
    """The rate law against the exponential, in log and in population coordinates."""
    fig, axes = newfig(th, 1, 2, figsize=(10.0, 4.0))
    c0, c1, _ = th["series"]
    s = np.linspace(-4.0, 1.6, 600)
    axes[0].plot(s, np.expm1(s), color=c0, linewidth=2.0, label="Lotka-Volterra  $e^s - 1$", zorder=3)
    axes[0].plot(s, phi(s), color=c1, linewidth=1.8, label="prototype  $\\phi(s) = \\max(s, -s_0)$", zorder=4)
    axes[0].axvline(-S0, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    axes[0].axhline(0, color=th["grid"], linewidth=0.8)
    axes[0].annotate("floor $-s_0$", xy=(-S0, 1.6), xytext=(6, -2), textcoords="offset points",
                     fontsize=8, color=th["ink2"])
    axes[0].set_ylim(-1.6, 3.4)
    style(axes[0], th, "$s$  (log deviation of the other population)", "rate law",
          "In log coordinates: one slope at the origin, one floor")
    legend(axes[0], th, loc="upper left")

    v = np.linspace(0.0, 3.0, 600)
    with np.errstate(divide="ignore"):
        proto = ALPHA*np.minimum(S0, -np.log(v))
    axes[1].plot(v, ALPHA*(1.0 - v), color=c0, linewidth=2.0, label="Lotka-Volterra  $\\alpha(1 - v/v^*)$", zorder=3)
    axes[1].plot(v, proto, color=c1, linewidth=1.8, label="prototype  $\\alpha\\min(s_0, \\ln(v^*/v))$", zorder=4)
    axes[1].axvline(np.exp(-S0), color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    axes[1].axhline(0, color=th["grid"], linewidth=0.8)
    axes[1].annotate("$v = v^* e^{-s_0}$", xy=(np.exp(-S0), 1.25), xytext=(6, 0), textcoords="offset points",
                     fontsize=8, color=th["ink2"])
    axes[1].set_ylim(-2.2, 1.4)
    style(axes[1], th, "predators  $v / v^*$", "prey per-capita growth rate  $/\\alpha$",
          "In populations: logarithmic and capped, against linear")
    legend(axes[1], th, loc="lower left")
    fig.suptitle("The prototype's rate law: Lotka-Volterra's exponential replaced by a line with a floor",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-law")


PHASE_TROUGHS = (-0.5, -1.5, -3.0)


def fig_phase(th, name):
    """Orbits at matched prey troughs, in log coordinates and in populations."""
    fig, axes = newfig(th, 1, 2, figsize=(10.4, 4.8))
    c0, c1, _ = th["series"]
    for k, xm in enumerate(PHASE_TROUGHS):
        _, x, y = orbit(H_from_trough(xm))
        _, _, _, xl, yl = lv_run(xm)
        lab_p = "prototype" if k == 0 else None
        lab_l = "Lotka-Volterra" if k == 0 else None
        axes[0].plot(xl, yl, color=c0, linewidth=1.9, zorder=3, label=lab_l)
        axes[0].plot(x, y, color=c1, linewidth=1.5, zorder=4, label=lab_p)
        axes[1].plot(np.exp(xl), np.exp(yl), color=c0, linewidth=1.9, zorder=3, label=lab_l)
        axes[1].plot(np.exp(x), np.exp(y), color=c1, linewidth=1.5, zorder=4, label=lab_p)
        axes[0].annotate(f"$e^{{{xm:g}}}$", xy=(x.min(), 0.0), xytext=(-3, (6, -13, 6)[k]),
                         textcoords="offset points", ha="right", fontsize=8, color=th["ink2"])
    _floors(axes[0], th)
    axes[0].plot([0], [0], "o", color=th["ink"], markersize=5, zorder=6)
    axes[1].plot([1], [1], "o", color=th["ink"], markersize=5, zorder=6)
    axes[1].axvline(np.exp(-S0), color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    axes[1].axhline(np.exp(-S0), color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    axes[0].annotate("prey scarce", xy=(-2.9, 2.2), fontsize=8, color=th["ink2"])
    axes[0].annotate("labels: prey trough $u_{min}/u^*$", xy=(0.03, 0.03), xycoords="axes fraction",
                     fontsize=7.5, color=th["ink2"])
    axes[0].annotate("predators\nscarce", xy=(1.0, -2.6), fontsize=8, color=th["ink2"])
    axes[0].annotate("both scarce", xy=(-2.9, -2.6), fontsize=8, color=th["ink2"])
    axes[0].set_xlim(-3.8, 2.7)
    axes[0].set_ylim(-3.6, 2.7)
    axes[1].set_xlim(0, 10)
    axes[1].set_ylim(0, 10)
    style(axes[0], th, "$\\xi = \\ln(u/u^*)$", "$\\eta = \\ln(v/v^*)$",
          "Log coordinates: the floors cut the plane into four regions")
    style(axes[1], th, "prey  $u/u^*$", "predators  $v/v^*$", "The same orbits in populations")
    legend(axes[0], th, loc="upper right")
    legend(axes[1], th, loc="upper right")
    fig.suptitle("Closed orbits at the same prey trough, $\\alpha = \\gamma = 1$, $s_0 = 1$",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-phase")


TIME_TROUGH = -4.0
REGION_LABEL = {1: "predators scarce", 2: "prey scarce", 3: "both scarce"}


def fig_time(th, name):
    """Populations against time over one period, prototype beside Lotka-Volterra."""
    fig, axes = newfig(th, 1, 2, figsize=(10.4, 4.0), sharey=True)
    c0, c1, _ = th["series"]
    H = H_from_trough(TIME_TROUGH)
    T, pieces = circuit(H)
    t, x, y = orbit(H)
    Tl, _, tl, xl, yl = lv_run(TIME_TROUGH)
    T0 = TWO_PI/np.sqrt(ALPHA*GAMMA)
    for ax, tt, xx, yy, TT, title in ((axes[0], t, x, y, T, "Prototype"),
                                      (axes[1], tl, xl, yl, Tl, "Lotka-Volterra")):
        # start each trace at the prey trough so the two panels line up
        i0 = int(np.argmin(xx))
        tt = np.concatenate([tt[i0:] - tt[i0], tt[:i0] + TT - tt[i0]])
        xx = np.concatenate([xx[i0:], xx[:i0]])
        yy = np.concatenate([yy[i0:], yy[:i0]])
        ax.plot(tt/T0, np.exp(xx), color=c0, linewidth=1.9, label="prey  $u/u^*$", zorder=4)
        ax.plot(tt/T0, np.exp(yy), color=c1, linewidth=1.9, label="predators  $v/v^*$", zorder=4)
        ax.axhline(np.exp(-S0), color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
        ax.set_yscale("log")
        style(ax, th, "time  $/\\,T_0$", "population / equilibrium" if ax is axes[0] else "",
              f"{title}: period $T = {TT/T0:.2f}\\,T_0$")
    # shade the regions of the prototype's circuit, shifted to start at the trough
    tt0 = 0.0
    spans = []
    for reg, dt, _, _ in pieces:
        spans.append((region_name(reg), tt0, tt0 + dt))
        tt0 += dt
    t_trough = t[int(np.argmin(x))]
    shown = []
    for reg, a, b in spans:
        if reg == 0:
            continue
        for shift in (-T, 0.0, T):
            lo, hi = a - t_trough + shift, b - t_trough + shift
            lo, hi = max(lo, 0.0), min(hi, T)
            if hi > lo:
                axes[0].axvspan(lo/T0, hi/T0, color=th["grid"], alpha=0.55, zorder=1)
                shown.append((reg, lo, hi))
    for reg in set(r for r, _, _ in shown):
        lo, hi = max(((lo, hi) for r, lo, hi in shown if r == reg), key=lambda p: p[1] - p[0])
        axes[0].annotate(REGION_LABEL[reg], xy=((lo + hi)/(2*T0), 25), ha="center",
                         fontsize=7.5, color=th["ink2"])
    axes[0].set_ylim(np.exp(-9.5), 60)
    legend(axes[0], th, loc="lower right")
    fig.suptitle(f"One cycle from a prey trough of $e^{{{TIME_TROUGH:g}}}u^*$: straight lines "
                 "on a log axis are the floors at work", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-time")


def fig_period(th, name, data):
    """Period, prey peak and predator lag against the prey trough, both systems."""
    fig, axes = newfig(th, 1, 3, figsize=(12.6, 4.0))
    c0, c1, c2 = th["series"]
    T0 = TWO_PI/np.sqrt(ALPHA*GAMMA)
    xm = data["xm_dense"]
    axes[0].plot(-xm, data["lv_T"]/T0, color=c0, linewidth=2.0, label="Lotka-Volterra (integrated)", zorder=3)
    axes[0].plot(-xm, data["pr_T"]/T0, color=c1, linewidth=1.8, label="prototype, closed form", zorder=4)
    axes[0].plot(-np.array(data["xm_pts"]), np.array(data["pr_Tn"])/T0, "o", color=c1, markersize=4.5,
                 markerfacecolor=th["surface"], label="prototype, integrated", zorder=5)
    axes[0].plot(-xm, (-xm/ALPHA)/T0, color=th["ink2"], linewidth=1.0, linestyle=(0, (4, 3)),
                 label="$|\\xi_{min}|/\\alpha$, the common asymptote", zorder=2)
    axes[0].axvline(S0, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    axes[0].annotate("corner", xy=(S0, 0.98), xytext=(4, 0), textcoords="offset points",
                     fontsize=8, color=th["ink2"])
    axes[0].set_ylim(0.9, None)
    style(axes[0], th, "prey trough  $-\\xi_{min} = \\ln(u^*/u_{min})$", "period  $T / T_0$",
          "Flat to the corner, then linear in the trough")
    legend(axes[0], th, loc="upper left")

    axes[1].plot(-xm, data["lv_xmax"], color=c0, linewidth=2.0, label="Lotka-Volterra: $\\sim \\ln|\\xi_{min}|$", zorder=3)
    axes[1].plot(-xm, data["pr_xmax"], color=c1, linewidth=1.8,
                 label="prototype: $\\sqrt{2 s_0 |\\xi_{min}| - s_0^2}$", zorder=4)
    axes[1].axvline(S0, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    style(axes[1], th, "prey trough  $-\\xi_{min}$", "prey peak  $\\xi_{max} = \\ln(u_{max}/u^*)$",
          "Where the two systems part: the prey peak")
    legend(axes[1], th, loc="upper left")

    axes[2].plot(-xm, data["lv_lag"]/T0, color=c0, linewidth=2.0, label="Lotka-Volterra", zorder=3)
    axes[2].plot(-xm, data["pr_lag"]/T0, color=c1, linewidth=1.8, label="prototype: exactly $T_0/4$", zorder=4)
    axes[2].plot(-xm, data["lv_lag"]/data["lv_T"], color=c0, linewidth=2.0, linestyle=(0, (2, 2)), zorder=3)
    axes[2].plot(-xm, data["pr_lag"]/data["pr_T"], color=c1, linewidth=1.8, linestyle=(0, (2, 2)), zorder=4)
    axes[2].plot([], [], color=th["ink2"], linewidth=1.6, label="lag $/\\,T_0$")
    axes[2].plot([], [], color=th["ink2"], linewidth=1.6, linestyle=(0, (2, 2)), label="lag $/\\,T$, as a fraction of the cycle")
    axes[2].set_ylim(0, 0.3)
    style(axes[2], th, "prey trough  $-\\xi_{min}$", "prey peak to predator peak",
          "The predator's lag behind the prey")
    legend(axes[2], th, loc="upper right")
    fig.suptitle("Period, shape and phase against amplitude, $\\alpha = \\gamma = 1$, $s_0 = 1$",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-period")


def _traj(f, y0, T, n=6000):
    s = solve_ivp(f, (0.0, T), y0, t_eval=np.linspace(0.0, T, n), rtol=1e-10, atol=1e-12,
                  method="DOP853", max_step=0.05)
    return s.y[0], s.y[1]


def fig_damped(th, name, cyc_small, cyc_large):
    """Density dependence: the spiral in, and the two limit cycles of the tent."""
    fig, axes = newfig(th, 1, 3, figsize=(12.6, 4.2))
    c0, c1, c2 = th["series"]
    # left: logistic prey, prototype against Lotka-Volterra, same c, same start
    y0 = [-3.0, 0.0]
    xl, yl = _traj(lv_field(c=C_DAMP), y0, 60.0)
    xp, yp = _traj(field(c=C_DAMP), y0, 60.0)
    axes[0].plot(xl, yl, color=c0, linewidth=1.5, label="Lotka-Volterra with logistic prey", zorder=3)
    axes[0].plot(xp, yp, color=c1, linewidth=1.3, label="prototype", zorder=4)
    _floors(axes[0], th)
    axes[0].plot([0], [0], "o", color=th["ink"], markersize=5, zorder=6)
    axes[0].set_xlim(-3.6, 2.4)
    axes[0].set_ylim(-3.2, 2.2)
    style(axes[0], th, "$\\xi$", "$\\eta$",
          f"Logistic prey, $c = {C_DAMP}$: a damper of ratio $\\zeta = c/2\\omega_0$")
    legend(axes[0], th, loc="lower right")

    for ax, xi1, cyc, title in ((axes[1], XI1_SMALL, cyc_small, "inside the inner region"),
                                (axes[2], XI1_LARGE, cyc_large, "touching the floors")):
        hump = (ZP, ZM, xi1)
        f = field(hump=hump)
        ye = eta_eq(hump)
        rs, T = cyc[0], cyc[1]
        for r0, lab in ((0.15*rs, "from inside"), (3.0*rs, "from outside")):
            x, y = _traj(f, [0.0, ye + r0], 12*T)
            ax.plot(x, y, color=c2 if lab == "from inside" else c0, linewidth=1.0, alpha=0.9,
                    label=lab, zorder=3)
        x, y = _traj(f, [0.0, ye + rs], T)
        ax.plot(x, y, color=c1, linewidth=2.2, label="limit cycle", zorder=5)
        _floors(ax, th)
        ax.axvline(xi1, color=th["ink2"], linewidth=1.0, linestyle=(0, (2, 2)), zorder=2)
        ax.annotate("$\\xi_1$", xy=(xi1, ax.get_ylim()[0]), xytext=(4, 6), textcoords="offset points",
                    fontsize=8, color=th["ink2"])
        ax.plot([0], [ye], "o", color=th["ink"], markersize=5, markerfacecolor=th["surface"], zorder=6)
        style(ax, th, "$\\xi$", "$\\eta$", f"The tent at $\\xi_1 = {xi1}$: cycle {title}")
        legend(ax, th, loc="lower right")
    axes[1].set_xlim(-1.4, 1.4)
    axes[1].set_ylim(-1.4, 1.4)
    axes[2].set_xlim(-4.6, 4.6)
    axes[2].set_ylim(-5.5, 4.6)
    fig.suptitle("Prey density dependence: a damper when it is logistic, a self-excited cycle when it "
                 "is a tent", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-damped")


EXAMPLE_RUNS = (
    # (index into EXAMPLES, start (xi, eta), duration, time unit, unit label, title)
    (1, (0.0, -3.0), 4.0*TWO_PI, 1.0, "time  $\\cdot k_3$", "Lotka's mechanism: product $Y$"),
    (2, (0.0, -3.0), 25.0, 1.0, "years", "SIR with births: infecteds $I$"),
    (3, (0.0, np.log(1.0/(1.0/230e-6*1.0/1.0))), 120e-6, 1e-6, "microseconds",
     "Laser switch-on: photon number $n$"),
)


def fig_examples(th, name):
    """The observable of each field, native equations beside the prototype."""
    fig, axes = newfig(th, 1, 3, figsize=(12.6, 4.0))
    c0, c1, c2 = th["series"]
    for ax, (k, y0, T, unit, ulab, title) in zip(axes, EXAMPLE_RUNS):
        label, (f, eq, mapping) = EXAMPLES[k]
        w0, z, Td, dec = predicted(mapping)
        t, en, ep = native_versus_prototype(f, eq, mapping, list(y0), T)
        _, _, e3 = native_versus_prototype(f, eq, mapping, list(y0), T, s1=S1, k=K)
        ax.plot(t/unit, np.exp(en), color=c0, linewidth=2.0, label="native equations", zorder=3)
        ax.plot(t/unit, np.exp(ep), color=c1, linewidth=1.5, label="prototype, two pieces", zorder=4)
        ax.plot(t/unit, np.exp(e3), color=c2, linewidth=1.5, label="prototype, three pieces", zorder=5)
        ax.axhline(np.exp(-S0), color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
        ax.set_yscale("log")
        style(ax, th, ulab, "observable / its equilibrium value", title)
        ax.annotate(f"$T_0 = {Td/unit:.3g}$, $\\zeta = {z:.3g}$", xy=(0.03, 0.04), xycoords="axes fraction",
                    fontsize=8, color=th["ink2"],
                    bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))
        legend(ax, th, loc="upper right")
    fig.suptitle("The same equation in three fields, integrated as each field writes it, "
                 "against the prototype fitted through the dictionary", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-examples")


def fig_three(th, name, data):
    """The three piece law, both fits, against Lotka-Volterra and the two piece law."""
    fig, axes = newfig(th, 2, 2, figsize=(10.4, 8.0))
    c0, c1, c2 = th["series"]
    T0 = TWO_PI/np.sqrt(ALPHA*GAMMA)
    ax = axes[0, 0]
    s = np.linspace(-3.0, 3.0, 800)
    ax.plot(s, np.expm1(s), color=c0, linewidth=2.0, label="Lotka-Volterra  $e^s - 1$", zorder=3)
    ax.plot(s, phi(s), color=c1, linewidth=1.6, label="two pieces", zorder=4)
    ax.plot(s, phi(s, S0, S1, K), color=c2, linewidth=1.8,
            label=f"three pieces, period matched: $s_1 = {S1:.2f}$, $k = {K:.1f}$", zorder=5)
    ax.plot(s, phi(s, S0, S1_PEAK, K_PEAK), color=c2, linewidth=1.8, linestyle=(0, (3, 2)),
            label=f"three pieces, peak matched: $s_1 = {S1_PEAK:.2f}$, $k = {K_PEAK:.1f}$", zorder=5)
    for x in (-S0, S1_PEAK, S1):
        ax.axvline(x, color=th["ink2"], linewidth=0.9, linestyle=(0, (5, 3)), zorder=2)
    ax.set_ylim(-1.5, 12)
    style(ax, th, "$s$", "rate law", "The law with its knee")
    legend(ax, th, loc="upper left")

    xm = data["xm_dense"]
    three = data["three"]
    ax = axes[0, 1]
    ax.plot(-xm, data["lv_T"]/T0, color=c0, linewidth=2.0, label="Lotka-Volterra", zorder=3)
    ax.plot(-xm, data["pr_T"]/T0, color=c1, linewidth=1.6, label="two pieces", zorder=4)
    ax.plot(-xm, three["period_T"]/T0, color=c2, linewidth=1.8, label="three pieces, period matched", zorder=5)
    ax.plot(-xm, three["peak_T"]/T0, color=c2, linewidth=1.8, linestyle=(0, (3, 2)),
            label="three pieces, peak matched", zorder=5)
    ax.axhline(1.0, color=th["grid"], linewidth=0.8)
    style(ax, th, "prey trough  $-\\xi_{min}$", "period  $T / T_0$", "Period: the knee can undo the floors' delay")
    legend(ax, th, loc="upper left")

    ax = axes[1, 0]
    ax.plot(-xm, data["lv_xmax"], color=c0, linewidth=2.0, label="Lotka-Volterra", zorder=3)
    ax.plot(-xm, data["pr_xmax"], color=c1, linewidth=1.6, label="two pieces", zorder=4)
    ax.plot(-xm, three["period_xmax"], color=c2, linewidth=1.8, label="three pieces, period matched", zorder=5)
    ax.plot(-xm, three["peak_xmax"], color=c2, linewidth=1.8, linestyle=(0, (3, 2)),
            label="three pieces, peak matched", zorder=5)
    style(ax, th, "prey trough  $-\\xi_{min}$", "prey peak  $\\xi_{max}$", "Prey peak: the gap the third piece closes")
    legend(ax, th, loc="upper left")

    ax = axes[1, 1]
    ax.plot(-xm, data["lv_lag"]/T0, color=c0, linewidth=2.0, label="Lotka-Volterra", zorder=3)
    ax.plot(-xm, data["pr_lag"]/T0, color=c1, linewidth=1.6, label="two pieces: exactly $T_0/4$", zorder=4)
    ax.plot(-xm, three["period_lag"]/T0, color=c2, linewidth=1.8, label="three pieces, period matched", zorder=5)
    ax.plot(-xm, three["peak_lag"]/T0, color=c2, linewidth=1.8, linestyle=(0, (3, 2)),
            label="three pieces, peak matched", zorder=5)
    ax.set_ylim(0, 0.3)
    style(ax, th, "prey trough  $-\\xi_{min}$", "lag  $/\\,T_0$", "Predator lag: shrinking once the knee is reached")
    legend(ax, th, loc="upper right")
    fig.suptitle("The third piece: a steeper slope above a knee, fitted two ways, $\\alpha = \\gamma = 1$, $s_0 = 1$",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "lotka-three")


# ----------------------------------------------------------------- driver
def _table(rows, header, fmt):
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join("---" for _ in header) + " |")
    for r in rows:
        print("| " + " | ".join(f.format(v) for f, v in zip(fmt, r)) + " |")


TABLE_TROUGHS = (-0.5, -1.0, -1.5, -2.0, -3.0, -5.0, -8.0)


def checks():
    """Print every number LOTKA.md quotes; return what the figures reuse."""
    out = {}
    w0 = np.sqrt(ALPHA*GAMMA)
    T0 = TWO_PI/w0
    f = field()

    print("## Equilibrium, continuity, conservation (alpha = gamma = 1, s0 = 1)")
    print(f"  field at the origin: {f(0, [0.0, 0.0])}")
    J = jacobian(f, [0.0, 0.0])
    print(f"  Jacobian eigenvalues at the origin: {np.round(np.linalg.eigvals(J), 9)}"
          f"  (formula +-i w0 = +-{w0}i)")
    for pt in ([-S0, 0.7], [0.4, -S0], [-S0, -S0]):
        lo = f(0, [pt[0] - 1e-9 if pt[0] == -S0 else pt[0], pt[1] - 1e-9 if pt[1] == -S0 else pt[1]])
        hi = f(0, [pt[0] + 1e-9 if pt[0] == -S0 else pt[0], pt[1] + 1e-9 if pt[1] == -S0 else pt[1]])
        print(f"  field across the floor at {pt}: jump {np.abs(np.subtract(hi, lo)).max():.1e}")
    for xm in (-1.5, -5.0):
        H = H_from_trough(xm)
        s = solve_ivp(f, (0.0, 3*period(H)), [np.sqrt(2*H/GAMMA), 0.0],
                      t_eval=np.linspace(0, 3*period(H), 600), rtol=RTOL, atol=ATOL,
                      method="DOP853", max_step=0.02)
        Hs = energy(s.y[0], s.y[1])
        print(f"  H drift over three periods from prey trough {xm}: {np.abs(Hs - H).max():.1e} on H = {H:.4f}")

    print("\n## Exact period: assembled formula, piecewise walk, direct integration")
    rows, pts, Tn = [], [], []
    Hc = 0.5*(ALPHA + GAMMA)*S0*S0
    for xm in TABLE_TROUGHS:
        H = H_from_trough(xm)
        Tf, Tw, Ti = period_formula(H), period(H), period_integrated(H)
        regs = sorted({region_name(p[0]) for p in circuit(H)[1]})
        rows.append((xm, H, "".join(str(r) for r in regs), Tf, Tw, Ti, abs(Tf - Ti)))
        pts.append(xm)
        Tn.append(Ti)
    _table(rows, ("prey trough", "H", "regions", "formula", "walk", "integrated", "formula - integrated"),
           ("{:+.1f}", "{:.4f}", "{}", "{:.9f}", "{:.9f}", "{:.9f}", "{:.1e}"))
    print(f"  corner energy Hc = {Hc}; the first floor is touched at H = {0.5*min(ALPHA, GAMMA)*S0**2}")
    out["xm_pts"], out["pr_Tn"] = pts, Tn
    print("  the same with alpha = 2, gamma = 0.5 (w0 = 1 still), s0 = 0.7:")
    rows = []
    for xm in (-0.4, -0.9, -1.5, -3.0, -6.0):
        H = H_from_trough(xm, 2.0, 0.5, 0.7)
        Tf, Ti = period_formula(H, 2.0, 0.5, 0.7), period_integrated(H, 2.0, 0.5, 0.7)
        regs = sorted({region_name(p[0]) for p in circuit(H, 2.0, 0.5, 0.7)[1]})
        rows.append((xm, H, "".join(str(r) for r in regs), Tf, Ti, abs(Tf - Ti)))
    _table(rows, ("prey trough", "H", "regions", "formula", "integrated", "difference"),
           ("{:+.1f}", "{:.4f}", "{}", "{:.9f}", "{:.9f}", "{:.1e}"))

    print("\n## Against Lotka-Volterra at the same prey trough")
    rows, lv_rows = [], {}
    for xm in TABLE_TROUGHS:
        H = H_from_trough(xm)
        _, xmax, ymin, ymax = extremes(H)
        Hl, xmaxl, yminl, ymaxl = lv_extremes(xm)
        Tl, lagl, *_ = lv_run(xm)
        lv_rows[xm] = (Tl, xmaxl, lagl)
        Tp = period(H)
        rows.append((xm, np.exp(xm), Tl/T0, Tp/T0, xmaxl, xmax, lagl/T0, 0.25))
    _table(rows, ("prey trough", "u_min/u*", "T/T0 LV", "T/T0 proto", "prey peak LV", "prey peak proto",
                  "lag/T0 LV", "lag/T0 proto"),
           ("{:+.1f}", "{:.4f}", "{:.4f}", "{:.4f}", "{:+.3f}", "{:+.3f}", "{:.4f}", "{:.4f}"))
    print("  at alpha = gamma both systems have eta_min = xi_min and eta_max = xi_max exactly "
          "(H is the same function of each coordinate); the predator columns would repeat the prey ones")
    print("  prototype lag by integration, prey trough -3 and -8: "
          f"{lag_integrated(H_from_trough(-3.0)):.9f}, {lag_integrated(H_from_trough(-8.0)):.9f}; "
          f"pi/(2 w0) = {0.5*np.pi/w0:.9f}")
    print("  large amplitude law T ~ |xi_min|/alpha, both systems, as T alpha / |xi_min|:")
    for xm in (-5.0, -10.0, -20.0, -40.0):
        Tp = period(H_from_trough(xm))
        Tl = lv_run(xm)[0]
        print(f"    trough {xm:+.0f}: LV {Tl*ALPHA/abs(xm):.4f}  prototype {Tp*ALPHA/abs(xm):.4f}")

    # dense curves for the figure
    xm_dense = -np.linspace(0.05, 8.0, 160)
    out["xm_dense"] = xm_dense
    out["pr_T"] = np.array([period(H_from_trough(x)) for x in xm_dense])
    out["pr_xmax"] = np.array([extremes(H_from_trough(x))[1] for x in xm_dense])
    out["pr_ymin"] = np.array([extremes(H_from_trough(x))[2] for x in xm_dense])
    out["pr_lag"] = np.full_like(xm_dense, 0.5*np.pi/w0)
    lv = [lv_run(x)[:2] for x in xm_dense]
    out["lv_T"] = np.array([a for a, _ in lv])
    out["lv_lag"] = np.array([b for _, b in lv])
    ext = [lv_extremes(x) for x in xm_dense]
    out["lv_xmax"] = np.array([e[1] for e in ext])
    out["lv_ymin"] = np.array([e[2] for e in ext])

    print("\n## The third piece: slope k above a knee s1")
    s1k, kk = fit_third_piece()
    s1p, kp = fit_period_piece(lv_T={xm: lv_rows[xm][0] for xm in PERIOD_TROUGHS})
    print(f"  peak matched at troughs {FIT_TROUGHS}: s1 = {s1k:.6f}, k = {kk:.6f} (constants {S1_PEAK}, {K_PEAK})")
    print(f"  period matched at troughs {PERIOD_TROUGHS}: s1 = {s1p:.6f}, k = {kp:.6f} (constants {S1}, {K})")
    print(f"  exponential's slope e^s at the two knees: {np.exp(S1_PEAK):.3f} and {np.exp(S1):.3f}")
    fk = field(s1=S1, k=K)
    for xm in (-1.5, -5.0):
        H = H_from_trough(xm, s1=S1, k=K)
        Tk = period(H, s1=S1, k=K)
        s = solve_ivp(fk, (0.0, 3*Tk), [Phi_inv(H/GAMMA, +1, S0, S1, K), 0.0],
                      t_eval=np.linspace(0, 3*Tk, 600), rtol=RTOL, atol=ATOL, method="DOP853", max_step=0.01)
        Hs = energy(s.y[0], s.y[1], s1=S1, k=K)
        print(f"  H drift over three periods from prey trough {xm}, period fit: {np.abs(Hs - H).max():.1e} on H = {H:.4f}")
    rows = []
    for xm in TABLE_TROUGHS:
        H = H_from_trough(xm, s1=S1, k=K)
        Tw, Ti = period(H, s1=S1, k=K), period_integrated(H, s1=S1, k=K)
        rows.append((xm, len(circuit(H, s1=S1, k=K)[1]), Tw, Ti, abs(Tw - Ti)))
    _table(rows, ("prey trough", "arcs", "walk", "integrated", "difference"),
           ("{:+.1f}", "{}", "{:.9f}", "{:.9f}", "{:.1e}"))
    print("  (period matched fit, alpha = gamma = 1; the two floors and two knees make nine regions)")
    rows = []
    for xm in (-0.9, -3.0, -6.0):
        args = (2.0, 0.5, 0.7, 0.3, 2.5)
        H = H_from_trough(xm, *args)
        Tw, Ti = period(H, *args), period_integrated(H, *args)
        rows.append((xm, len(circuit(H, *args)[1]), Tw, Ti, abs(Tw - Ti)))
    _table(rows, ("prey trough", "arcs", "walk", "integrated", "difference"),
           ("{:+.1f}", "{}", "{:.9f}", "{:.9f}", "{:.1e}"))
    print("  (alpha = 2, gamma = 0.5, s0 = 0.7, s1 = 0.3, k = 2.5)")
    rows = []
    for xm in TABLE_TROUGHS:
        Tl, xmaxl, lagl = lv_rows[xm]
        row = [xm, Tl/T0, xmaxl, lagl/T0]
        for s1, k in ((S1_PEAK, K_PEAK), (S1, K)):
            H = H_from_trough(xm, s1=s1, k=k)
            row += [period(H, s1=s1, k=k)/T0, extremes(H, s1=s1, k=k)[1], lag_integrated(H, s1=s1, k=k)/T0]
        rows.append(row)
    _table(rows, ("prey trough", "T/T0 LV", "peak LV", "lag/T0 LV", "T/T0 peak fit", "peak", "lag/T0",
                  "T/T0 period fit", "peak", "lag/T0"),
           ("{:+.1f}", "{:.4f}", "{:.3f}", "{:.4f}", "{:.4f}", "{:.3f}", "{:.4f}", "{:.4f}", "{:.3f}", "{:.4f}"))
    print("  large amplitude law, T alpha / |xi_min|, period fit:")
    for xm in (-10.0, -20.0, -40.0):
        Tk = period(H_from_trough(xm, s1=S1, k=K), s1=S1, k=K)
        print(f"    trough {xm:+.0f}: LV {lv_run(xm)[0]*ALPHA/abs(xm):.4f}  three pieces {Tk*ALPHA/abs(xm):.4f}"
              f"  two pieces {period(H_from_trough(xm))*ALPHA/abs(xm):.4f}")
    xs = -np.linspace(0.3, 3.0, 55)
    Ts = np.array([period(H_from_trough(x, s1=S1_PEAK, k=K_PEAK), s1=S1_PEAK, k=K_PEAK) for x in xs])/T0
    print(f"  peak fit: the period dips to {Ts.min():.4f} T0 at a trough of {xs[np.argmin(Ts)]:+.2f}")
    out["three"] = dict(xm_dense=xm_dense)
    for tag, (s1, k) in (("peak", (S1_PEAK, K_PEAK)), ("period", (S1, K))):
        Hs = [H_from_trough(x, s1=s1, k=k) for x in xm_dense]
        out["three"][tag + "_T"] = np.array([period(H, s1=s1, k=k) for H in Hs])
        out["three"][tag + "_xmax"] = np.array([extremes(H, s1=s1, k=k)[1] for H in Hs])
        out["three"][tag + "_lag"] = np.array([lag_integrated(H, s1=s1, k=k) for H in Hs])

    print(f"\n## Logistic prey as damping, c = {C_DAMP}")
    zeta = C_DAMP/(2*w0)
    fd = field(c=C_DAMP)
    print(f"  zeta = c/(2 w0) = {zeta}; Jacobian eigenvalues at the origin "
          f"{np.round(np.linalg.eigvals(jacobian(fd, [0.0, 0.0])), 6)}, formula "
          f"{np.round(-zeta*w0 + 1j*w0*np.sqrt(1 - zeta**2), 6)}")
    _, ys = peaks(fd, [0.0, 0.3], 4, 60.0)
    ratios = ys[1:, 1]/ys[:-1, 1]
    print(f"  inside the inner region, predator peak ratios per cycle: {np.round(ratios, 9)}; "
          f"exp(-2 pi zeta / sqrt(1 - zeta^2)) = {np.exp(-TWO_PI*zeta/np.sqrt(1 - zeta**2)):.9f}")
    n = 6
    _, xp = troughs(fd, [-6.0, 0.0], n, 400.0)
    _, xl = troughs(lv_field(c=C_DAMP), [-6.0, 0.0], n, 400.0)
    xe = trough_map(-6.0, C_DAMP, n)
    _table([(k + 1, xp[k], xe[k], xl[k]) for k in range(n)],
           ("cycle", "prototype, integrated", "prototype, energy map", "Lotka-Volterra, integrated"),
           ("{}", "{:+.4f}", "{:+.4f}", "{:+.4f}"))
    print("  (successive prey troughs after a first trough at -6, as log deviations; "
          f"first cycle loss by the energy map {energy_loss(H_from_trough(-6.0), C_DAMP):.4f} "
          f"on H = {H_from_trough(-6.0):.4f})")

    print(f"\n## The tent: zeta+ = {ZP}, zeta- = {ZM}")
    for xi1 in (XI1_SMALL, XI1_LARGE):
        hump = (ZP, ZM, xi1)
        ye = eta_eq(hump)
        fh = field(hump=hump)
        ev = np.linalg.eigvals(jacobian(fh, [0.0, ye]))
        rs, T, m, xmin, ymin = limit_cycle(hump, 2.0*xi1)
        pred_r, pred_T = README_R*GAMMA*xi1/w0, README_T/w0
        print(f"  xi1 = {xi1}: predator equilibrium eta_eq = {ye:+.4f}, eigenvalues there {np.round(ev, 6)}")
        print(f"    cycle: r* = {rs:.9f}, T = {T:.6f}, multiplier {m:.4f}; "
              f"xi_min = {xmin:+.4f}, eta_min = {ymin:+.4f} on the cycle")
        print(f"    README offset prototype scaled by v0 = gamma xi1 = {GAMMA*xi1}: r* = {pred_r:.9f}, "
              f"T = {pred_T:.6f}, multiplier {README_M}")
        out[f"cyc_{xi1}"] = (rs, T, m)
    print("  amplitude against xi1 (r*/xi1; constant while the cycle stays inside the inner region):")
    for xi1 in (0.1, 0.3, 0.5, 0.7, 1.0, 1.5):
        rs, T, m, xmin, ymin = limit_cycle((ZP, ZM, xi1), 2.0*xi1)
        print(f"    xi1 = {xi1:.1f}: r*/xi1 = {rs/xi1:.6f}, T = {T:.6f}, multiplier {m:.4f}, "
              f"min xi {xmin:+.3f}, min eta {ymin:+.3f}")
    print("  existence: zeta- < 0 < mean, as in the README:")
    for zp, zm in ((0.30, -0.10), (0.05, -0.10), (0.30, 0.10), (0.10, -0.30)):
        rs, T, m, _, _ = limit_cycle((zp, zm, XI1_SMALL), 1.0)
        got = ("grows unbounded" if (not np.isfinite(rs) or rs > 1e5)
               else "decays to the equilibrium" if rs < 1e-6
               else f"limit cycle r* = {rs:.6f}, T = {T:.4f}, multiplier {m:.4f}")
        print(f"    zeta+ = {zp:+.2f}, zeta- = {zm:+.2f}, mean {0.5*(zp + zm):+.3f}: {got}")

    print("\n## The same equation in other fields: dictionary checked by integration")
    rows = []
    for label, (f, eq, mapping) in EXAMPLES:
        w0, z, Td, dec = predicted(mapping)
        Tm, decm = measure(f, eq, mapping)
        rows.append((label, eq[0], eq[1], mapping["alpha"], mapping["gamma"], mapping["c"], w0, z,
                     Td, Tm, dec, decm))
    _table(rows, ("field", "u*", "v*", "alpha", "gamma", "c", "w0", "zeta", "period, formula",
                  "period, measured", "decrement, formula", "decrement, measured"),
           ("{}", "{:.4g}", "{:.4g}", "{:.4g}", "{:.4g}", "{:.4g}", "{:.4g}", "{:.4g}", "{:.6g}",
            "{:.6g}", "{:.6f}", "{:.6f}"))
    print("  (SIR and Goodwin in years, the laser in seconds; period is the damped one 2 pi / w0 sqrt(1 - zeta^2))")
    print("  large excursions, native against prototype, as the peak of the observable over its equilibrium:")
    for k, y0, T, unit, ulab, title in EXAMPLE_RUNS:
        label, (f, eq, mapping) = EXAMPLES[k]
        t, en, ep = native_versus_prototype(f, eq, mapping, list(y0), T)
        _, _, e3 = native_versus_prototype(f, eq, mapping, list(y0), T, s1=S1, k=K)
        i, j, m = np.argmax(en), np.argmax(ep), np.argmax(e3)
        print(f"    {label}: start eta = {y0[1]:.3f}; first peak native e^{en[i]:.3f} at t = {t[i]/unit:.4g}, "
              f"two pieces e^{ep[j]:.3f} at t = {t[j]/unit:.4g}, three pieces e^{e3[m]:.3f} at "
              f"t = {t[m]/unit:.4g} {ulab.split()[0]}")
    return out


def figures_(out=None):
    if out is None:
        out = checks()
    for name, th in THEMES.items():
        fig_law(th, name)
        fig_phase(th, name)
        fig_time(th, name)
        fig_period(th, name, out)
        fig_damped(th, name, out[f"cyc_{XI1_SMALL}"], out[f"cyc_{XI1_LARGE}"])
        fig_examples(th, name)
        fig_three(th, name, out)


if __name__ == "__main__":
    what = sys.argv[1:]
    out = None
    if not what or "checks" in what:
        out = checks()
    if not what or "figures" in what:
        figures_(out)
