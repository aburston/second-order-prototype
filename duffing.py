"""The piecewise Duffing prototype: a switched stiffness, not a switched damping.

Run ``python3 duffing.py`` to print every number ``DUFFING.md`` quotes and
to write its figures into ``figures/`` as ``duffing-*.png`` in both themes.
``python3 duffing.py checks`` prints the numbers only, ``python3 duffing.py
figures`` writes the figures only, and ``quick`` after either shrinks the
basin grids and the drive scans so a run takes a few minutes rather than a
quarter of an hour.
The drive scans and basin grids are cached in ``figures/.duffing-*.json``
and ``figures/.duffing-*.npy``; delete them, or pass ``fresh``, to recompute.

The model
---------

Every prototype in ``README.md`` and ``THREELEVEL.md`` switches the
*damping* ratio and keeps the stiffness linear, so none of them can have
more than one equilibrium. This one does the opposite: the damping is one
ordinary ``zeta`` everywhere and the *stiffness* is switched at a
displacement threshold, negative inside a band and positive outside it,
with the restoring force kept continuous::

    x'' + 2 zeta wn x' + g(x) = A cos(Om t)

    g(x) = -kappa wn^2 x                 |x| < x0        (the saddle band)
         =  wn^2 (x - sign(x) xe)        |x| > x0        (the wells)

    xe = (1 + kappa) x0                  so that g is continuous at +-x0

That is Duffing's double well ``-alpha x + beta x^3`` with the cubic
replaced by two straight lines: a saddle at the origin, a well either side
of it at ``+-xe``, and small oscillations in a well that are *exactly* the
linear prototype at ``wn`` and ``zeta``. Four parameters, each read off one
measurement:

``wn``      the frequency of small oscillations about a well
``zeta``    their decrement
``kappa``   the saddle stiffness ratio: the growth rate at the saddle is
            ``sqrt(kappa) wn``, and the well depth is ``kappa wn^2 x0 xe / 2``
``x0``      the half-width of the saddle band, which sets the amplitude
            scale and puts the wells at ``xe = (1 + kappa) x0``

Two variants share the force law and differ only in whether ``x`` wraps:

``periodic=False``  the **beam**: a buckled Euler strut has two buckled
                    states (the wells) and the straight configuration
                    between them (the saddle), and the stiffness keeps
                    rising beyond the wells.
``periodic=True``   the **pendulum**: ``x`` is the angle measured from the
                    *inverted* position, so the saddle is the top and the
                    hanging position appears twice, at ``x = +-pi``. Those
                    two wells are one point of the circle, so ``g`` is
                    tiled with period ``2 xe = 2 pi`` and the state lives
                    on a cylinder. Full rotations are then ordinary
                    trajectories that advance ``x`` by ``2 pi``.

The pendulum fixes ``xe = pi`` by geometry, so it has one shape parameter,
``kappa``, and the beam has two, ``kappa`` and ``x0``. Everything below is
in units ``wn = 1`` unless a section says otherwise; the scaling rule of
``scaling.py`` moves the results to any frequency and amplitude.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.colors import ListedColormap, to_rgb
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import ellipk

from figures import THEMES, style, legend, save, newfig, OUT

WN = 1.0

#: Shape parameters that reproduce named smooth systems, see ``fits()``.
KAPPA_PEND_SLOPE = 1.0                   # matches sin's slope at both ends
KAPPA_PEND_DEPTH = 4.0/(np.pi**2 - 4.0)  # matches the well depth 2 wn^2
KAPPA_DUFF_SLOPE = 0.5                   # Duffing's saddle rate is wn/sqrt2
KAPPA_DUFF_DEPTH = 1.0/3.0               # Duffing's well depth alpha^2/4beta

#: Damping used for the basin figures and the capture table.
ZETA_BEAM, ZETA_PEND = 0.10, 0.05

#: Holmes' forced Duffing, ``x'' + 0.25 x' - x + x^3 = A cos t``, and the
#: forced pendulum of Baker and Gollub, ``th'' + th'/2 + sin th = A cos(2t/3)``.
HOLMES_DELTA, HOLMES_OM = 0.25, 1.0
BG_Q, BG_OM = 2.0, 2.0/3.0
BEAM_AMPS = np.round(np.arange(0.10, 0.52, 0.02), 3)
PEND_AMPS = np.round(np.arange(0.60, 1.62, 0.05), 3)
BEAM_SHOW, PEND_SHOW = 0.40, 1.25

RTOL, ATOL = 1e-11, 1e-13
CACHE = os.path.join(OUT, ".duffing-{}.json")


# ------------------------------------------------------------------ model
def wells(kappa, x0):
    """Well position ``xe = (1 + kappa) x0``, forced by continuity at ``x0``."""
    return (1.0 + kappa)*x0


def depth(kappa, x0, wn=WN):
    """Well depth ``V(0) - V(xe) = kappa (1 + kappa) wn^2 x0^2 / 2``."""
    return 0.5*kappa*(1.0 + kappa)*wn**2*x0**2


def reduce(x, xe, periodic):
    """Bring ``x`` into the fundamental domain ``[-xe, xe)`` when periodic."""
    if not periodic:
        return x
    return (x + xe) % (2.0*xe) - xe


def force(x, kappa, x0, periodic=False, wn=WN):
    """Restoring acceleration ``g(x)``: saddle at 0, wells at ``+-xe``."""
    xe = wells(kappa, x0)
    x = reduce(np.asarray(x, float), xe, periodic)
    return np.where(np.abs(x) < x0, -kappa*wn**2*x,
                    wn**2*(x - np.sign(x)*xe))


def stiffness(x, kappa, x0, periodic=False, wn=WN):
    """``dg/dx``: ``-kappa wn^2`` in the band, ``+wn^2`` in the wells."""
    xe = wells(kappa, x0)
    x = reduce(np.asarray(x, float), xe, periodic)
    return np.where(np.abs(x) < x0, -kappa*wn**2, wn**2)


def potential(x, kappa, x0, periodic=False, wn=WN):
    """``V(x)`` with ``V(0) = 0`` at the saddle, so the wells sit at ``-depth``."""
    xe = wells(kappa, x0)
    x = reduce(np.asarray(x, float), xe, periodic)
    return np.where(np.abs(x) < x0, -0.5*kappa*wn**2*x**2,
                    0.5*wn**2*(np.abs(x) - xe)**2 - depth(kappa, x0, wn))


def energy(x, v, kappa, x0, periodic=False, wn=WN):
    """Total energy ``v^2/2 + V(x)``; zero on the separatrix."""
    return 0.5*np.asarray(v, float)**2 + potential(x, kappa, x0, periodic, wn)


def field(zeta, kappa, x0, periodic=False, amp=0.0, om=1.0, wn=WN):
    """Vector field ``f(t, y)`` in ``y = [x, xdot]`` for ``solve_ivp``.

    The force is continuous with a corner at ``+-x0``, exactly as the
    README's deadzone model was, so the field is Lipschitz: solutions are
    unique, nothing slides, and a general purpose integrator crosses the
    corners without chatter. Only the Jacobian jumps there.
    """
    xe = wells(kappa, x0)

    def f(t, y):
        x = reduce(y[0], xe, periodic)
        g = -kappa*wn**2*x if abs(x) < x0 else wn**2*(x - np.copysign(xe, x))
        return [y[1], -2.0*zeta*wn*y[1] - g + amp*np.cos(om*t)]
    return f


def variational(zeta, kappa, x0, periodic=False, amp=0.0, om=1.0, wn=WN):
    """The field with its 2x2 fundamental matrix, ``z = [x, xdot, J.flat]``.

    Because the field is continuous the saltation matrix at a corner is the
    identity, so integrating the variational equation through the corner
    with the jumping Jacobian is the whole story: no matching is needed.
    """
    xe = wells(kappa, x0)
    d = -2.0*zeta*wn

    def f(t, z):
        x = reduce(z[0], xe, periodic)
        if abs(x) < x0:
            g, k = -kappa*wn**2*x, -kappa*wn**2
        else:
            g, k = wn**2*(x - np.copysign(xe, x)), wn**2
        a, b, c, e = z[2], z[3], z[4], z[5]
        return [z[1], d*z[1] - g + amp*np.cos(om*t),
                c, e, -k*a + d*c, -k*b + d*e]
    return f


def eigenvalues(zeta, kappa, wn=WN):
    """Eigenvalues at the saddle and at a well, from the two linear pieces."""
    saddle = wn*(-zeta + np.sqrt(zeta**2 + kappa)), wn*(-zeta - np.sqrt(zeta**2 + kappa))
    well = wn*(-zeta + np.sqrt(complex(zeta**2 - 1.0))), wn*(-zeta - np.sqrt(complex(zeta**2 - 1.0)))
    return saddle, well


def check_equilibria(zeta, kappa, x0, wn=WN):
    """Evaluate the field at the claimed equilibria and difference the Jacobian.

    The README's working notes record a sign error that survived a reading
    and was caught only by evaluating the field at the claimed point, so
    that is done here rather than trusting the algebra.
    """
    f = field(zeta, kappa, x0, wn=wn)
    xe = wells(kappa, x0)
    rows = []
    for name, xs in (("saddle", 0.0), ("well +", xe), ("well -", -xe)):
        fx = f(0.0, [xs, 0.0])
        h = 1e-6
        J = np.array([[(f(0, [xs + h, 0])[i] - f(0, [xs - h, 0])[i])/(2*h),
                       (f(0, [xs, h])[i] - f(0, [xs, -h])[i])/(2*h)]
                      for i in range(2)])
        rows.append((name, xs, max(abs(fx[0]), abs(fx[1])), np.linalg.eigvals(J)))
    return rows


# ---------------------------------------------------- exact periods, undamped
def _radius(E, kappa, x0, wn):
    """Amplitude of the well arc about ``xe`` at energy ``E``."""
    return np.sqrt(2.0*(E + depth(kappa, x0, wn)))/wn


def period_well(E, kappa, x0, wn=WN):
    """Beam, one well, ``-depth < E < 0``.

    Below the corner the orbit never leaves the well and the period is the
    linear ``2 pi / wn``. Above it the orbit spends ``2 t_w`` in the well,
    an arc of the linear oscillator from the far turning point to the
    corner, and ``2 t_b`` in the band, a hyperbolic arc ``x = x_- cosh(mu t)``
    from the corner to the near turning point ``x_-`` and back.
    """
    R = _radius(E, kappa, x0, wn)
    if R <= kappa*x0:
        return 2.0*np.pi/wn
    mu = np.sqrt(kappa)*wn
    tw = np.arccos(-kappa*x0/R)/wn
    xm = np.sqrt(-2.0*E/kappa)/wn
    tb = np.arccosh(x0/xm)/mu
    return 2.0*(tw + tb)


def period_cross(E, kappa, x0, wn=WN):
    """Beam, both wells, ``E > 0``: four well arcs and four band half-crossings."""
    R = _radius(E, kappa, x0, wn)
    mu = np.sqrt(kappa)*wn
    tw = np.arccos(-kappa*x0/R)/wn
    tb = np.arcsinh(mu*x0/np.sqrt(2.0*E))/mu
    return 4.0*(tw + tb)


def period_libration(E, kappa, x0, wn=WN):
    """Pendulum swinging in its well, ``-depth < E < 0``.

    The well is bounded by a saddle band on *both* sides, so a large swing
    enters a band at each end and the orbit is symmetric about the well:
    four quarter arcs of the linear oscillator plus four band arcs.
    """
    R = _radius(E, kappa, x0, wn)
    if R <= kappa*x0:
        return 2.0*np.pi/wn
    mu = np.sqrt(kappa)*wn
    tw = np.arcsin(kappa*x0/R)/wn
    xm = np.sqrt(-2.0*E/kappa)/wn
    tb = np.arccosh(x0/xm)/mu
    return 4.0*(tw + tb)


def period_rotation(E, kappa, x0, wn=WN):
    """Pendulum going over the top, ``E > 0``: time to advance by ``2 pi``."""
    R = _radius(E, kappa, x0, wn)
    mu = np.sqrt(kappa)*wn
    tw = 2.0*np.arcsin(kappa*x0/R)/wn
    tb = 2.0*np.arcsinh(mu*x0/np.sqrt(2.0*E))/mu
    return tw + tb


def period(E, kappa, x0, periodic, wn=WN):
    """The right closed form for the orbit type at energy ``E``."""
    if periodic:
        return period_rotation(E, kappa, x0, wn) if E > 0 else period_libration(E, kappa, x0, wn)
    return period_cross(E, kappa, x0, wn) if E > 0 else period_well(E, kappa, x0, wn)


def period_integrated(E, kappa, x0, periodic, wn=WN):
    """Period by direct integration from the well bottom, for checking."""
    xe = wells(kappa, x0)
    v0 = np.sqrt(2.0*(E + depth(kappa, x0, wn)))
    f = field(0.0, kappa, x0, periodic, wn=wn)
    if periodic and E > 0:
        ev = lambda t, y: y[0] - 3.0*xe
        ev.terminal, ev.direction = True, 1
        s = solve_ivp(f, (0, 1e4), [xe, v0], events=ev, rtol=RTOL, atol=ATOL)
        return s.t_events[0][0]
    ev = lambda t, y: y[0] - xe
    ev.direction = 1
    s = solve_ivp(f, (0, 1e4), [xe, v0], events=ev, rtol=RTOL, atol=ATOL)
    te = s.t_events[0]
    return te[te > 1e-6][0]


def action(E, kappa, x0, wn=WN):
    """``oint xdot dx`` over one rotation of the pendulum, in closed form.

    The band contributes ``int sqrt(2E + kappa wn^2 x^2)`` over ``|x| < x0``
    and the well ``int sqrt(2(E + depth) - wn^2 u^2)`` over ``|u| < kappa x0``,
    both elementary. Multiplied by ``2 zeta wn`` this is the energy the
    damping removes per revolution to first order in ``zeta``.
    """
    p2, q = 2.0*E, np.sqrt(kappa)*wn
    a = x0
    band = a*np.sqrt(p2 + q*q*a*a) + (p2/q)*np.arcsinh(q*a/np.sqrt(p2))
    r2, a = 2.0*(E + depth(kappa, x0, wn)), kappa*x0
    well = a*np.sqrt(r2 - wn*wn*a*a) + (r2/wn)*np.arcsin(wn*a/np.sqrt(r2))
    return band + well


def action_quad(E, kappa, x0, wn=WN):
    """The same integral by quadrature, to check the closed form."""
    xe = wells(kappa, x0)
    return quad(lambda x: np.sqrt(2.0*(E - potential(x, kappa, x0, True, wn))),
                -xe, xe, limit=200)[0]


# ------------------------------------------------------- the smooth targets
def pendulum_period(amplitude_rad, wn=WN):
    """Exact period of ``th'' + wn^2 sin th = 0`` swinging to ``+-amplitude``."""
    return 4.0*ellipk(np.sin(amplitude_rad/2.0)**2)/wn


def pendulum_rotation_period(v, wn=WN):
    """Exact time per revolution with speed ``v`` at the bottom, ``v > 2 wn``."""
    k2 = 4.0*wn*wn/(v*v)
    return 4.0*ellipk(k2)/v


def duffing_field(delta=0.0, amp=0.0, om=1.0):
    """``x'' + delta x' - x + x^3 = amp cos(om t)``, the double well at alpha = beta = 1."""
    def f(t, y):
        return [y[1], -delta*y[1] + y[0] - y[0]**3 + amp*np.cos(om*t)]
    return f


def duffing_variational(delta, amp, om):
    def f(t, z):
        x, v = z[0], z[1]
        k = 1.0 - 3.0*x*x
        a, b, c, e = z[2], z[3], z[4], z[5]
        return [v, -delta*v + x - x**3 + amp*np.cos(om*t),
                c, e, k*a - delta*c, k*b - delta*e]
    return f


def pendulum_field(zeta, amp=0.0, om=1.0, wn=WN):
    """The true pendulum in the same coordinates, ``x`` from the inverted position."""
    def f(t, y):
        return [y[1], -2.0*zeta*wn*y[1] + wn*wn*np.sin(y[0]) + amp*np.cos(om*t)]
    return f


def pendulum_variational(zeta, amp, om, wn=WN):
    def f(t, z):
        x, v = z[0], z[1]
        k, d = wn*wn*np.cos(x), -2.0*zeta*wn
        a, b, c, e = z[2], z[3], z[4], z[5]
        return [v, d*v + wn*wn*np.sin(x) + amp*np.cos(om*t),
                c, e, k*a + d*c, k*b + d*e]
    return f


def duffing_well_period(a):
    """Duffing's in-well period at inner amplitude ``a`` (turning point at ``1 - a``)."""
    f = duffing_field()
    V = lambda x: -0.5*x*x + 0.25*x**4
    E = V(1.0 - a)
    v0 = np.sqrt(2.0*(E - V(1.0)))
    ev = lambda t, y: y[0] - 1.0
    ev.direction = 1
    s = solve_ivp(f, (0, 1e3), [1.0, v0], events=ev, rtol=RTOL, atol=ATOL)
    te = s.t_events[0]
    return te[te > 1e-6][0]


def fits():
    """The shape parameter ``kappa`` that matches each smooth target.

    With ``wn`` set by the well and ``xe`` by the well position, one number
    is left, and it can match the saddle's growth rate *or* the well depth
    (equivalently the escape speed from the well bottom), not both.
    """
    return {
        "pendulum": dict(slope=KAPPA_PEND_SLOPE, depth=KAPPA_PEND_DEPTH, xe=np.pi),
        "duffing": dict(slope=KAPPA_DUFF_SLOPE, depth=KAPPA_DUFF_DEPTH, xe=1.0),
    }


# ---------------------------------------------------------------- damped
def _settle_event(kappa, x0, periodic, wn):
    """Terminal event: the energy drops through zero, trapping the orbit."""
    def ev(t, y):
        return energy(y[0], y[1], kappa, x0, periodic, wn)
    ev.terminal, ev.direction = True, -1
    return ev


def settle(y0, zeta, kappa, x0, periodic=False, wn=WN, tmax=None):
    """Integrate until trapped and return ``(x, xdot, t)`` at trapping.

    Once ``E < 0`` the orbit cannot reach the saddle again and the damping
    only lowers ``E`` further, so the well it is in is the well it ends in.
    """
    if energy(y0[0], y0[1], kappa, x0, periodic, wn) < 0:
        return y0[0], y0[1], 0.0
    f = field(zeta, kappa, x0, periodic, wn=wn)
    tmax = tmax or 40.0/(zeta*wn)
    s = solve_ivp(f, (0.0, tmax), y0, events=_settle_event(kappa, x0, periodic, wn),
                  rtol=1e-9, atol=1e-11)
    if s.t_events[0].size:
        y = s.y_events[0][0]
        return y[0], y[1], s.t_events[0][0]
    return s.y[0, -1], s.y[1, -1], tmax


def cell(x, xe):
    """Index of the cell between saddles that ``x`` lies in, saddles at ``2 xe m``."""
    return int(np.floor(x/(2.0*xe)))


def _basin_job(args):
    y0, zeta, kappa, x0, periodic = args
    x, v, t = settle(list(y0), zeta, kappa, x0, periodic)
    xe = wells(kappa, x0)
    if periodic:
        # net number of times over the top: saddle crossings forwards less backwards
        return cell(x, xe) - cell(y0[0], xe)
    return 1 if x > 0 else -1


def basins(zeta, kappa, x0, periodic, xs, vs, workers=4, cache=None):
    """Which well each initial condition on the grid settles into."""
    path = cache and cache.replace(".json", ".npy")
    if path and os.path.exists(path):
        L = np.load(path)
        if L.shape == (len(vs), len(xs)):
            return L
    jobs = [((x, v), zeta, kappa, x0, periodic) for v in vs for x in xs]
    with Pool(workers) as p:
        out = p.map(_basin_job, jobs, chunksize=64)
    L = np.array(out).reshape(len(vs), len(xs))
    if path:
        np.save(path, L)
    return L


def stable_manifold(zeta, kappa, x0, wn=WN, T=40.0, eps=1e-6):
    """The saddle's stable manifold, integrated backwards from its eigenvector."""
    (lp, lm), _ = eigenvalues(zeta, kappa, wn)
    ev = np.array([1.0, lm])/np.hypot(1.0, lm)
    f = field(zeta, kappa, x0, wn=wn)
    fb = lambda t, y: [-c for c in f(t, y)]
    out = []
    for sgn in (1.0, -1.0):
        s = solve_ivp(fb, (0.0, T), sgn*eps*ev, rtol=RTOL, atol=ATOL,
                      t_eval=np.linspace(0, T, 4000))
        out.append((s.y[0], s.y[1]))
    return out


def turns(v, zeta, kappa, x0, wn=WN):
    """Full revolutions before capture, starting at the bottom with speed ``v``."""
    xe = wells(kappa, x0)
    x, _, _ = settle([xe, v], zeta, kappa, x0, True, wn)
    return cell(x, xe) - cell(xe, xe)


def capture_speeds(zeta, kappa, x0, n_max=5, wn=WN):
    """Smallest bottom speed ``v_n`` that completes ``n`` revolutions, by bisection."""
    out = []
    lo = np.sqrt(2.0*depth(kappa, x0, wn))*0.5
    hi = lo
    for n in range(1, n_max + 1):
        while turns(hi, zeta, kappa, x0, wn) < n:
            hi *= 1.3
        a, b = lo, hi
        for _ in range(32):
            m = 0.5*(a + b)
            if turns(m, zeta, kappa, x0, wn) < n:
                a = m
            else:
                b = m
        out.append(b)
        lo = b
    return out


def capture_speeds_averaged(zeta, kappa, x0, n_max=5, wn=WN):
    """The same thresholds from the energy map ``E -> E - 2 zeta wn J(E)``.

    First order in ``zeta``: each revolution removes the work the damping
    does over the undamped orbit at that energy. Exact in the limit
    ``zeta -> 0`` and a few per cent out at ``zeta = 0.05``.
    """
    d = depth(kappa, x0, wn)

    def n_of(v):
        E, n = 0.5*v*v - d, 0
        # a turn completes if the energy is still positive at the top,
        # half a revolution's dissipation later; the rest is lost coming down
        while E > 0 and E - zeta*wn*action(E, kappa, x0, wn) > 0:
            E -= 2.0*zeta*wn*action(E, kappa, x0, wn)
            n += 1
        return n

    out, lo = [], np.sqrt(2.0*d)
    for n in range(1, n_max + 1):
        hi = lo
        while n_of(hi) < n:
            hi *= 1.3
        out.append(brentq(lambda v: n_of(v) - n + 0.5, lo, hi, xtol=1e-10))
        lo = out[-1]
    return out


# ---------------------------------------------------------------- forced
def lyapunov(fvar, om, y0, n_skip=100, n=400):
    """Largest Lyapunov exponent from the fundamental matrix, once per drive period."""
    td = 2.0*np.pi/om
    y = np.array(y0, float)
    for i in range(n_skip):
        s = solve_ivp(fvar, (i*td, (i + 1)*td), [*y, 1, 0, 0, 1],
                      rtol=1e-9, atol=1e-11, method="DOP853")
        y = s.y[:2, -1]
    acc, w = 0.0, np.array([1.0, 0.0])
    for i in range(n):
        t0 = (n_skip + i)*td
        s = solve_ivp(fvar, (t0, t0 + td), [*y, 1, 0, 0, 1],
                      rtol=1e-9, atol=1e-11, method="DOP853")
        y, J = s.y[:2, -1], s.y[2:, -1].reshape(2, 2)
        w = J @ w
        nw = np.linalg.norm(w)
        acc += np.log(nw)
        w /= nw
    return acc/(n*td)


def forced_system(kind, amp):
    """``(variational field, drive frequency, start)`` for each system compared.

    The beam family is Holmes' Duffing at ``alpha = beta = 1``: ``wn = sqrt 2``,
    ``xe = 1``, ``2 zeta wn = 0.25``. The pendulum family is Baker and
    Gollub's driven pendulum: ``wn = 1``, ``2 zeta = 1/q``, ``xe = pi``.
    """
    if kind == "duffing":
        return duffing_variational(HOLMES_DELTA, amp, HOLMES_OM), HOLMES_OM, [1.0, 0.0]
    if kind in ("beam-slope", "beam-depth"):
        k = KAPPA_DUFF_SLOPE if kind == "beam-slope" else KAPPA_DUFF_DEPTH
        wn = np.sqrt(2.0)
        return (variational(HOLMES_DELTA/(2*wn), k, 1.0/(1 + k), False, amp, HOLMES_OM, wn),
                HOLMES_OM, [1.0, 0.0])
    if kind == "pendulum":
        return pendulum_variational(0.5/BG_Q, amp, BG_OM), BG_OM, [np.pi, 0.0]
    if kind in ("pend-slope", "pend-depth"):
        k = KAPPA_PEND_SLOPE if kind == "pend-slope" else KAPPA_PEND_DEPTH
        return (variational(0.5/BG_Q, k, np.pi/(1 + k), True, amp, BG_OM),
                BG_OM, [np.pi, 0.0])
    raise ValueError(kind)


def _lyap_job(args):
    kind, amp = args
    fvar, om, y0 = forced_system(kind, amp)
    return kind, float(amp), lyapunov(fvar, om, y0)


def scan(kinds, amps, cache, workers=4, fresh=False):
    """Largest Lyapunov exponent of each system across drive amplitude."""
    if not fresh and os.path.exists(cache):
        with open(cache) as fh:
            got = json.load(fh)
        if all(k in got and len(got[k]) == len(amps) for k in kinds):
            return {k: dict(zip(map(float, amps), got[k])) for k in kinds}
    jobs = [(k, a) for k in kinds for a in amps]
    with Pool(workers) as p:
        res = p.map(_lyap_job, jobs)
    out = {k: {} for k in kinds}
    for k, a, lam in res:
        out[k][a] = lam
    os.makedirs(OUT, exist_ok=True)
    with open(cache, "w") as fh:
        json.dump({k: [out[k][float(a)] for a in amps] for k in kinds}, fh)
    return out


def strobe(kind, amp, n_skip=300, n=4000):
    """Stroboscopic samples ``(x, xdot)`` once per drive period, transient dropped."""
    fvar, om, y0 = forced_system(kind, amp)
    f = lambda t, z: fvar(t, [*z, 1, 0, 0, 1])[:2]
    td = 2.0*np.pi/om
    y = np.array(y0, float)
    Y = np.empty((n, 2))
    for i in range(n_skip + n):
        s = solve_ivp(f, (i*td, (i + 1)*td), y, rtol=1e-9, atol=1e-11, method="DOP853")
        y = s.y[:, -1]
        if i >= n_skip:
            Y[i - n_skip] = y
    return Y


# ---------------------------------------------------------------- figures
def _ink_line(ax, th, x, y, **kw):
    ax.plot(x, y, color=th["ink2"], linewidth=1.0, linestyle=(0, (4, 3)), zorder=5, **kw)


def fig_force(th, name):
    """Force and potential: the piecewise law against Duffing and against sin."""
    fig, axes = newfig(th, 2, 2, figsize=(9.6, 6.4))
    c0, c1, c2 = th["series"]

    # beam against Duffing, alpha = beta = 1: wn = sqrt2, xe = 1
    wn = np.sqrt(2.0)
    x = np.linspace(-1.9, 1.9, 800)
    axes[0, 0].plot(x, x**3 - x, color=c0, linewidth=2.0, label="Duffing  $x^3 - x$", zorder=3)
    axes[1, 0].plot(x, -x*x/2 + x**4/4, color=c0, linewidth=2.0, label="Duffing", zorder=3)
    for k, c, lab in ((KAPPA_DUFF_SLOPE, c1, "piecewise, $\\kappa = 1/2$ (slopes)"),
                      (KAPPA_DUFF_DEPTH, c2, "piecewise, $\\kappa = 1/3$ (depth)")):
        x0 = 1.0/(1 + k)
        axes[0, 0].plot(x, force(x, k, x0, wn=wn), color=c, linewidth=1.8, label=lab, zorder=4)
        axes[1, 0].plot(x, potential(x, k, x0, wn=wn), color=c, linewidth=1.8, label=lab, zorder=4)
    for ax in axes[:, 0]:
        ax.axvline(0, color=th["grid"], linewidth=0.8)
    style(axes[0, 0], th, "$x$", "$g(x)$", "Beam: restoring force")
    style(axes[1, 0], th, "$x$", "$V(x)$", "Beam: potential")

    # pendulum against sin, over two periods
    x = np.linspace(-2*np.pi, 2*np.pi, 1200)
    axes[0, 1].plot(x, -np.sin(x), color=c0, linewidth=2.0, label="pendulum  $-\\sin x$", zorder=3)
    axes[1, 1].plot(x, np.cos(x) - 1, color=c0, linewidth=2.0, label="pendulum", zorder=3)
    for k, c, lab in ((KAPPA_PEND_SLOPE, c1, "piecewise, $\\kappa = 1$ (slopes)"),
                      (KAPPA_PEND_DEPTH, c2, "piecewise, $\\kappa = 0.68$ (depth)")):
        x0 = np.pi/(1 + k)
        axes[0, 1].plot(x, force(x, k, x0, True), color=c, linewidth=1.8, label=lab, zorder=4)
        axes[1, 1].plot(x, potential(x, k, x0, True), color=c, linewidth=1.8, label=lab, zorder=4)
    for ax in axes[:, 1]:
        for xs in (-np.pi, np.pi):
            ax.axvline(xs, color=th["grid"], linewidth=0.8)
        ax.set_xticks([-2*np.pi, -np.pi, 0, np.pi, 2*np.pi])
        ax.set_xticklabels(["$-2\\pi$", "$-\\pi$", "0", "$\\pi$", "$2\\pi$"])
    style(axes[0, 1], th, "$x$  (angle from the inverted position)", "$g(x)$",
          "Pendulum: restoring force, tiled")
    style(axes[1, 1], th, "$x$", "$V(x)$",
          "Pendulum: potential — the wells at $\\pm\\pi$ are one point of the circle")
    axes[0, 0].annotate("saddle", xy=(0, 0), xytext=(0, 14), textcoords="offset points",
                        ha="center", fontsize=8, color=th["ink2"])
    axes[1, 0].annotate("wells at $\\pm x_e$", xy=(1, -0.25), xytext=(20, -10),
                        textcoords="offset points", fontsize=8, color=th["ink2"])
    legend(axes[0, 0], th, loc="lower right")
    legend(axes[0, 1], th, loc="lower right")
    legend(axes[1, 0], th, loc="upper center")
    legend(axes[1, 1], th, loc="lower center")
    fig.suptitle("The piecewise Duffing force law: a saddle band of slope $-\\kappa\\omega_n^2$ "
                 "between wells of slope $+\\omega_n^2$", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "duffing-force")


def _portrait(ax, th, kappa, x0, periodic, energies, xlim, cross_label):
    """Undamped phase portrait at the listed energies plus the separatrix."""
    xe = wells(kappa, x0)
    f = field(0.0, kappa, x0, periodic)
    c0, c1, c2 = th["series"]
    d = depth(kappa, x0)
    for E in energies:
        v0 = np.sqrt(2.0*(E + d))
        T = period(E, kappa, x0, periodic)
        for sgn in ((1,) if periodic else (1, -1)):
            s = solve_ivp(f, (0, T), [sgn*xe, v0], t_eval=np.linspace(0, T, 1200),
                          rtol=RTOL, atol=ATOL)
            x = s.y[0]
            col = c1 if E < 0 else c2
            if periodic and E > 0:
                # rotation: draw both senses, unwrapped across the window
                for shift in (-4*np.pi, -2*np.pi, 0, 2*np.pi):
                    ax.plot(x + shift, s.y[1], color=col, linewidth=1.4, zorder=3)
                    ax.plot(-(x + shift), -s.y[1], color=col, linewidth=1.4, zorder=3)
            elif periodic:
                for shift in (-2*np.pi, 0.0):
                    ax.plot(x + shift, s.y[1], color=col, linewidth=1.4, zorder=3)
            else:
                ax.plot(x, s.y[1], color=col, linewidth=1.4, zorder=3)
    # separatrix: E = 0, drawn from the potential
    xs = np.linspace(xlim[0], xlim[1], 3000)
    V = potential(xs, kappa, x0, periodic)
    vs = np.where(V <= 0, np.sqrt(np.maximum(-2.0*V, 0.0)), np.nan)
    ax.plot(xs, vs, color=c0, linewidth=2.4, zorder=4, label="separatrix  $E = 0$")
    ax.plot(xs, -vs, color=c0, linewidth=2.4, zorder=4)
    ax.plot([0], [0], "o", color=th["ink"], markersize=6, zorder=6)
    ax.plot([xe, -xe], [0, 0], "o", color=th["ink"], markersize=6,
            markerfacecolor=th["surface"], zorder=6)
    for xs_ in (-x0, x0):
        ax.axvline(xs_, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    ax.plot([], [], color=c1, linewidth=1.4, label="in a well  $E < 0$")
    ax.plot([], [], color=c2, linewidth=1.4, label=cross_label)
    ax.set_xlim(*xlim)


def fig_phase(th, name):
    """Undamped portraits of the beam and the pendulum at the same kappa and xe."""
    k, x0 = KAPPA_PEND_SLOPE, np.pi/2
    fig, axes = newfig(th, 1, 2, figsize=(10.4, 4.6))
    d = depth(k, x0)
    _portrait(axes[0], th, k, x0, False, [-0.8*d, -0.3*d, -0.05*d, 0.4*d, 1.5*d],
              (-2.05*np.pi, 2.05*np.pi), "across both wells  $E > 0$")
    _portrait(axes[1], th, k, x0, True, [-0.8*d, -0.3*d, -0.05*d, 0.4*d, 1.5*d],
              (-2.05*np.pi, 2.05*np.pi), "rotation  $E > 0$")
    for ax in axes:
        ax.set_xticks([-2*np.pi, -np.pi, 0, np.pi, 2*np.pi])
        ax.set_xticklabels(["$-2\\pi$", "$-\\pi$", "0", "$\\pi$", "$2\\pi$"])
        ax.set_ylim(-5.2, 5.2)
    axes[1].axvspan(-np.pi, np.pi, color=th["grid"], alpha=0.35, zorder=1)
    axes[1].annotate("one turn of the circle", xy=(0, -4.7), ha="center", fontsize=8,
                     color=th["ink2"])
    axes[0].annotate("saddle band $|x| < x_0$", xy=(0, 4.7), ha="center", fontsize=8,
                     color=th["ink2"], bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"],
                                                 ec="none"), zorder=7)
    style(axes[0], th, "$x$", "$\\dot{x}$",
          "Beam: $\\kappa = 1$, $x_e = \\pi$ — two wells, one saddle")
    style(axes[1], th, "$x$  (angle from the inverted position)", "$\\dot{x}$",
          "Pendulum: the same law tiled — the wells at $\\pm\\pi$ are one point")
    legend(axes[0], th, loc="lower left")
    legend(axes[1], th, loc="lower left")
    fig.tight_layout()
    save(fig, name, "duffing-phase")


def fig_period(th, name, checks):
    """Period against energy and against amplitude, exact and against the targets."""
    fig, axes = newfig(th, 1, 3, figsize=(12.6, 4.0))
    c0, c1, c2 = th["series"]

    # 1. beam: period against energy, both branches, the log divergence
    k, x0 = KAPPA_DUFF_SLOPE, 1.0/(1 + KAPPA_DUFF_SLOPE)
    d = depth(k, x0)
    Ew = -d*np.logspace(0, -6, 400)
    Ec = d*np.logspace(-6, 1, 400)
    axes[0].plot(Ew/d, [period_well(E, k, x0)/(2*np.pi) for E in Ew], color=c1, linewidth=2.0,
                 label="in one well", zorder=3)
    axes[0].plot(Ec/d, [period_cross(E, k, x0)/(2*np.pi) for E in Ec], color=c2, linewidth=2.0,
                 label="across both wells", zorder=3)
    ex, ey = checks["beam_pts"]
    axes[0].plot(np.array(ex)/d, np.array(ey)/(2*np.pi), "o", color=th["ink"], markersize=5,
                 markerfacecolor=th["surface"], zorder=5, label="integrated")
    axes[0].axvline(0, color=th["ink2"], linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    axes[0].set_xscale("symlog", linthresh=1e-3)
    axes[0].set_xlim(-1.05, 10)
    axes[0].set_ylim(0.8, 5.5)
    axes[0].annotate("separatrix", xy=(0, 5.2), xytext=(8, 0), textcoords="offset points",
                     fontsize=8, color=th["ink2"])
    style(axes[0], th, "$E$ / well depth  (symlog)", "$T\\,\\omega_n / 2\\pi$",
          "Beam, $\\kappa = 1/2$: the log divergence at $E = 0$")
    legend(axes[0], th, loc="upper left")

    # 2. pendulum swing: period ratio against amplitude
    A = np.radians(np.linspace(1, 179.5, 400))
    axes[1].plot(np.degrees(A), pendulum_period(A)/(2*np.pi), color=c0, linewidth=2.2,
                 label="pendulum, exact", zorder=3)
    for kk, c, lab in ((KAPPA_PEND_SLOPE, c1, "piecewise, $\\kappa = 1$"),
                       (KAPPA_PEND_DEPTH, c2, "piecewise, $\\kappa = 0.68$")):
        xx0 = np.pi/(1 + kk)
        T = [period_libration(float(potential(np.pi - a, kk, xx0, True)), kk, xx0)/(2*np.pi)
             for a in A]
        axes[1].plot(np.degrees(A), T, color=c, linewidth=1.8, label=lab, zorder=4)
    axes[1].set_ylim(0.9, 3.0)
    style(axes[1], th, "swing amplitude, degrees", "$T / T_0$",
          "Pendulum swinging: flat until the corner")
    legend(axes[1], th, loc="upper left")

    # 3. pendulum rotating: time per revolution against speed at the bottom
    v = np.linspace(2.02, 8.0, 400)
    axes[2].plot(v, pendulum_rotation_period(v)/(2*np.pi), color=c0, linewidth=2.2,
                 label="pendulum, exact", zorder=3)
    for kk, c, lab in ((KAPPA_PEND_SLOPE, c1, "piecewise, $\\kappa = 1$"),
                       (KAPPA_PEND_DEPTH, c2, "piecewise, $\\kappa = 0.68$")):
        xx0 = np.pi/(1 + kk)
        dd = depth(kk, xx0)
        vv = v[v*v/2 > dd]
        T = [period_rotation(0.5*s*s - dd, kk, xx0)/(2*np.pi) for s in vv]
        axes[2].plot(vv, T, color=c, linewidth=1.8, label=lab, zorder=4)
    style(axes[2], th, "$\\dot{x}$ at the bottom, units of $\\omega_n$", "$T_{rot} / T_0$",
          "Pendulum rotating: time per revolution")
    legend(axes[2], th, loc="upper right")
    fig.tight_layout()
    save(fig, name, "duffing-period")


def fig_basins(th, name, beam, pend):
    """Basins of the damped beam, and winding number of the damped pendulum."""
    (xs_b, vs_b, L_b, man), (xs_p, vs_p, L_p) = beam, pend
    fig, axes = newfig(th, 1, 2, figsize=(11.0, 5.2))
    c0, c1, c2 = th["series"]

    def mix(a, b, f):
        return tuple((1 - f)*np.array(to_rgb(a)) + f*np.array(to_rgb(b)))

    # beam: two categorical fills, the manifold in ink
    cmap = ListedColormap([mix(c1, th["surface"], 0.55), mix(c2, th["surface"], 0.55)])
    axes[0].imshow(L_b, origin="lower", extent=(xs_b[0], xs_b[-1], vs_b[0], vs_b[-1]),
                   cmap=cmap, vmin=-1, vmax=1, aspect="auto", interpolation="nearest", zorder=1)
    for (x, v) in man:
        axes[0].plot(x, v, color=th["ink"], linewidth=1.2, zorder=4)
    k, x0 = KAPPA_DUFF_SLOPE, 1.0/(1 + KAPPA_DUFF_SLOPE)
    xe = wells(k, x0)
    axes[0].plot([0], [0], "o", color=th["ink"], markersize=6, zorder=6)
    axes[0].plot([xe, -xe], [0, 0], "o", color=th["ink"], markersize=6,
                 markerfacecolor=th["surface"], zorder=6)
    axes[0].plot([], [], color=th["ink"], linewidth=1.2, label="saddle's stable manifold")
    axes[0].fill([], [], color=mix(c1, th["surface"], 0.55), label="settles in $-x_e$")
    axes[0].fill([], [], color=mix(c2, th["surface"], 0.55), label="settles in $+x_e$")
    axes[0].set_xlim(xs_b[0], xs_b[-1])
    axes[0].set_ylim(vs_b[0], vs_b[-1])
    style(axes[0], th, "$x$", "$\\dot{x}$",
          f"Beam, $\\zeta = {ZETA_BEAM}$: which well")
    legend(axes[0], th, loc="upper left", bbox_to_anchor=(0.0, -0.14), ncol=1)

    # pendulum: diverging by net turns, neutral grey at zero
    nmax = int(np.max(np.abs(L_p)))
    nmax = max(nmax, 1)
    cols = []
    for n in range(-nmax, nmax + 1):
        if n == 0:
            cols.append(mix(th["axis"], th["surface"], 0.3))
        elif n < 0:
            cols.append(mix(th["div_neg"], th["surface"], 0.5 - 0.25*(min(abs(n), 3) - 1)))
        else:
            cols.append(mix(th["div_pos"], th["surface"], 0.5 - 0.25*(min(n, 3) - 1)))
    axes[1].imshow(L_p, origin="lower", extent=(xs_p[0], xs_p[-1], vs_p[0], vs_p[-1]),
                   cmap=ListedColormap(cols), vmin=-nmax - 0.5, vmax=nmax + 0.5,
                   aspect="auto", interpolation="nearest", zorder=1)
    for shift in (-np.pi, np.pi):
        axes[1].plot([shift], [0], "o", color=th["ink"], markersize=6,
                     markerfacecolor=th["surface"], zorder=6)
    axes[1].plot([0], [0], "o", color=th["ink"], markersize=6, zorder=6)
    for n, lab in ((-3, "three or more backwards"), (-2, "two turns backwards"),
                   (-1, "one turn backwards"), (0, "no full turn"), (1, "one turn forwards"),
                   (2, "two turns forwards"), (3, "three or more forwards")):
        if abs(n) <= nmax:
            axes[1].fill([], [], color=cols[n + nmax], label=lab)
    axes[1].set_xticks([-np.pi, 0, np.pi])
    axes[1].set_xticklabels(["$-\\pi$", "0", "$\\pi$"])
    axes[1].set_xlim(xs_p[0], xs_p[-1])
    axes[1].set_ylim(vs_p[0], vs_p[-1])
    style(axes[1], th, "$x$  (angle from the inverted position)", "$\\dot{x}$",
          f"Pendulum, $\\zeta = {ZETA_PEND}$: full turns before capture")
    legend(axes[1], th, loc="upper left", bbox_to_anchor=(0.0, -0.14), ncol=3)
    fig.tight_layout()
    save(fig, name, "duffing-basins")


def fig_forced(th, name, clouds):
    """Stroboscopic sections of the forced beam and pendulum beside their targets."""
    fig, axes = newfig(th, 2, 2, figsize=(10.0, 8.2))
    c0, c1, _ = th["series"]
    panels = [("duffing", axes[0, 0], c0, f"Duffing, $A = {BEAM_SHOW}$"),
              ("beam-depth", axes[0, 1], c1, f"piecewise beam, $\\kappa = 1/3$, $A = {BEAM_SHOW}$"),
              ("pendulum", axes[1, 0], c0, f"pendulum, $A = {PEND_SHOW}$"),
              ("pend-depth", axes[1, 1], c1, f"piecewise pendulum, $\\kappa = 0.68$, $A = {PEND_SHOW}$")]
    for kind, ax, c, title in panels:
        Y = clouds[kind]
        x = Y[:, 0]
        if kind.startswith("pend"):
            x = reduce(x, np.pi, True)
            ax.set_xlim(-np.pi, np.pi)
            ax.set_xticks([-np.pi, 0, np.pi])
            ax.set_xticklabels(["$-\\pi$", "0", "$\\pi$"])
            xlabel = "$x$  (angle from the inverted position)"
        else:
            xlabel = "$x$"
        ax.plot(x, Y[:, 1], ".", color=c, markersize=1.6, alpha=0.8, zorder=3)
        style(ax, th, xlabel, "$\\dot{x}$", title)
    fig.suptitle("Stroboscopic sections, once per drive period, after the transient",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "duffing-forced")


# ----------------------------------------------------------------- driver
def _table(rows, header, fmt):
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join("---" for _ in header) + " |")
    for r in rows:
        print("| " + " | ".join(f.format(v) for f, v in zip(fmt, r)) + " |")


def beam_points(k=0.5, x0=2.0/3.0):
    """Integrated beam periods at a spread of energies, for the table and the figure."""
    d = depth(k, x0)
    Es = [-0.9*d, -0.5*d, -0.3*d, -0.05*d, -1e-3*d, 1e-3*d, 0.05*d, 0.3*d, 1.5*d, 6.0*d]
    return Es, [period_integrated(E, k, x0, False) for E in Es]


def checks(quick=False):
    """Print every number DUFFING.md quotes, returning what the figures reuse."""
    out = {}
    print("## Equilibria and eigenvalues (wn = 1, zeta = 0.1, kappa = 1/2, x0 = 2/3)")
    for name, xs, res, eig in check_equilibria(0.1, 0.5, 2.0/3.0):
        print(f"  {name:7s} at x = {xs:+.6f}: |f| = {res:.1e}, eigenvalues {np.round(eig, 6)}")
    (lp, lm), (wp, wm) = eigenvalues(0.1, 0.5)
    print(f"  formula: saddle {lp:+.6f}, {lm:+.6f}; well {wp:.6f}, {wm:.6f}")
    print(f"  well depth kappa(1+kappa)x0^2/2 = {depth(0.5, 2/3):.6f}; "
          f"V(xe) = {float(potential(1.0, 0.5, 2/3)):.6f}")

    print("\n## Exact periods against integration")
    k, x0 = 0.5, 2.0/3.0
    d = depth(k, x0)
    rows = []
    pts = out["beam_pts"] = beam_points()
    for E, Tn in zip(*pts):
        Te = period(E, k, x0, False)
        rows.append(("beam", E/d, "well" if E < 0 else "cross", Te, Tn, abs(Te - Tn)))
    k, x0 = 1.0, np.pi/2
    d = depth(k, x0)
    for E in (-0.9*d, -0.3*d, -0.05*d, -1e-3*d, 1e-3*d, 0.3*d, 1.5*d):
        Te, Tn = period(E, k, x0, True), period_integrated(E, k, x0, True)
        rows.append(("pendulum", E/d, "swing" if E < 0 else "rotation", Te, Tn, abs(Te - Tn)))
    _table(rows, ("system", "E / depth", "orbit", "closed form", "integrated", "difference"),
           ("{}", "{:+.3f}", "{}", "{:.9f}", "{:.9f}", "{:.1e}"))

    print("\n## Log divergence at the separatrix: slope of T against -ln|E|, times mu")
    for lab, fn, k, x0, sgn in (("beam, well", period_well, 0.5, 2/3, -1),
                                ("beam, cross", period_cross, 0.5, 2/3, 1),
                                ("pendulum, swing", period_libration, 1.0, np.pi/2, -1),
                                ("pendulum, rotation", period_rotation, 1.0, np.pi/2, 1)):
        mu = np.sqrt(k)
        E1, E2 = sgn*1e-6, sgn*1e-8
        slope = (fn(E2, k, x0) - fn(E1, k, x0))/(np.log(abs(E1)) - np.log(abs(E2)))
        print(f"  {lab:20s} slope * mu = {slope*mu:.5f}")

    print("\n## Pendulum: period ratio against swing amplitude")
    rows = []
    for deg in (10, 30, 60, 90, 120, 150, 170, 179):
        a = np.radians(deg)
        r = [pendulum_period(a)/(2*np.pi)]
        for kk in (KAPPA_PEND_SLOPE, KAPPA_PEND_DEPTH):
            xx0 = np.pi/(1 + kk)
            r.append(period_libration(float(potential(np.pi - a, kk, xx0, True)), kk, xx0)/(2*np.pi))
        rows.append((deg, *r))
    _table(rows, ("amplitude", "pendulum", "kappa = 1", "kappa = 0.68"),
           ("{}", "{:.4f}", "{:.4f}", "{:.4f}"))
    print(f"  corner amplitudes: kappa=1 at {np.degrees(np.pi - np.pi/2):.1f} deg, "
          f"kappa=0.68 at {np.degrees(KAPPA_PEND_DEPTH*np.pi/(1+KAPPA_PEND_DEPTH)):.1f} deg")
    print(f"  well depth: pendulum 2, kappa=1 {depth(1.0, np.pi/2):.4f}, "
          f"kappa=0.68 {depth(KAPPA_PEND_DEPTH, np.pi/(1+KAPPA_PEND_DEPTH)):.4f}")
    print(f"  escape speed from the bottom: pendulum 2, kappa=1 {np.sqrt(2*depth(1.0, np.pi/2)):.4f}, "
          f"kappa=0.68 {np.sqrt(2*depth(KAPPA_PEND_DEPTH, np.pi/(1+KAPPA_PEND_DEPTH))):.4f}")

    print("\n## Pendulum: time per revolution against speed at the bottom")
    rows = []
    for v in (2.05, 2.2, 2.5, 3.0, 4.0, 6.0):
        r = [pendulum_rotation_period(v)/(2*np.pi)]
        for kk in (KAPPA_PEND_SLOPE, KAPPA_PEND_DEPTH):
            xx0 = np.pi/(1 + kk)
            dd = depth(kk, xx0)
            r.append(period_rotation(0.5*v*v - dd, kk, xx0)/(2*np.pi) if 0.5*v*v > dd else np.nan)
        rows.append((v, *r))
    _table(rows, ("bottom speed", "pendulum", "kappa = 1", "kappa = 0.68"),
           ("{:.2f}", "{:.4f}", "{:.4f}", "{:.4f}"))

    print("\n## Duffing (alpha = beta = 1): in-well period ratio against inner amplitude")
    wn = np.sqrt(2.0)
    rows = []
    for a in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99):
        r = [duffing_well_period(a)/(2*np.pi/wn)]
        for kk in (KAPPA_DUFF_SLOPE, KAPPA_DUFF_DEPTH):
            xx0 = 1.0/(1 + kk)
            E = float(potential(1.0 - a, kk, xx0, wn=wn))
            r.append(period_well(E, kk, xx0, wn)/(2*np.pi/wn))
        rows.append((a, *r))
    _table(rows, ("inner amplitude", "Duffing", "kappa = 1/2", "kappa = 1/3"),
           ("{:.2f}", "{:.4f}", "{:.4f}", "{:.4f}"))
    print(f"  well depth: Duffing 0.25, kappa=1/2 {depth(0.5, 2/3, wn):.4f}, "
          f"kappa=1/3 {depth(1/3, 3/4, wn):.4f}")

    print("\n## Action per revolution, closed form against quadrature (pendulum, kappa = 1)")
    for E in (0.01, 0.5, 3.0):
        print(f"  E = {E}: {action(E, 1.0, np.pi/2):.9f}  {action_quad(E, 1.0, np.pi/2):.9f}")

    print(f"\n## Pendulum capture speeds at zeta = {ZETA_PEND}, kappa = 1")
    k, x0 = 1.0, np.pi/2
    vn = capture_speeds(ZETA_PEND, k, x0, n_max=5)
    va = capture_speeds_averaged(ZETA_PEND, k, x0, n_max=5)
    rows = [(n + 1, vn[n], va[n], (vn[n] - vn[n - 1]) if n else np.nan) for n in range(5)]
    _table(rows, ("turns", "integrated v_n", "energy map", "step from previous"),
           ("{}", "{:.4f}", "{:.4f}", "{:.4f}"))
    print(f"  escape speed {np.sqrt(2*depth(k, x0)):.4f}; fast-rotation step 4 zeta wn xe = "
          f"{4*ZETA_PEND*np.pi:.4f}")
    out["capture"] = (vn, va)
    return out


def figures_(quick=False, fresh=False, checks_out=None):
    """Write every figure; the expensive inputs are computed once and cached."""
    if checks_out is None:
        checks_out = {"beam_pts": beam_points()}
    n = 90 if quick else 260
    # beam basins
    k, x0 = KAPPA_DUFF_SLOPE, 1.0/(1 + KAPPA_DUFF_SLOPE)
    xs = np.linspace(-2.2, 2.2, n)
    vs = np.linspace(-1.6, 1.6, int(n*0.72))
    L_b = basins(ZETA_BEAM, k, x0, False, xs, vs,
                 cache=CACHE.format("basin-beam" + ("-quick" if quick else "")))
    man = stable_manifold(ZETA_BEAM, k, x0, T=60.0)
    beam = (xs, vs, L_b, man)
    frac = np.mean(L_b > 0)
    print(f"  beam basins: {frac:.3f} of the window settles in +xe, {1-frac:.3f} in -xe")
    # pendulum winding
    xp = np.linspace(-np.pi, np.pi, n)
    vp = np.linspace(-4.5, 4.5, int(n*0.72))
    L_p = basins(ZETA_PEND, 1.0, np.pi/2, True, xp, vp,
                 cache=CACHE.format("basin-pend" + ("-quick" if quick else "")))
    pend = (xp, vp, L_p)
    print(f"  pendulum window: turns range {L_p.min()} to {L_p.max()}")

    # forced
    amps_b = BEAM_AMPS[::3] if quick else BEAM_AMPS
    amps_p = PEND_AMPS[::3] if quick else PEND_AMPS
    sb = scan(["duffing", "beam-slope", "beam-depth"], amps_b,
              CACHE.format("lyap-beam" + ("-quick" if quick else "")), fresh=fresh)
    sp = scan(["pendulum", "pend-slope", "pend-depth"], amps_p,
              CACHE.format("lyap-pend" + ("-quick" if quick else "")), fresh=fresh)
    print("\n## Largest Lyapunov exponent, forced beam family (Holmes: delta 0.25, Om 1)")
    _table([(a, sb["duffing"][a], sb["beam-slope"][a], sb["beam-depth"][a]) for a in map(float, amps_b)],
           ("A", "Duffing", "kappa = 1/2", "kappa = 1/3"), ("{:.2f}", "{:+.4f}", "{:+.4f}", "{:+.4f}"))
    print("\n## Largest Lyapunov exponent, forced pendulum family (q = 2, Om = 2/3)")
    _table([(a, sp["pendulum"][a], sp["pend-slope"][a], sp["pend-depth"][a]) for a in map(float, amps_p)],
           ("A", "pendulum", "kappa = 1", "kappa = 0.68"), ("{:.2f}", "{:+.4f}", "{:+.4f}", "{:+.4f}"))
    nc = 800 if quick else 4000
    with Pool(4) as p:
        got = p.starmap(strobe, [("duffing", BEAM_SHOW, 300, nc), ("beam-depth", BEAM_SHOW, 300, nc),
                                 ("pendulum", PEND_SHOW, 300, nc), ("pend-depth", PEND_SHOW, 300, nc)])
    clouds = dict(zip(["duffing", "beam-depth", "pendulum", "pend-depth"], got))

    for name, th in THEMES.items():
        fig_force(th, name)
        fig_phase(th, name)
        fig_period(th, name, checks_out)
        fig_basins(th, name, beam, pend)
        fig_forced(th, name, clouds)


if __name__ == "__main__":
    args = sys.argv[1:]
    quick, fresh = "quick" in args, "fresh" in args
    what = [a for a in args if a not in ("quick", "fresh")]
    out = None
    if not what or "checks" in what:
        out = checks(quick)
    if not what or "figures" in what:
        figures_(quick, fresh, out)
