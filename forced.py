"""Sinusoidal forcing of the symmetric deadzone prototype.

Run ``python3 forced.py`` for a self check.

The four unforced prototypes are planar and autonomous, so Poincare-Bendixson
rules out anything more complicated than an equilibrium or a limit cycle.
Adding a drive term supplies the third state (the drive phase) and with it the
possibility of entrainment, quasi-periodicity and chaos::

    xddot + 2 wn [ zm xdot + (zp - zm) dz(xdot) ] + wn^2 x = A cos(Om t)

    dz(v) = v - v0   for v >  v0
            0        for |v| <= v0
            v + v0   for v < -v0

Everything switched is unchanged from ``symmetric.py``; only the right hand
side is new. That matters, because it means the whole unforced analysis
carries over as the ``A = 0`` edge of every result here, and the reference
limit cycle supplies the frequency the drive competes with.

Two natural scales set the axes
-------------------------------

*Frequency.* The drive is compared against the **limit cycle** frequency
``w_lc = 2 pi / T``, not against ``wn``. ``T`` comes from
``symmetric.period_exact`` and depends only on ``wn``, ``zp`` and ``zm``.
For the reference pair ``zp = 0.3``, ``zm = -0.1`` this gives
``T = 6.319387``, ``w_lc = 0.994271``.

*Amplitude.* The only other quantity with dimensions in the unforced problem
is the deadzone half width ``v0``, which fixes the cycle size and nothing
else. The drive is an acceleration, so the dimensionless drive strength is
``a = A / (wn v0)``. Because the unforced system is positively homogeneous
apart from ``v0``, the pair ``(Om/w_lc, a)`` is the complete parameter plane:
scaling ``v0`` and ``A`` together moves the orbit but not its character.
``scaling_check`` verifies this.

What the routines here measure
------------------------------

``strobe`` samples the state once per drive period after a transient. A
periodic response locked to the drive shows up as a finite set of points; the
number of them is the order ``q`` of the lock (the response repeats after
``q`` drive periods). Quasi-periodic motion fills a closed curve, chaos a
fractal cloud.

``rotation_number`` measures how far the orbit winds around the origin per
drive period. On a ``p:q`` lock it sits exactly at ``p/q`` and stays there
over an interval of ``Om`` — the plateau of a devil's staircase, an Arnold
tongue in the ``(Om, A)`` plane.

``lyapunov`` estimates the largest Lyapunov exponent by evolving a nearby
trajectory and renormalising once per drive period. It is negative on a lock,
near zero on a torus, positive on a chaotic attractor.

``classify`` combines the three into one label. The strobe count decides a
lock; failing that the exponent decides between torus and chaos. The
threshold is deliberately loose (see ``LAM_TOL``): near a tongue boundary the
distinction between a long lock and a torus is not numerically decidable, and
the code says ``torus`` where it cannot tell.

A caution on the exponent
-------------------------

The field is continuous but not differentiable across the deadzone edges, so
there is no variational equation to integrate and the exponent has to come
from a finite separation. The value therefore depends mildly on ``d0`` and on
the renormalisation interval. Signs are robust and are what the
classification uses; the magnitudes should be read as indicative.
"""
import numpy as np
from scipy.integrate import solve_ivp

import symmetric
from symmetric import deadzone

WN = 1.0

#: Reference damping pair and deadzone half width, as used throughout the
#: README. ``zm < 0 < zp`` is the existence condition for the unforced cycle.
ZP, ZM, V0 = 0.3, -0.1, 0.25

#: A response counts as locked at order ``q`` when the strobe point repeats
#: after ``q`` drive periods to within this distance.
#:
#: Getting this test right took three goes. Clustering the strobe points at a
#: fixed radius of 1e-4 was wrong twice over. A weakly contracting case
#: settles onto an invariant curve only a few thousandths across, which the
#: fixed radius chopped into around twenty phantom clusters and reported as a
#: lock of order twenty; the give-away was that the count drifted with the
#: length of the discarded transient, 23 points after 200 drive periods and
#: 21 after 2000, which a real lock never does. Making the radius relative to
#: the section's own spread fixed that but left two subtler failures: a
#: *chain of small invariant curves*, one per island, clusters exactly like a
#: lock of the same order, and an orbit still spiralling in splits each of
#: its points into several tight clusters.
#:
#: Testing recurrence directly kills all three. A genuine period ``q``
#: response repeats to integration accuracy, around 1e-9. An island chain
#: sits at 1e-4 to 1e-5 and stays there however long the transient is run. A
#: drifting orbit fails outright. The margin between 1e-9 and 1e-5 is four
#: orders, so the threshold is not delicate.
LOCK_SCATTER = 1e-7

#: Largest lock order the strobe count will report. Beyond this a response is
#: indistinguishable from a torus over any affordable run.
QMAX = 24

#: Exponent magnitudes below this are read as zero. Chosen from the spread of
#: ``lyapunov`` over repeated runs on known locks and known tori.
LAM_TOL = 5e-3


def w_lc(zp=ZP, zm=ZM):
    """Angular frequency of the unforced limit cycle, ``2 pi / T``.

    Delegates the period to ``symmetric.period_exact``, which solves the
    exact half cycle reduction rather than integrating. Independent of
    ``v0``: the amplitude scales with the deadzone width, the period does
    not.
    """
    return 2.0*np.pi/symmetric.period_exact(zp, zm, 1.0)[0]


def field(zp, zm, v0, amp, om):
    """Right hand side of the forced deadzone oscillator.

    Returns ``f(t, y)`` for ``y = [x, xdot]``. Identical to
    ``symmetric.field`` with ``amp cos(om t)`` added to the acceleration,
    so ``amp = 0`` reproduces the unforced prototype exactly.
    """
    d = zp - zm

    def f(t, y):
        return [y[1],
                -WN**2*y[0] - 2.0*WN*(zm*y[1] + d*deadzone(y[1], v0))
                + amp*np.cos(om*t)]
    return f


#: Integrator settings. ``LSODA`` is used rather than the explicit default
#: because the deadzone kinks make an explicit stepper reject steps
#: constantly: on a representative run it agreed with ``RK45`` at
#: ``rtol = 1e-9`` to eight significant figures while taking a fifth of the
#: function evaluations. ``_tolerance_check`` re-confirms that agreement.
METHOD, RTOL, ATOL = "LSODA", 1e-9, 1e-11


def _run(zp, zm, v0, amp, om, t_eval, y0):
    """Integrate the forced field and return the state at ``t_eval``."""
    sol = solve_ivp(field(zp, zm, v0, amp, om), (0.0, t_eval[-1]), y0,
                    t_eval=t_eval, method=METHOD, rtol=RTOL, atol=ATOL)
    return sol.y.T


def strobe(amp, om, zp=ZP, zm=ZM, v0=V0, n_skip=500, n_keep=150, y0=None):
    """Stroboscopic section: the state sampled once per drive period.

    ``n_skip`` drive periods are discarded as transient and the next
    ``n_keep`` returned as an ``(n_keep, 2)`` array of ``(x, xdot)``. The
    default start ``y0`` is on the unforced cycle, which keeps the transient
    short for weak drives.
    """
    td = 2.0*np.pi/om
    if y0 is None:
        y0 = [0.0, 2.0*v0]
    t_eval = td*np.arange(n_skip, n_skip + n_keep + 1)
    return _run(zp, zm, v0, amp, om, np.concatenate(([0.0], t_eval)), y0)[1:]


def lock_order(pts, qmax=QMAX):
    """Order of the lock, or ``None`` if the response is not periodic.

    Returns the smallest ``q`` for which the strobe point repeats after
    ``q`` drive periods to within ``LOCK_SCATTER``, which is the definition
    of a period ``q`` orbit of the stroboscopic map. ``None`` means the
    section is a curve, an island chain or a cloud rather than a finite set
    — or that a genuine lock has not settled yet, which is answered by a
    longer ``n_skip`` rather than a looser threshold.
    """
    for k in range(1, qmax + 1):
        d = np.linalg.norm(pts[k:] - pts[:-k], axis=1)
        if d.size and np.max(d) < LOCK_SCATTER:
            return k
    return None


def rotation_number(amp, om, zp=ZP, zm=ZM, v0=V0, n_skip=150, n_keep=300,
                    y0=None):
    """Orbit windings around the origin per drive period.

    The phase is taken in the ``(x, xdot/wn)`` plane, where the unforced
    linear flow is a clockwise rotation, and unwrapped over ``n_keep`` drive
    periods after a transient. The result is the winding number ``w``: on a
    ``p:q`` lock it equals ``p/q``.

    Only meaningful while the orbit still encircles the origin. Under a
    strong drive it need not, and the number then loses its reading as a
    frequency ratio; ``classify`` reports it but does not rely on it.
    """
    td = 2.0*np.pi/om
    if y0 is None:
        y0 = [0.0, 2.0*v0]
    # sample finely inside each period so the unwrap cannot skip a turn
    per = 8
    t = td*np.arange(0, (n_skip + n_keep)*per + 1)/per
    y = _run(zp, zm, v0, amp, om, t, y0)
    k0 = n_skip*per
    th = np.unwrap(np.arctan2(-y[k0:, 1]/WN, y[k0:, 0]))
    return (th[-1] - th[0])/(2.0*np.pi)/n_keep


def lyapunov(amp, om, zp=ZP, zm=ZM, v0=V0, n_skip=100, n=300, d0=1e-8,
             y0=None):
    """Largest Lyapunov exponent, per unit time, by two-trajectory tracking.

    A twin trajectory is started ``d0`` away and rescaled back to ``d0``
    once per drive period; the exponent is the mean log stretch divided by
    the period. ``n_skip`` periods are run before accumulating so the
    separation aligns with the local expanding direction.

    The field is only piecewise smooth, so this is a finite-difference
    estimate rather than a variational one: read the sign, not the digits.
    """
    td = 2.0*np.pi/om
    if y0 is None:
        y0 = [0.0, 2.0*v0]
    f = field(zp, zm, v0, amp, om)

    def step(state, t0):
        sol = solve_ivp(f, (t0, t0 + td), state, method=METHOD,
                        rtol=RTOL, atol=ATOL)
        return sol.y[:, -1]

    a = np.array(y0, float)
    t0 = 0.0
    for _ in range(n_skip):
        a = step(a, t0)
        t0 += td
    b = a + np.array([d0, 0.0])
    total = 0.0
    for _ in range(n):
        a2, b2 = step(a, t0), step(b, t0)
        t0 += td
        d = np.hypot(*(b2 - a2))
        if d == 0.0:
            a = a2
            b = a2 + np.array([d0, 0.0])
            continue
        total += np.log(d/d0)
        a = a2
        b = a2 + (b2 - a2)*(d0/d)
    return total/(n*td)


def multipliers(amp, om, zp=ZP, zm=ZM, v0=V0, h=1e-6, n_skip=400):
    """Floquet multipliers of an entrained periodic response.

    Settles onto the attractor, takes the strobe point as a fixed point of
    the ``q``-fold stroboscopic map, and differences that map to get its
    ``2 x 2`` Jacobian. The eigenvalues are the multipliers.

    Their product is fixed by the trace formula: the flow's divergence is
    ``-2 zeta wn``, so ``mu1 mu2 = exp(2 Lambda)`` with ``Lambda`` the dwell
    weighted mean of the pole real parts times the period. That is the same
    ``exp(2 Lambda)`` the unforced cycle has, which is why the unforced
    Floquet result survives forcing unchanged.

    What the pair does individually is the part that is not forced by the
    algebra, and it is what decides whether chaos is available: a
    complex conjugate pair contracts every direction by the same factor and
    cannot fold, while a real pair splits the plane into a stretching and a
    contracting direction and can.

    Returns ``(mu, q)``: the multiplier pair as a complex array, and the
    lock order used. ``(None, None)`` if the response is not a lock.
    """
    pts = strobe(amp, om, zp, zm, v0, n_skip=n_skip, n_keep=QMAX + 4)
    q = lock_order(pts)
    if q is None:
        return None, None
    td = 2.0*np.pi/om
    f = field(zp, zm, v0, amp, om)
    t0 = n_skip*td                     # phase the strobe point was taken at
    p = pts[0]

    def qmap(state):
        sol = solve_ivp(f, (t0, t0 + q*td), state, method=METHOD,
                        rtol=RTOL, atol=ATOL)
        return sol.y[:, -1]

    j = np.empty((2, 2))
    for k in range(2):
        e = np.zeros(2)
        e[k] = h
        j[:, k] = (qmap(p + e) - qmap(p - e))/(2.0*h)
    return np.linalg.eigvals(j), q


def dwell_exponent(amp, om, zp=ZP, zm=ZM, v0=V0, n_skip=100, n=300,
                   y0=None):
    """Dwell weighted mean of the pole real parts, measured on the attractor.

    The unforced analysis reduced stability to ``Lambda = s+ t+ + s- t-``,
    the dwell weighted sum of the two half planes' pole real parts, with
    ``s = -zeta wn``. This measures the same quantity on a forced orbit by
    carrying a third state that accumulates the time spent outside the
    deadzone, so no sampling is involved and the fractions are as accurate
    as the integration.

    Returns ``(exponent, f_out)``: the mean of ``s`` per unit time and the
    fraction of the time spent outside the band.
    """
    td = 2.0*np.pi/om
    if y0 is None:
        y0 = [0.0, 2.0*v0]
    d = zp - zm

    def f(t, y):
        return [y[1],
                -WN**2*y[0] - 2.0*WN*(zm*y[1] + d*deadzone(y[1], v0))
                + amp*np.cos(om*t),
                1.0 if abs(y[1]) > v0 else 0.0]

    warm = solve_ivp(f, (0.0, n_skip*td), list(y0) + [0.0], method=METHOD,
                     rtol=RTOL, atol=ATOL)
    y1 = warm.y[:, -1].copy()
    y1[2] = 0.0
    run = solve_ivp(f, (n_skip*td, (n_skip + n)*td), y1, method=METHOD,
                    rtol=RTOL, atol=ATOL)
    f_out = run.y[2, -1]/(n*td)
    return -WN*(zp*f_out + zm*(1.0 - f_out)), f_out


def classify(amp, om, zp=ZP, zm=ZM, v0=V0, quick=False):
    """Label the forced response and return the evidence with it.

    Returns ``(label, q, w, lam)``: the label, the strobe count (``None`` if
    it exceeded ``QMAX``), the rotation number and the Lyapunov exponent.
    Labels are ``"lock p:q"``, ``"torus"`` and ``"chaos"``.

    The strobe count is trusted first because it is a direct observation of
    the section; the exponent is only consulted when the section is not a
    finite set. With ``quick`` the exponent is skipped and anything that is
    not a lock is returned as ``"torus"`` — useful for sweeps where only the
    tongues matter.
    """
    pts = strobe(amp, om, zp, zm, v0)
    q = lock_order(pts)
    w = rotation_number(amp, om, zp, zm, v0)
    if q is not None:
        p = int(round(w*q))
        return ("lock %d:%d" % (p, q), q, w, None)
    if quick:
        return ("torus", None, w, None)
    lam = lyapunov(amp, om, zp, zm, v0)
    return ("chaos" if lam > LAM_TOL else "torus", None, w, lam)


def _stair_point(args):
    """Worker for :func:`staircase` (module level so it can be pickled)."""
    amp, om, zp, zm, v0 = args
    return rotation_number(amp, om, zp, zm, v0)


def _map_point(args):
    """Worker for :func:`regime_map` (module level so it can be pickled)."""
    amp, om, zp, zm, v0 = args
    lab, q, w, lam = classify(amp, om, zp, zm, v0)
    return lab, (q if q is not None else 0), w, (0.0 if lam is None else lam)


def _pool_map(fn, args, workers=None):
    """Map over a process pool, falling back to a serial map."""
    import multiprocessing as mp
    n = workers or mp.cpu_count()
    if n <= 1:
        return [fn(a) for a in args]
    with mp.Pool(n) as pool:
        return pool.map(fn, args, chunksize=1)


def staircase(a, r_lo=0.4, r_hi=2.2, n=181, zp=ZP, zm=ZM, v0=V0,
              workers=None):
    """Rotation number against drive frequency, at one drive strength.

    ``a`` is the dimensionless drive ``A / (wn v0)`` and the sweep runs over
    ``Om / w_lc`` from ``r_lo`` to ``r_hi``. Returns ``(ratios, w)``.

    Read the result as a devil's staircase: flat runs are locks, and the
    width of the ``1:1`` plateau is the ``1:1`` Arnold tongue cut at this
    drive strength. At ``a = 0`` the curve is the smooth unforced ratio
    ``w_lc / Om`` with no plateaus at all, because an autonomous cycle has
    nothing to lock to.
    """
    ratios = np.linspace(r_lo, r_hi, n)
    wl = w_lc(zp, zm)
    args = [(a*WN*v0, r*wl, zp, zm, v0) for r in ratios]
    return ratios, np.array(_pool_map(_stair_point, args, workers))


def regime_map(ratios, amps, zp=ZP, zm=ZM, v0=V0, workers=None):
    """Classify the response over a grid of drive frequency and strength.

    ``ratios`` are ``Om / w_lc`` and ``amps`` are ``A / (wn v0)``. Returns
    ``(labels, q, w, lam)``, each of shape ``(len(amps), len(ratios))``:
    the label strings, the lock orders (0 where the section was not finite),
    the rotation numbers and the Lyapunov exponents (0 where not computed,
    which is exactly where a lock was found).
    """
    args = [(a*WN*v0, r*w_lc(zp, zm), zp, zm, v0)
            for a in amps for r in ratios]
    out = _pool_map(_map_point, args, workers)
    sh = (len(amps), len(ratios))
    lab = np.array([o[0] for o in out], dtype=object).reshape(sh)
    q = np.array([o[1] for o in out]).reshape(sh)
    w = np.array([o[2] for o in out]).reshape(sh)
    lam = np.array([o[3] for o in out]).reshape(sh)
    return lab, q, w, lam


def _lock_point(args):
    """Worker for :func:`tongue_edges`: is this point locked at ``p/q``?"""
    amp, om, zp, zm, v0, target = args
    return abs(rotation_number(amp, om, zp, zm, v0) - target) < W_TOL


#: A rotation number this close to ``p/q`` counts as locked. On a genuine
#: lock the measured number is exact to twelve digits or better, so the
#: threshold is far above the noise and far below any torus.
W_TOL = 1e-6


def tongue_edges(a, p=1, q=1, zp=ZP, zm=ZM, v0=V0, span=1.4, steps=14):
    """Locate the edges of the ``p:q`` Arnold tongue at one drive strength.

    The tongue is the set of drive frequencies for which the rotation number
    is pinned at ``p / q``. It is a single interval containing
    ``Om / w_lc = q / p``, so each edge can be bisected for independently.

    Args:
        a: dimensionless drive ``A / (wn v0)``.
        p, q: the lock to trace. ``1, 1`` is the main tongue.
        zp, zm, v0: the prototype's parameters.
        span: how far either side of ``q/p`` to look for an unlocked point,
            as a multiplier on ``q/p``.
        steps: bisection steps per edge; ``14`` puts the edges inside
            ``5e-5`` of ``Om / w_lc``.

    Returns:
        ``(r_left, r_right)`` in units of ``Om / w_lc``, or ``(nan, nan)``
        if the centre itself is not locked — which is what a tongue that has
        not yet opened at this drive strength looks like.
    """
    amp = a*WN*v0
    wl = w_lc(zp, zm)
    target = p/q
    centre = q/p

    def locked(r):
        return _lock_point((amp, r*wl, zp, zm, v0, target))

    if not locked(centre):
        return float("nan"), float("nan")
    out = []
    for direction in (-1.0, 1.0):
        far = centre*(1.0 + direction*(span - 1.0)) if direction > 0 \
            else centre/span
        if locked(far):
            out.append(far)
            continue
        lo, hi = centre, far
        for _ in range(steps):
            mid = 0.5*(lo + hi)
            if locked(mid):
                lo = mid
            else:
                hi = mid
        out.append(lo)
    return min(out), max(out)


def _edge_point(args):
    """Worker for :func:`tongue_width` (module level so it can be pickled)."""
    return tongue_edges(*args)


def tongue_width(amps, p=1, q=1, zp=ZP, zm=ZM, v0=V0, workers=None):
    """Trace one tongue's edges over a list of drive strengths.

    Returns ``(left, right)`` arrays the same length as ``amps``. Each edge
    search is a serial bisection, so the parallelism is across ``amps``.
    Rows where the tongue has not opened yet come back as ``nan``.
    """
    out = _pool_map(_edge_point,
                    [(a, p, q, zp, zm, v0) for a in amps], workers)
    return np.array([o[0] for o in out]), np.array([o[1] for o in out])


def scaling_check(amp=0.3, om=None, factor=3.0):
    """Confirm that ``(Om/w_lc, A/(wn v0))`` is the whole parameter plane.

    Scaling ``v0`` and ``A`` by the same factor should scale the orbit by
    that factor and leave the rotation number untouched. Returns
    ``(w_ref, w_scaled, ratio_of_orbit_extents)``.
    """
    if om is None:
        om = w_lc()
    w1 = rotation_number(amp, om, v0=V0)
    w2 = rotation_number(amp*factor, om, v0=V0*factor)
    p1 = strobe(amp, om, v0=V0)
    p2 = strobe(amp*factor, om, v0=V0*factor)
    return w1, w2, np.max(np.abs(p2))/np.max(np.abs(p1))


def derivative_check(amp=0.5, om=1.0, zp=ZP, zm=ZM, v0=V0, n=400):
    """Check that the derivative correspondence survives forcing.

    Differentiating the velocity switched equation turns it into the
    displacement switched one with the drive amplitude scaled by ``Om``:
    if ``x`` solves the deadzone model under ``A cos(Om t)`` then ``X = xdot``
    solves the displacement switched model under ``-A Om sin(Om t)``. So the
    peak velocity of the first equals the peak displacement of the second.
    Returns ``(max|xdot|, max|X|)``.
    """
    d = zp - zm
    td = 2.0*np.pi/om
    t = np.linspace(0.0, n*td, n*400)

    def g(tt, y):                      # displacement switched, derivative form
        z = zp if abs(y[0]) > v0 else zm
        return [y[1],
                -WN**2*y[0] - 2.0*WN*z*y[1] - amp*om*np.sin(om*tt)]

    ya = _run(zp, zm, v0, amp, om, t, [0.0, 2.0*v0])
    sol = solve_ivp(g, (0.0, t[-1]), [2.0*v0, ya[0, 1]], t_eval=t,
                    method=METHOD, rtol=RTOL, atol=ATOL)
    half = len(t)//2
    return (np.max(np.abs(ya[half:, 1])), np.max(np.abs(sol.y[0, half:])))


if __name__ == "__main__":
    wl = w_lc()
    print("unforced reference: zp=%.2f zm=%.2f v0=%.2f" % (ZP, ZM, V0))
    print("  T = %.6f   w_lc = %.6f" % (2.0*np.pi/wl, wl))

    print("\nweak drive should entrain 1:1")
    for r in (0.98, 1.00, 1.02):
        lab, q, w, lam = classify(0.05, r*wl, quick=True)
        print("  Om/w_lc = %.2f  ->  %-10s w = %.6f" % (r, lab, w))

    print("\nscaling: (Om/w_lc, A/(wn v0)) is the whole plane")
    w1, w2, ratio = scaling_check()
    print("  w = %.6f and %.6f, orbit extent ratio %.6f (expect 3)"
          % (w1, w2, ratio))

    print("\nderivative correspondence under forcing")
    a, b = derivative_check()
    print("  max|xdot| = %.6f   max|X| = %.6f   diff = %.2e"
          % (a, b, abs(a - b)))
