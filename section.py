"""Model agnostic stroboscopic analysis of a periodically driven oscillator.

Run ``python3 section.py`` for a self check, which is a cross validation
against ``forced.py``.

``forced.py`` built this machinery for one specific system, the deadzone
prototype. Comparing that prototype against a *smooth* oscillator needs the
same measurements applied to a different field, and re-implementing them
would invite the two copies to drift apart — which matters more than usual
here, because getting the lock test right took four attempts and every wrong
version produced confident false results.

So the engine lives here and takes the flow as an argument. Everything is
expressed in terms of a callable ``flow(t, y)`` for a planar state and a
drive period ``td``; nothing below knows what system it is looking at.

What is measured
----------------

``strobe`` samples the state once per drive period. A response locked to the
drive shows a finite set of points, quasi-periodic motion a closed curve,
chaos a fractal cloud.

``lock_order`` returns the smallest ``q`` for which the strobe point repeats
after ``q`` drive periods, which is the definition of a period ``q`` orbit
of the stroboscopic map.

``rotation_number`` counts orbit windings around a centre per drive period.
On a ``p:q`` lock it sits at ``p/q``.

``lyapunov`` estimates the largest Lyapunov exponent from a twin trajectory
renormalised once per drive period.

``multipliers`` differences the ``q``-fold stroboscopic map at its fixed
point, giving the Floquet multipliers of a locked response.

Why every threshold is relative
-------------------------------

Both the orbit's size and its settling time vary by orders of magnitude
across a parameter plane, so a fixed tolerance is a claim about a region
rather than a method. Two failures in ``forced.py`` came from exactly this:
an absolute recurrence threshold reported four cells as chaotic, with
Lyapunov exponents up to ``+0.30``, purely because the resonant orbit at the
centre of a tongue is ten times larger than off resonance and carries ten
times the integration error — all four are periodic; and a fixed transient
length reported an entire Arnold tongue as absent, when every cell in it
locks once the transient is given time to die.

The thresholds below are the ones that survived. They are stated relative to
the orbit, and the settling time is computed from the contraction rather
than fixed.
"""
import numpy as np
from scipy.integrate import solve_ivp

#: Integrator. ``LSODA`` handles both the deadzone prototype's kinks and a
#: smooth field; on a representative run it agreed with ``RK45`` at
#: ``rtol = 1e-9`` to eight significant figures for a fifth of the function
#: evaluations.
METHOD, RTOL, ATOL = "LSODA", 1e-9, 1e-11

#: A response is locked at order ``q`` when the strobe point repeats after
#: ``q`` drive periods to within this fraction of the orbit's own size. A
#: genuine lock repeats to about 1e-8 relative, an island chain holds at
#: 1e-4 however long the transient runs, so there are two orders of margin
#: either side.
LOCK_SCATTER_REL = 1e-6
LOCK_SCATTER_FLOOR = 1e-12

#: Twin trajectory separation, as a fraction of the orbit's size. It must sit
#: well above the integrator's relative error and well below the attractor's
#: curvature.
D0_REL = 1e-6

#: Largest lock order reported. Beyond this a response is indistinguishable
#: from a torus over any affordable run.
QMAX = 24

#: Exponent magnitudes below this read as zero.
#:
#: Set from the measured spread, not chosen. Re-running the exponent on the
#: cells a first pass called chaotic, at five times the run length and a
#: hundred times the separation, moved marginal values by up to 0.008 and
#: flipped four of seven verdicts: +0.0071 became -0.0003, +0.0079 became
#: +0.0004, +0.0123 became -0.0001, +0.0041 became -0.0010. The three that
#: held -- +0.0938, +0.1005, +0.0362 -- barely moved. So the noise floor on
#: these estimates is about 0.008 and a threshold of 5e-3 sat inside it.
#: At 2e-2 the confirmed and rejected cells separate with a factor of two
#: below and a factor of two above, and re-thresholding the stored maps at
#: 2e-2 reproduces the convergence test's verdicts exactly.
#:
#: A cell near the threshold still deserves :func:`confirm_chaos` rather
#: than trust.
LAM_TOL = 2e-2

#: A rotation number this close to ``p/q`` counts as locked.
W_TOL = 1e-6


def run(flow, t_eval, y0):
    """Integrate ``flow`` and return the state at each time in ``t_eval``."""
    sol = solve_ivp(flow, (t_eval[0], t_eval[-1]), y0, t_eval=t_eval,
                    method=METHOD, rtol=RTOL, atol=ATOL)
    return sol.y.T


def strobe(flow, td, y0, n_skip, n_keep=150):
    """State sampled once per drive period, after discarding a transient.

    Args:
        flow: ``f(t, y)`` for the planar state.
        td: drive period.
        y0: initial state.
        n_skip: drive periods to discard.
        n_keep: samples to return.

    Returns:
        ``(n_keep, 2)`` array of states.
    """
    t = td*np.arange(n_skip, n_skip + n_keep + 1)
    return run(flow, np.concatenate(([0.0], t)), y0)[1:]


def lock_order(pts, qmax=QMAX):
    """Smallest ``q`` with the strobe point repeating after ``q`` periods.

    ``None`` if the section is a curve, an island chain or a cloud rather
    than a finite set — or if a genuine lock has not settled yet, which is
    answered by a longer transient rather than a looser threshold.
    """
    tol = max(LOCK_SCATTER_FLOOR,
              LOCK_SCATTER_REL*float(np.max(np.abs(pts))))
    for k in range(1, qmax + 1):
        d = np.linalg.norm(pts[k:] - pts[:-k], axis=1)
        if d.size and np.max(d) < tol:
            return k
    return None


def lock_margin(pts, qmax=QMAX):
    """Smallest recurrence residual over ``q``, relative to the orbit's size.

    ``lock_order`` answers yes or no against ``LOCK_SCATTER_REL``; this
    reports how close the call was, which matters because the threshold is
    not always above the noise it is competing with.

    Integrating a strongly damped piecewise system accumulates error faster
    than the threshold allows for. On a 17 level staircase driven at the
    chaotic point, a genuine **period 4 orbit** gives a residual of
    1.5e-6 against a threshold of 1e-6, at every integrator tolerance from
    1e-9 down to 1e-12 — the error floor, not the dynamics. The orbit was
    therefore reported as unlocked and then, its exponent being negative, as
    a torus. The exact map in ``maps.py`` resolves the same orbit at 2.4e-9,
    three orders inside the threshold.

    So a "torus" verdict from this module means *no lock was detected*, not
    *no lock exists*. Where the margin is within an order or so of the
    threshold, believe the exact map instead.

    Returns ``(q_best, residual)``.
    """
    scale = max(LOCK_SCATTER_FLOOR, float(np.max(np.abs(pts))))
    best_q, best = None, np.inf
    for k in range(1, qmax + 1):
        d = np.linalg.norm(pts[k:] - pts[:-k], axis=1)
        if d.size:
            r = float(np.max(d))/scale
            if r < best:
                best_q, best = k, r
    return best_q, best


def rotation_number(flow, td, y0, n_skip, n_keep=300, centre=(0.0, 0.0),
                    per=8):
    """Windings of the orbit about ``centre`` per drive period.

    The phase is taken in the state plane and unwrapped over ``n_keep`` drive
    periods, sampled ``per`` times within each so the unwrap cannot skip a
    turn. On a ``p:q`` lock the result is ``p/q``.

    Only meaningful while the orbit encircles ``centre``; under a strong
    drive it need not, and the number then loses its reading as a frequency
    ratio.
    """
    t = td*np.arange(0, (n_skip + n_keep)*per + 1)/per
    y = run(flow, t, y0)
    k0 = n_skip*per
    th = np.unwrap(np.arctan2(-(y[k0:, 1] - centre[1]),
                              y[k0:, 0] - centre[0]))
    return (th[-1] - th[0])/(2.0*np.pi)/n_keep


def lyapunov(flow, td, y0, n_skip, n=300, d0=None):
    """Largest Lyapunov exponent per unit time, by twin trajectory tracking.

    A twin is started ``d0`` away and rescaled once per drive period; the
    exponent is the mean log stretch over the period. ``d0`` defaults to
    ``D0_REL`` times the size the orbit actually reaches, measured after the
    transient rather than assumed.

    For a field that is only piecewise smooth this is a finite difference
    estimate rather than a variational one: read the sign, not the digits.
    """
    def step(state, t0):
        sol = solve_ivp(flow, (t0, t0 + td), state, method=METHOD,
                        rtol=RTOL, atol=ATOL)
        return sol.y[:, -1]

    a = np.array(y0, float)
    t0 = 0.0
    scale = 0.0
    for _ in range(n_skip):
        a = step(a, t0)
        t0 += td
        scale = max(scale, float(np.max(np.abs(a))))
    if d0 is None:
        d0 = max(LOCK_SCATTER_FLOOR, D0_REL*scale)
    b = a + np.array([d0, 0.0])
    total = 0.0
    for _ in range(n):
        a2, b2 = step(a, t0), step(b, t0)
        t0 += td
        d = np.hypot(*(b2 - a2))
        if d == 0.0:
            a, b = a2, a2 + np.array([d0, 0.0])
            continue
        total += np.log(d/d0)
        a = a2
        b = a2 + (b2 - a2)*(d0/d)
    return total/(n*td)


def multipliers(flow, td, y0, n_skip, h_rel=1e-6):
    """Floquet multipliers of a locked response.

    Settles onto the attractor, takes the strobe point as a fixed point of
    the ``q``-fold stroboscopic map, and differences that map for its
    Jacobian. A complex conjugate pair contracts every direction alike; a
    real pair splits the plane into a stretching and a contracting
    direction, which is what a fold needs.

    Returns ``(mu, q)``, or ``(None, None)`` if the response is not a lock.
    """
    pts = strobe(flow, td, y0, n_skip, n_keep=QMAX + 4)
    q = lock_order(pts)
    if q is None:
        return None, None
    p = pts[0]
    t0 = n_skip*td
    h = max(1e-12, h_rel*float(np.max(np.abs(pts))))

    def qmap(state):
        sol = solve_ivp(flow, (t0, t0 + q*td), state, method=METHOD,
                        rtol=RTOL, atol=ATOL)
        return sol.y[:, -1]

    j = np.empty((2, 2))
    for k in range(2):
        e = np.zeros(2)
        e[k] = h
        j[:, k] = (qmap(p + e) - qmap(p - e))/(2.0*h)
    return np.linalg.eigvals(j), q


def confirm_chaos(flow, td, y0, n_skip, n_short=300, n_long=1500):
    """Re-test a positive exponent for convergence before believing it.

    A single Lyapunov estimate near the threshold is not evidence: the
    estimate is a finite difference on a finite run, and both the run length
    and the separation ``d0`` bias it. This recomputes the exponent at five
    times the run length and again at a hundred times the separation, and
    calls the cell chaotic only if every estimate clears ``LAM_TOL``.

    On the cells this was built for, the genuinely chaotic ones moved by
    less than 0.008 across all three estimates while four marginal ones
    changed sign.

    Returns ``(is_chaos, (lam_short, lam_long, lam_wide))``.
    """
    scale = float(np.max(np.abs(strobe(flow, td, y0, n_skip, n_keep=200))))
    ls = lyapunov(flow, td, y0, n_skip//2, n=n_short)
    ll = lyapunov(flow, td, y0, n_skip, n=n_long)
    lw = lyapunov(flow, td, y0, n_skip, n=n_long,
                  d0=max(LOCK_SCATTER_FLOOR, 1e-4*scale))
    return (min(ls, ll, lw) > LAM_TOL), (ls, ll, lw)


def settle_periods(mu_cycle, r=1.0, tol=1e-9, floor=500, cap=8000):
    """Drive periods to discard, from the unforced cycle's contraction.

    A transient decays by ``mu_cycle`` per cycle, so reaching ``tol`` takes
    ``log(tol)/log(mu_cycle)`` cycles, and a cycle is ``r`` drive periods.
    This varies by more than an order of magnitude across the cases of
    interest, which is why it cannot be a constant.
    """
    if not (0.0 < mu_cycle < 1.0):
        return floor
    return int(min(cap, max(floor,
                            np.ceil(np.log(tol)/np.log(mu_cycle)*r))))


def classify(flow, td, y0, n_skip, centre=(0.0, 0.0), quick=False):
    """Label the response and return the evidence with it.

    Returns ``(label, q, w, lam)`` with labels ``"lock p:q"``, ``"torus"``
    and ``"chaos"``. The strobe count is trusted first because it is a direct
    observation of the section; the exponent is consulted only when the
    section is not a finite set.
    """
    pts = strobe(flow, td, y0, n_skip)
    q = lock_order(pts)
    w = rotation_number(flow, td, y0, max(150, n_skip//3), centre=centre)
    if q is not None:
        return ("lock %d:%d" % (int(round(w*q)), q), q, w, None)
    # A near miss is not a torus. The threshold competes with the
    # integrator's own error, so a residual within an order of it means the
    # test could not decide; say so rather than defaulting to "torus".
    qm, margin = lock_margin(pts)
    if margin < 30.0*LOCK_SCATTER_REL:
        return ("undecided(q~%d)" % qm, None, w, None)
    if quick:
        return ("torus", None, w, None)
    lam = lyapunov(flow, td, y0, max(100, n_skip//5))
    return ("chaos" if lam > LAM_TOL else "torus", None, w, lam)


if __name__ == "__main__":
    # Cross validation: drive this engine with the deadzone field and check
    # it reproduces forced.py's own answers on the cases that module's lock
    # test was debugged against. Two independent code paths agreeing is
    # better evidence than one path re-run.
    import forced

    cases = ((0.3, -0.1, 0.45, 1.00, 1), (0.3, -0.1, 0.60, 1.35, 4),
             (0.3, -0.1, 1.00, 2.00, 2), (0.3, -0.1, 0.30, 0.60, None),
             (0.01, -0.003, 0.60, 0.50, None),
             (0.01, -0.003, 0.10, 1.40, None),
             (0.02, -0.006, 0.10, 1.80, None),
             (0.01, -0.003, 1.00, 1.10, 1))
    print("cross validation against forced.py")
    print("%7s %9s %6s %6s %8s %8s %8s"
          % ("zp", "zm", "a", "r", "forced", "section", "expect"))
    ok = True
    for zp, zm, a, r in ((c[0], c[1], c[2], c[3]) for c in cases):
        pass
    for zp, zm, a, r, want in cases:
        wl = forced.w_lc(zp, zm)
        om = r*wl
        qf = forced.lock_order(forced.strobe(a*forced.V0, om, zp, zm,
                                             forced.V0))
        ns = settle_periods(forced.contraction(zp, zm), r)
        qs = lock_order(strobe(forced.field(zp, zm, forced.V0, a*forced.V0,
                                            om),
                               2.0*np.pi/om, [0.0, 2.0*forced.V0], ns))
        good = (qf == qs == want)
        ok &= good
        print("%7.3f %9.4f %6.2f %6.2f %8s %8s %8s  %s"
              % (zp, zm, a, r, qf, qs, want, "ok" if good else "MISMATCH"))
    print("\nboth paths agree with the expected classification:", ok)
