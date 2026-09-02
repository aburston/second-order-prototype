"""Van der Pol under the same forcing analysis, as a control.

Run ``python3 vanderpol.py`` for a self check.

The prototypes in this repository put the nonlinearity in a *switch*: the
damping ratio takes one of two values either side of a boundary. Van der Pol
puts it in a *polynomial*::

    xddot - mu (1 - x^2) xdot + wn^2 x = A cos(Om t)

Both are second order oscillators with nonlinear damping and a limit cycle,
and forcing either one is the same experiment. Running the identical
measurements on both is the point of this module: it says which conclusions
were about nonlinear damping in general and which were about the switch.

Three structural differences, which is what the comparison turns on
---------------------------------------------------------------------

*The damping is unbounded rather than saturating.* Outside the deadzone the
prototype's damping ratio is exactly ``zp``, however large the orbit; Van der
Pol's is ``-mu(1 - x^2)``, which grows without limit. So the prototype is
asymptotically linear at large amplitude and Van der Pol is never linear
anywhere.

*There is no free amplitude scale.* The prototype's cycle amplitude is
proportional to the deadzone half width ``v0``, so ``v0`` can be scaled out
and the forced problem has exactly two parameters, ``Om/w_lc`` and
``A/(wn v0)``. Van der Pol's cycle sits at ``x ~ 2`` whatever ``mu`` is,
fixed by the polynomial; nothing can be scaled out, and ``mu`` stays a third
parameter. ``scale_check`` demonstrates the failure directly.

*The switch is on displacement, not velocity.* ``(1 - x^2)`` depends on
``x``, so the closer prototype is the displacement switched member of the
four, not the deadzone. The deadzone is still the one that carries the
forcing results, and the two are related by differentiation.

What is measured, and by what
-----------------------------

Everything comes from ``section.py``, the same engine the deadzone results
came from, driven with this field instead. Nothing in the classification is
re-implemented here, which is deliberate: that lock test took four attempts
to get right and a second copy would drift from the first.
"""
import numpy as np
from scipy.integrate import solve_ivp

import section

WN = 1.0

#: Relaxation parameter values spanning the range of behaviour: nearly
#: harmonic, moderately nonlinear, and a relaxation oscillator.
MU_NEAR_HARMONIC, MU_MODERATE, MU_RELAXATION = 0.1, 1.0, 5.0


def field(mu, amp=0.0, om=1.0):
    """Right hand side of the forced Van der Pol oscillator.

    Returns ``f(t, y)`` for ``y = [x, xdot]``. With ``amp = 0`` this is the
    autonomous oscillator, whose limit cycle the forcing competes with.
    """
    def f(t, y):
        return [y[1],
                mu*(1.0 - y[0]**2)*y[1] - WN**2*y[0] + amp*np.cos(om*t)]
    return f


def cycle(mu, t_settle=200.0, t_scan=200.0):
    """Period and a point on the unforced limit cycle.

    Settles onto the cycle, then times successive maxima of ``x`` — the
    crossings of ``xdot = 0`` with ``x`` positive — which is a clean section
    because the cycle encircles the origin once per period.

    Returns:
        ``(T, y_on_cycle)``.
    """
    f = field(mu)
    warm = solve_ivp(f, (0.0, t_settle), [2.0, 0.0], method=section.METHOD,
                     rtol=section.RTOL, atol=section.ATOL)
    y0 = warm.y[:, -1]

    def top(t, y):
        return y[1]
    top.direction = -1.0                       # xdot falling through zero

    sol = solve_ivp(f, (0.0, t_scan), y0, method=section.METHOD,
                    rtol=section.RTOL, atol=section.ATOL, events=top,
                    dense_output=True)
    ts = [t for t in sol.t_events[0] if t > 1e-6]
    if len(ts) < 3:
        return float("nan"), y0
    # average several periods rather than trust one crossing
    T = float((ts[-1] - ts[1])/(len(ts) - 2))
    return T, np.array(sol.sol(ts[1]))


def w_lc(mu):
    """Angular frequency of the unforced limit cycle."""
    T, _ = cycle(mu)
    return 2.0*np.pi/T


def amplitude(mu):
    """Peak displacement of the unforced limit cycle."""
    T, y = cycle(mu)
    if not np.isfinite(T):
        return float("nan")
    sol = solve_ivp(field(mu), (0.0, T), y, method=section.METHOD,
                    rtol=section.RTOL, atol=section.ATOL,
                    t_eval=np.linspace(0.0, T, 4000))
    return float(np.max(np.abs(sol.y[0])))


def contraction(mu):
    """Non-trivial Floquet multiplier of the unforced limit cycle.

    The autonomous period map has one multiplier at exactly 1, the neutral
    direction along the cycle; the other is how fast a transient collapses
    onto it, and is the quantity the deadzone analysis called
    ``exp(2 Lambda)``. Measured here by differencing the period map, since
    there is no closed form to appeal to.

    Returns the multiplier furthest from 1, or ``nan`` when the value is not
    resolvable — see :func:`contraction_resolved`, which this calls. A
    relaxation oscillator contracts so hard that its transverse multiplier
    falls below anything a finite difference in double precision can see,
    and returning a number there would be reporting noise.
    """
    val, _ = contraction_resolved(mu)
    return val


def contraction_resolved(mu, steps=(1e-7, 1e-6, 1e-5, 1e-4, 1e-3),
                         spread=0.25):
    """Measure the cycle's multiplier and say whether it is resolved at all.

    A single finite difference cannot tell a genuinely tiny multiplier from
    its own rounding error, and for a relaxation oscillator the difference
    matters: at ``mu = 5`` sweeping the step over five orders gives -0.746,
    -6.99e-4, +6.12e-4, +2.38e-4, -4.96e-5, -2.32e-3 and -0.249 — the sign
    flips and the magnitude moves four orders, so no digit of it is real.
    At ``mu = 0.1`` the same sweep holds 0.53307 across six orders, and at
    ``mu = 1`` it holds 8.6e-4 to about a tenth. So the test is agreement
    across steps, not the value from any one of them.

    Returns ``(value, resolved)``: the median estimate, and whether the
    estimates agree to within ``spread`` in relative terms and share a sign.
    An unresolved case gets ``(nan, False)`` — its true multiplier is
    smaller than double precision can express through this map, which is
    itself the finding.
    """
    T, y = cycle(mu)
    if not np.isfinite(T):
        return float("nan"), False
    f = field(mu)

    def pmap(state):
        sol = solve_ivp(f, (0.0, T), state, method=section.METHOD,
                        rtol=section.RTOL, atol=section.ATOL)
        return sol.y[:, -1]

    vals = []
    for h in steps:
        j = np.empty((2, 2))
        for k in range(2):
            e = np.zeros(2)
            e[k] = h
            j[:, k] = (pmap(y + e) - pmap(y - e))/(2.0*h)
        ev = np.linalg.eigvals(j)
        vals.append(float(np.real(ev[np.argmax(np.abs(ev - 1.0))])))
    v = np.array(vals)
    med = float(np.median(v))
    ok = bool(np.all(np.sign(v) == np.sign(med)) and med != 0.0
              and np.max(np.abs(v - med))/abs(med) < spread)
    return (med if ok else float("nan")), ok


def _settle(mu, r, floor=200):
    """Drive periods to discard, from this oscillator's own contraction."""
    return section.settle_periods(abs(contraction(mu)), r, floor=floor)


def classify(mu, a, r, n_skip=None, quick=False):
    """Label the forced response at drive strength ``a``, ratio ``r``.

    ``a`` is the drive amplitude ``A`` directly, because Van der Pol has no
    free scale to divide it by; ``r`` is ``Om / w_lc``.
    """
    wl = w_lc(mu)
    om = r*wl
    if n_skip is None:
        n_skip = _settle(mu, r)
    _, y = cycle(mu)
    return section.classify(field(mu, a, om), 2.0*np.pi/om, list(y),
                            n_skip, quick=quick)


def _cls_point(args):
    """Worker for :func:`regime_map` (module level so it can be pickled)."""
    mu, a, r = args
    lab, q, w, lam = classify(mu, a, r)
    return lab, (q if q is not None else 0), w, (0.0 if lam is None else lam)


def regime_map(mu, ratios, amps, workers=None):
    """Classify the response over a grid of drive frequency and amplitude.

    Returns ``(labels, q, w, lam)``, each of shape ``(len(amps),
    len(ratios))``.
    """
    import multiprocessing as mp
    args = [(mu, a, r) for a in amps for r in ratios]
    n = workers or mp.cpu_count()
    with mp.Pool(n) as pool:
        out = pool.map(_cls_point, args, chunksize=1)
    sh = (len(amps), len(ratios))
    return (np.array([o[0] for o in out], dtype=object).reshape(sh),
            np.array([o[1] for o in out]).reshape(sh),
            np.array([o[2] for o in out]).reshape(sh),
            np.array([o[3] for o in out]).reshape(sh))


def scale_check(mu=1.0, factor=3.0):
    """Show that Van der Pol has no free amplitude scale to divide out.

    The deadzone prototype's cycle amplitude is proportional to ``v0``, so
    scaling ``v0`` and ``A`` together rescales the orbit and leaves the
    dynamics alone — which is why its forced problem has two parameters.
    Van der Pol has no such parameter: the cycle sits where the polynomial
    puts it. Scaling the drive by ``factor`` therefore does *not* scale the
    orbit by ``factor``.

    Returns ``(ratio_of_orbit_extents, factor)``; the prototype returns the
    factor exactly, Van der Pol does not.
    """
    wl = w_lc(mu)
    _, y = cycle(mu)
    out = []
    for a in (0.2, 0.2*factor):
        pts = section.strobe(field(mu, a, wl), 2.0*np.pi/wl, list(y),
                             _settle(mu, 1.0), n_keep=200)
        out.append(float(np.max(np.abs(pts))))
    return out[1]/out[0], factor


if __name__ == "__main__":
    print("unforced Van der Pol")
    print("%8s %10s %10s %14s" % ("mu", "T", "amplitude", "contraction"))
    for mu in (0.1, 0.5, 1.0, 2.0, 5.0):
        T, _ = cycle(mu)
        print("%8.2f %10.5f %10.5f %14.3e"
              % (mu, T, amplitude(mu), contraction(mu)))

    print("\namplitude stays near 2 whatever mu: that is the fixed scale")
    print("no free scale to divide the drive by:")
    got, want = scale_check()
    print("  orbit extent ratio %.4f for a drive scaled by %.1f"
          " (the prototype gives %.1f exactly)" % (got, want, want))
