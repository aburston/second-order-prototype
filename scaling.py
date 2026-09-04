"""Moving the three level prototype to another frequency range.

Run ``python3 scaling.py`` for the self check, which is the evidence behind
the scaling section of ``THREELEVEL.md``; ``python3 scaling.py table``
prints the worked examples that section quotes.

The rule
--------

Every number in ``THREELEVEL.md`` is for the reference model, ``wn = 1``
and edges of order 2::

    x'' + 2 zeta(x) x' + x = A cos(Om t),   zeta switched at |x| = a, b

To place the same behaviour at natural frequency ``wn`` and amplitude
scale ``lam``, keep the three ratios, multiply the edges by ``lam``, the
drive frequency by ``wn`` and the drive acceleration by ``lam wn^2``::

    zeta_k -> zeta_k          a, b -> lam a, lam b
    Om     -> wn Om           A    -> lam wn^2 A

Then ``y(t) = lam x(wn t)`` solves the scaled equation exactly, so the
scaled model is the reference model with its clock run ``wn`` times faster
and its ruler ``lam`` times coarser: every period divides by ``wn``, every
displacement multiplies by ``lam``, every velocity by ``lam wn``, and
everything dimensionless -- lock orders, plateau edges in ``Om/w_lc``,
map multipliers, chaotic verdicts, ``mu`` -- is unchanged. A Lyapunov
exponent, a rate, multiplies by ``wn``.

The proof is one substitution. With ``y = lam x(wn t)``, ``y' = lam wn x'``
and ``y'' = lam wn^2 x''``, so::

    y'' + 2 zeta(y/lam) wn y' + wn^2 y
      = lam wn^2 [x'' + 2 zeta(x) x' + x]
      = lam wn^2 A cos(Om wn t)

and ``zeta(y/lam)`` is the reference staircase because the edges were
scaled with ``y``. The one non-obvious factor is ``wn^2`` on the drive: an
acceleration amplitude must scale with the spring's ``wn^2 x``, not with
``wn``, which is why the dimensionless drive strength is ``A/(wn^2 b)``.

What is checked
---------------

1. Trajectories: the scaled model integrated in its own units, divided by
   ``lam`` and read at ``t/wn``, is the reference trajectory, unforced and
   under a drive, at two scales three orders of magnitude apart.
2. The free cycle: the scaled model's amplitude is ``lam`` times the exact
   reference radius and its period ``1/wn`` times the exact period.
3. The classifier: driven at scaled frequencies, the scaled model gets the
   same verdict -- lock 3, chaos, lock 4, chaos -- as the reference, and its
   Lyapunov exponent divided by ``wn`` is the reference exponent. The
   check also shows where the classifier's own floors are not scale free.
4. The drive strength: ``A/(wn^2 b)`` is invariant, ``A/(wn b)`` is not.
5. A Van der Pol oscillator in physical units, ``x'' - eps (1 - x^2/X^2) x'
   + wn^2 x = 0``, is the dimensionless one at ``mu = eps/wn`` with
   ``lam = X``: integrated at 50 Hz it has Van der Pol's period over
   ``wn`` and amplitude times ``X``.
"""
import numpy as np
from scipy.integrate import solve_ivp

import section
import staircase

#: The reference model: the three level fit to Van der Pol at mu = 5.
REF_LEVELS, REF_EDGES = staircase.THREE_FITTED

#: The two physical scales the check runs at: a 50 Hz oscillator with a
#: millimetre per model unit, and a kilohertz one with a micrometre.
SCALES = ((2.0*np.pi*50.0, 1.0e-3), (2.0*np.pi*1000.0, 1.0e-6))

#: Drive used for the classifier check, the one the model was fitted at.
DRIVE_AMP = staircase.CMP_AMP


# --------------------------------------------------------------- the rule
def scale_model(levels, edges, wn, lam):
    """Levels and edges of the reference model at ``wn`` and ``lam``."""
    return tuple(levels), tuple(lam*e for e in edges)


def scale_drive(amp, om, wn, lam):
    """Drive acceleration and frequency: ``A -> lam wn^2 A``, ``Om -> wn Om``."""
    return lam*wn**2*amp, wn*om


def scale_state(y, wn, lam):
    """A phase plane point: ``x -> lam x``, ``xdot -> lam wn xdot``."""
    return lam*y[0], lam*wn*y[1]


def field(levels, edges, wn=1.0, amp=0.0, om=1.0):
    """The staircase oscillator at any natural frequency.

    `staircase.field` is the ``wn = 1`` case; this is the same right hand
    side with ``wn`` written in, so the scaled model can be integrated in
    its own units rather than by trusting the rule it is meant to test.
    """
    lv, ed = np.asarray(levels, float), np.asarray(edges, float)

    def f(t, y):
        z = lv[int(np.searchsorted(ed, abs(y[0]), "right"))]
        return [y[1], -wn**2*y[0] - 2.0*z*wn*y[1] + amp*np.cos(om*t)]
    return f


def drive_strength(amp, wn, b):
    """The dimensionless drive strength, ``A/(wn^2 b)``."""
    return amp/(wn**2*b)


# ------------------------------------------------------------- the checks
def _integrate(f, y0, t_eval, scale):
    """Integrate with tolerances relative to the model's own scale."""
    atol = np.array([1e-13*scale[0], 1e-13*scale[1]])
    sol = solve_ivp(f, (t_eval[0], t_eval[-1]), list(y0), method="LSODA",
                    rtol=1e-10, atol=atol, t_eval=t_eval)
    return sol.y


def check_trajectory(wn, lam, amp=0.0, om=1.0, y0=(2.5, 0.0), t_end=60.0,
                     n=3000):
    """Max deviation of the rescaled trajectory from the reference, relative
    to the reference amplitude, in ``x`` and in ``xdot``."""
    t = np.linspace(0.0, t_end, n)
    ref = _integrate(field(REF_LEVELS, REF_EDGES, 1.0, amp, om), y0, t,
                     (1.0, 1.0))
    lv, ed = scale_model(REF_LEVELS, REF_EDGES, wn, lam)
    amp_s, om_s = scale_drive(amp, om, wn, lam)
    sc = _integrate(field(lv, ed, wn, amp_s, om_s), scale_state(y0, wn, lam),
                    t/wn, (lam, lam*wn))
    back = np.vstack([sc[0]/lam, sc[1]/(lam*wn)])
    return (np.max(np.abs(back[0] - ref[0]))/np.max(np.abs(ref[0])),
            np.max(np.abs(back[1] - ref[1]))/np.max(np.abs(ref[1])))


def free_cycle_scaled(wn, lam, n_settle=60, n_time=40):
    """Amplitude and period of the scaled model's free cycle by integration,
    timing successive maxima of ``x`` after settling."""
    lv, ed = scale_model(REF_LEVELS, REF_EDGES, wn, lam)
    f = field(lv, ed, wn)
    T_guess = 2.0*np.pi/wn
    warm = solve_ivp(f, (0.0, 2.0*n_settle*T_guess), [2.5*lam, 0.0],
                     method="LSODA", rtol=1e-10,
                     atol=[1e-13*lam, 1e-13*lam*wn])

    def top(t, y):
        return y[1]
    top.direction = -1.0
    sol = solve_ivp(f, (0.0, 2.0*n_time*T_guess), warm.y[:, -1],
                    method="LSODA", rtol=1e-10,
                    atol=[1e-13*lam, 1e-13*lam*wn], events=top,
                    dense_output=True)
    te = sol.t_events[0]
    te = te[te > 1e-8/wn]
    T = float((te[-1] - te[1])/(len(te) - 2))
    return float(sol.sol(te[1])[0]), T


def _classify(args):
    """Worker: classify one driven case, reference or scaled."""
    wn, lam, om = args
    lv, ed = scale_model(REF_LEVELS, REF_EDGES, wn, lam)
    amp_s, om_s = scale_drive(DRIVE_AMP, om, wn, lam)
    flow = field(lv, ed, wn, amp_s, om_s)
    td = 2.0*np.pi/om_s
    y0 = list(scale_state((2.0, 0.0), wn, lam))
    q = section.lock_order(section.strobe(flow, td, y0, staircase.CMP_NSKIP))
    if q:
        return "lock%d" % q, 0.0
    lam_e = section.lyapunov(flow, td, y0, staircase.CMP_NSKIP//2, n=400)
    return ("chaos" if lam_e/wn > section.LAM_TOL else "torus"), lam_e


#: Drive frequencies for the classifier check, in reference units: inside
#: the 3:1 plateau, in the first chaotic band, on the 4:1 lock, in the
#: second chaotic band, on the 5:1 plateau.
CLASSIFY_OMS = (2.40, 2.435, 2.445, 2.470, 2.52)


def check_classifier(wn, lam, oms=CLASSIFY_OMS, workers=None):
    """Verdicts and exponents for the reference and the scaled model."""
    from multiprocessing import Pool
    jobs = [(1.0, 1.0, om) for om in oms] + [(wn, lam, om) for om in oms]
    with Pool(workers) as pool:
        out = pool.map(_classify, jobs)
    n = len(oms)
    return list(zip(oms, out[:n], out[n:]))


def vdp_physical_cycle(eps, wn, X, n_settle=60, n_time=40):
    """Free amplitude and period of ``x'' - eps (1 - x^2/X^2) x' + wn^2 x = 0``
    by integration, timing successive maxima after settling."""
    def f(t, y):
        return [y[1], eps*(1.0 - (y[0]/X)**2)*y[1] - wn**2*y[0]]
    T_guess = 2.0*np.pi/wn
    warm = solve_ivp(f, (0.0, 2.0*n_settle*T_guess), [2.5*X, 0.0],
                     method="LSODA", rtol=1e-10, atol=[1e-13*X, 1e-13*X*wn])

    def top(t, y):
        return y[1]
    top.direction = -1.0
    sol = solve_ivp(f, (0.0, 2.0*n_time*T_guess), warm.y[:, -1],
                    method="LSODA", rtol=1e-10, atol=[1e-13*X, 1e-13*X*wn],
                    events=top, dense_output=True)
    te = sol.t_events[0]
    te = te[te > 1e-8/wn]
    T = float((te[-1] - te[1])/(len(te) - 2))
    return float(sol.sol(te[1])[0]), T


# ---------------------------------------------------------- worked table
def worked_examples(scales=SCALES):
    """The reference model beside its scaled copies, in physical units."""
    R, T = staircase.free_cycle(REF_LEVELS, REF_EDGES)
    rows = []
    for wn, lam in ((1.0, 1.0),) + tuple(scales):
        lv, ed = scale_model(REF_LEVELS, REF_EDGES, wn, lam)
        amp_s, om_s = scale_drive(DRIVE_AMP, 2.47, wn, lam)
        rows.append(dict(wn=wn, lam=lam, edges=ed, amp=amp_s, om=om_s,
                         R=lam*R, T=T/wn, w_lc=2.0*np.pi*wn/T,
                         strength=drive_strength(amp_s, wn, ed[1]),
                         c_over_m=tuple(2.0*z*wn for z in lv)))
    return rows


def print_table(rows):
    print("%-14s %-9s %-9s %-11s %-11s %-11s %-11s %-9s %-11s"
          % ("wn [rad/s]", "lam", "b", "A", "Om [rad/s]", "free R", "free T [s]",
             "A/(wn^2 b)", "2 zeta_2 wn"))
    for r in rows:
        print("%-14.6g %-9.3g %-9.4g %-11.5g %-11.5g %-11.4g %-11.4g %-9.4f %-11.5g"
              % (r["wn"], r["lam"], r["edges"][1], r["amp"], r["om"], r["R"],
                 r["T"], r["strength"], r["c_over_m"][2]))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "table":
        print_table(worked_examples())
        sys.exit(0)

    R, T = staircase.free_cycle(REF_LEVELS, REF_EDGES)
    print("reference: mu = 5 fit, levels %s edges %s"
          % (tuple(round(z, 3) for z in REF_LEVELS),
             tuple(round(e, 3) for e in REF_EDGES)))
    print("  exact free cycle: R = %.6f  T = %.6f\n" % (R, T))

    print("1. trajectories: |scaled/lam - reference| / max|reference|")
    for wn, lam in SCALES:
        for amp, om, what in ((0.0, 1.0, "unforced"),
                              (DRIVE_AMP, 2.40, "driven, lock 3"),
                              (DRIVE_AMP, 2.47, "driven, chaotic band, 60 units")):
            dx, dv = check_trajectory(wn, lam, amp, om)
            print("  wn = %8.2f  lam = %6.0e  %-34s x %.1e  xdot %.1e"
                  % (wn, lam, what, dx, dv))
    print()

    print("2. free cycle by integration of the scaled model")
    for wn, lam in SCALES:
        Rs, Ts = free_cycle_scaled(wn, lam)
        print("  wn = %8.2f  lam = %6.0e  R/lam = %.6f (%.1e)  T*wn = %.6f (%.1e)"
              % (wn, lam, Rs/lam, abs(Rs/lam - R)/R, Ts*wn, abs(Ts*wn - T)/T))
    print()

    print("3. drive strength: A = %g at wn = 1 with b = %.4f" % (DRIVE_AMP, REF_EDGES[1]))
    for wn, lam in SCALES:
        lv, ed = scale_model(REF_LEVELS, REF_EDGES, wn, lam)
        amp_s, _ = scale_drive(DRIVE_AMP, 2.47, wn, lam)
        print("  wn = %8.2f  lam = %6.0e  A/(wn^2 b) = %.6f   A/(wn b) = %.6g"
              % (wn, lam, drive_strength(amp_s, wn, ed[1]), amp_s/(wn*ed[1])))
    print("  reference                  A/(wn^2 b) = %.6f   A/(wn b) = %.6g\n"
          % (drive_strength(DRIVE_AMP, 1.0, REF_EDGES[1]),
             DRIVE_AMP/REF_EDGES[1]))

    wn, lam = SCALES[0]
    print("4. classifier at A = %g, wn = %.2f, lam = %.0e  (verdict, exponent)"
          % (DRIVE_AMP, wn, lam))
    print("   %-7s %-22s %-22s %s" % ("Om", "reference", "scaled", "scaled/wn"))
    for om, (lr, er), (ls, es) in check_classifier(wn, lam):
        print("   %-7.3f %-7s %-14.5f %-7s %-14.5f %.5f"
              % (om, lr, er, ls, es, es/wn))

    wn, X = SCALES[0]
    mu = staircase.CMP_MU
    Rp, Tp = vdp_physical_cycle(mu*wn, wn, X)
    print("5. physical Van der Pol, eps = mu wn = %.1f, wn = %.2f, X = %.0e"
          % (mu*wn, wn, X))
    print("  R/X = %.4f  T*wn = %.4f   Van der Pol at mu = %g: R = %.4f  T = %.4f\n"
          % (Rp/X, Tp*wn, mu, staircase.VDP_R, staircase.VDP_T))
