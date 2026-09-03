"""Two thresholds instead of one: a staircase of damping levels.

Run ``python3 staircase.py`` for a self check.

Every prototype so far switches the damping ratio **once**, between two
values, across a single boundary or a symmetric pair of them. This module
switches it more than once. Keeping the symmetric displacement case — the
one `displacement.py` calls a piecewise constant Van der Pol — and adding a
second threshold gives five zones in the phase plane::

    |x| > b        outer      zeta_2
    a < |x| < b    middle     zeta_1
    |x| < a        core       zeta_0
    a < |x| < b    middle     zeta_1     (mirror)
    |x| > b        outer      zeta_2     (mirror)

so three damping levels and two thresholds ``0 < a < b``. The code takes an
arbitrary number of levels, because the two interesting questions both need
that generality: whether a second threshold buys new behaviour, and how
closely a staircase of ``n`` levels can imitate a genuinely smooth
nonlinearity.

    xddot + 2 zeta(x) wn xdot + wn^2 x = A cos(Om t)

    zeta(x) = levels[k],  k = how many thresholds |x| exceeds

What the second threshold buys
------------------------------

**A second limit cycle, and with it bistability.** The single threshold
model has an effective damping that runs monotonically from ``zeta_0`` at
small amplitude to ``zeta_1`` at large amplitude, so it can cross zero only
once and has at most one limit cycle. Three levels let the effective damping
turn twice. With ``zeta_0 > 0``, ``zeta_1 < 0``, ``zeta_2 > 0`` the origin
attracts, a band of negative damping sits above it, and heavy damping caps
it: an unstable cycle nested inside a stable one, with the origin and the
outer cycle both attracting and the inner cycle the boundary between their
basins.

That is hard excitation — a system that sits quietly until something knocks
it past a threshold and then runs away to a large oscillation, and stays
there. It cannot be represented by any of the four earlier prototypes, and
it is a common failure mode in the field.

The effective damping, which is where the analysis lives
-------------------------------------------------------

For a near circular orbit of radius ``R`` the damping force does work at a
rate proportional to ``zeta(x) xdot^2``, so the cycle averaged damping
weights each zone by its share of ``xdot^2``. Writing ``x = R cos(theta)``,
the zone ``|x| > e`` occupies the angles where ``|cos theta| > e/R``, and
its share of the weight is

    w(e/R) = (2 phi - sin 2 phi) / pi,    phi = arccos(e/R)

which is the same ``2 phi - sin 2 phi`` the single threshold models already
use. Stacking the levels,

    <zeta>(R) = zeta_0 + sum_k (zeta_k - zeta_{k-1}) w(e_k / R)

and a limit cycle sits wherever this crosses zero. For two levels it
reduces to ``2 phi - sin 2 phi = pi rho`` exactly as before, so the earlier
amplitude equation is the ``n = 2`` case of this one. ``check_reduces``
verifies that against `displacement.py`.

This is an averaging result, exact only as the damping goes to zero. The
exact answer comes from ``half_map``, which composes the arcs analytically,
and the two are compared in the self check rather than assumed to agree.
"""
import numpy as np
from scipy.optimize import brentq

from frequency import kernels, TMAX_STEPS

WN = 1.0
NGRID = 24000

#: The worked three level example: quiet at the origin, self exciting in a
#: band, heavily damped beyond it. Bistable.
BISTABLE_LEVELS, BISTABLE_EDGES = (0.15, -0.25, 0.40), (0.6, 1.6)


def zeta_at(x, levels, edges):
    """Damping ratio at displacement ``x``."""
    return levels[int(np.searchsorted(np.asarray(edges), abs(x), "right"))]


def field(levels, edges, amp=0.0, om=1.0):
    """Right hand side of the staircase oscillator, forced or not."""
    lv, ed = np.asarray(levels, float), np.asarray(edges, float)

    def f(t, y):
        z = lv[int(np.searchsorted(ed, abs(y[0]), "right"))]
        return [y[1], -WN**2*y[0] - 2.0*z*WN*y[1] + amp*np.cos(om*t)]
    return f


def _arc_to_x(zeta, xi, vi, xtarget):
    """Advance within one zone until the displacement reaches ``xtarget``.

    Every zone is an oscillation about the origin — the damping term
    vanishes at ``xdot = 0`` whatever ``x`` is — so no centre offset is
    needed and the kernels apply directly.

    Returns ``(t, v)`` on arrival, or ``(nan, nan)`` if the target is never
    reached, which is how the caller learns the orbit turned around first.
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
            t = brentq(xt, g[k - 1], g[k], xtol=1e-15, rtol=8.9e-16)
            c, s = kernels(zeta, t)
            return t, vi*c - (WN**2*xi + zeta*WN*vi)*s
    return np.nan, np.nan


def _arc_to_turn(zeta, xi, vi):
    """Advance within one zone until the velocity returns to zero.

    Returns ``(t, x)`` at the turning point. Used for the last arc of a half
    cycle, and for any arc where the orbit turns before reaching the next
    threshold.
    """
    def vt(t):
        c, s = kernels(zeta, t)
        return vi*c - (WN**2*xi + zeta*WN*vi)*s

    for tmax in TMAX_STEPS:
        g = np.linspace(1e-9, tmax, NGRID)
        val = vt(g)
        ok = np.isfinite(val)
        idx = np.nonzero(ok & (np.sign(val) == -np.sign(val[0])))[0]
        if len(idx):
            k = idx[0]
            t = brentq(vt, g[k - 1], g[k], xtol=1e-15, rtol=8.9e-16)
            c, s = kernels(zeta, t)
            return t, xi*c + (vi + zeta*WN*xi)*s
    return np.nan, np.nan


def half_map(amp0, levels, edges):
    """Exact half cycle map: peak ``amp0`` to the next peak, and its duration.

    Starts at ``(amp0, 0)`` and runs to the next zero of the velocity, which
    by oddness is at ``-amp1``. Displacement decreases monotonically over
    the whole half cycle — within a zone the flow is a linear oscillation
    about the origin, whose only velocity zeros are its own extremes — so
    the thresholds are met in a fixed order and each arc is solved exactly.

    Returns ``(amp1, t_half)``, or ``(nan, nan)`` if an arc fails.
    """
    walls = sorted([-e for e in edges] + list(edges), reverse=True)
    x, v, t_tot = float(amp0), 0.0, 0.0
    for _ in range(2*len(edges) + 2):
        eps = 1e-12*max(1.0, abs(x))
        z = zeta_at(x - eps, levels, edges)
        nxt = next((w for w in walls if w < x - eps), None)
        # Both events have to be timed and the earlier one taken. Solving
        # only for the wall is wrong: the arc formula is the linear solution
        # for this zone, which goes on oscillating for ever, so its first
        # arrival at a distant wall can be several oscillations after the
        # physical trajectory turned round and left the zone. That produced
        # half cycles of 9.06 against a half period of pi, an amplitude
        # jumping from 0.59 to 2.33 across a step of 0.01 in the initial
        # peak, and a Floquet multiplier of 1e12.
        t_turn, xturn = _arc_to_turn(z, x, v)
        t_wall, vwall = (_arc_to_x(z, x, v, nxt) if nxt is not None
                         else (np.nan, np.nan))
        if np.isfinite(t_wall) and (not np.isfinite(t_turn)
                                    or t_wall < t_turn):
            x, v, t_tot = nxt, vwall, t_tot + t_wall
            continue
        if not np.isfinite(t_turn):
            return np.nan, np.nan
        return -xturn, t_tot + t_turn
    return np.nan, np.nan


def w_share(c):
    """Share of the cycle's ``xdot^2`` weight lying beyond ``|x| = c R``."""
    if c >= 1.0:
        return 0.0
    if c <= 0.0:
        return 1.0
    phi = np.arccos(c)
    return (2.0*phi - np.sin(2.0*phi))/np.pi


def mean_damping(r, levels, edges):
    """Cycle averaged damping ratio of a near circular orbit of radius ``r``.

    Zero where a limit cycle sits, positive where the amplitude decays and
    negative where it grows, so its sign is the local stability of that
    amplitude and its zero crossings are the cycles.
    """
    z = levels[0]
    for k, e in enumerate(edges):
        z += (levels[k + 1] - levels[k])*w_share(e/r)
    return z


def cycles_predicted(levels, edges, rmax=40.0, n=4000):
    """Amplitudes where the averaged damping changes sign.

    Returns a list of ``(radius, stable)``. Stability is the sign of the
    slope: the averaged damping rising through zero means larger orbits are
    damped and smaller ones are driven, which is a stable cycle.
    """
    grid = np.linspace(min(edges)*0.2, rmax, n)
    val = np.array([mean_damping(r, levels, edges) for r in grid])
    out = []
    for k in range(len(grid) - 1):
        if np.isfinite(val[k]) and np.isfinite(val[k + 1]) \
                and val[k]*val[k + 1] < 0:
            r = brentq(mean_damping, grid[k], grid[k + 1],
                       args=(levels, edges), xtol=1e-13)
            out.append((r, val[k + 1] > val[k]))
    return out


def cycles_exact(levels, edges, rmax=40.0, n=600):
    """Amplitudes fixed by the exact half cycle map, with their multipliers.

    A symmetric cycle is a fixed point of the half cycle map, because the
    field is odd. Returns a list of ``(radius, multiplier, stable)``, the
    multiplier being that of the **full** cycle, ``(dH/dA)^2``.

    The multiplier comes back as ``nan`` where it is not resolvable. That is
    not a failure but a measurement: fitted to Van der Pol at ``mu = 5`` the
    outer levels reach ``zeta = 20``, and one pass through them destroys all
    memory of the starting amplitude — ``H(r+h)`` and ``H(r-h)`` come back
    bit identical even for ``h`` three per cent of ``r``. The true
    multiplier is then smaller than double precision can express through
    this map, and any number printed for it would be rounding noise. Van
    der Pol itself does the same thing at that ``mu``; see
    ``vanderpol.contraction_resolved``.
    """
    grid = np.linspace(min(edges)*0.2, rmax, n)
    res = np.array([half_map(a, levels, edges)[0] - a for a in grid])
    out = []
    for k in range(len(grid) - 1):
        if np.isfinite(res[k]) and np.isfinite(res[k + 1]) \
                and res[k]*res[k + 1] < 0:
            r = brentq(lambda a: half_map(a, levels, edges)[0] - a,
                       grid[k], grid[k + 1], xtol=1e-13)
            ds = []
            for hrel in (1e-5, 1e-4, 1e-3):
                h = hrel*r
                ds.append((half_map(r + h, levels, edges)[0]
                           - half_map(r - h, levels, edges)[0])/(2.0*h))
            d = float(np.median(ds))
            resolved = d != 0.0 and all(np.sign(x) == np.sign(d) for x in ds)
            out.append((r, d*d if resolved else float("nan"),
                        abs(d) < 1.0))
    return out


def period(r, levels, edges):
    """Exact period of the cycle at radius ``r`` — twice the half cycle."""
    return 2.0*half_map(r, levels, edges)[1]


# ----------------------------------------------------- imitating Van der Pol
def vdp_zeta(x, mu):
    """Van der Pol's damping ratio, ``-mu (1 - x^2) / 2``.

    Van der Pol is written ``xddot - mu(1 - x^2) xdot + x = 0``; matching it
    to ``xddot + 2 zeta wn xdot + wn^2 x = 0`` with ``wn = 1`` gives this.
    """
    return -0.5*mu*(1.0 - x**2)


def vdp_staircase(mu, n, xmax=3.0):
    """Staircase of ``n`` levels approximating Van der Pol out to ``xmax``.

    Thresholds are placed evenly over ``(0, xmax)`` and each level is the
    mean of Van der Pol's damping ratio across its zone, which is the choice
    that makes the staircase's averaged damping match at large radius rather
    than at a point.

    The outermost zone is unbounded and its level is therefore constant: a
    staircase saturates however many levels it has, and Van der Pol does
    not. Raising ``n`` and ``xmax`` pushes the saturation outwards without
    ever removing it, and that is the structural difference the comparison
    is about.

    Returns ``(levels, edges)``.
    """
    edges = tuple(np.linspace(0.0, xmax, n + 1)[1:-1])
    bounds = (0.0,) + edges + (xmax,)
    levels = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        xs = np.linspace(lo, hi, 400)
        levels.append(float(np.mean(vdp_zeta(xs, mu))))
    return tuple(levels), edges


# --------------------------------------------- driving it beside Van der Pol
#: Drive and fit used for the side by side. ``A`` and ``Om`` near 2.466 are
#: the classic chaotic forced Van der Pol case, and ``XMAX`` comfortably
#: contains the driven orbit, which never leaves ``|x| ~ 2.15`` — so the
#: staircase's outer plateau is never visited and the only thing the level
#: count changes is how finely the nonlinearity is resolved where the orbit
#: actually goes.
CMP_MU, CMP_XMAX, CMP_AMP = 5.0, 3.0, 5.0
CMP_NSKIP = 400

_SCAN_CACHE = {}


def forced_label(flow, om, n_skip=CMP_NSKIP, y0=(2.0, 0.0)):
    """Classify a driven response as a lock, a torus or chaos.

    Uses `section.py`, the same engine the earlier forcing work used, so the
    staircase and Van der Pol are measured identically and the only
    difference between them is the field.

    Returns ``(label, lam)``; ``lam`` is zero on a lock, where it is not
    computed.
    """
    import section
    td = 2.0*np.pi/om
    q = section.lock_order(section.strobe(flow, td, list(y0), n_skip))
    if q:
        return "lock%d" % q, 0.0
    lam = section.lyapunov(flow, td, list(y0), n_skip//2, n=400)
    return ("chaos" if lam > section.LAM_TOL else "torus"), lam


def _scan_one(args):
    """Worker for :func:`window_scan` (module level so it can be pickled)."""
    tag, om, mu, xmax, amp = args
    if tag == "vdp":
        import vanderpol
        flow = vanderpol.field(mu, amp, om)
    else:
        lv, ed = vdp_staircase(mu, int(tag), xmax)
        flow = field(lv, ed, amp, om)
    return tag, om, forced_label(flow, om)


def window_scan(oms, level_counts, mu=CMP_MU, xmax=CMP_XMAX, amp=CMP_AMP,
                workers=None):
    """Classify staircase and Van der Pol responses across drive frequency.

    A coarse sweep here is worse than none. The chaotic bands sit in the
    narrow transitions between one lock and the next, and a grid stepping
    0.05 in ``Om`` jumped straight over Van der Pol's, reporting zero
    chaotic points for a system confirmed chaotic at 2.466 by the same code
    minutes earlier. Sample finely enough to resolve the transition, not the
    range.

    Returns ``{tag: [(om, label, lam), ...]}`` with ``"vdp"`` among the tags.
    Results are cached per argument set so both themes of a figure share one
    scan.
    """
    import multiprocessing as mp
    key = (tuple(oms), tuple(level_counts), mu, xmax, amp)
    if key in _SCAN_CACHE:
        return _SCAN_CACHE[key]
    tags = [str(n) for n in level_counts] + ["vdp"]
    args = [(t, om, mu, xmax, amp) for t in tags for om in oms]
    with mp.Pool(workers or mp.cpu_count()) as pool:
        out = pool.map(_scan_one, args, chunksize=1)
    res = {t: [] for t in tags}
    for tag, om, (lab, lam) in out:
        res[tag].append((om, lab, lam))
    for t in tags:
        res[t].sort()
    _SCAN_CACHE[key] = res
    return res


def window_agreement(scan):
    """Share of Van der Pol's chaotic frequencies each staircase reproduces.

    Returns ``{tag: (n_chaotic, n_shared, jaccard)}``. The Jaccard index —
    shared divided by combined — is the honest measure here because a
    staircase can be chaotic where Van der Pol is not, and that counts
    against agreement just as a miss does.
    """
    ref = {om for om, lab, _ in scan["vdp"] if lab == "chaos"}
    out = {}
    for tag, rows in scan.items():
        s = {om for om, lab, _ in rows if lab == "chaos"}
        union = len(s | ref)
        out[tag] = (len(s), len(s & ref), (len(s & ref)/union) if union else
                    float("nan"))
    return out


#: The comparison window of ``VANDERPOL.md``, and a wide one for the level floor.
NARROW_OMS = tuple(np.round(np.linspace(2.40, 2.56, 33), 4))
WIDE_OMS = tuple(np.round(np.arange(1.80, 3.2001, 0.005), 4))
FLOOR_LEVELS = (2, 3, 5)


def confirm_scan(scan, mu=CMP_MU, xmax=CMP_XMAX, amp=CMP_AMP):
    """Re-test every chaotic verdict in a scan with `section.confirm_chaos`.

    Returns ``{tag: [(om, confirmed, (lam_short, lam_long, lam_wide))]}``
    over the chaotic cells only.
    """
    import section
    out = {}
    for tag, rows in scan.items():
        out[tag] = []
        for om, lab, lam in rows:
            if lab != "chaos":
                continue
            if tag == "vdp":
                import vanderpol
                flow = vanderpol.field(mu, amp, om)
            else:
                lv, ed = vdp_staircase(mu, int(tag), xmax)
                flow = field(lv, ed, amp, om)
            ok, lams = section.confirm_chaos(flow, 2.0*np.pi/om, [2.0, 0.0],
                                             CMP_NSKIP)
            out[tag].append((om, ok, lams))
    return out


def runs(rows):
    """Compress ``[(om, label, lam), ...]`` into ``[(label, om_lo, om_hi, n)]``."""
    out = []
    for om, lab, _ in rows:
        if out and out[-1][0] == lab:
            out[-1] = (lab, out[-1][1], om, out[-1][3] + 1)
        else:
            out.append((lab, om, om, 1))
    return out


def level_floor(level_counts=FLOOR_LEVELS, workers=None):
    """Sweep the coarsest staircases for chaos, to find the level floor.

    The level count comparison established chaos at 9 levels and above at the
    classic point, and found three chaotic frequencies at 5 levels in the
    narrow window, but 2 and 3 levels were only ever tested at the single
    frequency 2.466. This sweeps them: first over that comparison's window, then
    over a wide one, ``WIDE_OMS``, in case their period-adding transitions
    sit elsewhere — the coarser the staircase the further its free cycle is
    from Van der Pol's, so its locks need not be where Van der Pol's are.
    Every chaotic verdict is then re-tested with `section.confirm_chaos`.

    Prints the tables ``VANDERPOL.md`` quotes and returns
    ``(narrow, wide, confirmed_narrow, confirmed_wide)``.
    """
    import time
    t0 = time.time()
    narrow = window_scan(NARROW_OMS, level_counts, workers=workers)
    print("narrow window %.3f to %.3f, %d points, %.0fs"
          % (NARROW_OMS[0], NARROW_OMS[-1], len(NARROW_OMS), time.time() - t0),
          flush=True)
    _print_scan(narrow)
    t0 = time.time()
    wide = window_scan(WIDE_OMS, level_counts, workers=workers)
    print("\nwide window %.3f to %.3f, %d points, %.0fs"
          % (WIDE_OMS[0], WIDE_OMS[-1], len(WIDE_OMS), time.time() - t0),
          flush=True)
    _print_scan(wide)
    cn, cw = confirm_scan(narrow), confirm_scan(wide)
    print("\nconfirmation of every chaotic verdict, wide window")
    for tag, rows in cw.items():
        for om, ok, (ls, ll, lw) in rows:
            print("  %5s  Om = %.3f  %s  lam short %+.4f long %+.4f wide %+.4f"
                  % (tag, om, "confirmed" if ok else "REJECTED", ls, ll, lw),
                  flush=True)
    return narrow, wide, cn, cw


def floor_crosschecks(mu=CMP_MU, xmax=CMP_XMAX, amp=CMP_AMP):
    """Three independent checks on the coarse staircases' chaotic cells.

    The exact-Jacobian exponent from `maps.py`, which has no noise floor;
    the same cells from several initial conditions, which is what tells a
    coexisting lock from a fragile verdict; and the forcing chapter's deadzone
    prototype at the corresponding drive. The two-level staircase is the
    displacement-switched model, and the deadzone model is its
    derivative: a drive ``A cos(Om t)`` on the staircase is ``(A/Om) sin``
    on the deadzone, with the same damping ratios and ``v0 = x0``.
    """
    import maps
    import section
    import forced
    print("exact-Jacobian Lyapunov exponent, per unit time")
    for n, om in ((2, 1.890), (2, 3.000), (2, 2.500), (3, 2.585), (3, 2.300)):
        lv, ed = vdp_staircase(mu, n, xmax)
        lam = maps.forced_lyapunov(maps.staircase_model(lv, ed), [2.0, 0.0],
                                   amp, om, n_skip=400, n=1500)
        print("  %d levels  Om = %.3f  %+.4f" % (n, om, lam), flush=True)

    print("two levels from several initial conditions")
    lv, ed = vdp_staircase(mu, 2, xmax)
    for om in (1.890, 3.000):
        flow = field(lv, ed, amp, om)
        for y0 in ([2.0, 0.0], [0.5, 0.0], [-1.0, 3.0], [1.7, -2.0]):
            lab, q, w, lam = section.classify(flow, 2.0*np.pi/om, y0, CMP_NSKIP)
            print("  Om = %.3f  y0 = %s  %s  w = %.4f" % (om, y0, lab, w),
                  flush=True)

    zm, zp = lv
    v0 = float(ed[0])
    print("the deadzone prototype, zp = %.3f zm = %.3f v0 = %.2f, drive A/Om"
          % (zp, zm, v0))
    for om in (1.890, 3.000, 2.500):
        flow = forced.field(zp, zm, v0, amp/om, om)
        print("  Om = %.3f  r = %.2f  a = %.2f"
              % (om, om/forced.w_lc(zp, zm), amp/om/v0))
        for y0 in ([2.0, 0.0], [0.0, 2.0], [0.5, 0.5], [-10.0, 0.0]):
            lab, q, w, lam = section.classify(flow, 2.0*np.pi/om, y0, CMP_NSKIP)
            print("    y0 = %s  %s  w = %.4f  lam = %s"
                  % (y0, lab, w, "-" if lam is None else "%+.4f" % lam),
                  flush=True)


# ------------------------------------------ normalising the coarse models
#: Van der Pol's free cycle at ``CMP_MU``: the amplitude and period every
#: normalised model below is made to match.
VDP_R, VDP_T = 2.0215, 11.612
NORM_OMS = tuple(np.round(np.arange(2.30, 2.7001, 0.005), 4))


def free_cycle(levels, edges):
    """Radius on the section and period of the outermost stable cycle."""
    ex = cycles_exact(levels, edges)
    r = ex[-1][0]
    return r, period(r, levels, edges)


def free_cycle_num(levels, edges, n_settle=60, rtol=1e-9):
    """Radius and period of the free cycle by integration.

    `free_cycle` composes the arcs exactly but its cycle search can stall on
    strongly overdamped three level shapes; this integrates instead, settles
    for ``n_settle`` estimated cycles and times successive maxima of ``x``.
    Agrees with `free_cycle` to about ``1e-6`` where both work.
    """
    from scipy.integrate import solve_ivp
    import section
    f = field(levels, edges)
    warm = solve_ivp(f, (0.0, 12.0*n_settle), [2.5, 0.0], method=section.METHOD,
                     rtol=rtol, atol=1e-11)

    def top(t, y):
        return y[1]
    top.direction = -1.0
    sol = solve_ivp(f, (0.0, 400.0), warm.y[:, -1], method=section.METHOD,
                    rtol=rtol, atol=1e-11, events=top, dense_output=True)
    te = sol.t_events[0]
    te = te[te > 1e-8]
    T = float((te[-1] - te[1])/(len(te) - 2))
    return float(sol.sol(te[1])[0]), T


def two_level_matched(z1, T=VDP_T, R=VDP_R):
    """The two level model with outer ratio ``z1`` matched to a free cycle.

    The period of the two level model depends on its damping ratios alone
    and the edge sets the amplitude exactly proportionally, so matching a
    period fixes ``zeta_0`` given ``zeta_1``, and matching the amplitude
    then fixes the edge. What is left is the one dimensional family in
    ``z1``: every member has the same free cycle and a different damping
    *shape*. Returns ``(levels, edges)``.
    """
    z0 = brentq(lambda z: free_cycle((z, z1), (1.0,))[1] - T, -4.0, -0.05,
                xtol=1e-6)
    r, _ = free_cycle((z0, z1), (1.0,))
    return (z0, z1), (R/r,)


def uniform_matched(n, mu=CMP_MU, xmax=CMP_XMAX, T=VDP_T, R=VDP_R):
    """The fitted ``n`` level staircase normalised the naive way.

    Every damping ratio is scaled by one factor to hit the period, then
    every edge by one factor to hit the amplitude — the two knobs a
    person would reach for. Returns ``(levels, edges, s, k)``.
    """
    lv, ed = vdp_staircase(mu, n, xmax)
    s = brentq(lambda f: free_cycle_num(tuple(f*z for z in lv), ed)[1] - T,
               0.3, 4.0, xtol=1e-6)
    lv2 = tuple(s*z for z in lv)
    r, _ = free_cycle_num(lv2, ed)
    k = R/r
    return lv2, tuple(k*e for e in ed), s, k


def three_level_matched(z2, mu=CMP_MU, xmax=CMP_XMAX, T=VDP_T, R=VDP_R):
    """A three level model with outer ratio ``z2`` matched to a free cycle.

    Keeps the fitted middle level and the edge ratio, solves the core
    level for the period and scales the edges for the amplitude. One of
    the two shape freedoms a three level model has after its free cycle
    is fixed. Returns ``(levels, edges)``.
    """
    lv, ed = vdp_staircase(mu, 3, xmax)
    z1 = lv[1]
    z0 = brentq(lambda z: free_cycle_num((z, z1, z2), ed)[1] - T, -5.0, -0.05,
                xtol=1e-6)
    r, _ = free_cycle_num((z0, z1, z2), ed)
    k = R/r
    return (z0, z1, z2), tuple(k*e for e in ed)


#: What `fit_bands` lands on with 20% leeway on the free cycle, every
#: parameter free, from the matched two level model and the uniformly
#: scaled three level one. About half an hour and three quarters of an hour
#: respectively; ``python3 staircase.py fit`` re-runs them.
TWO_FITTED = ((-1.2419249232482135, 8.328703618461768), (1.43646639597176,))
THREE_FITTED = ((-1.7350655010015574, 3.8359503860037094, 15.047130340757839),
                (1.0750014202405769, 1.9812442858436983))

#: The same fit to Van der Pol at ``mu = 1``: free cycle 2.0086 and 6.6633,
#: plateau targets the end of lock 1 at ratio 2.195 and the start of lock 3
#: at 2.555 at ``A = 5``, window ratio 1.8 to 3.3 at 0.05. About an hour.
MU1 = 1.0
VDP_R_MU1, VDP_T_MU1 = 2.0086, 6.66329
VDP_END1_MU1, VDP_START3_MU1 = 2.195, 2.555
THREE_FITTED_MU1 = ((-0.35783384494211046, 0.8653252168440303, 3.57309031296195),
                    (1.1597467884015344, 1.9836384057971845))


def fit_mu1(maxfev=60, workers=None, log=print):
    """Re-run the ``mu = 1`` fit from the uniformly matched sampled staircase."""
    import vanderpol
    wl = vanderpol.w_lc(MU1)
    lv, ed, _, _ = uniform_matched(3, mu=MU1, T=VDP_T_MU1, R=VDP_R_MU1)
    oms = tuple(np.round(np.arange(1.80, 3.301, 0.05)*wl, 5))
    return fit_bands(lv, ed, maxfev=maxfev, workers=workers, log=log,
                     targets=(VDP_END1_MU1*wl, VDP_START3_MU1*wl),
                     free=(VDP_R_MU1, VDP_T_MU1), amp=CMP_AMP, oms=oms,
                     locks=("lock1", "lock3"))

#: Van der Pol's period adding region at ``CMP_AMP``: the last frequency
#: locked 3:1 and the first locked 5:1 from there on, from the fine scan.
VDP_END3, VDP_START5 = 2.4275, 2.4975
COARSE_OMS = tuple(np.round(np.arange(2.30, 2.7001, 0.01), 4))


def _label_at(levels, edges, om, amp=CMP_AMP):
    return forced_label(field(levels, edges, amp, om), om)[0]


def band_edges(levels, edges, amp=CMP_AMP, workers=None, oms=COARSE_OMS,
               locks=("lock3", "lock5")):
    """Where the lower lock plateau ends and the upper one begins.

    A coarse sweep over ``oms``, then two bisections on each edge, so the
    edges are located to a quarter of the grid step. The lock 5 plateau is
    allowed quasi-periodic blips — Van der Pol's is unbroken, but a coarse
    model's can carry an isolated torus inside it — so its start is the
    first frequency from which nothing but lock 5 and tori follows. Returns
    ``(end3, start5, labels)``; an edge outside the window is reported one
    step beyond it. ``locks`` names the two plateaus, lock 3 and lock 5 by
    default.
    """
    import multiprocessing as mp
    lo_lock, hi_lock = locks
    args = [("s", om, (levels, edges), CMP_MU, amp) for om in oms]
    with mp.Pool(workers or mp.cpu_count()) as pool:
        out = pool.map(_scan_system, args, chunksize=1)
    labels = [lab for _, _, (lab, _) in sorted(out, key=lambda o: o[1])]
    oms = list(oms)
    step = oms[1] - oms[0]
    ok5 = (hi_lock, "torus")
    i = 0
    while i < len(labels) and labels[i] == lo_lock:
        i += 1
    j = len(labels)
    while j > 0 and labels[j - 1] in ok5:
        j -= 1

    def bisect(lo, hi, good):
        """``lo`` satisfies ``good``, ``hi`` does not; narrow twice."""
        for _ in range(2):
            mid = 0.5*(lo + hi)
            if good(_label_at(levels, edges, mid, amp)):
                lo = mid
            else:
                hi = mid
        return 0.5*(lo + hi)

    if i == 0:
        end3 = oms[0] - step
    elif i == len(labels):
        end3 = oms[-1] + step
    else:
        end3 = bisect(oms[i - 1], oms[i], lambda lab: lab == lo_lock)
    if j == len(labels):
        start5 = oms[-1] + step
    elif j == 0:
        start5 = oms[0] - step
    else:
        start5 = bisect(oms[j], oms[j - 1], lambda lab: lab in ok5)
    return end3, start5, labels


def fit_bands(levels, edges, leeway=0.2, maxfev=40, workers=None,
              log=print, targets=(VDP_END3, VDP_START5), free=(VDP_R, VDP_T),
              amp=CMP_AMP, oms=COARSE_OMS, locks=("lock3", "lock5")):
    """Move a coarse model's period adding region onto Van der Pol's.

    Every damping ratio and every edge is free. The objective is the
    squared distance of the two plateau edges from `VDP_END3` and
    `VDP_START5`, in units of the coarse step, plus a penalty once the
    free amplitude or period leaves ``leeway`` of Van der Pol's. Nelder-Mead
    on the parameters, with the positive ratios and the edges in log form.
    Each evaluation is a coarse sweep, about a minute.

    ``targets`` are the two plateau edges to hit, ``free`` the amplitude and
    period the leeway is measured against, ``amp`` the drive, ``oms`` the
    coarse window and ``locks`` the two plateaus; the defaults are Van der
    Pol at ``mu = 5``. Returns ``(levels, edges, end3, start5, r, T)`` at
    the best point seen.
    """
    from scipy.optimize import minimize
    t_end, t_start = targets
    R_free, T_free = free
    step = oms[1] - oms[0]
    n = len(levels)
    z0 = levels[0]
    pos = [z > 0 for z in levels[1:]]

    def unpack(p):
        lv = [p[0]]
        for k, ispos in enumerate(pos):
            lv.append(float(np.exp(p[1 + k])) if ispos else float(p[1 + k]))
        ed = tuple(float(np.exp(v)) for v in p[n:])
        return tuple(lv), ed

    p0 = [z0] + [np.log(z) if ispos else z for z, ispos in zip(levels[1:], pos)]
    p0 += [np.log(e) for e in edges]
    best = {"J": np.inf}

    def J(p):
        lv, ed = unpack(p)
        if any(ed[k] >= ed[k + 1] for k in range(len(ed) - 1)) or lv[0] >= 0:
            return 1e4
        try:
            r, T = free_cycle_num(lv, ed)
        except Exception:
            return 1e4
        dr, dT = abs(r/R_free - 1.0), abs(T/T_free - 1.0)
        pen = 1e3*(max(0.0, dr - leeway)**2 + max(0.0, dT - leeway)**2)
        end3, start5, labels = band_edges(lv, ed, amp=amp, workers=workers,
                                          oms=oms, locks=locks)
        val = ((end3 - t_end)/step)**2 + ((start5 - t_start)/step)**2 + pen
        log("  J=%9.3f  levels %s edges %s  r %.3f T %.3f  end3 %.4f start5 %.4f"
            % (val, tuple(round(z, 3) for z in lv),
               tuple(round(e, 3) for e in ed), r, T, end3, start5))
        if val < best["J"]:
            best.update(J=val, lv=lv, ed=ed, end3=end3, start5=start5, r=r, T=T)
        return val

    steps = [0.4] + [0.3]*(n - 1) + [0.15]*len(edges)
    simplex = [np.array(p0, float)]
    for k, st in enumerate(steps):
        q = np.array(p0, float)
        q[k] += st
        simplex.append(q)
    minimize(J, np.array(p0, float), method="Nelder-Mead",
             options=dict(initial_simplex=np.array(simplex), maxfev=maxfev,
                          xatol=1e-3, fatol=0.05))
    return best["lv"], best["ed"], best["end3"], best["start5"], best["r"], best["T"]


def _scan_system(args):
    """Worker for :func:`system_scan` (module level so it can be pickled)."""
    tag, om, sysd, mu, amp = args
    if sysd is None:
        import vanderpol
        flow = vanderpol.field(mu, amp, om)
    else:
        flow = field(sysd[0], sysd[1], amp, om)
    return tag, om, forced_label(flow, om)


def system_scan(oms, systems, mu=CMP_MU, amp=CMP_AMP, workers=None):
    """`window_scan` for explicit ``{tag: (levels, edges)}`` systems.

    A tag mapped to ``None`` is Van der Pol itself. Returns the same
    ``{tag: [(om, label, lam), ...]}`` shape as `window_scan`.
    """
    import multiprocessing as mp
    args = [(t, om, sysd, mu, amp) for t, sysd in systems.items()
            for om in oms]
    with mp.Pool(workers or mp.cpu_count()) as pool:
        out = pool.map(_scan_system, args, chunksize=1)
    res = {t: [] for t in systems}
    for tag, om, (lab, lam) in out:
        res[tag].append((om, lab, lam))
    for t in res:
        res[t].sort()
    return res


def normalise(workers=None):
    """Can a coarse staircase be normalised onto Van der Pol's regime map?

    Three families, each member matched to Van der Pol's free amplitude and
    period first: the two level model along its remaining shape freedom
    ``z1``; the fitted 2, 3 and 5 level staircases scaled uniformly; and
    the three level model along its outer ratio. Then the two models
    `fit_bands` produced with 20% leeway, `TWO_FITTED` and `THREE_FITTED`,
    and the plain fitted staircases for reference. Each is swept over
    ``NORM_OMS`` beside Van der Pol and scored by the Jaccard agreement of
    `window_agreement`. Prints the tables ``VANDERPOL.md`` quotes and returns
    ``(systems, scan)``.
    """
    systems = {"vdp": None}
    for n in (2, 3):
        systems["fitted %d" % n] = vdp_staircase(CMP_MU, n, CMP_XMAX)
    for z1 in (3.0, 5.0, 6.5, 7.25, 8.0, 9.0, 10.63, 12.0, 15.0, 25.0, 40.0):
        systems["two z1=%.2f" % z1] = two_level_matched(z1)
    for n in (2, 3, 5):
        lv, ed, s, k = uniform_matched(n)
        systems["uniform %d" % n] = (lv, ed)
        print("uniform %d: zeta scale %.4f, edge scale %.4f" % (n, s, k))
    for z2 in (6.0, 8.0, 10.0, 20.0):
        systems["three z2=%.0f" % z2] = three_level_matched(z2)
    systems["two fitted bands"] = TWO_FITTED
    systems["three fitted bands"] = THREE_FITTED
    for tag, sysd in systems.items():
        if sysd is not None:
            r, T = free_cycle(*sysd)
            print("  %-16s levels %s edges %s  r %.4f T %.3f"
                  % (tag, tuple(round(z, 3) for z in sysd[0]),
                     tuple(round(float(e), 3) for e in sysd[1]), r, T))
    scan = system_scan(NORM_OMS, systems, workers=workers)
    agree = window_agreement(scan)
    print("\n%-16s %8s %7s %8s" % ("system", "chaotic", "shared", "jaccard"))
    for tag in systems:
        n, sh, j = agree[tag]
        print("%-16s %8d %7d %8.3f" % (tag, n, sh, j))
    print()
    _print_scan(scan)
    return systems, scan


# --------------------------------------------------- the regime map
#: The control section's drive grid, refined: ratios of the drive frequency
#: to Van der Pol's free cycle frequency, and absolute drive amplitudes
#: with an unforced row.
MAP_R = tuple(np.round(np.arange(0.5, 8.001, 0.1), 2))
MAP_A = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0)


def _regime_point(args):
    """Worker for :func:`regime_map` (module level so it can be pickled)."""
    import section
    lv, ed, om, amp, n_skip = args
    lab, q, w, lam = section.classify(field(lv, ed, amp, om), 2.0*np.pi/om,
                                      [2.0, 0.0], n_skip)
    return lab, (q if q is not None else 0), w, (0.0 if lam is None else lam)


def regime_map(levels, edges, oms, amps, n_skip=CMP_NSKIP, workers=None):
    """Classify a staircase's driven response over a grid of ``(om, amp)``.

    The same measurement `vanderpol.regime_map` makes, with absolute drive
    frequencies so the two can share an axis. Returns
    ``(labels, q, w, lam)``, each of shape ``(len(amps), len(oms))``.
    """
    import multiprocessing as mp
    args = [(levels, edges, om, a, n_skip) for a in amps for om in oms]
    with mp.Pool(workers or mp.cpu_count()) as pool:
        out = pool.map(_regime_point, args, chunksize=1)
    sh = (len(amps), len(oms))
    return (np.array([o[0] for o in out], dtype=object).reshape(sh),
            np.array([o[1] for o in out]).reshape(sh),
            np.array([o[2] for o in out]).reshape(sh),
            np.array([o[3] for o in out]).reshape(sh))


def _confirm_point(args):
    """Worker for :func:`confirm_map`."""
    import section
    flow_kind, om, amp, n_skip, mu = args
    if flow_kind == "vdp":
        import vanderpol
        flow = vanderpol.field(mu, amp, om)
        _, y0 = vanderpol.cycle(mu)
        y0 = list(y0)
    else:
        flow = field(flow_kind[0], flow_kind[1], amp, om)
        y0 = [2.0, 0.0]
    ok, lams = section.confirm_chaos(flow, 2.0*np.pi/om, y0, n_skip)
    return ok, lams


def confirm_map(flow_kind, oms, amps, lab, lam, n_skip=CMP_NSKIP, workers=None,
                mu=CMP_MU):
    """Re-test every chaotic cell of a regime map with `section.confirm_chaos`.

    ``flow_kind`` is ``"vdp"`` or ``(levels, edges)``. Returns a boolean
    array the shape of the map, true where a chaotic verdict survived.
    """
    import multiprocessing as mp
    cells = [(i, j) for i in range(len(amps)) for j in range(len(oms))
             if lab[i, j] == "chaos"]
    args = [(flow_kind, oms[j], amps[i], n_skip, mu) for i, j in cells]
    with mp.Pool(workers or mp.cpu_count()) as pool:
        out = pool.map(_confirm_point, args, chunksize=1)
    keep = np.zeros(lab.shape, dtype=bool)
    for (i, j), (ok, _) in zip(cells, out):
        keep[i, j] = ok
    return keep


def regime_compare(workers=None, log=print, mu=CMP_MU, model=None):
    """The fitted three level model and Van der Pol over the whole drive grid.

    Absolute drive frequencies ``MAP_R`` times Van der Pol's free cycle
    frequency at ``mu``, amplitudes ``MAP_A`` including no drive at all.
    ``model`` is the ``(levels, edges)`` to compare, `THREE_FITTED` by
    default. Every chaotic verdict is re-tested. Returns a dict with both
    raw maps and the confirmation masks; ``figures.fig_regime_three`` draws
    it.
    """
    import time
    import vanderpol
    wl = vanderpol.w_lc(mu)
    oms = tuple(float(r*wl) for r in MAP_R)
    lv, ed = THREE_FITTED if model is None else model
    out = dict(oms=oms, ratios=MAP_R, amps=MAP_A, w_lc=wl, mu=mu, model=(lv, ed))
    t0 = time.time()
    out["three"] = regime_map(lv, ed, oms, MAP_A, workers=workers)
    log("three level map: %d cells, %.0fs, chaotic %d"
        % (len(oms)*len(MAP_A), time.time() - t0,
           int(np.sum(out["three"][0] == "chaos"))))
    t0 = time.time()
    out["vdp"] = vanderpol.regime_map(mu, np.array(MAP_R), np.array(MAP_A),
                                      workers=workers)
    log("Van der Pol map: %.0fs, chaotic %d"
        % (time.time() - t0, int(np.sum(out["vdp"][0] == "chaos"))))
    for tag, kind in (("three", (lv, ed)), ("vdp", "vdp")):
        t0 = time.time()
        lab, q, w, lam = out[tag]
        out[tag + "_ok"] = confirm_map(kind, oms, MAP_A, lab, lam,
                                       workers=workers, mu=mu)
        log("%s: %d of %d chaotic cells confirmed, %.0fs"
            % (tag, int(out[tag + "_ok"].sum()), int(np.sum(lab == "chaos")),
               time.time() - t0))
    return out


#: Transitions the coarse regime map cannot resolve, swept at 0.01 in the
#: ratio: the two chaotic transitions Van der Pol has at ``A = 10`` and its
#: second band at ``A = 5``.
TRANSITION_WINDOWS = ((10.0, 3.90, 4.40), (10.0, 7.00, 7.70), (5.0, 6.10, 6.70))


def _transition_point(args):
    """Worker for :func:`regime_transitions`: classify, confirm if chaotic."""
    import section
    import vanderpol
    tag, r, amp, wl = args
    om = r*wl
    lv, ed = THREE_FITTED
    if tag == "vdp":
        flow = vanderpol.field(CMP_MU, amp, om)
        y0 = list(vanderpol.cycle(CMP_MU)[1])
    else:
        flow = field(lv, ed, amp, om)
        y0 = [2.0, 0.0]
    lab, q, w, lam = section.classify(flow, 2.0*np.pi/om, y0, CMP_NSKIP)
    if lab == "chaos":
        ok, _ = section.confirm_chaos(flow, 2.0*np.pi/om, y0, CMP_NSKIP)
        lab = "chaos" if ok else "torus"
    return tag, amp, r, lab, lam


def regime_transitions(windows=TRANSITION_WINDOWS, workers=None):
    """Sweep the fitted three level model and Van der Pol across the
    transitions the coarse map steps over, at 0.01 in the ratio, with
    every chaotic verdict confirmed. Returns
    ``{(tag, amp): [(r, label, lam), ...]}`` and prints the runs.
    """
    import multiprocessing as mp
    import vanderpol
    wl = vanderpol.w_lc(CMP_MU)
    args = [(tag, float(r), amp, wl) for amp, lo, hi in windows
            for r in np.round(np.arange(lo, hi + 0.001, 0.01), 3)
            for tag in ("three", "vdp")]
    with mp.Pool(workers or mp.cpu_count()) as pool:
        out = pool.map(_transition_point, args, chunksize=1)
    res = {}
    for tag, amp, r, lab, lam in out:
        res.setdefault((tag, amp), []).append((r, lab, lam))
    for key in sorted(res):
        rows = sorted(res[key])
        res[key] = rows
        print("%s at A = %g: %d chaotic, confirmed" %
              (key[0], key[1], sum(1 for _, l, _ in rows if l == "chaos")))
        for lab, lo, hi, n in runs(rows):
            if n > 1 or lab == "chaos":
                print("    %-12s r %.2f - %.2f (%d)" % (lab, lo, hi, n))
    return res


def _print_scan(scan):
    for tag in scan:
        rows = scan[tag]
        nch = sum(1 for _, lab, _ in rows if lab == "chaos")
        print("  %5s: %d chaotic of %d; runs:" % (tag, nch, len(rows)))
        for lab, lo, hi, n in runs(rows):
            print("         %-6s %.3f - %.3f  (%d)" % (lab, lo, hi, n))


def check_reduces(zp=0.4, zm=-0.15, x0=1.0):
    """Confirm the averaging here is the two level equation already in use.

    `displacement.py` solves ``2 phi - sin 2 phi = pi rho`` for the
    symmetric displacement model. That is this module's ``mean_damping``
    with two levels, so both should give the same cycle radius.
    """
    import displacement
    here = cycles_predicted((zm, zp), (x0,))
    return (here[0][0] if here else np.nan,
            displacement.amplitude(zp, zm, x0, symmetric=True))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "sweep":
        level_floor()
        print()
        floor_crosschecks()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "normalise":
        normalise()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "regime":
        import pickle
        pickle.dump(regime_compare(), open("regime.pkl", "wb"))
        print()
        regime_transitions()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "fit1":
        print("landed on", fit_mu1())
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "regime1":
        import pickle
        pickle.dump(regime_compare(mu=MU1, model=THREE_FITTED_MU1),
                    open("regime_mu1.pkl", "wb"))
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "fit":
        for lv, ed in (two_level_matched(7.25), uniform_matched(3)[:2]):
            print("fitting from", lv, ed)
            print("  landed on", fit_bands(lv, ed))
        sys.exit(0)
    print("the two level case must reproduce displacement.py")
    a, b = check_reduces()
    print("  staircase %.10f   displacement.py %.10f   diff %.2e\n"
          % (a, b, abs(a - b)))

    print("three levels, bistable: levels %s edges %s"
          % (BISTABLE_LEVELS, BISTABLE_EDGES))
    print("  averaged prediction:")
    for r, st in cycles_predicted(BISTABLE_LEVELS, BISTABLE_EDGES):
        print("    r = %8.5f   %s" % (r, "stable" if st else "unstable"))
    print("  exact half cycle map:")
    for r, m, st in cycles_exact(BISTABLE_LEVELS, BISTABLE_EDGES):
        print("    r = %8.5f   multiplier %9.6f   %s   T = %.6f"
              % (r, m, "stable" if st else "unstable",
                 period(r, BISTABLE_LEVELS, BISTABLE_EDGES)))

    print("\nstaircase approximations to Van der Pol, mu = 1")
    print("  %3s %10s %10s" % ("n", "levels", "cycle r"))
    for n in (2, 3, 5, 9, 17):
        lv, ed = vdp_staircase(1.0, n)
        ex = cycles_exact(lv, ed)
        r = ex[-1][0] if ex else float("nan")
        print("  %3d %10d %10.5f" % (n, len(lv), r))
