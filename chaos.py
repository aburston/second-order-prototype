"""The chaotic region of the three level prototype in its nearly harmonic mode.

Run ``python3 chaos.py`` to print every number ``CHAOS.md`` quotes and to
write its figures into ``figures/`` as ``chaos-mu1-*.png``, light and dark.
``python3 chaos.py quick`` skips the drive sweep altogether and draws only
the figures that do not need it. The sweep is cached in ``figures/`` after
its first run, so re-drawing a figure costs a minute rather than a quarter
of an hour; ``python3 chaos.py fresh`` recomputes it.

``THREELEVEL.md`` fits the prototype to Van der Pol at two values of the
relaxation parameter. At ``mu = 5`` the oscillator is a relaxation
oscillator and its chaos sits in the transitions of a period adding
sequence; ``STROBOSCOPIC.md`` maps that case. At ``mu = 1`` -- Van der Pol
in its **ordinary oscillating mode**, a free cycle that is nearly a
sinusoid -- there is no period adding at all: two tongues, tori almost
everywhere else, and one narrow band of chaos reached by period doubling of
the 1:1 lock at a drive near **half** the cycle frequency. That band is
what this module is about.

What it establishes, in order:

1. **Where the chaos is.** A sweep of drive ratio at ``A = 1``, the drive
   strength at which the band exists, classifying every cell as a lock, a
   torus or chaos by ``section.py``, and re-testing every chaotic verdict
   with ``section.confirm_chaos``. The prototype's band and Van der Pol's
   are both found, and they do not sit on top of one another.

2. **What it looks like in time and in the plane.** The time series, the
   phase plane orbit, and the two together saying why the response is not
   periodic when the drive is.

3. **The stroboscopic map.** The phase plane sampled once per drive period
   is the section, exactly as in ``STROBOSCOPIC.md``, and the samples are
   the attractor. Each sample is coloured by the stretch its own step
   applies -- the largest singular value of the one period Jacobian, exact
   for the prototype from ``maps.strobe_step`` and from the variational
   equation for Van der Pol -- so the figure shows not just where the
   attractor is but where the map pulls it apart.

4. **The fold.** A short segment of initial conditions laid across the
   attractor and carried forward three drive periods, coloured by position
   along the segment. It stretches and folds back onto the attractor,
   which is the mechanism the positive exponent measures.

5. **The comparison.** Every one of those on Van der Pol at the same
   ``mu``, the real oscillator the prototype was fitted to, with the
   attractor's Lyapunov exponent, its area contraction per drive period,
   its Kaplan-Yorke dimension and a box counting dimension for both.

Nothing here is fitted. The parameters are ``staircase.THREE_FITTED_MU1``,
fitted in ``THREELEVEL.md`` to two lock plateau edges at a different drive
strength (``A = 5``) and a different part of the frequency range; the
chaotic band is a prediction of that fit, and this module measures how good
a prediction it is.
"""
import json
import os
import sys
import time

import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp

import maps
import section
import staircase
import vanderpol

# ------------------------------------------------------------- the systems
#: Van der Pol in its ordinary oscillating mode, and the three level
#: prototype fitted to it there.
MU = 1.0
LEVELS, EDGES = staircase.THREE_FITTED_MU1

#: The drive ratio is measured against Van der Pol's own free cycle
#: frequency for both systems, as everywhere else in the repository, so a
#: ratio means the same drive frequency on each.
W_LC = vanderpol.w_lc(MU)

#: The drive strength at which this mode has chaos. ``THREELEVEL.md``'s
#: grid found Van der Pol's one chaotic cell at ``A = 1``; at ``A = 5`` the
#: 1:1 tongue is wide enough to swallow the subharmonic region on both
#: systems.
AMP = 1.0

#: The sweep. Fine enough to resolve bands two or three cells wide, which
#: is what these are: the 0.1 grid of ``THREELEVEL.md``'s regime maps steps
#: straight over both of them. The subharmonic window is swept five times
#: finer again, because the cascade that ends in the chaos happens inside
#: two cells of the coarse grid.
BAND_RS = tuple(np.unique(np.round(np.concatenate(
    (np.arange(0.40, 0.8001, 0.005), np.arange(0.46, 0.6001, 0.001))), 4)))

#: Windows the cascade is drawn on, one per system, chosen from the sweep.
CASCADE = {"proto": (0.550, 0.580), "vdp": (0.480, 0.515)}

#: Drive ratios worked in detail, one inside each system's own band. They
#: are not the same number, which is the point of the comparison.
R_PROTO, R_VDP = 0.565, 0.490

#: Stroboscopic samples: a transient discarded, then the attractor.
N_SKIP, N_ATTRACTOR = 400, 12000

#: Orbit drawn in the phase plane and in time, in drive periods, and how
#: finely it is sampled within each.
N_TRACE, N_TIME, TRACE_PER = 40, 20, 400

MODEL = maps.staircase_model(LEVELS, EDGES)


def om_of(r):
    """Drive frequency at drive ratio ``r``."""
    return r*W_LC


def flow(kind, r):
    """Vector field of either system at drive ratio ``r``."""
    om = om_of(r)
    if kind == "vdp":
        return vanderpol.field(MU, AMP, om)
    return staircase.field(LEVELS, EDGES, AMP, om)


def td_of(r):
    """Drive period at drive ratio ``r``."""
    return 2.0*np.pi/om_of(r)


# --------------------------------------------------- 1. where the chaos is
def _band_point(args):
    """Worker for :func:`band_scan`: classify one cell and keep its section."""
    kind, r = args
    td = td_of(r)
    f = flow(kind, r)
    pts = section.strobe(f, td, [2.0, 0.0], N_SKIP, n_keep=80)
    q = section.lock_order(pts)
    if q is not None:
        return kind, r, "lock%d" % q, 0.0, pts[:, 0]
    lam = section.lyapunov(f, td, [2.0, 0.0], N_SKIP//2, n=400)
    lab = "chaos" if lam > section.LAM_TOL else "torus"
    return kind, r, lab, float(lam), pts[:, 0]


#: Where the sweep is cached. It is eleven minutes of arithmetic that does
#: not change when a figure does, so it is written once and reused; delete
#: the file, or pass ``fresh``, to recompute it.
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "figures", ".chaos-scan.json")


def band_scan(rs=BAND_RS, workers=4, cache=CACHE, fresh=False):
    """Classify both systems across the drive sweep.

    Returns ``{kind: [(r, label, lam, xs), ...]}`` with ``xs`` the
    stroboscopic displacements, which are what the orbit diagram draws.
    The cache is keyed on the grid, so changing ``BAND_RS`` recomputes.
    """
    key = "%.4f:%.4f:%d" % (rs[0], rs[-1], len(rs))
    if cache and not fresh and os.path.exists(cache):
        with open(cache) as f:
            got = json.load(f)
        if got.get("key") == key:
            print("  (sweep read from %s)" % os.path.basename(cache))
            return {k: [(r, lab, lam, np.array(xs))
                        for r, lab, lam, xs in rows]
                    for k, rows in got["scan"].items()}
    jobs = [(k, float(r)) for k in ("proto", "vdp") for r in rs]
    with Pool(workers) as p:
        res = p.map(_band_point, jobs)
    out = {"proto": [], "vdp": []}
    for kind, r, lab, lam, xs in res:
        out[kind].append((r, lab, lam, xs))
    for v in out.values():
        v.sort(key=lambda row: row[0])
    if cache:
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        with open(cache, "w") as f:
            json.dump({"key": key,
                       "scan": {k: [(r, lab, lam, list(map(float, xs)))
                                    for r, lab, lam, xs in rows]
                                for k, rows in out.items()}}, f)
    return out


def _confirm_point(args):
    """Worker for :func:`confirm_band`."""
    kind, r = args
    ok, lams = section.confirm_chaos(flow(kind, r), td_of(r), [2.0, 0.0],
                                     N_SKIP)
    return kind, r, bool(ok), tuple(float(l) for l in lams)


def confirm_band(scan, workers=4, cache=None, fresh=False):
    """Re-test every chaotic cell at five times the run and a wider twin.

    ``section.confirm_chaos``'s own docstring is the reason this is not
    optional: four of seven marginal cells changed sign under it when the
    threshold was first set.
    """
    jobs = [(k, r) for k, rows in scan.items()
            for r, lab, _, _ in rows if lab == "chaos"]
    if not jobs:
        return []
    if cache and not fresh and os.path.exists(cache):
        with open(cache) as f:
            got = json.load(f)
        if got.get("jobs") == [[k, r] for k, r in jobs]:
            return [(k, r, ok, tuple(l)) for k, r, ok, l in got["res"]]
    with Pool(workers) as p:
        res = p.map(_confirm_point, jobs)
    if cache:
        with open(cache, "w") as f:
            json.dump({"jobs": [[k, r] for k, r in jobs],
                       "res": [(k, r, ok, list(l)) for k, r, ok, l in res]},
                      f)
    return res


def apply_confirmation(scan, confirmed):
    """Downgrade every chaotic cell that fails the convergence test.

    The first pass is one exponent on one run length. A cell that does not
    hold up under ``section.confirm_chaos`` is not chaotic here, and the
    figures and the band edges are drawn from the labels *after* this, not
    before, so nothing in the document rests on an unconfirmed cell.

    Returns the number of cells downgraded.
    """
    bad = {(k, r) for k, r, ok, _ in confirmed if not ok}
    n = 0
    for kind, rows in scan.items():
        for i, (r, lab, lam, xs) in enumerate(rows):
            if (kind, r) in bad:
                rows[i] = (r, "torus", lam, xs)
                n += 1
    return n


def bands(rows):
    """Contiguous runs of chaotic cells in one system's sweep."""
    out, run = [], []
    for r, lab, lam, _ in rows:
        if lab == "chaos":
            run.append(r)
        elif run:
            out.append((run[0], run[-1]))
            run = []
    if run:
        out.append((run[0], run[-1]))
    return out


# ----------------------------------------------- 3. the stroboscopic map
def proto_attractor(r, n=N_ATTRACTOR, n_skip=N_SKIP):
    """Stroboscopic samples of the prototype, with the map's own Jacobian.

    Uses the exact map of ``maps.py`` -- arcs composed analytically with a
    saltation matrix at every wall the orbit crosses -- rather than
    integrating the discontinuous field, so the Jacobian is the map's and
    not a difference of trajectories.

    Returns ``(Y, stretch, logdet)``: the samples, the largest singular
    value of each step's Jacobian, and ``log|det J|`` of each step.
    """
    om = om_of(r)
    y = np.array([2.0, 0.0])
    for _ in range(n_skip):
        y, _, _ = maps.strobe_step(MODEL, y, 0.0, AMP, om)
    Y = np.empty((n, 2))
    s = np.empty(n)
    ld = np.empty(n)
    for i in range(n):
        Y[i] = y
        y, _, J = maps.strobe_step(MODEL, y, 0.0, AMP, om)
        s[i] = np.linalg.svd(J, compute_uv=False)[0]
        ld[i] = np.log(abs(np.linalg.det(J)))
    return Y, s, ld


def _vdp_variational(r):
    """Right hand side of Van der Pol with its 2x2 fundamental matrix."""
    om = om_of(r)

    def f(t, z):
        x, v = z[0], z[1]
        g = MU*(1.0 - x*x)
        a, b, c, d = z[2], z[3], z[4], z[5]
        # d/dt [dx/dy0] = [[0, 1], [-1 - 2 mu x v, g]] [dx/dy0]
        k = -1.0 - 2.0*MU*x*v
        return [v, -x + g*v + AMP*np.cos(om*t),
                c, d, k*a + g*c, k*b + g*d]
    return f


def vdp_attractor(r, n=N_ATTRACTOR, n_skip=N_SKIP):
    """The same three things for Van der Pol, from the variational equation.

    The field is smooth, so the Jacobian is integrated alongside the state
    rather than composed: same quantity, the only method each system admits.
    """
    f = _vdp_variational(r)
    td = td_of(r)
    y = np.array([2.0, 0.0])
    z = np.empty(6)

    def step(y):
        z[:2], z[2:] = y, (1.0, 0.0, 0.0, 1.0)
        sol = solve_ivp(f, (0.0, td), z, method=section.METHOD,
                        rtol=section.RTOL, atol=section.ATOL)
        w = sol.y[:, -1]
        return w[:2], w[2:].reshape(2, 2)

    for _ in range(n_skip):
        y, _ = step(y)
    Y = np.empty((n, 2))
    s = np.empty(n)
    ld = np.empty(n)
    for i in range(n):
        Y[i] = y
        y, J = step(y)
        s[i] = np.linalg.svd(J, compute_uv=False)[0]
        ld[i] = np.log(abs(np.linalg.det(J)))
    return Y, s, ld


def _attractor_job(args):
    kind, r = args
    return (proto_attractor if kind == "proto" else vdp_attractor)(r)


def attractors(workers=2):
    """Both attractors, each in its own process."""
    jobs = [("proto", R_PROTO), ("vdp", R_VDP)]
    with Pool(workers) as p:
        got = p.map(_attractor_job, jobs)
    return dict(zip(("proto", "vdp"), got))


# -------------------------------------------------------- what they measure
def box_dimension(Y, ks=range(2, 9)):
    """Box counting dimension of a point set, with the range it was fitted on.

    The points are scaled to their own bounding square first, so the answer
    does not depend on the units of ``x`` against ``xdot``, and boxes are
    counted at ``2^-k`` of that square. The slope is fitted only where the
    count is neither saturated by the box grid nor by the sample -- between
    ten boxes and a twentieth of the sample size -- because outside that
    range a finite point cloud reports the dimension of its own sampling.

    Returns ``(dim, ks_used, counts)``.
    """
    p = Y - Y.min(axis=0)
    p = p/max(p.max(), 1e-300)
    ks = list(ks)
    counts = []
    for k in ks:
        n = 2**k
        idx = np.minimum((p*n).astype(int), n - 1)
        counts.append(len(set(map(tuple, idx))))
    counts = np.array(counts, float)
    ok = (counts > 10) & (counts < len(Y)/20.0)
    if ok.sum() < 2:
        ok = np.ones(len(ks), bool)
    a = np.polyfit(np.log(2.0**np.array(ks, float)[ok]), np.log(counts[ok]), 1)
    return float(a[0]), [k for k, o in zip(ks, ok) if o], counts


def kaplan_yorke(lam1, logdet, td):
    """Dimension of the attractor from the two exponents.

    ``lam1`` is measured by twin trajectories; the *sum* of the exponents is
    the mean area contraction per unit time, which is known exactly from
    the map's own determinant, so the second exponent follows without a
    second twin. For a planar map with ``lam1 > 0 > lam2`` the Kaplan-Yorke
    dimension is ``1 + lam1/|lam2|``.

    Returns ``(lam2, dky, lam_sum)``.
    """
    lam_sum = float(np.mean(logdet))/td
    lam2 = lam_sum - lam1
    dky = 1.0 + lam1/abs(lam2) if lam2 < 0.0 < lam1 else float("nan")
    return lam2, dky, lam_sum


def edge_jump(r=R_PROTO, v=None, eps=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6), Y=None):
    """How the map behaves across the inner zone edge, either side of ``a``.

    A point starting just inside the core and one starting just outside it
    are in different zones from the first instant, so the two orbits differ
    for the whole period. The question is what survives as the two starts
    close up: the *state* after a period does -- the map is continuous --
    but the *derivative* need not, and here it does not.

    This is a different edge from the grazing edges ``STROBOSCOPIC.md``
    measures, which are where an orbit just touches a wall on its way round
    and across which the map is continuously differentiable. This one is
    the wall itself, in the section.

    Returns ``[(eps, |dy|, |dJ|/|J|), ...]``. ``v`` defaults to the
    attractor's own strand where it crosses the edge.
    """
    a = EDGES[0]
    om = om_of(r)
    if v is None:
        Y = proto_attractor(r, n=2000)[0] if Y is None else Y
        m = ((Y[:, 1] > -0.30) & (Y[:, 1] < -0.16)
             & (np.abs(Y[:, 0] - a) < 0.01))
        v = float(np.median(Y[m, 1]))
    out = []
    for e in eps:
        ym, _, Jm = maps.strobe_step(MODEL, [a - e, v], 0.0, AMP, om)
        yp, _, Jp = maps.strobe_step(MODEL, [a + e, v], 0.0, AMP, om)
        out.append((e, float(np.hypot(*(yp - ym))),
                    float(np.linalg.norm(Jp - Jm)/np.linalg.norm(Jm))))
    return v, out


# ------------------------------------------------------------ 4. the fold
def fold_images(kind, r, y0, y1, n_pts=600, n_img=3):
    """A segment of initial conditions and its first few images.

    Returns a list of ``(n_pts, 2)`` arrays: the segment, then its image
    after one drive period, two, and so on. Colouring each by position
    along the original segment shows the stretch and the fold directly.
    """
    s = np.linspace(0.0, 1.0, n_pts)[:, None]
    P = np.asarray(y0, float) + s*(np.asarray(y1, float) - np.asarray(y0))
    out = [P.copy()]
    om = om_of(r)
    td = td_of(r)
    if kind == "proto":
        for _ in range(n_img):
            P = np.array([maps.strobe_step(MODEL, y, 0.0, AMP, om)[0]
                          for y in P])
            out.append(P.copy())
    else:
        f = vanderpol.field(MU, AMP, om)
        for _ in range(n_img):
            P = np.array([solve_ivp(f, (0.0, td), y, method=section.METHOD,
                                    rtol=section.RTOL,
                                    atol=section.ATOL).y[:, -1] for y in P])
            out.append(P.copy())
    return out


# ------------------------------------------------------------ orbits to draw
def trace(kind, r, n_periods, per=400, y0=None):
    """A continuous stretch of orbit on the attractor.

    Returns ``(t, x, v)`` with ``t`` in drive periods.
    """
    td = td_of(r)
    f = flow(kind, r)
    if y0 is None:
        y0 = section.strobe(f, td, [2.0, 0.0], N_SKIP, n_keep=1)[0]
    t = np.linspace(0.0, n_periods*td, per*n_periods + 1)
    y = section.run(f, t, list(y0))
    return t/td, y[:, 0], y[:, 1]


def strobe_pts(kind, r, n, n_skip=N_SKIP):
    """Stroboscopic samples by integration, the section engine's own route."""
    return section.strobe(flow(kind, r), td_of(r), [2.0, 0.0], n_skip,
                          n_keep=n)


def check_map_against_integration(r=R_PROTO, n=40):
    """The exact map and the integrator, on the same points.

    ``STROBOSCOPIC.md`` makes this check at ``mu = 5``; it is repeated here
    because the parameters are different ones and the map is what every
    Jacobian in this module comes from.

    Returns the largest disagreement over ``n`` steps, relative to the
    orbit's size.
    """
    om = om_of(r)
    td = td_of(r)
    f = flow("proto", r)
    y = np.array([2.0, 0.0])
    for _ in range(N_SKIP):
        y, _, _ = maps.strobe_step(MODEL, y, 0.0, AMP, om)
    worst, scale = 0.0, 0.0
    for _ in range(n):
        ym, _, _ = maps.strobe_step(MODEL, y, 0.0, AMP, om)
        yi = solve_ivp(f, (0.0, td), y, method=section.METHOD,
                       rtol=1e-11, atol=1e-13).y[:, -1]
        worst = max(worst, float(np.hypot(*(ym - yi))))
        scale = max(scale, float(np.max(np.abs(ym))))
        y = ym
    return worst/scale


# ---------------------------------------------------------------- figures
SYS_LABEL = {
    "proto": "three level prototype, $\\zeta = (%.2f, %.2f, %.2f)$, "
             "edges $(%.2f, %.2f)$" % (LEVELS + EDGES),
    "vdp": "Van der Pol, $\\mu = %g$" % MU,
}
SYS_SHORT = {"proto": "prototype", "vdp": "Van der Pol"}


def _ramp(th, slot=0):
    """A sequential ramp in one hue: the series colour, lightened and darkened.

    Magnitude gets a single hue running light to dark, never a rainbow and
    never two hues. The pale end is the series colour mixed most of the way
    into the surface, so it still reads against it, and the dark end is the
    same colour mixed towards ink.
    """
    import matplotlib
    from matplotlib.colors import to_rgb

    def mix(a, b, f):
        return tuple((1.0 - f)*np.array(to_rgb(a)) + f*np.array(to_rgb(b)))

    c = th["series"][slot]
    return matplotlib.colors.LinearSegmentedColormap.from_list(
        "ramp", [mix(c, th["surface"], 0.78), c, mix(c, th["ink"], 0.55)])


def _zones(ax, th, kind):
    """The prototype's zone edges as chrome; nothing for Van der Pol."""
    if kind != "proto":
        return
    for e in EDGES:
        for xv in (e, -e):
            ax.axvline(xv, color=th["ink2"], linewidth=1.0,
                       linestyle=(0, (5, 3)), zorder=4)


def _square(ax, pts_list, margin=0.12, half=None, centre=None):
    """Set square limits around the given point sets.

    Panels that are compared with one another are given the same width and
    each is centred on its own data, so a shape read off one panel is the
    same size as the shape read off the other. Returns the half width used.
    """
    P = np.vstack([np.asarray(p, float).reshape(-1, 2) for p in pts_list])
    c = np.asarray(centre, float) if centre is not None else \
        0.5*(P.max(axis=0) + P.min(axis=0))
    if half is None:
        half = 0.5*float(np.max(P.max(axis=0) - P.min(axis=0)))*(1.0 + margin)
    ax.set_xlim(c[0] - half, c[0] + half)
    ax.set_ylim(c[1] - half, c[1] + half)
    ax.set_aspect("equal", adjustable="box")
    return half


def fig_band(th, name, scan):
    """Where the chaos is: the orbit diagram, the exponent, and the cascade.

    Top row: every stroboscopic displacement of the settled response
    against drive ratio -- a finite set of points on a lock, a smear on a
    torus or in chaos. Middle row: the largest Lyapunov exponent, with
    ``section.LAM_TOL`` drawn as chrome. Bottom row: the same orbit diagram
    over each system's own subharmonic window, where the period doubling
    that ends in the chaos is resolved. Chaotic cells are shaded in every
    row that shows them.
    """
    from figures import newfig, style, legend, save
    fig, axes = newfig(th, 3, 2, figsize=(12.4, 9.4),
                       gridspec_kw=dict(height_ratios=(2.0, 1.1, 1.6)))
    for j, kind in enumerate(("proto", "vdp")):
        rows = scan[kind]
        top, mid, bot = axes[0][j], axes[1][j], axes[2][j]
        chaotic = bands(rows)
        for lo, hi in chaotic:
            for ax in (top, mid, bot):
                ax.axvspan(lo - 0.0005, hi + 0.0005, color=th["series"][1],
                           alpha=0.18, zorder=1, linewidth=0)
        for ax, size in ((top, 1.0), (bot, 1.8)):
            for r, lab, lam, xs in rows:
                col = th["series"][2] if lab == "chaos" else th["ink2"]
                ax.plot(np.full(len(xs), r), xs, ".", color=col,
                        markersize=size, alpha=0.8, zorder=3)
        rs = [r for r, _, _, _ in rows]
        lams = [lam for _, _, lam, _ in rows]
        mid.plot(rs, lams, "-", color=th["series"][0], linewidth=1.0,
                 zorder=3)
        mid.axhline(section.LAM_TOL, color=th["ink2"], linewidth=1.0,
                    linestyle=(0, (5, 3)), zorder=4)
        mid.text(0.995, section.LAM_TOL, "chaos threshold  ",
                 transform=mid.get_yaxis_transform(), fontsize=8,
                 color=th["ink2"], va="bottom", ha="right", zorder=5)
        mid.axhline(0.0, color=th["axis"], linewidth=0.8, zorder=2)
        mid.text(0.01, 0.06, "not computed on a lock, drawn as zero there",
                 transform=mid.transAxes, fontsize=7.5, color=th["ink2"],
                 zorder=5)
        r_mark = R_PROTO if kind == "proto" else R_VDP
        for ax in (top, bot):
            ax.axvline(r_mark, color=th["ink"], linewidth=1.0, zorder=5)
        bot.annotate("worked below, $r = %.3f$" % r_mark,
                     xy=(r_mark, 0.03), xycoords=("data", "axes fraction"),
                     xytext=(5, 0), textcoords="offset points", fontsize=8,
                     color=th["ink"], zorder=6)
        top.plot([], [], ".", color=th["series"][2], markersize=6,
                 label="chaotic cell, confirmed")
        top.plot([], [], ".", color=th["ink2"], markersize=6,
                 label="lock or torus")
        style(top, th, "", "$x$ at the strobe", SYS_LABEL[kind])
        style(mid, th, "drive ratio $r = \\Omega/\\omega_{lc}$",
              "$\\lambda$")
        style(bot, th, "drive ratio $r$", "$x$ at the strobe",
              "the subharmonic window, swept at 0.001")
        legend(top, th, loc="upper right", markerscale=2)
        top.set_xlim(BAND_RS[0], BAND_RS[-1])
        mid.set_xlim(BAND_RS[0], BAND_RS[-1])
        bot.set_xlim(*CASCADE[kind])
        top.set_ylim(-3.0, 3.0)
        mid.set_ylim(-0.03, 0.07)
        lo, hi = CASCADE[kind]
        inside = [x for r, _, _, xs in rows if lo <= r <= hi for x in xs]
        pad = 0.06*(max(inside) - min(inside))
        bot.set_ylim(min(inside) - pad, max(inside) + pad)
    fig.suptitle("The chaotic region in the nearly harmonic mode: "
                 "$\\mu = 1$, $A = %g$, the drive swept through half the "
                 "free cycle frequency" % AMP, color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "chaos-mu1-band")


def fig_time(th, name, series):
    """The time series, with the instants the section is cut at marked.

    The drive is periodic and drawn as chrome underneath; the response is
    not, and the stroboscopic samples are the response read at the drive's
    own period, which is what makes the section a section.
    """
    from figures import newfig, style, legend, save
    fig, axes = newfig(th, 2, 1, figsize=(12.0, 6.2), sharex=True)
    for ax, kind in zip(axes, ("proto", "vdp")):
        t, x, _, per = series[kind]
        n = int(np.searchsorted(t, N_TIME)) + 1
        t, x = t[:n], x[:n]
        # The samples are read off this same trajectory rather than
        # recomputed: two integrations of a chaotic orbit from the same
        # start diverge, so a separately strobed run would put the dots on
        # the attractor but not on this curve.
        cut = np.arange(0, n, per)
        r = R_PROTO if kind == "proto" else R_VDP
        ax.plot(t, x, color=th["series"][0], linewidth=1.1, zorder=3,
                label="$x(t)$")
        ax.plot(t[cut], x[cut], "o",
                color=th["series"][1], markersize=4.0, zorder=4,
                label="stroboscopic samples, one per drive period")
        xm = float(np.max(np.abs(x)))
        drive = np.cos(2.0*np.pi*t)
        ax.plot(t, 0.35*xm*drive - 1.45*xm, color=th["ink2"], linewidth=0.8,
                alpha=0.7, zorder=2, label="drive $A\\cos\\Omega t$, offset")
        for e in EDGES if kind == "proto" else ():
            for xv in (e, -e):
                ax.axhline(xv, color=th["ink2"], linewidth=0.9,
                           linestyle=(0, (5, 3)), zorder=2)
        style(ax, th, "", "$x$",
              "%s at $r = %.3f$" % (SYS_LABEL[kind], r))
        ax.set_ylim(-1.95*xm, 1.5*xm)
        legend(ax, th, loc="upper right", ncol=3)
    style(axes[1], th, "time, in drive periods", "$x$")
    fig.suptitle("Aperiodic response to a periodic drive: $A = %g$, the "
                 "drive near half the free cycle frequency" % AMP,
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "chaos-mu1-time")


def fig_phase(th, name, series, att):
    """The phase plane: the orbit, and the section it is cut on."""
    from figures import newfig, style, legend, save
    fig, axes = newfig(th, 1, 2, figsize=(12.0, 5.6))
    half = 1.1*max(float(np.max(np.abs(np.column_stack(
        (series[k][1], series[k][2]))))) for k in series)
    for ax, kind in zip(axes, ("proto", "vdp")):
        _, x, v, _ = series[kind]
        Y = att[kind][0]
        r = R_PROTO if kind == "proto" else R_VDP
        ax.plot(x, v, color=th["series"][0], linewidth=0.45, alpha=0.6,
                zorder=2, label="orbit, %d drive periods" % N_TRACE)
        ax.plot(Y[:, 0], Y[:, 1], ".", color=th["series"][1], markersize=1.3,
                zorder=3, label="stroboscopic samples, %d" % len(Y))
        _zones(ax, th, kind)
        style(ax, th, "$x$", "$\\dot{x}$",
              "%s, $r = %.3f$" % (SYS_LABEL[kind], r))
        _square(ax, [np.column_stack((x, v))], half=half,
                centre=(0.0, 0.0))
        legend(ax, th, loc="upper left", markerscale=4)
    fig.suptitle("Phase plane and the section cut through it, $A = %g$"
                 % AMP, color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "chaos-mu1-phase")


def fig_section(th, name, att, zooms):
    """The stroboscopic map, coloured by the stretch each step applies.

    Colour is a magnitude on one sequential ramp -- the largest singular
    value of the one period Jacobian at that sample, in logarithm -- and the
    two systems share the scale, so a colour means the same stretch on
    both. Identity is carried by the panels, not by hue. The lower row is
    the boxed window of the upper, at which scale the attractor is still
    made of strands.
    """
    import matplotlib
    import matplotlib.patches
    from figures import newfig, style, legend, save
    lo = min(np.log(att[k][1]).min() for k in att)
    hi = max(np.percentile(np.log(att[k][1]), 99.5) for k in att)
    norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
    ramp = _ramp(th)
    fig, axes = newfig(th, 2, 2, figsize=(12.2, 9.4))
    half = 0.5*max(float(np.max(np.ptp(att[k][0], axis=0))) for k in att)*1.35
    sc = None
    for j, kind in enumerate(("proto", "vdp")):
        Y, s, _ = att[kind]
        r = R_PROTO if kind == "proto" else R_VDP
        box = zooms[kind]
        for i, ax in enumerate((axes[0][j], axes[1][j])):
            sc = ax.scatter(Y[:, 0], Y[:, 1], c=np.log(s), cmap=ramp,
                            norm=norm, s=2.6 if i else 1.6, linewidths=0,
                            zorder=3, rasterized=True)
            _zones(ax, th, kind)
            if i == 0:
                ax.add_patch(matplotlib.patches.Rectangle(
                    (box[0], box[2]), box[1] - box[0], box[3] - box[2],
                    fill=False, edgecolor=th["ink"], linewidth=1.0,
                    zorder=5))
                _square(ax, [Y], half=half,
                        centre=0.5*(Y.max(axis=0) + Y.min(axis=0)))
                style(ax, th, "$x$", "$\\dot{x}$",
                      "%s, $r = %.3f$" % (SYS_SHORT[kind], r))
                ax.plot([], [], "s", color=th["ink"], markerfacecolor="none",
                        markersize=6, label="window below")
                legend(ax, th, loc="upper left")
            else:
                ax.set_xlim(box[0], box[1])
                ax.set_ylim(box[2], box[3])
                ax.set_aspect("equal", adjustable="box")
                style(ax, th, "$x$", "$\\dot{x}$",
                      "the same samples, %.3f wide" % (box[1] - box[0]))
    cb = fig.colorbar(sc, ax=axes.ravel().tolist(), pad=0.02, fraction=0.03)
    cb.set_label("$\\log$ of the step's largest singular value",
                 color=th["ink2"], fontsize=9)
    cb.ax.tick_params(colors=th["ink2"], labelsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("The stroboscopic map: %d samples, each coloured by the "
                 "stretch its own step applies" % len(att["proto"][0]),
                 color=th["ink"], fontsize=11)
    save(fig, name, "chaos-mu1-section")


def fig_fold(th, name, folds, att):
    """Stretch and fold: one segment of the plane, carried three periods.

    One panel per image rather than all four on top of one another, so
    which image a curve is comes from where the panel sits and not from a
    colour: colour inside a panel is position along the *original*
    segment, a magnitude on one sequential ramp, and it is what shows the
    order being preserved as the segment is drawn out and doubled back.
    The attractor is drawn behind each panel as chrome.
    """
    from figures import newfig, style, legend, save
    ramp = _ramp(th, 1)
    n_img = len(folds["proto"])
    fig, axes = newfig(th, 2, n_img, figsize=(3.1*n_img, 6.8))
    for i, kind in enumerate(("proto", "vdp")):
        imgs = folds[kind]
        Y = att[kind][0]
        s_col = np.linspace(0.0, 1.0, len(imgs[0]))
        half = 0.5*float(np.max(np.vstack(imgs + [Y]).max(axis=0)
                                - np.vstack(imgs + [Y]).min(axis=0)))*1.15
        centre = 0.5*(np.vstack(imgs + [Y]).max(axis=0)
                      + np.vstack(imgs + [Y]).min(axis=0))
        for k, P in enumerate(imgs):
            ax = axes[i][k]
            ax.plot(Y[:, 0], Y[:, 1], ".", color=th["axis"], markersize=0.9,
                    zorder=2, label="the attractor")
            ax.scatter(P[:, 0], P[:, 1], c=s_col, cmap=ramp, s=5.0,
                       linewidths=0, zorder=4)
            _zones(ax, th, kind)
            _square(ax, [P, Y], half=half, centre=centre)
            style(ax, th, "$x$" if i else "", "$\\dot{x}$" if k == 0 else "",
                  "%s, %d drive period%s" % (SYS_SHORT[kind], k,
                                             "" if k == 1 else "s"))
            if k == 0:
                legend(ax, th, loc="upper left", markerscale=6)
    fig.suptitle("Why the exponent is positive: a segment of initial "
                 "conditions, stretched and folded back onto the attractor. "
                 "Colour is position along the original segment.",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "chaos-mu1-fold")


# ------------------------------------------------------------------- main
def zoom_box(Y, half):
    """A square window of the given half width on an attractor's middle."""
    c = np.median(Y, axis=0)
    return (c[0] - half, c[0] + half, c[1] - half, c[1] + half)


def main(quick=False, fresh=False):
    t0 = time.time()
    print("Nearly harmonic mode, mu = %g, drive strength A = %g" % (MU, AMP))
    print("  Van der Pol free cycle: w_lc = %.5f, T = %.4f"
          % (W_LC, 2.0*np.pi/W_LC))
    print("  prototype levels %s edges %s" % (LEVELS, EDGES))
    print("  exact map against the integrator, relative: %.2e"
          % check_map_against_integration())

    scan = None
    if not quick:
        print("\nDrive sweep, ratio %.2f to %.2f, %d cells per system"
              % (BAND_RS[0], BAND_RS[-1], len(BAND_RS)))
        scan = band_scan(fresh=fresh)
        for kind in ("proto", "vdp"):
            print("  %-11s first pass:  %s" % (
                SYS_SHORT[kind],
                ", ".join("%.3f-%.3f" % b for b in bands(scan[kind]))
                or "none"))
        confirmed = confirm_band(scan, cache=CACHE.replace(
            "scan", "confirm"), fresh=fresh)
        n_bad = apply_confirmation(scan, confirmed)
        print("  %d of %d chaotic cells failed confirmation"
              % (n_bad, len(confirmed)))
        for kind in ("proto", "vdp"):
            print("  %-11s confirmed:   %s" % (
                SYS_SHORT[kind],
                ", ".join("%.3f-%.3f" % b for b in bands(scan[kind]))
                or "none"))
        for kind in ("proto", "vdp"):
            lo, hi = CASCADE[kind]
            seq, last = [], None
            for r, lab, _, _ in scan[kind]:
                if lo <= r <= hi and lab != last:
                    seq.append("%s from %.3f" % (lab, r))
                    last = lab
            print("  %-11s through its window %.3f-%.3f: %s"
                  % (SYS_SHORT[kind], lo, hi, ", ".join(seq)))
        print("  every chaotic cell, at five times the run and a wider "
              "twin:")
        for kind, r, ok, lams in confirmed:
            print("    %-11s r = %.3f  %-9s lam = %+.4f %+.4f %+.4f"
                  % (SYS_SHORT[kind], r, "confirmed" if ok else "REJECTED",
                     *lams))

    print("\nThe worked points: prototype r = %.3f, Van der Pol r = %.3f"
          % (R_PROTO, R_VDP))
    att = attractors()
    series = {}
    for kind in ("proto", "vdp"):
        r = R_PROTO if kind == "proto" else R_VDP
        t, x, v = trace(kind, r, N_TRACE, per=TRACE_PER, y0=att[kind][0][-1])
        series[kind] = (t, x, v, TRACE_PER)
    for kind in ("proto", "vdp"):
        r = R_PROTO if kind == "proto" else R_VDP
        Y, s, ld = att[kind]
        td = td_of(r)
        lam1 = section.lyapunov(flow(kind, r), td, [2.0, 0.0], N_SKIP, n=1500)
        lam2, dky, lsum = kaplan_yorke(lam1, ld, td)
        dim, ks, counts = box_dimension(Y)
        print("  %s, r = %.3f, Om = %.4f, drive period %.3f"
              % (SYS_SHORT[kind], r, om_of(r), td))
        print("    samples %d, x in [%.3f, %.3f], xdot in [%.3f, %.3f]"
              % (len(Y), Y[:, 0].min(), Y[:, 0].max(),
                 Y[:, 1].min(), Y[:, 1].max()))
        print("    lam1 = %+.4f, lam2 = %+.4f, sum = %+.4f per unit time"
              % (lam1, lam2, lsum))
        print("    area factor per drive period: %.4f"
              % float(np.exp(np.mean(ld))))
        print("    stretch per step: min %.3f, median %.3f, max %.3f"
              % (float(s.min()), float(np.median(s)), float(s.max())))
        print("    Kaplan-Yorke dimension %.3f, box counting %.3f over "
              "k = %s" % (dky, dim, ks))

    for kind in ("proto", "vdp"):
        _, x, _, _ = series[kind]
        ax = np.abs(x)
        print("  %-11s orbit peaks at %.3f; time with |x| below %.2f: "
              "%.1f%%, between the edges: %.1f%%, beyond %.2f: %.1f%%"
              % (SYS_SHORT[kind], float(ax.max()), EDGES[0],
                 100.0*float(np.mean(ax < EDGES[0])),
                 100.0*float(np.mean((ax >= EDGES[0]) & (ax <= EDGES[1]))),
                 EDGES[1], 100.0*float(np.mean(ax > EDGES[1]))))

    v_edge, jumps = edge_jump(Y=att["proto"][0])
    print("\n  across the inner edge a = %.4f at xdot = %.4f:" % (EDGES[0],
                                                                 v_edge))
    for e, dy, dj in jumps:
        print("    eps %8.1e   |dy| %.3e   |dJ|/|J| %.3f" % (e, dy, dj))

    # Both zoom windows are the same size, so a strand seen in one panel is
    # the same width as a strand seen in the other.
    zhalf = 0.09*max(float(np.max(np.ptp(att[k][0], axis=0))) for k in att)
    zooms = {k: zoom_box(att[k][0], zhalf) for k in att}
    print("\n  zoom windows: " + ", ".join(
        "%s [%.2f, %.2f] x [%.2f, %.2f]" % ((SYS_SHORT[k],) + zooms[k])
        for k in ("proto", "vdp")))

    folds = {}
    for kind in ("proto", "vdp"):
        r = R_PROTO if kind == "proto" else R_VDP
        Y = att[kind][0]
        c = np.median(Y, axis=0)
        w = 0.18*max(np.ptp(Y[:, 0]), np.ptp(Y[:, 1]))
        folds[kind] = fold_images(kind, r, (c[0] - w, c[1]), (c[0] + w, c[1]))
        L = [float(np.sum(np.hypot(*np.diff(P, axis=0).T))) for P in folds[kind]]
        print("  %-11s segment length by image: %s"
              % (SYS_SHORT[kind], ", ".join("%.3f" % l for l in L)))

    from figures import THEMES
    print("\nFigures")
    for name, th in THEMES.items():
        if scan is not None:
            fig_band(th, name, scan)
        fig_time(th, name, series)
        fig_phase(th, name, series, att)
        fig_section(th, name, att, zooms)
        fig_fold(th, name, folds, att)
    print("\ndone in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main(quick="quick" in sys.argv[1:], fresh="fresh" in sys.argv[1:])
