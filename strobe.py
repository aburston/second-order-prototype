"""The stroboscopic map of the three level prototype as a two dimensional section.

Run ``python3 strobe.py`` to print every number ``STROBOSCOPIC.md`` quotes
and to write its figures into ``figures/`` as ``strobe-*.png``, light and
dark. The whole run takes a few minutes on four cores; ``python3 strobe.py
quick`` skips the phase plane partition and its figure.

``MAPS.md`` sampled the prototypes on a geometric section, the line
``xdot = 0``, and reduced each cycle to a scalar recurrence in the
amplitude. This module changes the section. The forced prototype is
sampled **once per drive period**, so the section is the whole phase plane
and the map sends a point ``(x, xdot)`` to where it is one drive period
later. Only the three level prototype of ``THREELEVEL.md`` is treated, at
the drive it was fitted at, and every learning result is repeated on the
Van der Pol oscillator it was fitted to.

What the module establishes, in order:

1. **The difference equation.** Within one drive period the orbit visits a
   finite sequence of zones -- its *itinerary* -- and dwells a known time
   in each. The map is then a product of the ``5x5`` forced arc matrices of
   ``maps.py``, exactly, and on the plane it reads ``y_{n+1} = A(tau) y_n +
   b(tau)`` with ``A`` a product of zone transition matrices. Checked
   against the event stepping to ``1e-14``.

2. **The partition.** The itinerary is a function of the starting point,
   so the plane is tiled into cells of constant itinerary, and inside each
   cell the map is smooth. The tiling is measured on a grid, and the
   attractor is found to visit a small number of the cells.

3. **The edges.** Across a grazing edge the map is continuous with a
   continuous derivative; the singularity sits one order higher, as a
   square root in the derivative's variation. That is why getting a
   point's cell wrong costs almost nothing.

4. **Multiple prediction problems.** The map is fitted from stroboscopic
   snapshots alone, the way it would be from a measurement, as one global
   regression, as a gated set of per cell regressions, as local linear
   models in the manner of Farmer and Sidorowich, and as a regression of
   the dwells feeding the exact matrix product. One step and several
   steps ahead, recursively and directly.

5. **Van der Pol.** The same fits on Van der Pol's own snapshots, and the
   prototype's map compared with Van der Pol's point by point.
"""
import sys
import time
from collections import Counter, defaultdict

import numpy as np
from scipy.integrate import solve_ivp
from multiprocessing import Pool

import maps
import staircase
import vanderpol
import section

# -------------------------------------------------------------- the system
#: The three level prototype fitted to Van der Pol at ``mu = 5``, and the
#: drive it was fitted at: ``A = 5`` at ``Om = 2.470``, inside the chaotic
#: band both systems have at the lock 3 to lock 5 transition.
LEVELS, EDGES = staircase.THREE_FITTED
MU, AMP, OM = staircase.CMP_MU, staircase.CMP_AMP, 2.470
TD = 2.0*np.pi/OM

#: A drive frequency on the lock 3 plateau of both systems, for the poles.
OM_LOCK = 2.30

#: Feature scaling for every fit: the attractor is about three times taller
#: in ``xdot`` than it is wide in ``x``, so distances are measured with
#: ``xdot`` divided by three. Without this a nearest neighbour is nearly
#: always a neighbour in velocity alone.
SX = np.array([1.0, 3.0])

#: Snapshots used for every fit: a transient discarded, then a training
#: record and a held out test record from the same free running orbit.
N_SKIP, N_TRAIN, N_TEST = 400, 4000, 2000

#: Grid for the phase plane partition.
X_RANGE, V_RANGE, N_GRID = (-2.4, 2.4), (-8.0, 8.0), 161

MODEL = maps.staircase_model(LEVELS, EDGES)
VDP = vanderpol.field(MU, AMP, OM)


# ------------------------------------------------------- the exact map
def step(y, model=MODEL, amp=AMP, om=OM):
    """One drive period of the prototype, from drive phase zero.

    Returns ``(y_next, J, itinerary, dwells)``. The itinerary is a tuple of
    ``(zone, wall)`` pairs, one per arc, with ``wall`` the level the arc
    ends on and ``None`` for the last arc, which ends on the clock; the
    dwells are the arc durations, which sum to the drive period.
    """
    rec = []
    y2, _, J = maps.strobe_step(model, np.asarray(y, float), 0.0, amp, om,
                                record=rec)
    # A wall met and left again into the same zone is not a change of
    # zone, so consecutive arcs in one zone are merged: the itinerary is
    # the sequence of zones, and the dwells the time in each.
    merged = []
    for k, t, w in rec:
        if merged and merged[-1][0] == k:
            merged[-1] = (k, merged[-1][1] + t, w)
        else:
            merged.append((k, t, w))
    lab = tuple((k, None if w is None else round(w, 6)) for k, t, w in merged)
    return y2, J, lab, np.array([t for _, t, _ in merged])


def _step_worker(y):
    try:
        y2, _, lab, dw = step(y)
    except RuntimeError:
        return None
    return y2, lab, dw


def orbit(y0, n_skip, n, model=MODEL, amp=AMP, om=OM):
    """Stroboscopic samples of the prototype, with itineraries and dwells.

    Returns ``(Y, labels, dwells)`` with ``Y`` of shape ``(n + 1, 2)``;
    ``labels[i]`` and ``dwells[i]`` describe the step from ``Y[i]`` to
    ``Y[i + 1]``.
    """
    y = np.asarray(y0, float)
    for _ in range(n_skip):
        y, _, _, _ = step(y, model, amp, om)
    Y = np.empty((n + 1, 2))
    Y[0] = y
    labels, dwells = [], []
    for i in range(n):
        y, _, lab, dw = step(y, model, amp, om)
        Y[i + 1] = y
        labels.append(lab)
        dwells.append(dw)
    return Y, labels, dwells


def cycle_matrix(label, dwells, model=MODEL, amp=AMP, om=OM):
    """The step as one ``5x5`` matrix, from its itinerary and dwells.

    The product of the forced arc matrices of ``maps.py`` in the order the
    arcs are taken. On the stroboscopic section the drive state is
    ``(1, 0)`` at every sample, so the step is affine in ``y``::

        y_{n+1} = A y_n + b,   A = M[:2, :2],   b = M[:2, 2:] @ (1, 0, 1)

    with ``A`` the product of the zone transition matrices alone.
    """
    M = np.eye(5)
    for (k, _), t in zip(label, dwells):
        M = maps.forced_arc_matrix(model.zones[k][0], 0.0, t, amp, om) @ M
    return M


def apply_cycle(M, y):
    """Apply a cycle matrix to a state on the section."""
    return (M @ np.array([y[0], y[1], 1.0, 0.0, 1.0]))[:2]


def affine_parts(M):
    """``(A, b)`` of the affine step a cycle matrix encodes."""
    return M[:2, :2].copy(), M[:2, 2:] @ np.array([1.0, 0.0, 1.0])


FLOW = staircase.field(LEVELS, EDGES, AMP, OM)


def _integrate_period(y, rtol=1e-10, atol=1e-12, method="LSODA"):
    sol = solve_ivp(FLOW, (0.0, TD), np.asarray(y, float), method=method,
                    rtol=rtol, atol=atol)
    return sol.y[:, -1]


def check_integration(Y, workers=None):
    """Every step of the orbit against direct integration of the field.

    Returns the per step error. The integrator, not the map, sets the
    floor here: crossing a damping discontinuity costs an explicit or
    multistep integrator an error that tightening its tolerance does not
    remove, which ``MAPS.md`` measured at about ``1e-6`` on this model.
    """
    with Pool(workers) as p:
        Z = np.array(p.map(_integrate_period, list(Y[:-1]), chunksize=50))
    return np.linalg.norm(Z - Y[1:], axis=1)


def check_product(Y, labels, dwells, n=50):
    """Worst disagreement between the matrix product and the event stepping.

    Also confirms that ``A`` is the product of the zone ``Phi`` blocks and
    that the dwells of every step sum to the drive period.
    """
    worst = 0.0
    for i in range(n):
        M = cycle_matrix(labels[i], dwells[i])
        worst = max(worst, float(np.max(np.abs(apply_cycle(M, Y[i])
                                                - Y[i + 1]))))
        P = np.eye(2)
        for (k, _), t in zip(labels[i], dwells[i]):
            P = maps.phi(MODEL.zones[k][0], t) @ P
        assert np.allclose(P, M[:2, :2], atol=1e-12)
        assert abs(dwells[i].sum() - TD) < 1e-12
    return worst


# ---------------------------------------------------------- the partition
def partition(n=N_GRID, workers=None):
    """Itinerary of every point of a grid over the phase plane.

    Returns ``(xs, vs, labels)`` with ``labels`` a ``(len(vs), len(xs))``
    object array of itineraries, ``None`` where the step failed.
    """
    xs = np.linspace(*X_RANGE, n)
    vs = np.linspace(*V_RANGE, n)
    grid = [np.array([x, v]) for v in vs for x in xs]
    with Pool(workers) as p:
        res = p.map(_step_worker, grid, chunksize=50)
    labels = np.empty((len(vs), len(xs)), dtype=object)
    for i, r in enumerate(res):
        labels[i // len(xs), i % len(xs)] = None if r is None else r[1]
    return xs, vs, labels


def describe(label):
    """An itinerary as text: zones visited, with the wall crossed between."""
    out = []
    for k, w in label:
        out.append(str(k))
        if w is not None:
            out.append("|%+.2f|" % w)
    return " ".join(out)


# ------------------------------------------------------------- the edges
def find_boundary(Y, rng, max_gap=0.8, tries=400):
    """A point on an edge between two cells, with the unit normal across it.

    Draws pairs of attractor points with different itineraries and bisects
    the segment between them. Returns ``(y_b, n, label_a, label_b)`` with
    ``label_a`` the itinerary on the ``-n`` side; ``None`` if no clean edge
    was found.
    """
    for _ in range(tries):
        i, j = rng.integers(0, len(Y), 2)
        a, b = Y[i], Y[j]
        if np.linalg.norm(a - b) > max_gap:
            continue
        la = step(a)[2]
        lb = step(b)[2]
        if la == lb:
            continue
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5*(lo + hi)
            if step(a + mid*(b - a))[2] == la:
                lo = mid
            else:
                hi = mid
        yb = a + 0.5*(lo + hi)*(b - a)
        nrm = (b - a)/np.linalg.norm(b - a)
        if step(yb - 1e-7*nrm)[2] == la and step(yb + 1e-7*nrm)[2] == lb:
            return yb, nrm, la, lb
    return None


def edge_scaling(yb, nrm, ss=(1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6)):
    """How the map varies across an edge, against distance from it.

    Symmetric differences straddling the edge, so that whatever is smooth
    cancels and only the edge's own contribution is measured: the second
    difference of the map, ``|P(y_b + s n) + P(y_b - s n) - 2 P(y_b)|``,
    and the change in the Jacobian, ``|J(y_b + s n) - J(y_b - s n)|``.
    Returns the rows, the two fitted power law exponents, and the jump in
    the map itself across the edge at ``s = 1e-9``.

    A map that is ``C^1`` across the edge with a square root in its
    derivative gives exponents ``3/2`` and ``1/2``; one with a derivative
    jump gives ``1`` and ``0``; a smooth one gives ``2`` and ``1``.
    """
    P0, _, _, _ = step(yb)
    Pm, Jm, _, _ = step(yb - 1e-9*nrm)
    Pp, Jp, _, _ = step(yb + 1e-9*nrm)
    jump = float(np.linalg.norm(Pp - Pm))
    rows = []
    for s in ss:
        Pa, Ja, _, _ = step(yb + s*nrm)
        Pb, Jb, _, _ = step(yb - s*nrm)
        rows.append((s, float(np.linalg.norm(Pa + Pb - 2.0*P0)),
                     float(np.max(np.abs(Ja - Jb)))))
    r = np.array(rows)
    slopes = []
    for col in (1, 2):
        ok = r[:, col] > 10.0*max(jump, 1e-14)
        slopes.append(float(np.polyfit(np.log(r[ok, 0]), np.log(r[ok, col]),
                                       1)[0]) if ok.sum() > 2 else np.nan)
    return rows, slopes, jump


# ------------------------------------------------------------- the locks
def lock_order(om, model=MODEL, n_skip=600, n_keep=64, qmax=11, tol=1e-8):
    """Lock order of the prototype at ``om`` from the exact map, 0 if none."""
    y = np.array([2.0, 0.0])
    for _ in range(n_skip):
        y, _, _ = maps.strobe_step(model, y, 0.0, AMP, om)
    pts = [y]
    for _ in range(n_keep):
        y, _, _ = maps.strobe_step(model, y, 0.0, AMP, om)
        pts.append(y)
    pts = np.array(pts)
    scale = float(np.max(np.abs(pts)))
    for q in range(1, qmax + 1):
        if np.max(np.abs(pts[q:] - pts[:-q]))/scale < tol:
            return q, pts
    return 0, pts


def locked_poles(om, q, model=MODEL, n_newton=30):
    """The period ``q`` orbit at ``om`` and the eigenvalues of its Jacobian.

    Newton on ``P^q(y) - y`` with the exact Jacobian product, from the
    iterated orbit. On a two dimensional section a locked orbit has two
    multipliers, and these are its poles in the z domain.
    """
    _, pts = lock_order(om, model)
    y = pts[-1].copy()
    for _ in range(n_newton):
        z, J = y.copy(), np.eye(2)
        for _ in range(q):
            z, _, j = maps.strobe_step(model, z, 0.0, AMP, om)
            J = j @ J
        dy = np.linalg.solve(J - np.eye(2), -(z - y))
        y = y + dy
        if np.linalg.norm(dy) < 1e-14:
            break
    return y, np.linalg.eigvals(J), float(np.linalg.norm(z - y))


# -------------------------------------------------------------- learning
def features(X, d):
    """Monomials of the scaled state up to degree ``d``."""
    Z = np.atleast_2d(X)/SX
    cols = [np.ones(len(Z))]
    for p in range(1, d + 1):
        for i in range(p + 1):
            cols.append(Z[:, 0]**(p - i)*Z[:, 1]**i)
    return np.array(cols).T


def lstsq(F, T, ridge=1e-10):
    """Ridge regularised least squares, ``F w = T``."""
    n = F.shape[1]
    return np.linalg.solve(F.T @ F + ridge*np.eye(n), F.T @ T)


def rel_errors(pred, true, scale):
    """Median, 90th percentile and maximum error, relative to ``scale``."""
    e = np.linalg.norm(pred - true, axis=1)/scale
    e = e[np.isfinite(e)]
    return float(np.median(e)), float(np.percentile(e, 90)), float(e.max())


def nearest_label(x, Xref, Lref):
    """The itinerary of the nearest training snapshot: the learned gate."""
    d2 = np.sum(((Xref - x)/SX)**2, axis=1)
    return Lref[int(np.argmin(d2))]


def local_linear(x, Xref, Tref, k, Lref=None, label=None):
    """Weighted local linear prediction from the ``k`` nearest snapshots.

    Farmer and Sidorowich's method: a separate affine model for every
    query, fitted to its neighbours with tricube weights. With ``Lref``
    and ``label`` given the neighbours are restricted to one cell. Returns
    ``None`` when fewer than four neighbours are available.
    """
    Z = Xref/SX
    d2 = np.sum((Z - x/SX)**2, axis=1)
    if Lref is not None:
        d2 = np.where(np.array([l == label for l in Lref]), d2, np.inf)
    idx = np.argsort(d2)[:k]
    idx = idx[np.isfinite(d2[idx])]
    if len(idx) < 4:
        return None
    F = np.hstack([np.ones((len(idx), 1)), Z[idx] - x/SX])
    h = np.sqrt(d2[idx].max()) + 1e-12
    w = (1.0 - (np.sqrt(d2[idx])/h)**3)**3 + 1e-6
    Fw = F*w[:, None]
    coef = np.linalg.solve(Fw.T @ F + 1e-10*np.eye(3), Fw.T @ Tref[idx])
    return coef[0]


def fit_cells(X, T, L, d):
    """One polynomial of degree ``d`` per cell, lowered where data is short."""
    cells = defaultdict(list)
    for i, l in enumerate(L):
        cells[l].append(i)
    W = {}
    for l, idx in cells.items():
        idx = np.array(idx)
        dd = d
        while dd > 0 and len(idx) < 2*features(X[:1], dd).shape[1]:
            dd -= 1
        W[l] = (dd, lstsq(features(X[idx], dd), T[idx], ridge=1e-8))
    return W


def predict_cells(Xq, gates, W):
    out = np.full((len(Xq), next(iter(W.values()))[1].shape[1]), np.nan)
    for i, (x, g) in enumerate(zip(Xq, gates)):
        if g in W:
            dd, w = W[g]
            out[i] = features(x, dd) @ w
    return out


def fit_dwells(X, D, L, d):
    """Per cell polynomials for all but the last dwell; the last is fixed.

    The dwells of a step sum to the drive period, so a cell with ``m`` arcs
    has ``m - 1`` free dwells and the last is what remains.
    """
    cells = defaultdict(list)
    for i, l in enumerate(L):
        cells[l].append(i)
    W = {}
    for l, idx in cells.items():
        idx = np.array(idx)
        m = len(l)
        if m == 1:
            W[l] = (0, None)
            continue
        T = np.array([D[i][:m - 1] for i in idx])
        dd = d
        while dd > 0 and len(idx) < 2*features(X[:1], dd).shape[1]:
            dd -= 1
        W[l] = (dd, lstsq(features(X[idx], dd), T, ridge=1e-8))
    return W


def predict_dwells(Xq, gates, W):
    """Predicted dwells through the exact matrix product."""
    out = np.full((len(Xq), 2), np.nan)
    for i, (x, g) in enumerate(zip(Xq, gates)):
        if g not in W:
            continue
        dd, w = W[g]
        if w is None:
            dw = np.array([TD])
        else:
            dw = features(x, dd) @ w
            dw = np.append(dw, TD - dw.sum())
        out[i] = apply_cycle(cycle_matrix(g, dw), x)
    return out


def predict_local(Xq, Xref, Tref, k, Lref=None):
    """Local linear prediction of every query, gated when ``Lref`` is given."""
    out = np.full((len(Xq), Tref.shape[1]), np.nan)
    for i, x in enumerate(Xq):
        g = nearest_label(x, Xref, Lref) if Lref is not None else None
        r = local_linear(x, Xref, Tref, k, Lref, g)
        if r is not None:
            out[i] = r
    return out


def horizons(Y, n_train, k=8, hmax=8, Lref=None):
    """Recursive against direct prediction ``h`` periods ahead.

    Recursive iterates the one step local linear model ``h`` times; direct
    fits a fresh local linear model from ``y_n`` to ``y_{n+h}``, one
    prediction problem per horizon. Returns rows of
    ``(h, recursive median, recursive p90, direct median, direct p90)``.
    """
    scale = np.sqrt(np.mean(Y**2))
    N = len(Y) - 1
    Xtr, Ytr = Y[:n_train], Y[1:n_train + 1]
    Ltr = None if Lref is None else Lref[:n_train]
    Xte = Y[n_train:N]
    rows = []
    for h in range(1, hmax + 1):
        Xc = Xte[:len(Xte) - h].copy()
        for _ in range(h):
            Xc = predict_local(Xc, Xtr, Ytr, k, Ltr)
        truth = Y[n_train + h:N + 1][:len(Xc)]
        rec = rel_errors(Xc, truth, scale)
        Xh, Yh = Y[:n_train - h + 1], Y[h:n_train + 1]
        Lh = None if Lref is None else Lref[:n_train - h + 1]
        direct = rel_errors(predict_local(Xte[:len(Xte) - h], Xh, Yh, k, Lh),
                            truth, scale)
        rows.append((h, rec[0], rec[1], direct[0], direct[1]))
    return rows


# ---------------------------------------------------------- Van der Pol
def vdp_step(y, rtol=1e-10, atol=1e-12):
    """One drive period of Van der Pol from drive phase zero."""
    sol = solve_ivp(VDP, (0.0, TD), np.asarray(y, float), method="LSODA",
                    rtol=rtol, atol=atol)
    return sol.y[:, -1]


def vdp_orbit(y0, n_skip, n, rtol=1e-10, atol=1e-12):
    """Stroboscopic samples of Van der Pol, transient discarded."""
    t = TD*np.arange(0, n_skip + n + 1)
    sol = solve_ivp(VDP, (0.0, t[-1]), np.asarray(y0, float), t_eval=t,
                    method="LSODA", rtol=rtol, atol=atol)
    return sol.y.T[n_skip:]


def _proto_map(y):
    return maps.strobe_step(MODEL, y, 0.0, AMP, OM)[0]


def maps_compared(pts, workers=None):
    """The two one period maps at the same points, relative to the scale."""
    with Pool(workers) as p:
        Pp = np.array(p.map(_proto_map, pts, chunksize=20))
        Pv = np.array(p.map(vdp_step, pts, chunksize=20))
    return Pp, Pv


def overlap(A, B, delta):
    """Fraction of the points of ``A`` within ``delta`` of a point of ``B``."""
    hits = 0
    for a in A:
        if np.min(np.sum(((B - a)/SX)**2, axis=1)) < delta**2:
            hits += 1
    return hits/len(A)


# -------------------------------------------------------------- figures
def fig_partition(th, name, xs, vs, labels, Y):
    """The phase plane tiled by itinerary, with the attractor on top.

    Cells are filled by how many arcs their itinerary has, a magnitude, on
    one sequential ramp; identity is carried by the drawn cell edges, not
    by colour. The zone edges of the prototype are chrome, since they are
    the coordinate system the tiling is built on.
    """
    import matplotlib
    from matplotlib.collections import LineCollection
    from figures import newfig, style, legend, save

    ids = {}
    idg = np.full(labels.shape, -1)
    arcs = np.full(labels.shape, np.nan)
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            l = labels[i, j]
            if l is None:
                continue
            idg[i, j] = ids.setdefault(l, len(ids))
            arcs[i, j] = len(l)
    ramp = matplotlib.colors.LinearSegmentedColormap.from_list(
        "arcs", [th["grid"], th["series"][0]])
    fig, ax = newfig(th, figsize=(7.8, 6.4))
    dx, dv = xs[1] - xs[0], vs[1] - vs[0]
    n_max = int(np.nanmax(arcs))
    mesh = ax.pcolormesh(np.append(xs - dx/2, xs[-1] + dx/2),
                         np.append(vs - dv/2, vs[-1] + dv/2), arcs,
                         cmap=ramp, vmin=0.5, vmax=n_max + 0.5,
                         shading="flat", zorder=1, rasterized=True)
    segs = []
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1] - 1):
            if idg[i, j] != idg[i, j + 1]:
                xm = 0.5*(xs[j] + xs[j + 1])
                segs.append([(xm, vs[i] - dv/2), (xm, vs[i] + dv/2)])
    for i in range(labels.shape[0] - 1):
        for j in range(labels.shape[1]):
            if idg[i, j] != idg[i + 1, j]:
                vm = 0.5*(vs[i] + vs[i + 1])
                segs.append([(xs[j] - dx/2, vm), (xs[j] + dx/2, vm)])
    ax.add_collection(LineCollection(segs, colors=th["ink2"], linewidths=0.5,
                                     zorder=2))
    for e in EDGES:
        for xv in (e, -e):
            ax.axvline(xv, color=th["ink"], linewidth=1.0,
                       linestyle=(0, (5, 3)), zorder=3)
    ax.plot(Y[:, 0], Y[:, 1], ".", color=th["series"][1], markersize=1.4,
            zorder=4, label="stroboscopic samples of the attractor, %d"
            % len(Y))
    ax.plot([], [], color=th["ink2"], linewidth=0.8, label="cell edges")
    ax.plot([], [], color=th["ink"], linewidth=1.0, linestyle=(0, (5, 3)),
            label="zone edges $\\pm a$, $\\pm b$")
    cb = fig.colorbar(mesh, ax=ax, ticks=range(1, n_max + 1), pad=0.02,
                      fraction=0.04)
    cb.set_label("arcs in the itinerary", color=th["ink2"], fontsize=9)
    cb.ax.tick_params(colors=th["ink2"], labelsize=8)
    cb.outline.set_visible(False)
    style(ax, th, "$x$", "$\\dot{x}$",
          "Itinerary cells of the stroboscopic map, %d on this window, "
          "filled by arc count" % len(ids))
    ax.set_xlim(*X_RANGE)
    ax.set_ylim(*V_RANGE)
    legend(ax, th, loc="upper left", markerscale=4)
    fig.tight_layout()
    save(fig, name, "strobe-partition")


def fig_horizon(th, name, rows_p, rows_v, lam_p, lam_v):
    """Prediction error against horizon, recursive and direct, both systems."""
    from figures import newfig, style, legend, save
    fig, ax = newfig(th, figsize=(7.2, 4.6))
    for rows, col, sysname in ((rows_p, th["series"][0], "prototype"),
                               (rows_v, th["series"][1], "Van der Pol")):
        h = [r[0] for r in rows]
        ax.plot(h, [r[2] for r in rows], "-o", color=col, linewidth=2,
                markersize=5, label="%s, recursive" % sysname)
        ax.plot(h, [r[4] for r in rows], "--s", color=col, linewidth=2,
                markersize=5, markerfacecolor=th["surface"],
                label="%s, direct" % sysname)
    h = np.arange(1, rows_p[-1][0] + 1)
    for rows, lam in ((rows_p, lam_p), (rows_v, lam_v)):
        ax.plot(h, rows[0][2]*np.exp(lam*TD*(h - 1)), ":", color=th["ink2"],
                linewidth=1.2)
    ax.text(h[-1], rows_v[0][2]*np.exp(lam_v*TD*(h[-1] - 1))*0.8,
            "one step error grown at $e^{\\lambda h T}$", color=th["ink2"],
            fontsize=8, ha="right", va="top")
    ax.set_yscale("log")
    style(ax, th, "horizon $h$, drive periods ahead",
          "90th percentile error, relative to the attractor's rms radius",
          "Predicting $h$ periods ahead from snapshots, local linear models")
    legend(ax, th, loc="upper left")
    fig.tight_layout()
    save(fig, name, "strobe-horizon")


# ----------------------------------------------------------------- main
def _table(rows, header, fmt):
    print("  " + header)
    for r in rows:
        print("  " + fmt % r)


def main(quick=False):
    from figures import THEMES
    t_start = time.time()
    rng = np.random.default_rng(1)

    print("three level prototype, zeta = %s, edges = %s" % (
        tuple(round(z, 4) for z in LEVELS), tuple(round(e, 4) for e in EDGES)))
    print("drive A = %g, Om = %.3f, period T = %.6f" % (AMP, OM, TD))

    # ---- the orbit, and the difference equation
    N = N_TRAIN + N_TEST
    Y, labels, dwells = orbit([2.0, 0.0], N_SKIP, N)
    scale = float(np.sqrt(np.mean(Y**2)))
    print("\nprototype: %d stroboscopic samples after %d discarded, rms "
          "radius %.4f, x in [%.3f, %.3f], xdot in [%.3f, %.3f]"
          % (N, N_SKIP, scale, Y[:, 0].min(), Y[:, 0].max(),
             Y[:, 1].min(), Y[:, 1].max()))
    print("1. matrix product against event stepping, worst over 50 steps: "
          "%.1e" % check_product(Y, labels, dwells))
    t0 = time.time()
    e = check_integration(Y)
    worst = int(np.argmax(e))
    print("   every step against integration (LSODA, rtol 1e-10): median "
          "%.1e, max %.1e, %d steps above 1e-8 and %d above 1e-6 (%.0fs)"
          % (np.median(e), e.max(), (e > 1e-8).sum(), (e > 1e-6).sum(),
             time.time() - t0))
    alt = [_integrate_period(Y[worst], 1e-12, 1e-14),
           _integrate_period(Y[worst], 1e-10, 1e-12, "DOP853")]
    print("   the worst step, itinerary [%s]: map against LSODA at 1e-12 "
          "%.1e, against DOP853 %.1e; the two integrators against each "
          "other %.1e" % (describe(labels[worst]),
                          np.linalg.norm(alt[0] - Y[worst + 1]),
                          np.linalg.norm(alt[1] - Y[worst + 1]),
                          np.linalg.norm(alt[0] - alt[1])))
    print("   crossings the step's scan missed and recovered over the orbit: "
          "%d; excursions below resolution dropped: %d"
          % (maps.STROBE_REPAIRS["recovered"], maps.STROBE_REPAIRS["dropped"]))
    A_, b_ = affine_parts(cycle_matrix(labels[0], dwells[0]))
    print("   first step: itinerary [%s], dwells %s" % (
        describe(labels[0]), np.round(dwells[0], 4)))
    print("   A = %s\n   b = %s" % (np.round(A_, 6).tolist(),
                                    np.round(b_, 6).tolist()))

    # ---- itineraries on the attractor
    count = Counter(labels)
    print("\n2. distinct itineraries on the attractor over %d steps: %d"
          % (N, len(count)))
    print("   share  arcs  itinerary")
    cum = 0.0
    for lab, c in count.most_common():
        cum += c/N
        print("   %5.1f%%  %d     %s" % (100.0*c/N, len(lab), describe(lab)))
    top8 = sum(c for _, c in count.most_common(8))/N
    print("   the eight most visited cover %.1f%%" % (100.0*top8))

    # ---- the partition
    if not quick:
        t0 = time.time()
        xs, vs, grid_labels = partition()
        flat = [l for l in grid_labels.ravel() if l is not None]
        gc = Counter(flat)
        print("\n   grid %dx%d over x in %s, xdot in %s: %d distinct "
              "itineraries, %d failed points (%.0fs)"
              % (N_GRID, N_GRID, X_RANGE, V_RANGE, len(gc),
                 sum(l is None for l in grid_labels.ravel()),
                 time.time() - t0))
        print("   arcs per itinerary on the grid: %s"
              % sorted(Counter(len(l) for l in flat).items()))
        on_attractor = set(count)
        print("   grid cells the attractor visits: %d of %d; attractor "
              "itineraries absent from the grid: %d"
              % (len(on_attractor & set(gc)), len(gc),
                 len(on_attractor - set(gc))))
        for th_name, th in THEMES.items():
            fig_partition(th, th_name, xs, vs, grid_labels, Y[:3000])

    # ---- the edges
    print("\n3. across cell edges: symmetric second difference of the map and "
          "difference of the Jacobian, against distance s either side")
    seen = set()
    n_edges = 0
    while n_edges < 3:
        found = find_boundary(Y[:1500], rng)
        if found is None:
            break
        yb, nrm, la, lb = found
        key = (la, lb)
        if key in seen:
            continue
        seen.add(key)
        n_edges += 1
        rows, slopes, jump = edge_scaling(yb, nrm)
        _, _, _, dwa = step(yb - 1e-7*nrm)
        _, _, _, dwb = step(yb + 1e-7*nrm)
        print("   edge %d at (%.4f, %.4f) between [%s] and [%s]"
              % (n_edges, yb[0], yb[1], describe(la), describe(lb)))
        print("     dwells either side: %s | %s"
              % (np.round(dwa, 5), np.round(dwb, 5)))
        _table(rows, "      s      second difference   Jacobian difference",
               "  %8.1e   %14.2e   %14.2e")
        print("     fitted exponents: second difference %.2f, Jacobian "
              "difference %.2f; jump in the map across the edge %.1e"
              % (slopes[0], slopes[1], jump))

    # ---- locks and poles
    print("\n4. lock order from the exact map, against Van der Pol")
    for om in (2.20, 2.30, 2.40, 2.45, 2.47, 2.50, 2.60, 3.00):
        q, _ = lock_order(om)
        pts = section.strobe(VDP if om == OM else vanderpol.field(MU, AMP, om),
                             2.0*np.pi/om, [2.0, 0.0], 400)
        qv = section.lock_order(pts) or 0
        print("   Om = %.2f  prototype lock %d   Van der Pol lock %d"
              % (om, q, qv))
    q, _ = lock_order(OM_LOCK)
    y3, ev, res = locked_poles(OM_LOCK, q)
    print("   period %d orbit at Om = %.2f through (%.6f, %.6f), residual "
          "%.1e" % (q, OM_LOCK, y3[0], y3[1], res))
    print("   its poles: %s, moduli %s"
          % (np.round(ev, 6).tolist(), np.round(np.abs(ev), 6).tolist()))
    tdl = 2.0*np.pi/OM_LOCK
    fl = vanderpol.field(MU, AMP, OM_LOCK)
    mv, qv = section.multipliers(fl, tdl, [2.0, 0.0], 400)
    print("   Van der Pol's period %s orbit, by finite differences "
          "(section.multipliers): %s, moduli %s"
          % (qv, None if mv is None else np.round(mv, 6).tolist(),
             None if mv is None else np.round(np.abs(mv), 6).tolist()))

    # ---- learning from snapshots
    Xtr, Ytr, Ltr, Dtr = Y[:N_TRAIN], Y[1:N_TRAIN + 1], labels[:N_TRAIN], \
        dwells[:N_TRAIN]
    Xte, Yte, Lte = Y[N_TRAIN:N], Y[N_TRAIN + 1:N + 1], labels[N_TRAIN:]
    print("\n5. fitting the map from %d training snapshots, tested on the "
          "next %d; errors relative to the rms radius as median / 90th "
          "percentile / max" % (N_TRAIN, N_TEST))
    print("   itineraries in training %d, in test %d, unseen in training %d"
          % (len(set(Ltr)), len(set(Lte)), len(set(Lte) - set(Ltr))))
    print("   one global polynomial")
    for d in (2, 4, 6, 8):
        W = lstsq(features(Xtr, d), Ytr)
        print("     degree %d: %.1e / %.1e / %.1e"
              % ((d,) + rel_errors(features(Xte, d) @ W, Yte, scale)))
    gates = [nearest_label(x, Xtr, Ltr) for x in Xte]
    wrong = np.array([g != l for g, l in zip(gates, Lte)])
    print("   learned gate, nearest snapshot's itinerary: %.1f%% right "
          "(%d wrong of %d)" % (100.0*(1 - wrong.mean()), wrong.sum(),
                                len(Lte)))
    print("   per cell polynomial, learned gate | true cell")
    for d in (1, 2, 3):
        W = fit_cells(Xtr, Ytr, Ltr, d)
        p1 = rel_errors(predict_cells(Xte, gates, W), Yte, scale)
        p2 = rel_errors(predict_cells(Xte, Lte, W), Yte, scale)
        print("     degree %d: %.1e / %.1e / %.1e | %.1e / %.1e / %.1e"
              % ((d,) + p1 + p2))
    print("   local linear, k nearest snapshots, ungated | gated")
    for k in (8, 16, 32):
        p1 = rel_errors(predict_local(Xte, Xtr, Ytr, k), Yte, scale)
        p2 = rel_errors(predict_local(Xte, Xtr, Ytr, k, Ltr), Yte, scale)
        print("     k = %2d: %.1e / %.1e / %.1e | %.1e / %.1e / %.1e"
              % ((k,) + p1 + p2))
    print("   dwells regressed per cell, then the exact product, true cell")
    for d in (1, 2, 3):
        W = fit_dwells(Xtr, Dtr, Ltr, d)
        print("     polynomial degree %d: %.1e / %.1e / %.1e"
              % ((d,) + rel_errors(predict_dwells(Xte, Lte, W), Yte, scale)))
    pd_ = np.full((len(Xte), 2), np.nan)
    for i, (x, g) in enumerate(zip(Xte, Lte)):
        m = len(g)
        if m == 1:
            pd_[i] = apply_cycle(cycle_matrix(g, np.array([TD])), x)
            continue
        D = np.array([d[:m - 1] if len(d) == m else np.full(m - 1, np.nan)
                      for d in Dtr])
        r = local_linear(x, Xtr, D, 8, Ltr, g)
        if r is not None:
            pd_[i] = apply_cycle(cycle_matrix(g, np.append(r, TD - r.sum())),
                                 x)
    print("     local linear, k = 8: %.1e / %.1e / %.1e"
          % rel_errors(pd_, Yte, scale))
    pl = predict_local(Xte, Xtr, Ytr, 8, Ltr)
    e = np.linalg.norm(pl - Yte, axis=1)/scale
    print("   local linear k = 8, gated: error where the gate is right "
          "%.1e / %.1e / %.1e, where it is wrong %.1e / %.1e / %.1e"
          % (rel_errors(pl[~wrong], Yte[~wrong], scale)
             + rel_errors(pl[wrong], Yte[wrong], scale)))
    pe = np.array([apply_cycle(cycle_matrix(l, d), x)
                   for x, l, d in zip(Xte[:200], Lte[:200], dwells[N_TRAIN:
                                                                    N_TRAIN + 200])])
    print("   for reference, the true dwells through the product: "
          "%.1e / %.1e / %.1e" % rel_errors(pe, Yte[:200], scale))

    # ---- Van der Pol
    t0 = time.time()
    Yv = vdp_orbit([2.0, 0.0], N_SKIP, N)
    scale_v = float(np.sqrt(np.mean(Yv**2)))
    print("\n6. Van der Pol, mu = %g, same drive: %d samples, rms radius "
          "%.4f, x in [%.3f, %.3f], xdot in [%.3f, %.3f] (%.0fs)"
          % (MU, N, scale_v, Yv[:, 0].min(), Yv[:, 0].max(),
             Yv[:, 1].min(), Yv[:, 1].max(), time.time() - t0))
    floor = max(np.linalg.norm(vdp_step(y) - vdp_step(y, 1e-12, 1e-14))
                for y in Yv[:20])
    print("   integration floor of its one period map, rtol 1e-10 against "
          "1e-12: %.1e absolute, %.1e relative" % (floor, floor/scale_v))
    Xv, Yv1 = Yv[:N_TRAIN], Yv[1:N_TRAIN + 1]
    Xvt, Yvt = Yv[N_TRAIN:N], Yv[N_TRAIN + 1:N + 1]
    W = lstsq(features(Xv, 4), Yv1)
    print("   one global polynomial, degree 4: %.1e / %.1e / %.1e"
          % rel_errors(features(Xvt, 4) @ W, Yvt, scale_v))
    for k in (8, 16):
        print("   local linear, k = %d: %.1e / %.1e / %.1e"
              % ((k,) + rel_errors(predict_local(Xvt, Xv, Yv1, k), Yvt,
                                   scale_v)))

    # ---- horizons, both systems
    t0 = time.time()
    lam_p = maps.forced_lyapunov(MODEL, [2.0, 0.0], AMP, OM, n_skip=300, n=600)
    lam_v = section.lyapunov(VDP, TD, [2.0, 0.0], 300, n=600)
    print("\n7. Lyapunov exponents per unit time: prototype (exact Jacobian) "
          "%.4f, Van der Pol (twin trajectory) %.4f; per drive period %.4f "
          "and %.4f (%.0fs)" % (lam_p, lam_v, lam_p*TD, lam_v*TD,
                                time.time() - t0))
    rows_p = horizons(Y, N_TRAIN, Lref=labels)
    rows_v = horizons(Yv, N_TRAIN)
    print("   h periods ahead, local linear k = 8, median / p90:")
    print("   h   prototype recursive     prototype direct       "
          "Van der Pol recursive   Van der Pol direct    e^(lam h T)")
    for rp, rv in zip(rows_p, rows_v):
        print("   %d   %.1e / %.1e     %.1e / %.1e     %.1e / %.1e     "
              "%.1e / %.1e    %.1f / %.1f"
              % (rp[0], rp[1], rp[2], rp[3], rp[4], rv[1], rv[2], rv[3],
                 rv[4], np.exp(lam_p*TD*rp[0]), np.exp(lam_v*TD*rp[0])))
    for th_name, th in THEMES.items():
        fig_horizon(th, th_name, rows_p, rows_v, lam_p, lam_v)

    # ---- the two maps, point by point
    t0 = time.time()
    Pp, Pv = maps_compared(Yv[:1500])
    d = np.linalg.norm(Pp - Pv, axis=1)/scale_v
    print("\n8. the prototype's map against Van der Pol's at 1500 points of "
          "Van der Pol's attractor, relative to its rms radius: median %.3f, "
          "p90 %.3f, max %.3f; within 1%% at %.1f%% of points, within 10%% "
          "at %.1f%% (%.0fs)"
          % (np.median(d), np.percentile(d, 90), d.max(),
             100.0*np.mean(d < 0.01), 100.0*np.mean(d < 0.1),
             time.time() - t0))
    Pp2, Pv2 = maps_compared(Y[:1500])
    d2 = np.linalg.norm(Pp2 - Pv2, axis=1)/scale_v
    print("   at 1500 points of the prototype's attractor: median %.3f, p90 "
          "%.3f, max %.3f" % (np.median(d2), np.percentile(d2, 90), d2.max()))
    for delta in (0.02, 0.05, 0.10):
        print("   attractors as sets, delta = %.2f: Van der Pol points near "
              "the prototype's %.1f%%, prototype points near Van der Pol's "
              "%.1f%%" % (delta, 100.0*overlap(Yv[:1500], Y[:3000], delta),
                          100.0*overlap(Y[:1500], Yv[:3000], delta)))
    print("\ndone in %.0fs" % (time.time() - t_start))


if __name__ == "__main__":
    main(quick="quick" in sys.argv[1:])
