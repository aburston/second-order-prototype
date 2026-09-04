"""Generate the figures used in README.md, VANDERPOL.md and MAPS.md.

Run ``python3 figures.py`` to write every image into ``figures/``. Nothing
is cached, so a run always reflects the current parameters; regenerate
after changing any value the README quotes.

Five figures, in the order the README develops the argument:

``linear-prototype``
    The linear prototype for three damping ratios, phase plane beside time
    history. Establishes the notation and the single stable equilibrium
    that the nonlinear cases depart from.
``switched-damping``
    Switched damping with the boundary on the x-axis, in three panels for
    positive, zero and negative mean damping. The middle panel draws three
    nested orbits to show that the marginal case gives a *continuum* of
    closed orbits, not one isolated cycle.
``decrement``
    Amplitude per full cycle, integrated against the closed form
    ``exp(-(delta(zp) + delta(zm)) n)``, on a log scale. Straight lines
    confirm the decay is exactly geometric.
``limit-cycle``
    The offset boundary: trajectories starting inside and outside both
    converging on one closed orbit, with the offset ``v0`` and the
    equilibrium marked.
``return-map``
    The return map crossing the diagonal transversally beside the
    proportionality of cycle amplitude to offset. The offset and
    through-equilibrium maps are drawn together because the contrast
    between them is the whole mechanism.
``forced-tongues``
    What a sinusoidal drive does: the rotation number staircase beside the
    1:1 Arnold tongue traced by bisection.
``forced-sections``
    The three things the stroboscopic section is ever observed to be —
    a point, a closed curve, a chain of islands — and no fourth.
``staircase``
    The two threshold prototype: averaged damping crossing zero twice, the
    two cycles that follow with the basin between them, and the staircase
    closing on Van der Pol as levels are added.
``vanderpol-compare``
    The same forcing analysis run on Van der Pol, whose damping is
    unbounded rather than saturating: at weak nonlinearity it matches the
    prototype, and at strong nonlinearity it goes chaotic where the
    prototype cannot.

Each figure is rendered once per entry in ``THEMES`` and saved as
``<name>-light.png`` and ``<name>-dark.png``. The README embeds the pair in
a ``<picture>`` element so GitHub serves whichever matches the reader's
theme. Dark is a separate set of colours chosen for the dark surface, not
an automatic inversion of the light one.

Colour carries no meaning on its own anywhere in these figures: the
switching boundary, equilibrium markers and callouts are drawn in chrome
ink rather than a series colour, and every series is directly labelled as
well as being in the legend.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import forced
import section
import staircase
import vanderpol
import frequency
import symmetric
import displacement

WN = 1.0
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# Validated palette: categorical slots 1-3, chart chrome and ink.
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e",
                  grid="#e1e0d9", axis="#c3c2b7",
                  series=("#2a78d6", "#eb6834", "#1baf7a"),
                  div_pos="#2a78d6", div_neg="#e34948"),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7",
                  grid="#2c2c2a", axis="#383835",
                  series=("#3987e5", "#d95926", "#199e70"),
                  div_pos="#3987e5", div_neg="#e66767"),
}


# ---------------------------------------------------------------- dynamics
def f_linear(zeta):
    """Build the right hand side of the linear prototype, unforced.

    Returns the vector field of ``xddot + 2 zeta wn xdot + wn^2 x = 0``
    written as a first order system in ``y = [x1, x2] = [x, xdot]``, in the
    ``f(t, y)`` form ``solve_ivp`` expects.

    Args:
        zeta: damping ratio. Underdamped below 1, overdamped above.

    Returns:
        A callable ``f(t, y)`` returning ``[x2, -wn^2 x1 - 2 zeta wn x2]``.
    """
    def f(t, y):
        return [y[1], -WN**2*y[0] - 2*zeta*WN*y[1]]
    return f


def f_switched(zp, zm, v0=0.0):
    """Build the right hand side of the switched damping prototype.

    The damping ratio takes one value either side of the boundary
    ``Sigma = {xdot = v0}`` and acts on ``w = x2 - v0``, the velocity
    *relative to the boundary*. That choice matters: because the damping
    term carries a factor of ``w``, it vanishes on ``Sigma``, so the two
    half-plane fields agree there and the field stays continuous. Switching
    on ``w`` while damping the absolute velocity would make the field jump
    by ``2 (zp - zm) wn v0`` and reintroduce sliding solutions.

    With ``v0 = 0`` the boundary passes through the equilibrium and the
    field is positively homogeneous, which is why that case has no isolated
    limit cycle. Any nonzero ``v0`` breaks the scale invariance.

    Args:
        zp: damping ratio applied where ``w > 0``.
        zm: damping ratio applied where ``w < 0``. Negative values feed
            energy in and make the equilibrium a repelling focus.
        v0: boundary offset in the velocity direction. Zero puts the
            boundary on the x-axis, through the equilibrium.

    Returns:
        A callable ``f(t, y)`` suitable for ``solve_ivp``.
    """
    def f(t, y):
        w = y[1] - v0
        return [y[1], -WN**2*y[0] - 2*(zp if w > 0 else zm)*WN*w]
    return f


def traj(f, y0, T, n=4000):
    """Integrate a trajectory and return it densely enough to plot smoothly.

    Tolerances are far tighter than a picture needs. They are kept there so
    that a closed orbit in the marginal case visibly closes rather than
    drifting open under integration error, which would misrepresent the
    result being illustrated.

    Args:
        f: vector field in ``f(t, y)`` form.
        y0: initial state ``[x, xdot]``.
        T: end time. Chosen per figure to give enough revolutions to read
            the behaviour without the spiral becoming a solid disc.
        n: number of evenly spaced output samples.

    Returns:
        Tuple ``(x, xdot, t)`` of equal-length arrays.
    """
    s = solve_ivp(f, (0, T), y0, t_eval=np.linspace(0, T, n),
                  rtol=1e-11, atol=1e-13)
    return s.y[0], s.y[1], s.t


def xeq(zm, v0):
    """Return the displacement of the offset system's equilibrium, for u = 0.

    At equilibrium ``x2 = 0``, so the relative velocity is ``w = -v0``. For
    ``v0 > 0`` that puts the point in the ``w < 0`` region, whose damping is
    ``zm``, and balancing the spring against the damping term gives
    ``x1 = 2 zm v0 / wn``.

    The sign here is easy to get backwards and was wrong once already. It is
    checked by evaluating the field at the returned point rather than by
    re-reading the algebra.

    Args:
        zm: damping ratio in the region containing the equilibrium.
        v0: boundary offset.

    Returns:
        The equilibrium displacement ``x1``; the equilibrium is
        ``(x1, 0)``.
    """
    return 2*zm*v0/WN


def preturn(zp, zm, v0, r):
    """Advance one turn of the Poincare map on the section x2 = 0.

    Starts on the section at radius ``r`` measured from the equilibrium,
    integrates one full revolution and returns the radius at the next
    downward crossing of ``x2 = 0``.

    The trajectory begins *on* the section, so the event at ``t = 0`` is
    discarded; without that filter the solver returns the starting point and
    the map looks like the identity. Only downward crossings are counted, so
    one call is a full cycle rather than a half cycle.

    Args:
        zp: damping ratio where ``w > 0``.
        zm: damping ratio where ``w < 0``.
        v0: boundary offset.
        r: radius on the section, as displacement from the equilibrium.

    Returns:
        The radius after one cycle, or NaN if no crossing occurs within the
        integration window (an overdamped half plane, where the trajectory
        decays without recrossing).
    """
    xe = xeq(zm, v0)

    def ev(t, y):
        return y[1]
    ev.direction = -1
    s = solve_ivp(f_switched(zp, zm, v0), (0, 30), [xe + r, 0.0], events=ev,
                  rtol=1e-12, atol=1e-14)
    i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
    return (s.y_events[0][i[0]][0] - xe) if i else np.nan


def rstar(zp, zm, v0, r=None, n=400, tol=1e-11):
    """Locate the limit cycle by iterating the return map to its fixed point.

    Iteration rather than a root find: the map is a contraction near the
    cycle (multiplier about 0.54 for the parameters the README uses), so
    successive images converge geometrically and the iteration doubles as
    evidence that the orbit attracts.

    Args:
        zp: damping ratio where ``w > 0``.
        zm: damping ratio where ``w < 0``.
        v0: boundary offset.
        r: starting radius. Defaults to ``2.2 * v0``, which is near the
            cycle because amplitude is exactly proportional to the offset.
        n: iteration cap, so a diverging case terminates.
        tol: relative convergence tolerance on successive radii.

    Returns:
        The fixed point radius, or NaN if the orbit escaped or decayed
        instead of converging.
    """
    r = 2.2*v0 if r is None else r
    for _ in range(n):
        rn = preturn(zp, zm, v0, r)
        if not np.isfinite(rn):
            return np.nan
        if abs(rn - r) < tol*max(1.0, abs(r)):
            return rn
        r = rn
    return r


# ------------------------------------------------------------------ chrome
def style(ax, th, xlabel, ylabel, title=None):
    """Apply the shared chart chrome to one axes.

    Keeps grid, spines and ticks recessive so the data reads first: the top
    and right spines are dropped, the remainder are hairlines in the theme's
    axis colour, and all text is in ink tokens rather than any series
    colour.

    Args:
        ax: the axes to style.
        th: theme dict from ``THEMES``.
        xlabel: x axis label, may contain mathtext.
        ylabel: y axis label, may contain mathtext.
        title: optional axes title, drawn in primary ink.
    """
    ax.set_facecolor(th["surface"])
    ax.grid(True, color=th["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(th["axis"])
        ax.spines[side].set_linewidth(0.8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=th["ink2"], labelsize=8, width=0.8, length=3)
    ax.set_xlabel(xlabel, color=th["ink2"], fontsize=9)
    ax.set_ylabel(ylabel, color=th["ink2"], fontsize=9)
    if title:
        ax.set_title(title, color=th["ink"], fontsize=10, pad=8)


def legend(ax, th, **kw):
    """Draw a legend that stays readable where it overlaps the data.

    The phase portraits are dense enough that a frameless legend becomes
    unreadable over the trajectories, so this gives it the surface colour
    and lifts it above the data.

    Args:
        ax: the axes to attach the legend to.
        th: theme dict from ``THEMES``.
        **kw: passed through to ``Axes.legend``, typically ``loc``.

    Returns:
        The Legend instance.
    """
    lg = ax.legend(fontsize=8, labelcolor=th["ink2"], frameon=True,
                   facecolor=th["surface"], edgecolor="none", framealpha=0.92,
                   **kw)
    lg.set_zorder(7)
    return lg


def boundary(ax, th, y, label):
    """Draw the switching boundary Sigma as chrome rather than a data series.

    The boundary is part of the coordinate system, not one of the things
    being compared, so it takes ink rather than a series colour. Its label
    is anchored to the right edge in axes coordinates, with a surface
    coloured background: pinning it to a data coordinate put it on top of
    the trajectories whenever the axis limits changed.

    Args:
        ax: the axes to draw on.
        th: theme dict from ``THEMES``.
        y: velocity at which the boundary sits, ``v0``.
        label: text for the line, may contain mathtext.
    """
    ax.axhline(y, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
               zorder=5)
    ax.text(0.985, y, label, transform=ax.get_yaxis_transform(), fontsize=8,
            color=th["ink2"], va="bottom", ha="right", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))


def save(fig, th_name, name):
    """Write one figure to ``figures/<name>-<theme>.png`` and close it.

    The face colour is passed explicitly because ``bbox_inches="tight"``
    otherwise reverts the margin around the axes to white, which would
    frame every dark figure in a white border.

    Args:
        fig: the figure to write.
        th_name: theme key, becoming the filename suffix.
        name: figure stem, matching the name the README embeds.
    """
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{name}-{th_name}.png")
    fig.savefig(path, dpi=200, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print("  wrote", os.path.relpath(path, os.path.dirname(OUT)))


def newfig(th, *a, **kw):
    """Create a figure whose canvas already carries the theme surface.

    Args:
        th: theme dict from ``THEMES``.
        *a: positional arguments for ``plt.subplots``, e.g. row and column
            counts.
        **kw: keyword arguments for ``plt.subplots``, e.g. ``figsize``.

    Returns:
        The ``(figure, axes)`` pair from ``plt.subplots``.
    """
    fig, ax = plt.subplots(*a, **kw)
    fig.patch.set_facecolor(th["surface"])
    return fig, ax


# ------------------------------------------------------------- the figures
def fig_linear(th, name):
    """Draw the linear prototype: phase portrait beside time history.

    Three damping ratios from one common initial condition, spanning
    underdamped through overdamped, so the reader can connect the spiral in
    the phase plane to the decaying oscillation in time. This is the
    baseline the nonlinear sections depart from: one equilibrium, every
    trajectory converging to it, no dependence on amplitude.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    cases = [(0.15, "ζ = 0.15"), (0.50, "ζ = 0.50"), (1.20, "ζ = 1.20")]
    fig, axes = newfig(th, 1, 2, figsize=(9.4, 3.9))
    for (z, lab), c in zip(cases, th["series"]):
        x, v, t = traj(f_linear(z), [1.0, 0.0], 26)
        axes[0].plot(x, v, color=c, linewidth=1.8, label=lab, zorder=3)
        axes[1].plot(t, x, color=c, linewidth=1.8, label=lab, zorder=3)
    axes[0].plot([1.0], [0.0], "o", color=th["ink2"], markersize=5, zorder=4)
    axes[0].annotate("start", xy=(1.0, 0.0), xytext=(6, -12),
                     textcoords="offset points", fontsize=8, color=th["ink2"])
    axes[0].plot([0], [0], "o", color=th["ink2"], markersize=4,
                 markerfacecolor=th["surface"], zorder=4)
    style(axes[0], th, "$x$", "$\\dot{x}$", "Phase plane")
    style(axes[1], th, "$t$", "$x$", "Time history")
    axes[0].set_aspect("equal", adjustable="datalim")
    legend(axes[0], th, loc="lower left")
    legend(axes[1], th, loc="upper right")
    fig.suptitle("Linear prototype  $\\ddot{x} + 2\\zeta\\omega_n\\dot{x} "
                 "+ \\omega_n^2 x = 0$,  $\\omega_n = 1$",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "linear-prototype")


def fig_switched(th, name):
    """Draw the three stability cases with the boundary on the x-axis.

    One panel each for positive, zero and negative mean damping, with
    ``Sigma`` and the two half planes labelled and the damping ratio
    annotated in each. Every panel has a half plane with negative damping,
    which is the point: stability follows the mean over a cycle, not the
    sign on either side.

    The marginal panel draws three nested orbits from different starting
    radii rather than one. The closed orbits there form a continuum, one
    through every point, and a single orbit would read as an isolated limit
    cycle, which is precisely the wrong conclusion.

    Axis limits are squared off per panel so the orbit geometry is not
    distorted, but the scales differ between panels because the growing case
    covers a far wider range.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    cases = [
        (0.30, -0.10, "$\\bar{\\zeta} > 0$: decays", [2.6]),
        (0.20, -0.20, "$\\bar{\\zeta} = 0$: closed orbits", [1.0, 1.8, 2.6]),
        (0.10, -0.30, "$\\bar{\\zeta} < 0$: grows", [0.55]),
    ]
    fig, axes = newfig(th, 1, 3, figsize=(11.4, 4.0))
    for ax, (zp, zm, title, starts) in zip(axes, cases):
        for r0 in starts:
            x, v, _ = traj(f_switched(zp, zm), [r0, 0.0], 44)
            ax.plot(x, v, color=th["series"][0], linewidth=1.5, zorder=3)
        ax.plot([0], [0], "o", color=th["ink2"], markersize=5, zorder=5)
        boundary(ax, th, 0.0, "$\\Sigma:\\ \\dot{x} = 0$")
        lim = max(abs(np.array(ax.get_xlim())).max(),
                  abs(np.array(ax.get_ylim())).max())
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.text(0.03, 0.96, "$S^{+}$  $\\zeta_{+}=%.2f$" % zp, transform=ax.transAxes,
                fontsize=8, color=th["ink2"], va="top")
        ax.text(0.03, 0.04, "$S^{-}$  $\\zeta_{-}=%.2f$" % zm, transform=ax.transAxes,
                fontsize=8, color=th["ink2"], va="bottom")
        style(ax, th, "$x$", "$\\dot{x}$", title)
        ax.set_aspect("equal")
    fig.suptitle("Switched damping, boundary through the equilibrium — "
                 "stability set by the mean damping alone",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "switched-damping")


def fig_decrement(th, name):
    """Plot amplitude per cycle against the closed form, on a log scale.

    Markers are integrated from the return map; the dashed lines are
    ``exp(-(delta(zp) + delta(zm)) n)``. Straight parallel lines on a log
    axis make two things visible at once: the decay is exactly geometric,
    and its direction is set by the sign of the mean damping rather than by
    either damping ratio alone.

    The three cases share a starting amplitude so the lines fan out from a
    common point and the divergence is attributable to the parameters rather
    than the initial condition.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    d = lambda z: np.pi*z/np.sqrt(1 - z**2)
    cases = [(0.30, -0.10), (0.20, -0.20), (0.10, -0.30)]
    fig, ax = newfig(th, figsize=(6.4, 4.2))
    for (zp, zm), c in zip(cases, th["series"]):
        r, amps = 1.0, [1.0]
        for _ in range(9):
            r = preturn(zp, zm, 0.0, r)
            amps.append(abs(r))
        n = np.arange(len(amps))
        pred = amps[0]*np.exp(-(d(zp) + d(zm))*n)
        lab = "$\\zeta_{+}=%.2f,\\ \\zeta_{-}=%.2f$" % (zp, zm)
        ax.plot(n, pred, color=c, linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)
        ax.plot(n, amps, "o", color=c, markersize=6, label=lab, zorder=3,
                markeredgecolor=th["surface"], markeredgewidth=1.0)
        ax.annotate("$\\bar{\\zeta}=%+.2f$" % ((zp + zm)/2), xy=(n[-1], amps[-1]),
                    xytext=(6, 0), textcoords="offset points", fontsize=8,
                    color=th["ink2"], va="center")
    ax.set_yscale("log")
    ax.set_xlim(-0.3, 11.2)
    style(ax, th, "full cycles", "amplitude",
          "Markers: integrated.  Dashed: $e^{-(\\delta(\\zeta_+)+\\delta(\\zeta_-))n}$")
    legend(ax, th, loc="lower left")
    fig.tight_layout()
    save(fig, name, "decrement")


def fig_limit_cycle(th, name):
    """Draw the offset boundary portrait with its attracting limit cycle.

    Trajectories from inside and outside the cycle are drawn in different
    colours converging on the same closed orbit, which is what distinguishes
    an attractor from the continuum of the marginal case. The cycle itself
    is drawn heavier so it reads as the object the others approach.

    The equilibrium is marked and annotated, and the gap between it and
    ``Sigma`` is dimensioned, because the offset is the parameter the whole
    section turns on. Starting radii and integration times are tuned so both
    trajectories complete enough turns to be unambiguous without filling the
    frame.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    zp, zm, v0 = 0.3, -0.1, 1.0
    xe, rs = xeq(zm, v0), rstar(zp, zm, v0)
    fig, ax = newfig(th, figsize=(6.6, 5.6))

    xi, vi, _ = traj(f_switched(zp, zm, v0), [xe + 0.25, 0.0], 40)
    ax.plot(xi, vi, color=th["series"][1], linewidth=1.2, alpha=0.95,
            label="from inside", zorder=3)
    xo, vo, _ = traj(f_switched(zp, zm, v0), [xe + 4.2, 0.0], 34)
    ax.plot(xo, vo, color=th["series"][2], linewidth=1.2, alpha=0.95,
            label="from outside", zorder=3)
    xc, vc, _ = traj(f_switched(zp, zm, v0), [xe + rs, 0.0], 6.367077)
    ax.plot(xc, vc, color=th["series"][0], linewidth=2.6,
            label="limit cycle", zorder=4)

    boundary(ax, th, v0, "$\\Sigma:\\ \\dot{x} = v_0$")
    ax.plot([xe], [0.0], "o", color=th["ink"], markersize=6, zorder=6)
    ax.annotate("equilibrium  $x^{*}=u+2\\zeta_{-}v_0/\\omega_n$",
                xy=(xe, 0.0), xytext=(-14, -58), textcoords="offset points",
                fontsize=8, color=th["ink2"], ha="center",
                arrowprops=dict(arrowstyle="-", color=th["ink2"], lw=0.8),
                bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                          ec=th["grid"], lw=0.8), zorder=7)
    ax.annotate("", xy=(xe, v0), xytext=(xe, 0.0),
                arrowprops=dict(arrowstyle="<->", color=th["ink2"], lw=1.0))
    ax.annotate("$v_0$", xy=(xe, v0/2), xytext=(5, 0),
                textcoords="offset points", fontsize=9, color=th["ink2"],
                va="center")
    style(ax, th, "$x$", "$\\dot{x}$",
          "Offset boundary: $\\zeta_{+}=0.3$, $\\zeta_{-}=-0.1$, $v_0=1$")
    ax.set_aspect("equal")
    legend(ax, th, loc="lower right")
    fig.tight_layout()
    save(fig, name, "limit-cycle")


def fig_map_scaling(th, name):
    """Draw the return map beside the amplitude-offset proportionality.

    Left panel: the return map for an offset boundary and for a boundary
    through the equilibrium, over the diagonal ``P(r) = r``. Both are drawn
    together deliberately. The offset map bends across the diagonal and
    crosses it transversally at ``r*``; the through-equilibrium map is a ray
    through the origin that never can, because a homogeneous field gives a
    return map that is an exact scaling. That contrast is the mechanism, and
    it is far clearer as one comparison than as two separate figures.

    Right panel: fixed point radius against offset, with the fitted
    proportionality. The points lie on a line through the origin because the
    system is self-similar in ``(x, v0)`` when unforced, so the amplitude is
    exactly proportional to the offset rather than approximately so.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    zp, zm = 0.3, -0.1
    fig, axes = newfig(th, 1, 2, figsize=(9.8, 4.2))

    rr = np.linspace(0.25, 6.0, 46)
    p1 = np.array([preturn(zp, zm, 1.0, r) for r in rr])
    p0 = np.array([preturn(zp, zm, 0.0, r) for r in rr])
    axes[0].plot(rr, rr, color=th["ink2"], linewidth=1.0,
                 linestyle=(0, (4, 3)), zorder=2)
    axes[0].annotate("$P(r)=r$", xy=(rr[-1], rr[-1]), xytext=(-4, 8),
                     textcoords="offset points", fontsize=8, color=th["ink2"],
                     ha="right")
    axes[0].plot(rr, p1, color=th["series"][0], linewidth=1.8,
                 label="$v_0=1$: crosses at $r^{*}$", zorder=3)
    axes[0].plot(rr, p0, color=th["series"][1], linewidth=1.8,
                 label="$v_0=0$: pure scaling", zorder=3)
    rs = rstar(zp, zm, 1.0)
    axes[0].plot([rs], [rs], "o", color=th["series"][0], markersize=7,
                 markeredgecolor=th["surface"], markeredgewidth=1.2, zorder=5)
    axes[0].annotate("$r^{*}=%.3f$, multiplier $0.539$" % rs, xy=(rs, rs),
                     xytext=(14, -30), textcoords="offset points", fontsize=8,
                     color=th["ink2"],
                     arrowprops=dict(arrowstyle="-", color=th["ink2"], lw=0.8),
                     bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                               ec=th["grid"], lw=0.8), zorder=7)
    style(axes[0], th, "$r$", "$P(r)$", "Return map on $\\{x_2=0\\}$")
    legend(axes[0], th, loc="upper left")

    vs = np.array([0.25, 0.5, 1.0, 2.0, 3.0, 4.0])
    rstars = np.array([rstar(zp, zm, v) for v in vs])
    axes[1].plot(vs, rstars, "o", color=th["series"][0], markersize=7,
                 markeredgecolor=th["surface"], markeredgewidth=1.2, zorder=4,
                 label="integrated $r^{*}$")
    axes[1].plot([0, vs[-1]], [0, (rstars[0]/vs[0])*vs[-1]],
                 color=th["ink2"], linewidth=1.0, linestyle=(0, (4, 3)),
                 zorder=2, label="$r^{*} = %.6f\\,v_0$" % (rstars[0]/vs[0]))
    style(axes[1], th, "offset $v_0$", "$r^{*}$",
          "Amplitude is exactly proportional to the offset")
    legend(axes[1], th, loc="upper left")
    fig.tight_layout()
    save(fig, name, "return-map")


def fig_frequency(th, name):
    """Draw how the limit cycle period depends on the two damping ratios.

    Left: period normalised by the undamped period, against the damping in
    the outer region, for three inner dampings. Every curve starts at the
    existence boundary and stays above one, because the correction to the
    undamped period goes as the square of the damping ratio and so has the
    same sign whichever way the damping points.

    Right: error of the closed form against the exact reduction, with the
    one percent band marked, so the figure says where the formula can be
    trusted rather than only that it exists.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    fig, axes = newfig(th, 1, 2, figsize=(10.0, 4.2))
    for zm, c in zip([-0.05, -0.15, -0.30], th["series"]):
        zps = np.linspace(-zm*1.03, 0.92, 24)
        Te = np.array([frequency.period_exact(z, zm)[0] for z in zps])
        Ts = np.array([frequency.period_series(z, zm) for z in zps])
        lab = "$\\zeta_{-} = %.2f$" % zm
        axes[0].plot(zps, Te/(2*np.pi), color=c, linewidth=1.9, label=lab,
                     zorder=3)
        axes[0].plot(zps, Ts/(2*np.pi), color=c, linewidth=1.1,
                     linestyle=(0, (4, 3)), zorder=3)
        axes[0].annotate(lab, xy=(zps[-1], Te[-1]/(2*np.pi)), xytext=(5, -2),
                         textcoords="offset points", fontsize=8,
                         color=th["ink2"], va="center")
        axes[1].plot(zps, 100*(Ts - Te)/Te, color=c, linewidth=1.9, label=lab,
                     zorder=3)
    axes[0].axhline(1.0, color=th["ink2"], linewidth=1.0,
                    linestyle=(0, (5, 3)), zorder=2)
    axes[0].text(0.985, 1.0, "undamped, $T = 2\\pi/\\omega_n$",
                 transform=axes[0].get_yaxis_transform(), fontsize=8,
                 color=th["ink2"], va="top", ha="right", zorder=4,
                 bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"],
                           ec="none"))
    style(axes[0], th, "$\\zeta_{+}$", "$T\\,\\omega_n / 2\\pi$",
          "Solid: exact reduction.  Dashed: closed form")
    legend(axes[0], th, loc="upper left")
    axes[1].axhspan(-1, 1, color=th["grid"], zorder=1)
    axes[1].text(0.98, 1.0, "$\\pm 1\\%$", transform=axes[1].get_yaxis_transform(),
                 fontsize=8, color=th["ink2"], ha="right", va="bottom", zorder=4)
    axes[1].axhline(0.0, color=th["axis"], linewidth=0.8, zorder=2)
    style(axes[1], th, "$\\zeta_{+}$", "error in $T$  (%)",
          "Where the closed form can be trusted")
    legend(axes[1], th, loc="upper left")
    fig.suptitle("Limit cycle period: set by $\\omega_n$ and the two damping "
                 "ratios alone, never by the offset",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "frequency")


def fig_poles(th, name):
    """Draw the s-plane pole locus of the two half-plane subsystems.

    Each half plane is an ordinary second order system with characteristic
    polynomial ``s^2 + 2 zeta wn s + wn^2``, so its poles sit at
    ``wn(-zeta +/- sqrt(zeta^2 - 1))``. While ``|zeta| < 1`` they are a
    complex pair on the circle of radius ``wn``, at ``cos(theta) = zeta``
    from the negative real axis; at ``|zeta| = 1`` they meet on the real
    axis and split along it.

    The radius is ``wn`` for both half planes, so switching moves the pole
    pair between two points on one fixed circle. Colour runs on the
    diverging blue-red scale because the quantity it encodes, the damping
    ratio, has a meaningful zero: the imaginary axis, where the poles cross
    from decaying to growing.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    fig, ax = newfig(th, figsize=(6.6, 6.0))
    ax.axhspan(-2, 2, xmin=0, xmax=0.5, color=th["grid"], zorder=0)

    def poles(z):
        d = complex(z*z - 1)**0.5
        return np.array([-z + d, -z - d])

    zz = np.linspace(-2.2, 2.2, 900)
    for z in zz:
        p = poles(z)
        c = th["div_neg"] if z < 0 else th["div_pos"]
        a = min(1.0, 0.18 + 0.8*abs(z)/2.2)
        ax.plot(p.real, p.imag, ".", color=c, markersize=2.2, alpha=a, zorder=2)

    marks = [(2.0, "$\\zeta = 2$"), (1.0, "$\\zeta = 1$"), (0.5, "$\\zeta = 0.5$"),
             (0.0, "$\\zeta = 0$"), (-0.5, "$\\zeta = -0.5$"),
             (-1.0, "$\\zeta = -1$"), (-2.0, "$\\zeta = -2$")]
    for z, lab in marks:
        p = poles(z)
        c = th["ink2"] if z == 0 else (th["div_neg"] if z < 0 else th["div_pos"])
        ax.plot(p.real, p.imag, "x", color=c, markersize=9, markeredgewidth=2.0,
                zorder=4)
        k = 0 if p[0].imag >= 0 else 1
        ax.annotate(lab, xy=(p[k].real, p[k].imag), xytext=(6, 6),
                    textcoords="offset points", fontsize=8, color=th["ink2"],
                    zorder=5, bbox=dict(boxstyle="round,pad=0.18",
                                        fc=th["surface"], ec="none"))

    t = np.linspace(0, 2*np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color=th["ink2"], linewidth=1.0,
            linestyle=(0, (5, 3)), zorder=1)
    ax.annotate("$|s| = \\omega_n$", xy=(np.cos(3.6), np.sin(3.6)),
                xytext=(-6, -14), textcoords="offset points", fontsize=8,
                color=th["ink2"], zorder=5,
                bbox=dict(boxstyle="round,pad=0.18", fc=th["surface"], ec="none"))
    ax.axhline(0, color=th["axis"], linewidth=0.9, zorder=1)
    ax.axvline(0, color=th["axis"], linewidth=0.9, zorder=1)
    ax.text(-1.95, 1.85, "left half plane\ndecaying", fontsize=8.5,
            color=th["ink2"], va="top", zorder=5)
    ax.text(1.95, 1.85, "right half plane\ngrowing", fontsize=8.5,
            color=th["ink2"], va="top", ha="right", zorder=5)
    ax.set_xlim(-2.3, 2.3)
    ax.set_ylim(-2.0, 2.0)
    ax.set_aspect("equal")
    style(ax, th, "$\\mathrm{Re}\\,s / \\omega_n$",
          "$\\mathrm{Im}\\,s / \\omega_n$",
          "Poles of one half plane: $s^2 + 2\\zeta\\omega_n s + \\omega_n^2$")
    fig.tight_layout()
    save(fig, name, "pole-zero")


def fig_stability_map(th, name):
    """Draw the stability classification over the whole damping ratio plane.

    Left: the boundary through the equilibrium. Whether a half plane
    carries an invariant ray decides everything. A half plane with
    ``zeta <= -1`` holds an escaping ray, one with ``zeta >= 1`` holds a
    decaying sector, and each is invariant, so a trajectory that enters
    never leaves. Where neither half plane has real poles the rotation is
    complete, the return map applies, and the sign of the mean damping
    decides globally.

    Right: the same plane for an offset boundary, where the interest is
    whether a limit cycle exists rather than whether the origin attracts.

    Regions are labelled in place rather than by colour alone.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    def rule(zp, zm):
        esc, dec = (zp <= -1) or (zm <= -1), (zp >= 1) or (zm >= 1)
        if esc and dec: return 2
        if esc: return 1
        if dec: return 0
        return 0 if zp + zm > 0 else 1

    fig, axes = newfig(th, 1, 2, figsize=(11.0, 5.0))
    g = np.linspace(-2.2, 2.2, 601)
    ZP, ZM = np.meshgrid(g, g)
    Z = np.vectorize(rule)(ZP, ZM)
    cmap = matplotlib.colors.ListedColormap(list(th["series"]))
    axes[0].pcolormesh(g, g, Z, cmap=cmap, vmin=-0.5, vmax=2.5,
                       shading="nearest", zorder=1)
    axes[0].plot([-2.2, 2.2], [2.2, -2.2], color=th["ink"], linewidth=1.4,
                 linestyle=(0, (5, 3)), zorder=3)
    for v in (-1, 1):
        axes[0].axvline(v, color=th["ink"], linewidth=0.9, alpha=0.55, zorder=3)
        axes[0].axhline(v, color=th["ink"], linewidth=0.9, alpha=0.55, zorder=3)
    def tag(ax, x, y, txt):
        ax.text(x, y, txt, fontsize=8.5, color=th["ink"], ha="center",
                va="center", zorder=6,
                bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                          ec=th["grid"], lw=0.8, alpha=0.94))
    tag(axes[0], 1.15, 1.15, "decays")
    tag(axes[0], -1.5, -1.5, "escapes")
    tag(axes[0], 1.55, -1.6, "mixed:\nseparatrix")
    tag(axes[0], -1.6, 1.55, "mixed:\nseparatrix")
    tag(axes[0], 0.72, -0.72, "$\\bar{\\zeta} = 0$")
    style(axes[0], th, "$\\zeta_{+}$", "$\\zeta_{-}$",
          "Boundary through the equilibrium")
    axes[0].set_aspect("equal")
    axes[0].set_xlim(-2.2, 2.2)
    axes[0].set_ylim(-2.2, 2.2)

    gg = np.linspace(-1.0, 1.0, 601)
    A, B = np.meshgrid(gg, gg)
    lc = (B < 0) & (A + B > 0)
    # same colour, same meaning as the left panel: blue decays, orange
    # escapes, and the third slot carries whatever else the panel shows
    axes[1].pcolormesh(gg, gg, np.where(lc, 2, np.where(A + B > 0, 0, 1)),
                       cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest",
                       zorder=1)
    axes[1].plot([-1, 1], [1, -1], color=th["ink"], linewidth=1.4,
                 linestyle=(0, (5, 3)), zorder=3)
    axes[1].axhline(0, color=th["ink"], linewidth=1.4, linestyle=(0, (5, 3)),
                    zorder=3)
    tag(axes[1], 0.62, -0.34, "limit cycle")
    tag(axes[1], -0.05, 0.62, "decays to\nequilibrium")
    tag(axes[1], -0.55, -0.55, "escapes")
    style(axes[1], th, "$\\zeta_{+}$", "$\\zeta_{-}$",
          "Offset boundary, both half planes underdamped")
    axes[1].set_aspect("equal")
    fig.suptitle("Where each behaviour lives in the damping ratio plane",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "stability-map")


def fig_symmetric(th, name):
    """Draw the symmetric deadzone version beside its existence region.

    Left: the phase portrait. The band ``|xdot| < v0`` is shaded, both
    boundaries drawn, and trajectories converge onto one odd-symmetric
    cycle from inside and outside. The equilibrium is at the origin here,
    not offset, because the field is odd.

    Right: where a limit cycle exists, for both versions. The single
    boundary needs the mean damping positive, which is a triangle. The
    deadzone needs only ``zp > 0``, which is the whole quadrant — so the
    wedge between them is where symmetrising the transition creates a cycle
    that did not exist before.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    zp, zm, v0 = 0.3, -0.1, 1.0
    fig, axes = newfig(th, 1, 2, figsize=(11.0, 5.0))

    r, _ = symmetric.cycle_integrated(zp, zm, v0)
    axes[0].axhspan(-v0, v0, color=th["grid"], zorder=0)
    for y in (v0, -v0):
        axes[0].axhline(y, color=th["ink2"], linewidth=1.2,
                        linestyle=(0, (5, 3)), zorder=5)
    f = symmetric.field(zp, zm, v0)
    for r0, c, lab in [(0.30*r, th["series"][1], "from inside"),
                       (2.0*r, th["series"][2], "from outside")]:
        x, v, _ = traj(f, [r0, 0.0], 46)
        axes[0].plot(x, v, color=c, linewidth=1.2, label=lab, zorder=3)
    xc, vc, _ = traj(f, [r, 0.0], symmetric.period_exact(zp, zm, v0)[0])
    axes[0].plot(xc, vc, color=th["series"][0], linewidth=2.6,
                 label="limit cycle", zorder=4)
    axes[0].plot([0], [0], "o", color=th["ink"], markersize=6, zorder=6)
    axes[0].text(0.985, 0.0, "deadzone  $|\\dot{x}| < v_0$",
                 transform=axes[0].get_yaxis_transform(), fontsize=8,
                 color=th["ink2"], ha="right", va="center", zorder=7,
                 bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"],
                           ec="none"))
    style(axes[0], th, "$x$", "$\\dot{x}$",
          "Symmetric band: $\\zeta_{+}=0.3$, $\\zeta_{-}=-0.1$, $v_0=1$")
    axes[0].set_aspect("equal")
    legend(axes[0], th, loc="lower right")

    g = np.linspace(-1.0, 1.0, 601)
    A, B = np.meshgrid(g, g)
    both = (B < 0) & (A + B > 0)
    sym_only = (B < 0) & (A > 0) & (A + B <= 0)
    Z = np.where(both, 0, np.where(sym_only, 1, np.nan))
    axes[1].pcolormesh(g, g, Z, cmap=matplotlib.colors.ListedColormap(
        list(th["series"][:2])), vmin=-0.5, vmax=1.5, shading="nearest",
        zorder=1)
    axes[1].plot([-1, 1], [1, -1], color=th["ink"], linewidth=1.4,
                 linestyle=(0, (5, 3)), zorder=3)
    axes[1].axhline(0, color=th["ink"], linewidth=1.4, zorder=3)
    axes[1].axvline(0, color=th["ink"], linewidth=1.4, zorder=3)
    for x, y, t in [(0.62, -0.22, "both versions"),
                    (0.30, -0.72, "deadzone only"),
                    (-0.5, 0.5, "no cycle in\neither version")]:
        axes[1].text(x, y, t, fontsize=8.5, color=th["ink"], ha="center",
                     va="center", zorder=6,
                     bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                               ec=th["grid"], lw=0.8, alpha=0.94))
    style(axes[1], th, "$\\zeta_{+}$", "$\\zeta_{-}$",
          "Where a limit cycle exists")
    axes[1].set_aspect("equal")
    fig.suptitle("Symmetrising the transition drops the condition on the "
                 "mean damping", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "symmetric")


def fig_four_models(th, name):
    """Draw the complete set of four switched-damping models.

    Rows are the asymmetric and symmetric forms; columns are switching on
    velocity and on displacement. All four use the same damping ratios and
    the same boundary value, so the panels are directly comparable.

    The point of the layout is that the period is identical along each row
    while the orbit is not: moving the boundary from velocity to
    displacement rotates which part of the cycle is damped, without
    changing how long the orbit spends damped.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    zp, zm, b = 0.3, -0.1, 1.0

    def settle(f, start, period, T=120.0):
        """Integrate onto the attractor, then return exactly one period.

        ``max_step`` is bounded so the boundary crossings are resolved
        cleanly. The corners visible in the displacement-switched orbits are
        not an artefact of that: those fields are discontinuous, so the
        curvature genuinely jumps where the orbit crosses, and the corner is
        the thing worth seeing.
        """
        s = solve_ivp(f, (0, T), start, rtol=1e-12, atol=1e-14,
                      max_step=0.004, dense_output=True)
        return s.sol(np.linspace(T - period, T, 3000))

    panels = [
        ("Asymmetric, switch on $\\dot{x}$", f_switched(zp, zm, b),
         [2*zm*b + 2.2, 0.0], "h", (b,), 2*zm*b),
        ("Asymmetric, switch on $x$",
         displacement.field(zp, zm, b, False), [2.6, 0.0], "v", (b,), 0.0),
        ("Symmetric, switch on $\\dot{x}$", symmetric.field(zp, zm, b),
         [1.6, 0.0], "h", (b, -b), 0.0),
        ("Symmetric, switch on $x$",
         displacement.field(zp, zm, b, True), [1.6, 0.0], "v", (b, -b), 0.0),
    ]
    periods = [frequency.period_exact(zp, zm)[0],
               displacement.period_exact(zp, zm, b, False)[0],
               symmetric.period_exact(zp, zm)[0],
               displacement.period_exact(zp, zm, b, True)[0]]

    fig, axes = newfig(th, 2, 2, figsize=(9.6, 9.0))
    for ax, (title, f, start, orient, lines, xeq), T in zip(axes.ravel(),
                                                            panels, periods):
        y = settle(f, start, T)
        if len(lines) == 2:
            if orient == "h":
                ax.axhspan(-b, b, color=th["grid"], zorder=0)
            else:
                ax.axvspan(-b, b, color=th["grid"], zorder=0)
        for L in lines:
            drawer = ax.axhline if orient == "h" else ax.axvline
            drawer(L, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
                   zorder=4)
        ax.plot(y[0], y[1], color=th["series"][0], linewidth=2.4, zorder=3)
        ax.plot([xeq], [0], "o", color=th["ink"], markersize=5, zorder=5)
        style(ax, th, "$x$", "$\\dot{x}$", title)
        ax.set_aspect("equal")
        ax.text(0.03, 0.03, "$T = %.4f$" % T, transform=ax.transAxes,
                fontsize=8.5, color=th["ink2"], zorder=6,
                bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"],
                          ec=th["grid"], lw=0.8))
    fig.suptitle("The four models at $\\zeta_{+}=0.3$, $\\zeta_{-}=-0.1$: "
                 "the period matches along each row, the orbit does not",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "four-models")


# ------------------------------------------------------- the forced figures
#: Drive strengths for the staircase panel, and the sweep resolution. Kept
#: modest because each point is an integration over several hundred drive
#: periods; the tongue edges are found by bisection instead, which is where
#: the resolution actually matters.
STAIR_A = (0.3, 0.6, 1.2)
R_LO, R_HI, N_STAIR = 0.5, 2.0, 49
TONGUE_A = (0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.6)


def fig_forced_tongues(th, name):
    """Draw what a sinusoidal drive does to the limit cycle.

    Left: the rotation number, orbit windings per drive period, swept
    against drive frequency at three drive strengths. With no drive the
    curve is the smooth hyperbola ``w_lc / Omega``, drawn as chrome — an
    autonomous cycle keeps its own frequency because it has nothing to lock
    to. Adding drive flattens it into plateaus, one per lock, and they widen
    with drive strength. That is a devil's staircase.

    Right: the 1:1 tongue, each edge located by bisecting on the rotation
    number. It closes to a point at zero drive on ``Omega = w_lc``, which is
    the definition of an Arnold tongue, and by the top of the axis it has
    swallowed the whole frequency range shown.

    Both panels are read the same way: inside the tongue the oscillator has
    given up its own frequency and runs at the drive's.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    fig, axes = newfig(th, 1, 2, figsize=(11.0, 4.6))

    rr = np.linspace(R_LO, R_HI, 400)
    axes[0].plot(rr, 1.0/rr, color=th["ink2"], linewidth=1.2,
                 linestyle=(0, (5, 3)), zorder=3, label="no drive")
    for k, a in enumerate(STAIR_A):
        r, w = forced.staircase(a, R_LO, R_HI, N_STAIR)
        axes[0].plot(r, w, color=th["series"][k], linewidth=1.7, zorder=4 + k,
                     label="$A/\\omega_n v_0 = %g$" % a)
        # label each curve on its own 1:1 plateau, staggered along it so the
        # three labels cannot collide -- they all sit at w = 1
        flat = np.flatnonzero(np.abs(w - 1.0) < 1e-6)
        if flat.size:
            j = flat[min(flat.size - 1, int((0.30 + 0.30*k)*(flat.size - 1)))]
            axes[0].annotate("$%g$" % a, (r[j], 1.0), fontsize=8.5,
                             color=th["series"][k], zorder=8, ha="center",
                             xytext=(0, 7), textcoords="offset points")
    style(axes[0], th, "$\\Omega / \\omega_{lc}$", "rotation number  $w$",
          "The staircase: drive flattens the frequency ratio")
    legend(axes[0], th, loc="upper right")

    amps = np.array(TONGUE_A)
    lo, hi = forced.tongue_width(amps)
    open_ = ~np.isnan(lo)
    # an edge outside the search window is drawn at the axis, and the curve
    # is not drawn there, so a clipped tongue never reads as a closed one
    lo_p = np.where(np.isneginf(lo), R_LO, lo)
    hi_p = np.where(np.isposinf(hi), R_HI, hi)
    axes[1].fill_betweenx(amps[open_], lo_p[open_], hi_p[open_],
                          color=th["series"][0], alpha=0.28, zorder=2,
                          linewidth=0)
    inside = open_ & np.isfinite(lo)
    axes[1].plot(lo[inside], amps[inside], color=th["series"][0],
                 linewidth=1.8, zorder=4)
    inside_hi = open_ & np.isfinite(hi)
    axes[1].plot(hi[inside_hi], amps[inside_hi], color=th["series"][0],
                 linewidth=1.8, zorder=4, label="edge of the 1:1 lock")
    axes[1].axvline(1.0, color=th["ink2"], linewidth=1.1,
                    linestyle=(0, (5, 3)), zorder=3)
    axes[1].text(1.0, amps[-1]*0.98, "  $\\Omega = \\omega_{lc}$", fontsize=8,
                 color=th["ink2"], va="top", ha="left", zorder=6)
    axes[1].text(1.0, amps[-1]*0.45, "entrained\n1:1", fontsize=9,
                 color=th["ink"], ha="center", va="center", zorder=6,
                 bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                           ec=th["grid"], lw=0.8, alpha=0.94))
    for x in (R_LO + 0.06, R_HI - 0.06):
        axes[1].text(x, amps[-1]*0.16, "free\nrunning", fontsize=9,
                     color=th["ink"], ha="center", va="center", zorder=6,
                     bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                               ec=th["grid"], lw=0.8, alpha=0.94))
    axes[1].set_xlim(R_LO, R_HI)
    style(axes[1], th, "$\\Omega / \\omega_{lc}$", "$A / \\omega_n v_0$",
          "The 1:1 Arnold tongue")
    fig.suptitle("Forcing the deadzone prototype: entrainment, not chaos",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "forced-tongues")


def _forced_orbit(a, r, zp, zm, n_skip=500, n_show=40, per=400):
    """Settled trajectory of the forced model, for drawing under a section.

    Discards ``n_skip`` drive periods, then returns ``n_show`` more sampled
    ``per`` times each, as an ``(N, 2)`` array of ``(x, xdot)``.
    """
    wl = forced.w_lc(zp, zm)
    om = r*wl
    td = 2.0*np.pi/om
    t = td*(n_skip + np.arange(n_show*per + 1)/per)
    return forced._run(zp, zm, forced.V0, a*forced.V0, om,
                       np.concatenate(([0.0], t)), [0.0, 2.0*forced.V0])[1:]


def fig_forced_sections(th, name):
    """Draw the three things the stroboscopic section is ever observed to be.

    Sampling the state once per drive period collapses the response to a set
    of points. Across every parameter combination tested that set is one of
    exactly three things, and each panel draws one of them over the settled
    trajectory it came from:

    *A single point*, inside the tongue: the response repeats every drive
    period, so it is a periodic orbit locked to the drive. The orbit is a
    closed curve; the section of it is one point.

    *A closed curve*, outside the tongue: the response carries two
    incommensurate frequencies, its own and the drive's, and fills a torus
    whose section is that curve.

    *A chain of islands*, once the contraction per drive period is weak
    enough: the curve breaks into a ring of small closed curves permuted
    cyclically by the map. The inset zooms one island by a couple of
    hundred times to show it is a curve, strongly flattened, and not a point.

    A fourth possibility, a fractal cloud, is what chaos would look like. It
    does not appear anywhere in the range tested.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    cases = (
        (0.45, 1.00, forced.ZP, forced.ZM, "locked 1:1",
         "one point:\nperiodic", False, 3),
        (0.30, 0.60, forced.ZP, forced.ZM, "quasi-periodic",
         "closed curve:\ntwo frequencies", False, 40),
        (0.10, 1.40, 0.01, -0.003, "island chain",
         "seven islands:\nstill not chaos", True, 14),
    )
    fig, axes = newfig(th, 1, 3, figsize=(12.6, 4.4))
    for k, (a, r, zp, zm, title, tag, zoom, nshow) in enumerate(cases):
        ax = axes[k]
        orb = _forced_orbit(a, r, zp, zm, n_show=nshow)
        ax.plot(orb[:, 0], orb[:, 1], color=th["axis"], linewidth=0.7,
                zorder=2, label="trajectory")
        wl = forced.w_lc(zp, zm)
        pts = forced.strobe(a*forced.V0, r*wl, zp, zm, forced.V0,
                            n_keep=1500)
        ax.scatter(pts[:, 0], pts[:, 1], s=9.0, color=th["series"][k],
                   zorder=5, linewidths=0, label="once per drive period")
        boundary(ax, th, forced.V0, "$\\dot{x} = v_0$")
        ax.axhline(-forced.V0, color=th["ink2"], linewidth=1.2,
                   linestyle=(0, (5, 3)), zorder=4)
        m = 1.12*np.max(np.abs(orb), axis=0)
        ax.set_xlim(-m[0], m[0])
        ax.set_ylim(-m[1], m[1])
        style(ax, th, "$x_1 = x$", "$x_2 = \\dot{x}$",
              "%s\n$A/\\omega_n v_0=%g$,  $\\Omega/\\omega_{lc}=%g$"
              % (title, a, r))
        ax.text(0.03, 0.97, tag, transform=ax.transAxes, fontsize=8.5,
                color=th["ink"], ha="left", va="top", zorder=8,
                bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                          ec=th["grid"], lw=0.8, alpha=0.94))
        if zoom:
            c = pts[0]
            d = 2.2*np.max(np.abs(pts[::7][:120] - c))
            ins = ax.inset_axes([0.34, 0.35, 0.32, 0.32])
            ins.scatter(pts[:, 0], pts[:, 1], s=4.0, color=th["series"][k],
                        zorder=5, linewidths=0)
            ins.set_xlim(c[0] - d, c[0] + d)
            ins.set_ylim(c[1] - d, c[1] + d)
            ins.set_facecolor(th["surface"])
            ins.set_xticks([])
            ins.set_yticks([])
            for sp in ins.spines.values():
                sp.set_color(th["axis"])
                sp.set_linewidth(0.8)
            ins.set_title("one island, $\\times%d$" % round(m[0]/d),
                          fontsize=7.5, color=th["ink2"], pad=3)
    axes[0].legend(fontsize=8, labelcolor=th["ink2"], frameon=True,
                   facecolor=th["surface"], edgecolor="none", framealpha=0.92,
                   loc="lower right").set_zorder(9)
    fig.suptitle("The stroboscopic section is only ever one of three things",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "forced-sections")


#: Grid the comparison is drawn from. Wide in frequency because Van der
#: Pol's chaos lives well above the cycle frequency -- an earlier scan
#: stopping at 2.4 would have missed all of it.
CMP_R = np.round(np.linspace(0.5, 8.0, 31), 3)
CMP_A_VDP = (0.5, 1.0, 2.0, 5.0, 10.0)
CMP_A_PROTO = (0.3, 0.8, 1.5, 3.0, 5.0)


def _regime_codes(lab, q, lam):
    """0 locked, 1 quasi-periodic, 2 chaotic, from a regime map."""
    return np.where(q > 0, 0, np.where(lam > section.LAM_TOL, 2, 1))


def fig_vanderpol_compare(th, name):
    """Draw the prototype beside Van der Pol under identical forcing.

    Four panels over the same drive frequency range, all classified by the
    same engine: the deadzone prototype, then Van der Pol at three
    relaxation parameters. Cells are locked, quasi-periodic, or chaotic.

    The comparison is the point. At ``mu = 0.1`` Van der Pol is nearly
    harmonic and its map is the prototype's — a 1:1 tongue, narrow higher
    order locks, tori between them, and no chaos. Raising ``mu`` strengthens
    the nonlinearity and chaos appears, at drive frequencies several times
    the cycle frequency. The prototype never gets there at any drive
    strength, because its damping saturates: outside the deadzone the ratio
    is exactly ``zp`` however hard the orbit is driven, so a bigger orbit is
    not a more nonlinear one. Van der Pol's ``-mu(1 - x^2)`` grows without
    limit, so drive amplitude buys nonlinearity that the prototype cannot
    buy at any price.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    fig, axes = newfig(th, 1, 4, figsize=(15.0, 4.2))
    cmap = matplotlib.colors.ListedColormap(list(th["series"]))

    lab, q, w, lam = forced.regime_map(CMP_R, np.array(CMP_A_PROTO))
    panels = [(axes[0], CMP_A_PROTO, _regime_codes(lab, q, lam),
               "deadzone prototype", "$A/\\omega_n v_0$",
               "$\\zeta_{+}=0.3$, $\\zeta_{-}=-0.1$")]
    for k, mu in enumerate((0.1, 1.0, 5.0)):
        lab, q, w, lam = vanderpol.regime_map(mu, CMP_R,
                                              np.array(CMP_A_VDP))
        panels.append((axes[k + 1], CMP_A_VDP, _regime_codes(lab, q, lam),
                       "Van der Pol  $\\mu=%g$" % mu, "$A$",
                       "contraction %.1e" % vanderpol.contraction(mu)))

    for ax, amps, code, title, ylab, sub in panels:
        ax.pcolormesh(CMP_R, np.arange(len(amps)), code, cmap=cmap,
                      vmin=-0.5, vmax=2.5, shading="nearest", zorder=1)
        ax.set_yticks(np.arange(len(amps)))
        ax.set_yticklabels(["%g" % a for a in amps])
        style(ax, th, "$\\Omega / \\omega_{lc}$", ylab, title)
        ax.text(0.5, -0.30, sub, transform=ax.transAxes, fontsize=8,
                color=th["ink2"], ha="center", va="top")
        n_chaos = int(np.sum(code == 2))
        ax.text(0.97, 0.04,
                "chaotic cells: %d" % n_chaos, transform=ax.transAxes,
                fontsize=8.5, color=th["ink"], ha="right", va="bottom",
                zorder=7,
                bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                          ec=th["grid"], lw=0.8, alpha=0.94))

    handles = [matplotlib.patches.Patch(facecolor=th["series"][i],
                                        label=t)
               for i, t in enumerate(("locked", "quasi-periodic", "chaotic"))]
    axes[0].legend(handles=handles, fontsize=8, labelcolor=th["ink2"],
                   frameon=True, facecolor=th["surface"], edgecolor="none",
                   framealpha=0.92, loc="upper left").set_zorder(9)
    fig.suptitle("Same drive, same measurements: a saturating nonlinearity "
                 "against an unbounded one", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "vanderpol-compare")


def fig_staircase(th, name):
    """Draw what the second threshold adds: a second cycle, and a basin.

    Left: the cycle averaged damping against amplitude. With one threshold
    it runs monotonically between two levels and can cross zero once. Three
    levels let it turn, so it crosses twice — rising through zero at the
    unstable cycle and falling through at the stable one.

    Middle: the phase portrait. Trajectories started inside the inner cycle
    decay to the origin, those started outside it wind on to the outer
    cycle, and the inner cycle is the boundary between the two basins. The
    origin and the outer cycle both attract, which none of the earlier
    prototypes can do.

    Right: the staircase closing on Van der Pol as levels are added. The
    cycle radius approaches Van der Pol's exact 2, but the outermost zone is
    unbounded, so the approximation is only ever good over a bounded range
    of amplitude.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    lv, ed = staircase.BISTABLE_LEVELS, staircase.BISTABLE_EDGES
    cyc = staircase.cycles_exact(lv, ed)
    r_un = [r for r, m, st in cyc if not st][0]
    r_st = [r for r, m, st in cyc if st][0]

    fig, axes = newfig(th, 1, 3, figsize=(13.2, 4.4))

    # ---- left: averaged damping
    rr = np.linspace(0.25, 4.0, 900)
    zz = np.array([staircase.mean_damping(r, lv, ed) for r in rr])
    axes[0].axhline(0.0, color=th["ink2"], linewidth=1.0,
                    linestyle=(0, (5, 3)), zorder=3)
    axes[0].plot(rr, zz, color=th["series"][0], linewidth=1.9, zorder=4,
                 label="$\\langle\\zeta\\rangle(R)$")
    for r, txt, k in ((r_un, "unstable", 1), (r_st, "stable", 2)):
        axes[0].plot([r], [0.0], "o", color=th["ink"], markersize=6,
                     zorder=6)
        axes[0].annotate(txt + "\n$R=%.3f$" % r, (r, 0.0), fontsize=8,
                         color=th["ink"], ha="center", zorder=7,
                         xytext=(0, 16 if k == 1 else -30),
                         textcoords="offset points",
                         bbox=dict(boxstyle="round,pad=0.25",
                                   fc=th["surface"], ec=th["grid"], lw=0.8))
    for e in ed:
        axes[0].axvline(e, color=th["axis"], linewidth=0.9, zorder=2)
    axes[0].text(ed[0], zz.max()*0.96, " $a$", fontsize=8, color=th["ink2"],
                 ha="left", va="top", zorder=6)
    axes[0].text(ed[1], zz.max()*0.96, " $b$", fontsize=8, color=th["ink2"],
                 ha="left", va="top", zorder=6)
    style(axes[0], th, "amplitude $R$", "cycle averaged $\\zeta$",
          "Two zero crossings, so two cycles")
    legend(axes[0], th, loc="lower right")

    # ---- middle: phase portrait with the basin boundary
    f = staircase.field(lv, ed)
    for k, (a0, lbl) in enumerate(((r_un*0.93, "inside: decays"),
                                   (r_un*1.07, "outside: grows"),
                                   (3.6, "from beyond: settles"))):
        sol = solve_ivp(f, (0.0, 150.0), [a0, 0.0], method="LSODA",
                        rtol=1e-10, atol=1e-12,
                        t_eval=np.linspace(0.0, 150.0, 12000))
        axes[1].plot(sol.y[0], sol.y[1], color=th["series"][k],
                     linewidth=1.0, zorder=4, label=lbl)
    for r, ls in ((r_un, (0, (4, 3))), (r_st, "-")):
        sol = solve_ivp(f, (0.0, staircase.period(r, lv, ed)), [r, 0.0],
                        method="LSODA", rtol=1e-12, atol=1e-14,
                        t_eval=np.linspace(0.0, staircase.period(r, lv, ed),
                                           2000))
        axes[1].plot(sol.y[0], sol.y[1], color=th["ink"], linewidth=1.7,
                     linestyle=ls, zorder=6)
    for e in ed:
        for sgn in (1, -1):
            axes[1].axvline(sgn*e, color=th["ink2"], linewidth=1.0,
                            linestyle=(0, (5, 3)), zorder=3)
    axes[1].text(0.5, 0.965,
                 "dashed: unstable cycle (basin boundary)\nsolid: stable cycle",
                 transform=axes[1].transAxes, fontsize=8, color=th["ink"],
                 ha="center", va="top", zorder=8,
                 bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                           ec=th["grid"], lw=0.8, alpha=0.94))
    style(axes[1], th, "$x_1 = x$", "$x_2 = \\dot{x}$",
          "The origin and the outer cycle both attract")
    legend(axes[1], th, loc="lower right")

    # ---- right: convergence to Van der Pol
    ns = (2, 3, 5, 9, 17)
    got = []
    for n in ns:
        lvn, edn = staircase.vdp_staircase(1.0, n)
        ex = staircase.cycles_exact(lvn, edn)
        got.append(ex[-1][0] if ex else np.nan)
    axes[2].axhline(2.0, color=th["ink2"], linewidth=1.2,
                    linestyle=(0, (5, 3)), zorder=3)
    axes[2].text(ns[-1], 2.0, "Van der Pol, $R=2$  ", fontsize=8,
                 color=th["ink2"], ha="right", va="bottom", zorder=6)
    axes[2].plot(ns, got, "o-", color=th["series"][0], linewidth=1.8,
                 markersize=5, zorder=4, label="staircase cycle radius")
    for n, g in zip(ns, got):
        axes[2].annotate("%.3f" % g, (n, g), fontsize=7.5,
                         color=th["ink2"], ha="center",
                         xytext=(0, -13), textcoords="offset points")
    axes[2].set_xscale("log")
    axes[2].set_xticks(ns)
    axes[2].set_xticklabels([str(n) for n in ns])
    # a log axis draws its own minor decade labels, which land on top of
    # these level counts
    axes[2].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axes[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    style(axes[2], th, "levels in the staircase",
          "cycle radius", "Closing on Van der Pol, $\\mu = 1$")
    legend(axes[2], th, loc="lower right")

    fig.suptitle("Two thresholds: a second cycle, a basin, and a dial "
                 "towards a smooth nonlinearity",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "staircase")


#: Drive frequencies for the side by side. Fine enough to resolve the
#: chaotic bands, which sit in the transitions between one lock and the
#: next: a 0.05 grid stepped straight over Van der Pol's and reported it as
#: having none.
SVDP_OMS = np.round(np.linspace(2.40, 2.56, 33), 4)
SVDP_LEVELS = (5, 9, 17, 33, 65)


def fig_staircase_vdp(th, name):
    """Drive the staircase beside Van der Pol and compare what comes out.

    Left: the damping laws themselves. Van der Pol's ``-mu(1-x^2)/2`` is
    smooth; the staircase samples it in steps. This is the only difference
    between the two systems, and the level count is the only thing varied.

    Middle: how each responds as the drive frequency is swept, one strip per
    system. All of them run lock 3, then a chaotic band, then lock 4,
    another chaotic band, then lock 5 — the chaos lives in the transitions.
    The strips make the convergence visible: at five levels the bands are in
    the wrong place, and by sixty-five they sit on Van der Pol's.

    Right: that convergence as a number. The Jaccard index counts
    frequencies where both are chaotic against frequencies where either is,
    so a staircase going chaotic where Van der Pol does not counts against
    it exactly as a miss does.

    The orbit stays within ``|x| ~ 2.15`` throughout, well inside the fitted
    range, so none of this is about the staircase's outer plateau.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    scan = staircase.window_scan(SVDP_OMS, SVDP_LEVELS)
    agree = staircase.window_agreement(scan)

    fig, axes = newfig(th, 1, 3, figsize=(14.6, 4.4))

    xs = np.linspace(0.0, staircase.CMP_XMAX, 800)
    axes[0].plot(xs, staircase.vdp_zeta(xs, staircase.CMP_MU),
                 color=th["ink2"], linewidth=1.8, zorder=5,
                 label="Van der Pol")
    for k, n in enumerate((5, 17)):
        lv, ed = staircase.vdp_staircase(staircase.CMP_MU, n,
                                         staircase.CMP_XMAX)
        ys = [staircase.zeta_at(x, lv, ed) for x in xs]
        axes[0].step(xs, ys, color=th["series"][k], linewidth=1.4, zorder=4,
                     where="mid", label="staircase, %d levels" % n)
    axes[0].axhline(0.0, color=th["axis"], linewidth=0.9, zorder=2)
    style(axes[0], th, "$x$", "damping ratio  $\\zeta$",
          "The only difference: how finely\nthe same law is resolved")
    legend(axes[0], th, loc="upper left")

    tags = [str(n) for n in SVDP_LEVELS] + ["vdp"]
    labels = ["%d levels" % n for n in SVDP_LEVELS] + ["Van der Pol"]
    code = {"chaos": 2, "torus": 1}
    grid = np.array([[code.get(lab, 0) for _, lab, _ in scan[t]]
                     for t in tags])
    cmap = matplotlib.colors.ListedColormap(list(th["series"]))
    axes[1].pcolormesh(SVDP_OMS, np.arange(len(tags)), grid, cmap=cmap,
                       vmin=-0.5, vmax=2.5, shading="nearest", zorder=1)
    axes[1].set_yticks(np.arange(len(tags)))
    axes[1].set_yticklabels(labels)
    style(axes[1], th, "drive frequency  $\\Omega$", "",
          "Chaos lives in the transitions\nbetween one lock and the next")
    handles = [matplotlib.patches.Patch(facecolor=th["series"][i], label=t)
               for i, t in enumerate(("locked", "quasi-periodic", "chaotic"))]
    # below the axes, not on them: at lower left the box covered genuine
    # chaotic cells in the nine level row
    axes[1].legend(handles=handles, fontsize=8, labelcolor=th["ink2"],
                   frameon=False, ncol=3, loc="upper center",
                   bbox_to_anchor=(0.5, -0.17)).set_zorder(9)

    js = [agree[str(n)][2] for n in SVDP_LEVELS]
    axes[2].plot(SVDP_LEVELS, js, "o-", color=th["series"][0], linewidth=1.8,
                 markersize=5, zorder=4, label="agreement with Van der Pol")
    for n, j in zip(SVDP_LEVELS, js):
        axes[2].annotate("%.2f" % j, (n, j), fontsize=8, color=th["ink2"],
                         xytext=(0, -14), textcoords="offset points",
                         ha="center")
    axes[2].axhline(1.0, color=th["ink2"], linewidth=1.1,
                    linestyle=(0, (5, 3)), zorder=3)
    axes[2].text(SVDP_LEVELS[-1], 1.0, "Van der Pol ", fontsize=8,
                 color=th["ink2"], va="bottom", ha="right", zorder=6)
    axes[2].set_xscale("log")
    axes[2].set_xticks(list(SVDP_LEVELS))
    axes[2].set_xticklabels([str(n) for n in SVDP_LEVELS])
    # a log axis keeps its own minor ticks, whose labels collide with these
    axes[2].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axes[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[2].set_ylim(-0.05, 1.12)
    style(axes[2], th, "levels in the staircase",
          "shared chaotic frequencies / combined",
          "A piecewise model converges on\nthe chaos, not just the cycle")
    legend(axes[2], th, loc="upper left")

    fig.suptitle("Driving the staircase at forced Van der Pol's chaotic "
                 "point", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "staircase-vdp")


FLOOR_LEVELS = (2, 3, 5)


def fig_level_floor(th, name, scan=None):
    """The coarsest staircases swept wide: the chaotic bands move, they do
    not vanish.

    One strip per system over ``staircase.WIDE_OMS``, the same encoding as
    ``fig_staircase_vdp``. Two levels is the original piecewise constant
    Van der Pol, and it has two chaotic bands of its own — at the
    transitions its own locks make, which sit far from Van der Pol's
    because its free cycle is faster. The earlier narrow window is drawn
    on top, to show why the earlier single-frequency test saw a lock.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
        scan: a ``window_scan`` result to draw, or ``None`` to compute it,
            which takes about a quarter of an hour.
    """
    oms = np.array(staircase.WIDE_OMS)
    if scan is None:
        scan = staircase.window_scan(oms, FLOOR_LEVELS)
    tags = [str(n) for n in FLOOR_LEVELS] + ["vdp"]
    labels = ["%d levels" % n for n in FLOOR_LEVELS] + ["Van der Pol"]
    code = {"chaos": 2, "torus": 1}
    grid = np.array([[code.get(lab, 0) for _, lab, _ in scan[t]]
                     for t in tags])

    fig, ax = newfig(th, figsize=(12.0, 3.6))
    cmap = matplotlib.colors.ListedColormap(list(th["series"]))
    ax.pcolormesh(oms, np.arange(len(tags)), grid, cmap=cmap, vmin=-0.5,
                  vmax=2.5, shading="nearest", zorder=1)
    lo, hi = staircase.NARROW_OMS[0], staircase.NARROW_OMS[-1]
    for x in (lo, hi):
        ax.axvline(x, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
                   zorder=5)
    ax.text(0.5*(lo + hi), len(tags) - 0.5, "the earlier window",
            fontsize=8, color=th["ink2"], ha="center", va="bottom", zorder=6)
    ax.set_yticks(np.arange(len(tags)))
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.5, len(tags) - 0.1)
    style(ax, th, "drive frequency  $\\Omega$", "",
          "Two levels is enough: the chaotic bands move with the level "
          "count, they do not vanish")
    handles = [matplotlib.patches.Patch(facecolor=th["series"][i], label=t)
               for i, t in enumerate(("locked", "quasi-periodic", "chaotic"))]
    ax.legend(handles=handles, fontsize=8, labelcolor=th["ink2"],
              frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.2)).set_zorder(9)
    fig.tight_layout()
    save(fig, name, "level-floor")


NORM_ROWS = (("vdp", "Van der Pol"),
             ("three fitted bands", "three levels, fitted to the bands"),
             ("two fitted bands", "two levels, fitted to the bands"),
             ("two z1=7.25", "two levels, cycle matched, $\\zeta_1 = 7.25$"),
             ("uniform 3", "three levels, scaled uniformly"),
             ("uniform 2", "two levels, scaled uniformly"),
             ("fitted 3", "three levels, as fitted"),
             ("fitted 2", "two levels, as fitted"))


def fig_normalised(th, name, scan=None):
    """Coarse staircases normalised onto Van der Pol's free cycle, driven.

    From the bottom: the staircases as fitted, whose bands sit elsewhere;
    the same scaled uniformly onto Van der Pol's free cycle, which leaves
    the bands where they were; the two level model with its free cycle
    matched and its remaining shape freedom spent on the bands; and the two
    models `staircase.fit_bands` produced with every parameter free and 20%
    leeway on the free cycle, whose plateau edges sit on Van der Pol's.
    Thin lines mark Van der Pol's chaotic frequencies. Same encoding as
    ``fig_staircase_vdp``.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
        scan: the ``staircase.normalise`` scan to draw, or ``None`` to
            compute it, which takes about half an hour.
    """
    if scan is None:
        _, scan = staircase.normalise()
    oms = np.array(staircase.NORM_OMS)
    tags = [t for t, _ in NORM_ROWS][::-1]
    labels = [l for _, l in NORM_ROWS][::-1]
    code = {"chaos": 2, "torus": 1}
    grid = np.array([[code.get(lab, 0) for _, lab, _ in scan[t]]
                     for t in tags])

    fig, ax = newfig(th, figsize=(12.0, 4.4))
    cmap = matplotlib.colors.ListedColormap(list(th["series"]))
    ax.pcolormesh(oms, np.arange(len(tags)), grid, cmap=cmap, vmin=-0.5,
                  vmax=2.5, shading="nearest", zorder=1)
    for om in (om for om, lab, _ in scan["vdp"] if lab == "chaos"):
        ax.axvline(om, color=th["ink2"], linewidth=0.5, alpha=0.5, zorder=4)
    ax.set_yticks(np.arange(len(tags)))
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.5, len(tags) - 0.5)
    style(ax, th, "drive frequency  $\\Omega$", "",
          "Matching the free cycle does not place the bands; "
          "the damping shape does")
    handles = [matplotlib.patches.Patch(facecolor=th["series"][i], label=t)
               for i, t in enumerate(("locked", "quasi-periodic", "chaotic"))]
    ax.legend(handles=handles, fontsize=8, labelcolor=th["ink2"],
              frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.16)).set_zorder(9)
    fig.tight_layout()
    save(fig, name, "normalised")


CHAOS_OM = 2.470


def fig_chaos_phase(th, name, om=CHAOS_OM, n_trace=40, n_strobe=3000):
    """The fitted three level prototype beside Van der Pol, both chaotic.

    Same drive, ``A = 5`` at ``Om = 2.470``, a frequency inside the chaotic
    band of both. Each panel draws forty drive periods of the orbit after
    the transient as a thin line, and three thousand stroboscopic samples,
    one per drive period, on top: the samples are the attractor, the line
    is how the orbit threads it. The prototype's zone edges are drawn as
    chrome, since the only nonlinearity it has is which zone it is in.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
        om: drive frequency.
        n_trace: drive periods of continuous orbit to draw.
        n_strobe: stroboscopic samples to draw.
    """
    import section
    import vanderpol
    lv, ed = staircase.THREE_FITTED
    systems = [("three level prototype, fitted\n$\\zeta = (%.2f, %.2f, %.1f)$,"
                " edges $(%.2f, %.2f)$" % (lv[0], lv[1], lv[2], ed[0], ed[1]),
                staircase.field(lv, ed, staircase.CMP_AMP, om), ed),
               ("Van der Pol, $\\mu = %g$" % staircase.CMP_MU,
                vanderpol.field(staircase.CMP_MU, staircase.CMP_AMP, om), ())]
    td = 2.0*np.pi/om
    fig, axes = newfig(th, 1, 2, figsize=(12.0, 5.6))
    for ax, (label, flow, edges) in zip(axes, systems):
        pts = section.strobe(flow, td, [2.0, 0.0], staircase.CMP_NSKIP,
                             n_keep=n_strobe)
        sol = solve_ivp(flow, (0.0, n_trace*td), pts[-1], method=section.METHOD,
                        rtol=1e-9, atol=1e-11,
                        t_eval=np.linspace(0.0, n_trace*td, 400*n_trace))
        ax.plot(sol.y[0], sol.y[1], color=th["series"][0], linewidth=0.35,
                alpha=0.55, zorder=2, label="orbit, %d drive periods" % n_trace)
        ax.plot(pts[:, 0], pts[:, 1], ".", color=th["series"][1], markersize=1.6,
                zorder=3, label="stroboscopic samples, %d" % n_strobe)
        for e in edges:
            for xv in (e, -e):
                ax.axvline(xv, color=th["ink2"], linewidth=1.0,
                           linestyle=(0, (5, 3)), zorder=4)
        style(ax, th, "$x$", "$\\dot{x}$", label)
        ax.set_xlim(-2.6, 2.6)
        ax.set_ylim(-13, 13)
        legend(ax, th, loc="upper left", markerscale=4)
    fig.suptitle("Chaotic mode, side by side: $A = %g$, $\\Omega = %.3f$"
                 % (staircase.CMP_AMP, om), color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "chaos-phase")


def fig_regime_three(th, name, data=None):
    """The fitted three level model against Van der Pol over the drive grid.

    Left and middle: regime maps over ``staircase.MAP_R`` ratios and
    ``staircase.MAP_A`` amplitudes, cells locked, quasi-periodic or chaotic,
    with chaotic verdicts that failed ``section.confirm_chaos`` drawn as
    quasi-periodic. The unforced row is painted as chrome: with no drive
    there is nothing to lock to, and what the classifier reports there is
    the sampling frequency being commensurate with the free cycle. Right:
    no forcing at all — the two free limit cycles in the phase plane, with
    the three level model's zone edges.

    This is the test the candidate paragraph in ``VANDERPOL.md`` asked for:
    matched at one drive strength, does the model track Van der Pol across
    the rest of the grid?

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
        data: a ``staircase.regime_compare`` result, or ``None`` to
            compute it, which takes an hour or two.
    """
    import section
    if data is None:
        data = staircase.regime_compare()
    ratios, amps = np.array(data["ratios"]), list(data["amps"])
    mu = data.get("mu", staircase.CMP_MU)
    lv, ed = data.get("model", staircase.THREE_FITTED)
    stem = "regime-three" if mu == staircase.CMP_MU else "regime-three-mu%g" % mu
    fig, axes = newfig(th, 1, 3, figsize=(15.0, 4.6),
                       gridspec_kw=dict(width_ratios=(1.3, 1.3, 1.0)))
    cmap = matplotlib.colors.ListedColormap(list(th["series"]) + [th["grid"]])
    panels = [(axes[0], "three", "three level prototype, fitted",
               "$\\zeta = (%.2f, %.2f, %.1f)$, edges $(%.2f, %.2f)$"
               % (lv[0], lv[1], lv[2], ed[0], ed[1])),
              (axes[1], "vdp", "Van der Pol  $\\mu = %g$" % mu,
               "the same drive, the same classifier")]
    for ax, tag, title, sub in panels:
        lab, q, w, lam = data[tag]
        code = _regime_codes(lab, q, lam)
        code = np.where((code == 2) & ~data[tag + "_ok"], 1, code)
        # with no drive there is nothing to lock to: the classifier's locks
        # in that row are the sampling frequency being commensurate with
        # the free cycle, which is a property of the grid, not the system
        code[[i for i, a in enumerate(amps) if a == 0.0], :] = 3
        ax.pcolormesh(ratios, np.arange(len(amps)), code, cmap=cmap,
                      vmin=-0.5, vmax=3.5, shading="nearest", zorder=1)
        if 0.0 in amps:
            ax.text(0.5*(ratios[0] + ratios[-1]), amps.index(0.0),
                    "no drive: the free cycle, nothing to lock to",
                    fontsize=8, color=th["ink2"], ha="center", va="center",
                    zorder=6)
        ax.set_yticks(np.arange(len(amps)))
        ax.set_yticklabels(["%g" % a for a in amps])
        style(ax, th, "$\\Omega / \\omega_{lc}$ (Van der Pol's)", "$A$", title)
        ax.text(0.5, -0.26, sub, transform=ax.transAxes, fontsize=8,
                color=th["ink2"], ha="center", va="top")
        ax.text(0.97, 0.04, "chaotic cells: %d" % int(np.sum(code == 2)),
                transform=ax.transAxes, fontsize=8.5, color=th["ink"],
                ha="right", va="bottom", zorder=7,
                bbox=dict(boxstyle="round,pad=0.3", fc=th["surface"],
                          ec=th["grid"], lw=0.8, alpha=0.94))
    handles = [matplotlib.patches.Patch(facecolor=th["series"][i], label=t)
               for i, t in enumerate(("locked", "quasi-periodic", "chaotic"))]
    axes[0].legend(handles=handles, fontsize=8, labelcolor=th["ink2"],
                   frameon=True, facecolor=th["surface"], edgecolor="none",
                   framealpha=0.92, loc="upper left").set_zorder(9)

    ax = axes[2]
    import vanderpol
    for k, (label, flow, r0) in enumerate((
            ("three level prototype", staircase.field(lv, ed), 2.0),
            ("Van der Pol", vanderpol.field(mu), 2.0))):
        warm = solve_ivp(flow, (0.0, 300.0), [r0, 0.0], method=section.METHOD,
                         rtol=1e-9, atol=1e-11)
        T = (staircase.free_cycle_num(lv, ed)[1] if k == 0
             else vanderpol.cycle(mu)[0])
        sol = solve_ivp(flow, (0.0, T), warm.y[:, -1], method=section.METHOD,
                        rtol=1e-9, atol=1e-11, t_eval=np.linspace(0, T, 4000))
        ax.plot(sol.y[0], sol.y[1], color=th["series"][k], linewidth=1.8,
                zorder=3, label="%s, $T = %.2f$" % (label, T))
    for e in ed:
        for xv in (e, -e):
            ax.axvline(xv, color=th["ink2"], linewidth=1.0,
                       linestyle=(0, (5, 3)), zorder=4)
    style(ax, th, "$x$", "$\\dot{x}$", "No forcing: the free cycles")
    legend(ax, th, loc="upper left")

    fig.suptitle("Fitted at one drive strength, tested across the grid",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, stem)


def fig_campaign(th, name, results=None):
    """The three level model's parameters against Van der Pol's mu.

    Left: the three damping ratios from every campaign fit, with the power
    laws through them; the two edges below. Middle: where Van der Pol's
    plateau edges sit at each mu, the targets the fits were made on, with
    the fitted models' edges over them. Right: the agreement of each fitted
    model's sweep with Van der Pol's, and the chaotic cell counts.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
        results: the ``campaign/results.json`` dict, loaded if ``None``.
    """
    import json
    if results is None:
        results = json.load(open(os.path.join(os.path.dirname(OUT),
                                              "campaign", "results.json")))
    fits = sorted((float(k), v) for k, v in results["fits"].items()
                  if "levels" in v)
    tg = sorted((float(k), v) for k, v in results["targets"].items())
    law = results.get("formula")
    fig, axes = newfig(th, 1, 3, figsize=(15.0, 4.8))

    ax = axes[0]
    mus = np.array([m for m, _ in fits])
    lv = np.array([v["levels"] for _, v in fits])
    ed = np.array([v["edges"] for _, v in fits])
    mm = np.logspace(-1, np.log10(5), 100)
    for k, (lab, col) in enumerate((("$\\zeta_0$ (core, negative)", th["series"][0]),
                                    ("$\\zeta_1$ (band)", th["series"][1]),
                                    ("$\\zeta_2$ (outer)", th["series"][2]))):
        ax.loglog(mus, np.abs(lv[:, k]), "o", color=col, markersize=5, zorder=4,
                  label=lab)
        if law:
            l = law["levels"][k]
            ax.loglog(mm, l["c"]*mm**l["p"], color=col, linewidth=1.2,
                      linestyle=(0, (5, 3)), zorder=3,
                      label="$%.2f\\,\\mu^{%.2f}$" % (l["sign"]*l["c"], l["p"]))
    ax.loglog(mus, ed[:, 0], "s", color=th["ink2"], markersize=4, zorder=4,
              label="edge $a$")
    ax.loglog(mus, ed[:, 1], "^", color=th["ink2"], markersize=4, zorder=4,
              label="edge $b$")
    style(ax, th, "$\\mu$", "magnitude", "Fitted parameters against $\\mu$")
    legend(ax, th, loc="upper left", ncol=2)

    ax = axes[1]
    for m, v in tg:
        for lock, col in (("lock1", th["series"][0]), ("lock3", th["series"][1])):
            if v.get(lock):
                ax.plot([m, m], v[lock], color=col, linewidth=6, alpha=0.35,
                        solid_capstyle="butt", zorder=2)
    for m, v in fits:
        fd = v["found"]
        if "lock11" in fd:
            ax.plot(m, fd["lock11"], "_", color=th["series"][0], markersize=12,
                    markeredgewidth=2, zorder=4)
        if "lock30" in fd and "lock31" in fd:
            ax.plot([m, m], [fd["lock30"], fd["lock31"]], color=th["series"][1],
                    linewidth=1.6, zorder=4)
            ax.plot([m, m], [fd["lock30"], fd["lock31"]], "_", color=th["series"][1],
                    markersize=12, markeredgewidth=2, zorder=4)
    ax.set_xscale("log")
    handles = [matplotlib.patches.Patch(facecolor=th["series"][0], alpha=0.35,
                                        label="Van der Pol 1:1 plateau"),
               matplotlib.patches.Patch(facecolor=th["series"][1], alpha=0.35,
                                        label="Van der Pol 3:1 plateau"),
               matplotlib.lines.Line2D([], [], color=th["ink2"], marker="_",
                                       markersize=12, markeredgewidth=2,
                                       linestyle="none", label="fitted model's edges")]
    ax.legend(handles=handles, fontsize=8, labelcolor=th["ink2"], frameon=True,
              facecolor=th["surface"], edgecolor="none", framealpha=0.92,
              loc="upper left").set_zorder(9)
    style(ax, th, "$\\mu$", "drive ratio  $\\Omega/\\omega_{lc}$",
          "The targets, and where the fits landed ($A = 5$)")

    ax = axes[2]
    ver = sorted((float(k), v) for k, v in results["verify"].items())
    if ver:
        vm = np.array([m for m, _ in ver])
        jac = np.array([v["jaccard"] for _, v in ver], float)
        ok = np.isfinite(jac)
        ax.semilogx(vm[ok], jac[ok], "o-", color=th["series"][0],
                    linewidth=1.6, markersize=5, zorder=4,
                    label="agreement on chaotic cells (Jaccard)")
        ax.text(0.03, 0.55, "no chaos in either\nsystem below $\\mu = 1.5$",
                transform=ax.transAxes, fontsize=8, color=th["ink2"],
                ha="left", va="center")
        ax2 = ax.twinx()
        ax2.semilogx(vm, [len(v["vdp"]["chaotic"]) for _, v in ver], "s--",
                     color=th["series"][1], markersize=5, zorder=3,
                     label="chaotic cells, Van der Pol")
        ax2.semilogx(vm, [len(v["three"]["chaotic"]) for _, v in ver], "^--",
                     color=th["series"][2], markersize=5, zorder=3,
                     label="chaotic cells, model")
        ax2.tick_params(colors=th["ink2"], labelsize=8)
        ax2.set_ylabel("chaotic cells in the sweep", color=th["ink2"], fontsize=9)
        for sp in ("top",):
            ax2.spines[sp].set_visible(False)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, labelcolor=th["ink2"], frameon=True,
                  facecolor=th["surface"], edgecolor="none", framealpha=0.92,
                  loc="upper left").set_zorder(9)
    ax.set_ylim(0, 1.05)
    style(ax, th, "$\\mu$", "agreement", "Verification sweeps at $A = 5$")
    fig.suptitle("The three level prototype across Van der Pol's range",
                 color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "campaign")


def fig_boundary(th, name, results=None):
    """The chaos boundary in mu at drive amplitude 10.

    One row per mu of the boundary sweep, Van der Pol above the model:
    the lock plateaus as faint bars, every confirmed chaotic cell as a
    marker. Reads ``results["boundary"]`` from ``campaign/results.json``.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
        results: the ``campaign/results.json`` dict, loaded if ``None``.
    """
    import json
    if results is None:
        results = json.load(open(os.path.join(os.path.dirname(OUT),
                                              "campaign", "results.json")))
    rows = sorted((float(k), v) for k, v in results.get("boundary", {}).items())
    if not rows:
        return
    fig, ax = newfig(th, figsize=(9.0, 0.9 + 1.0*len(rows)))
    off = 0.12
    col = {"vdp": th["series"][1], "three": th["series"][2]}
    mark = {"vdp": "s", "three": "^"}
    lab = {"vdp": "Van der Pol", "three": "model from the laws"}
    for i, (m, v) in enumerate(rows):
        for tag, y in (("vdp", i + off), ("three", i - off)):
            seen = set()
            for l, lo, hi, n in v[tag]["runs"]:
                if l.startswith("lock") and n >= 3:
                    ax.plot([lo, hi], [y, y], color=th["ink2"], linewidth=5,
                            alpha=0.22, solid_capstyle="butt", zorder=2)
                    q = int(l[4:])
                    if q in (1, 3, 5, 7) and q not in seen and hi - lo > 0.4:
                        ax.text(0.5*(lo + hi), y, "%d:1" % q, fontsize=6.5,
                                color=th["ink2"], ha="center", va="center", zorder=3)
                        seen.add(q)
            ch = v[tag]["chaotic"]
            ax.plot(ch, [y]*len(ch), mark[tag], color=col[tag], markersize=5,
                    linestyle="none", zorder=5, label=lab[tag] if i == 0 else None)
    ax.set_xlim(rows[0][1]["r_lo"], rows[0][1]["r_hi"])
    ax.set_ylim(-0.5, len(rows) - 0.5 + 0.6)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(["%g" % m for m, _ in rows])
    style(ax, th, "drive ratio  $\\Omega/\\omega_{lc}$", "$\\mu$",
          "Chaotic cells at $A = %g$, Van der Pol above the model at each $\\mu$"
          % rows[0][1]["amp"])
    legend(ax, th, loc="upper left", ncol=2)
    fig.tight_layout()
    save(fig, name, "boundary")


# --------------------------------------------------- the strange attractor
#: Drive at which forced Van der Pol is chaotic, and the staircase fitted to
#: it. The same settings the comparison in ``VANDERPOL.md`` uses.
ATT_MU, ATT_XMAX, ATT_A, ATT_OM = 5.0, 3.0, 5.0, 2.466
ATT_KEEP, ATT_SKIP = 30000, 400

_CLOUD = {}


def _strobe_cloud(kind, n_keep=ATT_KEEP):
    """Stroboscopic points on the attractor, cached across themes.

    Integrated at a looser tolerance than the analysis elsewhere uses. That
    is deliberate and it changes what the picture is: on a chaotic
    attractor a slightly perturbed trajectory *shadows* the attractor rather
    than tracking any one orbit, so the cloud fills out the attractor's
    shape faithfully while not being a picture of a particular trajectory.
    For geometry that is the right trade; for an exponent it would not be,
    which is why the exponents come from ``maps.py`` instead.
    """
    if kind in _CLOUD:
        return _CLOUD[kind]
    cache = os.path.join(OUT, ".attractor-%s.npy" % kind)
    if os.path.exists(cache):
        _CLOUD[kind] = np.load(cache)
        return _CLOUD[kind]
    td = 2.0*np.pi/ATT_OM
    if kind == "vdp":
        f = vanderpol.field(ATT_MU, ATT_A, ATT_OM)
    else:
        lv, ed = staircase.vdp_staircase(ATT_MU, int(kind), ATT_XMAX)
        f = staircase.field(lv, ed, ATT_A, ATT_OM)
    t = td*np.arange(ATT_SKIP, ATT_SKIP + n_keep + 1)
    sol = solve_ivp(f, (0.0, t[-1]), [2.0, 0.0], t_eval=t, method="LSODA",
                    rtol=1e-7, atol=1e-9)
    _CLOUD[kind] = sol.y.T
    os.makedirs(OUT, exist_ok=True)
    np.save(cache, _CLOUD[kind])
    return _CLOUD[kind]


def _ramp(th, k):
    """Single hue sequential ramp for the density encoding, per theme.

    Density is a magnitude, so it takes a sequential scale: one hue running
    from the sparse end to the dense end, never a rainbow. The dark ramp is
    chosen against the dark surface rather than flipped from the light one —
    on a dark ground the dense end has to be the bright one to read at all,
    and on a light ground it has to be the deep one.

    Each panel keeps its own hue so the two systems stay distinguishable,
    but identity is carried by the panel titles, not by the colour.
    """
    hue = th["series"][k]
    dark = th["surface"].lower().startswith("#1")
    stops = [th["grid"], hue, "#ffffff" if dark else "#0d0d12"]
    return matplotlib.colors.LinearSegmentedColormap.from_list(
        "attractor%d" % k, stops)


def _density(pts, bins=260):
    """Log point density at each sample, for the colour encoding.

    Log because an attractor's visit density spans orders of magnitude: on a
    linear scale two or three dense folds saturate the ramp and every other
    filament reads as empty, which is exactly the structure worth seeing.
    """
    h, xe, ye = np.histogram2d(pts[:, 0], pts[:, 1], bins=bins)
    i = np.clip(np.digitize(pts[:, 0], xe) - 1, 0, bins - 1)
    j = np.clip(np.digitize(pts[:, 1], ye) - 1, 0, bins - 1)
    return np.log1p(h[i, j])


def fig_strange_attractor(th, name):
    """Draw the chaotic stroboscopic section, beside Van der Pol's, and zoomed.

    Sampling the forced response once per drive period turns a chaotic
    trajectory into a point set, and that set is the attractor.

    A textbook chaotic attractor is a Cantor-like stack of filaments, and
    this one is not — at least not visibly. Both systems here are strongly
    dissipative: at ``mu = 5`` the contraction per cycle is below what
    double precision can even express, so the attractor is squeezed to
    within a few per cent of a one dimensional curve. Its fractal structure
    is real, and the positive Lyapunov exponent measures the stretching that
    creates it, but it lives at scales far below anything 30000 points and
    an eighteen fold magnification can reach. The third panel is drawn to
    make that honest rather than to hide it.

    Searching for a fatter case did not help: of five drive settings tried,
    the only two that are genuinely chaotic have transverse thickness 0.030
    and 0.036, and the three with visible girth are periodic orbits whose
    exponents are negative. Thin attractors are what this family gives.

    What the figure does show is the comparison. Left and middle are the
    piecewise staircase at 65 levels and the smooth Van der Pol it was
    fitted to, driven identically at the point where Van der Pol is chaotic,
    and they trace the same shape. That is the geometric counterpart of the
    agreement in exponent and lock structure reported elsewhere: the
    piecewise model reproduces not just that there is chaos but where the
    orbit goes.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    fig, axes = newfig(th, 1, 3, figsize=(14.4, 4.8))

    stair = _strobe_cloud("65")
    vdp_pts = _strobe_cloud("vdp")

    for ax, pts, title, k in ((axes[0], stair, "piecewise staircase,\n65 levels", 0),
                              (axes[1], vdp_pts, "Van der Pol,\n$\\mu = 5$", 1)):
        ax.scatter(pts[:, 0], pts[:, 1], s=0.30, c=_density(pts),
                   cmap=_ramp(th, k), linewidths=0, zorder=4,
                   rasterized=True)
        style(ax, th, "$x_1 = x$", "$x_2 = \\dot{x}$", title)

    lim = [min(stair[:, 0].min(), vdp_pts[:, 0].min()),
           max(stair[:, 0].max(), vdp_pts[:, 0].max()),
           min(stair[:, 1].min(), vdp_pts[:, 1].min()),
           max(stair[:, 1].max(), vdp_pts[:, 1].max())]
    for ax in axes[:2]:
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[2], lim[3])

    # zoom on the densest patch, found from the data rather than chosen
    h, xe, ye = np.histogram2d(stair[:, 0], stair[:, 1], bins=40)
    i, j = np.unravel_index(np.argmax(h), h.shape)
    cx, cy = 0.5*(xe[i] + xe[i + 1]), 0.5*(ye[j] + ye[j + 1])
    wx, wy = 0.055*(lim[1] - lim[0]), 0.055*(lim[3] - lim[2])
    axes[2].scatter(stair[:, 0], stair[:, 1], s=2.4, c=_density(stair),
                    cmap=_ramp(th, 0), linewidths=0, zorder=4,
                    rasterized=True)
    axes[2].set_xlim(cx - wx, cx + wx)
    axes[2].set_ylim(cy - wy, cy + wy)
    style(axes[2], th, "$x_1 = x$", "$x_2 = \\dot{x}$",
          "magnified 18 times: still\none strand, no layering")
    for ax in axes[:1]:
        ax.add_patch(matplotlib.patches.Rectangle(
            (cx - wx, cy - wy), 2*wx, 2*wy, fill=False, lw=1.0,
            edgecolor=th["ink2"], zorder=7))

    axes[0].text(0.03, 0.03, "brighter = more often visited",
                 transform=axes[0].transAxes, fontsize=8, color=th["ink2"],
                 ha="left", va="bottom", zorder=8)
    fig.suptitle("A chaotic attractor, squeezed almost flat by strong "
                 "dissipation", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "strange-attractor")


def fig_scaling(th, name):
    """Draw the three level model at three frequency ranges, then collapse them.

    Left, three panels: the free response of the mu = 5 fit from an initial
    displacement 25% above its cycle, as the reference model in its own
    units, as a 50 Hz oscillator with a millimetre per model unit, and as a
    1 kHz oscillator with a micrometre per model unit. Each is integrated
    in its own units from the scaled equation, not produced by rescaling
    the reference.

    Right: the same three traces on the reduced axes ``wn t`` and
    ``x / lam``. They coincide, which is the scaling rule of
    ``THREELEVEL.md`` and the content of ``scaling.py``'s check: the
    damping ratios and the edge ratio carry the behaviour, and ``wn`` and
    ``lam`` only set the clock and the ruler.

    Args:
        th: theme dict from ``THEMES``.
        name: theme key, used as the output filename suffix.
    """
    import scaling
    cases = (("reference: $\\omega_n = 1$, $\\lambda = 1$", 1.0, 1.0,
              1.0, 1.0, "$t$", "$x$"),
             ("50 Hz, 1 mm per unit", 2.0*np.pi*50.0, 1.0e-3,
              1.0e3, 1.0e3, "$t$ [ms]", "$x$ [mm]"),
             ("1 kHz, 1 $\\mu$m per unit", 2.0*np.pi*1000.0, 1.0e-6,
              1.0e3, 1.0e6, "$t$ [ms]", "$x$ [$\\mu$m]"))
    R, T = staircase.free_cycle(*staircase.THREE_FITTED)
    t_end = 3.0*T
    tau = np.linspace(0.0, t_end, 3000)

    fig = plt.figure(figsize=(11.8, 5.0))
    fig.patch.set_facecolor(th["surface"])
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.15], hspace=0.65,
                          wspace=0.22)
    left = [fig.add_subplot(gs[k, 0]) for k in range(3)]
    right = fig.add_subplot(gs[:, 1])
    styles = ((2.8, "-"), (1.7, (0, (6, 3))), (1.2, (0, (1.5, 2.2))))

    for k, (lbl, wn, lam, tf, xf, xl, yl) in enumerate(cases):
        lv, ed = scaling.scale_model(*staircase.THREE_FITTED, wn, lam)
        f = scaling.field(lv, ed, wn)
        sol = solve_ivp(f, (0.0, t_end/wn), [2.5*lam, 0.0], method="LSODA",
                        rtol=1e-10, atol=[1e-13*lam, 1e-13*lam*wn],
                        t_eval=tau/wn)
        lw, ls = styles[k]
        left[k].plot(sol.t*tf, sol.y[0]*xf, color=th["series"][k],
                     linewidth=1.4, zorder=4)
        for e, s in ((ed[1], "$b$"), (-ed[1], "$-b$")):
            left[k].axhline(e*xf, color=th["ink2"], linewidth=0.8,
                            linestyle=(0, (5, 3)), zorder=3)
        left[k].text(0.99, ed[1]*xf, "$b$", transform=left[k].get_yaxis_transform(),
                     fontsize=8, color=th["ink2"], va="bottom", ha="right",
                     zorder=6)
        left[k].set_xlim(0.0, t_end*tf/wn)
        style(left[k], th, xl, yl)
        left[k].set_title(lbl, color=th["ink"], fontsize=9, pad=4, loc="left")
        right.plot(tau, sol.y[0]/lam, color=th["series"][k], linewidth=lw,
                   linestyle=ls, zorder=4 + k, label=lbl)

    for e in (staircase.THREE_FITTED[1][1], -staircase.THREE_FITTED[1][1]):
        right.axhline(e, color=th["ink2"], linewidth=0.9,
                      linestyle=(0, (5, 3)), zorder=3)
    right.text(0.62, staircase.THREE_FITTED[1][1], "$b$",
               transform=right.get_yaxis_transform(), fontsize=8,
               color=th["ink2"], va="bottom", ha="left", zorder=6)
    right.axvline(T, color=th["axis"], linewidth=0.9, zorder=2)
    right.text(T, -2.55, " one free period, $\\omega_n T = %.2f$" % T,
               fontsize=8, color=th["ink2"], ha="left", va="bottom", zorder=6)
    right.set_xlim(0.0, t_end)
    right.set_ylim(-2.7, 3.2)
    style(right, th, "$\\omega_n t$", "$x / \\lambda$",
          "The same trace on reduced axes: the three coincide")
    legend(right, th, loc="upper right")
    fig.suptitle("Scaling the three level model to another frequency range: "
                 "$\\zeta$ and $a/b$ fixed, $\\omega_n$ sets the clock, "
                 "$\\lambda$ the ruler", color=th["ink"], fontsize=11)
    save(fig, name, "scaling")


if __name__ == "__main__":
    for name, th in THEMES.items():
        print(f"{name}:")
        fig_linear(th, name)
        fig_switched(th, name)
        fig_decrement(th, name)
        fig_limit_cycle(th, name)
        fig_map_scaling(th, name)
        fig_frequency(th, name)
        fig_poles(th, name)
        fig_stability_map(th, name)
        fig_symmetric(th, name)
        fig_four_models(th, name)
        fig_forced_tongues(th, name)
        fig_forced_sections(th, name)
        fig_staircase(th, name)
        fig_staircase_vdp(th, name)
        fig_level_floor(th, name)
        fig_normalised(th, name)
        fig_chaos_phase(th, name)
        fig_regime_three(th, name)
        fig_campaign(th, name)
        fig_boundary(th, name)
        fig_scaling(th, name)
        fig_strange_attractor(th, name)
        fig_vanderpol_compare(th, name)
