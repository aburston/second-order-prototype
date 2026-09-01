"""Generate the phase plane figures used in README.md.

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
