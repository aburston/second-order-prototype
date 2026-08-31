"""Generate the phase plane figures used in README.md.

Run:  python3 figures.py          (writes figures/*.png, light and dark)

Every figure is rendered twice, once per theme, and embedded in the README
through a <picture> element so GitHub picks the right one.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

WN = 1.0
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# Validated palette: categorical slots 1-3, chart chrome and ink.
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e",
                  grid="#e1e0d9", axis="#c3c2b7",
                  series=("#2a78d6", "#eb6834", "#1baf7a")),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7",
                  grid="#2c2c2a", axis="#383835",
                  series=("#3987e5", "#d95926", "#199e70")),
}


# ---------------------------------------------------------------- dynamics
def f_linear(zeta):
    def f(t, y):
        return [y[1], -WN**2*y[0] - 2*zeta*WN*y[1]]
    return f


def f_switched(zp, zm, v0=0.0):
    """Switched damping acting on the velocity relative to the boundary."""
    def f(t, y):
        w = y[1] - v0
        return [y[1], -WN**2*y[0] - 2*(zp if w > 0 else zm)*WN*w]
    return f


def traj(f, y0, T, n=4000):
    s = solve_ivp(f, (0, T), y0, t_eval=np.linspace(0, T, n),
                  rtol=1e-11, atol=1e-13)
    return s.y[0], s.y[1], s.t


def xeq(zm, v0):
    """Equilibrium of the offset system (u = 0)."""
    return 2*zm*v0/WN


def preturn(zp, zm, v0, r):
    """One return to the section {x2 = 0, x1 > xeq}."""
    xe = xeq(zm, v0)

    def ev(t, y):
        return y[1]
    ev.direction = -1
    s = solve_ivp(f_switched(zp, zm, v0), (0, 30), [xe + r, 0.0], events=ev,
                  rtol=1e-12, atol=1e-14)
    i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
    return (s.y_events[0][i[0]][0] - xe) if i else np.nan


def rstar(zp, zm, v0, r=None, n=400, tol=1e-11):
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
    lg = ax.legend(fontsize=8, labelcolor=th["ink2"], frameon=True,
                   facecolor=th["surface"], edgecolor="none", framealpha=0.92,
                   **kw)
    lg.set_zorder(7)
    return lg


def boundary(ax, th, y, label):
    """Draw the switching boundary as chrome, not as a data series."""
    ax.axhline(y, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
               zorder=5)
    ax.text(0.985, y, label, transform=ax.get_yaxis_transform(), fontsize=8,
            color=th["ink2"], va="bottom", ha="right", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))


def save(fig, th_name, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{name}-{th_name}.png")
    fig.savefig(path, dpi=200, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print("  wrote", os.path.relpath(path, os.path.dirname(OUT)))


def newfig(th, *a, **kw):
    fig, ax = plt.subplots(*a, **kw)
    fig.patch.set_facecolor(th["surface"])
    return fig, ax


# ------------------------------------------------------------- the figures
def fig_linear(th, name):
    """Linear prototype: phase portrait and time history for three zetas."""
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
    """Boundary through the equilibrium: the three stability cases."""
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
    """Measured amplitude per cycle against the closed form decrement."""
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
    """Offset boundary: the hyperbolic limit cycle and its basin."""
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
    """Return map, and the exact proportionality of amplitude to offset."""
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


if __name__ == "__main__":
    for name, th in THEMES.items():
        print(f"{name}:")
        fig_linear(th, name)
        fig_switched(th, name)
        fig_decrement(th, name)
        fig_limit_cycle(th, name)
        fig_map_scaling(th, name)
