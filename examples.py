"""Physical systems modelled by the prototypes in README.md.

Run ``python3 examples.py`` to print every number quoted in ``EXAMPLES.md``
and to write the figures it embeds into ``figures/`` as ``example-*.png``,
one light and one dark rendering of each, in the same style as
``figures.py``. ``python3 examples.py checks`` prints the numbers only and
``python3 examples.py figures`` writes the figures only; the full run takes
a while because every table is re-integrated at tight tolerance.

Two systems, each mapped onto the prototypes they fit:

``governor``
    Maxwell's governor (1868), reduced to second order. A true governor
    with a one-sided flyball brake puts the switching boundary through the
    equilibrium; a moderator brake set above synchronous speed on a
    synchronised machine gives the offset boundary and a hunting limit
    cycle; a governor deadband gives the symmetric velocity band. The
    brake force ``F (xdot - V)`` acting on the *excess* speed is Maxwell's
    own form, and is exactly the relative-velocity damping the README needs
    to keep the field continuous.

``oscillator``
    A single-transistor LC oscillator with inductive feedback, the
    transistor's gain collapsing outside its linear range. This is Van der
    Pol's triode with the characteristic replaced by a clip, and it is the
    README's piecewise Van der Pol exactly. A smooth ``tanh`` saturation
    with the same small-signal gain and the same saturation current is
    integrated alongside to measure what the hard switch costs. A stage
    that clips on one side only is the asymmetric displacement model,
    bounded only while the loop gain is below two.

Everything is in the README's normalised units, ``wn = 1`` and the band
half-width equal to one, with the physical parameters mapped onto the two
damping ratios in the docstring of each section. The prototypes are never
altered: every physical system is one of them exactly, and where a real
nonlinearity is smoother than the switch (the transistor's saturation) it
is integrated separately as the reference the prototype is compared with.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import figures
from figures import THEMES, style, legend, save, newfig, traj
import frequency
import symmetric
import displacement
import stability

WN = 1.0


def delta(zeta):
    """Logarithmic decrement per half cycle, ``pi zeta / sqrt(1 - zeta^2)``."""
    return np.pi*zeta/np.sqrt(1 - zeta**2)


# ================================================================ governor
# Maxwell's equation for an engine of inertia M, angle x, driving torque P,
# resistance R, with a brake applying liquid friction F (xdot - V) above the
# set speed V and reducing the driving power by G y through an accumulated
# motion y. Reduced to second order here by taking the accumulation to be
# the speed error itself, ydot = xdot - V (a true governor), and the net
# torque to vary with speed as P - R = T0 + c (xdot - V). With M = G = 1:
#
#     zeta_plus  = (F - c) / 2      braked, xdot > V
#     zeta_minus = -c / 2           unbraked, xdot < V
#
# In the frame moving at V the equilibrium is on the boundary. The same
# brake on a machine synchronised to a grid at speed ws (synchronising
# torque K delta, K = 1) with the brake set at ws + v0 gives the offset
# model, and a deadband governor at ws gives the band model.

def zetas(F, c):
    """Map a brake coefficient and a self-excitation slope to damping ratios.

    Normalised so that ``2 sqrt(MG) = 2 sqrt(MK) = 2``.
    """
    return (F - c)/2.0, -c/2.0


def governor_field(F, c, mode, v0=1.0):
    """Vector field of one governor configuration, in physical coordinates.

    ``mode`` is ``"integral"`` (true governor, frame moving at ``V``),
    ``"offset"`` (synchronised machine, brake set ``v0`` above synchronous)
    or ``"deadband"`` (synchronised machine, no governor action within
    ``v0`` of synchronous). The state is ``[angle error, speed error]`` and
    the equilibrium is the origin in every mode.
    """
    def f(t, y):
        x, v = y
        if mode == "integral":
            brake = F*v if v > 0 else 0.0
        elif mode == "offset":
            brake = F*(v - v0) if v > v0 else 0.0
        else:
            brake = F*symmetric.deadzone(v, v0)
        return [v, -x + c*v - brake]
    return f


def offset_rstar(zp, zm, v0=1.0, r=6.0, n=400, tol=1e-11):
    """Fixed point of the README's return map, giving up once it escapes.

    ``figures.rstar`` iterates to its cap however large the radius grows;
    on a diverging case that means integrating astronomically large orbits
    for minutes. Stopping at a radius of a million costs nothing that is
    reported and turns those cases from minutes into seconds.
    """
    for _ in range(n):
        rn = figures.preturn(zp, zm, v0, r)
        if not np.isfinite(rn) or rn > 1e6:
            return np.nan
        if abs(rn - r) < tol*max(1.0, abs(r)):
            return rn
        r = rn
    return r


def rstar_bounded(zp, zm, v0=1.0, r=6.0, n=400, tol=1e-11):
    """``figures.rstar`` with an escape cap, so a diverging case ends early."""
    for _ in range(n):
        rn = figures.preturn(zp, zm, v0, r)
        if not np.isfinite(rn) or rn > 1e6:
            return np.nan
        if abs(rn - r) < tol*max(1.0, abs(r)):
            return rn
        r = rn
    return r


def cycle_ratio(zp, zm):
    """Amplitude ratio over one full cycle by integration, and the period."""
    def f(t, y):
        return [y[1], -y[0] - 2*(zp if y[1] > 0 else zm)*y[1]]

    def ev(t, y):
        return y[1]
    ev.direction = -1
    s = solve_ivp(f, (0, 40), [1.0, 0.0], events=ev, rtol=1e-12, atol=1e-14)
    k = [i for i, t in enumerate(s.t_events[0]) if t > 1e-6][0]
    return s.y_events[0][k][0], s.t_events[0][k]


def brake_intervals(t, v, thresh):
    """Start and end times of the intervals where the brake is engaged."""
    on = v > thresh
    edges = np.flatnonzero(np.diff(on.astype(int)))
    starts = list(t[edges[on[edges + 1]] + 1])
    ends = list(t[edges[~on[edges + 1]] + 1])
    if on[0]:
        starts.insert(0, t[0])
    if on[-1]:
        ends.append(t[-1])
    return list(zip(starts, ends))


def check_governor():
    print("Maxwell's governor")
    print("  true governor, one-sided brake: amplitude ratio per cycle")
    print("     F     c   zeta+  zeta-   integrated   e^-(d+ + d-)   verdict")
    for F, c in [(0.6, 0.0), (0.3, 0.2), (0.4, 0.2), (0.8, 0.2), (0.6, 0.6),
                 (1.0, 0.4)]:
        zp, zm = zetas(F, c)
        r, T = cycle_ratio(zp, zm)
        v = ("neutral" if abs(r - 1) < 1e-9 else
             "decays" if r < 1 else "grows")
        print(f"  {F:>4.1f}  {c:>4.1f}  {zp:>+5.2f}  {zm:>+5.2f}   {r:.6f}"
              f"     {np.exp(-(delta(zp) + delta(zm))):.6f}     {v}")
    print("  classification when a region is overdamped (rule / observed):")
    for F, c in [(2.4, 0.4), (2.8, 0.4), (3.0, 2.4), (4.8, 2.4)]:
        zp, zm = zetas(F, c)
        print(f"     F={F:.1f} c={c:.1f}  zeta+={zp:+.1f} zeta-={zm:+.1f}:  "
              f"{stability.classify(zp, zm)} / "
              f"{stability.observe(zp, zm, n=16)}")

    print("  overspeed brake set v0 above synchronous (offset boundary):")
    print("     F     c   zeta+  zeta-   mean     r*/v0      T")
    for F, c in [(0.8, 0.2), (0.5, 0.2), (1.2, 0.2), (0.3, 0.2), (0.8, 0.6)]:
        zp, zm = zetas(F, c)
        r = offset_rstar(zp, zm)
        ok = np.isfinite(r)
        T = frequency.period_exact(zp, zm)[0] if ok else np.nan
        got = f"{r:8.4f}  {T:8.4f}" if ok else "grows unbounded"
        print(f"  {F:>4.1f}  {c:>4.1f}  {zp:>+5.2f}  {zm:>+5.2f}  {(zp+zm)/2:+.3f}"
              f"   {got}")
    print("  governor deadband at synchronous (symmetric band):")
    print("     F     c   zeta+  zeta-    r*/v0      T")
    for F, c in [(0.8, 0.2), (0.8, 0.6), (0.7, 0.6)]:
        zp, zm = zetas(F, c)
        r, T = symmetric.cycle_integrated(zp, zm)
        got = f"{r:8.4f}  {T:8.4f}" if np.isfinite(r) else "grows unbounded"
        print(f"  {F:>4.1f}  {c:>4.1f}  {zp:>+5.2f}  {zm:>+5.2f}   {got}")


# ============================================================== oscillator
# Shared cycle finder, then the LC oscillator.

def driven_cycle(f, r0, n=400, tol=1e-11, cap=None):
    """Limit cycle by iterating the return map on ``{thetadot = 0}``.

    Returns ``(R+, R-, T, status)``: the positive and negative extremes,
    the period, and ``"cycle"``, ``"over the top"`` (the displacement
    reached ``cap``) or ``"no return"``.
    """
    def ev(t, y):
        return y[1]
    ev.direction = -1

    def top(t, y):
        return abs(y[0]) - (np.inf if cap is None else cap)
    top.terminal = True
    r, T, neg = r0, np.nan, np.nan
    for _ in range(n):
        s = solve_ivp(f, (0, 60), [r, 0.0], events=[ev, top], rtol=1e-12,
                      atol=1e-14, dense_output=True)
        if len(s.t_events[1]):
            return np.nan, np.nan, np.nan, "over the top"
        i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
        if not i:
            return np.nan, np.nan, np.nan, "no return"
        rn, T = s.y_events[0][i[0]][0], s.t_events[0][i[0]]
        neg = s.sol(np.linspace(0, T, 4000))[0].min()
        if abs(rn - r) < tol*max(1.0, abs(r)):
            return rn, neg, T, "cycle"
        r = rn
    return r, neg, T, "cycle"


# ============================================================== oscillator
# LC tank, loss conductance G, tickler feedback current f(v) from the
# transistor. C v'' + (G - f'(v)) v' + v/L = 0. With w0 = 1/sqrt(LC),
# Q = w0 C / G, small-signal loop gain A = gm / G, and v measured in units
# of the linear-range half-width I0/gm:
#
#     zeta(v) = (1/2Q) (1 - A g(v)),   g = 1 inside |v| < 1 else 0  (clip)
#                                       g = sech^2(v)               (tanh)

def oscillator_field(Q, A, smooth=False):
    zp = 1/(2*Q)

    def f(t, y):
        v, dv = y
        g = 1/np.cosh(v)**2 if smooth else (1.0 if abs(v) < 1 else 0.0)
        return [dv, -v - 2*zp*(1 - A*g)*dv]
    return f


def oscillator_cycle(Q, A, smooth=False, tol=1e-11):
    """Amplitude and period of the oscillator's cycle, from the energy seed.

    Near the start-up threshold the return map barely contracts, so a
    tight tolerance costs hundreds of turns; the figure sweep passes a
    looser one, still far below plotting resolution.
    """
    zp = 1/(2*Q)
    R0 = displacement.amplitude(zp, zp*(1 - A), 1.0, symmetric=True)
    R, _, T, st = driven_cycle(oscillator_field(Q, A, smooth), R0, cap=None,
                               tol=tol)
    return R, T


def check_oscillator():
    print("\nTransistor LC oscillator: clip (the prototype) against tanh")
    print("    Q    A   zeta+   zeta-    R clip   R energy   R tanh   "
          "T clip    T tanh")
    for Q, A in [(10, 1.5), (10, 2.0), (10, 3.0), (10, 5.0), (3, 2.0),
                 (3, 5.0)]:
        zp, zm = 1/(2*Q), (1 - A)/(2*Q)
        Rc, Tc = oscillator_cycle(Q, A)
        Rt, Tt = oscillator_cycle(Q, A, smooth=True)
        Re = displacement.amplitude(zp, zm, 1.0, symmetric=True)
        Tx = displacement.period_exact(zp, zm, symmetric=True)[0]
        print(f"   {Q:>2}  {A:>3.1f}  {zp:+.3f}  {zm:+.3f}   {Rc:.4f}   "
              f"{Re:.4f}   {Rt:.4f}   {Tc:.5f}  {Tt:.5f}   "
              f"(exact reduction {Tx:.5f})")
    print("  clipping on one side only (asymmetric model), Q = 10:")
    print("    A    zeta-    mean     R integrated   energy balance   T")
    for A in [1.2, 1.5, 1.8, 1.95, 2.5]:
        zp, zm = 0.05, (1 - A)*0.05
        R, T = displacement.cycle_integrated(zp, zm, 1.0, symmetric=False)
        Re = displacement.amplitude(zp, zm, 1.0, symmetric=False)
        got = (f"{R:12.4f}   {Re:12.4f}   {T:.5f}" if np.isfinite(Re)
               else "grows unbounded")
        print(f"   {A:>4.2f}  {zm:+.3f}  {(zp+zm)/2:+.4f}   {got}")
    print("  large loop gain: R / (A v0) against the hard limiter's 4/pi "
          f"= {4/np.pi:.4f}")
    for A in [2, 5, 10, 20, 50]:
        print(f"     A={A:>3}:  energy balance {displacement.amplitude(0.05, 0.05*(1 - A), 1.0, True)/A:.4f}")


# ================================================================= figures
def fig_governor(th, name):
    """Speed histories for the three governor configurations."""
    fig, axes = newfig(th, 1, 3, figsize=(13.5, 4.2))
    c, v0, T = 0.2, 1.0, 46.0

    ax = axes[0]
    for F, col, lab in [(0.3, th["series"][1], "$F=0.3$, grows"),
                        (0.4, th["series"][2], "$F=0.4=2c$, neutral"),
                        (0.8, th["series"][0], "$F=0.8$, decays")]:
        x, v, t = traj(governor_field(F, c, "integral"), [0.0, 0.6], T)
        ax.plot(t, v, color=col, linewidth=1.3, label=lab, zorder=3)
    ax.axhline(0, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
               zorder=5)
    ax.text(0.985, 0.0, "brake on above $V$",
            transform=ax.get_yaxis_transform(), fontsize=8, color=th["ink2"],
            ha="right", va="bottom", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))
    style(ax, th, "$t\\,\\omega_n$", "$\\dot{x}-V$",
          "True governor, one-sided brake ($c=0.2$)")
    ax.set_ylim(-2.2, 2.2)
    legend(ax, th, loc="lower left")

    ax = axes[1]
    F = 0.8
    x, v, t = traj(governor_field(F, c, "offset", v0), [0.0, 0.05], T)
    ax.plot(t, v, color=th["series"][0], linewidth=1.3,
            label="speed error, from a small hunt", zorder=3)
    for a, b in brake_intervals(t, v, v0):
        ax.axvspan(a, b, color=th["grid"], zorder=0)
    ax.axhline(v0, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
               zorder=5)
    ax.text(0.985, v0, "brake set at $\\omega_s+v_0$",
            transform=ax.get_yaxis_transform(), fontsize=8, color=th["ink2"],
            ha="right", va="bottom", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))
    ax.axhline(0, color=th["axis"], linewidth=0.8, zorder=1)
    style(ax, th, "$t\\,\\omega_n$", "$\\dot{\\delta}$",
          "Overspeed brake on a synchronised machine ($F=0.8$)")
    ax.set_ylim(-3.2, 3.2)
    legend(ax, th, loc="lower left")

    ax = axes[2]
    x, v, t = traj(governor_field(F, c, "deadband", v0), [0.0, 0.05], T)
    ax.plot(t, v, color=th["series"][0], linewidth=1.3,
            label="speed error, from a small hunt", zorder=3)
    ax.axhspan(-v0, v0, color=th["grid"], zorder=0)
    for y in (v0, -v0):
        ax.axhline(y, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
                   zorder=5)
    ax.text(0.985, 0.0, "deadband  $|\\dot{\\delta}| < v_0$",
            transform=ax.get_yaxis_transform(), fontsize=8, color=th["ink2"],
            ha="right", va="center", zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))
    style(ax, th, "$t\\,\\omega_n$", "$\\dot{\\delta}$",
          "Governor deadband at synchronous speed ($F=0.8$)")
    ax.set_ylim(-3.2, 3.2)
    legend(ax, th, loc="lower left")

    fig.suptitle("Maxwell's governor: where the brake threshold sits relative "
                 "to the equilibrium decides everything", color=th["ink"],
                 fontsize=11)
    fig.tight_layout()
    save(fig, name, "example-governor")


_CURVES = {}


def _amplitude_curve(Q, As, smooth):
    """Cycle amplitudes over loop gain, computed once and reused per theme."""
    key = (Q, smooth, len(As))
    if key not in _CURVES:
        _CURVES[key] = [oscillator_cycle(Q, a, smooth)[0] for a in As]
    return _CURVES[key]


def fig_oscillator(th, name):
    """Clip against tanh saturation in a transistor LC oscillator."""
    Q = 10
    fig, axes = newfig(th, 1, 2, figsize=(11.0, 4.8))
    ax = axes[0]
    A = 3.0
    ax.axvspan(-1, 1, color=th["grid"], zorder=0)
    for xv in (1, -1):
        ax.axvline(xv, color=th["ink2"], linewidth=1.2, linestyle=(0, (5, 3)),
                   zorder=5)
    fc = oscillator_field(Q, A)
    x, v, _ = traj(fc, [0.05, 0.0], 90)
    ax.plot(x, v, color=th["series"][1], linewidth=0.9,
            label="start-up from noise", zorder=3)
    R, T = oscillator_cycle(Q, A)
    x, v, _ = traj(fc, [R, 0.0], T)
    ax.plot(x, v, color=th["series"][0], linewidth=2.4,
            label="clip: the prototype", zorder=4)
    R, T = oscillator_cycle(Q, A, smooth=True)
    x, v, _ = traj(oscillator_field(Q, A, True), [R, 0.0], T)
    ax.plot(x, v, color=th["series"][2], linewidth=1.6,
            linestyle=(0, (5, 3)), label="tanh saturation", zorder=4)
    ax.text(0.0, 0.03, "transistor linear\n$|v| < v_0$", fontsize=8,
            color=th["ink2"], ha="center", va="bottom", zorder=6,
            transform=ax.get_xaxis_transform(),
            bbox=dict(boxstyle="round,pad=0.25", fc=th["surface"], ec="none"))
    style(ax, th, "$v / v_0$", "$\\dot{v} / (\\omega_0 v_0)$",
          f"$Q = {Q}$, loop gain $A = {A:.0f}$")
    ax.set_aspect("equal")
    legend(ax, th, loc="upper right")

    ax = axes[1]
    As = np.linspace(1.25, 8.0, 24)
    for smooth, col, lab, lw in [(False, th["series"][0], "clip", 2.0),
                                 (True, th["series"][2], "tanh", 1.6)]:
        ax.plot(As, [oscillator_cycle(Q, a, smooth, tol=1e-7)[0] for a in As],
                color=col, linewidth=lw, label=lab, zorder=3)
    ax.plot(As, [displacement.amplitude(1/(2*Q), (1 - a)/(2*Q), 1.0, True)
                 for a in As], color=th["series"][1], linewidth=1.4,
            linestyle=(0, (5, 3)), label="energy balance $v_0/\\cos\\phi$",
            zorder=3)
    A1 = np.linspace(1.05, 1.97, 60)
    ax.plot(A1, [displacement.amplitude(1/(2*Q), (1 - a)/(2*Q), 1.0, False)
                 for a in A1], color=th["ink2"], linewidth=1.2,
            linestyle=(0, (2, 2)), label="clip on one side only", zorder=3)
    ax.axvline(2.0, color=th["axis"], linewidth=0.8, zorder=1)
    ax.set_ylim(0, 11)
    style(ax, th, "small-signal loop gain $A = g_m/G$", "$R / v_0$",
          "Amplitude against loop gain")
    legend(ax, th, loc="upper left")
    fig.suptitle("Transistor LC oscillator: the gain switch is the piecewise "
                 "Van der Pol", color=th["ink"], fontsize=11)
    fig.tight_layout()
    save(fig, name, "example-oscillator")


if __name__ == "__main__":
    import sys
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "checks"):
        check_governor()
        check_oscillator()
        print()
    if what not in ("all", "figures"):
        sys.exit(0)
    for name, th in THEMES.items():
        print(f"{name}:")
        fig_governor(th, name)
        fig_oscillator(th, name)
