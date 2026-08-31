# Second order nonlinear prototype 

This project covers the second order nonlinear prototype. The starting point is
the second order **linear** prototype; the nonlinear prototype builds on it by
introducing a switching boundary in the phase plane.

## Second order linear prototype

As a single second order ordinary differential equation:

```math
\ddot{x} + 2\zeta\omega_n\dot{x} + \omega_n^2 x = \omega_n^2 u(t)
```

where $\omega_n$ is the undamped natural frequency, $\zeta$ is the damping
ratio, and $u(t)$ is the input. In raw coefficient form this is

```math
m\ddot{x} + c\dot{x} + kx = f(t)
```

with $\omega_n = \sqrt{k/m}$ and $\zeta = c / (2\sqrt{km})$.

## As a set of first order ordinary differential equations

Taking the states as position and velocity, $x_1 = x$ and $x_2 = \dot{x}$:

```math
\begin{aligned}
\dot{x}_1 &= x_2 \\
\dot{x}_2 &= -\omega_n^2 x_1 - 2\zeta\omega_n x_2 + \omega_n^2 u
\end{aligned}
```

The second derivative $\ddot{x}$ no longer appears: one second order equation
has become two coupled first order equations, with $\dot{x}_2$ carrying what
was $\ddot{x}$.

## State space form

```math
\begin{bmatrix}\dot{x}_1 \\ \dot{x}_2\end{bmatrix}
=
\begin{bmatrix}0 & 1 \\ -\omega_n^2 & -2\zeta\omega_n\end{bmatrix}
\begin{bmatrix}x_1 \\ x_2\end{bmatrix}
+
\begin{bmatrix}0 \\ \omega_n^2\end{bmatrix} u
```

```math
y = \begin{bmatrix}1 & 0\end{bmatrix}
\begin{bmatrix}x_1 \\ x_2\end{bmatrix}
```

## Nonlinear prototype: switched damping across a boundary on the x-axis

The first nonlinear prototype keeps the dynamics linear on either side of a
single straight switching boundary, and puts the nonlinearity in the **damping
coefficient**. Taking the phase plane with $x$ on the horizontal axis and
$\dot{x}$ on the vertical axis, the boundary is the x-axis itself:

```math
\Sigma = \{ (x_1, x_2) \in \mathbb{R}^2 : h(x_1, x_2) = x_2 = 0 \}
```

splitting the plane into the two open half planes

```math
S^{+} = \{ x_2 > 0 \}, \qquad S^{-} = \{ x_2 < 0 \}
```

Since $x_2 = \dot{x}$, the boundary is crossed whenever the velocity changes
sign. The damping ratio takes a different value on each side:

```math
\ddot{x} + 2\zeta(\dot{x})\,\omega_n\dot{x} + \omega_n^2 x = \omega_n^2 u(t),
\qquad
\zeta(\dot{x}) =
\begin{cases}
\zeta_{+} & \dot{x} > 0 \\
\zeta_{-} & \dot{x} < 0
\end{cases}
```

Writing the mean and half difference of the two damping ratios as
$\bar{\zeta} = (\zeta_{+} + \zeta_{-})/2$ and
$\Delta\zeta = (\zeta_{+} - \zeta_{-})/2$, the switch can be folded into a
single term and the nonlinearity shows up as an absolute value:

```math
\ddot{x} + 2\omega_n\bar{\zeta}\,\dot{x} + 2\omega_n\Delta\zeta\,\lvert\dot{x}\rvert + \omega_n^2 x = \omega_n^2 u(t)
```

Setting $\Delta\zeta = 0$ recovers the linear prototype exactly.

### As a set of first order ordinary differential equations

With $x_1 = x$ and $x_2 = \dot{x}$ as before:

```math
\begin{aligned}
\dot{x}_1 &= x_2 \\
\dot{x}_2 &= -\omega_n^2 x_1 - 2\omega_n\left(\bar{\zeta}x_2 + \Delta\zeta\lvert x_2 \rvert\right) + \omega_n^2 u
\end{aligned}
```

Equivalently, as one linear system per half plane, $\dot{x} = A^{\pm}x + b\,u$:

```math
A^{\pm} =
\begin{bmatrix}
0 & 1 \\
-\omega_n^2 & -2\zeta_{\pm}\omega_n
\end{bmatrix},
\qquad
b = \begin{bmatrix} 0 \\ \omega_n^2 \end{bmatrix}
\qquad \text{on } S^{\pm}
```

with eigenvalues $\lambda_{\pm} = \omega_n\left(-\zeta_{\pm} \pm \sqrt{\zeta_{\pm}^2 - 1}\right)$.

### The field is continuous across the boundary

This is the structural difference from putting the switch in a *force* term.
The two fields differ only in the damping term $-2\zeta_{\pm}\omega_n x_2$,
and that term vanishes on $\Sigma$ where $x_2 = 0$. So

```math
f^{+}(x) = f^{-}(x) \quad \text{for all } x \in \Sigma
```

The vector field is continuous everywhere — only its Jacobian jumps across
$\Sigma$. The system is therefore piecewise smooth and continuous, and
Lipschitz, so ordinary solutions exist and are unique: no Filippov convex
combination is needed, there is no sliding or sticking set, and every
trajectory crosses $\Sigma$ transversally except at an equilibrium. The
boundary changes the *rate* at which the system gains or loses energy, not the
force acting on it.

### Stability

Take $u = 0$, so the equilibrium is the origin. The system is not
differentiable there, so stability is decided by the half-cycle map rather
than by linearisation.

Within $S^{-}$, a trajectory leaving the boundary at $(x_1, 0)$ with
$x_1 = A_n > 0$ follows the linear system with damping $\zeta_{-}$, and for
$\lvert\zeta_{-}\rvert < 1$ returns to the boundary after exactly half a
damped period at $(-A_n e^{-\delta(\zeta_{-})},\, 0)$, where the logarithmic
decrement per half cycle is

```math
\delta(\zeta) = \frac{\pi\zeta}{\sqrt{1-\zeta^2}}
```

The next half cycle runs through $S^{+}$ with $\zeta_{+}$, so the amplitude
after a full cycle is

```math
A_{n+1} = A_n\,e^{-\left(\delta(\zeta_{+}) + \delta(\zeta_{-})\right)}
```

The origin is asymptotically stable exactly when
$\delta(\zeta_{+}) + \delta(\zeta_{-}) > 0$. Since $\delta$ is odd and
strictly increasing on $(-1, 1)$, this collapses to a condition on the mean
damping alone:

```math
\zeta_{+} + \zeta_{-} > 0
\qquad \Longleftrightarrow \qquad
\bar{\zeta} > 0
```

So one half plane may have **negative** damping and the origin still be
asymptotically stable, provided the other half plane damps harder. What
matters is the average over a cycle, not the sign on either side.

Three cases follow:

| condition | behaviour |
| --- | --- |
| $\bar{\zeta} > 0$ | origin asymptotically stable |
| $\bar{\zeta} = 0$ | full cycle map is the identity — a continuum of closed orbits |
| $\bar{\zeta} < 0$ | origin unstable, oscillation grows without bound |

**There is no limit cycle.** With $u = 0$ the field is positively homogeneous
of degree one, $f(\lambda x) = \lambda f(x)$ for $\lambda > 0$, so the return
map is an exact scaling and its behaviour cannot depend on amplitude.
Stability is global, and the marginal case gives a continuum of periodic
orbits rather than an isolated one. Producing an isolated limit cycle requires
amplitude dependence — a $\zeta$ that varies with $x$, or a boundary that does
not pass through the equilibrium.

For constant $u \neq 0$ the equilibrium moves to $(u, 0)$, which still lies on
$\Sigma$, and the same analysis applies in the shifted coordinate $x_1 - u$.

### Notes for numerical work

- Both half planes must be underdamped, $\lvert\zeta_{\pm}\rvert < 1$, for the
  half cycle map above to be defined. If the trajectory enters a half plane
  with $\zeta_{\pm} \ge 1$ it decays monotonically to the equilibrium without
  recrossing $\Sigma$; with $\zeta_{\pm} \le -1$ it diverges monotonically.
- The field is continuous, so a general purpose integrator will not chatter
  the way it would across a discontinuous force. Accuracy still drops at each
  crossing because the Jacobian jumps, so use event detection on $x_2 = 0$ to
  place a step boundary at the switch.
- $\bar{\zeta} = 0$ is non-hyperbolic and structurally fragile: integration
  error alone will make the orbits drift in or out.
