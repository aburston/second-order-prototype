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

## Nonlinear prototype: a linear boundary along the x-axis

The first nonlinear prototype keeps the dynamics linear on either side of a
single straight switching boundary. Taking the phase plane with $x$ on the
horizontal axis and $\dot{x}$ on the vertical axis, that boundary is the
x-axis itself:

```math
\Sigma = \{ (x_1, x_2) \in \mathbb{R}^2 : h(x_1, x_2) = x_2 = 0 \}
```

It splits the plane into the two open half planes

```math
S^{+} = \{ x_2 > 0 \}, \qquad S^{-} = \{ x_2 < 0 \}
```

Since $x_2 = \dot{x}$, the boundary is crossed whenever the velocity changes
sign, and the natural nonlinearity to attach to it is a force that depends on
the direction of motion — Coulomb friction:

```math
\ddot{x} + 2\zeta\omega_n\dot{x} + \omega_n^2 x + \mu\,\mathrm{sign}(\dot{x}) = \omega_n^2 u(t)
```

where $\mu = F_c / m \ge 0$ is the Coulomb friction force per unit mass. Setting
$\mu = 0$ recovers the linear prototype exactly.

### As a set of first order ordinary differential equations

With $x_1 = x$ and $x_2 = \dot{x}$ as before:

```math
\begin{aligned}
\dot{x}_1 &= x_2 \\
\dot{x}_2 &= -\omega_n^2 x_1 - 2\zeta\omega_n x_2 - \mu\,\mathrm{sign}(x_2) + \omega_n^2 u
\end{aligned}
```

The right hand side is discontinuous across $\Sigma$, so this is not a single
smooth vector field but two of them, one per half plane:

```math
f^{\pm}(x) =
\begin{bmatrix}
x_2 \\
-\omega_n^2 x_1 - 2\zeta\omega_n x_2 \mp \mu + \omega_n^2 u
\end{bmatrix}
\quad \text{on } S^{\pm}
```

Each of $f^{+}$ and $f^{-}$ is affine, so the system is piecewise linear: the
prototype is linear everywhere except on the boundary itself.

### Behaviour on the boundary

Because $f^{+} \neq f^{-}$ on $\Sigma$, the solution there is taken in the
Filippov sense, as the convex combination of the two fields:

```math
\dot{x} \in F(x) = \{ (1-\alpha) f^{-}(x) + \alpha f^{+}(x) : \alpha \in [0,1] \},
\qquad x \in \Sigma
```

Which of the two behaviours occurs is decided by the components of $f^{\pm}$
normal to $\Sigma$. With $\nabla h = (0, 1)^{\top}$ and $x_2 = 0$:

```math
\nabla h \cdot f^{\pm} = \omega_n^2 (u - x_1) \mp \mu
```

**Sticking.** Both fields point at the boundary when
$\nabla h \cdot f^{+} < 0$ and $\nabla h \cdot f^{-} > 0$, which reduces to a
single condition on the net applied force:

```math
\lvert \omega_n^2 (u - x_1) \rvert < \mu
```

The set of such points is the stick band

```math
\hat{\Sigma} = \left\{ (x_1, 0) : u - \frac{\mu}{\omega_n^2} < x_1 < u + \frac{\mu}{\omega_n^2} \right\}
```

On $\hat{\Sigma}$ the Filippov combination gives $\dot{x}_2 = 0$, and
$\dot{x}_1 = x_2 = 0$ already holds, so the sliding vector field vanishes: the
whole band is a continuum of equilibria. The mass stops and stays stopped
until $u(t)$ moves far enough to break it out. This is the qualitative
departure from the linear prototype, which has the single equilibrium
$x_1 = u$.

**Crossing.** Where the two normal components share a sign, that is
$\lvert \omega_n^2 (u - x_1) \rvert > \mu$, the trajectory passes straight
through $\Sigma$: the velocity changes sign and the friction term switches
with it, so the mass reverses direction without ever stopping.

### Notes for numerical work

- $\mathrm{sign}(\cdot)$ is discontinuous, so a fixed step integrator will
  chatter around $\Sigma$. Either use event detection to locate the crossing
  and restart the integration, or regularise the switch, for example with
  $\tanh(x_2 / \varepsilon)$ for a small $\varepsilon > 0$.
- The model above takes the static and kinetic friction coefficients to be
  equal. Allowing $\mu_s > \mu_k$ widens the stick band relative to the
  sliding force and is the usual next refinement.
