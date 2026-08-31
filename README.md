# Second order nonlinear prototype 

This project covers the second order nonlinear prototype. The starting point is
the second order **linear** prototype, described below; the nonlinear form
builds on it.

## Second order linear prototype

As a single second order ordinary differential equation:

$$\ddot{x} + 2\zeta\omega_n\dot{x} + \omega_n^2 x = \omega_n^2 u(t)$$

where $\omega_n$ is the undamped natural frequency, $\zeta$ is the damping
ratio, and $u(t)$ is the input. In raw coefficient form this is

$$m\ddot{x} + c\dot{x} + kx = f(t)$$

with $\omega_n = \sqrt{k/m}$ and $\zeta = c / (2\sqrt{km})$.

## As a set of first order ordinary differential equations

Taking the states as position and velocity, $x_1 = x$ and $x_2 = \dot{x}$:

$$
\begin{aligned}
\dot{x}_1 &= x_2 \\
\dot{x}_2 &= -\omega_n^2 x_1 - 2\zeta\omega_n x_2 + \omega_n^2 u
\end{aligned}
$$

The second derivative $\ddot{x}$ no longer appears: one second order equation
has become two coupled first order equations, with $\dot{x}_2$ carrying what
was $\ddot{x}$.

## State space form

$$
\begin{bmatrix}\dot{x}_1 \\ \dot{x}_2\end{bmatrix}
=
\begin{bmatrix}0 & 1 \\ -\omega_n^2 & -2\zeta\omega_n\end{bmatrix}
\begin{bmatrix}x_1 \\ x_2\end{bmatrix}
+
\begin{bmatrix}0 \\ \omega_n^2\end{bmatrix} u,
\qquad
y = \begin{bmatrix}1 & 0\end{bmatrix}
\begin{bmatrix}x_1 \\ x_2\end{bmatrix}
$$
