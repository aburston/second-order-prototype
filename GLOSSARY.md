# Glossary

The documents in this repository borrow vocabulary from three fields that do
not usually meet: **control engineering** (damping ratios, poles, the
s-plane), **nonlinear dynamics** (limit cycles, attractors, chaos), and the
smaller literature on **non-smooth systems** (Filippov solutions, saltation,
grazing). A term that is everyday in one of those is often unfamiliar in the
others.

Each entry below gives a plain explanation first, then says where the term is
used here and why it matters, then points at a standard text. The
explanations assume first-year university mathematics — calculus, complex
numbers, matrices and eigenvalues — and nothing else.

Reference keys in brackets are resolved in [References](#references) at the
end.

Throughout, $`x`$ is a displacement, $`\dot{x}`$ its velocity, $`\omega_n`$ the
natural frequency and $`\zeta`$ the damping ratio.

---

## 1. The linear oscillator

Everything in this repository is built from one equation, so its vocabulary
comes first.

```math
\ddot{x} + 2\zeta\omega_n\dot{x} + \omega_n^2 x = 0
```

### Natural frequency, $`\omega_n`$

How fast the system would oscillate with no damping at all, in radians per
second. For a mass on a spring, $`\omega_n = \sqrt{k/m}`$. Divide by
$`2\pi`$ to get cycles per second.

*Reference:* [Åström & Murray], ch. 6; [Ogata], ch. 5.

### Damping ratio, $`\zeta`$ ("zeta")

A dimensionless number saying how much the oscillation is being bled away
relative to how fast it swings. It is the single most important quantity in
this repository because it is the thing the prototypes **switch**.

- $`\zeta = 0`$ — no damping; the oscillation continues for ever.
- $`0 \lt \zeta \lt 1`$ — **underdamped**; it oscillates while decaying.
- $`\zeta = 1`$ — critically damped; returns without overshoot, as fast as
  possible.
- $`\zeta \gt 1`$ — **overdamped**; returns slowly without oscillating.
- $`\zeta \lt 0`$ — *negative* damping. Energy is being *added*, so the
  oscillation grows. Physically this means something is pumping the system —
  an engine, an amplifier, a belt dragging past a block. Several prototypes
  here use a negative $`\zeta`$ on one side of a boundary, which is what makes
  them self-exciting.

*Reference:* [Åström & Murray], ch. 6; [Ogata], ch. 5.

### Damped natural frequency, $`\omega_d`$

The frequency actually observed when damping is present:
$`\omega_d = \omega_n\sqrt{1-\zeta^2}`$. Slightly slower than $`\omega_n`$.
When $`\lvert\zeta\rvert \gt 1`$ the square root turns imaginary, the cosine
becomes a **cosh**, and the motion stops oscillating — which is why the code
here computes with $`\omega_d^2`$ rather than $`\omega_d`$, so one formula
covers both cases.

*Reference:* [Ogata], ch. 5.

### Logarithmic decrement, $`\delta`$

How much the amplitude shrinks in one half-swing, measured as a logarithm:

```math
\delta(\zeta) = \frac{\pi\zeta}{\sqrt{1-\zeta^2}}
```

If a peak has height $`A`$, the next peak on the other side has height
$`A e^{-\delta}`$. It is the classic way to measure damping from a ringdown:
count how fast the peaks decay and you can solve back for $`\zeta`$. In this
repository it is the building block for every stability result, because the
switched models are just two different decrements alternating.

*Reference:* [Ogata], ch. 5; [Åström & Murray], ch. 6.

### Ringdown

The free decaying oscillation you get by disturbing a system and letting go —
struck bell, plucked string, nudged governor. `speculation.md` uses it as the
basic field measurement, because the decrement is readable straight off it.

### Pole, s-plane, pole–zero diagram

Write the solution as $`e^{st}`$; the values of $`s`$ that work are the
**poles**. For our oscillator they are

```math
s = \omega_n\left(-\zeta \pm \sqrt{\zeta^2-1}\right)
```

Plotting them on the complex plane (the **s-plane**) makes stability visual:
poles in the left half decay, poles in the right half grow, and how far they
sit from the imaginary axis says how fast. `README.md` uses this to classify
every combination of the two damping ratios.

*Reference:* [Franklin, Powell & Emami-Naeini], ch. 3; [Åström & Murray],
ch. 6.

### Unit circle, z-domain

The discrete-time counterpart of the s-plane. When a system is sampled — or,
as here, observed once per cycle — stability is judged by whether the
relevant multipliers lie **inside the unit circle** in the complex plane
rather than in the left half-plane. `MAPS.md` uses this to classify cycles.

*Reference:* [Franklin, Powell & Emami-Naeini], ch. 8.

### State space, phase plane

Instead of one second-order equation, write two first-order ones by treating
velocity as a second variable: $`x_1 = x`$, $`x_2 = \dot{x}`$. The
**state space** is the space of all $`(x_1, x_2)`$ pairs; for a
second-order system it is a plane, so it is called the **phase plane**. A
solution is then a *curve* in that plane rather than a graph against time,
and questions about behaviour become questions about geometry.

*Reference:* [Strogatz], ch. 5–6.

---

## 2. Geometry in the phase plane

### Trajectory / orbit

The curve a solution traces in the phase plane. "Orbit" is the same thing;
dynamicists prefer it.

### Equilibrium (fixed point)

A state where nothing changes — the system sits there for ever. For our
oscillator it is the origin. **Stable** means nearby states are drawn in,
**unstable** means they are pushed away.

*Reference:* [Strogatz], ch. 5.

### Limit cycle

A closed loop in the phase plane that the system settles onto: a
self-sustaining oscillation with its **own** amplitude, set by the system
rather than by how you started it. This is the central object in the
repository.

The crucial contrast is with a linear oscillator, which either decays to
nothing or grows without bound, and whose closed orbits (in the marginal
case) come in a continuous family — nudge it and you get a different orbit,
which stays. A limit cycle is **isolated**: nudge it and it comes back.

*Reference:* [Strogatz], ch. 7; [Andronov, Vitt & Khaikin].

### Isolated / hyperbolic cycle

**Isolated** means no other closed orbit sits arbitrarily close to it.
**Hyperbolic** means it attracts or repels at an exponential rate — the
multiplier (below) is strictly inside or outside the unit circle, not on it.
Hyperbolic cycles are robust: a small change to the equations moves them
slightly but does not destroy them, which is why the repository cares whether
a cycle is hyperbolic or merely marginal.

*Reference:* [Guckenheimer & Holmes], ch. 1.

### Basin of attraction

The set of starting states that end up on a particular attractor. When a
system has two attractors — as the two-threshold prototype does, with the
origin and an outer cycle — the phase plane splits into two basins, and the
boundary between them is often an unstable cycle.

*Reference:* [Strogatz], ch. 7–8.

### Bistability, hard excitation

**Bistable**: two stable states coexist for the *same* parameter values, and
which one you get depends on where you started. **Hard excitation** is the
practical consequence: the system sits quietly at rest, and a small
disturbance dies away, but a large enough one pushes it past the basin
boundary and it jumps to a big oscillation — and stays there. A common and
unpleasant failure mode in machinery.

*Reference:* [Strogatz], ch. 8; [Nayfeh & Mook], ch. 3.

### Separatrix

A curve dividing the phase plane into regions with different fates. A basin
boundary is a separatrix.

### Invariant ray

A straight line through the origin that trajectories stay on once they reach
it. In `README.md` these appear when a region is overdamped: an overdamped
linear system has real poles, and each real pole direction is such a ray.
Because a trajectory cannot leave one, the ray's own stability decides the
fate of everything that lands on it.

### Self-excited oscillation

An oscillation that starts itself, drawing energy from a steady source, with
no oscillating input. Clock escapements, violin strings, and the transistor
oscillator in `EXAMPLES.md` are all self-excited. The signature is negative
damping at small amplitude and positive damping at large.

*Reference:* [Andronov, Vitt & Khaikin]; [Strogatz], ch. 7.

---

## 3. Piecewise and non-smooth systems

This is the least familiar vocabulary, and it is where the repository does
most of its careful work.

### Piecewise linear

The system is an ordinary linear system inside each region of the phase
plane, but *which* linear system depends on where you are. All the prototypes
here are of this kind: the damping ratio takes one value on one side of a
boundary and another value on the other side. The nonlinearity is entirely in
the switching, never in the dynamics once you have switched.

*Reference:* [di Bernardo et al.], ch. 1–2.

### Switching boundary, $`\Sigma`$ ("Sigma")

The surface in the phase plane where the rule changes. Here it is always a
straight line — either a line of constant velocity ($`\dot{x} = v_0`$) or a
line of constant displacement ($`x = x_0`$).

### Zone / region

One of the areas the boundaries divide the phase plane into. The
two-threshold prototype has five: a core, two middle zones, two outer ones.

### Deadzone (deadband)

A band around zero where a nonlinearity does nothing, switching on only
outside it. Here the damping takes its inner value throughout
$`\lvert\dot{x}\rvert \lt v_0`$ and its outer value beyond. Physically this
is a governor with slack, a valve that does not respond to small errors, or
an amplifier that ignores small signals.

*Reference:* [Åström & Murray], ch. 4 (as a common actuator nonlinearity).

### Continuous field, discontinuous field

The **field** is the arrow at each point saying which way the state moves. A
field is *continuous* across a boundary if the arrows match on both sides,
even though the rule generating them changed. This matters enormously: with a
continuous field, solutions behave normally. With a discontinuous one — the
arrows jump — you need extra theory (below).

The repository's velocity-switched models are deliberately built so the
switched term vanishes on the boundary, keeping the field continuous. The
displacement-switched ones cannot do this, and their field jumps.

*Reference:* [di Bernardo et al.], ch. 2.

### Filippov solution, sliding

When the field is discontinuous *and* the arrows on the two sides point
towards each other, no ordinary solution exists — the state cannot go either
way. **Filippov's** construction resolves this by allowing the state to
travel *along* the boundary, taking a blend of the two fields. This is called
**sliding**, and it is what happens physically when a block sticks to a belt
instead of slipping.

Sliding is a legitimate phenomenon but a complication, and the prototypes
here are designed to avoid it. Where the field is continuous it cannot occur
at all; where it jumps (the displacement models) the geometry happens to
prevent it, because the jump is parallel to the boundary rather than across
it.

*Reference:* [Filippov]; [di Bernardo et al.], ch. 2.

### Transversal crossing

The trajectory crosses the boundary cleanly, at an angle, and carries on. The
ordinary, well-behaved case — as opposed to sliding, or grazing.

### Grazing

The trajectory touches a boundary tangentially — reaching it exactly, with no
component of motion across it — and returns to the side it came from. Grazing
is the characteristic difficulty of non-smooth systems: on one side of a
grazing the orbit crosses into the next zone, on the other side it does not,
so the *sequence of zones visited* changes discontinuously. Quantities that
depend on that sequence are therefore not differentiable there.

`MAPS.md` records grazing as the reason the mapping approach was abandoned.

*Reference:* [di Bernardo et al.], ch. 6–7.

### Saltation matrix

Latin *saltus*, a leap. When you ask how a *small perturbation* evolves
across a switching boundary, it is not enough to apply each side's dynamics
in turn: a neighbouring trajectory reaches the boundary at a slightly
*different time*, and that time difference shears the perturbation. The
saltation matrix is the correction:

```math
S = I + \frac{\bigl(f^{+}-f^{-}\bigr) g^{\mathsf T}}{g^{\mathsf T} f^{-}}
```

with $`f^{-}, f^{+}`$ the field before and after and $`g`$ the boundary's
normal. When the field is continuous, $`f^{+} = f^{-}`$ and $`S`$ is the
identity — no correction needed.

*Reference:* [Leine & Nijmeijer], ch. 5; [di Bernardo et al.], ch. 2.

### Virtual centre

When damping acts on a *relative* velocity — for instance $`\dot{x} - v_0`$,
the speed relative to a moving belt — the oscillation in that zone is centred
not on the origin but on a shifted point. That point may lie outside the zone
where the rule applies, so the system never actually sits there: it is
"virtual", a centre the motion curves around without reaching.

### Positive homogeneity

A field is positively homogeneous if $`f(\lambda x) = \lambda f(x)`$ for
$`\lambda \gt 0`$ — scaling the state scales the velocity by the same factor.
Consequence: the picture looks identical at every magnification, so behaviour
cannot depend on amplitude, so **no isolated limit cycle is possible**. A
switching boundary through the equilibrium gives exactly this, which is why
the repository's first prototype has no limit cycle and why offsetting the
boundary is what creates one.

### Unfolding parameter

A parameter that, when turned on, breaks a degeneracy and reveals the
structure hiding inside it. Here the boundary offset $`v_0`$ is the unfolding
parameter: at $`v_0 = 0`$ there is a whole continuum of closed orbits, and any
$`v_0 \neq 0`$ collapses it to one isolated cycle.

*Reference:* [Kuznetsov], ch. 2–3.

---

## 4. Periodic orbits and their stability

### Poincaré section, return map

Rather than follow a trajectory continuously, watch it only when it crosses a
chosen surface — the **section** — and record where it lands each time. This
turns a continuous flow into a sequence of points, and a closed orbit into a
**fixed point** of the resulting **return map** (or Poincaré map). Much
easier to analyse: the question "is this cycle stable?" becomes "does this
map pull nearby points towards its fixed point?"

*Reference:* [Strogatz], ch. 8; [Guckenheimer & Holmes], ch. 1.

### Stroboscopic section

The same idea for a *driven* system, but sampling at a fixed rhythm — once
per drive period — rather than at a geometric surface. Like photographing a
spinning wheel with a strobe light: if it is locked to the flash you see a
frozen picture, if not you see it creep.

*Reference:* [Ott], ch. 1; [Pikovsky, Rosenblum & Kurths], ch. 3.

### Floquet multiplier

How much a small disturbance to a periodic orbit is multiplied over one full
cycle. $`\lvert\mu\rvert \lt 1`$ means disturbances shrink and the cycle is
stable; $`\lvert\mu\rvert \gt 1`$ means it is unstable. It is the discrete-time
analogue of a pole, which is why the unit circle replaces the imaginary axis.
For a planar autonomous system one multiplier is always exactly 1 — nudging
the state *along* the orbit just shifts its timing — and the other one is the
informative one.

*Reference:* [Guckenheimer & Holmes], ch. 1; [Kuznetsov], ch. 1.

### Monodromy matrix

The matrix that maps a perturbation at the start of a cycle to the
perturbation one cycle later. Its eigenvalues are the Floquet multipliers.
In `MAPS.md` it is built as a product of per-zone matrices and saltation
factors.

*Reference:* [Leine & Nijmeijer], ch. 5.

### Dwell time

How long the trajectory spends in one zone before crossing into the next. The
repository's central stability result is that a cycle's growth or decay is
the **dwell-weighted** sum of the pole real parts — that is, each zone
contributes in proportion to the time spent in it.

### Transition matrix, $`\Phi`$

The matrix that advances a *linear* system by a given time:
$`y(t) = \Phi(t)\,y(0)`$. Because each zone here is linear, one matrix
advances the state across a whole arc with no integration.

*Reference:* [Åström & Murray], ch. 6 ("matrix exponential").

### Amplitude equation

An equation whose solution is the cycle's amplitude, obtained by requiring
that the energy put in over a cycle balances the energy taken out. In this
repository it takes the recurring form
$`2\varphi - \sin 2\varphi = \pi\rho`$.

### Averaging, energy balance

A technique for weakly nonlinear oscillators: over one cycle the system is
nearly sinusoidal, so replace the exact equations with their average over
that cycle. The averaged damping is what decides whether the amplitude grows
or shrinks, and where it crosses zero is a limit cycle. Exact only in the
limit of weak damping, so the repository always checks it against an exact
calculation.

*Reference:* [Nayfeh & Mook], ch. 3; [Strogatz], ch. 7.

### Poincaré–Bendixson theorem

A powerful result about *planar* systems: a trajectory that stays in a
bounded region and does not approach an equilibrium must approach a closed
orbit. The consequence used repeatedly here is the negative one — **a planar
autonomous system cannot be chaotic**. Chaos needs a third dimension, which
is why the repository has to add forcing to look for it.

*Reference:* [Strogatz], ch. 7.

---

## 5. Forcing, entrainment and chaos

### Forcing (drive)

An externally imposed oscillation added to the equation, $`A\cos\Omega t`$.
Because the drive's phase is a third variable that keeps advancing, a forced
planar system is effectively three-dimensional — which is what opens the door
to chaos.

### Entrainment (phase locking, mode locking)

The oscillator abandons its own frequency and runs at the drive's, or at a
simple ratio of it. A pendulum clock nudged by a slightly-off rhythm ends up
keeping the rhythm's time.

*Reference:* [Pikovsky, Rosenblum & Kurths], ch. 3.

### Rotation (winding) number

How many turns the oscillator makes per drive period, on average. On a
$`p{:}q`$ lock it is exactly the rational $`p/q`$ and *stays* there over a
range of drive frequencies — a plateau. Between plateaus it varies smoothly
and is irrational.

*Reference:* [Pikovsky, Rosenblum & Kurths], ch. 3; [Arnol'd].

### Arnold tongue

Plot drive frequency horizontally and drive strength vertically, and shade
where the system is locked at a given ratio. The locked region is a wedge
that narrows to a point at zero drive and widens as drive grows — a "tongue".
Each rational ratio has its own.

*Reference:* [Pikovsky, Rosenblum & Kurths], ch. 3; [Arnol'd].

### Devil's staircase

The graph of rotation number against drive frequency: flat wherever the
system is locked, rising in between. Since locks occur at *every* rational
ratio, the result is a staircase with infinitely many steps — continuous,
non-decreasing, and constant almost everywhere.

*Reference:* [Pikovsky, Rosenblum & Kurths], ch. 3.

### Quasi-periodic motion, torus

Motion with two frequencies whose ratio is irrational — the oscillator's own
and the drive's — so it never exactly repeats. In the state space the
trajectory winds around the surface of a doughnut (a **torus**), and its
stroboscopic section is a closed curve rather than a finite set of points.

*Reference:* [Ott], ch. 1; [Strogatz], ch. 8.

### Chaos

Motion that is bounded, aperiodic, and **sensitively dependent on initial
conditions**: two starts differing by any tiny amount end up completely
different. Not randomness — the equations are deterministic and repeatable —
but unpredictability in practice, because you can never specify the start
precisely enough.

*Reference:* [Strogatz], ch. 9; [Ott], ch. 1.

### Sensitive dependence

The defining property above, made quantitative by the Lyapunov exponent. This
repository measures it directly: two starts $`10^{-10}`$ apart diverge to full
scale, while still tracing the same attractor to within a few per cent.

### Lyapunov exponent, $`\lambda`$

The average exponential rate at which nearby trajectories separate. Positive
means chaos, negative means they converge (a lock), zero means neither (a
torus). It is an *average over the attractor*, which is why it is
well-defined even though no individual trajectory is predictable.

*Reference:* [Ott], ch. 4; [Benettin et al.] for the standard algorithm.

### Attractor, strange attractor

The set a system settles onto after transients die: a point, a closed curve,
or something more complicated. A **strange** attractor is one with fractal
structure — layers within layers at every magnification — which is what
chaotic dissipative systems produce. Strong damping squeezes it very thin,
which is why the attractors in this repository look nearly one-dimensional.

*Reference:* [Ott], ch. 2–3.

### Period doubling

A cycle loses stability and is replaced by one taking *twice* as long, then
four times, and so on; the cascade accumulates and ends in chaos. In the
z-domain it is a multiplier leaving the unit circle through $`-1`$.

*Reference:* [Strogatz], ch. 10; [Kuznetsov], ch. 4.

### Fold (saddle-node) of cycles

Two limit cycles, one stable and one unstable, approach each other as a
parameter changes, collide, and both vanish. In the z-domain a multiplier
leaves through $`+1`$. This is how the two-threshold prototype's bistability
would be destroyed.

*Reference:* [Kuznetsov], ch. 3–5.

### Neimark–Sacker bifurcation

A cycle loses stability to a *torus*: a complex-conjugate pair of multipliers
crosses the unit circle, and the locked response becomes quasi-periodic. The
discrete-time analogue of a Hopf bifurcation.

*Reference:* [Kuznetsov], ch. 4.

### Subharmonic, harmonic

A response component at a *lower* frequency than the drive (subharmonic —
e.g. one third of it) or at a multiple of it (harmonic). `speculation.md`
uses harmonic content in an FFT to tell the models apart.

### Relaxation oscillator

An oscillator whose cycle is not sinusoidal but consists of slow build-ups
punctuated by fast jumps — a dripping tap, a neon flasher, van der Pol at
large $`\mu`$. Strongly nonlinear, and very strongly damped transverse to its
cycle.

*Reference:* [van der Pol]; [Strogatz], ch. 7.

---

## 6. Named systems

### Van der Pol oscillator

```math
\ddot{x} - \mu\left(1-x^2\right)\dot{x} + \omega_n^2 x = 0
```

The standard self-exciting oscillator: damping is *negative* for
$`\lvert x\rvert \lt 1`$ and positive outside, so small oscillations grow and
large ones decay, and every start converges on one limit cycle of amplitude
about 2. Introduced for vacuum-tube circuits and now the reference example
for self-sustained oscillation. This repository uses it as its control:
the prototypes are piecewise approximations to it, and forced Van der Pol is
the case where chaos is known to occur.

*Reference:* [van der Pol]; [Strogatz], ch. 7. For chaos in the forced case:
[Cartwright & Littlewood], [Levinson], [Ueda & Akamatsu].

### Duffing oscillator

$`\ddot{x} + \delta\dot{x} + \alpha x + \beta x^3 = \gamma\cos\omega t`$ — a
nonlinear *stiffness* rather than nonlinear damping. Its signature is that
the resonant frequency shifts with amplitude. With $`\alpha \lt 0`$ it is
the **double well**: a saddle at the origin between two stable wells,
the model of a buckled beam, and under a drive one of the classic routes
to chaos. `DUFFING.md` builds the piecewise version, with the cubic
replaced by a band of negative stiffness between two linear wells, and
uses the same double well for the pendulum seen from its inverted
position.

*Reference:* [Nayfeh & Mook], ch. 4; [Strogatz], ch. 12; [Guckenheimer &
Holmes], ch. 2 for the forced double well; [Baker & Gollub] for the
driven pendulum.

### Maxwell's governor

The flyball governor that regulates a steam engine's speed, analysed by
Maxwell in the paper that founded control theory. `EXAMPLES.md` fits several
prototypes to governor variants — a one-sided brake, an overspeed trip, a
deadband.

*Reference:* [Maxwell].

### LC tank, negative resistance

An inductor and capacitor together oscillate at
$`\omega_n = 1/\sqrt{LC}`$ — the **tank** circuit. Real components have
resistance, which damps it; an active device (transistor, tunnel diode)
arranged to supply energy behaves like a **negative resistance**, cancelling
the loss and sustaining the oscillation. Saturation of that device is what
limits the amplitude, and it is the switch modelled in `EXAMPLES.md`.

*Reference:* [Andronov, Vitt & Khaikin]; [Ueda & Akamatsu].

---

## 7. Numerical and computational vocabulary

### Stiff equations, LSODA

An equation is **stiff** when it contains processes on wildly different
timescales, so an ordinary method must take absurdly small steps for
stability rather than accuracy. **LSODA** is a solver that detects this and
switches method automatically. Used here because the strongly damped zones
are stiff.

*Reference:* [Petzold]; [Hairer & Wanner].

### Absolute and relative tolerance (`atol`, `rtol`)

The accuracy a numerical solver is asked to achieve per step — `rtol` as a
fraction of the value, `atol` as an absolute floor. Requesting more accuracy
than the arithmetic can deliver does not help, which is a recurring theme in
this repository's numerical notes.

### Event detection

Stopping an integration exactly when some condition is met — here, when the
trajectory reaches a switching boundary. Requires solving a scalar equation
for the crossing time.

*Reference:* [Hairer, Nørsett & Wanner], ch. II.

### Transcendental equation

An equation not solvable by algebra — one mixing polynomials with
exponentials or trigonometric functions, like
$`\alpha\cos t + \beta\sin t = \gamma`$ with a decaying factor. It must be
solved numerically. Every boundary-crossing time here is of this kind, which
is the reason the "exact" maps are not closed-form.

### Root finding, bracketing, Brent's method

Finding where a function crosses zero. **Bracketing** means finding two
points where it has opposite signs, guaranteeing a root between them;
**Brent's method** then narrows the bracket quickly and reliably.

*Reference:* [Press et al.], ch. 9.

### Affine map, augmented state

A **linear** map is $`y \mapsto My`$; an **affine** one is $`y \mapsto My + b`$,
with a constant offset. An affine map can be made linear by carrying an extra
component fixed at 1 — the **augmented state** — so the offset becomes part
of the matrix. `MAPS.md` uses this to turn each arc into a matrix.

### Jacobian

The matrix of partial derivatives — how each output responds to each input.
For a map it tells you how a small perturbation is transformed, so its
eigenvalues give stability.

### Shadowing

A deep and reassuring result about chaotic systems: although a computed
trajectory diverges from the true one starting at the same point, there
usually exists a *different* true trajectory that stays close to the computed
one throughout. So a numerically computed chaotic orbit is not "the" orbit,
but it is *an* orbit — which is why a computed attractor is trustworthy as a
**shape** even when no individual point on it is right.

*Reference:* [Hammel, Yorke & Grebogi].

### Fast Fourier Transform (FFT)

An efficient algorithm for decomposing a signal into its constituent
frequencies. `speculation.md` proposes using the pattern of harmonics in a
measured signal to decide which prototype fits.

*Reference:* [Press et al.], ch. 12.

---

## References

- **[Andronov, Vitt & Khaikin]** A. A. Andronov, A. A. Vitt & S. E. Khaikin,
  *Theory of Oscillators*, Pergamon Press, 1966. The classic treatment of
  self-excited and relaxation oscillators, including piecewise models.
- **[Arnol'd]** V. I. Arnol'd, "Small denominators I: On the mapping of a
  circle into itself", *American Mathematical Society Translations* 46
  (1965), 213–284. Origin of the tongues.
- **[Åström & Murray]** K. J. Åström & R. M. Murray, *Feedback Systems: An
  Introduction for Scientists and Engineers*, 2nd ed., Princeton University
  Press, 2021. Freely available online; the friendliest starting point for
  the control-theory vocabulary.
- **[Baker & Gollub]** G. L. Baker & J. P. Gollub, *Chaotic Dynamics: An
  Introduction*, 2nd ed., Cambridge University Press, 1996. The driven
  damped pendulum as the worked example, with the drive parameters
  `DUFFING.md` borrows.
- **[Benettin et al.]** G. Benettin, L. Galgani, A. Giorgilli &
  J.-M. Strelcyn, "Lyapunov characteristic exponents for smooth dynamical
  systems and for Hamiltonian systems; a method for computing all of them",
  *Meccanica* 15 (1980), parts I and II. The standard algorithm.
- **[Cartwright & Littlewood]** M. L. Cartwright & J. E. Littlewood, "On
  non-linear differential equations of the second order I", *Journal of the
  London Mathematical Society* 20 (1945), 180–189. The first identification
  of chaotic behaviour in the forced van der Pol equation.
- **[di Bernardo et al.]** M. di Bernardo, C. J. Budd, A. R. Champneys &
  P. Kowalczyk, *Piecewise-smooth Dynamical Systems: Theory and
  Applications*, Springer, 2008. The reference for everything in section 3.
- **[Filippov]** A. F. Filippov, *Differential Equations with Discontinuous
  Righthand Sides*, Kluwer, 1988.
- **[Franklin, Powell & Emami-Naeini]** G. F. Franklin, J. D. Powell &
  A. Emami-Naeini, *Feedback Control of Dynamic Systems*, Pearson. For the
  s-plane and the z-plane.
- **[Guckenheimer & Holmes]** J. Guckenheimer & P. Holmes, *Nonlinear
  Oscillations, Dynamical Systems, and Bifurcations of Vector Fields*,
  Springer, 1983. Graduate level; the standard reference for Floquet theory
  and Poincaré maps.
- **[Hairer, Nørsett & Wanner]** E. Hairer, S. P. Nørsett & G. Wanner,
  *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed.,
  Springer, 1993.
- **[Hairer & Wanner]** E. Hairer & G. Wanner, *Solving Ordinary Differential
  Equations II: Stiff and Differential-Algebraic Problems*, 2nd ed.,
  Springer, 1996.
- **[Hammel, Yorke & Grebogi]** S. M. Hammel, J. A. Yorke & C. Grebogi, "Do
  numerical orbits of chaotic dynamical processes represent true orbits?",
  *Journal of Complexity* 3 (1987), 136–145.
- **[Kuznetsov]** Y. A. Kuznetsov, *Elements of Applied Bifurcation Theory*,
  3rd ed., Springer, 2004. The catalogue of bifurcations.
- **[Leine & Nijmeijer]** R. I. Leine & H. Nijmeijer, *Dynamics and
  Bifurcations of Non-Smooth Mechanical Systems*, Springer, 2004. Clearest
  treatment of the saltation matrix.
- **[Levinson]** N. Levinson, "A second order differential equation with
  singular solutions", *Annals of Mathematics* 50 (1949), 127–153.
- **[Maxwell]** J. C. Maxwell, "On governors", *Proceedings of the Royal
  Society of London* 16 (1868), 270–283.
- **[Nayfeh & Mook]** A. H. Nayfeh & D. T. Mook, *Nonlinear Oscillations*,
  Wiley, 1979. For averaging and the Duffing equation.
- **[Ogata]** K. Ogata, *Modern Control Engineering*, Pearson. Second-order
  response, damping ratio, logarithmic decrement.
- **[Ott]** E. Ott, *Chaos in Dynamical Systems*, 2nd ed., Cambridge
  University Press, 2002. Attractors, Lyapunov exponents, fractal dimension.
- **[Petzold]** L. Petzold, "Automatic selection of methods for solving stiff
  and nonstiff systems of ordinary differential equations", *SIAM Journal on
  Scientific and Statistical Computing* 4 (1983), 136–148. The LSODA method.
- **[Pikovsky, Rosenblum & Kurths]** A. Pikovsky, M. Rosenblum & J. Kurths,
  *Synchronization: A Universal Concept in Nonlinear Sciences*, Cambridge
  University Press, 2001. Entrainment, tongues, rotation numbers.
- **[Press et al.]** W. H. Press, S. A. Teukolsky, W. T. Vetterling &
  B. P. Flannery, *Numerical Recipes*, 3rd ed., Cambridge University Press,
  2007. Root finding and the FFT.
- **[Strogatz]** S. H. Strogatz, *Nonlinear Dynamics and Chaos*, 2nd ed.,
  CRC Press, 2015. **Start here.** Written for exactly this level, and covers
  most of sections 2, 4 and 5.
- **[Ueda & Akamatsu]** Y. Ueda & N. Akamatsu, "Chaotically transitional
  phenomena in the forced negative-resistance oscillator", *IEEE Transactions
  on Circuits and Systems* 28 (1981), 217–224. Source of the drive parameters
  used for the chaotic case in this repository.
- **[van der Pol]** B. van der Pol, "On relaxation-oscillations", *The London,
  Edinburgh and Dublin Philosophical Magazine* 2 (1926), 978–992.

### If you are reading one thing

[Strogatz] covers most of this glossary at exactly this level and is the
standard undergraduate introduction. Follow it with [di Bernardo et al.] for
the piecewise material, which Strogatz does not treat.
