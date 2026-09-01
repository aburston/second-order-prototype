# Speculation: from these prototypes to a field method

Forward-looking notes, kept separate from `README.md` because most of what
follows is **not** established. Three labels are used throughout:

- **Established** — proved or verified numerically, and written up in the README.
- **Implied** — follows from established results but has not been worked out
  or checked.
- **Speculative** — a direction, with the obstacles named.

## The target

An engineer measures a vibration, fits one of these prototypes quickly,
reads off whether to expect a limit cycle or chaos, and knows which
parameter to move — with a controller or otherwise — to get the behaviour
they want. The analogy is PID tuning rules: a caricature model, but one
whose parameters map onto knobs, with charts that make the mapping usable
without solving anything.

## What the existing results already contribute

**Established.** The four prototypes have an unusually clean parameter
separation, which is the part that makes a chart plausible at all:

| parameter | what it sets | what it does not touch |
| --- | --- | --- |
| $`\omega_n`$ | the timescale, nothing else | shape, stability |
| boundary $`v_0`$ or $`x_0`$ | the amplitude, exactly proportionally | period, stability, shape |
| $`\zeta_{+}, \zeta_{-}$ | stability, period, orbit shape | — |

So amplitude and frequency are independently adjustable, and only two
numbers govern behaviour. **Established:** the existence boundaries are
sharp and simple — $`\zeta_{-} \lt 0 \lt \bar{\zeta}`$ for the asymmetric
models, $`\zeta_{-} \lt 0 \lt \zeta_{+}`$ for the symmetric ones.

## Identification: which measurements pin down which parameters

**Implied, not yet worked out.** Two measurable quantities are
*dimensionless invariants* — they depend on $`\zeta_{+}, \zeta_{-}`$ alone,
and not on $`\omega_n`$ or the boundary:

1. **The Floquet multiplier** $`\mu = e^{2\Lambda}`$, observable as the rate
   at which a disturbed system settles back onto its cycle. $`\Lambda`$ is a
   dwell-weighted sum of pole real parts, and dwell times scale as
   $`1/\omega_n`$, so $`\Lambda`$ carries no units.
2. **Harmonic ratios** of the steady waveform. Amplitude scales out exactly
   and time scales with $`\omega_n`$, so $`h_2/h_1`$, $`h_3/h_1`$ and so on
   are functions of the two damping ratios only.

That suggests an identification order which needs no prior knowledge of the
system's stiffness or mass:

```
harmonics + settling rate  ->  zeta+, zeta-
period                     ->  wn
amplitude                  ->  boundary
```

Two dimensionless measurements for two dimensionless unknowns, then two
scale factors read off directly. **This has not been tested**: whether the
map from $`(\zeta_{+}, \zeta_{-})`$ to those observables is injective, and
how it behaves under noise, is unknown.

## Model selection from an FFT

**Verified here, not yet in the README.** Harmonic content separates the
four models, at $`\zeta_{+} = 0.3`$, $`\zeta_{-} = -0.1`$:

| model | $`h_2/h_1`$ | $`h_3/h_1`$ | $`h_4/h_1`$ |
| --- | --- | --- | --- |
| asymmetric, switch on $`\dot{x}`$ | 0.0440 | 0.0066 | 0.0006 |
| asymmetric, switch on $`x`$ | 0.0880 | 0.0199 | 0.0024 |
| symmetric, switch on $`\dot{x}`$ | 0.0000 | 0.0124 | 0.0000 |
| symmetric, switch on $`x`$ | 0.0000 | 0.0373 | 0.0000 |

Two things fall out, both usable at a screen:

**Even harmonics decide symmetric against asymmetric.** The symmetric fields
are odd, so their cycles carry odd harmonics only — exactly zero at
$`h_2`$ and $`h_4`$. A visible second harmonic rules the symmetric models
out. This is rigorous, not empirical: odd symmetry forbids even harmonics.

**The $`n`$-th harmonic ratio doubles, triples, quadruples between the
pairs.** $`0.0440 \to 0.0880`$ is exactly $`\times 2`$, and
$`0.0066 \to 0.0199`$ and $`0.0124 \to 0.0373`$ are exactly $`\times 3`$.
That is the differentiation relation showing up in the spectrum: if one
model's displacement is the other's velocity, its $`n`$-th harmonic is
scaled by $`n\omega`$. So the velocity- and displacement-switched pair,
which are *identical* in period, dwell times and multiplier, are cleanly
separable in the spectrum after all — provided you know whether the sensor
measured displacement or velocity. Confusing the two swaps the model.

## Chaos needs a third state

**Established, by Poincaré–Bendixson.** None of the four can be chaotic.
They are planar and autonomous, so the only attractors available are
equilibria and limit cycles. Any chaotic behaviour must come from
somewhere else:

- periodic forcing, making the system effectively three dimensional;
- a third state — an actuator lag, a compliant mount, a thermal drift;
- coupling two prototypes together.

**Speculative, and the natural next step.** Forcing is the cheapest of the
three and the closest to the field case, where a machine runs at a drive
frequency. A forced piecewise-linear oscillator entrains over a range of
drive frequencies and amplitudes — the Arnold tongues — and the classical
route to chaos is tongue overlap. What sets tongue width is precisely what
is already computed exactly here: the unforced cycle's frequency, and how
strongly it attracts, which is $`\mu`$. So the existing results are the
inputs to that analysis rather than a detour from it.

The rule of thumb that would come out, if it works: **drive frequency
relative to the natural cycle frequency, and drive amplitude relative to
the boundary parameter** are the two axes on which locking, quasi-periodic
and chaotic behaviour separate. That chart is the deliverable the target
needs, and it is one forcing term away.

## What would have to be true, and what probably is not

Worth stating plainly, since the ambition invites over-claiming.

- **Two states is often not enough.** Real vibration frequently needs more.
  The prototypes are a caricature in the same way a first order plus dead
  time model is for PID — useful because it is crude, not despite it.
- **Piecewise-constant damping is not physical.** Real nonlinearity is
  smooth. The piecewise version is chosen for exact solvability. Whether
  the fitted parameters mean anything mechanically, or are only descriptive,
  is untested. A comparison against a genuine Van der Pol would settle how
  much the sharp switch distorts things.
- **Identification assumes a visible limit cycle and a measurable
  transient.** A system sitting quietly on its cycle gives the period and
  amplitude but not $`\mu`$; getting that needs a disturbance and a
  recording of the settling.
- **Noise is unaddressed.** Harmonic ratios at the $`10^{-3}`$ level, as in
  the $`h_4`$ column above, will not survive a real measurement. Which
  observables are robust is an open question and probably the deciding one
  for whether any of this is practical.

## Suggested order of work

1. **Forcing.** Add a drive term to the symmetric models and map the
   locking regions against drive frequency and amplitude. This is where
   chaos first becomes possible and where the target chart lives.
2. **Identifiability.** Check whether $`(\zeta_{+}, \zeta_{-})`$ is
   recoverable from harmonics plus settling rate, and how that degrades
   with noise.
3. **Smooth comparison.** Fit the piecewise prototypes to Van der Pol data
   and see what the fitted parameters do.
4. **Only then, the chart.** A single figure with the parameter plane, the
   existence boundaries, and the forced behaviour regions marked.
