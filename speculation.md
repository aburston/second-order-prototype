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
scale factors read off directly.

### Tested, and it half works

**$`\zeta_{+}`$ is recoverable. $`\zeta_{-}`$ saturates.** Both invariants
stop responding to $`\zeta_{-}`$ once it passes roughly $`-0.3`$, for the
symmetric velocity model:

| $`\zeta_{-}`$ | $`-0.05`$ | $`-0.1`$ | $`-0.2`$ | $`-0.3`$ | $`-0.45`$ | $`-0.8`$ |
| --- | --- | --- | --- | --- | --- | --- |
| $`h_3/h_1`$ at $`\zeta_{+}=0.05`$ | 0.00328 | 0.00378 | 0.00403 | 0.00410 | 0.00414 | 0.00416 |
| multiplier at $`\zeta_{+}=0.05`$ | 0.7439 | 0.7358 | 0.7321 | 0.7311 | 0.7306 | 0.7303 |
| $`h_3/h_1`$ at $`\zeta_{+}=0.8`$ | 0.00995 | 0.01738 | 0.02811 | 0.03558 | 0.04335 | 0.05389 |
| multiplier at $`\zeta_{+}=0.8`$ | 0.0482 | 0.0261 | 0.0133 | 0.0087 | 0.0056 | 0.0030 |

The reason is simple, and it is a property of the system rather than of the
method. As $`\zeta_{-}`$ becomes more negative the cycle grows, so it
spends proportionally less of each revolution inside the destabilising
region. Past a point the orbit stops noticing how negative $`\zeta_{-}`$
is, and no measurement of the orbit can recover it.

What that means in practice:

- **At larger $`\zeta_{+}`$ the multiplier stays useful.** At
  $`\zeta_{+} = 0.8`$ it moves by a factor of sixteen across the same
  $`\zeta_{-}`$ range where it barely moves at $`\zeta_{+} = 0.05`$. So a
  strongly damped outer region makes the inner one identifiable.
- **At small $`\zeta_{+}`$, report $`\zeta_{-}`$ as a bound, not a value.**
  Anything past about $`-0.3`$ fits the data equally well.
- **The higher harmonics do not rescue it.** $`h_5/h_1`$ is an order of
  magnitude smaller than $`h_3/h_1`$ and not monotonic in $`\zeta_{-}`$;
  $`h_7/h_1`$ sits at $`10^{-4}`$, below any realistic noise floor.

This is good news for the engineering goal rather than bad. A parameter
that the data cannot determine is also a parameter the behaviour does not
depend on, over that range — so a fit that pins $`\zeta_{+}`$ and bounds
$`\zeta_{-}`$ is enough to predict what the system will do.

Still untested: behaviour under noise, and whether the same saturation
appears in the asymmetric models.

## First check the class, before fitting anything

**Established, and the cheapest test available.** All four prototypes have
**amplitude-independent frequency**. The period is fixed by $`\omega_n`$ and
the two damping ratios; the boundary parameter sets amplitude and nothing
else. Integrating at boundary values of $`0.25, 1, 4, 16`$ returns the same
period to nine decimal places.

A stiffness nonlinearity does the opposite. Duffing and its relatives put
the nonlinearity in the restoring force, so their frequency moves with
amplitude — that is where their jump and hysteresis behaviour comes from.

So excite the system at two different amplitudes and watch the frequency:

| observation | conclusion |
| --- | --- |
| frequency unchanged | damping nonlinearity — one of these four, fitting applies |
| frequency shifts with amplitude | stiffness nonlinearity — none of these four will fit, at any parameters |

This is a better class test than the harmonic signature, because it is one
robust frequency measurement rather than a ratio at the $`10^{-3}`$ level.
Use it first; use the harmonics below only to choose *which* of the four.

**Speculative.** The counterpart family, when it is wanted, is the obvious
one: switch the **stiffness** at a displacement boundary instead of the
damping — the bilinear or clearance oscillator. The same solvable-arcs
machinery would apply, and it would give amplitude-dependent frequency by
construction. That is a separate build, and none of the results here carry
over to it.

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

## Where to look next: other disciplines that have solved parts of this

**Speculative throughout, and of a different kind to the rest of this file.**
Everything above is a direction with an experiment attached. This section is
a reading list: the claim is only that these fields appear to have built,
for their own reasons, the machinery that the work here kept running into.
None of it has been checked against the sources from inside this repository,
and the summaries below are from memory rather than from the papers. Treat
each entry as a lead to verify, not as a result.

The context for asking is that the first three items of the order of work
above are now done, and two routes onward have closed rather than opened:

- the discrete-map route, in `MAPS.md`, because the zone sequence a
  trajectory takes cannot be known ahead of stepping it, and because
  grazing makes the map's derivative undefined on a dense set of
  parameters;
- the closed-form route for the smooth control, in `VANDERPOL.md`, because
  one revolution of Van der Pol is an Abel equation of the second kind and
  only its expansion in $`\mu`$ integrates.

Both obstacles are old, and neither is peculiar to this repository. That is
the reason to go looking sideways.

### Power electronics: the same system, with a design language attached

A switched-mode converter is a second order piecewise-linear system with a
switching boundary in the state plane — structurally the same object as the
prototypes here, with duty cycle in the role of the boundary parameter.
Because these ship in volume, the field was forced to produce engineering
answers rather than phase portraits, and the two that look most transferable
are:

- **Border-collision bifurcation theory.** The normal form for what happens
  when a fixed point of the return map crosses a switching boundary. This is
  exactly the grazing event that made the Jacobian in `MAPS.md` undefined,
  and the treatment is the opposite of the one attempted there: rather than
  chain arcs globally, work with the *local* map at the boundary, which is
  piecewise-linear with two branches. Its outcomes — period doubling, direct
  transition to chaos, robust chaos — are then classified as regions in the
  plane of the two branch slopes. Two numbers and a chart of regions is very
  close to the deliverable this file has been calling the target.
- **Fast-scale against slow-scale instability** as vocabulary, with slope
  compensation as the standard corrective term. That is the PID analogy
  done properly: measure, identify the region, apply the known fix.

### Impact mechanics: what the grazing singularity actually is

Vibro-impact dynamics — gear rattle, rotor-stator rub, machining chatter,
loose components — met the "cannot guarantee which side the trajectory
leaves by" problem decades ago. The result to look up is **Nordmark's**:
near grazing the return map is not smooth but carries a **square root**
singularity, so the local map behaves like $`\sqrt{\cdot}`$ rather than
linearly. If that holds for these prototypes too, it says the derivative in
`MAPS.md` did not merely fail to converge — there was no derivative there to
find, and a linear difference equation was the wrong object rather than a
badly conditioned one. The map itself remains tractable; it is only
differentiability that is lost.

This is the most direct diagnosis of the closed route, and the cheapest
thing on this list to check: take one of the staircase models, put an orbit
onto a grazing tangency by tuning the boundary, and see whether the
amplitude response goes as the square root of the distance past it.

### Hybrid systems verification: propagate sets, not trajectories

Reachability analysis for hybrid automata exists precisely because a
verifier cannot know which discrete mode a trajectory will take. Its answer
to that is not to choose: propagate a **set** of states through every branch
at once, represented as zonotopes or support functions, and return a
guaranteed over-approximation of where the state can go. The output is
"this region cannot be left" or "this region is reachable", with a proof,
instead of a single trajectory that is wrong as soon as the zone sequence
differs.

Given that the goal here is stability *boundaries* rather than particular
trajectories, this may be the right shape of tool, and the tooling is mature
and open source (CORA, SpaceEx, Flow\*, JuliaReach). The same mathematics
underlies verification of ReLU neural networks, which is piecewise-linear
with hyperplane boundaries for the same reason, and which is why the
computational side has had money spent on it.

### Threshold time series: the statistics of finding the boundary

The most unwitting of the parallels. **Threshold autoregressive** models
(Tong's TAR and SETAR) and **Markov-switching** models are piecewise-linear
difference equations with a switching threshold — very nearly the object
`MAPS.md` set out to build, arrived at independently in econometrics. What
that literature has and dynamics does not is the **estimation theory**:

- how to locate the threshold from a noisy record;
- how to *test* whether a threshold exists at all, against the null of a
  single linear regime;
- how to put a confidence interval on the threshold once found.

For an engineer holding a measured vibration and no prior idea where the
switch sits, that is the missing half of the identification section above.
It also speaks directly to the open question at the end of this file, which
is whether any of these observables survive noise: threshold estimation is
built for noisy data from the start, where the fitting sketched here assumed
clean records.

### The measurement front end: two things that already exist

- **Relay autotuning.** Drive the plant with a relay, let it settle into a
  limit cycle, read amplitude and frequency, back out the model, set the
  controller. That is literally the PID workflow named as the target at the
  top of this file, already productised, and its theory — describing
  functions for relay feedback, Tsypkin's locus — gives exact limit cycle
  existence and stability conditions for switched systems.
- **Hilbert transform instantaneous modal analysis.** From a single free
  decay record, Feldman's method extracts instantaneous frequency and
  **instantaneous damping as a function of amplitude**. That curve is the
  $`\zeta(A)`$ staircase, measured directly. If it works as advertised it
  would tell an engineer which prototype they have from one hammer test,
  with no fitting search at all — and it measures the same amplitude
  dependence that the class test above uses, so the two cross-check.

### What each field supplies

| obstacle met here | field | what it appears to supply |
| --- | --- | --- |
| grazing makes the map non-differentiable | impact mechanics | the square root normal form; the derivative was never there |
| zone sequence unknown in advance | hybrid systems verification | set propagation over all branches, with guarantees |
| want regions, not trajectories | power electronics | border-collision classification in the branch-slope plane |
| boundary location unknown from data | threshold time series | estimation, testing and confidence intervals for the threshold |
| noise robustness untested | threshold time series | the same, since it assumes noise from the start |
| fitting from a field measurement | relay autotuning, Hilbert methods | amplitude, frequency and $`\zeta(A)`$ from one test |

### The order I would take them in

1. **Nordmark's square root**, because it is a diagnosis of a specific
   failure already recorded here, and can be tested against the existing
   staircase code in an afternoon.
2. **Feldman's instantaneous damping** on a synthetic decay from each
   prototype, because it turns identification from a fit into a
   measurement, and the prototypes provide their own ground truth.
3. **Border-collision normal form**, once the boundary behaviour is
   understood, because it is the format the final chart should be in.
4. **Reachability tooling**, only if the chart needs guarantees rather than
   observations — it is the largest investment of the four.

### References

These are pointers to locate the work, not sources consulted here.

1. M. di Bernardo, C. Budd, A. Champneys, P. Kowalczyk, *Piecewise-smooth
   Dynamical Systems: Theory and Applications*, Springer, 2008. The
   canonical text for Filippov systems, sliding, grazing and the saltation
   matrix.
2. S. Banerjee, G. Verghese (eds), *Nonlinear Phenomena in Power
   Electronics*, IEEE Press, 2001. Border collision and robust chaos in
   switched converters.
3. A. Nordmark, "Non-periodic motion caused by grazing incidence in an
   impact oscillator", *Journal of Sound and Vibration* 145(2), 1991.
4. H. Nusse, J. Yorke, "Border-collision bifurcations including period two
   to period three for piecewise smooth systems", *Physica D* 57, 1992.
5. M. Wiercigroch, B. de Kraker (eds), *Applied Nonlinear Dynamics and Chaos
   of Mechanical Systems with Discontinuities*, World Scientific, 2000.
6. M. Althoff, "An introduction to CORA", ARCH 2015; G. Frehse et al.,
   "SpaceEx: scalable verification of hybrid systems", CAV 2011; X. Chen,
   E. Abraham, S. Sankaranarayanan, "Flow\*: an analyzer for non-linear
   hybrid systems", CAV 2013.
7. H. Tong, *Non-linear Time Series: A Dynamical System Approach*, Oxford,
   1990. SETAR models.
8. B. Hansen, "Sample splitting and threshold estimation", *Econometrica*
   68(3), 2000. Confidence intervals for an estimated threshold.
9. J. Hamilton, "A new approach to the economic analysis of nonstationary
   time series and the business cycle", *Econometrica* 57(2), 1989.
   Markov switching.
10. K. Astrom, T. Hagglund, "Automatic tuning of simple regulators with
    specifications on phase and amplitude margins", *Automatica* 20(5),
    1984. The relay autotuner.
11. Y. Tsypkin, *Relay Control Systems*, Cambridge, 1984.
12. M. Feldman, "Non-linear system vibration analysis using Hilbert
    transform — I. Free vibration analysis method FREEVIB", *Mechanical
    Systems and Signal Processing* 8(2), 1994; and "II. Forced vibration
    analysis method FORCEVIB", 8(3), 1994.
13. S. Masri, T. Caughey, "A nonparametric identification technique for
    nonlinear dynamic problems", *Journal of Applied Mechanics* 46, 1979.
    The restoring force surface method.
14. D. Barton, "Control-based continuation: bifurcation and stability
    analysis for physical experiments", *Mechanical Systems and Signal
    Processing* 84, 2017. Tracing unstable branches on a rig.
15. L. Glass, M. Mackey, *From Clocks to Chaos: The Rhythms of Life*,
    Princeton, 1988. Circle maps, Arnold tongues and entrainment measured
    in a physical system.
