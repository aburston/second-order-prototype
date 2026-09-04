# Three level prototype — data sheet

Behavioural model of a self-excited oscillator of the Van der Pol class.
Second order, damping ratio switched on displacement at two thresholds.
Piecewise linear and exact by pieces. Under sinusoidal drive it predicts
lock, beat or chaos, and the drive frequencies at which those change. It
also covers hard excitation, free response only.

**Status.** Proven against Van der Pol at $`\mu = 1`$, 2 and 5 across a
drive grid; fitted at eleven values of $`\mu`$ from 0.1 to 5 and checked at
8; closed form parameter laws. Fitted to exact integrations, never to
measured data.

Formulas below are stated, not derived. Proofs and evidence: §7.

---

## 1. Model

```math
\ddot{x} + 2\zeta(x)\,\omega_n\dot{x} + \omega_n^2 x = A\cos\Omega t,
\qquad
\zeta(x) = \begin{cases} \zeta_{2} & |x| \gt b \\ \zeta_{1} & a \lt |x| \lt b \\ \zeta_{0} & |x| \lt a \end{cases}
\qquad 0 \lt a \lt b
```

| symbol | meaning | units |
| --- | --- | --- |
| $`\omega_n`$ | natural frequency, $`\sqrt{k/m}`$; sets every time | rad/s |
| $`a, b`$ | inner and outer thresholds; set the amplitude scale | units of $`x`$ |
| $`\zeta_{0}`$ | core damping ratio, $`\lvert x \rvert \lt a`$; negative for self-excitation | — |
| $`\zeta_{1}`$ | band damping ratio | — |
| $`\zeta_{2}`$ | outer damping ratio, $`\lvert x \rvert \gt b`$ | — |
| $`A, \Omega`$ | drive acceleration amplitude and angular frequency | $`x`$/s², rad/s |

Behaviour depends only on the three ratios and on $`a/b`$. $`\omega_n`$ is a
timescale, $`b`$ an amplitude scale; §6 moves both. From
$`m\ddot{x} + c(x)\dot{x} + kx = F\cos\Omega t`$: $`c_{k} = 2\zeta_{k} m \omega_n`$,
$`A = F/m`$.

**Drive coordinates.** A drive enters through two numbers only:

```math
r = \frac{\Omega}{\omega_{lc}}, \qquad \text{strength} = \frac{A}{\omega_n^{2} b}
```

$`\omega_{lc} = 2\pi/T`$ is the **free cycle** frequency, not $`\omega_n`$.
All drive ratios in this sheet are against Van der Pol's $`\omega_{lc}`$ at
the same $`\mu`$.

**Reference units.** Every characteristic is quoted at $`\omega_n = 1`$,
$`b \approx 2`$. Fit, sweep and classify there; scale afterwards (§6, L12).

---

## 2. Design data

### 2.1 Parameter laws

For a Van der Pol class oscillator of relaxation parameter $`\mu`$, valid
$`0.1 \le \mu \le 8`$:

```math
\zeta_{0} = -0.365\,\mu^{0.93}, \qquad
\zeta_{1} = 0.824\,\mu^{0.96}, \qquad
\zeta_{2} = 3.29\,\mu^{1.17}, \qquad
a = 1.20\lambda, \quad b = 2.13\lambda
```

$`\lambda = R/2`$, with $`R`$ the measured free amplitude; equivalently
$`\lambda`$ is the displacement at which the damping changes sign. Edge
tolerance $`a \pm 0.05`$, $`b \pm 0.09`$.

### 2.2 Design table

$`T\omega_n`$ is the free cycle period in reference units: **target** is Van
der Pol's, which identifies $`\mu`$; **model** is this prototype's, which
sets $`\omega_n`$ from a measured period (§3.2 step 4).

| $`\mu`$ | $`T\omega_n`$ target | $`T\omega_n`$ model | $`\zeta_{0}`$ | $`\zeta_{1}`$ | $`\zeta_{2}`$ |
| --- | --- | --- | --- | --- | --- |
| 0.1 | 6.287 | 6.288 | $`-0.043`$ | 0.090 | 0.22 |
| 0.2 | 6.299 | 6.305 | $`-0.082`$ | 0.176 | 0.50 |
| 0.3 | 6.318 | 6.316 | $`-0.119`$ | 0.259 | 0.80 |
| 0.5 | 6.381 | 6.393 | $`-0.192`$ | 0.424 | 1.46 |
| 0.7 | 6.473 | 6.440 | $`-0.262`$ | 0.585 | 2.17 |
| 1 | 6.663 | 6.701 | $`-0.365`$ | 0.824 | 3.29 |
| 1.5 | 7.096 | 7.138 | $`-0.532`$ | 1.216 | 5.29 |
| 2 | 7.630 | 7.756 | $`-0.695`$ | 1.603 | 7.40 |
| 3 | 8.859 | 9.252 | $`-1.014`$ | 2.366 | 11.9 |
| 4 | 10.204 | 10.469 | $`-1.325`$ | 3.118 | 16.7 |
| 5 | 11.612 | 11.961 | $`-1.631`$ | 3.863 | 21.6 |
| 8 | 16.0 | 16.686 | $`-2.524`$ | 6.066 | 37.5 |

$`\zeta_{2}`$ is weakly determined: at $`\mu = 5`$ values of 15 and 22 fit
the lock plateaus equally well, because the driven orbit barely enters the
outer zone. Take it from the law and leave it.

### 2.3 Reference builds

At $`\omega_n = 1`$, in `staircase.py`.

| build | $`\zeta_{0}, \zeta_{1}, \zeta_{2}`$ | $`a, b`$ | free amplitude | free $`T\omega_n`$ |
| --- | --- | --- | --- | --- |
| $`\mu = 5`$ proof fit (`THREE_FITTED`) | $`-1.735,\ 3.836,\ 15.05`$ | 1.075, 1.981 | 1.996 | 12.44 |
| $`\mu = 1`$ proof fit (`THREE_FITTED_MU1`) | $`-0.358,\ 0.865,\ 3.573`$ | 1.160, 1.984 | 1.981 | 6.670 |
| bistable case (`BISTABLE_*`) | $`0.15,\ -0.25,\ 0.40`$ | 0.6, 1.6 | 1.168 unstable / 2.253 stable | — |

### 2.4 Operating modes

$`\langle\zeta\rangle(R)`$ is the cycle averaged damping ratio at amplitude
$`R`$; it runs monotonically from $`\zeta_{0}`$ to $`\zeta_{2}`$ and a limit
cycle sits at each zero crossing.

| sign pattern | free response | under drive |
| --- | --- | --- |
| $`\zeta_{0} \gt 0`$, $`\langle\zeta\rangle`$ never negative | decays to rest | linear, no locking |
| $`\zeta_{0} \lt 0 \lt \zeta_{1} \le \zeta_{2}`$, all below ~1 | one near sinusoidal cycle | 1:1 and 3:1 tongues, tori between, one narrow period doubled chaotic band |
| $`\zeta_{0} \lt 0`$, $`\zeta_{2} \gg 1`$ | one cycle, crawl and jump | 1:1, 3:1, 5:1, 7:1 period adding, chaos at the transitions |
| $`\zeta_{0} \gt 0`$, $`\zeta_{1} \lt 0`$, $`\zeta_{2} \gt 0`$ | stable origin, unstable inner cycle, stable outer cycle | unmapped (L10) |

---

## 3. Setting the model

### 3.1 Gate — is this the right model?

Drive at two amplitudes and watch the response.

- Response **frequency** shifts with amplitude → stiffness nonlinearity.
  Stop; nothing in this repository fits it (L1).
- Nearly sinusoidal free cycle, chaos not of interest → two level model.
  Ringdown only → linear prototype.
- Lock plateaus and the transitions between them **move with drive
  amplitude** → the damping shape matters. Use this model, set from §3.2 or
  §3.3.

### 3.2 Known Van der Pol class damping — no fitting

Damping known or believed to follow $`-\varepsilon(1 - x^{2}/X^{2})`$.

| step | measure | sets |
| --- | --- | --- |
| 1 | free period $`T`$, free amplitude $`R`$ | — |
| 2 | $`\mu = \varepsilon/\omega_n`$, or read $`T\omega_n`$ off §2.2 | $`\mu`$ |
| 3 | — | $`\zeta_{0}, \zeta_{1}, \zeta_{2}`$ from §2.1 |
| 4 | — | $`\lambda = R/2`$, then $`a, b`$; $`\omega_n = (T\omega_n)_{\text{model}}/T`$ |
| check | one lock plateau at one drive strength | edges within a few per cent → done, else §3.3 |

Cost: one sweep for the check. Accuracy: plateau edges within 0.02 in
drive ratio, chaos at the same transitions — tested at $`\mu = 2.5`$ and 3.5
with no fitting at all.

### 3.3 Unknown damping law — fit to the driven response

A ringdown will not do it: nothing measured without a drive depends on the
three ratios beyond a few per cent. Two drive strengths are needed.

| step | measure | sets |
| --- | --- | --- |
| 1 | free $`T`$, free $`R`$ | held within 20% as a constraint |
| 2 | at one drive strength, the last frequency locked $`p{:}1`$ and the first locked $`(p+2){:}1`$ | all five of $`\zeta_{0}, \zeta_{1}, \zeta_{2}, a, b`$ |
| 3 | a lock plateau at a **second** drive strength | check only — edges right → the fit holds across drive |
| 4 | free $`T`$, free $`R`$ | $`\omega_n`$ and $`\lambda`$ as §3.2 step 4 |

Sweep drive frequency at fixed strength, spacing finer than a tongue edge.
`staircase.fit_bands` moves the five parameters by Nelder–Mead onto the two
measured edges. Start from §2.1 at the nearest $`\mu`$ — from $`\mu = 1.5`$
up, every campaign fit started that way landed within a coarse step before
the first evaluation. Budget: ~1 min per evaluation on four cores, 40 to 60
evaluations.

### 3.4 Hard excitation

For a system that sits quiet until knocked past a threshold, then sustains
a large oscillation. Free response only (L10).

| measure | sets |
| --- | --- |
| small signal logarithmic decrement at the origin | $`\zeta_{0} \gt 0`$ |
| smallest kick that starts the large oscillation | inner (unstable) cycle radius |
| sustained amplitude | outer (stable) cycle radius |

The two radii are the zero crossings of $`\langle\zeta\rangle(R)`$
(`staircase.cycles_predicted`, exact by `staircase.cycles_exact`). With
$`\zeta_{0}`$ measured they give two equations in four unknowns: choose
$`\zeta_{2}`$ and $`a/b`$, solve for $`\zeta_{1}`$ and the edge scale.

### 3.5 Simpler model instead

| situation | use |
| --- | --- |
| $`\mu \lesssim 0.3`$, or a nearly sinusoidal cycle | two level model (L11) |
| chaos wanted with fewest parameters, window interiors not needed | two level, $`\zeta \approx (-1.24,\ 8.33)`$, $`x_0 = 1.44`$ at $`\mu = 5`$ |
| free response in closed form, no drive | Hopf normal form |
| ringdown | linear prototype |

---

## 4. Driven response

### 4.1 Lock plateaus at strength 2.5 ($`A = 5`$, $`b \approx 2`$)

Van der Pol's edges in drive ratio; the fits reproduce all of these to
within 0.04.

| $`\mu`$ | 0.3 | 0.5 | 1 | 1.5 | 2 | 3 | 4 | 5 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:1 ends | 2.11 | 2.13 | 2.21 | 2.29 | 2.38 | 2.49 | 2.59 | 2.68 | 2.81 |
| 3:1 plateau | 2.87–3.12 | 2.77–3.21 | 2.56–3.41 | 2.47–3.61 | 2.47–3.79 | 2.53–4.11 | 2.59–4.33 | 2.68–4.48 | 2.81–4.74 |

From $`\mu = 4`$ the 1:1 and 3:1 plateaus meet directly. Below
$`\mu = 0.3`$ there is no 3:1 plateau at this drive (L11). Across the full
grid at $`\mu = 5`$ the model places the 1:1, 3:1, 5:1 and 7:1 plateaus
within one or two cells of 0.1 in drive ratio at every drive strength from
0.25 to 5, from a fit at one strength.

### 4.2 Chaos

**Rule: where the model's lock plateaus meet, chaos is coming; where they
meet cleanly, it is not.**

At strength 2.5: none up to $`\mu = 1`$; from $`\mu = 1.5`$ at the 1:1 → 3:1
transition, within 0.05 in ratio; from $`\mu = 4`$ above the 3:1 plateau,
within 0.05 to 0.1. At $`\mu = 3`$ the model gets two of Van der Pol's three
transitions and switches through a torus at the third. At strength 5
($`A = 10`$) chaos first appears between $`\mu = 2`$ and 2.125 in both
systems.

Working procedure, with a model from §3.2 or §3.3:

1. Sweep drive frequency at the working strength, in reference units,
   0.01 to 0.02 in ratio.
2. Read where the plateaus meet. Report those frequencies to within 0.05
   in drive ratio (0.1 at strength 5), and the **width** of any chaotic
   window as unknown (L4).
3. Confirm every chaotic verdict with `section.confirm_chaos`. Isolated
   narrow locks at integer ratios are suspect (L6).

---

## 5. Limits

Inside the tested envelope the behaviour is established. Outside it is
untested — which is not the same as wrong, but nothing here vouches for it.

**Tested envelope:** $`\mu`$ 0.1 to 5 fitted, 8 checked; drive strength 0 to
5 and drive ratio $`r`$ 0.5 to 8 at $`\mu = 1`$, 2 and 5, elsewhere to 2.5
and 6; any $`\omega_n`$ and amplitude scale; exact integrations only.

- **L1 — Stiffness nonlinearity is out of scope.** Response frequency
  shifting with drive amplitude rules the model out at every parameter set.
- **L2 — Beyond the envelope, nothing is known.** Above $`\mu = 8`$,
  strength 5 or ratio 8 the laws are extrapolation; at $`\mu = 8`$ they hold
  to 5%, at 10 they are a guess.
- **L3 — Plateau edges are good to a few per cent, not better.** Free
  period runs up to 7% long ($`\mu = 5`$ proof fit; campaign fits within
  4.4%), and every lock sits about 2% low in ratio.
- **L4 — Chaos: where, not how wide.** Window widths are not reproduced —
  at $`\mu = 4`$, strength 5, eight cells against Van der Pol's four. Cell
  by cell agreement is 0.2 to 0.5.
- **L5 — Outer saturation.** Past $`b`$ the model stays linear at
  $`\zeta_{2}`$ for ever. Low ratio with high strength pushes the orbit well
  past $`b`$ and onto a plateau the real system does not have. On the tested
  grid the orbit never exceeded $`\lvert x \rvert \approx 2.1`$.
- **L6 — Narrow integer locks are artefacts.** The model locks more readily
  than a smooth system: orders 4 to 11 at integer drive ratios where Van der
  Pol has tori, at every $`\mu`$ tested. Confirm before believing.
- **L7 — Multistability.** Where chaos coexists with a lock, which one is
  reported depends on the initial state. No basins have been computed under
  drive.
- **L8 — The levels are a shape, not a sampling.** The core sits 1.3 times
  and the outer 1.2 to 1.5 times the underlying damping law. Do not read a
  fitted $`\zeta_{k}`$ as physical damping at that displacement.
- **L9 — Not validated on measured data.** All fits used exact
  integrations; the classifier's thresholds are set from integrator noise
  and have never seen a noisy sweep.
- **L10 — Hard excitation is free response only.** No regime map, no chaos
  verdicts, and no recipe for the two free choices in §3.4 beyond the
  worked case.
- **L11 — Below $`\mu \approx 0.3`$ the third level does nothing.** No 3:1
  plateau at strength 2.5, and the two level model is already
  behaviourally Van der Pol. The laws still give parameters; the data
  cannot check them.
- **L12 — Analyse in reference units.** `section.py` has absolute floors
  (integration tolerance $`10^{-11}`$ in $`x`$, a Lyapunov threshold per
  unit time) set for $`\omega_n = 1`$, $`b \approx 2`$. Never run the
  classifier on a model in millimetres and kilohertz.
- **L13 — Exact, but not a predictive recurrence.** The sequence of zones a
  trajectory visits can only be found by stepping, and at a grazing the
  Jacobian is undefined. Use the pieces as an event driven integrator.

---

## 6. Scaling

To place the reference model at natural frequency $`\omega_n`$ with
$`\lambda`$ units of displacement per model unit. Exact — no re-fitting,
behaviour preserved.

```math
\zeta_{k} \to \zeta_{k}, \qquad
(a, b) \to (\lambda a,\ \lambda b), \qquad
\Omega \to \omega_n\Omega, \qquad
A \to \lambda\,\omega_n^{2} A
```

| quantity | multiply by |
| --- | --- |
| times: period, settling, transient | $`1/\omega_n`$ |
| frequencies, rates, Lyapunov exponents | $`\omega_n`$ |
| displacements: thresholds, amplitudes, basin boundary | $`\lambda`$ |
| velocities | $`\lambda\omega_n`$ |
| accelerations, drive amplitude $`A`$ | $`\lambda\omega_n^{2}`$ |
| damping $`c_{k}`$ per unit mass | $`\omega_n`$ |
| stiffness per unit mass | $`\omega_n^{2}`$ |
| $`a/b`$, $`\mu`$, $`r`$, strength, lock orders, multipliers, verdicts | 1 |

Two traps:

- The drive scales as $`\omega_n^{2}`$. Doubling the frequency at the same
  displacement scale needs **four times** the drive acceleration to stay at
  the same point of the regime map.
- $`\omega_n`$ is not the cycle frequency. A relaxation oscillator running
  at $`f`$ hertz needs $`\omega_n = (T\omega_n)_{\text{model}} f`$, well
  above $`2\pi f`$.

---

## 7. Sources

Every number above comes from one of these; that is where the derivations,
the proofs and the evidence are.

| for | see |
| --- | --- |
| definition, fits, laws, proof against Van der Pol, scaling | `THREELEVEL.md`; `staircase.py`, `campaign.py`, `scaling.py` |
| model selection across the repository, why drive not free cycle sets the ratios | `VANDERPOL.md` |
| exact arcs, crossings, Floquet multipliers, event driven stepping | `MAPS.md`; `maps.py` |
| the two level and linear prototypes | `README.md` |
| lock and chaos classifier and its thresholds | `section.py` |
| physical systems by prototype | `EXAMPLES.md` |
| terms | `GLOSSARY.md` |

Figures: `python3 figures.py`. Fitted constants live in `staircase.py`, so
nothing here needs an optimiser re-run.
