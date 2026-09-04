# Three level prototype — data sheet

**Device.** Second order oscillator with switched damping: three damping
ratios, two displacement thresholds, five parameters plus a timescale and
an amplitude scale. Piecewise linear, exact by pieces.

**Function.** Behavioural model of a self-excited oscillator of the Van
der Pol class — and of hard excitation, which no simpler prototype in this
repository has. Predicts, for a periodic drive of given strength and
frequency, whether the response locks, beats or goes chaotic, and at which
drive frequencies the transitions between those fall.

**Status.** Proven against Van der Pol at $`\mu = 1`$, 2 and 5 across a
drive grid, fitted at eleven values of $`\mu`$ from 0.1 to 5 and checked
at 8, with closed form parameter laws. Not yet fitted to measured data.

This sheet states what the model is, what it can be used for, how to set
it for each scenario, and where it stops being valid. It contains no
proofs. Every number is taken from `THREELEVEL.md`, `VANDERPOL.md`,
`MAPS.md` and the scripts they name, which is where the evidence is; a
reference like *THREELEVEL §Formulas* points at the section that
establishes the claim.

---

## 1. Definition

```math
\ddot{x} + 2\zeta(x)\,\omega_n\dot{x} + \omega_n^2 x = A\cos\Omega t,
\qquad
\zeta(x) = \begin{cases} \zeta_{2} & |x| \gt b \\ \zeta_{1} & a \lt |x| \lt b \\ \zeta_{0} & |x| \lt a \end{cases}
\qquad 0 \lt a \lt b
```

| symbol | meaning | units | carries behaviour? |
| --- | --- | --- | --- |
| $`\omega_n`$ | natural frequency of the spring; sets every time | rad/s | no — timescale only |
| $`a, b`$ | inner and outer thresholds; set the amplitude scale | units of $`x`$ | only through $`a/b`$ |
| $`\zeta_{0}`$ | damping ratio in the core, $`\lvert x \rvert \lt a`$; negative for self-excitation | — | yes |
| $`\zeta_{1}`$ | damping ratio in the band | — | yes |
| $`\zeta_{2}`$ | damping ratio outside, $`\lvert x \rvert \gt b`$ | — | yes |
| $`A, \Omega`$ | drive acceleration amplitude and angular frequency | $`x`$/s², rad/s | through $`r`$ and the strength below |

**Reference units.** Every characteristic below is stated at
$`\omega_n = 1`$ with $`b`$ of order 2. Section 6 moves it anywhere else.

**Dimensionless drive.** A drive enters through exactly two numbers:

```math
r = \frac{\Omega}{\omega_{lc}}, \qquad \text{strength} = \frac{A}{\omega_n^2\, b}
```

with $`\omega_{lc} = 2\pi/T`$ the *free cycle* frequency, not
$`\omega_n`$. Drive ratios in this sheet are quoted against Van der Pol's
$`\omega_{lc}`$ at the same $`\mu`$, as the regime maps are.

**Physical form.** For $`m\ddot{x} + c(x)\dot{x} + kx = F\cos\Omega t`$:
$`\omega_n = \sqrt{k/m}`$, $`c_{k} = 2\zeta_{k}\,m\,\omega_n`$,
$`A = F/m`$.

---

## 2. Applications

What the model is established to do (THREELEVEL §The proof at μ = 5,
§Van der Pol at μ = 1, §Beyond the fitted range, §The formulas alone):

- **Lock structure under drive.** For a Van der Pol class oscillator it
  places the 1:1, 3:1, 5:1 and 7:1 lock plateaus within one or two cells
  of 0.1 in drive ratio at every drive strength on the tested grid, from
  a fit at a single strength.
- **Chaos prediction.** It has chaos at every transition Van der Pol has
  it at and at no transition Van der Pol does not, from $`\mu = 1.5`$ at
  $`A = 5`$ and from $`\mu = 2.25`$ at $`A = 10`$, the bands within a
  twentieth to a tenth of a unit of drive ratio.
- **Parameters without fitting.** For Van der Pol at any $`\mu`$ from 0.1
  to 8 the three ratios come from closed form laws (§4.3); models built
  from the laws alone at $`\mu`$ never fitted put plateau edges within
  0.02 of Van der Pol's and chaos at the same transitions.
- **Hard excitation.** Quiet origin, threshold, large sustained
  oscillation: the sign pattern $`\zeta_{0} \gt 0`$, $`\zeta_{1} \lt 0`$,
  $`\zeta_{2} \gt 0`$ gives two nested cycles with the inner one the basin
  boundary. Free behaviour only; see warning W6.
- **Exact pieces.** Every arc is a closed form kernel, every crossing a
  scalar root find, and a locked cycle's Floquet multipliers come from a
  product of matrices without differencing (MAPS.md). Useful as a fast,
  exact integrator between events; not as a predictive difference
  equation (W9).

What it is *not* for: linear ringdown (use the linear prototype), a
nearly sinusoidal self-excited cycle with no interest in chaos (the two
level model does the same with fewer numbers), any stiffness
nonlinearity (W1).

---

## 3. Validity envelope and warnings

The model is only ever as good as the grid it was proved on. Inside the
envelope its behaviour is established; outside, it is untested, which is
not the same as wrong — but nothing here vouches for it.

### 3.1 Tested envelope

| quantity | tested range | where |
| --- | --- | --- |
| relaxation parameter $`\mu`$ | 0.1 to 5 fitted; 8 checked from the laws | THREELEVEL §The fits, §Beyond the fitted range |
| drive strength $`A/(\omega_n^2 b)`$ | 0 to about 5 ($`A`$ to 10 with $`b \approx 2`$) at $`\mu = 1`$, 2, 5; 2.5 ($`A = 5`$) elsewhere | §The proof at μ = 5, §Van der Pol at μ = 1, §μ = 2 across drive strength |
| drive ratio $`r`$ | 0.5 to 8 at $`\mu = 1`$, 2, 5; 0.5 to 6 elsewhere | same |
| $`\omega_n`$, amplitude scale | any; exact similarity, checked at 50 Hz/mm and 1 kHz/µm | §Moving it to another frequency range |
| chaos boundary in $`\mu`$ at $`A = 10`$ | located between 2 and 2.125 on both systems | §The chaos boundary in μ at A = 10 |
| data | exact Van der Pol integrations only | §Gaps |

### 3.2 Warnings

**W1 — Stiffness nonlinearity: the model does not apply at any
parameters.** Drive the system at two amplitudes. If the response
frequency shifts with amplitude, the nonlinearity is in the stiffness,
not the damping, and no prototype in this repository fits it
(VANDERPOL §Setting the parameters).

**W2 — Outside the envelope of 3.1 nothing is known.** Above $`\mu = 8`$,
above strength 5, above ratio 8: untested. The laws hold at 8 to 5%
on the ratios; at 10 they are a guess.

**W3 — Outer saturation.** Beyond $`b`$ the model is linear with damping
$`\zeta_{2}`$ for ever, and a smooth law keeps steepening. On the tested
grid the driven orbit never leaves $`\lvert x \rvert \approx 2.1`$ with
$`b \approx 2.0`$ to 2.1 because a stronger drive at these ratios grows
velocity, not displacement; a drive that does push the orbit well past
$`b`$ — low ratio, high strength — takes the model onto a plateau the real
system does not have (THREELEVEL §Why the saturation worry did not bite).

**W4 — Frequency offset.** The $`\mu = 5`$ fit's free period is 7% long
and every lock of it sits 2% low in ratio. The campaign fits are within
4.4% in period. Read plateau edges as accurate to a few per cent,
not better.

**W5 — Band widths are not reproduced.** The model says *whether* chaos
comes and *at which transition*, to a tenth of a unit of ratio. It does
not give how wide the chaotic window is: at $`\mu = 4`$, $`A = 10`$ its
first band is eight cells against Van der Pol's four, at $`\mu = 3`$ seven
against ten. Cell by cell agreement is 0.2 to 0.5.

**W6 — Hard excitation under drive is unmapped.** The bistable sign
pattern has its free cycles, multipliers and basin boundary and nothing
else. No regime map, no chaos verdicts.

**W7 — Extra narrow locks.** The model locks slightly more readily than a
smooth system: narrow locks of order 4 to 11 at integer drive ratios where
Van der Pol has tori, at every $`\mu`$ tested. Treat an isolated narrow
lock at an integer ratio as an artefact until confirmed.

**W8 — Multistability.** Where chaos coexists with a lock, which one the
classifier reports depends on the initial state. No basins have been
computed for the driven model.

**W9 — The exact map is not predictive.** The chain of zones a
trajectory visits can only be discovered by stepping, and at a grazing
(tangency with a threshold) the map's Jacobian is undefined. Use the
pieces as an event-driven integrator; do not expect a closed
recurrence (MAPS §The zone sequence cannot be known in advance).

**W10 — Below $`\mu \approx 0.3`$ the third level does nothing.** There is
no 3:1 plateau at $`A = 5`$, the driven targets barely depend on shape,
and the two level prototype is already behaviourally Van der Pol. The
laws still give parameters there, but the data cannot check them.

**W11 — Measured data untested.** All fits used exact integrations. The
classifier's thresholds are set from integrator noise; on a noisy sweep
they have not been exercised. The fit's own leeway is 20% on the free
cycle.

**W12 — The parameters are a shape, not a sampling.** The fitted levels
are not Van der Pol's damping law averaged over the zones: the core is
1.3 times the law's, the outer 1.2 to 1.5 times. Do not read a fitted
$`\zeta_{k}`$ as the physical damping at that displacement.

**W13 — Analyse in reference units.** The classifier in `section.py` has
absolute floors (integration tolerance $`10^{-11}`$ in $`x`$, a Lyapunov
threshold per unit time) set for $`\omega_n = 1`$, $`b \approx 2`$. Fit,
sweep and classify at the reference scale, then scale the results with
§6; do not run the classifier on a model in millimetres and kilohertz.

---

## 4. Characteristics

### 4.1 Operating modes by sign pattern

$`\langle\zeta\rangle(R)`$ is the cycle averaged damping at amplitude
$`R`$; it runs from $`\zeta_{0}`$ through $`\zeta_{1}`$ to $`\zeta_{2}`$
and a limit cycle sits at each zero crossing (THREELEVEL §The reduction
extends rather than restarts).

| mode | sign pattern | free response | under drive |
| --- | --- | --- | --- |
| damped | $`\zeta_{0} \gt 0`$, $`\langle\zeta\rangle`$ never negative | decays to rest | linear response, no locking |
| soft excitation, nearly harmonic | $`\zeta_{0} \lt 0 \lt \zeta_{1} \le \zeta_{2}`$, all ratios below about 1 | one cycle, near sinusoidal | two tongues (1:1, 3:1), tori elsewhere, one narrow period-doubled chaotic band near $`r = 0.5`$ at $`A = 1`$ |
| soft excitation, relaxation | $`\zeta_{0} \lt 0`$, $`\zeta_{2} \gg 1`$ | one cycle, crawl and jump | period adding 1:1, 3:1, 5:1, 7:1 with chaotic bands at the transitions once the outer levels are heavy |
| hard excitation, bistable | $`\zeta_{0} \gt 0`$, $`\zeta_{1} \lt 0`$, $`\zeta_{2} \gt 0`$ | origin stable, unstable inner cycle, stable outer cycle | unmapped (W6) |

### 4.2 Reference parameter sets

All at $`\omega_n = 1`$. Van der Pol's own cycle for comparison.

| set | $`\zeta_{0}, \zeta_{1}, \zeta_{2}`$ | $`a, b`$ | free amplitude | free period $`T\omega_n`$ | source |
| --- | --- | --- | --- | --- | --- |
| $`\mu = 5`$ proof fit (`staircase.THREE_FITTED`) | $`-1.735,\ 3.836,\ 15.05`$ | 1.075, 1.981 | 1.996 | 12.44 | fitted to plateau edges 2.4275, 2.4975 at $`A = 5`$ |
| Van der Pol $`\mu = 5`$ | — | — | 2.0215 | 11.612 | |
| $`\mu = 1`$ proof fit (`staircase.THREE_FITTED_MU1`) | $`-0.358,\ 0.865,\ 3.573`$ | 1.160, 1.984 | 1.981 | 6.670 | fitted to 2.195, 2.555 at $`A = 5`$ |
| Van der Pol $`\mu = 1`$ | — | — | 2.0086 | 6.663 | |
| bistable worked case (`staircase.BISTABLE_*`) | $`0.15,\ -0.25,\ 0.40`$ | 0.6, 1.6 | inner 1.16784 (unstable, multiplier 4.188), outer 2.25340 (stable, multiplier 0.1632) | | THREELEVEL §What the second threshold buys |

Campaign fits, one objective throughout (1:1 plateau end and 3:1 plateau
at $`A = 5`$, free cycle within 20%), from THREELEVEL §The fits:

| $`\mu`$ | $`\zeta_{0}`$ | $`\zeta_{1}`$ | $`\zeta_{2}`$ | $`a`$ | $`b`$ | free amplitude | free $`T\omega_n`$ | Van der Pol $`T`$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.1 | $`-0.043`$ | 0.090 | 0.22 | 1.203 | 2.124 | 2.132 | 6.288 | 6.287 |
| 0.2 | $`-0.093`$ | 0.178 | 0.52 | 1.226 | 2.177 | 2.233 | 6.305 | 6.299 |
| 0.3 | $`-0.108`$ | 0.246 | 0.69 | 1.278 | 2.238 | 2.203 | 6.316 | 6.318 |
| 0.5 | $`-0.190`$ | 0.444 | 1.47 | 1.255 | 2.081 | 2.127 | 6.393 | 6.381 |
| 0.7 | $`-0.237`$ | 0.537 | 2.48 | 1.246 | 2.286 | 2.160 | 6.440 | 6.473 |
| 1 | $`-0.362`$ | 0.924 | 3.37 | 1.209 | 2.159 | 2.030 | 6.701 | 6.663 |
| 1.5 | $`-0.543`$ | 1.214 | 5.43 | 1.163 | 2.025 | 2.048 | 7.138 | 7.096 |
| 2 | $`-0.730`$ | 1.544 | 7.30 | 1.165 | 2.054 | 2.091 | 7.756 | 7.630 |
| 3 | $`-1.074`$ | 2.353 | 10.72 | 1.131 | 2.018 | 2.047 | 9.252 | 8.859 |
| 4 | $`-1.319`$ | 3.020 | 16.59 | 1.161 | 2.109 | 2.118 | 10.469 | 10.204 |
| 5 | $`-1.563`$ | 4.025 | 22.42 | 1.189 | 2.179 | 2.094 | 11.961 | 11.612 |
| 8 (outside the law fit) | $`-2.475`$ | 6.136 | 39.13 | 1.178 | 2.303 | 2.119 | 16.686 | 16.0 |

### 4.3 Parameter laws

Power laws through the eleven campaign fits, $`\mu`$ from 0.1 to 5, edges
constant (THREELEVEL §Formulas):

```math
\zeta_{0} = -0.365\,\mu^{0.93}, \qquad
\zeta_{1} = 0.824\,\mu^{0.96}, \qquad
\zeta_{2} = 3.29\,\mu^{1.17}, \qquad
a = 1.20 \pm 0.05, \quad b = 2.13 \pm 0.09
```

in units where Van der Pol's free amplitude is 2 and $`\omega_n = 1`$.
Evaluated:

| $`\mu`$ | 0.3 | 0.5 | 1 | 2 | 3 | 5 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $`\zeta_{0}`$ | $`-0.12`$ | $`-0.19`$ | $`-0.37`$ | $`-0.70`$ | $`-1.01`$ | $`-1.63`$ | $`-2.52`$ |
| $`\zeta_{1}`$ | 0.26 | 0.42 | 0.82 | 1.61 | 2.37 | 3.88 | 6.10 |
| $`\zeta_{2}`$ | 0.80 | 1.46 | 3.29 | 7.40 | 11.9 | 21.6 | 37.5 |

Notes. The outer level is weakly determined and does not matter: fits at
$`\mu = 5`$ put it at 15 and at 22 with the plateau edges matched equally
well, because the driven orbit barely enters the outer zone. Set it from
the law and leave it. The core and band are proportional to $`\mu`$ to
within their exponents' distance from one, which is how Van der Pol's own
law scales.

### 4.4 Driven response — lock plateaus

The $`\mu = 5`$ proof fit against Van der Pol, drive ratio units, Van der
Pol's $`\omega_{lc}`$ (THREELEVEL §The proof at μ = 5):

| $`A`$ | lock 1 | lock 3 | lock 5 | lock 7 |
| --- | --- | --- | --- | --- |
| 0.5 | 0.8–1.1 vs 0.9–1.1 | 2.7–2.9 vs 2.9–3.1 | 4.5–4.8 vs 4.9–5.1 | 6.4–6.6 vs 7.0 |
| 1 | 0.6–1.3 vs 0.7–1.3 | 2.4–3.1 vs 2.7–3.3 | 4.3–5.0 vs 4.8–5.2 | 6.3–6.8 vs 6.9–7.1 |
| 2 | 0.5–1.6 vs 0.5–1.7 | 2.2–3.5 vs 2.2–3.6 | 4.1–5.3 vs 4.5–5.5 | 6.1–7.0 vs 6.8–7.3 |
| 5 | to 2.7 vs to 2.6 | 2.8–4.4 vs 2.7–4.4 | 4.7–6.1 vs 4.7–6.1 | 6.6–7.8 vs 6.7–7.8 |
| 10 | to 4.1 vs to 4.0 | 4.2–5.8 vs 4.2–5.6 | 5.9–7.3 vs 5.7–7.1 | 7.6–8 vs 7.4–8 |

Van der Pol's plateau edges at $`A = 5`$ against $`\mu`$, the campaign's
targets, all matched by the fits to within 0.04 (§The fits):

| $`\mu`$ | 0.3 | 0.5 | 1 | 1.5 | 2 | 3 | 4 | 5 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:1 ends | 2.11 | 2.13 | 2.21 | 2.29 | 2.38 | 2.49 | 2.59 | 2.68 | 2.81 |
| 3:1 plateau | 2.87–3.12 | 2.77–3.21 | 2.56–3.41 | 2.47–3.61 | 2.47–3.79 | 2.53–4.11 | 2.59–4.33 | 2.68–4.48 | 2.81–4.74 |

Below $`\mu = 0.3`$ there is no 3:1 plateau at this drive (W10). From
$`\mu = 4`$ the 1:1 and 3:1 plateaus meet directly.

### 4.5 Driven response — where chaos is

Rule (THREELEVEL §Predicting chaos): **where the model's lock plateaus
meet, chaos is coming; where they meet cleanly, it is not.** At $`A = 5`$:

| $`\mu`$ | Van der Pol | model |
| --- | --- | --- |
| ≤ 1 | none | none |
| 1.5, 2 | one or two cells at the 1:1 → 3:1 transition | same transition, within 0.05 |
| 3 | 1:1 → 3:1; 3:1 → 4:1; 4:1 → 5:1 | 1:1 → 3:1; 4:1 → 5:1 (misses 3:1 → 4:1, switches through a torus) |
| 4, 5 | transitions above the 3:1 plateau | same region, cells within 0.05 to 0.1 |

At $`A = 10`$ (§The chaos boundary in μ at A = 10): first chaos between
$`\mu = 2`$ and 2.125 in both systems, at the end of the 3:1 plateau;
above the 5:1 plateau from between 3 and 3.5 in both; below the 7:1
plateau by 4 in Van der Pol and 3.5 in the model. At $`\mu = 5`$,
$`A = 10`$ the 5:1 → 7:1 band has fifteen confirmed cells on both, the
model's 2% high in ratio.

Nearly harmonic mode, $`\mu = 1`$: the only chaos is a band reached by
period doubling of the 1:1 lock at $`A = 1`$, Van der Pol at
$`r = 0.48`$–0.50, the model at 0.56–0.57.

---

## 5. Tuning procedures by scenario

Read §3 first. Each procedure lists what to measure, what it sets, and
the check that says the model holds.

### 5.1 Is this the right model? (gate, before anything else)

1. Drive at two amplitudes. Frequency shifts → stiffness nonlinearity,
   stop (W1).
2. Free waveform: nearly sinusoidal with no interest in chaos → two level
   model, `README.md` §Switching on displacement. Ringdown only → linear
   prototype.
3. Drive at two amplitudes and watch the *driven* response: if the lock
   plateaus and the transitions between them move with drive amplitude,
   the damping shape across the visited range matters and this model is
   needed, fitted to those transitions.

### 5.2 Van der Pol class oscillator, no fitting

For a system whose damping is known or believed to follow
$`-\varepsilon(1 - x^2/X^2)`$.

| step | measure | sets | how |
| --- | --- | --- | --- |
| 1 | free period $`T`$ and free amplitude $`R`$ | — | ringdown to the cycle, or steady state |
| 2 | relaxation parameter | $`\mu`$ | $`\mu = \varepsilon/\omega_n`$ if $`\varepsilon`$ is known; else from the waveform: Van der Pol's $`T\omega_n`$ runs 6.29 at $`\mu = 0.1`$, 6.66 at 1, 7.63 at 2, 11.6 at 5, 16.0 at 8 (THREELEVEL §Where the targets sit), so $`T\omega_n`$ reads $`\mu`$ once $`\omega_n`$ is known, and the third harmonic ratio is about 0.1 at $`\mu = 1`$ |
| 3 | — | $`\zeta_{0}, \zeta_{1}, \zeta_{2}`$ | the laws of §4.3 at that $`\mu`$ |
| 4 | — | $`a, b`$ | $`(1.20, 2.13)\,\lambda`$ with $`\lambda = R/2`$; equivalently $`\lambda = X`$, the displacement at which the damping changes sign |
| 5 | — | $`\omega_n`$ | $`\omega_n = (T\omega_n)_{\text{model}}/T`$, with $`(T\omega_n)_{\text{model}}`$ from the campaign table at that $`\mu`$ (not $`2\pi/T`$: W4, §6) |
| check | one lock plateau at one drive strength | — | if its edges are where the model puts them to a few per cent, done; if not, go to 5.3 |

Cost: none beyond one sweep for the check. Accuracy: plateau edges within
0.02 in ratio and chaos at the same transitions, tested at $`\mu = 2.5`$
and 3.5 without fitting (THREELEVEL §The formulas alone).

### 5.3 Self-excited oscillator of unknown damping law — fit to the driven response

The free cycle cannot set the three ratios: nothing measured without a
drive depends on them beyond a few per cent. A ringdown is not enough;
two drive strengths are (THREELEVEL §What to measure).

| step | measure | sets | how |
| --- | --- | --- | --- |
| 1 | free $`T`$, free $`R`$ | leeway targets | held within 20% during the fit |
| 2 | at one drive strength, the last frequency locked $`p{:}1`$ and the first locked $`(p+2){:}1`$ — two plateau edges of adjacent locks | $`\zeta_{0}, \zeta_{1}, \zeta_{2}, a, b`$ | sweep drive frequency at fixed strength, spacing finer than a tongue edge; `staircase.fit_bands` moves all five parameters by Nelder–Mead until the model's edges sit on them |
| 3 | a lock plateau at a second drive strength | check, not a parameter | edges where the model puts them → the fit holds across drive; if not, the model is a fit at one strength only |
| 4 | free $`T`$ | $`\omega_n`$ | as 5.2 step 5 |
| 5 | free $`R`$ | $`\lambda`$ | scale the fitted edges by the measured amplitude over the model's |

Cost: each evaluation is a coarse sweep with bisection on the edges,
about a minute on four cores; a fit takes forty to sixty. Start from the
laws of §4.3 at the nearest $`\mu`$; from $`\mu = 1.5`$ up every campaign
fit started that way landed within a coarse step before the first
evaluation. Do not expect the fitted levels to be physical damping (W12).

### 5.4 Predicting chaos at a working drive

With a model from 5.2 or 5.3:

1. Sweep its drive frequency at the working drive strength, in reference
   units (W13), at 0.01 to 0.02 in ratio across the range of interest.
2. Read where the lock plateaus meet. Chaos comes at the first transition
   above the 1:1 lock once the ratios are of order one, and at the
   transitions above the 3:1 lock once they are of order three; a clean
   meeting means none.
3. Report the transition frequencies as the model's to within a twentieth
   of a unit of drive ratio (a tenth at $`A = 10`$), and the *width* of any
   chaotic window as unknown (W5).
4. Confirm any chaotic verdict with `section.confirm_chaos`; isolated
   narrow locks at integer ratios are suspect (W7).

### 5.5 Hard excitation

For a system that sits quietly until knocked past a threshold, then
sustains a large oscillation.

| measure | sets | how |
| --- | --- | --- |
| small-signal decrement at the origin | $`\zeta_{0} \gt 0`$ | logarithmic decrement, as for the linear prototype |
| threshold amplitude | inner unstable cycle radius | the smallest kick that leads to the large oscillation |
| sustained amplitude | outer stable cycle radius | steady state |

The two radii are the two zero crossings of $`\langle\zeta\rangle(R)`$
(`staircase.cycles_predicted`, exact by `staircase.cycles_exact`). With
$`\zeta_{0}`$ measured they constrain the remaining four numbers
$`\zeta_{1}, \zeta_{2}, a, b`$ by two equations: choose the outer level
and the edge ratio, solve for the band level and the edge scale; from
the worked case's averaged radii that solve recovers its $`-0.25`$ and
1.6 from any reasonable start. The worked case $`(0.15, -0.25, 0.40)`$,
edges $`(0.6, 1.6)`$, has exact radii 1.168 and 2.253. **Free behaviour only** — the driven response of this
pattern is unmapped (W6), and there is no recipe for the two free
choices beyond the worked case.

### 5.6 Moving an existing model to another frequency or amplitude

See §6. No re-fitting; the behaviour is exactly preserved.

### 5.7 Choosing a simpler model

| situation | use instead |
| --- | --- |
| $`\mu \lesssim 0.3`$, or nearly sinusoidal cycle | two level model, displacement switched; the third level does nothing (W10) |
| chaos under drive wanted with the fewest parameters, inside of the chaotic region not needed | two level model, $`\zeta \approx (-1.24, 8.33)`$, $`x_0 = 1.44`$ at $`\mu = 5`$ (VANDERPOL §Setting the parameters) |
| free response in closed form, no drive | Hopf normal form |
| ringdown | linear prototype |

---

## 6. Frequency and amplitude scaling

To place the reference model at natural frequency $`\omega_n`$ with
$`\lambda`$ units of displacement per model unit (THREELEVEL §Moving it
to another frequency range, checked by `scaling.py`):

```math
\zeta_{k} \to \zeta_{k}, \qquad
(a, b) \to (\lambda a,\ \lambda b), \qquad
\Omega \to \omega_n\Omega, \qquad
A \to \lambda\,\omega_n^2 A
```

| quantity | multiply by |
| --- | --- |
| times: period, settling, transient | $`1/\omega_n`$ |
| frequencies: free cycle, drive, plateau edges in rad/s | $`\omega_n`$ |
| displacements: edges, amplitudes, basin boundary | $`\lambda`$ |
| velocities | $`\lambda\omega_n`$ |
| accelerations, the drive amplitude $`A`$ included | $`\lambda\omega_n^2`$ |
| damping coefficients $`c_{k}`$ per unit mass | $`\omega_n`$ |
| stiffness per unit mass | $`\omega_n^2`$ |
| ratios, $`a/b`$, $`\mu`$, $`r`$, strength, lock orders, multipliers, verdicts | 1 |
| Lyapunov exponents and any rate | $`\omega_n`$ |

Worked, the $`\mu = 5`$ fit driven in its chaotic band at $`\Omega = 2.47`$,
$`A = 5`$:

| | reference | 50 Hz, 1 mm per unit | 1 kHz, 1 µm per unit |
| --- | --- | --- | --- |
| $`\omega_n`$ | 1 rad/s | 314.2 rad/s | 6283 rad/s |
| edges | 1.075, 1.981 | 1.075, 1.981 mm | 1.075, 1.981 µm |
| free amplitude, period | 1.996, 12.44 s | 1.996 mm, 39.6 ms (25.3 Hz) | 1.996 µm, 1.98 ms (505 Hz) |
| drive $`\Omega`$, $`A`$ | 2.47 rad/s, 5 | 776 rad/s (123.5 Hz), 493 m/s² | 15 519 rad/s (2.47 kHz), 197 m/s² |
| strength $`A/(\omega_n^2 b)`$ | 2.524 | 2.524 | 2.524 |

Two traps. The drive scales with $`\omega_n^2`$: doubling the frequency at
the same amplitude needs four times the drive acceleration to stay at the
same point of the regime map. And $`\omega_n`$ is not the free cycle
frequency: a relaxation oscillator measured at $`f`$ hertz needs
$`\omega_n = (T\omega_n)_{\text{model}}\,f`$, well above $`2\pi f`$.

---

## 7. Reproduction and where the evidence is

| claim | document section | script |
| --- | --- | --- |
| definition, averaged damping, two cycles, bistable case | THREELEVEL §Definition | `python3 staircase.py` |
| behaviour map by sign pattern | THREELEVEL §Behaviours | — |
| $`\mu = 5`$ fit and driven proof across the grid | THREELEVEL §Fitting it, §The proof at μ = 5 | `staircase.py fit`, `staircase.py regime` |
| $`\mu = 1`$ fit and proof | THREELEVEL §Van der Pol at μ = 1 | `staircase.py fit1`, `regime1` |
| campaign fits, laws, verification, $`\mu = 8`$, $`\mu = 2`$ grid | THREELEVEL §Parameters against the relaxation parameter | `campaign.py survey`, `fit`, `verify`, `check`, `formula` |
| chaos boundary in $`\mu`$ at $`A = 10`$ | THREELEVEL §The chaos boundary | `campaign.py boundary` |
| scaling rule and its check | THREELEVEL §Moving it to another frequency range | `scaling.py`, `scaling.py table` |
| exact arcs, crossings, Floquet multipliers | MAPS.md | `maps.py` |
| classifier and its thresholds | `section.py` docstring; MAPS §The exponent, measured | `section.py` |
| why drive, not free cycle, sets the ratios | VANDERPOL §normalisation chapter | `staircase.py normalise` |
| model selection across the repository | VANDERPOL §Setting the parameters | — |

Figures: `python3 figures.py` regenerates all of them. Fitted constants
are stored in `staircase.py`, so the tables here do not depend on
re-running any optimiser.
