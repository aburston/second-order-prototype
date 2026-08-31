# Working notes for Claude

## The project

A second order nonlinear prototype. `README.md` develops it in stages:
the linear prototype, then a piecewise linear system with a switching
boundary in the phase plane where the **damping coefficient** is the
switched quantity, then the offset boundary case that produces a
hyperbolic limit cycle. `limit_cycle.py` is the Poincare map analysis
behind the numbers quoted in the last section.

Read the README before extending it — each section builds on the previous
one and reuses its notation ($x_1 = x$, $x_2 = \dot{x}$, $\zeta_{\pm}$,
$\bar{\zeta}$, $\Sigma$).

## Workflow: use pull requests

Do **not** merge into `main` directly. For every change:

1. branch off the current `main`
2. commit and push the branch
3. open a pull request against `main` and hand over the link
4. leave it for the repository owner to review and merge

Ref deletion is blocked for remote sessions, so merged branches cannot be
cleaned up from a session — leave them and let GitHub's "automatically
delete head branches" setting handle it.

## Maths in markdown

Write display equations as fenced ` ```math ` blocks, not `$$ ... $$`.
GitHub's markdown pass strips the `\\` line breaks inside `$$` blocks,
which silently breaks any `aligned` or `bmatrix` environment. Inline
`$ ... $` on a single line is fine.

Fenced math blocks are a GitHub extension: they render on github.com but
show as plain code in editor previews without a math extension.

## Figures

`figures.py` generates every image in the README into `figures/`, each one
rendered twice — `-light.png` and `-dark.png` — and embedded through a
`<picture>` element so GitHub serves the right one for the reader's theme.
Regenerate with `python3 figures.py` after changing any parameter the README
quotes, and commit the PNGs alongside the code.

Colours come from a validated categorical palette (slots 1-3: blue, orange,
aqua) with the switching boundary, equilibrium markers and callouts drawn in
chrome ink rather than a series colour, so colour never carries identity on
its own. Every series is also directly labelled.

## Verify before committing

The README states quantitative results — decay factors, stability
conditions, multipliers, scaling laws. Check them numerically before
committing rather than asserting them from the algebra, and say plainly in
the text which claims are verified and which are only observed over a
tested range. A sign error in the equilibrium formula survived a first
reading and was caught only by evaluating the vector field at the claimed
point; do that kind of check routinely.
