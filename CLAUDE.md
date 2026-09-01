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

1. branch off the current `main`, onto a **new** branch named for that
   change — never reuse a branch across pull requests, not even one whose
   pull request has already merged
2. commit and push the branch
3. open a pull request against `main` and hand over the link
4. leave it for the repository owner to review and merge

One branch per pull request keeps each change reviewable on its own and
allows more than one to be open at a time.

### Branches cannot be deleted from a session

Deleting a remote branch is blocked by the agent proxy, by both routes:

- `git push origin --delete <branch>` fails with HTTP 403
- `DELETE /repos/.../git/refs/heads/<branch>` returns *"Write access to
  this GitHub API path is not permitted through this proxy"*

So never promise to delete a head branch after closing a pull request —
it cannot be done from here. Say plainly that it needs doing outside the
session.

That does not leave branches lying around, because the repository has
**Settings → General → Automatically delete head branches** switched on.
GitHub removes each head branch when its pull request merges, with no
session involvement at all. Nothing needs doing about branch cleanup.

If merged branches ever start accumulating again, that setting is what to
check. It is not writable from a session (*"Repository settings writes are
not permitted through this proxy"*), so it has to be set by hand. The
reliable way to confirm it is on is to merge a pull request and see whether
its branch disappears — the setting once read as off while a person
believed they had enabled it.

## Maths in markdown

Three rules, each of which was a real rendering bug before it was one:

1. **Display equations** go in fenced ` ```math ` blocks, never `$$ ... $$`.
   GitHub's markdown pass strips the `\\` line breaks inside `$$`, which
   silently breaks any `aligned` or `bmatrix` environment.
2. **Inline maths** uses the `` $`...`$ `` form, never bare `$ ... $`. The
   backticks shield the content from the markdown pass. With bare `$`, a
   mis-paired delimiter desynchronises matching for the rest of the
   paragraph and takes neighbouring expressions down with it — 17 of 131
   inline expressions were silently not rendering.
3. **Never put a literal `<` or `>` inside inline maths.** Use `\lt` and
   `\gt`. GitHub HTML-escapes the characters to `&lt;` and `&gt;` inside
   the maths, and MathJax then draws the entity text literally. Display
   blocks escape correctly and do not need this.

Fenced math blocks and `` $`...`$ `` are GitHub extensions: they render on
github.com but show as plain code in editor previews without a math
extension.

### Checking that the maths actually renders

Do not trust a local preview — GitHub's pipeline is what matters, and it
fails silently. Push the branch and read back GitHub's own server-side
render:

```
curl -sS -H "User-Agent: Mozilla/5.0" \
  "https://github.com/aburston/second-order-prototype/tree/<branch>" > page.html
```

Count `js-inline-math` and `js-display-math` elements in the result and
compare against the expressions in the file. Any expression missing an
element is not rendering; any element containing `&lt;`, `&gt;` or `&amp;`
will show entity text to the reader. Both should be zero.

## Figures

`figures.py` generates every image in the README into `figures/`, each one
rendered twice — `-light.png` and `-dark.png` — and embedded through a
`<picture>` element so GitHub serves the right one for the reader's theme.
Regenerate with `python3 figures.py` after changing any parameter the README
quotes, and commit the PNGs alongside the code.

`examples.py` does the same for `EXAMPLES.md`, writing `figures/example-*.png`
and printing every number that file quotes. It imports the plotting chrome
from `figures.py`, so the two stay in one style. The prototypes are never
altered in `EXAMPLES.md`: a physical system is either one of them exactly or
is integrated separately as the reference a tuned prototype is compared with.

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
