# Learn You a Tropical Geometry for Great Good!

## Chapter 13 — Trees Are Invisible

> **DRAFT — Phase 1, 2026-11-17.** Converted from Day 9 of
> the study: prose converted, program assembled
> (`hidden_files/code/ch13_chapter_program_draft.py` —
> the parameterized Baker–Norine machinery of
> `chip_firing.py` inlined into Day 9's
> `day09_full_cubic.py`), figures staged from the study.
> **NOT yet re-verified end-to-end from this draft text** —
> that is Phase 2, together with ClaudeDuet's referee pass
> in #minplus-book. Phase 3 compresses and ships. This is
> the last chapter of the main trail: Parts I–VI complete
> here; the Deep Dives follow.

Chapter 12 cheated, honestly and on purpose. It took the
smooth tropical cubic of Chapter 3, threw away everything
except the hexagon at its centre, and computed the group
law on that cycle alone: `Jac(C₆) ≅ ℤ/6`, Abel–Jacobi a
bijection, addition = rotation. But the actual cubic is
bigger than its cycle. It has nine vertices, not six: the
hexagon, plus three extra vertices hanging off it. This
chapter puts them back — and the headline is that the
Jacobian refuses to notice:

> **Graft any trees you like onto the cycle. Every divisor
> pushes onto the cycle, Abel–Jacobi collapses each tree
> to its attachment point, the canonical divisor stays
> principal, Riemann–Roch still holds — and
> `Jac(G) ≅ ℤ/6`, exactly as before. Trees are invisible.**

The one thing that *does* break is injectivity of
Abel–Jacobi — and it breaks in exactly the place
Chapter 12 flagged: §12.5's bijection rested on
Baker–Norine Theorem 1.8, whose hypothesis is
2-edge-connectedness. A tree hangs by a single bridge, so
the hypothesis fails, and the computation shows the
failure is *precisely* the trees: nothing else collapses.

Sources are Chapter 12's: Baker–Norine chip-firing theory
(M. Baker & S. Norine, arXiv:math/0608360), with the
metric coda from Mikhalkin–Zharkov (arXiv:math/0612267).

### 13.1 The full graph, rebuilt — not remembered

The program does not hardcode the graph. It re-runs
Chapter 3's pipeline — heights `c(i,j) = i² + j² + ij` on
the degree-3 lattice triangle, lower subdivision, dual
curve — and reads the bounded graph off the result:
9 vertices (numbered by their dual faces), 9 bounded
edges, 9 rays, balancing asserted at every vertex. The
bounded edges are

```
cycle:   (3,4) (4,8) (8,7) (7,6) (6,5) (5,3)
bridges: (0,4) (1,5) (2,7)
```

so the genus is unchanged,

```
g = E − V + 1 = 9 − 9 + 1 = 1.
```

**The 2-core.** Strip leaves repeatedly — a vertex of
degree ≤ 1 in what remains is peeled — and what survives
is the **2-core**, the graph's minimal skeleton. Here
leaf-stripping removes exactly vertices `0, 1, 2` and
stops: the core is `{3, 4, 5, 6, 7, 8}`, every core vertex
has exactly 2 core neighbours, and walking it from
`q = 3` gives the cycle in the order

```
3 → 4 → 8 → 7 → 6 → 5 → (back to 3).
```

Each peeled vertex hangs off the core by one **bridge** —
an edge whose removal disconnects the graph — and the
program asserts each off-core component touches the core
at exactly one vertex (two attachment points would close
a second cycle, contradicting the core's maximality):

```
pendant 0 attaches at 4,   pendant 1 attaches at 5,
pendant 2 attaches at 7.
```

In this cubic each "tree" is a single pendant vertex, the
smallest tree there is. The theorems below are stated and
checked for these trees, but nothing in the arguments
uses their smallness — only that they meet the core in
one vertex and contain no cycle. Vertex degrees on `G`:
the attachments `4, 5, 7` have degree 3, the other core
vertices degree 2, the pendants degree 1. Keep the base
point on the cycle: `q = 3`, as in the metric data of
Day 9.

### 13.2 Predict first: the spanning-tree bet

Chapter 12's first route to `|Jac| = 6` was Kirchhoff:
the order of the Jacobian is the number of spanning
trees. Before computing anything on `G`, predict.

A spanning tree must connect every vertex, so it must
contain **every bridge** — delete a bridge and its tree
is cut off, with no other route back. The three bridges
are therefore in all spanning trees, forced. What remains
free is the cycle part, and a spanning tree of a 6-cycle
is the cycle minus one edge: 6 choices. Prediction:

```
# spanning trees of G = 6,   hence   |Jac(G)| = 6
```

— the trees contribute *choices of nothing*: every
spanning tree contains them entirely. If the group also
comes out cyclic, Chapter 12's group survives grafting
untouched.

### 13.3 The Jacobian is still `ℤ/6`

The program collects the bet two ways, as in §12.4:

1. **Matrix-Tree.** `det` of the reduced Laplacian of `G`
   is `6` — the prediction of §13.2, on the nose.
2. **Smith normal form.** The reduced Laplacian is now
   8×8, and its SNF is `diag(1, 1, 1, 1, 1, 1, 1, 6)`:
   seven trivial factors and the same single `6`. So

   ```
   Jac(G) ≅ ℤ/6,
   ```

   cyclic, not just of order 6 — and *identical* to
   Chapter 12's group, not merely the same size. The
   grafted trees add no torsion, no new classes, nothing.

A third count arrives in §13.5: the Abel–Jacobi image has
exactly 6 distinct classes, the same six as Day 8's.

### 13.4 The push-in lemma: every divisor lives on the cycle

Why is the Jacobian unchanged? Because, up to linear
equivalence, there was never anything on the trees to
begin with. Take a divisor `D` and a pendant leaf `p` with
neighbour `u`. Firing `p` exactly `D(p)` times (borrowing
if `D(p)` is negative) moves all of `p`'s chips to `u`:

```
D(u) += D(p),   D(p) = 0,
```

and firing is principal by definition, so the class
`[D]` does not change. Peel the trees leaf by leaf —
a leaf's neighbour is always still present, since it lies
on the path toward the core — and every chip ends up on
the 2-core. Hence:

> **Push-in lemma.** Every divisor on `G` is linearly
> equivalent to one supported on the cycle.

The program verifies this in bulk, not by example: for
**all `3⁹ = 19,683` divisors with entries in `{-1, 0, 1}`**
it pushes the divisor in and asserts the q-reduced
representative — the class fingerprint of Chapter 12's
stabiliser — is unchanged. A worked instance, drawn in
Figure 13.1B: the degree-0 divisor with 2 chips on
pendant 0, −1 on each of pendants 1 and 2, +1 on core
vertex 4 and −1 at `q` pushes to a divisor supported
entirely on `{3,…,8}`, same class, trees empty.

This is the computational shadow of a classical move:
contracting a bridge does not change the Jacobian
(Baker–Norine, Cor. 4.7, in the study's notes from
Day 8). Retract the whole graph onto its minimal skeleton
and divisor theory retracts with it.

### 13.5 Abel–Jacobi collapses the trees

Compute `S_q(v) = [(v) − (q)]` for all nine vertices, as
q-reduced representatives in vertex order `0,…,8`. Every
one is already reduced as written:

```
v₀ ↦ [0, 0, 0, -1, 1, 0, 0, 0, 0]
v₁ ↦ [0, 0, 0, -1, 0, 1, 0, 0, 0]
v₂ ↦ [0, 0, 0, -1, 0, 0, 0, 1, 0]
v₃ ↦ [0, 0, 0,  0, 0, 0, 0, 0, 0]      (q itself)
v₄ ↦ [0, 0, 0, -1, 1, 0, 0, 0, 0]
v₅ ↦ [0, 0, 0, -1, 0, 1, 0, 0, 0]
v₆ ↦ [0, 0, 0, -1, 0, 0, 1, 0, 0]
v₇ ↦ [0, 0, 0, -1, 0, 0, 0, 1, 0]
v₈ ↦ [0, 0, 0, -1, 0, 0, 0, 0, 1]
```

Read the collisions off the table:

```
S_q(0) = S_q(4),   S_q(1) = S_q(5),   S_q(2) = S_q(7).
```

**Each pendant has the Abel–Jacobi class of its
attachment point.** It could not be otherwise — `(0) −
(4)` is the divisor of firing the leaf once, hence
principal — but the table makes the factorisation
visible: the Abel–Jacobi map on `G` is the Abel–Jacobi
map on the cycle, precomposed with "walk to your
attachment point". Its image is exactly Chapter 12's six
classes.

And so injectivity fails, by exactly three collisions.
Chapter 12's bijection used Baker–Norine Theorem 1.8:
`S^(1)` is injective iff the graph is 2-edge-connected.
`G` is not — each bridge is a 1-edge cut — and the
computation shows the theorem's failure mode is not
diffuse: the *only* identifications are pendant =
attachment. On the cycle itself the map is still
injective; the trees are the whole story.

The group law agrees. Number the cycle positions from
`q`: position of `3, 4, 8, 7, 6, 5` is `0, 1, 2, 3, 4, 5`,
and Chapter 12's rotation law says a vertex's class is
its position in `ℤ/6`. The pendants inherit their
attachments' positions — `0 ↦ 1`, `1 ↦ 5`, `2 ↦ 3` — the
program checks the generator `S_q(4)` (position 1) still
has order exactly 6, and the full consistency law

```
[pendant] + [x] = [attachment] + [x]   for every vertex x,
```

so a pendant does not merely *sit* in the same class as
its attachment; it *acts* identically in the group.

### 13.6 The canonical divisor is principal; Riemann–Roch survives

On the cycle, Chapter 12 had `K = 0` for free — every
vertex had degree 2. On `G` the canonical divisor
`K = Σ (deg(v) − 2)·(v)` is honestly nonzero:

```
K = [-1, -1, -1, 0, +1, +1, 0, +1, 0]
```

`−1` at each pendant (degree 1), `+1` at each attachment
(degree 3), `0` at the other core vertices — and
`deg K = 0 = 2g − 2`, as it must at genus 1. Is `K`
principal? Yes, and the witness is one line: `K` is
*minus* the sum of the three leaf columns of the
Laplacian,

```
K = −( Q[:,0] + Q[:,1] + Q[:,2] ),
```

so borrowing once at each leaf (the firing script
`f = (−1, −1, −1, 0, 0, 0, 0, 0, 0)`) sends `0` to `K`.
(Day 9's log said "sum of the Laplacian columns, one
firing per leaf" — the class statement was right, the
sign was loose: it is the negative sum, i.e. one
*borrowing* per leaf. The draft states it with the sign;
Phase 2 asserts the equality as written here.)

With `K` principal, `r(K − D) = r(−D)`, and tropical
Riemann–Roch on `G` again reads

```
r(D) − r(K − D) = deg(D) + 1 − g = deg(D).
```

The program verifies it on **8,953 divisors** — all
entries in `{-1, 0, 1}` with degree in `{-1, 0, 1}` —
rank computed from Baker–Norine's definition, exactly as
in Chapter 12. Grafting trees changed the canonical
divisor from `0` to a nonzero principal divisor, and the
theorem did not blink: another sense in which divisor
theory sees only the skeleton.

> **Honest Box — two small catches from Day 9, kept.**
> (1) sympy's `smith_normal_form` first returned
> `diag(1×7, −6)` for the reduced Laplacian — invariant
> factors are defined up to sign, so the code takes
> absolute values; a negative "group order" is a
> convention, not a discovery. (2) The first cycle-walk
> for the metric coda asserted every step had exactly
> one unvisited core neighbour, which is false at the
> *first* step, where the start vertex has two. Both
> were caught by the program's own assertions, and both
> fixes are in the assembled chapter program.

### 13.7 Metric coda: the retraction, drawn

Give the cycle its tropical metric — Chapter 12's
§12.8 correction applied from the start: every cycle edge
has **lattice length 1** (the two diagonal edges have
displacement `(1,1)`, exactly one primitive step), so the
total length is

```
L = 6,   Jac_metric(G) = ℝ / 6ℤ,
```

not the Euclidean `4 + 2√2` that Day 9's figure still
uses for its spacing (the Phase 3 regeneration redraws
it, as it did for Chapter 12). The metric Abel–Jacobi map
sends a point at lattice arclength `x` from `q` to `x`
around that circle; the cycle vertices land at the
integers, in cycle order:

```
3 ↦ 0,   4 ↦ 1,   8 ↦ 2,   7 ↦ 3,   6 ↦ 4,   5 ↦ 5.
```

What does it do to a point *inside* a tree — here, to a
pendant vertex? The tree is a bridge segment hanging off
the cycle; moving along it changes no divisor class
(§13.5), and metrically the map is constant on it: the
whole tree lands on its attachment point,

```
0 ↦ 1,   1 ↦ 5,   2 ↦ 3,
```

the same positions as §13.5's group computation, because
they are the same statement — combinatorially, the
retraction of `G` onto its 2-core induces an isomorphism
of Jacobians; metrically, Abel–Jacobi factors through
that retraction. This is the Baker–Payne–Rabinoff
picture in miniature, the one the study flagged on
Day 9: tropicalisation retracts a curve onto its
skeleton, and everything divisor-theoretic — group,
group law, Riemann–Roch — lives on the skeleton. The
trees were decoration. Invisible, twice over: once to
the group, once to the circle.

That closes Part VI, and with it the main trail of this
book. The cubic of Chapter 3 has now given up everything
it was carrying: a cycle (Ch. 3), a group law (Ch. 12),
and the proof that the group law never needed the rest
of the curve (this chapter). The Deep Dives take the
same machinery to bigger graphs — the Petersen graph of
Chapter 11 gets a Jacobian of its own there — where the
skeletons stop being single cycles and the groups stop
being single cyclic factors.

![Chip-firing on the full tropical cubic: the graph and its AJ classes, a push-in demonstration, the metric Jacobian, and K principal](figures/day09-trees-invisible.png)

*Figure 13.1 — A: the full cubic: red 6-cycle, grey
bridges, blue rays; pendants labelled by the Abel–Jacobi
class of their attachment (§13.5). B: the push-in demo of
§13.4 — a divisor with chips on the leaves, and its
equivalent core-supported divisor. C: the metric Jacobian
circle with pendants collapsed onto their roots (§13.7).
D: `K` principal — `+1` at the attachments, `−1` at the
leaves, killed by one borrowing per leaf (§13.6).
(The draft PNG's panel C still spaces the vertices by
Euclidean arclength — the pre-correction drawing; the
Phase 3 SVG regeneration redraws it with the equal
lattice spacing of §13.7, and the picture-per-table pass
happens then too, as with earlier chapters.)*

### 13.8 The whole chapter, as a program

*(Pointer in the published chapter, per house rule; the
draft keeps the listing.)* The chapter program is
`hidden_files/code/ch13_chapter_program_draft.py`,
assembled for this chapter from Day 9's two files with
the Baker–Norine machinery inlined: the cubic rebuilt
from its subdivision (balancing asserted); leaf-stripping
to the 2-core and tree roots; Matrix-Tree and Smith
normal form; the push-in lemma on all 19,683 divisors;
the Abel–Jacobi table, its three collisions, the order-6
generator and the pendant group law; `K` principal with
the signed leaf-column witness; Riemann–Roch on 8,953
divisors; and the metric export, with lattice length
`L = 6` asserted alongside the Euclidean total. Phase 2
will re-run it from this draft text and assert the
Try It answers below.

### Try It

1. Using the positions of §13.7, compute
   `S_q(0) + S_q(6)` as a point class: which vertex (or
   vertices) carry it?
2. A spanning tree of `G` must contain all three bridges
   (§13.2). How many edges of the cycle does it contain,
   and why does that force the count 6?
3. Compute `K` at vertices 0 and 4 directly from the
   definition `K = Σ (deg(v) − 2)·(v)`, and check the
   firing script of §13.6 really produces `K` from `0`
   at those two coordinates.
4. True or false, with a reason from the text: "Adding
   a pendant vertex to a graph always leaves its
   Jacobian unchanged — and Abel–Jacobi stays injective."

*Answers (to be asserted explicitly in Phase 2):*

1. `S_q(0)` is position 1 (its attachment 4's class);
   `S_q(6)` is position 4. Sum: position `1 + 4 = 5`,
   the class of vertex **5** — and equally of pendant
   **1**, its tree-mate.
2. A spanning tree on 9 vertices has 8 edges; 3 are the
   forced bridges, leaving 5 cycle edges — a 6-cycle
   minus one edge, 6 choices. Hence 6 spanning trees,
   matching Matrix-Tree.
3. `deg(0) = 1`, so `K(0) = −1`; `deg(4) = 3`, so
   `K(4) = +1`. Borrowing once at leaf 0 subtracts its
   Laplacian column's negative: the script
   `f = (−1,−1,−1,0,…,0)` gives `(Qf)(0) = −1` (only
   leaf 0's own firing touches vertex 0, contributing
   `−deg(0) = −1`) and `(Qf)(4) = +1` (only leaf 0's
   borrowing touches vertex 4 among the leaves, sending
   it one chip). Matches `K` at both coordinates.
4. **False as stated — half of it.** The Jacobian *is*
   unchanged (push-in + Matrix-Tree, §§13.3–13.4), but
   Abel–Jacobi does **not** stay injective: the new
   pendant and its attachment have the same class
   (§13.5). Theorem 1.8's 2-edge-connected hypothesis is
   exactly what a bridge violates.

---

*Next: the Deep Dives — beginning with the defect
ledger: the Kähler package off the uniform locus
(Days 12, 14–15 of the study).*

*Source & provenance (draft): Day 9 (2026-10-19) of the
tropical geometry study — code
`hidden_files/code/day09_full_cubic.py` and
`hidden_files/code/chip_firing.py`, assembled chapter
program `hidden_files/code/ch13_chapter_program_draft.py`,
figure `files/figures/day09-trees-invisible.png`
(generator `hidden_files/code/day09_figures.py`), study
log Day 9 entry, run log `hidden_files/day09-run.log`.
Sources: M. Baker & S. Norine, "Riemann–Roch and
Abel–Jacobi theory on a finite graph",
arXiv:math/0608360 (stabilisation Prop. 3.1, Thms
1.7/1.8/1.12; bridge contraction Cor. 4.7 as cited in
the study notes); Mikhalkin–Zharkov,
arXiv:math/0612267, for the metric Jacobian; the
retraction remark of §13.7 names Baker–Payne–Rabinoff
as the study's flagged source for the Berkovich picture.
Graph data asserted by the program: core
`{3,4,5,6,7,8}`, cycle order `3,4,8,7,6,5`, bridges
`(0,4), (1,5), (2,7)`, `q = 3`. Metric correction from
Chapter 12 §12.8: lattice length `L = 6`, Euclidean total
`4 + 2√2`.*
