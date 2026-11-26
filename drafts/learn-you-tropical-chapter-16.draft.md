# Learn You a Tropical Geometry for Great Good!

## Chapter 16 — The Petersen Graph Has a Jacobian
### Deep Dive 3

> **DRAFT — Phase 1, 2026-11-26.** Prose converted, program
> assembled, figure carried over from Day 17 (PNG; Phase 3
> regenerates as SVG). NOT yet re-verified end-to-end from
> this text — that is Phase 2, with the referee notes.
> Two items are flagged for Phase 2 in the checklist at the
> end: the new §16.4 pairwise-difference claim (computed in
> this draft's assembly run, beyond Day 17's checks) and
> Figure 16.1's text panels.

Chapter 11 built a famous graph and called it a moduli
space. Chapter 12 built a machine — Baker–Norine
chip-firing — and ran it on a hexagon. This chapter points
the machine at the graph:

> **The Petersen graph, read as a tropical curve of genus
> 6, has Jacobian `Jac ≅ ℤ/2 × (ℤ/10)³`, of order 2000 —
> computed three ways that land on the same number, with a
> fourth, spectral, confirmation — and every vertex class
> `[(v) − (q)]` has order exactly 10, the exponent of the
> group. Stronger, and new in this draft: every one of the
> 45 pairwise vertex differences `[(a) − (b)]` has order
> exactly 10 as well. Riemann–Roch holds at genus 6, with
> canonical divisor `K = (1, …, 1)`.**

The Deep Dive door warning applies: Chapter 11 is assumed
for the graph (splits, compatibility, the Petersen
invariants), and Chapters 12–13 are assumed for the
machinery (divisors, the Laplacian, q-reduced
representatives, Abel–Jacobi, Riemann–Roch). This chapter
does not re-teach them; it runs them at a scale where the
answers stop being checkable by eye — 2000 group elements,
none of them listed — and have to be checkable by
structure instead.

Sources: M. Baker & S. Norine, "Riemann–Roch and
Abel–Jacobi theory on a finite graph"
(arXiv:math/0608360) — Jacobian as degree-0 divisors
modulo principal ones, Kirchhoff's theorem
`|Jac| = #spanning trees`, Theorem 1.8 (`S_q¹` injective
iff 2-edge-connected), Theorem 1.12 (Riemann–Roch),
Proposition 3.1 (stabilisation). The graph: D. Speyer &
B. Sturmfels, "The tropical Grassmannian"
(arXiv:math/0304218), Example 3.3, via Chapter 11.

## 16.1 The same graph, a second life

Build the graph exactly as Chapter 11 did. The vertices
are the 2-subsets of `{0, 1, 2, 3, 4}` — the ten splits of
a 5-leaf tree — and two vertices are joined when the
2-subsets are disjoint:

```
V = 10,   E = 15,   every vertex has degree 3.
```

In Chapter 11 those vertices were *rays of a fan* and the
edges were *trees*. Here the graph itself is the curve.
Read as a combinatorial tropical curve — a finite graph,
the skeleton of the tropical moduli space `M_{0,5}` — its
genus is the cycle rank Chapter 11 already computed for
the homotopy count:

```
g = E − V + 1 = 15 − 10 + 1 = 6.
```

**Predict first.** Two predictions, one safe and one
honest about its limits. Safe: Kirchhoff's classical count
says the Petersen graph has **2000 spanning trees**, so if
the Baker–Norine theory of Chapters 12–13 applies,
`|Jac|` should be 2000 `= 2⁴ · 5³`. Not guessed: the
*group structure* — which abelian group of order 2000 —
is whatever the Smith normal form says. Day 17 recorded
the order prediction before running and deliberately did
not predict the decomposition.

One structural fact decides in advance which of
Chapter 12's phenomena survive and which of Chapter 13's
failures cannot happen here. The program checks,
edge by edge, that removing any single edge leaves the
graph connected:

```
bridgeless (2-edge-connected): removing any of the 15
  edges keeps the graph connected.
```

So there are no bridges, no grafted trees, nothing for
Chapter 13's push-in lemma to push: **the whole graph is
its own 2-core.** By Baker–Norine Theorem 1.8, Abel–Jacobi
`S_q¹` should therefore be *injective* — the exact
opposite of Chapter 13, where each pendant vertex
collapsed onto its attachment point. §16.5 collects that
bet.

Fix the base vertex `q = {0,1}` for the whole chapter
(the Day 17 figure draws it as `{1,2}` in 1-based labels;
the program and this text use 0-based splits throughout).

## 16.2 The group: order 2000, three ways — and a fourth

Chapter 12's discipline, at Chapter 12's standard: no
number is believed until independent routes agree.

**Route 1 — Kirchhoff.** The determinant of any reduced
Laplacian (delete `q`'s row and column) counts spanning
trees:

```
det(reduced Laplacian) = 2000.
```

The prediction of §16.1 lands.

**Route 2 — Smith normal form.** The same reduced
Laplacian, over `ℤ`, has Smith diagonal

```
diag(1, 1, 1, 1, 1, 2, 10, 10, 10),     product 2000.
```

> **Honest Box — the diagonal is not the group.** The
> nine diagonal entries are the SNF of a `9 × 9` matrix
> (the reduced Laplacian has one row per vertex ≠ `q`).
> The Jacobian's invariant factors are what remains after
> dropping the entries equal to 1 — a factor `ℤ/1` is the
> trivial group and contributes nothing:
>
> ```
> Jac ≅ ℤ/2 × ℤ/10 × ℤ/10 × ℤ/10,
> ```
>
> with the divisibility chain `2 | 10 | 10 | 10` the
> program asserts, as invariant factors must satisfy.
> Conflating the two lists — reporting a group with five
> extra trivial factors, or worse, reading the order off
> only the nontrivial factors of a *different* matrix —
> is the standard slip here, so the draft keeps both
> lists in view and asserts each against the other.

**Route 3 — G-parking functions.** A q-reduced degree-0
divisor is determined by its values off `q`, and
Chapter 12's singleton-set argument bounds them: P1 gives
`D(v) ≥ 0`, and P2 applied to `A = {v}` gives
`D(v) < outdeg_A(v) = deg(v) = 3`, so `D(v) ∈ {0, 1, 2}`.
Enumerate all `3⁹ = 19,683` candidates, keep those
passing Dhar's full burning condition (P2 for every
`A ⊆ V ∖ {q}`), and set `D(q)` to minus the sum:

```
G-parking functions (q-reduced degree-0 reps) = 2000.
```

One representative per class, 2000 classes — the count
agrees, and this route never computes a determinant: it
enumerates reduced divisors directly.

**Route 4 — the spectrum (a confirmation, honestly
labelled).** The full Laplacian's eigenvalues, computed
exactly (sympy, integer arithmetic):

```
spectrum = {0: ×1,   2: ×5,   5: ×4}.
```

Kirchhoff's theorem in spectral form says
`|Jac| = (1/V) · Π (nonzero eigenvalues)`, and indeed
`(2⁵ · 5⁴)/10 = 20000/10 = 2000`. This is a genuine
fourth *computation* — eigenvalues, not a determinant of
the reduced matrix — but it is not a fourth *theorem*:
routes 1, 2 and 4 all read the Laplacian, by three
different algorithms. Route 3 is the structurally
independent count. The draft says so rather than
inflating "three ways" into "four ways": the headline
remains the three routes of the outline — determinant,
normal form, parking functions — with the spectrum as a
corroborating aside the Phase 1 assembly added.

```
|Jac| = 2000:  det  =  SNF product  =  parking count
               =  (2⁵ · 5⁴)/10        (spectral check)
Jac ≅ ℤ/2 × (ℤ/10)³,   exponent 10.
```

## 16.3 Every vertex class has order exactly 10

The exponent of `ℤ/2 × (ℤ/10)³` is 10: no element has
order larger than 10, and some element attains it. The
chapter's sharpest fact is *which* elements attain it:
all of the obvious ones.

```
orders of [(v) − (q)], v ≠ q:
  all nine classes have order exactly 10.
```

**How much of this is symmetry?** Less than it first
looks, and the draft is careful about the difference.
Permuting the five leaves (`S₅`) is a symmetry of the
graph, but the class `[(v) − (q)]` is defined relative to
the *fixed* base point `q`, so only the stabiliser of
`q` acts. Its orbits off `q` are two (asserted in the
program): the **3 splits disjoint from `q`** — the
2-subsets of `{2,3,4}` — and the **6 splits meeting `q`
in exactly one leaf**. Symmetry therefore forces *at
most two* order values among the nine classes, one per
orbit. That the two values coincide at 10 is a
computation, not a symmetry argument — and Chapter 18's
quartic skeleton is the control showing the distinction
matters: there the vertex classes take orders
`{7, 14, 28}`, plural, exactly because no such coincidence
occurs. (That chapter is Deep Dive 4; the contrast was
already visible on Day 18.)

**The draft's strengthening.** Day 17 checked the nine
classes above. This draft's assembled program checks
every pair: for all `45 = C(10,2)` unordered pairs
`{a, b}`,

```
[(a) − (b)]  has order exactly 10.
```

Same orbit bookkeeping applies — unordered pairs fall
into two `S₅`-orbits, adjacent (disjoint splits, distance
1) and non-adjacent (meeting splits, distance 2) — so
symmetry again reduces the work to two cases, and the
program verifies both, pair by pair. The statement is
flagged here and in the checklist: it is *new in
Phase 1*, computed in the assembly run
(`hidden_files/ch16-draft-run.log`), and Phase 2 must
re-verify it from this text like every other claim
before the chapter ships.

Why it is worth a section and not a footnote: a Jacobian
of order 2000 whose every difference of two points sits
at the *maximum possible order* is rigid in a way the
bare group decomposition does not state. None of the
point-differences collapses early — there is no pair of
vertices whose difference is torsion of order 1, 2 or 5.
The points of this curve, as elements of its Jacobian,
are as far from each other as the group allows.

## 16.4 Abel–Jacobi, the group law, and no collapse

The Abel–Jacobi classes are the ten divisors
`S_q(v) = [(v) − (q)]`, each represented by its unique
q-reduced representative (Baker–Norine Proposition 3.1,
the stabiliser of Chapters 12–13, reused unchanged):

```
the 10 classes S_q(v) are pairwise distinct.
```

Injectivity is §16.1's bet collected: the graph is
bridgeless, Theorem 1.8 applies, and the computation
agrees. Set this beside Chapter 13 once more, because
the pair of chapters is the theorem's two sides: graft
trees onto a cycle and Abel–Jacobi collapses each tree
to its attachment point (the Jacobian ignores the trees);
take a graph with no trees at all and nothing collapses
(the Jacobian sees every vertex separately).

The group law on representatives is addition followed by
re-reduction, and the program checks the defining
identity in the form Chapter 12 used:

```
[(a) − (q)] + [(b) − (q)] = [(a) + (b) − 2(q)]
  as q-reduced representatives.
```

Nothing about the law is special to this graph — that is
the point. The *same* code that added points on a hexagon
in Chapter 12 adds them here, on a genus-6 graph whose
group has 2000 elements; only the graph data changed.

## 16.5 Riemann–Roch at genus 6

The canonical divisor is `K = Σ_v (deg(v) − 2)·(v)`.
Three-regularity makes it the flattest possible
nonzero `K` in the book:

```
K = (1, 1, …, 1),   deg K = 10 = 2g − 2.
```

Hand checks first, as always:

```
r(0) = 0,   r(K) = g − 1 = 5,   r((q)) = 0.
```

Then the theorem itself — Baker–Norine Theorem 1.12,
`r(D) − r(K − D) = deg D − g + 1` — on random divisors,
rank computed from the definition exactly as in
Chapter 12 (via q-reduced representatives):

```
RR identity verified on 40 random divisors
  (entries in {−1, 0, 1}, degrees 0–3; seed 17).
```

Forty is a sample, not an exhaustion — at genus 6 the
Chapter 12 move (all `3⁶` divisors of a box) does not
scale to a meaningful box here, and the draft does not
pretend otherwise. What carries the weight is that the
rank function, the canonical divisor and the identity
are the same objects as in the genus-1 chapter, and the
identity holds on every divisor the program drew.

A stabilisation exhibit, for the figure: the divisor
`[4, −2, 3, 0, −1, 2, 0, 1, −3, 2]` (in the vertex order
of the program) stabilises under Proposition 3.1 to the
q-reduced `[0, 2, 0, 1, 1, 1, 1, 0, 0, 0]` in 10 fired
sets, degree preserved. Figure 16.1C draws the
before → after chip counts on the graph.

![Chip-firing on the Petersen graph: the labelled graph with base point, the group of order 2000, a stabilisation exhibit, and the Riemann–Roch data](figures/day17-petersen-jacobian.png)

*Figure 16.1 — A: the Petersen graph in the classic
pentagon/pentagram drawing, vertices labelled by their
splits (1-based in the drawing; §16.1's `q` is the red
`{1,2}`), genus `g = 6`. B: the group
`ℤ/2 × (ℤ/10)³` and its three routes to 2000 (§16.2).
C: the stabilisation exhibit of §16.5, chips
before → q-reduced after. D: Riemann–Roch data
(`K = (1,…,1)`, hand checks, 40 divisors). (Day 17 PNG;
Phase 3 regenerates as SVG. Two honest flags for that
regeneration, recorded now: panels B and D are text
receipts, not drawings — under the book's
figure house rule the SVG pass must give the drawn
counterparts (the three counts and the order-10
uniformity, drawn) or demote the panels from the figure;
and panel A's inner pentagram is a straight-line drawing
of the asserted edges, not the classic pentagram
geometry — Day 17 noted this rather than hiding it.)*

## 16.6 What the Jacobian of a moduli skeleton means

It is worth one section of scope discipline, because
the headline invites a slide. Two Jacobians are in the
room, and the chapter computes one of them.

*Not* computed: the Jacobian of the classical moduli
space `M_{0,5}`, or of any algebraic curve. *Computed:*
the Baker–Norine Jacobian of the finite graph that
Chapter 11 identified as the tropical Grassmannian's
skeleton — the combinatorial shadow of `M_{0,5}`'s
Berkovich skeleton. Metrically, a genus-6 tropical curve
has a Jacobian that is a 6-dimensional real torus; the
finite group of this chapter is what that torus's
shadow looks like when the metric is forgotten and only
the graph remains — exactly as Chapter 12's `ℤ/6` is the
combinatorial shadow of its circle `ℝ/6ℤ`.

That is the Deep Dive's real subject, and it is the
study's two oldest threads meeting: the graph arrived
as *moduli* (Part V) and the machinery arrived as
*divisor theory* (Part VI), and neither needed the other
to be re-derived. Feed the moduli skeleton into the
chip-firing machine and a clean finite group falls out:
`ℤ/2 × (ℤ/10)³`, order 2000, every point-difference at
full order. The retraction story that would make
"skeleton of `M_{0,5}`" more than a name — Baker–Payne–
Rabinoff — remains on the open list (Part VIII); this
chapter claims the graph, its group, and nothing
upstream of them.

## 16.7 The whole chapter, as a program

The draft program is
`hidden_files/code/ch16_chapter_program_draft.py`
(Day 17's program, assembled for this chapter, plus the
spectral route of §16.2, the orbit count of §16.3, and
the 45-pair order check of §16.3; it imports the
`ChipFiring` machinery from `chip_firing.py` in the same
directory, as Day 17 did — the Phase 3 bundle program
will be merged self-contained, as Chapters 14–15 were).
It asserts: the graph census `(V, E) = (10, 15)`,
3-regularity, genus 6; bridgelessness edge by edge;
Kirchhoff `det = 2000`; the SNF diagonal
`[1,1,1,1,1,2,10,10,10]`, its product, the invariant
factors `[2,10,10,10]` and their divisibility chain;
the exact Laplacian spectrum `{0:1, 2:5, 5:4}` and the
spectral Kirchhoff count; 2000 G-parking functions over
the `3⁹` candidates; Abel–Jacobi injectivity (10 distinct
classes); the `3 + 6` orbit split off `q`; order exactly
10 for all nine vertex classes and for all 45 pairwise
differences; the group-law identity; `K = (1,…,1)`,
`deg K = 10`, the three hand ranks; Riemann–Roch on 40
random divisors; and the stabilisation exhibit, target
included. Assembly run:
`hidden_files/ch16-draft-run.log` — all assertions
passed (≈40 s). NOT yet re-verified from this draft
text — that is Phase 2, per the cycle.

## 16.8 Try It

1. **Genus and `K` by hand.** From `V = 10`, `E = 15`,
   compute `g`, and then `K` for a 3-regular graph.
   Check `deg K = 2g − 2`.
2. **Read the Smith normal form.** The reduced
   Laplacian's SNF diagonal is
   `[1, 1, 1, 1, 1, 2, 10, 10, 10]`. Write down
   `Jac` as a product of cyclic groups, its order, and
   its exponent. Which entries did you drop, and why is
   dropping them legitimate?
3. **The parking bound.** Using only P1 and P2 with
   `A = {v}`, show a q-reduced degree-0 divisor on this
   graph has `D(v) ∈ {0, 1, 2}` at every `v ≠ q`. How
   many candidates does that leave, and how many pass?
4. **Two orbits.** With `q = {0,1}`, list the splits
   disjoint from `q` and count the splits meeting `q` in
   one element. Symmetry under the stabiliser of `q`
   then forces at most how many distinct values among
   the nine orders of §16.3? What settles the actual
   value(s)?

## Answers (to be asserted in Phase 2; items 3–4's
## counts are already asserted in the draft program)

- **1.** `g = 15 − 10 + 1 = 6`. `K(v) = deg(v) − 2 = 1`
  at every vertex, so `K = (1,…,1)` and
  `deg K = 10 = 2·6 − 2`. ✓
- **2.** Drop the five `1`s (`ℤ/1` is trivial):
  `Jac ≅ ℤ/2 × (ℤ/10)³`; order `2 · 10³ = 2000`;
  exponent `lcm(2, 10) = 10` — the largest invariant
  factor, since `2 | 10 | 10 | 10`.
- **3.** P1: `D(v) ≥ 0` off `q`. P2 on `A = {v}` needs
  some vertex of `A` with `D(v) < outdeg_A(v)`; `A` has
  only `v`, and all 3 of its edges leave `A`, so
  `outdeg_A(v) = 3` and `D(v) ≤ 2`. Candidates:
  `3⁹ = 19,683`; exactly 2000 pass the full P2 and are
  the q-reduced degree-0 representatives.
- **4.** Disjoint from `{0,1}`: the 2-subsets of
  `{2,3,4}` — `{2,3}, {2,4}, {3,4}`, three of them.
  Meeting `q` in one element: choose the shared element
  (2 ways) and the other element from the 3 outside —
  `2 · 3 = 6`. So at most **two** distinct order values
  are forced by symmetry; that both equal 10 is settled
  by the computation of §16.3, not by the symmetry.

## What we learned

- **A moduli skeleton has a divisor theory, and it is
  computable end-to-end.** The Petersen graph of
  Chapter 11, fed to Chapters 12–13's machinery
  unchanged, has `Jac ≅ ℤ/2 × (ℤ/10)³` — order 2000 by
  determinant, by Smith normal form, and by an
  enumeration of reduced divisors, with the Laplacian
  spectrum `{0, 2⁵, 5⁴}` as a fourth computation.
- **Order-10 uniformity is a computation, not a
  symmetry.** Stabiliser orbits reduce the vertex
  classes to two cases; both land at the exponent — and
  so does every one of the 45 pairwise differences.
  No point-difference on this curve is low-order
  torsion.
- **Bridgelessness is the whole Abel–Jacobi story.**
  Chapter 13's collapse and this chapter's injectivity
  are one theorem (Baker–Norine 1.8) read from its two
  sides: trees collapse, 2-cores do not.

## Next

Deep Dive 4 — quartic, quintic, and Dhar: three hexagons
sharing edges at genus 3 (`Jac = ℤ/7 × ℤ/28`, vertex
orders *not* uniform — §16.3's control), six hexagons at
genus 6 (`Jac = ℤ/82 × ℤ/410`), and Dhar's burning
algorithm carrying Riemann–Roch where the subset test
cannot (Days 18–20).

---

*Source & provenance (draft): Day 17 (2026-10-27) of the
tropical geometry study — code
`hidden_files/code/day17_petersen_jacobian.py` and
`day17_figures.py`, summary `day17_summary.json`, study
log entry for that day, run log
`hidden_files/day17-run.log`. Machinery:
`hidden_files/code/chip_firing.py` (Baker–Norine,
parameterised on Day 9 from Day 8's verified
implementation). Graph construction as in Chapter 11
(Speyer–Sturmfels, arXiv:math/0304218, Ex. 3.3).
Preflight 2026-11-25 re-ran Day 17 unmodified, clean
(`hidden_files/ch16-preflight-day17.log`), with an
independent sympy cross-check of the determinant, the
SNF diagonal, and the spectrum. Phase 1 assembly run
2026-11-26: `hidden_files/ch16-draft-run.log`, all
assertions passed, including the three draft additions
(spectral route, orbit split, 45 pairwise orders).*

## Phase 2 checklist (editorial — removed before shipping)

- [ ] Re-run the chapter program from this text, claim
      by claim; assert the four Try It answers in a
      companion verify program (Ch. 15 pattern). The
      §16.3 pairwise claim (all 45 differences of order
      exactly 10) is new in Phase 1 — it gets the full
      re-verification, not a free pass from the assembly
      run.
- [ ] Make the bundle program self-contained (inline
      `ChipFiring`, as Chs. 14–15's bundles inlined
      their machinery). PHASE 3, as planned.
- [ ] Figure 16.1: panels B and D are text receipts.
      Phase 3 SVG regeneration must draw the
      picture-per-table counterparts (three routes to
      2000; the order-10 uniformity) or demote the
      panels — house rule. Also unify the figure's
      1-based split labels with the text's 0-based
      convention, or state the offset on the figure.
- [ ] Check the §16.2 spectral aside's framing against
      the referee's independence standard: routes 1/2/4
      share the Laplacian; route 3 does not. Text
      already says so — referee to confirm the wording
      does not overclaim "four ways" anywhere.
- [ ] Referee notes from #minplus-book: check each, fix
      what holds, reply to what doesn't.
- [ ] Phase 3: SVG figure regeneration with asserted
      values, compression pass, ship.
