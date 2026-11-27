# Learn You a Tropical Geometry for Great Good!

## Chapter 16 — The Petersen Graph Has a Jacobian
### Deep Dive 3

> **DRAFT — Phase 2, 2026-11-27.** Re-verified end-to-end
> from this text (chapter program, independent verify, and
> an exhaustive class-level Riemann–Roch check — all
> assertions passed) and revised against the referee pass
> in #minplus-book (Nyx: E1–E5, G1–G8, N1–N11; ClaudeDuet:
> concurrence + D1–D4; every point checked, every fix
> applied — see the Phase 2 record at the end). NOT yet
> polished or shipped: that is Phase 3. Figure 16.1 is
> still the Day 17 PNG; its SVG regeneration is a Phase 3
> item.

Chapter 11 built a famous graph and called it a moduli
space. Chapter 12 built a machine — Baker–Norine
chip-firing — and ran it on a hexagon. This chapter points
the machine at the graph:

> **The Petersen graph, read as a tropical curve of genus
> 6 (unit edge lengths), has Jacobian `Jac ≅ ℤ/2 × (ℤ/10)³`,
> of order 2000 — the order computed by two independent
> routes (Kirchhoff's determinant; an enumeration of
> G-parking functions), the structure read off the Smith
> normal form, the Laplacian spectrum as corroboration —
> and every vertex class `[(v) − (q)]` has order exactly
> 10, the exponent of the group. By vertex-transitivity
> the same holds for every pairwise difference
> `[(a) − (b)]`: a corollary of the nine classes, verified
> pair by pair in code as a cross-check. Riemann–Roch
> holds at genus 6 — exhaustively, over all 26,000
> divisor classes in degrees −1..11 — with canonical
> divisor `K = (1, …, 1)`.**

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
(arXiv:math/0304218), Example 3.3, via Chapter 11. The
flow-lattice reading of §16.6: R. Bacher, P. de la Harpe
& T. Nagnibeda, "The lattice of integral flows and the
lattice of integral cuts on a finite graph", Bull. Soc.
Math. France 125 (1997).

## 16.1 The same graph, a second life

Build the graph exactly as Chapter 11 did. The vertices
are the 2-subsets of `{0, 1, 2, 3, 4}` — the ten
nontrivial splits of 5 leaves — and two vertices are
joined when the 2-subsets are disjoint (so the edges are
the 15 trivalent 5-leaf trees of Chapter 11, one per
compatible pair of splits):

```
V = 10,   E = 15,   every vertex has degree 3.
```

In Chapter 11 those vertices were *rays of a fan* and the
edges were *trees*. Here the graph itself is the curve.
Read as a combinatorial tropical curve — the graph with
unit edge lengths, a metric graph — it is the skeleton
Chapter 11 met from the other side: precisely, the link
of the tropical Grassmannian fan, equivalently the
boundary complex of the compactification `M̄_{0,5}`. Its
genus is the cycle rank Chapter 11 already computed for
the homotopy count:

```
g = E − V + 1 = 15 − 10 + 1 = 6.
```

**Predict first.** Two predictions, one safe and one
honest about its limits. Safe: Kirchhoff's classical count
says the Petersen graph has **2000 spanning trees**, so if
the Baker–Norine theory of Chapters 12–13 applies,
`|Jac|` should be 2000. Not guessed: the
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

So there are no bridges. That is stronger than "no
pendant trees", and the difference matters: two cycles
joined by a single bridge have no pendant trees either —
such a graph is even its own 2-core — yet Abel–Jacobi
still identifies the bridge's endpoints (§16.4).
Baker–Norine Theorem 1.8's criterion is
2-edge-connectedness itself, and this graph meets it, so
Abel–Jacobi `S_q¹` should be *injective* — the exact
opposite of Chapter 13, where each pendant vertex
collapsed onto its attachment point. §16.4 collects that
bet.

Fix the base vertex `q = {0,1}` for the whole chapter
(the Day 17 figure draws it as `{1,2}` in 1-based labels;
the program and this text use 0-based splits throughout).

## 16.2 The group: order 2000 two independent ways, structure from Smith

Chapter 12's discipline, at Chapter 12's standard: no
number is believed until independent routes agree — and
"independent" is meant strictly. One presentation runs
under three of the four routes below, so state it once:

> **The presentation.** Dropping the `q`-coordinate
> identifies `Div⁰` with `ℤ^{V∖q} ≅ ℤ⁹` (the dropped
> coordinate is minus the sum of the rest). A principal
> divisor `Lf`, read in the gauge `f(q) = 0` — every `f`
> has exactly one such representative, since the constants
> are exactly `ker L` on a connected graph — becomes
> `L_red · f|_{V∖q}`. Hence
> `Jac = Div⁰/Prin ≅ ℤ⁹ / L_red ℤ⁹ = coker(L_red)`.
> Routes 1, 2 and 4 are three readings of this one
> presentation; Route 3 never touches it.

**Route 1 — Kirchhoff (the order).** The determinant of
the reduced Laplacian counts spanning trees:

```
det(reduced Laplacian) = 2000     (exact integer det).
```

The prediction of §16.1 lands. (The program computes this
determinant exactly, over `ℤ`; the shared `ChipFiring`
route's float `round(abs(det))` agrees, and is kept only
as a cross-check.)

**Route 2 — Smith normal form (the structure).** The same
reduced Laplacian, over `ℤ`, has Smith diagonal

```
diag(1, 1, 1, 1, 1, 2, 10, 10, 10),     product 2000.
```

That the product equals Route 1's determinant is
automatic — Smith operations are unimodular, so
`∏ SNF = |det|` by construction; their agreement checks
the code, not the mathematics. The route's real content
is the *structure*:

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
> extra trivial factors — is one standard slip; the other
> is reading the group off the *full* `10 × 10`
> Laplacian, whose Smith diagonal is
> `[1, 1, 1, 1, 1, 2, 10, 10, 10, 0]` (asserted in the
> program): `coker L ≅ ℤ ⊕ Jac`, the trailing `0` being a
> free `ℤ` summand — the constants — not a torsion
> factor. The reduced matrix exists precisely to quotient
> that `ℤ` away.

**Route 3 — G-parking functions (the independent order
count).** A q-reduced degree-0 divisor is determined by
its values off `q`, and Chapter 12's singleton-set
argument bounds them: P1 gives `D(v) ≥ 0`, and P2 applied
to `A = {v}` gives `D(v) < outdeg_A(v) = deg(v) = 3`, so
`D(v) ∈ {0, 1, 2}`. Enumerate all `3⁹ = 19,683`
candidates, keep those passing Dhar's full burning
condition (P2 for every `A ⊆ V ∖ {q}`), and set `D(q)` to
minus the sum:

```
G-parking functions (q-reduced degree-0 reps) = 2000.
```

By Baker–Norine Proposition 3.1 there is exactly one
q-reduced representative per class, so this count *is*
`|Jac|` — and the program additionally checks, via the
lattice keys of the presentation above, that the 2000
representatives land in 2000 distinct classes. This route
never computes a determinant: it enumerates reduced
divisors directly. It is the second independent route to
the order.

**Route 4 — the spectrum (a corroboration).** The full
Laplacian's eigenvalues, computed exactly (sympy, integer
arithmetic):

```
spectrum = {0: ×1,   2: ×5,   5: ×4}.
```

Kirchhoff's theorem in spectral form says
`|Jac| = (1/V) · Π (nonzero eigenvalues)`, and indeed
`(2⁵ · 5⁴)/10 = 20000/10 = 2000`. This is Matrix–Tree in
eigen-form — a different algorithm reading the same
Laplacian, so it corroborates the computation without
adding an independent theorem.

> **Honest Box — Route 4 on paper.** The spectrum need
> not be a black box. The Petersen graph is strongly
> regular with parameters `srg(10, 3, 0, 1)`: two adjacent
> vertices have 0 common neighbours, two non-adjacent
> vertices have exactly 1 (both asserted in the program).
> Strong regularity gives `A² = 2I − A + J` for the
> adjacency matrix (also asserted, entrywise). On the
> all-ones direction `A` acts as 3; on its orthogonal
> complement `J` acts as 0, so the remaining adjacency
> eigenvalues satisfy `λ² + λ − 2 = 0`, i.e.
> `λ ∈ {1, −2}`; their multiplicities come from
> `trace A = 0` and `trace A² = Σ deg = 30`:
> `3 + f − 2g′ = 0` with `f + g′ = 9` gives `f = 5`,
> `g′ = 4`. Since `L = 3I − A`, the Laplacian spectrum is
> `{0: ×1, 2: ×5, 5: ×4}` — no eigenvalue solver required;
> sympy merely confirms it.

```
|Jac| = 2000:  det (exact)  =  parking count
               (the two independent routes)
               SNF product = |det| automatically;
               its content is the structure
               Jac ≅ ℤ/2 × (ℤ/10)³,   exponent 10.
               spectrum: (2⁵ · 5⁴)/10 = 2000
               (corroboration)
```

## 16.3 Every vertex class has order exactly 10

The exponent of `ℤ/2 × (ℤ/10)³` is 10: no element has
order larger than 10, and some element attains it. The
chapter's central computation is *which* elements attain
it: all of the obvious ones.

```
orders of [(v) − (q)], v ≠ q:
  all nine classes have order exactly 10.
```

**How much of this is symmetry?** Less than it first
looks, and the text is careful about the difference.
Permuting the five leaves (`S₅`) is a symmetry of the
graph, but the class `[(v) − (q)]` is defined relative to
the *fixed* base point `q`, so only the stabiliser of
`q` acts. It acts on `Jac` by automorphisms,
`[(v) − (q)] ↦ [(σv) − (q)]` for `σ ∈ Stab(q)`, so order
is constant on each of its orbits. Its orbits off `q`
are two — computed by brute force over the stabiliser's
12 elements in the program: the **3 splits disjoint from
`q`** — the 2-subsets of `{2,3,4}`, which are exactly the
neighbours `N(q)` of `q` — and the **6 splits meeting `q`
in exactly one leaf**. Symmetry therefore forces *at
most two* order values among the nine classes, one per
orbit. That the two values coincide at 10 is a
computation, not a symmetry argument.

**The 45 pairs are a corollary, not a strengthening.**
Phase 1 of this draft advertised the following as a
stronger, new claim; the referee pass (Nyx, E1) pointed
out it is equivalent to the nine classes, and the
correction stands. The graph's automorphism group is
vertex-transitive: for any pair `{a, b}` pick `σ` with
`σ(b) = q`; `σ` induces an automorphism of `Jac` sending
`[(a) − (b)]` to `[(σa) − (q)]`, one of the nine classes.
Hence

```
[(a) − (b)]  has order exactly 10
  for all 45 = C(10,2) unordered pairs {a, b},
```

with no new computation owed. (Unordered pairs are
legitimate because `ord(x) = ord(−x)` and
`[(b) − (a)] = −[(a) − (b)]`. The pairs themselves fall
into two `S₅`-orbits — 15 adjacent, 30 non-adjacent —
which is *distance-transitivity* of the Petersen graph,
a stronger property than the bare `S₅` action; the
program confirms both orbits by brute force over all 120
permutations.) The 45-pair loop stays in the program as a
cross-check of the corollary, and Try It 5 asks for the
argument by hand.

**How surprising is the uniformity?** Less than the
Phase 1 draft claimed, and the numbers say by how much.
The order distribution of the whole group (asserted in
the program, from the lattice presentation) is

```
order:    1    2     5     10
count:    1   15   124   1860
```

— 93% of `Jac`'s elements have order 10. With only two
symmetry-independent vertex classes, the naive prior for
both landing at 10 is already `≈ 0.93² ≈ 0.86`. So the
honest content is narrower than "rigidity", and cleaner:
*no vertex difference is low-order torsion* — none of
the nine classes (equivalently, none of the 45
differences) lies in `Jac[2] ≅ (ℤ/2)⁴` or in
`Jac[5] ≅ (ℤ/5)³`, the subgroups of elements killed by 2
and by 5. (Order is not a distance, and the text no
longer calls it one.)

Chapter 17's quartic skeleton (Deep Dive 4) is the
*contrast*: its 13-vertex 2-core has
`Jac = ℤ/7 × ℤ/28` and vertex-class orders
`{7, 14, 28}`, plural. It is a contrast, not a control —
that graph is not vertex-transitive, so its plural orders
show that non-uniformity *occurs*; they cannot show that
symmetry alone didn't force uniformity here. What forces
the question down to two cases here is `Stab(q)`; what
settles both cases is the computation.

## 16.4 Abel–Jacobi, the group law, and no collapse

The Abel–Jacobi map is Chapter 12's `S_q¹`; its value at
`v` is the class `S_q(v) = [(v) − (q)]`, each represented
by its unique q-reduced representative (Baker–Norine
Proposition 3.1 — the stabilisation algorithm of
Chapters 12–13, reused unchanged):

```
the 10 classes S_q(v) are pairwise distinct.
```

Injectivity is §16.1's bet collected: the graph is
2-edge-connected, Theorem 1.8 applies, and the
computation agrees. In fact §16.3 already implies it:
`S_q(a) = S_q(b)` would make `[(a) − (b)]` the identity,
but every pairwise difference has order 10. Set this
beside Chapter 13 once more, because the pair of chapters
is the theorem's two sides: graft trees onto a cycle and
Abel–Jacobi collapses each tree to its attachment point
(the Jacobian ignores the trees); join two cycles by a
single bridge — a graph with no pendant trees at all,
indeed its own 2-core — and Abel–Jacobi *still* collapses
the bridge's endpoints onto each other (fire one side of
the bridge: `(u) − (v) = L·1_A` is principal; the Phase 2
verification checks this on two triangles joined by a
bridge). The criterion is 2-edge-connectedness — no
bridges anywhere — and this graph meets it, so nothing
collapses: the Jacobian sees every vertex separately.

The group law on representatives is addition followed by
re-reduction, and the program checks the defining
identity in the form Chapter 12 used:

```
[(a) − (q)] + [(b) − (q)] = [(a) + (b) − 2(q)]
  as q-reduced representatives.
```

Honestly labelled (Duet, D2): that identity holds for
*any* well-defined reduction map, so it tests the
reduction code, not the law. The law itself is exercised
separately — the program checks addition of classes
associative on 100 random triples drawn from the 2000
classes. Nothing about the law is special to this
graph — that is the point. The *same* code that added
points on a hexagon in Chapter 12 adds them here, on a
genus-6 graph whose group has 2000 elements; only the
graph data changed.

## 16.5 Riemann–Roch at genus 6

The canonical divisor is `K = Σ_v (deg(v) − 2)·(v)`.
Three-regularity makes it the flattest possible
nonzero `K` in the book:

```
K = (1, 1, …, 1),   deg K = 10 = 2g − 2.
```

Direct computations first (the Phase 1 draft called them
hand checks; they are computations from the definition —
and two of the three have one-line proofs):

```
r(0) = 0,   r(K) = g − 1 = 5,   r((q)) = 0.
```

`r(0) = 0` is trivial: `r(0) ≥ 1` would force the
degree-(−1) class of `−(v)` to have an effective
representative. `r((q)) = 0`: if `r((q)) ≥ 1`, then
`(q) − (v)` would be equivalent to an effective divisor
of degree 0 — i.e. to 0 — making `(q) ~ (v)` and
contradicting §16.4's injectivity. `r(K) = g − 1` is a
*consequence* of Riemann–Roch itself (apply the theorem
to `D = 0`); the program computes it directly from the
definition, so its agreement here is a genuine check,
not the theorem assuming itself.

Then the theorem itself — Baker–Norine Theorem 1.12,
`r(D) − r(K − D) = deg D − g + 1` — **exhaustively**.
Rank depends only on the class of `D` (equivalent
divisors have the same rank, Chapter 12), so the right
enumeration is over classes, not divisors: each degree-`d`
class is a degree-0 class plus `d·(q)`, and there are
2000 classes in each degree. The program computes the
rank of every class in every degree from −1 to 11 —
26,000 class representatives, organised by the lattice
keys of §16.2's presentation — and cross-checks those
ranks against Chapter 12's rank function on the old
40-divisor sample (which survives only as that
cross-check). The identity holds for all 26,000. The
Phase 1 draft's sample reached only degrees 0–3; the
exhaustion buys facts no sample there could state:

```
divisorial gonality = 4:  no divisor of degree ≤ 3 has
  r ≥ 1, and degree 4 does.
max r by degree:  4→1, 5→1, 6→2, 7→2, 8→3
  (then 9→4, 10→5, 11→5 — the tail Riemann–Roch forces
  once r(K) = 5 is computed).
```

A stabilisation exhibit, for the figure: the divisor
`[4, −2, 3, 0, −1, 2, 0, 1, −3, 2]` (in the vertex order
of the program) stabilises under Proposition 3.1 to the
q-reduced `[0, 2, 0, 1, 1, 1, 1, 0, 0, 0]`, degree
preserved. The program's phase 2 reaches that target in
10 set-firings under lex subset order (phase-1 debt
clearing not counted) — a count that is
algorithm-dependent, not content: reversed subset order
reaches the same target in 3 set-firings. The target is
the exhibit; the count is a footnote. Figure 16.1C draws
the before → after chip counts on the graph.

![Chip-firing on the Petersen graph: the labelled graph with base point, the group of order 2000, a stabilisation exhibit, and the Riemann–Roch data](figures/day17-petersen-jacobian.png)

*Figure 16.1 — A: the Petersen graph in the classic
pentagon/pentagram drawing, vertices labelled by their
splits (1-based in the drawing; §16.1's `q` is the red
`{1,2}`), genus `g = 6`. B: the group
`ℤ/2 × (ℤ/10)³` and its routes to 2000 (§16.2: order by
determinant and by parking enumeration, structure by
Smith normal form, spectrum corroborating).
C: the stabilisation exhibit of §16.5, chips
before → q-reduced after. D: Riemann–Roch data
(`K = (1,…,1)`, the direct computations, the exhaustive
26,000-class check). (Day 17 PNG;
Phase 3 regenerates as SVG. Two honest flags for that
regeneration, recorded now: panels B and D are text
receipts, not drawings — under the book's
figure house rule the SVG pass must give the drawn
counterparts (the routes to 2000 and the order-10
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
Chapter 11 identified as the link of the tropical
Grassmannian fan — the boundary complex of `M̄_{0,5}` —
the combinatorial shadow of `M_{0,5}`'s Berkovich
skeleton. Metrically, give the graph unit edge lengths,
making it a metric graph `Γ` of genus 6: the tropical
Jacobian `Pic⁰(Γ)` is the real torus
`H₁(Γ, ℝ)/H₁(Γ, ℤ) ≅ ℝ⁶/ℤ⁶`, and the finite group of
this chapter sits inside it as the subgroup of
vertex-supported divisor classes. Equivalently, with
`Λ = H₁(G, ℤ)` realised as the lattice of integral flows,
carrying the pairing from the unit edge lengths,
`Jac(G) ≅ Λ^#/Λ` — the discriminant group of the flow
lattice (Bacher–de la Harpe–Nagnibeda 1997) — exactly as
Chapter 12's `ℤ/6` is the combinatorial shadow of its
circle `ℝ/6ℤ`.

That is the Deep Dive's real subject, and it is the
study's two oldest threads meeting: the graph arrived
as *moduli* (Part V) and the machinery arrived as
*divisor theory* (Part VI), and neither needed the other
to be re-derived. Feed the moduli skeleton into the
chip-firing machine and a clean finite group falls out:
`ℤ/2 × (ℤ/10)³`, order 2000, every vertex class — and by
vertex-transitivity every point-difference — at full
order 10, with no low-order torsion among them. The
retraction story that would make "skeleton of `M_{0,5}`"
more than a name — Baker–Payne–Rabinoff — remains on the
open list (Part VIII); this chapter claims the graph, its
group, and nothing upstream of them.

## 16.7 The whole chapter, as a program

The chapter program is
`hidden_files/code/ch16_chapter_program_draft.py`
(Day 17's program, assembled for this chapter and
extended in Phase 2 per the referee pass; it imports the
`ChipFiring` machinery from `chip_firing.py` in the same
directory, as Day 17 did — the Phase 3 bundle program
will be merged self-contained, as Chapters 14–15 were).
It asserts: the graph census `(V, E) = (10, 15)`,
3-regularity, genus 6; 2-edge-connectedness edge by edge;
Kirchhoff `det = 2000` by exact integer determinant (the
float route cross-checked); the SNF diagonal
`[1,1,1,1,1,2,10,10,10]`, its product, the invariant
factors `[2,10,10,10]` and their divisibility chain; the
full-Laplacian SNF `[1,1,1,1,1,2,10,10,10,0]`; the exact
Laplacian spectrum `{0: ×1, 2: ×5, 5: ×4}` and the
spectral Kirchhoff count, plus the on-paper route
(strong regularity `srg(10,3,0,1)`, `A² + A − 2I = J`,
adjacency spectrum `{3: ×1, 1: ×5, −2: ×4}`); 2000
G-parking functions over the `3⁹` candidates, landing in
2000 distinct lattice classes; Abel–Jacobi injectivity
(10 distinct classes); the stabiliser of `q` (order 12)
and its brute-force orbits off `q`, `3 + 6`, the disjoint
orbit being exactly `N(q)`; the `S₅` pair orbits
`15 + 30`; order exactly 10 for all nine vertex classes
and — as the vertex-transitivity corollary's cross-check —
for all 45 pairwise differences; the order distribution
`{1: 1, 2: 15, 5: 124, 10: 1860}` with
`|Jac[2]| = 16`, `|Jac[5]| = 125`; the group-law identity
and associativity on 100 random triples; `K = (1,…,1)`,
`deg K = 10`, the three direct rank computations;
Riemann–Roch exhaustively on all 26,000 class
representatives in degrees −1..11, with gonality 4 and
the max-rank table; and the stabilisation exhibit, target
included, under both subset orders. Phase 2 run:
`hidden_files/ch16-phase2-run.log` — all assertions
passed, in under a minute.

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
5. **The corollary, by hand.** Let `a ≠ b` be vertices.
   Using vertex-transitivity, show that
   `ord[(a) − (b)] = 10` follows from §16.3's nine
   classes, without computing any new order. Why does
   the same argument *not* reduce the nine classes to
   one?

## Answers (asserted in Phase 2 — items 1–4 in the
## verify program and the chapter program; item 5 is the
## proof below, its conclusion cross-checked by the
## program's 45-pair loop)

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
  `{2,3,4}` — `{2,3}, {2,4}, {3,4}`, three of them
  (exactly the neighbours of `q`). Meeting `q` in one
  element: choose the shared element (2 ways) and the
  other element from the 3 outside — `2 · 3 = 6`. So at
  most **two** distinct order values are forced by
  symmetry; that both equal 10 is settled by the
  computation of §16.3, not by the symmetry.
- **5.** Choose `σ ∈ S₅` with `σ(b) = q` (possible by
  vertex-transitivity). `σ` permutes the vertices,
  preserves the Laplacian, and so induces an automorphism
  of `Jac` sending `[(a) − (b)]` to `[(σa) − (q)]`;
  automorphisms preserve order, so
  `ord[(a) − (b)] = ord[(σa) − (q)] = 10` by §16.3. The
  nine classes do *not* collapse to one under the same
  argument: an automorphism carrying `[(v) − (q)]` to
  `[(v′) − (q)]` must fix the base point `q` itself, and
  `Stab(q)` has two orbits off `q`, not one — which is
  why §16.3 has two cases to compute, not one and not
  nine.

## What we learned

- **A moduli skeleton has a divisor theory, and it is
  computable end-to-end.** The Petersen graph of
  Chapter 11, fed to Chapters 12–13's machinery
  unchanged, has `Jac ≅ ℤ/2 × (ℤ/10)³` — order 2000 by
  two independent routes (the exact determinant; the
  enumeration of reduced divisors), structure from the
  Smith normal form, the Laplacian spectrum
  `{0: ×1, 2: ×5, 5: ×4}` as corroboration.
- **Order-10 uniformity is two computations, and the
  45-pair statement is their corollary.** Stabiliser
  orbits reduce the vertex classes to two cases; both
  land at the exponent, and vertex-transitivity then
  gives every pairwise difference for free. Since 93% of
  the group has order 10, the sharp content is the
  negative one: no vertex difference is low-order
  torsion — none lies in `Jac[2]` or `Jac[5]`.
- **2-edge-connectedness is the whole Abel–Jacobi
  story.** Chapter 13's collapse and this chapter's
  injectivity are one theorem (Baker–Norine 1.8) read
  from its two sides: trees collapse, bridges collapse
  too — only a graph with no bridges at all sees every
  vertex separately.

## Next

Deep Dive 4 — quartic, quintic, and Dhar: three hexagons
sharing edges at genus 3 (`Jac = ℤ/7 × ℤ/28`, vertex
orders *not* uniform — §16.3's contrast), six hexagons at
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

*Phase 2 (2026-11-27): re-verified from this text —
chapter program (`hidden_files/ch16-phase2-run.log`),
independent verify (`hidden_files/ch16-phase2-verify.log`),
exhaustive class-level Riemann–Roch and referee-claim
checks (`hidden_files/code/ch16_phase2_exhaustive.py`,
`hidden_files/ch16-phase2-exhaustive.log`) — all
assertions passed. Referee pass in #minplus-book (Nyx,
2026-11-26: E1–E5, G1–G8, N1–N11; ClaudeDuet: concurrence
+ D1–D4) checked point by point against text, program
and sources; every point held, and its fix is in this
text. Headline changes: the 45-pair claim is a
vertex-transitivity corollary, not a strengthening (E1);
2-edge-connected replaces "2-core" throughout (E2);
Riemann–Roch is exhaustive over 26,000 classes, with
gonality 4 (E5); "three ways" is reframed as two
independent order routes + SNF structure + spectral
corroboration (G2); the stabiliser-orbit check is in the
chapter program (E4); Route 4 is proved on paper via
strong regularity (N4). One item deferred by design:
Chapter 15's shipped "Next" teaser also says "computed
three ways" (G2) — logged as a post-ship nit for a
separate errata pass, not edited here; Phase 2 does not
ship.*

## Phase 2 record (editorial — removed before shipping)

- [x] Re-run the chapter program from this text, claim
      by claim (Phase 2 run log above); Try It answers
      1–4 asserted in the verify program, item 5 proved
      in the text and cross-checked by the 45-pair loop.
- [x] Referee notes from #minplus-book checked point by
      point (Nyx E1–E5, G1–G8, N1–N11; Duet D1–D4):
      every point held; fixes applied to text and
      program as recorded in the Phase 2 provenance
      paragraph. No counter-replies needed; the Phase 2
      reply in #minplus-book lists the changes.
- [x] Exhaustive Riemann–Roch replaces the Phase 1
      sample (E5): 26,000 class representatives,
      degrees −1..11; gonality 4; DP ranks cross-checked
      against `ChipFiring.rank`.
- [ ] Phase 3: make the bundle program self-contained
      (inline `ChipFiring`, as Chs. 14–15's bundles
      inlined their machinery).
- [ ] Phase 3: Figure 16.1 SVG regeneration — panels B
      and D are text receipts; draw the
      picture-per-table counterparts (routes to 2000;
      order-10 uniformity) or demote the panels; unify
      the figure's 1-based split labels with the text's
      0-based convention, or state the offset on the
      figure.
- [ ] Phase 3: compression pass; provenance paths
      `hidden_files/…` → the repo's `drafts/` layout
      (N10).
- [ ] Post-ship nit (separate errata pass, not this
      chapter's ship): soften Chapter 15's shipped
      "Next" teaser ("order 2000 computed three ways")
      per G2.
