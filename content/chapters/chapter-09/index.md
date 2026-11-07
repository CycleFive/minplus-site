+++
title = "Hodge Theory With No Hodge Theory Underneath"
description = "The Chow ring of U(3,4), built from its flats alone: intersection numbers (1, 3, 3) recovered three independent ways, and the Hodge–Riemann relations verified exactly — fractions, not floats — on the primitive hyperplane."
weight = 9
[extra]
number = 9
part = "Part IV — Matroids are linear spaces"
lead = "Chapter 8's fan has an algebra built from flats alone, and the Kähler package holds on it anyway, exactly, in code — proved from the flats for every matroid, including ones no variety realizes."
+++


Chapter 8 ended with a fan: the Bergman fan `B(U(3,4))`, the
tropical plane. Fans have algebra. The **Chow ring** `A*(M)` of a
matroid is built from its flats alone — no fan geometry, no
vectors, no field — and on it, Adiprasito–Huh–Katz proved the
Kähler package (Poincaré duality, Hard Lefschetz, Hodge–Riemann)
that classical geometry gets from Hodge theory on a smooth
projective variety. No Hodge theory is needed underneath: the
proof never leaves the lattice of flats, and it covers matroids
that no variety realizes at all (the Fano plane, §9.3). One
honesty note on this chapter's example: `U(3,4)` itself *is*
realizable — four lines in general position in the plane — and
its Chow ring is the cohomology ring of the plane blown up at
the six points where those lines cross, so there the relations
*are* the classical Hodge index theorem in disguise. The
computation below never touches that variety, and the same
program, fed different flats, checks the package where there is
provably nothing to blow up. This chapter verifies the relations
*exactly* — fractions, not floats — on the smallest interesting
case, `M = U(3,4)`. This is the mechanism behind June Huh's
Fields Medal work, in miniature, in code.

Prediction, in the book's usual order: if the algebra is really
there, then the same numbers should be reachable by routes that
share nothing — pure ring theory, pure Möbius combinatorics, and
pure flag-counting should collide on one answer. They do, in §9.2.

## 9.1 The ring, from the flats

Presentation (Feichtner–Yuzvinsky, arXiv:math/0305142; via
Baker's survey arXiv:1705.07960 §4.1): one variable `x_F` per
nonempty proper flat, modulo

- **(CH1)** for all `a, b ∈ E`: `Σ_{F∋a} x_F = Σ_{F∋b} x_F`
  (linear relations);
- **(CH2)** `x_F · x_G = 0` whenever `F, G` are incomparable
  (monomial relations — only chains of flats survive).

For `U(3,4)`: 10 generators, 3 independent linear relations,
graded dimensions **(1, 7, 1)** — degree 2 is one-dimensional, as
Poincaré duality wants. The degree map `deg: A² → ℤ` is pinned by
`deg = 1` on every maximal-chain monomial, and the computation
finds it *constant* on all 12 maximal chains (one nonzero
constant — the program's raw functional gives −1; the sign is
the representative's, the constancy is the content — normalized
to 1) — the tropical reading: `deg` is the
unique integer weighting making the Bergman fan balanced
(Chapter 8's balancing, returning as algebra).

One subdivision note, so this chapter and Chapter 8 agree: the
presentation above is built from the lattice of flats alone — no
fan, and no subdivision of one, is chosen to define it. When the
ring *is* read as the Chow ring of the Bergman fan, it is the
fan's **fine subdivision** (Chapter 8, description 2) whose rays
are the variables `x_F` — one per nonempty proper flat, 10 here.
The coarsest fan structure on the same set has only 4 rays and 6
cones and carries no such variables. Definition: flats. Fan
reading: fine.

## 9.2 The numbers (1, 3, 3), three ways

With `α = Σ_{F∋0} x_F` and `β = Σ_{F∌0} x_F` (Huh–Katz's two
classes; Baker §4.4), set `m_k = deg(α^{2−k} β^k)`. The chapter
computes:

- `deg(α²) = 1`, `deg(α·β) = 3`, `deg(β²) = 3` — by intersection
  theory in the ring;
- the same (1, 3, 3) as the unsigned coefficients of the reduced
  characteristic polynomial `χ̄(t) = t² − 3t + 3`, itself
  Möbius-verified from `χ(t) = t³ − 4t² + 6t − 3`. (Identifying
  the coefficients as intersection numbers `deg(α^{2−k} β^k)` is
  Huh–Katz, arXiv:1104.2519, proved there for realizable
  matroids; Baker Thm 4.2 is AHK's all-matroids form, on this
  ring. AHK's own step is the Hodge–Riemann relations of §9.3,
  for every matroid at once.);
- the same (1, 3, 3) a third way, by Björner's count of initial
  descending flags (Baker Prop 4.3) — pure combinatorics, no ring
  at all.

Three routes, one answer. And the point of the numbers:
`3² ≥ 1·3` — log-concavity, an instance of the Rota–Welsh
conjecture, which is what the Hodge–Riemann relations below are
*for*.

| quantity | route 1: ring | route 2: Möbius | route 3: Björner |
|---|---|---|---|
| `m_0` | 1 | 1 | 1 |
| `m_1` | 3 | 3 | 3 |
| `m_2` | 3 | 3 | 3 |

![The numbers (1, 3, 3) arriving three ways: by intersection theory in the ring, as the unsigned coefficients of the reduced characteristic polynomial, and by Björner's descending-flag count](figures/ch09-three-routes.svg)

*Figure 9.2 — The table above, drawn. Three computations that
share no machinery land on the same three bars; every bar was
read off the chapter program's own summary, after the program
asserted the three routes pairwise equal.*

## 9.3 Hodge–Riemann, checked

An ample class comes from a strictly submodular function:
`l = Σ_F |F|(4−|F|) x_F`, with `c(S) = |S|(n−|S|)` verified
strictly submodular on all 110 ordered incomparable subset pairs
(55 unordered) by brute force. Then:

- `deg(l²) = 120 > 0` (the `k = 0` relation);
- on the 6-dimensional primitive hyperplane
  `P¹ = ker(l : A¹ → A²)`, the form `−deg(·,·)` is positive
  definite — Sylvester leading minors `97/18, 260/9, 925/6, 435,
  2133/2, 2430`, all positive, exact fractions, **in the
  program's basis of `P¹`**. Eigenvalues (float, for the figure,
  same basis): `1, 1, 10/3, 9, 9, 9`.

The minors and eigenvalues are coordinates, not invariants: they
belong to the basis the form's matrix is written in, and they
move when the basis does. The program checks one determinant-1
change of basis itself: the minors become `32/9, 116/3, 470/3,
237, 540, 2430` and the eigenvalues ≈ `0.55, 1.59, 3.49, 6.34,
7.14, 17.55` — only the last minor survives, because it is the
determinant and the change has determinant 1. What no basis can
move is the *signature*: by Sylvester's law of inertia the
pairing on `A¹` has signature `(1, 6)` in every basis — one
positive direction, `l`'s, and six negative — and
positive-definiteness on `P¹` is exactly that basis-free fact.
Sylvester's criterion certifies it from any basis; the numbers
above are the certificate our basis happened to write.

> **Honest Box.** The first ample-class candidate,
> `c(F) = rk(F)(3−rk(F))`, is NOT strictly submodular:
> `A = {0,1}`, `B = {0,2,3}` gives equality — both sides are
> exactly 2, re-asserted in the chapter program — caught by the
> brute-force check before it could poison the ample class. The
> fix, `c(S) = |S|(n−|S|)`, works because `|·|` is modular and
> `t(n−t)` is strictly concave — and was verified exactly anyway.
> Trust the check, not the lemma-shaped intuition.

![The Chow ring of U(3,4): the 10 flats coloured alpha/beta, the intersection numbers (1,3,3) with the log-concavity check, and the Hodge–Riemann eigenvalues on P¹ in the program's basis](figures/day07-chow-ring-huh.svg)

*Figure 9.1 — Three panels: the 10 flats coloured α/β with a
highlighted maximal chain; the intersection numbers (1, 3, 3)
with the log-concavity check; the Hodge–Riemann eigenvalues on
P¹ in the program's basis — all positive, as the exact Sylvester
check in the text requires (the signature is the invariant; the
spectrum is the basis's). Every value drawn was asserted in the
chapter program before rendering.*

**No variety, same package.** The program runs the identical
machinery on the Fano plane `F7` — 7 points, 7 lines of 3 points
each — a matroid realizable over no field of characteristic
≠ 2, hence not over ℂ: nothing to blow up, no cohomology to
borrow. Out come 14 generators, graded dimensions `(1, 8, 1)`,
`deg` constant on all 21 maximal chains, `m = (1, 6, 8)` — the
unsigned coefficients of `χ̄(t) = t² − 6t + 8` — and signature
`(1, 7)` on `A¹`: the Kähler package from a lattice of flats,
where there is provably no Hodge theory underneath because
there is no variety. That, not `U(3,4)`, is the theorem's honest
advertisement — and where this chapter's title earns its keep.

## 9.4 The whole chapter, as a program

Everything in this chapter is one runnable program: [ch09_chow_ring.py](./ch09_chow_ring.py) — self-contained apart from numpy (eigenvalues, display only; every asserted quantity is an exact Fraction). Its assertions *are* the chapter's claims: graded dimensions (1, 7, 1) from 10 generators and 3 independent (CH1) relations; `deg` constant (one nonzero constant before normalization, = 1 after) on all 12 maximal-chain monomials; the Poincaré pairing on `A¹` nondegenerate (Gram rank 7); `m = (1, 3, 3)` by all three routes of §9.2, asserted pairwise equal, with `χ(t)` and `χ̄(t)` Möbius-verified; strict submodularity of `c(S) = |S|(4−|S|)` on all 110 ordered incomparable pairs, and the failed rank-based candidate's equality 2 = 2 at `A = {0,1}`, `B = {0,2,3}` reproduced as the Honest Box's evidence; Hodge–Riemann with `deg(l²) = 120` exactly and `−deg(·,·)` positive definite on the 6-dimensional `P¹`, Sylvester minors asserted as the exact list above — in the program's basis, with a determinant-1 change of that basis asserted to move the minors and eigenvalues while the signature (1, 6) stays; the Fano plane `F7` run through the same machinery (14 generators, dimensions (1, 8, 1), `m = (1, 6, 8)`, signature (1, 7)); and, asserted as the Try It answer key, the comparable-monomial count (22 survive (CH2), 33 are killed, all 22 nonzero under `deg`) and the single log-concavity inequality for (1, 4, 6). Run it, and if it prints its final line, the claims are what we said they are. (This section used to print the listing inline; a screenful of code you can download is just a download with extra scrolling.)

## 9.5 Try It

1. `A*(U(3,4))` has 10 generators and graded dimensions (1, 7, 1).
   How many independent linear relations are there in degree 1?
2. Why does `x_F · x_G` vanish when `F` and `G` are incomparable —
   what does a surviving degree-2 monomial look like?
3. The log-concavity instance here is `3² ≥ 1·3`. What is the
   corresponding inequality for the numbers (1, 4, 6) of
   `U(3,5)` (a later computation), and does it hold? (A 3-term
   sequence has exactly one inequality of this shape.)

## Answers

- **1.** **3** — 10 generators minus the 7-dimensional `A¹`.
- **2.** By decree — (CH2) is a defining relation, not a
  theorem — but a forced decree: a degree-2 monomial records a
  chain of flats, and no chain contains two incomparable flats
  in either order, so the product has nothing to record. A
  surviving degree-2 monomial is a product `x_F · x_G` with
  `F ⊆ G` or `G ⊆ F` — the maximal ones are the 12 chain monomials
  `x_{ {i} } · x_{ {i,j} }`, and each `x_F²` survives too (a flat
  is comparable with itself). The program counts 22 comparable
  monomials surviving and 33 incomparable ones killed — and
  `deg` is nonzero on all 22, so here "survives (CH2)" really
  does mean nonzero in `A²`.
- **3.** `4² ≥ 1·6`, i.e. `16 ≥ 6` — and yes, it holds.

## What we learned

- **The Kähler package needs no variety in the proof.**
  Poincaré duality, Hard Lefschetz's shadow and the
  Hodge–Riemann relations hold on a ring defined from a finite
  lattice of flats — verified here exactly, Sylvester minor by
  Sylvester minor, not quoted — on `F7`, where no variety
  exists, just as on `U(3,4)`, where one secretly does.
- **Three routes, one answer is evidence, not redundancy.** The
  ring, the characteristic polynomial and the flag count share
  no code path; their agreement on (1, 3, 3) is the chapter's
  real theorem-check, and it is Huh–Katz's identification in
  AHK's all-matroids form (Baker Thm 4.2), seen at the smallest
  size where it says anything.
- **Log-concavity is the corollary that pays.** `3² ≥ 1·3` is
  the Rota–Welsh conjecture for `U(3,4)` — the statement Huh's
  programme proves in general by exactly this mechanism.

## Next

Part IV closes here. Part V asks what a *moduli space* looks like
tropically — and the answer starts with a surprise from
Chapter 1: the tropical Grassmannian `Gr(2,4)`, the space of all
2-planes in 4-space, *is* the space of phylogenetic trees, and
itself a tropical line. Chapter 10 draws it.

---

*Source & provenance: Day 7 (2026-10-17) of the tropical geometry
study — code `hidden_files/code/day07_chow_ring.py` and
`hidden_files/code/day07b_hodge_riemann.py` (eigenvalues in
`hidden_files/day07b-eigenvalues.txt`), figure
`files/figures/day07-chow-ring-huh.png`, study log Day 7 entry.
Sources: M. Baker, "Hodge Theory in Combinatorics"
(arXiv:1705.07960), §§4.1, 4.3, 4.4 (FY presentation, deg, α/β,
Thm 4.2, Thm 4.1, Prop 4.3), re-verified against the ar5iv HTML
during Day 7 itself. The presentation is Feichtner–Yuzvinsky,
arXiv:math/0305142; the α/β identification is Huh–Katz,
arXiv:1104.2519; the Hodge–Riemann step for all matroids is
Adiprasito–Huh–Katz, arXiv:1511.02888. The chapter program was
merged into one self-contained file and re-run end-to-end on
2026-11-06, with every claim asserted —
including `deg(l²) = 120` exactly, the Sylvester-minor list, and
the failed candidate's equality — and re-run again by the figure
script on 2026-11-07 before rendering. Figures regenerated as SVG
the same day: code `hidden_files/code/ch09_extra_figures.py`,
which re-runs the chapter program, asserts every value it draws,
and only then renders. Referee patch 2026-11-07 (ClaudeDuet,
#minplus-book), each point verified against the program and the
sources before landing: the variety framing corrected (`U(3,4)`
is realizable; `F7` added as the variety-free case,
program-asserted), the §9.3 minors and eigenvalues labelled as
basis coordinates with the signature as the invariant, the α/β
identification credited to Huh–Katz, and the degree map's raw
constant asserted as constancy rather than as the
representative's −1.*

## Chapter code

The full program behind this chapter, as a download: [ch09_chow_ring.py](./ch09_chow_ring.py) — the same file §9.4 points at. Its assertions are the chapter's claims; run it.
