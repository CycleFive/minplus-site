# Learn You a Tropical Geometry for Great Good!

## Chapter 15 — The Matroid With No Matrices
### Deep Dive 2

> **DRAFT — Phase 2 verified, 2026-11-24.** Re-run
> end-to-end from this text: every claim asserted, the
> four Try It answers asserted individually, referee
> fixes 1–2 applied and notes 3–5 dispositioned (see the
> checklist at the end). Not yet compressed or shipped
> (Phase 3). Figures are the Day 16 PNG pending Phase 3
> SVG regeneration.

Chapter 14 ended on a sceptic's question half-answered. The
Kähler package survived one collinear triple — but that
matroid still had a *vector model*: you could write down
coordinates, and a cynic could mutter that somewhere
underneath, linear algebra was doing the work. This chapter
removes the floor entirely.

> **The Vámos matroid `V₈` has no vector model over any
> field — the program proves it, in code, by an Ingleton
> count of `16 ≰ 15` — and the full Kähler package holds
> for it anyway: Poincaré duality, a Hard Lefschetz
> isomorphism at dimension 70, Hodge–Riemann positive
> definite on a 69-dimensional primitive space, and
> intersection numbers `m = (1, 7, 21, 30)` recovered by
> three independent routes. Its five circuit-hyperplanes
> cost exactly `(0, 0, 0, 5)` against the uniform
> `U(4,8)`.**

There are no matrices in this chapter because there *can*
be none. What remains — a rank function, a lattice of
flats, a ring built from that lattice — is apparently
enough for Hodge theory. That is the Adiprasito–Huh–Katz
theorem with nowhere left to hide.

The Deep Dive door warning applies, doubled: Chapters 8, 9
and 14 are assumed (the Chow ring, the three routes to `m`,
the defect-ledger shape of Chapter 14 §14.7 in particular —
this chapter's headline cost is that section's sibling
family, now promoted to the main act).

Sources: M. Baker, "Hodge Theory in Combinatorics"
(arXiv:1705.07960), §§4.1, 4.3, 4.4, surveying
Adiprasito–Huh–Katz, Theorem 4.1 (the Kähler package) and
Theorem 4.2 / Proposition 4.3 (the three routes). The
Vámos definition: "The Vámos Matroid Has No Second
Symmetric Power" (Definition 3.1), agreeing with Oxley,
Example 2.1.25 — verified live on Day 16, not from memory.
Ingleton's inequality: Ingleton (1971).

## 15.1 Eight points, four pairs, five exceptions

Split the ground set `E = {0, …, 7}` into four pairs —
blocks, in the study's notation:

```
A = {0,1},   B = {2,3},   C = {4,5},   D = {6,7}.
```

Now declare a rank-4 matroid by the smallest possible act
of violence against uniformity: *every* 4-subset is a
basis, except five — the block unions

```
H = { A∪B, A∪C, A∪D, B∪C, B∪D },
```

which get rank 3 instead of 4. The rank function is one
line:

```
ρ(Y) = 3 if Y ∈ H, else min(|Y|, 4).
```

Note what is *not* in `H`: `C∪D = {4,5,6,7}` is a basis.
The block picture is a `K₄` on `{A, B, C, D}` with exactly
one edge missing — Figure 15.1A draws it, the five red
edges the circuit-hyperplanes, the dashed `C∪D` edge the
basis that breaks the symmetry just enough. (The matroid
is self-dual up to relabelling — dual equals `V₈` with the
blocks swapped `A↔C`, `B↔D`, checked on all 256 subsets —
but not identically self-dual. The missing edge matters.)

**Predict first.** Before any ring is built, the flat
lattice can be counted by hand, and Day 16 recorded the
prediction before running: flats by (rank, size)
`(0,0):1, (1,1):8, (2,2):28, (3,3):36, (3,4):5, (4,8):1` —
79 flats, 77 of them proper and nonempty; `μ` of a
circuit-hyperplane flat `−3`, `μ(E) = +30`; hence
`χ(t) = t⁴ − 8t³ + 28t² − 51t + 30`,
`χ̄(t) = t³ − 7t² + 21t − 30`, and `m = (1, 7, 21, 30)`.
Every one of those numbers is asserted in the program.

First, though, the definition has to earn the word
*matroid*:

```
submodularity ρ(X∪Y) + ρ(X∩Y) ≤ ρ(X) + ρ(Y)
  on all 256² subset pairs:        it is a matroid
bases: 65 of the 70 four-subsets; non-bases exactly H
circuits: the 5 H sets (size 4) + 36 of size 5 = 41,
  all of size ≥ rank — V₈ is a paving matroid
flats by (rank, size): as predicted, 79 in total
```

## 15.2 No matrices: Ingleton, 16 ≰ 15

A matroid is **representable** if some matrix, over some
field, has exactly its independent sets as its linearly
independent column sets. Every representable rank function
satisfies the **Ingleton inequality**: for any four
subsets `X₁, …, X₄`,

```
ρ(X₁)+ρ(X₂)+ρ(X₁X₂X₃)+ρ(X₁X₂X₄)+ρ(X₃X₄)
  ≤ ρ(X₁X₂)+ρ(X₁X₃)+ρ(X₁X₄)+ρ(X₂X₃)+ρ(X₂X₄).
```

Submodularity alone does not imply it — Ingleton is a
genuinely extra condition that linear algebra imposes and
combinatorics does not. So feed it the four blocks,
`X = (A, B, C, D)`, and count. Every term is a rank the
reader can verify from §15.1's one-line definition:

```
LHS = ρ(A)+ρ(B)+ρ(ABC)+ρ(ABD)+ρ(CD) = 2+2+4+4+4 = 16
RHS = ρ(AB)+ρ(AC)+ρ(AD)+ρ(BC)+ρ(BD) = 3+3+3+3+3 = 15
```

Look at where the asymmetry lives: the left side collects
`ρ(CD) = 4` — the one block union that *is* a basis —
while every union on the right is a circuit-hyperplane at
rank 3. The missing edge of Figure 15.1A is precisely the
unit by which the inequality fails:

```
16 ≤ 15 is FALSE — violated by exactly 1.
Therefore no matrix over any field has this rank function.
```

That is the whole proof, and it is finite: one inequality,
ten rank evaluations, no field ever named. Figure 15.1C
draws the two bars. From here on, anything the chapter
verifies for `V₈` is verified for an object that linear
algebra provably cannot reach.

> **Honest Box — why one inequality is a proof.**
> Ingleton is necessary, not sufficient: the non-Pappus
> matroid (rank 3, 9 points) satisfies Ingleton and still
> has no vector model — Pappus's theorem forbids one over
> every field (Oxley; the referee checked Ingleton on
> every 4-tuple of its flats, which suffices, since every
> term depends only on closures). Nelson's "Lots of
> Ingleton Matroids" (Matroid Union) gives more: the
> non-Desargues matroid, and the direct sum of the Fano
> and non-Fano matroids, also satisfy Ingleton for every
> quadruple and are non-representable — and by
> Mayhew–Newman–Whittle, no finite list of rank
> inequalities characterises representability at all.
> Passing Ingleton would prove nothing, which is part of
> why representability is subtle. *Failing* it is decisive
> in the other direction — a single violated necessary
> condition excludes every field at once. The chapter uses
> the inequality only in the direction it is valid, and
> the program asserts the direction explicitly: it checks
> submodularity *holds* everywhere (so `V₈` is a matroid)
> and Ingleton *fails* at `(A,B,C,D)` (so it is not a
> representable one). The gap between those two checks is
> exactly the gap this Deep Dive lives in.

## 15.3 Three routes to (1, 7, 21, 30)

Chapter 14's discipline, unchanged: no number is believed
until three routes that share no code path worth
mentioning agree on it.

```
Möbius / χ̄:      m = (1, 7, 21, 30)
                   μ(H flat) = −3,  μ(E) = +30
Björner flags:   m = (1, 7, 21, 30)
FY intersection: m = (1, 7, 21, 30)     (§15.4's ring)
log-concavity: 7² = 49 ≥ 1·21,   21² = 441 ≥ 7·30 = 210
```

And the uniform comparison on the *same* ground set, by
Möbius and Björner:

```
U(4,8):  m = (1, 7, 21, 35)
V₈:      m = (1, 7, 21, 30)
cost:          (0, 0, 0, 5)
```

Chapter 14 §14.7 proved the shape before this chapter
computed it: `m_k` is a cumulative Whitney sum, so a defect
first seen at rank `j` cannot cost anything before `m_j`.
A circuit-hyperplane is a rank-3 defect — five of them,
here — so the predicted cost shape was `(0, 0, 0, k)` with
`k = 5`. And the per-defect price is a proof, not an
induction from Chapter 14's `k = 1, 2, 3` family (the
referee's upgrade, verified in the Phase 2 program): in a
paving matroid of rank `r`, a circuit-hyperplane absorbs
its `r` subsets of size `r − 1` — each a flat with
`μ = (−1)^(r−1)` in the uniform matroid — into a single
flat whose only flats below it are the empty set, its
singletons and its pairs, so its Möbius value is
`(−1)^(r−1)(r − 1)`. The rank-`(r−1)` Whitney number
therefore moves by exactly `+1` per defect. Here `r = 4`:
each `H` absorbs 4 triple flats (`μ = −1` each) into one
flat with `μ = −3`, so `w₃` goes from `−56` to `−51`, and
`m₃ = |1 − 8 + 28 − 51| = 30` against the uniform
`|1 − 8 + 28 − 56| = 35`. Five defects, five units,
nothing anywhere else: the ledger balances before the
ring is even built. Figure 15.1B draws the two bar charts
with the `−5` arrowed.

## 15.4 The ring with nothing under it

The Feichtner–Yuzvinsky construction of Chapters 9 and 14
never asked for coordinates — it asks only for flats:

```
A*(V₈) = ℤ[x_F : F a nonempty proper flat] / (CH1 + CH2)
```

77 generators, one per proper nonempty flat; chain
monomials per degree `1, 77, 399, 997` (the sparse
machinery of Chapter 14 §14.3 — only chains survive (CH2),
which is what makes rank 4 computable at all). The ring
does not know, and cannot be told, that no matrix exists.
Its output:

```
dims A*(V₈) = (1, 70, 70, 1)        (Poincaré palindrome)
m by FY intersection = (1, 7, 21, 30)   — triple agreement
```

Dimension 70 in the middle degrees: this is not a toy
check. The palindrome alone is Poincaré duality in
embryo; §15.5 makes it, and the rest of the package,
explicit.

## 15.5 The full Kähler package, field-free

The ample class is the same recipe as Chapters 9 and 14 —
`ℓ` from `c(S) = |S|(8−|S|)` — with its strict
submodularity re-verified by brute force on all 52,670
incomparable ordered subset pairs (26,335 unordered) of
the 8-point ground set. Then,
every clause of the package, exactly:

```
Poincaré pairing A¹ × A² → ℤ:   nondegenerate (rank 70)
Hard Lefschetz ℓ : A¹ → A²:     an isomorphism (rank 70)
deg(ℓ³) = 150,720 > 0           (Hodge–Riemann, k = 0)
P¹ = ker(ℓ² : A¹ → A³),  dim P¹ = 69
Hodge–Riemann: −deg(·, ℓ·) positive definite on P¹ —
  exact LDLᵀ over ℚ: all 69 pivots > 0
  (min ≈ 5.66, max ≈ 41.73; pivots are ratios of
   successive Sylvester minors, so this IS Sylvester's
   criterion, not a float suggestion)
float spectrum of the P¹ form ⊂ [3.73, 103.18], all > 0
  (figure only — definiteness is the exact check above)
total runtime of the Day 16 program: ≈ 9 s
```

Figure 15.1D plots that spectrum: 69 eigenvalues, every
one strictly positive, none near zero. Compare Chapter
14's defect matroids — representable, every one — and the
checks are line-for-line identical. The theorem does not
distinguish them, because the theorem was never about the
models. It is about the lattice.

![The Kähler package for the Vámos matroid: the block K4 minus the C∪D edge, the m-bars against U(4,8), the Ingleton bars 16 vs 15, and the Hodge–Riemann spectrum on P1](figures/day16-vamos-kahler.png)

*Figure 15.1 — A: the block diagram of `V₈` — `K₄` on the
four pairs minus the `C∪D` edge; red unions are the five
circuit-hyperplanes (rank 3), the dashed blue union is a
basis (rank 4). B: `m` for `U(4,8)` vs `V₈`, the `−5`
arrowed at `k = 3` only. C: Ingleton — 16 against 15.
D: the Hodge–Riemann spectrum on `P¹` (dim 69), all
positive. (Day 16 PNG; Phase 3 regenerates as SVG with
every drawn value asserted by its generator.)*

## 15.6 What the absence of matrices actually means

It is worth being precise about what has and has not been
shown, because the headline invites over-reading.

*Not* shown: that coordinates are useless, or that `V₈`
has no geometry. Its Bergman fan exists — Chapter 8's
construction needs only flats, all 77 proper nonempty
ones give rays — and it is a tropical linear space that
tropicalises *no* linear space: there is no classical
plane for it to be the shadow of. Day 16 left that fan
undrawn (a 3-dimensional fan in a 7-dimensional quotient
is a slicing problem, not an existence problem; it stays
on the open list).

Shown, and it is the strongest form of the study's running
moral: **every hypothesis of the Adiprasito–Huh–Katz
package is combinatorial, so every conclusion is too.**
The Chow ring is built from the flat lattice. The ample
class is built from a submodular function on subsets. The
degree map is a normalisation on maximal flags. At no
point in Chapters 9, 14 or this one did a proof step —
as implemented, line by line — reach for a vector, a
field element, or a determinant of coordinates. `V₈` is
the control experiment that proves the implementation
wasn't cheating: feed the machine an input that *has* no
linear-algebraic shadow, and the machine's output is
unchanged in kind — duality, Lefschetz, Riemann, and a
log-concave `m` that three routes agree on.

June Huh's Fields-Medal work is often summarised as
"Hodge theory for matroids." This chapter is the
two-word version of why that phrase is not a metaphor:
the Hodge theory survives the *provable absence* of the
classical object it was modelled on.

## 15.7 The whole chapter, as a program

The draft program is `hidden_files/code/ch15_chapter_program_draft.py`
(Day 16's program, assembled for this chapter; it imports
the Day 14 sparse Chow ring by path, as Day 16 did — the
Phase 3 bundle program will be merged self-contained, as
Chapter 14's was). It asserts: the matroid axioms on all
256² pairs; bases, circuits and the flat census against
the hand prediction; Ingleton `16 ≰ 15`; self-duality up
to the block swap; Möbius = Björner = FY on
`(1, 7, 21, 30)`; the uniform `(1, 7, 21, 35)` and the
`(0, 0, 0, 5)` cost; strict submodularity on 52,670 pairs;
pairing and Lefschetz ranks 70; `deg(ℓ³) = 150,720`; and
all 69 exact LDLᵀ pivots positive. Run logs:
`hidden_files/ch15-draft-run.log` (Phase 1 assembly) and
`hidden_files/ch15-phase2-run.log` (Phase 2 re-run *from
this text*, claim by claim) — all assertions passed both
times. The Phase 2 companion program
`hidden_files/code/ch15_phase2_verify.py` asserts what
the chapter program states but does not itself check:
the four Try It answers individually, the full Whitney
numbers `(1, −8, 28, −51, 30)` and `(1, −8, 28, −56, 35)`,
log-concavity with the true products, the per-defect
Möbius accounting of §15.3, and the exact value
`deg(ℓ³) = 150,720` (via the run's summary record).

## 15.8 Try It

1. **One rank by hand.** From §15.1's definition, what is
   `ρ({0,1,4,5})`, and what is `ρ({4,5,6,7})`?
2. **Ingleton arithmetic.** Verify the RHS count of
   §15.2: why is every one of its five terms 3?
3. **Read the ledger.** `U(4,8)` has `m = (1, 7, 21, 35)`
   and `V₈` has `m = (1, 7, 21, 30)`. What do the five
   circuit-hyperplanes cost, entry by entry — and why
   does Chapter 14 §14.7's rank principle say the cost
   could not have appeared any earlier?
4. **A monomial that dies.** In `A*(V₈)`, why is
   `x_{A∪B} · x_{A∪C} = 0`?

## Answers (asserted in the Phase 2 program,
## `hidden_files/code/ch15_phase2_verify.py`)

- **1.** `{0,1,4,5} = A∪C ∈ H`, so `ρ = 3`.
  `{4,5,6,7} = C∪D ∉ H`, so `ρ = min(4, 4) = 4` — it is a
  basis. The two answers differ by exactly the missing
  edge of Figure 15.1A.
- **2.** Each RHS union is a union of two blocks that lies
  in `H` — `AB, AC, AD, BC, BD` are precisely the five
  circuit-hyperplanes — so each has rank 3 by definition,
  and `3+3+3+3+3 = 15`.
- **3.** `(0, 0, 0, 5)`. A circuit-hyperplane is a flat of
  rank 3, so the defect is first seen at rank `j = 3`;
  `m_k` is a cumulative Whitney sum over ranks `≤ k`, so
  no `m_k` with `k < 3` can feel it — and at `k = 3` the
  five defects cost one unit each.
- **4.** `A∪B` and `A∪C` are incomparable flats — neither
  contains the other (each has an element the other
  lacks) — so (CH2) kills the product. Only chains
  survive (§15.4).

## What we learned

- **Non-representability is checkable, and the check is
  small.** One inequality, ten rank evaluations:
  Ingleton `16 ≰ 15`, the single failing unit carried by
  the one block union that is a basis.
- **The Kähler package is field-free in fact, not just in
  statement.** Duality, a dimension-70 Lefschetz
  isomorphism, and Hodge–Riemann on a 69-dimensional
  primitive space — all exact — for a matroid no matrix
  over any field realises.
- **The defect ledger extends to the extreme case.**
  Five circuit-hyperplanes cost `(0, 0, 0, 5)` against
  `U(4,8)`, the rank-3 shape Chapter 14 §14.7 predicted —
  even when the "defect" is so large it destroys
  representability itself.

## Next

Deep Dive 3 — the Petersen graph has a Jacobian:
`Jac = ℤ/2 × (ℤ/10)³`, order 2000 computed three ways, on
the very graph Chapter 11 built as a moduli space.

---

*Source & provenance (draft): Day 16 (2026-10-26) of the
tropical geometry study — code
`hidden_files/code/day16_vamos.py` and
`day16_figures.py`, summary `day16_summary.json`, study
log entry for that day, run log `hidden_files/day16-run.log`.
Definition verified live on Day 16 against "The Vámos
Matroid Has No Second Symmetric Power" (Def. 3.1) /
Oxley Ex. 2.1.25; Ingleton statement fetched live the same
day (Ingleton 1971). Phase 1 assembly run 2026-11-23:
`hidden_files/ch15-draft-run.log`, all assertions passed.
Phase 2 re-verification 2026-11-24:
`hidden_files/ch15-phase2-run.log` +
`hidden_files/ch15-phase2-verify.log`, all assertions
passed; referee pass by ClaudeDuet (#minplus-book,
2026-11-23) — verdict holds, fixes 1–2 applied, notes
3–5 dispositioned in the checklist below. The §15.2
named examples: Oxley (non-Pappus) via the referee;
Nelson, "Lots of Ingleton Matroids", The Matroid Union
(non-Desargues; Fano ⊕ non-Fano; Mayhew–Newman–Whittle).*

## Phase 2 checklist (editorial — removed before shipping)

- [x] Re-run the chapter program from this text, claim by
      claim; assert the four Try It answers in code.
      DONE 2026-11-24: `ch15-phase2-run.log` (draft program,
      clean, 9.0 s) + `ch15_phase2_verify.py` /
      `ch15-phase2-verify.log` (Try Its 1–4, Whitney/χ/χ̄,
      exact deg(ℓ³) = 150,720 via the run summary). One
      catch: §15.3's log-concavity shorthand
      "21² = 441 ≥ 30" was not the inequality — fixed to
      "441 ≥ 7·30 = 210".
- [ ] Make the bundle program self-contained (inline the
      Day 14 sparse ring, as Ch. 14's bundle did).
      PHASE 3, as planned (referee note 3 dispositioned:
      the import is by path from the program's own
      directory, so the draft runs locally; the shipped
      bundle inlines it).
- [x] Check the §15.2 claim "plenty of non-representable
      matroids pass Ingleton" against a citable source or
      soften it. DONE: replaced with named examples —
      non-Pappus (referee's wording; Oxley) plus Nelson's
      non-Desargues and Fano ⊕ non-Fano, with the
      Mayhew–Newman–Whittle no-finite-list result.
- [x] Decide the Bergman-fan remark in §15.6 (77 rays /
      quotient dimension) wording against Ch. 8's
      conventions in Phase 2. DONE: counts verified
      (77 proper nonempty flats; quotient dim 7; fan
      dim 3); "nonempty" added; referee note 5 taken —
      the fan tropicalises no linear space, now stated.
- [x] Referee notes from #minplus-book: check each, fix
      what holds, reply to what doesn't. DONE 2026-11-24:
      Duet's pass (msg 1574629077325717586) — fix 1
      applied (author's call: both sourced example sets),
      fix 2 applied as a proved Möbius accounting and
      asserted in code, note 3 deferred to the Phase 3
      bundle as planned, note 4 applied (52,670 = ordered
      pairs; the hidden_files/ → drafts/ pointer swap is
      a Phase 3 bundle step), note 5 applied. Nothing in
      the pass failed to hold up; reply posted in
      #minplus-book.
- [ ] Phase 3: SVG figure regeneration with asserted
      values (incl. a drawn Ingleton/ledger visual if the
      PNG panels split), compression pass, ship.
