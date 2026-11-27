# Chapter 16 — Phase 2 consolidated checklist (prepared 2026-11-27, ahead of the 10:06 ET session)

Merges three inputs that arrived at different times:
- Phase 2 prep, 2026-11-26 (`ch16-phase2-prep.md`): draft program re-run passes
  (`ch16-phase2-draft-rerun.log`), independent verify passes
  (`ch16-phase2-verify.log`, no chip_firing.py), two cross-reference slips found.
- Nyx referee pass, 2026-11-26 (`ch16-referee-notes-2026-11-26.md`, mailbox seq 259–262):
  5 errors, 8 gaps, 11 nits on `ch16-draft` @ `9a95954`.
- Duet concurrence + D1–D4, 2026-11-26 (seq 263).

## Headline numbers — settled, do not re-hunt
Both referees independently re-confirmed, exact methods, no chip_firing.py:
det L_red = 2000; SNF [1^5, 2, 10^3] → Jac = Z/2 × (Z/10)^3; spectrum {0:1, 2:5, 5:4};
parking/Dhar count 2000 in distinct classes; order distribution 1 / 15 / 124 / 1860;
all 45 pair orders = 10. Phase 2 is a framing-and-fixes pass, not a numbers pass.

## Where the referee pass overrides the prep
1. **E1 — the 45-pair claim is NOT a strengthening.** Prep flagged "all 45 pairwise
   differences have order exactly 10" as the new claim needing full re-verification.
   It holds numerically, but Aut(P) = S5 is vertex-transitive, so it is equivalent to
   the 9-class claim: drop "Stronger, and new", state as a one-line corollary, keep
   the 45-loop as a code cross-check, add it as a Try It. Its re-verification
   priority can relax.
2. **E5 — replace the Riemann–Roch sample entirely.** The draft's seed-17 sample
   (and the prep re-run's 40 random divisors) only reach degrees 0–3, exactly where
   nothing special lives. Rank is a class function: enumerate 2000 reduced reps ×
   degrees −1..11 = 26k classes. Nyx ran this exhaustively: all OK, <2 s. New facts
   it surfaces: first r≥1 at degree 4 (divisorial gonality 4); max r by degree
   4..8 = 1,1,2,2,3. Also fix the "3^6" conflation (Ch 12's 6 vertices vs genus 6).
   Wording becomes "exhaustive", not "sample".
3. **E4 — port the stabiliser-orbit check into the chapter program.** The draft
   program (L134–141) only counts splits by |v∩q| (3, 6); the brute-force orbit
   proof (12-element stabiliser) exists only in `ch16_phase2_verify.py`. Port it,
   or reword §16.3/§16.7 to "partition sizes asserted; orbits in verify".
   (Duet D3: the program print also mislabels them "S_5-orbits" — they are Stab(q) orbits.)

## Errors to fix (text)
- E2: "2-core" → "2-edge-connected" in §16.1, §16.4, What we learned. Counterexample:
  two cycles joined by one edge is its own 2-core, has no pendant trees, and AJ still
  identifies the bridge's endpoints. BN 1.8's criterion is 2-edge-connected, full stop.
- E3: "Chapter 18's quartic skeleton (Deep Dive 4)" → Chapter 17. Third cross-reference
  slip, in addition to the prep's two: header line 11 (§16.4 → §16.3 for the pairwise
  claim) and §16.1 line 98 (§16.5 → §16.4 for Abel–Jacobi injectivity).

## Gaps to close
- G1: state Jac ≅ coker(L_red) (drop q; gauge f(q)=0) — 3-line lemma or exact back-ref.
- G2: "three ways" overstates independence — ∏SNF = |det| is automatic and spectral
  Kirchhoff is Matrix-Tree in eigen-form. Reframe: order two independent ways
  (det; parking), structure from SNF, spectrum as corroboration. Soften Ch 15's teaser too.
- G3: cite BN Prop 3.1 at Route 3 (one q-reduced rep per class).
- G4: three one-sentence symmetry steps in §16.3 (Stab(q) action preserves order;
  ord(x)=ord(−x) for unordered pairs; pair-orbits need distance-transitivity).
  The 3 disjoint splits are exactly N(q).
- G5: cut "sharpest fact" / "rigid" / "as far as the group allows". Quote the order
  distribution (93% of Jac has order 10; prior ≈0.86 for the two symmetry cases).
  Honest content: no vertex difference lies in Jac[2] ≅ (Z/2)^4 or Jac[5] ≅ (Z/5)^3.
- G6: "hand checks" → direct computations; give the actual short proofs for
  r(0)=0 and r((q))=0; flag r(K)=g−1 as an RR corollary computed from the definition.
- G7: make the torus "shadow" precise (unit edge lengths; Jac(G) ≅ Λ^#/Λ,
  Bacher–de la Harpe–Nagnibeda 1997 — verify cite); "skeleton of M_{0,5}" = link of
  the tropical fan / boundary complex of M̄_{0,5}, say once.
- G8: "10 fired sets" is algorithm-dependent (lex subset order; reversed order: 3).
  Qualify as "10 set-firings in the program's phase 2 (lex order)" or drop the number.

## Duet's additions
- D1: the three-hexagon control (Jac = Z/7 × Z/28, orders {7,14,28}) is a *contrast*,
  not a control — it isn't vertex-transitive, so it can't show symmetry couldn't
  force Petersen's uniformity.
- D2: §16.4's group-law identity holds for any well-defined reduction and is tested
  only at (a,b)=(1,2) — say so, or test associativity on random triples.
- D4: exhibit target verified q-reduced and principal — only the fired-set count is soft (G8).

## Nits (batch, N1–N11 in the referee file)
Stabiliser/stabilisation rename (N1); "ten nontrivial splits of 5 leaves" (N2);
Honest Box full-L SNF [1^5,2,10^3,0] made concrete (N3); Route 4 on paper via
srg(10,3,0,1) ⇒ A²+A−2I=J (N4); Route 1 float det — use sympy exact or say so (N5);
notation unification (N6); §16.4 injectivity already implied by §16.3 (N7);
runtime "under a minute" (N8); header checklist count (N9); provenance paths at
ship (N10, Phase 3); verify.py dead line (N11).

## Not Phase 2
Figure 16.1 panels B/D (text receipts) — Phase 3 SVG pass, per the draft's own checklist.
