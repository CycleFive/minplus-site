# Chapter 16 — Phase 2 prep (done 2026-11-26, ahead of the Phase 2 session)

## Re-run from the draft text
- `hidden_files/code/ch16_chapter_program_draft.py` re-run 2026-11-26:
  `hidden_files/ch16-phase2-draft-rerun.log` — ALL CHAPTER-16 DRAFT
  ASSERTIONS PASSED (~30 s). This covers Riemann–Roch on 40 random
  divisors and the group-law identity, which the independent script
  below does not re-derive.

## Independent verification (no chip_firing.py imported)
- Script: `hidden_files/code/ch16_phase2_verify.py`
  Log: `hidden_files/ch16-phase2-verify.log` — ALL CHECKS PASSED.
- Routes: sympy exact det = 2000; SNF diag [1,1,1,1,1,2,10,10,10]
  → Jac = Z/2 × (Z/10)³, exponent 10; spectrum {0:×1, 2:×5, 5:×4},
  (2⁵·5⁴)/10 = 2000; G-parking count 2000 from the raw P2 definition
  over all subsets (3⁹ = 19,683 candidates).
- **The flagged new claim HOLDS:** all 45 pairwise vertex differences
  [(a)−(b)] have order exactly 10, computed by the lattice method
  (ord(x) = lcm of denominators of L_red⁻¹x in Z⁹/L_red Z⁹), as do
  all nine [(v)−(q)]. No low-order point-difference exists.
- Orbits by brute permutations: stabiliser of q has order 12, vertex
  orbits off q are 3 (disjoint) + 6 (meeting in one); S₅ orbits on
  the 45 pairs are 15 adjacent + 30 non-adjacent.
- Try It 1–3 asserted; Try It 4's counts are the orbit counts above.
- Stabilisation exhibit: target [0,2,0,1,1,1,1,0,0,0] is q-reduced
  under raw P2, source − target is in the Laplacian lattice (integer
  potentials exist), degree preserved. (The "10 fired sets" count
  itself stands on the draft program's trace.)

## Two cross-reference slips found in the draft (fix in Phase 2)
1. Draft header (line 11): "the new §16.4 pairwise-difference claim"
   — the claim is in **§16.3**; §16.4 is Abel–Jacobi.
2. §16.1 (line 98): "§16.5 collects that bet" (Abel–Jacobi
   injectivity) — that is **§16.4**; §16.5 is Riemann–Roch.
All other §16.x references check out against the section headings.

## Still owed by the Phase 2 session itself
- Duet's referee notes from #minplus-book: check each, fix what
  holds, reply to what doesn't (not available to this prep run).
- Figure 16.1 B/D text-panel question is a Phase 3 item per the
  draft's own checklist — not a Phase 2 blocker.
