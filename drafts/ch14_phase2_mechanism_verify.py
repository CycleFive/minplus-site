#!/usr/bin/env python3
"""Ch.14 Phase 2: GLM mechanism, restated precisely and asserted.

GLM (via Lothrop, 2026-11-20) proposed: a defect at rank j cannot
affect m_k for k < j, cost starts at m_j and propagates upward.
As stated it was too clean: the m_k are Chow-ring intersection
numbers, not rank-k Mobius sums. Restated against the actual
definition (Baker Thm 4.2: m_k are the unsigned coefficients of
chi-bar = chi/(t-1)), the shape claim becomes a theorem of the
synthetic division itself:

  chi(t) = sum_F mu(empty,F) t^{r-rk(F)}; write
  w_i = sum_{rk(F)=i} mu(empty,F)   (Whitney numbers, signed).
  Dividing by (t-1) (see mobius_chibar): quot[k] = w_0+...+w_k,
  so   m_k = |w_0 + ... + w_k|   -- m_k depends ONLY on w_i, i<=k.

  mu(empty,F) depends only on the interval [empty,F], i.e. on
  flats of rank < rk(F). Hence if two matroids on the same ground
  set have identical flats of every rank < j (so identical w_i
  for i < j), then m_k agree for all k < j. A defect whose first
  non-uniform rank is j can first cost at m_j.

This script asserts, for every defect pair in the chapter plus
the circuit-hyperplane family:
  (1) m_k == |cumsum(w)_k| for every k;
  (2) flats of rank < j are identical to the uniform matroid's,
      flats of rank j differ;
  (3) w_i equal for i < j, w_j differs;
  (4) m_k equal for k < j, and the observed costs at k >= j.
What remains open (the narrowed §14.7 problem): the MAGNITUDES
above j -- why the rank-2 line defect's cost grows 1 -> 3 from
m_2 to m_3 (equivalently why mu(E) drops 10 -> 7) is propagation
through the Mobius recursion, not forced by (1)-(4) alone.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itertools import combinations
from ch14_chapter_program_draft import (uniform_flats, m_line_rank3,
                                        m_line_rank4, mobius_chibar)

E6 = list(range(6))
_DEFECTS = [frozenset({0,1,2,3}), frozenset({0,1,4,5}), frozenset({2,3,4,5})]

def paving(k):
    """Circuit-hyperplane paving matroids, copied from
    ch14_phase2_circuit_hyperplane_verify.py (that module runs its
    checks at import, so it is not imported here)."""
    D = _DEFECTS[:k]
    def rk(S):
        S = frozenset(S)
        if S in D: return 3
        return min(4, len(S))
    allsubs = [frozenset(c) for s in range(7) for c in combinations(E6, s)]
    def closure(S):
        S = set(S); r = rk(S)
        return frozenset(e for e in E6 if e in S or rk(set(S)|{e}) == r)
    flats = sorted({closure(S) for S in allsubs},
                   key=lambda f: (rk(f), len(f), sorted(f)))
    return flats, rk

def whitney(E, flats, rk):
    m, mu = mobius_chibar(E, flats, rk)
    r = max(rk(f) for f in flats)
    w = [0]*(r+1)
    for f in flats:
        w[rk(f)] += int(mu[f])
    cum, s = [], 0
    for x in w:
        s += x; cum.append(abs(s))
    assert cum[:-0 or None] and cum[:len(m)] == m, (w, cum, m)
    # cum has r entries after dropping the final total (which is 0)
    assert cum[len(m)] == 0
    return w, m

def flats_by_rank(flats, rk, j):
    return {i: {f for f in flats if rk(f) == i} for i in range(j+1)}

def check_pair(name, pair_u, pair_d, j, expect_cost):
    Eu, fu, ru = pair_u; Ed, fd, rd = pair_d
    assert Eu == Ed
    wu, mu_ = whitney(Eu, fu, ru)
    wd, md = whitney(Ed, fd, rd)
    bu, bd = flats_by_rank(fu, ru, j), flats_by_rank(fd, rd, j)
    for i in range(j):
        assert bu[i] == bd[i], (name, "flats differ below j", i)
        assert wu[i] == wd[i], (name, "w differs below j", i)
    assert bu[j] != bd[j], (name, "rank-j flats should differ")
    assert wu[j] != wd[j], (name, "w_j should differ")
    for k in range(j):
        assert mu_[k] == md[k], (name, "m_k differs below j", k)
    cost = [a-b for a, b in zip(mu_, md)]
    assert cost == expect_cost, (name, cost, expect_cost)
    print(f"{name}: w_u={wu} w_d={wd}; m_u={mu_} m_d={md}; "
          f"cost={cost}; flats & w agree below rank {j}, differ at {j}: OK")
    return wu, wd

# rank-3 pair: defect at j=2
check_pair("U(3,6) vs M_line (rank 3)",
           (list(range(6)), *uniform_flats(list(range(6)), 3)),
           m_line_rank3(), 2, [0, 0, 1])

# rank-4 line pair: defect at j=2; w_2 net -1 (line mu +2 vs 3 pairs +1)
_E4, _f4, _r4, _byrank = m_line_rank4()
wu, wd = check_pair("U(4,6) vs M_line^(4)",
                    (list(range(6)), *uniform_flats(list(range(6)), 4)),
                    (_E4, _f4, _r4), 2, [0, 0, 1, 3])
assert wd[2] - wu[2] == -1, "rank-2 Whitney net: 12*1+2 - 15*1 = -1"
print("  w_2 net -1 = merged line (+2) replacing three pair flats (+1 each): OK")

# circuit-hyperplane family: defects at j=3, costs (0,0,0,k)
for k in (1, 2, 3):
    fd, rd = paving(k)
    check_pair(f"U(4,6) vs circuit-hyperplane k={k}",
               (E6, *uniform_flats(E6, 4)),
               (E6, fd, rd), 3, [0, 0, 0, k])

print("\nGLM MECHANISM (restated): shape PROVED by the chi-bar")
print("synthetic division + lower-interval Mobius dependency,")
print("asserted on all three defect families. Magnitudes above")
print("rank j (the propagation step) remain the open problem.")
