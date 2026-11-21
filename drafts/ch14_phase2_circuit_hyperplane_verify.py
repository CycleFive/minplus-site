#!/usr/bin/env python3
"""Ch.14 Phase 2 (pre-run, proactivity 2026-11-20 evening):
verify the circuit-hyperplane brief family ("The Complete Defect
Ledger", 2026-10-29) in code, per the draft's Phase 2 checklist.

Construction: sparse paving rank-4 matroids on {0..5} whose only
non-bases are k circuit-hyperplanes (4-sets dropped to rank 3):
  D1={0,1,2,3}, D2={0,1,4,5}, D3={2,3,4,5} (pairwise triple-disjoint).
Rank: r(S)=min(4,|S|), except r(D)=3 for a defect D; supersets of a
defect still have rank 4 (they contain other 4-subset bases).
Flats = closures of all subsets under that rank function.

Brief claims to check (k=1,2,3):
  dims A* = (1, 36-3k, 36-3k, 1);  m = (1,5,10,10-k) three ways
  (FY / Mobius chi-bar / Bjorner);  maximal flags = 120-12k;
  deg(l^3) for c(S)=|S|(6-|S|) = 20880-768k;
  proper flats = 41-3k.
Also GLM rank-principle check: each defect IS a rank-3 flat, so a
rank-3 defect predicts cost concentrated at m_3 only: (0,0,0,k).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itertools import combinations
from fractions import Fraction
from ch14_chapter_program_draft import ChowRing, mobius_chibar, bjorner_count

E = list(range(6))
DEFECTS = [frozenset({0,1,2,3}), frozenset({0,1,4,5}), frozenset({2,3,4,5})]

def paving(k):
    D = DEFECTS[:k]
    def rk(S):
        S = frozenset(S)
        if S in D: return 3
        return min(4, len(S))
    allsubs = [frozenset(c) for s in range(7) for c in combinations(E, s)]
    # sanity: rank axioms spot-check submodularity on all pairs
    for A in allsubs:
        for B in allsubs:
            assert rk(A|B) + rk(A&B) <= rk(A) + rk(B), (A,B)
    def closure(S):
        S = set(S); r = rk(S)
        return frozenset(e for e in E if e in S or rk(set(S)|{e}) == r)
    flats = sorted({closure(S) for S in allsubs},
                   key=lambda f: (rk(f), len(f), sorted(f)))
    return flats, rk

def maximal_flags(flats, rk):
    Fp = [f for f in flats if f and f != frozenset(E)]
    cnt = 0
    for a in Fp:
        if rk(a) != 1: continue
        for b in Fp:
            if rk(b) != 2 or not a < b: continue
            for c in Fp:
                if rk(c) == 3 and b < c: cnt += 1
    return cnt

for k in (1, 2, 3):
    flats, rk = paving(k)
    proper = [f for f in flats if f and f != frozenset(E)]
    print(f"== k={k} defects ==")
    print(f"  flats total={len(flats)} proper={len(proper)} (brief proper 41-3k={41-3*k})")
    assert len(proper) == 41 - 3*k
    for D in DEFECTS[:k]:
        assert D in flats and rk(D) == 3, "defect must be a rank-3 flat"
    CR = ChowRing(E, flats, 4)
    dims = tuple(len(CR.basis[d]) for d in range(4))
    print(f"  dims={dims} expect {(1, 36-3*k, 36-3*k, 1)}")
    assert dims == (1, 36-3*k, 36-3*k, 1)
    # FY route (same alpha/beta construction as analyse())
    a_lift, b_lift = {}, {}
    for f, i in CR.gidx.items():
        (a_lift if E[0] in f else b_lift)[i] = Fraction(1)
    m_fy = []
    for kk in range(4):
        cur = {0: Fraction(1)}; dd = 0
        for lf in [a_lift]*(3-kk) + [b_lift]*kk:
            cur = CR.mul(dd, cur, lf); dd += 1
        m_fy.append(int(CR.deg(cur)))
    m_mb, mu = mobius_chibar(E, flats, rk)
    bj, _ = bjorner_count(E, flats, rk, 3)
    expect = [1, 5, 10, 10-k]
    print(f"  m FY={m_fy} Mobius={m_mb} Bjorner={bj} expect={expect}")
    assert m_fy == m_mb == bj == expect
    fl = maximal_flags(flats, rk)
    print(f"  maximal flags={fl} expect {120-12*k}")
    assert fl == 120 - 12*k
    # deg(l^3), ample c(S)=|S|(6-|S|)
    l_lift = {i: Fraction(len(f)*(6-len(f))) for f, i in CR.gidx.items()}
    cur = {0: Fraction(1)}; dd = 0
    for _ in range(3):
        cur = CR.mul(dd, cur, l_lift); dd += 1
    deg = int(CR.deg(cur))
    print(f"  deg(l^3)={deg} expect {20880-768*k}")
    assert deg == 20880 - 768*k
    print(f"  k={k}: ALL BRIEF CLAIMS VERIFIED; cost vs U(4,6) = (0,0,0,{k})")

print("\nCIRCUIT-HYPERPLANE FAMILY: VERIFIED for k=1,2,3")
