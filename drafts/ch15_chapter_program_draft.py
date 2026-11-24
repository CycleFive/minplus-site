#!/usr/bin/env python3
"""Chapter 15 (Deep Dive 2) -- The Matroid With No Matrices. DRAFT program, Phase 1.
Assembled 2026-11-23 from Day 16 (2026-10-26): the Vamos matroid V8 -- the Kahler package where no vector
model exists (AHK is field-free).

Day 15's open problem #3, second option: a NON-REPRESENTABLE matroid
through the same Feichtner-Yuzvinsky machinery.  If A*(M), Poincare
duality, Hard Lefschetz and Hodge-Riemann are really combinatorial
(Baker, arXiv:1705.07960, sec 4.1/4.3/4.4; AHK Thm 4.1/4.2), they must
hold for V8 exactly as for the representable matroids of days 7/12/14/15
-- there is no field, no vectors, and no Bergman fan of a linear space
anywhere in sight, only a rank function.

DEFINITION (verified live this session, not from memory).
"The Vamos Matroid Has No Second Symmetric Power" (arXiv HTML,
Def. 3.1, presentation agreeing with Oxley, Example 2.1.25):
E = {1..8}, blocks A={1,2}, B={3,4}, C={5,6}, D={7,8},
H = {A+B, A+C, A+D, B+C, B+D} (block unions); V8 is the rank-4 matroid
whose bases are all 4-subsets EXCEPT the five in H; rank function
rho(Y) = 3 if Y in H, else min(|Y|, 4).
Relabeled 0-based here: A={0,1}, B={2,3}, C={4,5}, D={6,7};
H = {0123, 0145, 0167, 2345, 2367}; the missing block union C+D=4567
IS a basis.  Hand predictions made before running (Mobius over the
flat lattice, derived in the session notes): flats by rank/size
(0,0):1, (1,1):8, (2,2):28, (3,3):36, (3,4):5, (4,8):1;
mu = +1,-1,+1,-1 (rank<=2 and triple flats), mu(H flat) = -3,
mu(E) = +30; chi = t^4-8t^3+28t^2-51t+30;
chi-bar = t^3-7t^2+21t-30; m = (1,7,21,30) -- vs uniform U(4,8)
m = (1,7,21,35): the five circuit-hyperplanes should cost exactly 5,
at the top coefficient only.

NON-REPRESENTABILITY, proved in code: the rank function satisfies
submodularity on all 256^2 subset pairs (it IS a matroid) but violates
the INGLETON inequality (Ingleton 1971; statement as at
en.wikipedia.org/wiki/Ingleton's_inequality, fetched live):
for X1..X4, representable rho must satisfy
  rho(X1)+rho(X2)+rho(X1+X2+X3)+rho(X1+X2+X4)+rho(X3+X4)
  <= rho(X1+X2)+rho(X1+X3)+rho(X1+X4)+rho(X2+X3)+rho(X2+X4).
With X=(A,B,C,D): LHS = 2+2+4+4+4 = 16, RHS = 3+3+3+3+3 = 15.
16 <= 15 is false, so no matrix over ANY field has this rank function.

All arithmetic exact (Fractions; LDL^T over Q for HRR in place of
69 Bareiss minors -- for a symmetric matrix, completing LDL^T with no
pivoting gives pivots = ratios of successive leading principal minors,
so all pivots > 0 iff Sylvester's criterion holds).  Eigenvalues float,
for the figure only.  Writes day16_summary.json next to this file. NOTE (Phase 1): imports the Day 14 sparse ChowRing via importlib, as Day 16 did; the Phase 3 bundle program will be merged self-contained, as Ch. 14 was.
"""
import importlib.util
import itertools
import json
import os
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "day14", os.path.join(HERE, "day14_rank4_chow.py"))
day14 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(day14)
ChowRing, mobius_chibar, dense_rank = (
    day14.ChowRing, day14.mobius_chibar, day14.dense_rank)

SUMMARY = {}
t_start = time.time()

# ---------- the Vamos matroid, from the cited definition ----------
E = list(range(8))
BLOCKS = [frozenset({0, 1}), frozenset({2, 3}),
          frozenset({4, 5}), frozenset({6, 7})]
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]  # all but (2,3)=C,D
H = [BLOCKS[i] | BLOCKS[j] for i, j in PAIRS]
Hset = set(H)

def rho(S):
    S = frozenset(S)
    return 3 if S in Hset else min(len(S), 4)

allsubs = [frozenset(c) for s in range(9)
           for c in itertools.combinations(E, s)]
rk_of = {S: rho(S) for S in allsubs}

print("== V8: matroid axioms from the rank function ==")
for A in allsubs:
    for B in allsubs:
        assert rk_of[A | B] + rk_of[A & B] <= rk_of[A] + rk_of[B], (A, B)
print("  submodularity on all 256^2 subset pairs: it is a matroid")
bases = [S for S in allsubs if len(S) == 4 and rk_of[S] == 4]
nonbases = [S for S in allsubs if len(S) == 4 and rk_of[S] < 4]
assert len(bases) == 65 and sorted(nonbases) == sorted(H)
print(f"  bases: 65 of 70 four-subsets; nonbases exactly the 5 in H")

# circuits: minimal dependent sets
circuits = []
for S in allsubs:
    if S and rk_of[S] < len(S) and all(
            rk_of[S - {e}] == len(S) - 1 for e in S):
        circuits.append(S)
c4 = [c for c in circuits if len(c) == 4]
c5 = [c for c in circuits if len(c) == 5]
assert sorted(c4) == sorted(H) and len(c5) == 36 and len(circuits) == 41
assert all(len(c) >= 4 for c in circuits)
print(f"  circuits: 5 of size 4 (the H sets) + 36 of size 5 = 41; "
      f"all sizes >= rank 4: V8 is a paving matroid")

# flats by closure
def closure(S):
    S = set(S); r = rk_of[frozenset(S)]
    return frozenset(e for e in E
                     if e in S or rk_of[frozenset(S | {e})] == r)
flats = sorted({closure(S) for S in allsubs},
               key=lambda f: (rk_of[f], len(f), sorted(f)))
from collections import Counter
by_rs = Counter((rk_of[f], len(f)) for f in flats)
assert by_rs == Counter({(0, 0): 1, (1, 1): 8, (2, 2): 28,
                          (3, 3): 36, (3, 4): 5, (4, 8): 1}), by_rs
print(f"  flats by (rank,size): {dict(sorted(by_rs.items()))} "
      f"(79 flats, 77 proper nonempty) -- as predicted by hand")

# ---------- non-representability: Ingleton violation ----------
print("== Ingleton: V8 is not representable over any field ==")
A, B, C, D = BLOCKS
lhs = (rho(A) + rho(B) + rho(A | B | C) + rho(A | B | D) + rho(C | D))
rhs = (rho(A | B) + rho(A | C) + rho(A | D) + rho(B | C) + rho(B | D))
assert (lhs, rhs) == (16, 15) and lhs > rhs
print(f"  X=(A,B,C,D): LHS = 2+2+4+4+4 = {lhs}, "
      f"RHS = 3+3+3+3+3 = {rhs}; {lhs} <= {rhs} FALSE -- "
      f"Ingleton violated by exactly 1, so no vector model exists")
SUMMARY["ingleton"] = {"lhs": lhs, "rhs": rhs}

# ---------- self-duality up to isomorphism (bonus check) ----------
def rho_dual(S):
    S = frozenset(S)
    return len(S) + rho(frozenset(E) - S) - 4
perm = {0: 4, 1: 5, 4: 0, 5: 1, 2: 6, 3: 7, 6: 2, 7: 3}  # A<->C, B<->D
for S in allsubs:
    assert rho_dual(S) == rho(frozenset(perm[e] for e in S)), S
print("  dual matroid = V8 relabeled by block swap (A<->C, B<->D): "
      "self-dual up to isomorphism, not identically")

# ---------- Bjorner's flag count, rank-general (Baker Prop 4.3) -------
def bjorner_count(E_, flats_all, rank_of_, d):
    Fp = [f for f in flats_all if f and f != frozenset(E_)]
    counts = []
    for k in range(d + 1):
        if k == 0:
            counts.append(1); continue
        cnt = 0
        for chain in itertools.combinations(Fp, k):
            if not all(chain[i] < chain[i + 1] for i in range(k - 1)):
                continue
            if not all(rank_of_(chain[i]) == i + 1 for i in range(k)):
                continue
            mins = [min(f) for f in chain]
            if (all(mins[i] > mins[i + 1] for i in range(k - 1))
                    and mins[-1] > 0):
                cnt += 1
        counts.append(cnt)
    return counts

print("== Three routes to m ==")
mb, mu = mobius_chibar(E, flats, lambda f: rk_of[f])
assert mu[BLOCKS[0] | BLOCKS[1]] == -3 and mu[frozenset(E)] == 30
bj = bjorner_count(E, flats, lambda f: rk_of[f], 3)
assert mb == bj == [1, 7, 21, 30], (mb, bj)
print(f"  Mobius/chi-bar {mb} == Bjorner flags {bj} "
      f"(mu(H flat) = -3, mu(E) = +30, as predicted by hand)")

# uniform comparison on the same ground set: U(4,8)
uflats = [frozenset(c) for s in range(4)
          for c in itertools.combinations(E, s)] + [frozenset(E)]
urk = lambda f: min(len(f), 4)
umb, _ = mobius_chibar(E, uflats, urk)
ubj = bjorner_count(E, uflats, urk, 3)
assert umb == ubj == [1, 7, 21, 35], (umb, ubj)
print(f"  uniform U(4,8): Mobius {umb} == Bjorner {ubj}; "
      f"Vamos cost = {[u - v for u, v in zip(umb, mb)]} (top coeff only)")
SUMMARY["U48_m"] = umb

# ---------- the Chow ring and the full Kahler package ----------
CR, _a, _b, a_lift, b_lift, m = day14.analyse(
    "Vamos V8 via day14.analyse", E, flats, lambda f: rk_of[f],
    (1, 70, 70, 1), [1, 7, 21, 30])
assert m == mb == bj == [1, 7, 21, 30]
print("  V8: FY == Mobius == Bjorner == (1,7,21,30): TRIPLE AGREEMENT")
SUMMARY["V8"] = {"dims": [1, 70, 70, 1], "m": m,
                 "flats_by_rank_size": {f"{k}": v for k, v in
                                         sorted(by_rs.items())},
                 "circuits": {"size4": 5, "size5": 36},
                 "chi": "t^4-8t^3+28t^2-51t+30"}

print("== V8: Poincare duality, Hard Lefschetz, Hodge-Riemann ==")
n = 8
cfun = lambda S: len(S) * (n - len(S))
cnt = 0
for A_ in allsubs:
    for B_ in allsubs:
        if not (A_ <= B_ or B_ <= A_):
            assert cfun(A_ | B_) + cfun(A_ & B_) < cfun(A_) + cfun(B_)
            cnt += 1
print(f"  c(S)=|S|(8-|S|) strictly submodular: {cnt} incomparable pairs")
l_lift = {i: Fraction(cfun(f)) for f, i in CR.gidx.items()}

pair = [[Fraction(0)] * 70 for _ in range(70)]
for i in range(70):
    lift_i = {CR.basis[1][i]: Fraction(1)}
    for j in range(70):
        pair[i][j] = CR.deg(CR.mul(2, {j: Fraction(1)}, lift_i))
assert dense_rank(pair) == 70
print("  Poincare pairing A^1 x A^2 nondegenerate (rank 70): OK")

L12 = CR.op_matrix(1, l_lift)
assert dense_rank(L12) == 70
print("  Hard Lefschetz: l: A^1 -> A^2 is an isomorphism (rank 70): OK")
cur = {0: Fraction(1)}; dd = 0
for _ in range(3):
    cur = CR.mul(dd, cur, l_lift); dd += 1
deg_l3 = CR.deg(cur)
assert deg_l3 > 0
print(f"  deg(l^3) = {deg_l3} > 0 (HRR k=0): OK")

def l2_func(a_coords):
    step = CR.mul(1, a_coords, l_lift)
    return CR.deg(CR.mul(2, step, l_lift))
row = [l2_func({i: Fraction(1)}) for i in range(70)]
piv = next(i for i, v in enumerate(row) if v != 0)
prim = []
for j in range(70):
    if j == piv:
        continue
    v = {j: Fraction(1), piv: -row[j] / row[piv]}
    assert l2_func(v) == 0
    prim.append(v)
assert len(prim) == 69
G = [[Fraction(0)] * 69 for _ in range(69)]
for i in range(69):
    for j in range(69):
        lb = CR.mul(1, prim[j], l_lift)
        G[i][j] = -CR.deg(CR.mul(2, lb, CR.lift1(prim[i])))
assert all(G[i][j] == G[j][i] for i in range(69) for j in range(69))
print("  Gram of -deg(., l .) on P^1 (dim 69) symmetric: OK")
# exact LDL^T (no pivoting): pivots d_k = minor_k / minor_{k-1}
M = [r[:] for r in G]
pivots = []
for k in range(69):
    p = M[k][k]
    assert p != 0, f"LDL breaks at {k}"
    pivots.append(p)
    for i in range(k + 1, 69):
        f = M[i][k] / p
        if f:
            for j in range(k, 69):
                M[i][j] -= f * M[k][j]
assert all(p > 0 for p in pivots)
print(f"  Hodge-Riemann on P^1: LDL^T completes, all 69 pivots > 0 "
      f"(min {float(min(pivots)):.4g}, max {float(max(pivots)):.4g}) "
      f"-- Sylvester: positive definite")
try:
    import numpy as np
    ev = sorted(float(x) for x in np.linalg.eigvalsh(
        np.array([[float(v) for v in r] for r in G])))
    assert ev[0] > 0
    SUMMARY["V8"]["hrr_eigenvalues"] = [round(x, 4) for x in ev]
    print(f"  HRR eigenvalues (float): min {ev[0]:.3f}, max {ev[-1]:.3f}, "
          f"all > 0")
except ImportError:
    ev = None
SUMMARY["V8"].update({"deg_l3": int(deg_l3), "hl_rank": 70,
                      "pairing_rank": 70, "ldl_pivot_min": str(min(pivots)),
                      "ldl_pivot_max": str(max(pivots)),
                      "runtime_s": round(time.time() - t_start, 2)})

with open(os.path.join(HERE, "day16_summary.json"), "w") as fh:
    json.dump(SUMMARY, fh, indent=1)
print(f"\nALL DAY-16 ASSERTIONS PASSED in {time.time() - t_start:.1f}s")
