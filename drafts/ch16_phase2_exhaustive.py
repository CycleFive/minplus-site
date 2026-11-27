#!/usr/bin/env python3
"""Ch.16 Phase 2 — exhaustive verification of the referee-pass claims
(2026-11-27). Independent of the draft program's RR sample; uses the
lattice presentation Jac = Z^9 / L_red Z^9 throughout (exact integer
arithmetic mod det = 2000 via the adjugate), plus ChipFiring only for
cross-checks where noted.

Checks:
  A. Keys: the 2000 G-parking reps land in 2000 distinct lattice
     classes (Nyx G3).
  B. Exhaustive Riemann-Roch (Nyx E5): rank is a class function, so
     enumerate classes, not divisors. For every degree d in -1..11
     and every one of the 2000 classes (26,000 class representatives),
     r(D) is computed by the exact recursion
        r(D) >= 0  iff [D] has an effective representative,
        r(D) >= s+1 iff r(D-(v)) >= s for every vertex v
     (the second line is BN's definition with E = (v) + E'), organised
     by lattice keys. RR: r(D) - r(K-D) = d - g + 1 for all 26,000.
     Also: gonality (first d with some r >= 1) and max r by degree.
     DP ranks are cross-checked against ChipFiring.rank (BN definition,
     subset machinery) on the old seed-17 sample + random classes.
  C. Order distribution of Jac (Nyx G5): {1:1, 2:15, 5:124, 10:1860};
     |Jac[2]| = 16, |Jac[5]| = 125, matching (Z/2)^4 x (Z/5)^3.
  D. Route 4 on paper (Nyx N4): Petersen is srg(10,3,0,1);
     A^2 + A - 2I = J; adjacency spectrum {3:1, 1:5, -2:4}.
  E. Full-Laplacian SNF (Nyx N3): [1,1,1,1,1,2,10,10,10,0],
     so coker L = Z + Jac.
  F. E2 counterexample: two triangles joined by a bridge is its own
     2-core, has no pendant trees, yet AJ identifies the bridge's
     endpoints ((u)-(v) = L * 1_A is principal).
  G. Exhibit fired-set counts (Nyx G8): lex subset order 10 phase-2
     set-firings, reversed order 3, same q-reduced target.
  H. Group-law associativity on random triples (Duet D2): the draft's
     single identity is tautological for any well-defined reduction;
     associativity of add_classes is not.
  I. Exact det (Nyx N5): sympy det(L_red) = 2000 (the shared
     ChipFiring route is float round(abs(det)); it agrees here).
"""
import itertools
import random
from math import gcd

import numpy as np
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from chip_firing import ChipFiring

# ---------- graph ----------
verts = list(itertools.combinations(range(5), 2))
vidx = {v: i for i, v in enumerate(verts)}
edges = [(vidx[a], vidx[b]) for a, b in itertools.combinations(verts, 2)
         if not set(a) & set(b)]
n, m = len(verts), len(edges)
assert (n, m) == (10, 15)
q = vidx[(0, 1)]
others = [i for i in range(n) if i != q]
oi = {v: k for k, v in enumerate(others)}   # position in the 9-vector
g = m - n + 1
assert g == 6
cf = ChipFiring(n, edges, q)

L = sp.zeros(n)
deg = [0] * n
for i, j in edges:
    deg[i] += 1; deg[j] += 1
for i in range(n):
    L[i, i] = deg[i]
for i, j in edges:
    L[i, j] -= 1; L[j, i] -= 1
Lr = L.extract(others, others)
DET = int(Lr.det())
assert DET == 2000
print("I. exact sympy det(L_red) = 2000; float route agrees:",
      cf.spanning_tree_count() == 2000)

adj_mat = Lr.adjugate()          # Lr^{-1} = adj / DET, exact integers
adj_np = np.array(adj_mat.tolist(), dtype=np.int64)

def key_of_offvec(x9):
    """Lattice key of a degree-0 divisor given off q (9-vector):
    adj @ x mod DET in (Z/DET)^9. Equal keys iff same class."""
    y = adj_np @ np.array(x9, dtype=np.int64)
    return tuple(int(v) % DET for v in y)

def order_of_offvec(x9):
    """ord(x) = lcm of denominators of Lr^{-1} x = DET/gcd per coord."""
    y = adj_np @ np.array(x9, dtype=np.int64)
    o = 1
    for v in y:
        d = DET // gcd(DET, abs(int(v)))   # denom of v/DET, lowest terms
        o = o * d // gcd(o, d)
    return o

# ---------- parking reps (raw route, via ChipFiring's P1/P2) ----------
parking = []
for vals in itertools.product(range(3), repeat=n - 1):
    D = np.zeros(n, dtype=int)
    for v, val in zip(others, vals):
        D[v] = val
    if cf.is_q_reduced(D):
        D[q] = -int(D.sum())
        parking.append(tuple(int(x) for x in D))
assert len(parking) == 2000, len(parking)
keys = [key_of_offvec([D[v] for v in others]) for D in parking]
assert len(set(keys)) == 2000
ALLKEYS = set(keys)
print("A. 2000 G-parking reps -> 2000 distinct lattice classes: OK (G3)")

# ---------- C. order distribution ----------
from collections import Counter
dist = Counter(order_of_offvec([D[v] for v in others]) for D in parking)
assert dist == {1: 1, 2: 15, 5: 124, 10: 1860}, dist
jac2 = sum(1 for D in parking
           if order_of_offvec([D[v] for v in others]) in (1, 2))
jac5 = sum(1 for D in parking
           if order_of_offvec([D[v] for v in others]) in (1, 5))
assert jac2 == 16 and jac5 == 125, (jac2, jac5)
print(f"C. order distribution {dict(sorted(dist.items()))}; "
      f"|Jac[2]|=16=(Z/2)^4, |Jac[5]|=125=(Z/5)^3; "
      f"1860/2000 = {1860/2000:.2f} have order 10: OK (G5)")

# ---------- D. srg / spectrum on paper ----------
A = np.zeros((n, n), dtype=int)
for i, j in edges:
    A[i, j] = A[j, i] = 1
adjset = {v: set() for v in range(n)}
for i, j in edges:
    adjset[i].add(j); adjset[j].add(i)
for a, b in itertools.combinations(range(n), 2):
    common = len(adjset[a] & adjset[b])
    assert common == (0 if b in adjset[a] else 1), (a, b, common)
J = np.ones((n, n), dtype=int)
assert np.array_equal(A @ A + A - 2 * np.eye(n, dtype=int), J)
print("D. srg(10,3,0,1): common neighbours 0 (adjacent) / 1 (not); "
      "A^2 + A - 2I = J: OK")
evA = sp.Matrix(A.tolist()).eigenvals()
specA = {int(k): int(v) for k, v in evA.items()}
assert specA == {3: 1, 1: 5, -2: 4}, specA
# on paper: eigenvalues of A off the all-ones direction satisfy
# l^2 + l - 2 = 0 (J acts as 0 there), so l in {1, -2}; multiplicities
# from trace A = 0 and trace A^2 = sum deg = 30:
#   3 + f - 2g' = 0, f + g' = 9  ->  f = 5, g' = 4.
assert 3 + 5 - 2 * 4 == 0 and 5 + 4 == 9 and 9 + 5 * 1 + 4 * 4 == 30
print("D. adjacency spectrum {3:x1, 1:x5, -2:x4} (mults from traces); "
      "L = 3I - A -> {0:x1, 2:x5, 5:x4}: OK (N4)")

# ---------- E. full-Laplacian SNF ----------
Df = smith_normal_form(L, domain=ZZ)
diagf = sorted(abs(int(Df[i, i])) for i in range(n))  # signs are convention
assert diagf == [0, 1, 1, 1, 1, 1, 2, 10, 10, 10], diagf
print("E. SNF(full L) = [1x5, 2, 10, 10, 10, 0] -> coker L = Z + Jac: "
      "OK (N3)")

# ---------- B. exhaustive Riemann-Roch by classes ----------
# columns: key contribution of one chip at v (v != q); q contributes 0.
col = {}
for v in range(n):
    if v == q:
        col[v] = (0,) * 9
    else:
        x = [0] * 9; x[oi[v]] = 1
        col[v] = key_of_offvec(x)
def kadd(a, b):
    return tuple((x + y) % DET for x, y in zip(a, b))
def ksub(a, b):
    return tuple((x - y) % DET for x, y in zip(a, b))

MAXD = 11
Eff = [{(0,) * 9}]
for d in range(1, MAXD + 1):
    s = set()
    for k in Eff[-1]:
        for v in range(n):
            s.add(kadd(k, col[v]))
    Eff.append(s)
assert Eff[0] == {(0,) * 9}
# R[d] = list of sets: R[d][s] = keys of degree-d classes with r >= s.
R = {}
for d in range(0, MAXD + 1):
    levels = [set(Eff[d])]
    s = 0
    while s < d:                       # r(D) <= deg D, so s+1 <= d
        prev = levels[-1]
        if not prev:
            break
        # classes at degree d-1 with r >= s; empty beyond its max rank
        lower = R[d - 1][s] if s < len(R[d - 1]) else set()
        nxt = {k for k in prev
               if all(ksub(k, col[v]) in lower for v in range(n))}
        if not nxt:
            break
        levels.append(nxt)
        s += 1
    R[d] = levels

def rank_of(d, k):
    if d < 0:
        return -1
    for s in range(len(R[d]) - 1, -1, -1):
        if k in R[d][s]:
            return s
    return -1

Kdiv = np.array([dv - 2 for dv in deg], dtype=int)
assert Kdiv.tolist() == [1] * 10
kK = key_of_offvec([int(Kdiv[v]) for v in others])  # key of K - 10(q)
# (K - 10 e_q) off q is just K|others since e_q drops out. ✓

checks = 0
for d in range(-1, MAXD + 1):
    for k in ALLKEYS:
        lhs = rank_of(d, k) - rank_of(10 - d, ksub(kK, k))
        assert lhs == d - g + 1, (d, k, lhs)
        checks += 1
assert checks == 13 * 2000 == 26000, checks
print(f"B. EXHAUSTIVE RR: r(D)-r(K-D) = deg D - g + 1 for all {checks} "
      f"class representatives (degrees -1..11 x 2000 classes): OK (E5)")

maxr = {d: max((rank_of(d, k) for k in ALLKEYS), default=-1)
        for d in range(0, MAXD + 1)}
print("B. max r by degree:", {d: maxr[d] for d in range(0, 12)})
assert [maxr[d] for d in (4, 5, 6, 7, 8)] == [1, 1, 2, 2, 3], maxr
gon = min(d for d in range(1, MAXD + 1) if maxr[d] >= 1)
assert gon == 4, gon
assert rank_of(0, (0,) * 9) == 0
assert rank_of(1, (0,) * 9) == 0          # r((q)) = 0: class key 0
assert rank_of(10, kK) == g - 1 == 5      # r(K) = g-1
print("B. gonality 4 (first r>=1 at degree 4); max r at degrees "
      "4..8 = 1,1,2,2,3; r(0)=0, r((q))=0, r(K)=5: OK")

# cross-check DP ranks against ChipFiring.rank (BN definition)
rng = random.Random(17)
sample = []
while len(sample) < 40:   # replay the draft program's seed-17 sample
    D = np.array([rng.choice([-1, 0, 1]) for _ in range(n)], dtype=int)
    if 0 <= int(D.sum()) <= 3:
        sample.append(D)
rng2 = random.Random(99)
for _ in range(30):
    rep = np.array(parking[rng2.randrange(2000)], dtype=int)
    d = rng2.randrange(0, 9)
    D = rep.copy(); D[q] += d
    sample.append(D)
for D in sample:
    d = int(D.sum())
    # key of D - d*(q): off-q coordinates are D[v] for v != q.
    k = key_of_offvec([int(D[v]) for v in others])
    assert cf.rank(D) == rank_of(d, k), (D, cf.rank(D), rank_of(d, k))
print(f"B. DP ranks agree with ChipFiring.rank on {len(sample)} divisors "
      f"(seed-17 sample + random classes, degrees 0..8): OK")

# ---------- F. E2 counterexample ----------
e2_edges = [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)]
def two_core(nn, ee):
    alive = set(range(nn)); el = list(ee); changed = True
    while changed:
        changed = False
        dv = {v: 0 for v in alive}
        for i, j in el:
            if i in alive and j in alive:
                dv[i] += 1; dv[j] += 1
        for v in list(alive):
            if dv[v] < 2:
                alive.discard(v); changed = True
    return alive
assert two_core(6, e2_edges) == set(range(6)), "graph is its own 2-core"
cf2 = ChipFiring(6, e2_edges, 0)
# (2)-(3) is principal: L @ 1_{0,1,2} = e_2 - e_3
f = np.array([1, 1, 1, 0, 0, 0], dtype=int)
assert (cf2.lap @ f).tolist() == [0, 0, 1, -1, 0, 0]
assert np.array_equal(cf2.aj(2), cf2.aj(3)), "AJ must identify 2 and 3"
print("F. two triangles + bridge: own 2-core, no pendant trees, bridge "
      "exists, yet S_q(2) = S_q(3) ((2)-(3) = L*1_A principal): "
      "BN 1.8 needs 2-edge-connected, not 2-core: OK (E2)")

# ---------- G. exhibit fired-set counts ----------
D0 = np.array([4, -2, 3, 0, -1, 2, 0, 1, -3, 2], dtype=int)
tgt, fired_lex = cf.stabilize(D0, _trace=True)
tgt2, fired_rev = cf.stabilize(D0, subset_order="reversed", _trace=True)
assert tgt.tolist() == [0, 2, 0, 1, 1, 1, 1, 0, 0, 0]
assert np.array_equal(tgt, tgt2)
print(f"G. exhibit: lex order {len(fired_lex)} phase-2 set-firings, "
      f"reversed order {len(fired_rev)}, same q-reduced target: "
      f"the count is algorithm-dependent: OK (G8)")
assert len(fired_lex) == 10 and len(fired_rev) == 3

# ---------- H. associativity ----------
rng3 = random.Random(7)
reps = [np.array(D, dtype=int) for D in parking]
for _ in range(100):
    a, b, c = (reps[rng3.randrange(2000)] for _ in range(3))
    lhs = cf.add_classes(cf.add_classes(a, b), c)
    rhs = cf.add_classes(a, cf.add_classes(b, c))
    assert np.array_equal(lhs, rhs)
print("H. add_classes associative on 100 random triples of classes: "
      "OK (D2)")

print("\nALL CH.16 PHASE-2 EXHAUSTIVE / REFEREE CHECKS PASSED")
