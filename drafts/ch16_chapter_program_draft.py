#!/usr/bin/env python3
"""Chapter 16 (Deep Dive 3) -- The Petersen Graph Has a Jacobian.
Chapter program, Phase 2 (2026-11-27): re-verified end-to-end from the
draft text, with the #minplus-book referee fixes (Nyx E1-E5, G1-G8,
N1-N11; Duet D1-D4) applied to the code:

  * Route 1 determinant is now EXACT (sympy integer det); the shared
    ChipFiring float route is kept as a cross-check (N5).
  * Full-Laplacian SNF asserted: [1,1,1,1,1,2,10,10,10,0], so
    coker L = Z + Jac (N3).
  * Route 4 done on paper: Petersen is srg(10,3,0,1), so
    A^2 + A - 2I = J, adjacency spectrum {3,1^5,(-2)^4} (multiplicities
    from the two trace identities), L = 3I - A gives {0,2^5,5^4} (N4).
  * The 2000 G-parking reps land in 2000 DISTINCT lattice classes
    (key = adj(L_red) x mod 2000), as BN Prop 3.1 requires (G3).
  * Stabiliser orbits by brute force: Stab(q) has order 12, vertex
    orbits off q are 3 (disjoint = N(q)) + 6 (meeting in one); S_5
    orbits on the 45 pairs are 15 adjacent + 30 non-adjacent (E4, D3).
  * Order distribution of Jac asserted: {1:1, 2:15, 5:124, 10:1860},
    |Jac[2]| = 16, |Jac[5]| = 125 (G5).
  * Group law: the draft's single identity is joined by associativity
    on 100 random triples -- the identity alone holds for ANY
    well-defined reduction, so it tests code, not law (D2).
  * Riemann-Roch is EXHAUSTIVE at class level (E5): rank is a class
    function, so all 2000 classes x degrees -1..11 = 26,000 class
    representatives are checked, ranks computed by the exact key DP
    (r >= s+1 iff r(D-(v)) >= s for all v) and cross-checked against
    ChipFiring.rank (BN's definition) on the old seed-17 sample.
    Gonality 4; max r by degree 4..8 = 1,1,2,2,3.
  * The stabilisation exhibit reports both subset orders: 10 phase-2
    set-firings in lex order, 3 in reversed order, same target (G8).

Headline claims (all asserted below):
  Jac(Petersen) ~= Z/2 x (Z/10)^3, |Jac| = 2000 -- order by two
  independent routes (exact Kirchhoff det; G-parking enumeration),
  structure from the Smith normal form, spectrum as corroboration.
  Every vertex class [(v)-(q)] has order exactly 10 = the exponent;
  by vertex-transitivity the same holds for all 45 pairwise
  differences [(a)-(b)] (a corollary, kept as a code cross-check).
  Abel-Jacobi is injective (Petersen is 2-edge-connected; BN 1.8).
  Riemann-Roch holds at genus 6, exhaustively (K = (1,...,1)).

NOTE: imports ChipFiring from chip_firing.py in this directory, as
Day 17 did; the Phase 3 bundle program will be merged self-contained,
as Chs. 14-15 were. Writes ch16_summary.json next to this file.
"""
import itertools
import json
import os
import random
from collections import Counter
from math import gcd, prod

import numpy as np
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from chip_firing import ChipFiring

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------- the Petersen graph, built the Chapter 11 way ----------
verts = list(itertools.combinations(range(5), 2))   # 2-subsets = splits
vidx = {v: i for i, v in enumerate(verts)}
edges = []
for a, b in itertools.combinations(verts, 2):
    if not set(a) & set(b):
        edges.append((vidx[a], vidx[b]))
n, m = len(verts), len(edges)
assert (n, m) == (10, 15)
deg = [0] * n
for i, j in edges:
    deg[i] += 1
    deg[j] += 1
assert all(d == 3 for d in deg), "Petersen is 3-regular"
g = m - n + 1
assert g == 6
print(f"Petersen graph: V={n}, E={m}, 3-regular, genus g={g}")

# bridgeless (2-edge-connected): removing any single edge keeps it connected
def connected(edge_list):
    adj = {v: set() for v in range(n)}
    for i, j in edge_list:
        adj[i].add(j); adj[j].add(i)
    seen, stack = set(), [0]
    while stack:
        u = stack.pop()
        if u not in seen:
            seen.add(u); stack.extend(adj[u])
    return len(seen) == n
assert connected(edges)
for e in edges:
    assert connected([x for x in edges if x != e]), f"edge {e} is a bridge"
print("bridgeless (2-edge-connected): OK -- BN Thm 1.8 predicts S_q^1 injective")

q = 0                       # vertex {0,1} (0-based; {1,2} in 1-based labels)
cf = ChipFiring(n, edges, q)
assert cf.genus == 6
others = [v for v in range(n) if v != q]
oi = {v: k for k, v in enumerate(others)}

# ---------- the group: Jac = coker(L_red), order + structure ----------
# Div^0 ~= Z^{V\q} by dropping the q-coordinate; principal divisors map
# to L_red Z^9 under the gauge f(q) = 0 (constants are ker L). So
# Jac ~= Z^9 / L_red Z^9: |det| is the order, the SNF the structure.
Lfull = Matrix(cf.lap.tolist())
Lr = Lfull.extract(others, others)

tau_exact = int(Lr.det())
tau_float = cf.spanning_tree_count()   # shared route: round(abs(det))
print(f"spanning trees (Kirchhoff det, exact) = {tau_exact} "
      f"(float route agrees: {tau_float == tau_exact})")
assert tau_exact == 2000 and tau_float == 2000
tau = tau_exact

snf = cf.snf_diagonal()
print(f"SNF of reduced Laplacian = {snf}")
assert snf == [1, 1, 1, 1, 1, 2, 10, 10, 10]
assert prod(snf) == tau
invariants = [d for d in snf if d > 1]
print(f"Jac(Petersen) ~= {' x '.join('Z/' + str(d) + 'Z' for d in invariants)}")
assert invariants == [2, 10, 10, 10]
assert all(invariants[i + 1] % invariants[i] == 0
           for i in range(len(invariants) - 1)), "invariant factors must divide"
exponent = invariants[-1]
assert exponent == 10
print(f"invariant factors {invariants}, exponent = {exponent}")

# the full Laplacian's SNF has one extra entry, 0 (N3):
Df = smith_normal_form(Lfull, domain=ZZ)
diagf = sorted(abs(int(Df[i, i])) for i in range(n))
assert diagf == [0, 1, 1, 1, 1, 1, 2, 10, 10, 10], diagf
print("SNF(full L) = [1,1,1,1,1,2,10,10,10,0] -> coker L = Z + Jac: OK")

# spectral corroboration: the Laplacian spectrum, computed exactly.
eig = Matrix(cf.lap.tolist()).eigenvals()
eig = {int(k): int(v) for k, v in eig.items()}
print(f"Laplacian spectrum (exact, sympy) = {eig}")
assert eig == {0: 1, 2: 5, 5: 4}, eig
tau_spec = 1
for val, mult in eig.items():
    if val != 0:
        tau_spec *= val ** mult
assert tau_spec % n == 0 and tau_spec // n == tau
print(f"spectral Kirchhoff: (2^5 * 5^4)/10 = {tau_spec // n} == det route: OK")

# Route 4 on paper (N4): Petersen is strongly regular srg(10,3,0,1) --
# adjacent vertices have 0 common neighbours, non-adjacent have 1 --
# so A^2 + A - 2I = J, the adjacency eigenvalues off the all-ones
# direction satisfy l^2 + l - 2 = 0, and the multiplicities follow
# from trace A = 0 and trace A^2 = sum deg = 30.
A = np.zeros((n, n), dtype=int)
adjset = {v: set() for v in range(n)}
for i, j in edges:
    A[i, j] = A[j, i] = 1
    adjset[i].add(j); adjset[j].add(i)
for a, b in itertools.combinations(range(n), 2):
    assert len(adjset[a] & adjset[b]) == (0 if b in adjset[a] else 1)
assert np.array_equal(A @ A + A - 2 * np.eye(n, dtype=int),
                      np.ones((n, n), dtype=int))
specA = {int(k): int(v) for k, v in Matrix(A.tolist()).eigenvals().items()}
assert specA == {3: 1, 1: 5, -2: 4}, specA
print("srg(10,3,0,1): A^2+A-2I = J; adjacency spectrum {3:x1, 1:x5, "
      "-2:x4} -> L = 3I-A spectrum as above: OK")

# lattice keys: class of a degree-0 divisor x (off q) is keyed by
# adj(L_red) @ x mod 2000 in (Z/2000)^9; equal keys iff same class.
adj_np = np.array(Lr.adjugate().tolist(), dtype=np.int64)
def key_of_offvec(x9):
    y = adj_np @ np.array(x9, dtype=np.int64)
    return tuple(int(v) % tau for v in y)
def order_of_offvec(x9):
    y = adj_np @ np.array(x9, dtype=np.int64)
    o = 1
    for v in y:
        d = tau // gcd(tau, abs(int(v)))
        o = o * d // gcd(o, d)
    return o

# G-parking functions: q-reduced degree-0 reps are exactly the divisors
# with 0 <= D(v) <= deg(v)-1 = 2 off q satisfying Dhar's condition (P2).
# BN Prop 3.1: exactly one q-reduced rep per class, so the count is
# |Jac| -- and the reps must land in distinct classes (checked via keys).
parking = []
for vals in itertools.product(range(3), repeat=n - 1):
    D = np.zeros(n, dtype=int)
    for v, val in zip(others, vals):
        D[v] = val
    if cf.is_q_reduced(D):
        D[q] = -int(D.sum())
        parking.append(tuple(int(x) for x in D))
print(f"G-parking functions (q-reduced deg-0 reps) = {len(parking)}")
assert len(parking) == tau, "parking count must equal |Jac|"
pkeys = [key_of_offvec([D[v] for v in others]) for D in parking]
assert len(set(pkeys)) == tau, "reps must land in distinct classes (G3)"
ALLKEYS = set(pkeys)
print("the 2000 reps land in 2000 distinct lattice classes: OK")

# ---------- Abel-Jacobi ----------
aj = {i: cf.aj(i) for i in range(n)}
assert np.all(aj[q] == 0)
keys = {tuple(int(x) for x in aj[i]) for i in range(n)}
assert len(keys) == n, "AJ must be injective on vertices (bridgeless)"
print("Abel-Jacobi S_q(v) = [(v)-(q)]: 10 distinct classes: OK")

# ---------- orbits, by brute-force permutations (E4) ----------
def perm_vertex(p, v):
    return vidx[tuple(sorted(p[x] for x in verts[v]))]
perms = list(itertools.permutations(range(5)))
stab = [p for p in perms if perm_vertex(p, q) == q]
assert len(stab) == 12          # 2! * 3!
def orbits(items, group, act):
    seen, orbs = set(), []
    for it in items:
        if it in seen:
            continue
        orb = {act(p, it) for p in group}
        seen |= orb
        orbs.append(orb)
    return orbs
vo = orbits(others, stab, perm_vertex)
assert sorted(len(o) for o in vo) == [3, 6]
qv = set(verts[q])
disjoint_idx = {vidx[v] for v in verts if v != verts[q] and not set(v) & qv}
assert disjoint_idx in vo and len(disjoint_idx) == 3
assert disjoint_idx == adjset[q], "the 3 disjoint splits are exactly N(q)"
meeting_idx = {v for v in others if v not in disjoint_idx}
assert meeting_idx in vo and len(meeting_idx) == 6
print("Stab(q) (order 12) orbits off q: 3 disjoint splits (= N(q)) + "
      "6 splits meeting q in one leaf: OK")
po = orbits(list(itertools.combinations(range(n), 2)), perms,
            lambda p, pr: tuple(sorted((perm_vertex(p, pr[0]),
                                        perm_vertex(p, pr[1])))))
assert sorted(len(o) for o in po) == [15, 30]
adj_pairs = {pr for pr in itertools.combinations(range(n), 2)
             if not set(verts[pr[0]]) & set(verts[pr[1]])}
assert adj_pairs in po and len(adj_pairs) == 15
print("S_5 orbits on the 45 pairs: 15 adjacent + 30 non-adjacent: OK "
      "(distance-transitivity, not just the S_5 action)")

def class_order(cls):
    """Order of a q-reduced representative, by repeated addition."""
    cur, k = np.array(cls, dtype=int).copy(), 1
    while not np.all(cur == 0):
        cur = cf.add_classes(cur, cls); k += 1
        assert k <= exponent, "order cannot exceed the exponent"
    return k

# orders of the vertex classes; the group exponent must be attained
orders = {}
for i in range(n):
    if i == q:
        continue
    orders[i] = class_order(aj[i])
    assert exponent % orders[i] == 0
print(f"orders of [(v)-(q)]: { {verts[i]: orders[i] for i in orders} }")
assert all(v == 10 for v in orders.values()), "every vertex class has order 10"
assert max(orders.values()) == exponent, "exponent must be attained by an element"

# corollary of vertex-transitivity (E1), kept as a code cross-check:
# every pairwise difference [(a)-(b)] also has order exactly 10.
pair_orders = {}
for a, b in itertools.combinations(range(n), 2):
    D = np.zeros(n, dtype=int); D[a] += 1; D[b] -= 1
    pair_orders[(verts[a], verts[b])] = class_order(cf.reduced_class(D))
assert len(pair_orders) == 45
assert all(v == 10 for v in pair_orders.values()), pair_orders
print("all 45 pairwise differences [(a)-(b)] have order exactly 10: OK "
      "(corollary of the nine classes by vertex-transitivity)")

# the group's order distribution (G5): most elements have order 10,
# so the content is that NO vertex difference is low-order torsion.
dist = Counter(order_of_offvec([D[v] for v in others]) for D in parking)
assert dist == {1: 1, 2: 15, 5: 124, 10: 1860}, dist
jac2 = sum(1 for D in parking
           if order_of_offvec([D[v] for v in others]) in (1, 2))
jac5 = sum(1 for D in parking
           if order_of_offvec([D[v] for v in others]) in (1, 5))
assert jac2 == 16 and jac5 == 125
print(f"order distribution of Jac: {dict(sorted(dist.items()))} "
      f"(93% order 10); |Jac[2]| = 16, |Jac[5]| = 125: OK")

# a concrete group-law identity: [(a)-(q)] + [(b)-(q)] computed two ways.
# NOTE (D2): this identity holds for ANY well-defined reduction, so it
# tests the reduction code; the law itself is exercised by the
# associativity check below.
a, b = 1, 2
lhs = cf.add_classes(aj[a], aj[b])
D = np.zeros(n, dtype=int); D[a] += 1; D[b] += 1; D[q] -= 2
rhs = cf.reduced_class(D)
assert np.array_equal(lhs, rhs)
print("group law: [a-q]+[b-q] == [(a)+(b)-2(q)] as q-reduced reps: OK")
rngA = random.Random(7)
reps = [np.array(D, dtype=int) for D in parking]
for _ in range(100):
    x, y, z = (reps[rngA.randrange(tau)] for _ in range(3))
    assert np.array_equal(cf.add_classes(cf.add_classes(x, y), z),
                          cf.add_classes(x, cf.add_classes(y, z)))
print("group law: add_classes associative on 100 random triples: OK")

# ---------- Riemann-Roch, exhaustive at class level (E5) ----------
K = np.array([d - 2 for d in deg], dtype=int)     # 3-regular => all ones
assert np.all(K == 1) and int(K.sum()) == 2 * g - 2 == 10
assert cf.rank(np.zeros(n, dtype=int)) == 0
assert cf.rank(K) == g - 1 == 5
assert cf.rank(np.eye(n, dtype=int)[q]) == 0
print("direct computations: r(0)=0, r(K)=g-1=5, r((q))=0: OK")

# Rank is a class function: enumerate classes, not divisors. Keys as
# above; effective keys per degree by multiset sums of column keys;
# ranks by the recursion r(D) >= s+1 iff r(D-(v)) >= s for all v
# (BN's definition with E = (v) + E'). RR is then checked for every
# class in every degree -1..11: 2000 x 13 = 26,000 representatives.
col = {}
for v in range(n):
    if v == q:
        col[v] = (0,) * 9
    else:
        x = [0] * 9; x[oi[v]] = 1
        col[v] = key_of_offvec(x)
def kadd(aa, bb):
    return tuple((x + y) % tau for x, y in zip(aa, bb))
def ksub(aa, bb):
    return tuple((x - y) % tau for x, y in zip(aa, bb))
MAXD = 11
Eff = [{(0,) * 9}]
for d in range(1, MAXD + 1):
    Eff.append({kadd(k, col[v]) for k in Eff[-1] for v in range(n)})
R = {}
for d in range(0, MAXD + 1):
    levels = [set(Eff[d])]
    s = 0
    while s < d:
        if not levels[-1]:
            break
        lower = R[d - 1][s] if s < len(R[d - 1]) else set()
        nxt = {k for k in levels[-1]
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
kK = key_of_offvec([int(K[v]) for v in others])
rr_checks = 0
for d in range(-1, MAXD + 1):
    for k in ALLKEYS:
        assert rank_of(d, k) - rank_of(10 - d, ksub(kK, k)) == d - g + 1
        rr_checks += 1
assert rr_checks == 26000, rr_checks
print(f"Riemann-Roch EXHAUSTIVE: identity holds for all {rr_checks} "
      f"class representatives (degrees -1..11 x 2000 classes): OK")
maxr = {d: max(rank_of(d, k) for k in ALLKEYS) for d in range(MAXD + 1)}
assert [maxr[d] for d in (4, 5, 6, 7, 8)] == [1, 1, 2, 2, 3], maxr
assert min(d for d in range(1, MAXD + 1) if maxr[d] >= 1) == 4
print(f"max r by degree 0..11: {[maxr[d] for d in range(12)]}; "
      f"divisorial gonality 4 (first r>=1 at degree 4): OK")
# cross-check the DP ranks against ChipFiring.rank (BN's definition,
# subset machinery) on the Phase 1 seed-17 sample:
rng = random.Random(17)
sample = []
while len(sample) < 40:
    D = np.array([rng.choice([-1, 0, 1]) for _ in range(n)], dtype=int)
    if 0 <= int(D.sum()) <= 3:
        sample.append(D)
for D in sample:
    d = int(D.sum())
    assert cf.rank(D) == rank_of(d, key_of_offvec([int(D[v]) for v in others]))
print("DP ranks agree with ChipFiring.rank on the 40-divisor seed-17 "
      "sample (now a cross-check, not the RR evidence): OK")

# ---------- a stabilization exhibit for the figure (G8) ----------
D0 = np.array([4, -2, 3, 0, -1, 2, 0, 1, -3, 2], dtype=int)
Dq, fired_lex = cf.stabilize(D0, _trace=True)
Dq2, fired_rev = cf.stabilize(D0, subset_order="reversed", _trace=True)
print(f"stabilization exhibit: {D0.tolist()} -> {Dq.tolist()} "
      f"in {len(fired_lex)} phase-2 set-firings (lex order; "
      f"{len(fired_rev)} in reversed order -- the count is "
      f"algorithm-dependent, the target is not)")
assert int(Dq.sum()) == int(D0.sum())
assert Dq.tolist() == [0, 2, 0, 1, 1, 1, 1, 0, 0, 0]
assert np.array_equal(Dq, Dq2) and len(fired_lex) == 10 and len(fired_rev) == 3

summary = {
    "V": n, "E": m, "genus": g,
    "spanning_trees": tau,
    "snf_diagonal": snf,
    "invariant_factors": invariants,
    "exponent": exponent,
    "laplacian_spectrum": {str(k): v for k, v in sorted(eig.items())},
    "parking_functions": len(parking),
    "parking_distinct_classes": len(set(pkeys)),
    "aj_orders": {str(verts[i]): orders[i] for i in orders},
    "pairwise_difference_orders_all_10": True,
    "pairwise_differences_checked": len(pair_orders),
    "order_distribution": {str(k): v for k, v in sorted(dist.items())},
    "stabiliser_order": len(stab),
    "orbits_off_q": {"disjoint": 3, "meeting_in_one": 6},
    "pair_orbits": {"adjacent": 15, "non_adjacent": 30},
    "rr_exhaustive_classes": rr_checks,
    "gonality": 4,
    "max_rank_by_degree": {str(d): maxr[d] for d in range(12)},
    "stabilization": {"from": D0.tolist(), "to": [int(x) for x in Dq],
                      "fired_sets_lex": len(fired_lex),
                      "fired_sets_reversed": len(fired_rev)},
}
with open(os.path.join(HERE, "ch16_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
print("ALL CHAPTER-16 PHASE-2 ASSERTIONS PASSED")
