#!/usr/bin/env python3
"""Chapter 16 (Deep Dive 3) -- The Petersen Graph Has a Jacobian. DRAFT program, Phase 1.
Assembled 2026-11-26 from Day 17 (2026-10-27): Baker-Norine chip-firing
on the Petersen graph, the skeleton of tropical M_{0,5}.

Chapter 11 built the Petersen graph as G'''_{2,5} = T_5, the moduli of
5-leaf trees (vertices = 2-subsets of {0,...,4} = splits, edges =
disjoint pairs). Chapters 12-13 put divisor theory on the tropical
cubic's skeleton (the 6-cycle, then the full cubic graph). This chapter
runs the same machinery (chip_firing.py, Baker & Norine
arXiv:math/0608360, verified Day 8) on the Petersen graph itself,
viewed as a combinatorial tropical curve of genus
g = E - V + 1 = 15 - 10 + 1 = 6.

Headline claims (all asserted below):
  Jac(Petersen) ~= Z/2 x (Z/10)^3, |Jac| = 2000, computed three ways
  (Kirchhoff det; Smith normal form product; G-parking count), plus a
  fourth spectral confirmation added in this draft (Laplacian
  eigenvalues {0 x1, 2 x5, 5 x4}, so (2^5 * 5^4)/10 = 2000).
  Every vertex class [(v)-(q)] has order exactly 10 = the exponent --
  and, a Phase 1 strengthening beyond Day 17, every one of the 45
  pairwise vertex differences [(a)-(b)] has order exactly 10 too.
  Abel-Jacobi is injective (Petersen is bridgeless; BN Thm 1.8).
  Riemann-Roch holds at genus 6 (K = (1,...,1) by 3-regularity).

NOTE (Phase 1): imports ChipFiring from chip_firing.py in this
directory, as Day 17 did; the Phase 3 bundle program will be merged
self-contained, as Chs. 14-15 were. Writes ch16_summary.json next to
this file.
"""
import itertools
import json
import os
import random
from math import prod

import numpy as np
from sympy import Matrix

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

# ---------- the group: order three ways, structure from SNF ----------
tau = cf.spanning_tree_count()
print(f"spanning trees (Kirchhoff det) = {tau}")
assert tau == 2000, tau

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

# fourth route (draft addition): the Laplacian spectrum.
# Kirchhoff via eigenvalues: tau = (1/n) * product of nonzero eigenvalues.
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

# G-parking functions: q-reduced degree-0 reps are exactly the divisors
# with 0 <= D(v) <= deg(v)-1 = 2 off q satisfying Dhar's condition (P2).
others = [v for v in range(n) if v != q]
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

# ---------- Abel-Jacobi ----------
aj = {i: cf.aj(i) for i in range(n)}
assert np.all(aj[q] == 0)
keys = {tuple(int(x) for x in aj[i]) for i in range(n)}
assert len(keys) == n, "AJ must be injective on vertices (bridgeless)"
print("Abel-Jacobi S_q(v) = [(v)-(q)]: 10 distinct classes: OK")

# the stabiliser of q in S_5 has two orbits off q: the splits disjoint
# from q and the splits meeting q in one leaf (draft addition).
qv = set(verts[q])
disjoint = [v for v in verts if v != verts[q] and not set(v) & qv]
meeting = [v for v in verts if v != verts[q] and len(set(v) & qv) == 1]
assert len(disjoint) == 3 and len(meeting) == 6, (disjoint, meeting)
print(f"S_5-orbits off q: {len(disjoint)} disjoint splits + "
      f"{len(meeting)} splits meeting q in one leaf")

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

# draft strengthening (beyond Day 17): EVERY pairwise vertex difference
# [(a)-(b)], a != b, also has order exactly 10 -- all 45 unordered pairs.
pair_orders = {}
for a, b in itertools.combinations(range(n), 2):
    D = np.zeros(n, dtype=int); D[a] += 1; D[b] -= 1
    pair_orders[(verts[a], verts[b])] = class_order(cf.reduced_class(D))
assert len(pair_orders) == 45
assert all(v == 10 for v in pair_orders.values()), pair_orders
print("all 45 pairwise differences [(a)-(b)] have order exactly 10: OK")

# a concrete group-law identity: [(a)-(q)] + [(b)-(q)] computed two ways
a, b = 1, 2
lhs = cf.add_classes(aj[a], aj[b])
D = np.zeros(n, dtype=int); D[a] += 1; D[b] += 1; D[q] -= 2
rhs = cf.reduced_class(D)
assert np.array_equal(lhs, rhs)
print("group law: [a-q]+[b-q] == [(a)+(b)-2(q)] as q-reduced reps: OK")

# ---------- Riemann-Roch ----------
K = np.array([d - 2 for d in deg], dtype=int)     # 3-regular => all ones
assert np.all(K == 1) and int(K.sum()) == 2 * g - 2 == 10
assert cf.rank(np.zeros(n, dtype=int)) == 0
assert cf.rank(K) == g - 1 == 5
assert cf.rank(np.eye(n, dtype=int)[q]) == 0
print("hand checks: r(0)=0, r(K)=g-1=5, r((q))=0: OK")

rng = random.Random(17)
checked = 0
while checked < 40:
    D = np.array([rng.choice([-1, 0, 1]) for _ in range(n)], dtype=int)
    d = int(D.sum())
    if not (0 <= d <= 3):
        continue
    lhs = cf.rank(D) - cf.rank(K - D)
    assert lhs == d - g + 1, (D, lhs, d)
    checked += 1
print(f"Riemann-Roch r(D)-r(K-D) = deg D - g + 1 verified on {checked} "
      f"random divisors (entries in {{-1,0,1}}, deg 0..3): OK")

# ---------- a stabilization exhibit for the figure ----------
D0 = np.array([4, -2, 3, 0, -1, 2, 0, 1, -3, 2], dtype=int)
Dq, fired = cf.stabilize(D0, _trace=True)
print(f"stabilization exhibit: {D0.tolist()} -> {Dq.tolist()} "
      f"in {len(fired)} fired sets")
assert int(Dq.sum()) == int(D0.sum())
assert Dq.tolist() == [0, 2, 0, 1, 1, 1, 1, 0, 0, 0]

summary = {
    "V": n, "E": m, "genus": g,
    "spanning_trees": tau,
    "snf_diagonal": snf,
    "invariant_factors": invariants,
    "exponent": exponent,
    "laplacian_spectrum": {str(k): v for k, v in sorted(eig.items())},
    "parking_functions": len(parking),
    "aj_orders": {str(verts[i]): orders[i] for i in orders},
    "pairwise_difference_orders_all_10": True,
    "pairwise_differences_checked": len(pair_orders),
    "orbits_off_q": {"disjoint": len(disjoint), "meeting_in_one": len(meeting)},
    "rr_random_checked": checked,
    "stabilization": {"from": D0.tolist(), "to": [int(x) for x in Dq],
                      "fired_sets": len(fired)},
}
with open(os.path.join(HERE, "ch16_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
print("ALL CHAPTER-16 DRAFT ASSERTIONS PASSED")
