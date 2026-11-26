#!/usr/bin/env python3
"""Ch.16 Phase 2 independent verification (2026-11-26).

Checks the Chapter 16 draft's claims WITHOUT the book's chip_firing.py:
  * graph census, genus, bridgelessness (direct BFS);
  * Kirchhoff det, Smith diagonal, spectrum -- sympy exact integer
    arithmetic on the Laplacian built here;
  * Jacobian group and element orders via the lattice presentation
    Jac = Z^9 / L_red Z^9: k*x is in the lattice iff L_red^{-1}(k*x)
    is integral, so ord(x) is the lcm of the denominators of
    L_red^{-1} x -- no chip-firing, no reduction algorithm;
  * the draft's new claim: all 45 pairwise vertex differences
    [(a)-(b)] have order exactly 10; the nine [(v)-(q)] likewise;
  * stabiliser-of-q orbits on vertices (3+6) and S5 orbits on
    unordered pairs (adjacent / non-adjacent) by brute permutations;
  * G-parking count from the raw P2 definition over all subsets;
  * Try It 1-4 answers;
  * stabilisation exhibit: target is q-reduced (raw P2) and
    source - target lies in the full Laplacian lattice (integer
    solution of L x = source - target exists), degree preserved.
Riemann-Roch is NOT re-derived here; it stands on the draft
program re-run (ch16-phase2-draft-rerun.log, all assertions passed).
"""
import itertools
from math import gcd, lcm

import sympy as sp

# ---- graph, built from splits as in Ch. 11 ----
verts = list(itertools.combinations(range(5), 2))
vidx = {v: i for i, v in enumerate(verts)}
edges = [(vidx[a], vidx[b]) for a, b in itertools.combinations(verts, 2)
         if not set(a) & set(b)]
n, m = len(verts), len(edges)
assert (n, m) == (10, 15)
deg = [0] * n
adj = {i: set() for i in range(n)}
for i, j in edges:
    deg[i] += 1; deg[j] += 1; adj[i].add(j); adj[j].add(i)
assert all(d == 3 for d in deg)
g = m - n + 1
assert g == 6
print("census V=10 E=15, 3-regular, genus 6: OK")

def connected(edge_list):
    a = {i: set() for i in range(n)}
    for i, j in edge_list:
        a[i].add(j); a[j].add(i)
    seen, stack = set(), [0]
    while stack:
        u = stack.pop()
        if u not in seen:
            seen.add(u); stack.extend(a[u])
    return len(seen) == n
assert connected(edges) and all(
    connected([e for e in edges if e != e0]) for e0 in edges)
print("bridgeless: removing any of the 15 edges keeps it connected: OK")

q = vidx[(0, 1)]
others = [i for i in range(n) if i != q]

# ---- Laplacian routes (sympy, exact) ----
L = sp.zeros(n)
for i in range(n):
    L[i, i] = deg[i]
for i, j in edges:
    L[i, j] -= 1; L[j, i] -= 1
Lr = L.extract([i for i in range(n) if i != q],
               [i for i in range(n) if i != q])
det = int(Lr.det())
assert det == 2000, det
print("Route 1 Kirchhoff: det(reduced Laplacian) = 2000: OK")

from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
D = smith_normal_form(Lr, domain=ZZ)
diag = [int(D[i, i]) for i in range(9)]
assert diag == [1, 1, 1, 1, 1, 2, 10, 10, 10], diag
inv = [d for d in diag if d > 1]
assert inv == [2, 10, 10, 10]
assert all(inv[i + 1] % inv[i] == 0 for i in range(len(inv) - 1))
prod = 1
for d in diag: prod *= d
assert prod == 2000
print("Route 2 SNF: diag [1x5,2,10,10,10], Jac = Z/2 x (Z/10)^3, "
      "order 2000, exponent 10: OK")

ev = L.eigenvals()
spec = {int(k): v for k, v in ev.items()}
assert spec == {0: 1, 2: 5, 5: 4}, spec
assert (2**5 * 5**4) // 10 == 2000
print("Route 4 spectrum {0:x1, 2:x5, 5:x4}; (2^5*5^4)/10 = 2000: OK")

# ---- orders via the lattice presentation, exact rationals ----
# x in Z^9 represents a class in Z^9 / Lr Z^9. k*x lies in the lattice
# iff Lr^{-1}(k*x) is integral, so ord(x) = lcm of the denominators
# of Lr^{-1} x (Lr symmetric, invertible over Q).
Lr_inv = Lr.inv()

def order_of_full_divisor(vec):
    """vec: degree-0 divisor on all 10 vertices; drop q coordinate."""
    x = sp.Matrix([int(vec[i]) for i in others])
    y = Lr_inv * x
    ord_ = 1
    for v in y:
        ord_ = lcm(ord_, int(sp.denom(v)))
    return ord_

orders_q = {}
for v in others:
    vec = [0] * n; vec[v] = 1; vec[q] = -1
    orders_q[v] = order_of_full_divisor(vec)
assert all(o == 10 for o in orders_q.values()), orders_q
print(f"nine [(v)-(q)] orders, lattice method: all exactly 10: OK")

pair_orders = {}
for a, b in itertools.combinations(range(n), 2):
    vec = [0] * n; vec[a] = 1; vec[b] = -1
    pair_orders[(a, b)] = order_of_full_divisor(vec)
assert len(pair_orders) == 45
assert all(o == 10 for o in pair_orders.values()), pair_orders
print("NEW CLAIM: all 45 pairwise [(a)-(b)] orders = exactly 10: OK")

# ---- orbits by brute-force permutations ----
def perm_vertex(p, v):
    return vidx[tuple(sorted(p[x] for x in verts[v]))]
perms = list(itertools.permutations(range(5)))
stab = [p for p in perms if perm_vertex(p, q) == q]
assert len(stab) == 12  # 2! * 3!
def orbits(items, group, act):
    seen, orbs = set(), []
    for it in items:
        if it in seen: continue
        orb = {act(p, it) for p in group}
        seen |= orb; orbs.append(orb)
    return orbs
vo = orbits(others, stab, perm_vertex)
assert sorted(len(o) for o in vo) == [3, 6]
disj = {v for v in others if not set(verts[v]) & set(verts[q])}
assert disj in vo and len(disj) == 3
print("stabiliser of q (order 12): vertex orbits off q are 3 (disjoint) "
      "+ 6 (meeting in one leaf): OK")
po = orbits(list(itertools.combinations(range(n), 2)), perms,
            lambda p, pr: tuple(sorted((perm_vertex(p, pr[0]),
                                        perm_vertex(p, pr[1])))))
assert sorted(len(o) for o in po) == [15, 30]
adj_pairs = {pr for pr in itertools.combinations(range(n), 2)
             if not set(verts[pr[0]]) & set(verts[pr[1]])}
assert adj_pairs in po and len(adj_pairs) == 15
print("S5 orbits on 45 pairs: 15 adjacent (disjoint) + 30 non-adjacent: OK")

# ---- parking functions from raw P2 ----
def outdeg(A, v):
    return sum(1 for u in adj[v] if u not in A)
def is_q_reduced(off):
    # off: dict vertex->value for v != q; P2 for every A subset of V\{q}
    for r in range(1, len(others) + 1):
        for A in itertools.combinations(others, r):
            S = set(A)
            if not any(off[v] < outdeg(S, v) for v in A):
                return False
    return True
count = 0
for vals in itertools.product(range(3), repeat=9):
    off = dict(zip(others, vals))
    if all(v >= 0 for v in vals) and is_q_reduced(off):
        count += 1
assert count == 2000, count
print("Route 3 G-parking functions (raw P2, all subsets): 2000: OK")

# ---- Try It answers ----
assert 15 - 10 + 1 == 6 and all(dd - 2 == 1 for dd in deg) and 10 == 2 * 6 - 2
assert 2 * 10**3 == 2000 and lcm(2, 10) == 10
assert 3**9 == 19683
print("TryIt 1,2,3: genus/K, SNF reading (drop five 1s), 3^9=19683 "
      "candidates -> 2000: OK")
# TryIt 4 counts already asserted in orbit section (3 disjoint, 6 meeting)

# ---- stabilisation exhibit ----
src = [4, -2, 3, 0, -1, 2, 0, 1, -3, 2]
tgt = [0, 2, 0, 1, 1, 1, 1, 0, 0, 0]
assert sum(src) == sum(tgt)
off = {v: tgt[v] for v in others}
assert tgt[q] == 0 and is_q_reduced(off)
diff = sp.Matrix(src) - sp.Matrix(tgt)
sol = L.LUsolve(diff) if False else None
# integer solution of L x = diff with x_q = 0 gauge: solve reduced system
rhs = sp.Matrix([int(diff[i]) for i in others])
xred = Lr.LUsolve(rhs)
assert all(v == int(v) for v in xred), xred
print("stabilisation exhibit: target q-reduced (raw P2), source-target "
      "in Laplacian lattice (integer potentials), degree preserved: OK")

print("\nALL CH.16 PHASE-2 INDEPENDENT CHECKS PASSED")
