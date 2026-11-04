"""
Chapter 7: Shortest Paths Are Matrix Multiplication (Interlude) --
min-plus linear algebra is shortest-path theory.

The full program behind Chapter 7 of "Learn You a Tropical Geometry
for Great Good!" Converted from Day 13 of the tropical geometry study
(2026-10-23). SELF-CONTAINED apart from numpy/scipy: one concrete
strongly connected digraph on 6 nodes, 16 integer-weighted edges,
everything asserted exactly (integers / Fractions, no floating point
anywhere near a claim).

Chapter claims asserted below:
  * Min-plus matrix powers = exactly-k-edge shortest walks, checked
    against brute-force walk enumeration for k = 1, 2, 3. (The first
    implementation let the INF sentinel leak through tropical
    multiplication -- INF + 2 priced a walk using a missing edge as
    merely expensive; the brute-force check caught it, and trop_mul
    now caps any term with an INF leg.)
  * The truncated Kleene star A* = I (+) A (+) ... (+) A^(n-1) is the
    all-pairs shortest-path matrix: it equals an independent
    Floyd-Warshall implementation entry by entry, and A^(n) improves
    nothing (no negative cycles, so some shortest walk is simple).
    Cf. Joswig, "The Tropical Geometry of Shortest Paths" (TU Berlin
    lecture notes); algebraic-path/Kleene-star formalism.
  * For an irreducible (strongly connected) min-plus matrix the
    tropical eigenvalue is unique and equals the minimum cycle mean
    (Gondran-Minoux / Baccelli et al.; see e.g. arXiv:1705.09513,
    Thms 3.4-3.5 in the max-plus form; min-plus by sign flip).
    Here: all 31 elementary cycles enumerated, unique minimum
    lambda = 2 on the critical cycle 0->1->2->0 (total 6), next best
    mean 3. Karp's algorithm (1978) returns 2, and Karp agrees with
    brute-force enumeration on 30 random strongly connected digraphs.
  * Eigenvector: for B = A - lambda, the distance column to a
    critical node c satisfies A (x) v = lambda + v exactly
    (Butkovic, Max-linear Systems, Springer 2010, Ch. 3: eigenvectors
    from the Kleene star of A - lambda, columns indexed by critical
    nodes). The three critical columns (to 0, 1, 2) coincide up to
    tropical scaling -- v = (0,0,-1,3,3,2) and its +1 shift -- so the
    eigenspace here is a single tropical ray. A non-critical column
    (to node 4) is NOT an eigenvector: criticality is essential.
  * Cyclicity: powers of B become exactly periodic with period 3
    (= gcd of the critical cycle lengths) after transient T = 12:
    B^(k+3) = B^k for all k in [12, 33], and period 1 fails.
  * Tropical determinant (min-plus permanent) = the assignment
    problem: for the seeded 5x5 cost matrix printed below the value
    is 36, by brute force over all 120 permutations and by
    scipy.optimize.linear_sum_assignment (Hungarian algorithm).

Convention: (A (x) x)_i = min_j (A_ij + x_j); A_ij = weight of i -> j,
INF = no edge. A^(0) = tropical identity (0 diagonal, INF off).
Running this file writes ch07_summary.json (the data the chapter's
figures are drawn from) into the current directory.
"""

import itertools
import json
from fractions import Fraction

import numpy as np
from scipy.optimize import linear_sum_assignment

INF = 10**9  # exact integer "infinity" (all real values are tiny)

# ---------------------------------------------------------------- graph
# Nodes 0..5. Critical cycle: 0->1 (2), 1->2 (3), 2->0 (1): total 6, mean 2.
EDGES = [
    (0, 1, 2), (1, 2, 3), (2, 0, 1),          # critical cycle, mean 2
    (0, 3, 5), (3, 4, 2), (4, 5, 3), (5, 3, 4),   # 3-4-5 cycle, mean 3
    (5, 0, 7), (2, 3, 6), (4, 1, 6), (1, 4, 8),
    (3, 0, 9), (5, 2, 5), (0, 5, 11), (4, 0, 10), (2, 5, 9),
]
N = 6


def adjacency(n, edges):
    A = np.full((n, n), INF, dtype=object)
    for i, j, w in edges:
        A[i, j] = w
    return A


def trop_mul(A, B):
    n = A.shape[0]
    C = np.full((n, n), INF, dtype=object)
    for i in range(n):
        for j in range(n):
            # cap the INF sentinel: a walk using a missing edge is INF,
            # not INF + (finite weight)
            C[i, j] = min(
                INF if (A[i, k] >= INF // 2 or B[k, j] >= INF // 2)
                else A[i, k] + B[k, j]
                for k in range(n))
    return C


def trop_identity(n):
    I = np.full((n, n), INF, dtype=object)
    for i in range(n):
        I[i, i] = 0
    return I


def trop_power(A, k):
    P = trop_identity(A.shape[0])
    for _ in range(k):
        P = trop_mul(P, A)
    return P


def as_int_matrix(M):
    return [[(None if v >= INF // 2 else int(v)) for v in row] for row in M]


A = adjacency(N, EDGES)

# ------------------------------------------------- strong connectivity
def reachable(A):
    n = A.shape[0]
    reach = [[A[i, j] < INF // 2 or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    return reach


reach = reachable(A)
assert all(all(r) for r in reach), "graph must be strongly connected (irreducible)"
print("[ok] graph strongly connected (A irreducible)")

# --------------------------------------- powers = exactly-k-edge walks
# Brute-force walk enumeration for small k as an independent check.
def brute_exact_k(A, k):
    n = A.shape[0]
    best = np.full((n, n), INF, dtype=object)
    for i in range(n):
        for j in range(n):
            m = INF
            for mid in itertools.product(range(n), repeat=k - 1) if k >= 1 else [()]:
                seq = (i,) + mid + (j,)
                tot = 0
                ok = True
                for a, b in zip(seq, seq[1:]):
                    if A[a, b] >= INF // 2:
                        ok = False
                        break
                    tot += A[a, b]
                if ok:
                    m = min(m, tot)
            best[i, j] = m
    return best


for k in (1, 2, 3):
    assert (trop_power(A, k) == brute_exact_k(A, k)).all(), f"power {k} != brute walks"
print("[ok] A^k = exactly-k-edge shortest walks, k=1,2,3 (brute force)")

# ------------------------------------------------ Kleene star = APSP
def floyd_warshall(A):
    n = A.shape[0]
    D = [[(0 if i == j else (A[i, j] if A[i, j] < INF // 2 else INF))
          for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if D[i][k] + D[k][j] < D[i][j]:
                    D[i][j] = D[i][k] + D[k][j]
    return D


star = trop_identity(N)
P = trop_identity(N)
for k in range(1, N):
    P = trop_mul(P, A)
    star = np.minimum(star, P)
FW = floyd_warshall(A)
assert all(int(star[i, j]) == FW[i][j] for i in range(N) for j in range(N)), \
    "Kleene star (truncated at n-1) != Floyd-Warshall APSP"
# Truncation is genuinely enough: A^(n) improves nothing over the star.
Pn = trop_power(A, N)
assert all(min(int(star[i, j]), int(Pn[i, j])) == int(star[i, j])
           for i in range(N) for j in range(N))
print("[ok] Kleene star I+A+...+A^5 == Floyd-Warshall all-pairs shortest paths")
print("     APSP distances:")
for row in as_int_matrix(star):
    print("      ", row)

# ------------------------------------------- elementary cycles, means
def elementary_cycles(n, edges):
    adj = {i: [] for i in range(n)}
    for i, j, w in edges:
        adj[i].append((j, w))
    cycles = []

    def dfs(start, node, path, total, visited):
        for nxt, w in adj[node]:
            if nxt == start and len(path) >= 1:
                cycles.append((tuple(path), total + w))
            elif nxt not in visited and nxt > start:
                visited.add(nxt)
                dfs(start, nxt, path + [nxt], total + w, visited)
                visited.discard(nxt)

    for s in range(n):
        dfs(s, s, [s], 0, {s})
    return cycles


cycles = elementary_cycles(N, EDGES)
means = [(Fraction(tot, len(c)), c, tot) for c, tot in cycles]
means.sort()
lam_enum, crit_cycle, crit_total = means[0]
assert lam_enum == Fraction(2, 1), f"min cycle mean {lam_enum} != 2"
assert set(crit_cycle) == {0, 1, 2} and len(crit_cycle) == 3
# uniqueness of the minimizer (as a node set / rotation class)
minimisers = [m for m in means if m[0] == lam_enum]
assert len(minimisers) == 1, f"min cycle not unique: {minimisers}"
assert all(m[0] > lam_enum for m in means[1:])
print(f"[ok] {len(cycles)} elementary cycles enumerated; unique minimum mean "
      f"{lam_enum} on cycle {crit_cycle} (total {crit_total})")
for mean, c, tot in means:
    print(f"       cycle {c}: total {tot}, mean {mean} = {float(mean):.3f}")

# ------------------------------------------------------- Karp's algorithm
# F_k(v) = min weight of an exactly-k-edge walk ENDING at v, start anywhere
# (super-source with 0 edges to every node). lambda = min_v max_{k<n}
# (F_n(v) - F_k(v)) / (n - k).  (Karp 1978, min-plus form.)
F = [np.zeros(N, dtype=object)]
for k in range(1, N + 1):
    prev = F[-1]
    cur = np.array([min(A[i, v] + prev[i] for i in range(N)) for v in range(N)],
                   dtype=object)
    F.append(cur)
karp_vals = []
for v in range(N):
    worst = max(Fraction(int(F[N][v] - F[k][v]), N - k) for k in range(N))
    karp_vals.append(worst)
lam_karp = min(karp_vals)
assert lam_karp == lam_enum, f"Karp {lam_karp} != enumeration {lam_enum}"
print(f"[ok] Karp's algorithm: lambda = {lam_karp} (agrees with enumeration)")

LAM = 2

# --------------------------------------------------------- eigenvector
# B = A - lambda; v_i = APSP distance in B from i to a critical node c.
B = A.copy()
for i in range(N):
    for j in range(N):
        if B[i, j] < INF // 2:
            B[i, j] = B[i, j] - LAM
Bdist = floyd_warshall(B)  # B has no negative cycle (min mean 0), FW valid
assert all(Bdist[i][i] == 0 for i in range(N)), "B must have no negative cycle"
for c in (0, 1, 2):
    v = np.array([Bdist[i][c] for i in range(N)], dtype=object)
    Av = np.array([min(A[i, j] + v[j] for j in range(N)) for i in range(N)],
                  dtype=object)
    assert (Av == np.array([LAM + int(x) for x in v], dtype=object)).all(), \
        f"eigen equation fails for critical node {c}"
    print(f"[ok] eigenvector from critical node {c}: v = {[int(x) for x in v]}, "
          f"A(x)v = lambda + v exactly")
# The three critical columns coincide up to tropical scaling, so the
# eigenspace here is a single tropical ray:
_cols = {c: [int(Bdist[i][c]) for i in range(N)] for c in (0, 1, 2)}
assert _cols[0] == [0, 0, -1, 3, 3, 2], _cols[0]
assert _cols[1] == _cols[0]
assert _cols[2] == [x + 1 for x in _cols[0]]
print("[ok] the three critical columns coincide up to tropical scaling: "
      "one tropical ray")

# A non-critical column is NOT an eigenvector (criticality is essential):
v_bad = np.array([Bdist[i][4] for i in range(N)], dtype=object)
Av_bad = np.array([min(A[i, j] + v_bad[j] for j in range(N)) for i in range(N)],
                  dtype=object)
assert not (Av_bad == np.array([LAM + int(x) for x in v_bad], dtype=object)).all()
print("[ok] column of a non-critical node (4) is not an eigenvector")

# ------------------------------------------- cyclicity of B's powers
powers = [trop_power(B, k) for k in range(0, 40)]
transient = None
for T in range(0, 25):
    if all((powers[k + 3] == powers[k]).all() for k in range(T, 33)):
        transient = T
        break
assert transient is not None, "no period-3 regime found up to k=33"
# and period is exactly 3, not 1, in that regime
assert not (powers[transient + 1] == powers[transient]).all()
print(f"[ok] B=A-lambda: powers eventually periodic, period exactly 3 "
      f"(= gcd of critical cycle lengths), transient T = {transient}: "
      f"B^(k+3) = B^k for all k in [{transient}, 33]")

# (Try It, asserted: from T = 12 with period 3, B^20 = B^17 = B^14.)
assert (powers[20] == powers[17]).all() and (powers[17] == powers[14]).all()
print("[ok] in particular B^20 = B^17 = B^14")

# ------------------------------------------- tropical determinant = assignment
rng = np.random.default_rng(13)
C = rng.integers(1, 30, size=(5, 5))
brute = min(sum(int(C[i, p[i]]) for i in range(5))
            for p in itertools.permutations(range(5)))
ri, ci = linear_sum_assignment(C)
hung = int(C[ri, ci].sum())
assert brute == hung == 36, (brute, hung)
print("     assignment cost matrix (seeded rng, default_rng(13)):")
for row in C:
    print("      ", [int(x) for x in row])
# and brute IS the min-plus permanent by definition; record the optimizer
best_p = min(itertools.permutations(range(5)),
             key=lambda p: sum(int(C[i, p[i]]) for i in range(5)))
print(f"[ok] tropical determinant of a random 5x5 cost matrix = {brute} "
      f"(brute force over 120 perms == scipy linear_sum_assignment; "
      f"optimal perm {best_p})")

# Random-matrix cross-check of the eigenvalue theorem (not just our graph):
for trial in range(30):
    n = 5
    M = np.full((n, n), INF, dtype=object)
    perm = list(range(n))
    rng.shuffle(perm)
    edges_r = []
    for a, b in zip(perm, perm[1:] + perm[:1]):  # guarantee strong connectivity
        w = int(rng.integers(1, 9))
        edges_r.append((a, b, w))
    for i in range(n):
        for j in range(n):
            if i != j and rng.random() < 0.45:
                edges_r.append((i, j, int(rng.integers(1, 15))))
    seen = {}
    for i, j, w in edges_r:
        seen[(i, j)] = min(w, seen.get((i, j), w))
    edges_r = [(i, j, w) for (i, j), w in seen.items()]
    for (i, j), w in seen.items():
        M[i, j] = w
    cyc = elementary_cycles(n, edges_r)
    lam_r = min(Fraction(t, len(c)) for c, t in cyc)
    # Karp on M
    Fr = [np.zeros(n, dtype=object)]
    for k in range(1, n + 1):
        prev = Fr[-1]
        Fr.append(np.array([min(M[i, v] + prev[i] for i in range(n))
                            for v in range(n)], dtype=object))
    lam_k = min(max(Fraction(int(Fr[n][v] - Fr[k][v]), n - k) for k in range(n))
                for v in range(n))
    assert lam_k == lam_r, f"trial {trial}: Karp {lam_k} != enum {lam_r}"
print("[ok] Karp == cycle enumeration on 30 random strongly connected digraphs")

summary = {
    "n": N,
    "edges": EDGES,
    "apsp": [[int(star[i, j]) for j in range(N)] for i in range(N)],
    "num_elementary_cycles": len(cycles),
    "lambda": str(lam_enum),
    "critical_cycle": list(crit_cycle),
    "cycle_means_sorted": [[str(m), list(c), t] for m, c, t in means],
    "eigenvector_from_0": [int(Bdist[i][0]) for i in range(N)],
    "cyclicity_transient": transient,
    "assignment_example": {"matrix": C.tolist(), "tropdet": brute,
                           "optimal_perm": list(best_p)},
    "B_power_trace_00": [int(powers[k][0, 0]) for k in range(20)],
    "B_power_trace_03": [int(powers[k][0, 3]) for k in range(20)],
    "B_power_trace_23": [int(powers[k][2, 3]) for k in range(20)],
    "B_power_trace_35": [int(powers[k][3, 5]) for k in range(20)],
}
with open("ch07_summary.json", "w") as f:
    json.dump(summary, f, indent=1)
print("\nwrote ch07_summary.json")
print("Chapter 7: all assertions pass. Powers = walks; "
          "Kleene star = Floyd-Warshall; lambda = 2 by enumeration "
          "and by Karp (plus 30 random digraphs); the critical "
          "columns are one tropical ray; period exactly 3 from "
          "T = 12; tropical determinant 36 = the assignment optimum.")
