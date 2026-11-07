#!/usr/bin/env python3
"""
Chapter 9: Hodge Theory With No Hodge Theory Underneath -- the Chow
ring A*(U(3,4)) and the Hodge-Riemann relations, verified exactly.

The full program behind Chapter 9 of "Learn You a Tropical Geometry
for Great Good!" Converted from Day 7 of the tropical geometry study
(2026-10-17: day07_chow_ring.py + day07b_hodge_riemann.py, merged
into ONE self-contained program in the Chapter 9 Phase 2 verification
session, 2026-11-06). SELF-CONTAINED apart from numpy (eigenvalues,
display only; every asserted quantity is an exact Fraction).

Presentation (Feichtner-Yuzvinsky, via M. Baker, "Hodge Theory in
Combinatorics", arXiv:1705.07960, sec 4.1):
  A*(M) = Z[x_F : F nonempty proper flat] / (CH1 + CH2)
  (CH1) for all a, b in E: sum_{F ni a} x_F = sum_{F ni b} x_F
  (CH2) x_F * x_G = 0 whenever F, G are incomparable.
deg: A^2 -> Z is pinned by deg = 1 on maximal-chain monomials; its
uniqueness (up to scale) is the tropical reading that deg is the
unique integer weighting making the Bergman fan balanced (Baker
sec 4.1).  Classes alpha = sum_{F ni 0} x_F, beta = sum_{F not ni 0}
x_F (sec 4.4); m_k = deg(alpha^{2-k} beta^k) are the unsigned
coefficients of the reduced characteristic polynomial (Thm 4.2).
Ample classes come from strictly submodular functions (sec 4.3);
Hodge-Riemann is Thm 4.1; Bjorner's flag count is Prop 4.3.

Subdivision note (settled in Phase 2, answering the Ch. 8 referee
thread): this presentation is built from the lattice of flats
ALONE -- no fan and no subdivision is chosen to define it. When it
is read as the Chow ring OF the Bergman fan, it is the fan's FINE
subdivision (Ch. 8, description 2) whose rays are the x_F: one per
nonempty proper flat (10 for U(3,4)). The coarsest fan structure on
the same set has only 4 rays and 6 cones and carries no such
generators. Definition: flats. Fan reading: fine.

Chapter claims asserted below, no exceptions:
  * 10 generators, 3 independent (CH1) relations, graded dims
    (1, 7, 1); deg constant (= -1 pre-normalization, = 1 after) on
    all 12 maximal-chain monomials; Poincare Gram on A^1 has
    rank 7 (nondegenerate).
  * m = (1, 3, 3) THREE ways: ring intersection numbers; Moebius /
    chi(t) = t^3 - 4t^2 + 6t - 3, chi-bar(t) = t^2 - 3t + 3;
    Bjorner's initial-descending flag count. Log-concavity 3^2 >= 1*3.
  * Ample class: c(S) = |S|(4-|S|) strictly submodular on all 110
    ORDERED incomparable subset pairs (55 unordered); the failed
    candidate c(F) = rk(F)(3-rk(F)) gives EQUALITY (2 = 2) at
    A = {0,1}, B = {0,2,3} -- the Honest Box's evidence, reproduced.
  * Hodge-Riemann: deg(l^2) = 120 exactly (> 0, the k = 0 case);
    dim P^1 = 6; -deg(.,.) on P^1 has Sylvester leading minors
    exactly [97/18, 260/9, 925/6, 435, 2133/2, 2430], all > 0, and
    (float, display) eigenvalues [1, 1, 10/3, 9, 9, 9].
  * Try It key: (1) 3 independent linear relations in degree 1;
    (2) a degree-2 monomial survives iff its two flats are
    comparable -- maximal ones are the 12 chain monomials;
    (3) U(3,5)'s (1, 4, 6) must satisfy the SINGLE log-concavity
    inequality 4^2 >= 1*6 (a 3-term sequence has exactly one).

Running this file writes ch09_summary.json into the current
directory (the data the chapter's figures are drawn from).
"""

import json
from fractions import Fraction
from itertools import combinations

import numpy as np

E = (0, 1, 2, 3)
N = len(E)
RANK = 3

# Flats of U(3,4): every subset with |S| < 3 is a flat; |S| >= 3 closes to E.
flats = [frozenset([i]) for i in E] + [frozenset(c) for c in combinations(E, 2)]
NF = len(flats)
assert NF == 10, NF
fidx = {f: i for i, f in enumerate(flats)}


def comparable(f, g):
    return f <= g or g <= f


def mat_rank(rows):
    m = [list(r) for r in rows]
    r = 0
    for c in range(len(m[0])):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pivv = m[r][c]
        m[r] = [v / pivv for v in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
    return r


def nullspace(rows, ncols):
    m = [list(r) for r in rows]
    pivots = {}
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pivv = m[r][c]
        m[r] = [v / pivv for v in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        pivots[c] = r
        r += 1
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * ncols
        v[fc] = 1
        for c, ri in pivots.items():
            v[c] = -m[ri][fc]
        basis.append(v)
    return basis


def rref_pivots(rows, ncols):
    m = [list(r) for r in rows]
    pivots = {}
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pivv = m[r][c]
        m[r] = [v / pivv for v in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        pivots[c] = r
        r += 1
    return m, pivots


# ---------- A^1: generators mod (CH1) ----------
rel1 = []
for a in (1, 2, 3):
    row = [Fraction(0)] * NF
    for f, i in fidx.items():
        if a in f:
            row[i] += 1
        if 0 in f:
            row[i] -= 1
    rel1.append(row)

rank_rel1 = mat_rank(rel1)
assert rank_rel1 == 3, rank_rel1  # Try It 1: 3 independent relations
dimA1 = NF - rank_rel1
assert dimA1 == 7, dimA1
print(f"A^1: 10 generators - 3 independent (CH1) relations = dim {dimA1}")

# ---------- A^2: monomials mod (CH2) + lifted (CH1) ----------
mons = [(i, j) for i in range(NF) for j in range(i, NF)]
midx = {m: k for k, m in enumerate(mons)}
NM = len(mons)

rel2 = []
for i, j in mons:
    if not comparable(flats[i], flats[j]):
        row = [Fraction(0)] * NM
        row[midx[(i, j)]] = 1
        rel2.append(row)
for a in (1, 2, 3):
    lin = [Fraction(0)] * NF
    for f, i in fidx.items():
        if a in f:
            lin[i] += 1
        if 0 in f:
            lin[i] -= 1
    for k in range(NF):
        row = [Fraction(0)] * NM
        for i, c in enumerate(lin):
            if c == 0:
                continue
            m = (i, k) if i <= k else (k, i)
            row[midx[m]] += c
        if any(v != 0 for v in row):
            rel2.append(row)

dimA2 = NM - mat_rank(rel2)
assert dimA2 == 1, dimA2
print(f"A^2: dim {dimA2}; graded dims (1, 7, 1) with A^0 = Z")

# Try It 2: surviving degree-2 monomials are exactly the comparable pairs.
surviving = [(i, j) for (i, j) in mons if comparable(flats[i], flats[j])]
killed = [(i, j) for (i, j) in mons if not comparable(flats[i], flats[j])]
assert len(surviving) + len(killed) == NM
print(f"Try It 2: {len(surviving)} comparable monomials survive (CH2); "
      f"{len(killed)} incomparable ones vanish")

# ---------- degree functional ----------
ns = nullspace(rel2, NM)
assert len(ns) == 1, len(ns)
d_raw = ns[0]
chains = []
for i in E:
    for j, k in combinations(E, 2):
        if i in (j, k):
            a, b = fidx[frozenset([i])], fidx[frozenset([j, k])]
            m = (a, b) if a <= b else (b, a)
            chains.append(midx[m])
assert len(chains) == 12, len(chains)
vals = {d_raw[c] for c in chains}
assert vals == {Fraction(-1)}, vals  # constant -1 pre-normalization
scale = next(iter(vals))
d = [v / scale for v in d_raw]
assert all(d[c] == 1 for c in chains)
print("deg: constant -1 on all 12 maximal-chain monomials pre-normalization; = 1 after")


def deg_of_poly(poly):
    return sum(c * d[midx[m]] for m, c in poly.items())


def lmul(p1, p2):
    out = {}
    for i, c1 in p1.items():
        for j, c2 in p2.items():
            m = (i, j) if i <= j else (j, i)
            out[m] = out.get(m, Fraction(0)) + c1 * c2
    return out


# ---------- route 1: intersection numbers in the ring ----------
la = {fidx[f]: Fraction(1) for f in flats if 0 in f}
lb = {fidx[f]: Fraction(1) for f in flats if 0 not in f}
mu = []
for k in range(3):
    factors = [la] * (2 - k) + [lb] * k
    mu.append(deg_of_poly(lmul(factors[0], factors[1])))
assert mu == [Fraction(1), Fraction(3), Fraction(3)], mu
print(f"route 1 (ring): deg(alpha^2, alpha*beta, beta^2) = {[str(x) for x in mu]}")

# ---------- route 2: Moebius / characteristic polynomial ----------
all_flats = [frozenset()] + flats + [frozenset(E)]
mu_mob = {frozenset(): Fraction(1)}
for f in flats + [frozenset(E)]:
    mu_mob[f] = -sum(mu_mob[g] for g in all_flats if g < f)
chi = {}
for f in all_flats:
    rk = 3 if f == frozenset(E) else len(f)
    chi[3 - rk] = chi.get(3 - rk, Fraction(0)) + mu_mob[f]
a3, a2, a1, a0 = chi[3], chi[2], chi[1], chi[0]
assert (a3, a2, a1, a0) == (1, -4, 6, -3), (a3, a2, a1, a0)
b2, b1, b0 = a3, a2 + a3, a1 + a2 + a3
assert a0 + a1 + a2 + a3 == 0
assert (b2, b1, b0) == (1, -3, 3), (b2, b1, b0)
assert [abs(b2), abs(b1), abs(b0)] == mu
print("route 2 (Moebius): chi = t^3 - 4t^2 + 6t - 3; chi-bar = t^2 - 3t + 3 -> (1, 3, 3)")

# ---------- Poincare duality: pairing on A^1 nondegenerate ----------
m1, piv1 = rref_pivots(rel1, NF)
free1 = [c for c in range(NF) if c not in piv1]
assert len(free1) == 7
basis1 = []
for fc in free1:
    v = [Fraction(0)] * NF
    v[fc] = 1
    for c, ri in piv1.items():
        v[c] = -m1[ri][fc]
    basis1.append(v)


def degprod(u, v):
    return deg_of_poly(lmul({i: c for i, c in enumerate(u) if c},
                            {j: c for j, c in enumerate(v) if c}))


G_pd = [[degprod(u, v) for v in basis1] for u in basis1]
assert mat_rank(G_pd) == 7
print("Poincare pairing on A^1: Gram rank 7 (nondegenerate)")

# ---------- ample class: strict submodularity, brute force ----------
all_subsets = [frozenset(s) for r_ in range(N + 1) for s in combinations(E, r_)]


def cfun(S):
    return len(S) * (N - len(S))


assert cfun(frozenset()) == 0 and cfun(frozenset(E)) == 0
assert all(cfun(S) >= 0 for S in all_subsets)
n_ordered = 0
unordered_pairs = set()
for A in all_subsets:
    for B in all_subsets:
        if A <= B or B <= A:
            continue
        assert cfun(A | B) + cfun(A & B) < cfun(A) + cfun(B), (A, B)
        n_ordered += 1
        unordered_pairs.add(frozenset((A, B)))
assert n_ordered == 110, n_ordered          # ordered incomparable pairs
assert len(unordered_pairs) == 55, len(unordered_pairs)
print("c(S) = |S|(4-|S|): strictly submodular on all 110 ordered "
      "incomparable pairs (55 unordered)")

# Honest Box evidence: the rank-based candidate FAILS, at the stated pair.
def rk_of(S):
    return min(len(S), RANK)


def c_bad(S):
    return rk_of(S) * (RANK - rk_of(S))


A_bad, B_bad = frozenset({0, 1}), frozenset({0, 2, 3})
lhs_bad = c_bad(A_bad | B_bad) + c_bad(A_bad & B_bad)
rhs_bad = c_bad(A_bad) + c_bad(B_bad)
assert lhs_bad == rhs_bad == 2, (lhs_bad, rhs_bad)  # equality, NOT strict
print("Honest Box: c(F) = rk(F)(3-rk(F)) gives equality 2 = 2 at "
      "A = {0,1}, B = {0,2,3} -- not strictly submodular")

ell = [Fraction(cfun(f)) for f in flats]

# ---------- Hodge-Riemann k = 0 ----------
q0 = degprod(ell, ell)
assert q0 == 120, q0
print(f"HRR k=0: deg(l^2) = {q0} > 0 (exact)")

# ---------- Hodge-Riemann k = 1: -deg(.,.) positive definite on P^1 ----------
row = [[degprod(ell, b) for b in basis1]]
P1 = nullspace(row, 7)
assert len(P1) == 6, len(P1)


def lincomb(coeffs):
    v = [Fraction(0)] * NF
    for c, b in zip(coeffs, basis1):
        for i in range(NF):
            v[i] += c * b[i]
    return v


P1vecs = [lincomb(p) for p in P1]
G = [[-degprod(u, v) for v in P1vecs] for u in P1vecs]
assert all(G[i][j] == G[j][i] for i in range(6) for j in range(6))


def det(M):
    n = len(M)
    Mx = [list(r) for r in M]
    sgn = 1
    for c in range(n):
        piv = next((i for i in range(c, n) if Mx[i][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            Mx[c], Mx[piv] = Mx[piv], Mx[c]
            sgn = -sgn
        for i in range(c + 1, n):
            f = Mx[i][c] / Mx[c][c]
            for j in range(c, n):
                Mx[i][j] -= f * Mx[c][j]
    p = Fraction(1)
    for i in range(n):
        p *= Mx[i][i]
    return sgn * p


minors = [det([r_[:k] for r_ in G[:k]]) for k in range(1, 7)]
expected_minors = [Fraction(97, 18), Fraction(260, 9), Fraction(925, 6),
                   Fraction(435), Fraction(2133, 2), Fraction(2430)]
assert minors == expected_minors, minors
assert all(m > 0 for m in minors)
print(f"HRR k=1: Sylvester minors {[str(m) for m in minors]} -- all > 0, exact")

eigs = sorted(float(v) for v in np.linalg.eigvalsh(
    np.array([[float(x) for x in r_] for r_ in G])))
expected_eigs = [1.0, 1.0, 10.0 / 3.0, 9.0, 9.0, 9.0]
assert all(abs(e - x) < 1e-9 for e, x in zip(eigs, expected_eigs)), eigs
print(f"eigenvalues (float, display): {[round(e, 4) for e in eigs]} = [1, 1, 10/3, 9, 9, 9]")

# ---------- route 3: Bjorner flag count (Prop 4.3) ----------
def initial_descending(k):
    count = 0
    for chain in combinations(flats, k):
        if not all(chain[i] < chain[i + 1] for i in range(k - 1)):
            continue
        if not all(len(chain[i]) == i + 1 for i in range(k)):
            continue
        mins = [min(f) for f in chain]
        if all(mins[i] > mins[i + 1] for i in range(k - 1)) and (k == 0 or mins[-1] > 0):
            count += 1
    return count


mk = [initial_descending(k) for k in range(3)]
assert mk == [1, 3, 3], mk
print(f"route 3 (Bjorner): flag counts {mk}")

# ---------- log-concavity ----------
m0, m1_, m2 = mk
assert [Fraction(x) for x in mk] == mu
assert m1_ ** 2 >= m0 * m2 and m1_ ** 2 == 9 and m0 * m2 == 3
# Try It 3 (corrected): a 3-term sequence has exactly ONE such inequality.
assert 4 ** 2 >= 1 * 6
print("log-concavity: 3^2 = 9 >= 1*3 = 3; Try It 3 key: (1,4,6) needs "
      "the single inequality 4^2 = 16 >= 1*6 = 6 -- and it holds")

summary = {
    "matroid": "U(3,4)",
    "generators": NF,
    "independent_linear_relations": rank_rel1,
    "graded_dims": [1, dimA1, dimA2],
    "maximal_chains": len(chains),
    "deg_on_chains_prenormalization": str(scale),
    "m_ring": [str(x) for x in mu],
    "chi_coeffs_t3_to_t0": [str(x) for x in (a3, a2, a1, a0)],
    "chi_bar_coeffs": [str(x) for x in (b2, b1, b0)],
    "m_bjorner": mk,
    "submodular_ordered_incomparable_pairs": n_ordered,
    "submodular_unordered_incomparable_pairs": len(unordered_pairs),
    "failed_candidate_equality": {"A": [0, 1], "B": [0, 2, 3],
                                  "lhs": int(lhs_bad), "rhs": int(rhs_bad)},
    "deg_l_squared": str(q0),
    "dim_P1": len(P1),
    "sylvester_minors": [str(m) for m in minors],
    "eigenvalues_float": eigs,
    "log_concavity": {"m1_squared": m1_ ** 2, "m0_times_m2": m0 * m2},
}
with open("ch09_summary.json", "w") as f:
    json.dump(summary, f, indent=1)
print("\nALL CHECKS PASSED: A*(U(3,4)) dims (1,7,1); m = (1,3,3) three ways;")
print("deg(l^2) = 120; -deg(.,.) positive definite on the 6-dim P^1, exactly.")
print("The Kahler package for A*(U(3,4)) is real, verified exactly.")
