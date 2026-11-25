#!/usr/bin/env python3
"""Chapter 15 (Deep Dive 2) -- The Matroid With No Matrices. SITE BUNDLE program.

Merged self-contained from Day 16 (2026-10-26) plus the Phase 2
verification groups (2026-11-24). The Vamos matroid V8: the Kahler
package where no vector model exists (AHK is field-free).

V8 is the rank-4 matroid on E = {0..7} with blocks A={0,1}, B={2,3},
C={4,5}, D={6,7} whose bases are all 4-subsets EXCEPT the five block
unions H = {A+B, A+C, A+D, B+C, B+D}; rho(Y) = 3 if Y in H else
min(|Y|, 4). Definition: "The Vamos Matroid Has No Second Symmetric
Power" Def. 3.1, agreeing with Oxley Ex. 2.1.25 (verified live Day 16).
Non-representability is proved in code: submodularity holds on all
256^2 pairs (it IS a matroid) but Ingleton fails at (A,B,C,D):
16 <= 15 is false, so no matrix over ANY field has this rank function.

Headline (all asserted): flats 79 (77 proper nonempty);
m = (1,7,21,30) by Mobius/chi-bar, Bjorner flags, and FY intersection
theory; U(4,8) has (1,7,21,35), so five circuit-hyperplanes cost
exactly (0,0,0,5); dims A*(V8) = (1,70,70,1); Poincare pairing and
Hard Lefschetz l: A^1 -> A^2 both rank 70; deg(l^3) = 150720 exactly;
Hodge-Riemann positive definite on P^1 (dim 69) by exact LDL^T over Q
(pivots = ratios of successive Sylvester minors).

Presentation: Feichtner-Yuzvinsky A*(M) = Z[x_F]/(CH1+CH2), Baker
"Hodge Theory in Combinatorics" arXiv:1705.07960 sec 4.1/4.3/4.4,
Thm 4.1/4.2, Prop 4.3, surveying Adiprasito-Huh-Katz. The sparse
chain-monomial Chow ring below is Day 14's machinery, inlined from
the shipped Chapter 14 bundle (a monomial is zero unless its flats
form a chain, by CH2). All arithmetic exact (Fractions); eigenvalues
float, for the figure only. Writes ch15_summary.json next to itself.
"""
import itertools
import json
import os
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = {}
t_start = time.time()

# ---------- sparse exact linear algebra (dict rows of Fractions) ----------
def gauss_jordan(rows, ncols):
    """Fully reduced row echelon form of sparse rows. Returns
    (pivot_rows, pivots): pivot_rows[pivot_col] = row dict with 1 at the
    pivot column and support only on free columns otherwise."""
    rows = [dict(r) for r in rows if r]
    used = [False] * len(rows)
    pivot_idx = {}
    for col in range(ncols):
        cand = [i for i in range(len(rows))
                if not used[i] and col in rows[i]]
        if not cand:
            continue
        i = min(cand, key=lambda j: len(rows[j]))
        used[i] = True
        pv = rows[i][col]
        rows[i] = {c: v / pv for c, v in rows[i].items()}
        for j in range(len(rows)):
            if j != i and col in rows[j]:
                f = rows[j][col]
                new = dict(rows[j])
                for c, v in rows[i].items():
                    nv = new.get(c, Fraction(0)) - f * v
                    if nv == 0:
                        new.pop(c, None)
                    else:
                        new[c] = nv
                rows[j] = new
        pivot_idx[col] = i
    # NOTE: read the pivot rows only now -- elimination rebinds rows[j] to
    # fresh dicts, so references captured mid-loop go stale (day-14 bug).
    return {col: rows[i] for col, i in pivot_idx.items()}

def reduce_vec(vec, pivot_rows):
    """Reduce a monomial vector (dict col->coeff) modulo the relations:
    returns dict supported on free (quotient-basis) columns only."""
    v = dict(vec)
    while True:
        todo = [c for c in v if c in pivot_rows]
        if not todo:
            return v
        for col in todo:
            if col not in v:
                continue
            f = v.pop(col)
            for c, coeff in pivot_rows[col].items():
                if c == col:
                    continue
                # pivot row: x_col + sum coeff*x_c = 0  =>  x_col = -sum
                nv = v.get(c, Fraction(0)) - f * coeff
                if nv == 0:
                    v.pop(c, None)
                else:
                    v[c] = nv

def det_bareiss(M):
    n = len(M); A = [list(map(int, row)) for row in M]; sgn = 1; prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            assert piv is not None, "singular in Bareiss"
            A[k], A[piv] = A[piv], A[k]; sgn = -sgn
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            A[i][k] = 0
        prev = A[k][k]
    return sgn * A[n - 1][n - 1]

def dense_rank(mat):
    """Rank of a dense Fraction matrix (small)."""
    m = [row[:] for row in mat]; nr = len(m); nc = len(m[0]); r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if m[i][c] != 0), None)
        if piv is None: continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]; m[r] = [v / pv for v in m[r]]
        for i in range(nr):
            if i != r and m[i][c] != 0:
                f = m[i][c]; m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
    return r

# ---------- Chow ring of a matroid given by its flats ----------
class ChowRing:
    def __init__(self, E, flats_all, rank):
        """flats_all: all flats incl. empty and E. rank = matroid rank."""
        self.E = list(E); self.rank = rank; self.d = rank - 1
        self.gens = [f for f in flats_all if f and f != frozenset(E)]
        self.gidx = {f: i for i, f in enumerate(self.gens)}
        self.NF = len(self.gens)
        # CH1 linear forms L_a - L_0 as dicts var->coeff
        self.linforms = []
        for a in self.E[1:]:
            lf = {}
            for f, i in self.gidx.items():
                c = (1 if a in f else 0) - (1 if self.E[0] in f else 0)
                if c: lf[i] = Fraction(c)
            self.linforms.append(lf)
        self._build_degrees()

    def _chain_monomials(self, k):
        """All multisets of k generators forming a chain (sorted tuples)."""
        if k == 0: return [()]
        out = []
        def rec(start_min_size, last_flat, cur):
            if len(cur) == k:
                out.append(tuple(cur)); return
            for i, f in enumerate(self.gens):
                if len(f) < start_min_size: continue
                if last_flat is not None and not last_flat <= f: continue
                cur.append(i); rec(len(f), f, cur); cur.pop()
        rec(0, None, [])
        return out

    def _insert(self, mon, var):
        """Monomial (tuple) times x_var, or None if not a chain (CH2 = 0)."""
        fv = self.gens[var]
        for i in mon:
            f = self.gens[i]
            if not (f <= fv or fv <= f): return None
        new = list(mon) + [var]
        new.sort(key=lambda i: (len(self.gens[i]), i))
        return tuple(new)

    def _build_degrees(self):
        self.mons = {}; self.midx = {}; self.pivots = {}; self.basis = {}
        for k in range(self.d + 1):
            mons = self._chain_monomials(k)
            self.mons[k] = mons
            self.midx[k] = {m: i for i, m in enumerate(mons)}
            if k == 0:
                self.pivots[k] = {}; self.basis[k] = [0]; continue
            rows = []
            for lf in self.linforms:
                for pm in self.mons[k - 1]:
                    row = {}
                    for var, c in lf.items():
                        nm = self._insert(pm, var)
                        if nm is not None:
                            j = self.midx[k][nm]
                            row[j] = row.get(j, Fraction(0)) + c
                    rows.append({j: v for j, v in row.items() if v})
            pr = gauss_jordan(rows, len(mons))
            self.pivots[k] = pr
            self.basis[k] = [j for j in range(len(mons)) if j not in pr]
        # degree functional: A^d is 1-dim; normalize on maximal flags
        d = self.d
        assert len(self.basis[d]) == 1, f"A^{d} dim {len(self.basis[d])} != 1"
        self.deg_scale = None
        flags = [m for m in self.mons[d] if len(set(m)) == d]
        coords = set()
        for m in flags:
            red = reduce_vec({self.midx[d][m]: Fraction(1)}, self.pivots[d])
            assert len(red) == 1
            coords.add(next(iter(red.values())))
        assert len(coords) == 1 and coords != {Fraction(0)}, \
            f"flag monomials disagree: {coords}"
        self.deg_scale = 1 / coords.pop()   # deg(basis monomial of A^d)

    # quotient element = dict basis_col -> Fraction, per degree
    def elem_from_poly(self, k, poly):
        """poly: dict monomial-tuple -> coeff (in the monomial basis)."""
        vec = {self.midx[k][m]: c for m, c in poly.items() if c}
        red = reduce_vec(vec, self.pivots[k])
        return {self.basis[k].index(j): v for j, v in red.items()}

    def gen_elem(self, var):
        return self.elem_from_poly(1, {(var,): Fraction(1)})

    def mul_gen(self, k, elem, var):
        """Multiply an A^k element (basis coords) by generator x_var."""
        out = {}
        for bi, c in elem.items():
            mon = self.mons[k][self.basis[k][bi]]
            nm = self._insert(mon, var)
            if nm is None: continue
            red = reduce_vec({self.midx[k + 1][nm]: c}, self.pivots[k + 1])
            for j, v in red.items():
                bj = self.basis[k + 1].index(j)
                out[bj] = out.get(bj, Fraction(0)) + v
        return {b: v for b, v in out.items() if v}

    def mul(self, k, elem, elem1):
        """Multiply A^k element by an A^1 element (given by its lift:
        dict generator var -> coeff)."""
        out = {}
        for var, cv in elem1.items():
            part = self.mul_gen(k, elem, var)
            for b, v in part.items():
                out[b] = out.get(b, Fraction(0)) + cv * v
        return {b: v for b, v in out.items() if v}

    def lift1(self, elem1_coords):
        """A^1 element in basis coords -> representative poly over gens."""
        return {self.basis[1][bi]: c for bi, c in elem1_coords.items()}

    def deg(self, elem_top):
        assert len(elem_top) <= 1
        if not elem_top: return Fraction(0)
        (bi, v), = elem_top.items()
        assert bi == 0
        return v * self.deg_scale

    def op_matrix(self, k, elem1_lift):
        """Dense matrix of multiplication by an A^1 element A^k -> A^{k+1}
        (columns = images of basis elements), as list of rows."""
        nk, nk1 = len(self.basis[k]), len(self.basis[k + 1])
        M = [[Fraction(0)] * nk for _ in range(nk1)]
        for bi in range(nk):
            img = self.mul(k, {bi: Fraction(1)}, elem1_lift)
            for bj, v in img.items(): M[bj][bi] = v
        return M

def mobius_chibar(E, flats_all, rank_of):
    """chi(t) via Mobius over the flat lattice; chi-bar = chi/(t-1);
    returns unsigned coeffs m_k = |coeff of t^{d-k}|, d = rank-1."""
    flats = sorted(flats_all, key=lambda f: (len(f), sorted(f)))
    mu = {}
    for f in flats:
        s = Fraction(0)
        for g in flats:
            if g < f: s += mu[g]
        mu[f] = -s if f else Fraction(1)
    # chi(t) = sum_F mu(F) t^{rank - rk(F)}: coeffs low->high
    r = max(rank_of(f) for f in flats)
    coeff = [Fraction(0)] * (r + 1)
    for f in flats:
        coeff[r - rank_of(f)] += mu[f]
    # divide by (t-1): synthetic on high->low coefficients
    high = coeff[::-1]
    assert high[0] == 1
    quot = []
    carry = Fraction(0)
    for c in high:
        carry = c + carry
        quot.append(carry)
    assert quot.pop() == 0, "t=1 not a root of chi"
    m = [abs(c) for c in quot]           # quot[k] ~ t^{d-k}
    assert all(v.denominator == 1 for v in m)
    return [int(v) for v in m], mu

def uniform(E, r):
    E = list(E); Es = frozenset(E)
    flats = [frozenset(c) for s in range(r) for c in combinations(E, s)]
    flats.append(Es)
    return flats, (lambda f: min(len(f), r))

def analyse(name, E, flats, rank_of, expect_dims, expect_m):
    print(f"== {name} ==")
    CR = ChowRing(E, flats, max(rank_of(f) for f in flats))
    dims = tuple(len(CR.basis[k]) for k in range(CR.d + 1))
    print(f"  generators NF={CR.NF}; chain monomials per degree: "
          f"{[len(CR.mons[k]) for k in range(CR.d + 1)]}")
    print(f"  dims A* = {dims}")
    assert dims == expect_dims, (dims, expect_dims)
    assert dims == dims[::-1], "Poincare palindrome fails"
    # alpha, beta
    a_lift, b_lift = {}, {}
    for f, i in CR.gidx.items():
        (a_lift if E[0] in f else b_lift)[i] = Fraction(1)
    alpha = CR.elem_from_poly(1, {(i,): c for i, c in a_lift.items()})
    beta = CR.elem_from_poly(1, {(i,): c for i, c in b_lift.items()})
    # m_k by FY intersection theory
    m_fy = []
    for k in range(CR.d + 1):
        el = {0: Fraction(1)} if CR.d - k == 0 else None
        cur = {0: Fraction(1)}  # in A^0
        deg_done = 0
        seq = [a_lift] * (CR.d - k) + [b_lift] * k
        for lf in seq:
            cur = CR.mul(deg_done, cur, lf); deg_done += 1
        m_fy.append(CR.deg(cur))
    assert all(v.denominator == 1 for v in m_fy)
    m_fy = [int(v) for v in m_fy]
    m_mb, mu = mobius_chibar(E, flats, rank_of)
    print(f"  m (FY intersection) = {m_fy}")
    print(f"  m (Mobius/chi-bar)  = {m_mb}")
    assert m_fy == m_mb == expect_m
    for k in range(1, CR.d):
        assert m_fy[k] ** 2 >= m_fy[k - 1] * m_fy[k + 1], "log-concavity fails"
    print(f"  log-concavity: {[f'{m_fy[k]}^2>={m_fy[k-1]*m_fy[k+1]}' for k in range(1, CR.d)]}")
    return CR, alpha, beta, a_lift, b_lift, m_fy

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
CR, _a, _b, a_lift, b_lift, m = analyse(
    "Vamos V8 (bundle Chow ring)", E, flats, lambda f: rk_of[f],
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
assert deg_l3 == 150720
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

with open(os.path.join(HERE, "ch15_summary.json"), "w") as fh:
    json.dump(SUMMARY, fh, indent=1)
print(f"\nDAY-16 CORE ASSERTIONS PASSED in {time.time() - t_start:.1f}s")


# ---------- Phase 2 verification groups (merged at ship, 2026-11-25) ----------
def verify_phase2():
    """Independent checks from hidden_files/code/ch15_phase2_verify.py
    (drafts/ch15_phase2_verify.py on GitHub): Try It answers asserted
    individually, full Whitney/chi/chi-bar for both matroids, the
    per-defect Mobius proof of sec 15.3 (with referee item 6's sign:
    w_{r-1} moves by (-1)^r per defect -- +1 here since r=4, while
    |m_{r-1}| drops by one in general), sec 15.6 counts, and the
    ordered-pairs count of sec 15.5."""
    print("== Phase 2 verification groups ==")
    # Try It 1
    assert A | C in Hset and rho(A | C) == 3
    assert C | D not in Hset and rho(C | D) == 4
    # Try It 2
    rhs_sets = [A | B, A | C, A | D, B | C, B | D]
    assert rhs_sets == H and all(rho(s) == 3 for s in rhs_sets)
    assert sum(rho(s) for s in rhs_sets) == 15
    # Whitney numbers both matroids
    def whitney(flats_, rk_):
        m_, mu_ = mobius_chibar(E, flats_, rk_)
        w = [0] * 5
        for f in flats_:
            w[rk_(f)] += int(mu_[f])
        return w, m_, mu_
    w_v, m_v, mu_v = whitney(flats, lambda f: rk_of[f])
    assert w_v == [1, -8, 28, -51, 30] and m_v == [1, 7, 21, 30]
    uflats = [frozenset(c) for s in range(4)
              for c in itertools.combinations(E, s)] + [frozenset(E)]
    urk = lambda f: min(len(f), 4)
    w_u, m_u, mu_u = whitney(uflats, urk)
    assert w_u == [1, -8, 28, -56, 35] and m_u == [1, 7, 21, 35]
    cum, s = [], 0
    for x in w_v[:-1]:
        s += x; cum.append(abs(s))
    assert cum == m_v
    # Try It 3 + log-concavity with true products
    assert [u - v for u, v in zip(m_u, m_v)] == [0, 0, 0, 5]
    assert 7 ** 2 >= 1 * 21 and 21 ** 2 == 441 >= 7 * 30 == 210
    # Per-defect proof (sec 15.3): each H absorbs 4 triple flats
    flatset = set(flats)
    def closure(S):
        S = set(S); r = rk_of[frozenset(S)]
        return frozenset(e for e in E
                         if e in S or rk_of[frozenset(S | {e})] == r)
    from collections import Counter as _Counter
    absorbed = []
    for Hf in H:
        subs3 = [frozenset(c) for c in itertools.combinations(sorted(Hf), 3)]
        assert len(subs3) == 4
        for t in subs3:
            assert t in set(uflats) and int(mu_u[t]) == -1
            assert t not in flatset and closure(t) == Hf
        below = [f for f in flats if f < Hf]
        assert _Counter((rk_of[f], len(f)) for f in below) == _Counter(
            {(0, 0): 1, (1, 1): 4, (2, 2): 6})
        assert int(mu_v[Hf]) == -3 == (-1) ** 3 * (4 - 1)
        absorbed.extend(subs3)
    assert len(set(absorbed)) == 20
    assert w_v[3] == 36 * (-1) + 5 * (-3) == -51
    assert w_v[3] - w_u[3] == 5  # (-1)^4 = +1 per defect, r = 4
    print("  per-defect Mobius proof: w_3 -56 -> -51, m_3 35 -> 30: OK")
    # Try It 4
    X, Y = A | B, A | C
    assert X in flatset and Y in flatset and not (X <= Y or Y <= X)
    # sec 15.6 counts
    proper = [f for f in flats if f and f != frozenset(E)]
    assert len(proper) == 77 and len(flats) == 79
    # sec 15.5 ordered pairs
    cnt = sum(1 for X_ in allsubs for Y_ in allsubs
              if not (X_ <= Y_ or Y_ <= X_))
    assert cnt == 52670 and cnt // 2 == 26335
    print("  Try Its 1-4, Whitney/chi/chi-bar both matroids, sec 15.5/15.6")
    print("  counts: ALL PHASE-2 CHECKS PASSED")

verify_phase2()
print("\nALL CHAPTER-15 ASSERTIONS PASSED")

