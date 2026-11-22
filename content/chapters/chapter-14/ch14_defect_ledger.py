#!/usr/bin/env python3
"""Chapter 14 (Deep Dive 1) -- The Defect Ledger.  SITE BUNDLE program.

Merged from the study programs of Days 12, 14 and 15, plus the
Phase 2 verification groups (Try It answers, the exact Day-12
spectrum factorisation, the per-flat w_3 bookkeeping that closes
the 1 -> 3 jump, and the circuit-hyperplane sibling family).
Re-run end-to-end from this bundle at ship time (2026-11-22):

  Day 12 (2026-10-22): the Kahler package off the uniform locus at
    rank 3 -- M_line, one 3-point line {0,1,2} on 6 points, dims
    (1,14,1), m=(1,5,9) three ways, Hodge-Riemann on a 13-dim P^1.
  Day 14 (2026-10-24): the sparse chain-monomial Chow ring -- a
    monomial is zero unless its flats form a chain (CH2), so degree 3
    for U(4,5) needs 245 chain monomials, not 2925 multiset monomials
    (2300 is the squarefree count -- the distinction is Duet's
    referee fix, 2026-11-21);
    U(4,5) dims (1,21,21,1), m=(1,4,6,4), genuine Hard Lefschetz
    l: A^1 -> A^2 an isomorphism, HRR on a 20-dim P^1.
  Day 15 (2026-10-25): Bjorner's initial-descending flag count
    (Baker Prop 4.3) generalised to d=3 as a third route to m, and
    M_line^(4), the rank-4 matroid with one collinear triple:
    dims (1,27,27,1), m=(1,5,9,7) three ways, full Kahler package.

Headline ledger line (asserted at the end):
  U(4,6) m=(1,5,10,10)  ->  M_line^(4) m=(1,5,9,7):
  one collinear triple costs exactly (0,0,1,3);
  mu(line flat)=+2 (uniform: +1), mu(E)=7 (uniform: 10).
Rank-3 precursor (also asserted): U(3,6) (1,5,10) -> M_line (1,5,9),
cost (0,0,1).

Presentation: Feichtner-Yuzvinsky A*(M) = Z[x_F]/(CH1+CH2), Baker
"Hodge Theory in Combinatorics" arXiv:1705.07960 sec 4.1/4.3/4.4,
Thm 4.1/4.2 (Adiprasito-Huh-Katz).  All arithmetic exact
(Fractions / Bareiss integers); eigenvalues float, for figures only.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))


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


# ---------- Bjorner's initial descending flag count (rank-general, from Day 15) ----------'s initial descending flag count (rank-general) ----------
def bjorner_count(E, flats_all, rank_of, d):
    """Returns (counts, examples): counts[k] = m_k for k=0..d, per Baker
    Prop 4.3; examples[k] = first few maximal-length... (first few flags
    of that length, as tuples of sorted tuples)."""
    Fp = [f for f in flats_all if f and f != frozenset(E)]
    counts, examples = [], {}
    for k in range(d + 1):
        if k == 0:
            counts.append(1)
            examples[0] = [()]
            continue
        cnt, ex = 0, []
        for chain in combinations(Fp, k):
            if not all(chain[i] < chain[i + 1] for i in range(k - 1)):
                continue
            if not all(rank_of(chain[i]) == i + 1 for i in range(k)):
                continue
            mins = [min(f) for f in chain]
            if (all(mins[i] > mins[i + 1] for i in range(k - 1))
                    and mins[-1] > 0):
                cnt += 1
                if len(ex) < 6:
                    ex.append(tuple(tuple(sorted(f)) for f in chain))
        counts.append(cnt)
        examples[k] = ex
    return counts, examples


# ---------- matroid constructions ----------
def uniform_flats(E, r):
    flats = [frozenset(c) for s in range(r) for c in combinations(E, s)]
    flats.append(frozenset(E))
    return flats, (lambda f: min(len(f), r))

def m_line_rank3():
    """Day 12's M_line: rank 3 on {0..5}, one 3-point line {0,1,2}."""
    E = list(range(6)); line = frozenset({0, 1, 2})
    rk = {frozenset(): 0, frozenset(E): 3}
    flats = []
    for i in E:
        f = frozenset([i]); flats.append(f); rk[f] = 1
    flats.append(line); rk[line] = 2
    for c in combinations(E, 2):
        if set(c) <= set(line): continue
        f = frozenset(c); flats.append(f); rk[f] = 2
    # rank-2 flats partition the 15 pairs (linear-space axiom)
    seen = {}
    for f in flats:
        if rk[f] == 2:
            for p in combinations(sorted(f), 2):
                assert p not in seen; seen[p] = f
    assert len(seen) == 15
    return E, [frozenset()] + flats + [frozenset(E)], (lambda f: rk[f])

def m_line_rank4():
    """Day 15's M_line^(4): v0=e1, v1=e2, v2=e1+e2, v3=e3, v4=e4,
    v5=(2,3,5,7) over Q.  The vector model is PROVED to give
    r(S) = min(4, |S| - [{0,1,2} subset of S]) on all 64 subsets."""
    E = list(range(6)); LINE = {0, 1, 2}
    V = {0:(1,0,0,0), 1:(0,1,0,0), 2:(1,1,0,0), 3:(0,0,1,0),
         4:(0,0,0,1), 5:(2,3,5,7)}
    def vec_rank(S):
        m = [[Fraction(V[i][r]) for i in S] for r in range(4)]
        rk = 0; row = 0
        for c in range(len(S)):
            piv = next((i for i in range(row, 4) if m[i][c] != 0), None)
            if piv is None: continue
            m[row], m[piv] = m[piv], m[row]
            pv = m[row][c]; m[row] = [x / pv for x in m[row]]
            for i in range(4):
                if i != row and m[i][c] != 0:
                    f = m[i][c]; m[i] = [a - f*b for a, b in zip(m[i], m[row])]
            row += 1; rk += 1
        return rk
    allsubs = [frozenset(c) for s in range(7) for c in combinations(E, s)]
    rk_of = {S: vec_rank(sorted(S)) for S in allsubs}
    for S in allsubs:
        assert rk_of[S] == min(4, len(S) - (1 if LINE <= set(S) else 0)), S
    for A in allsubs:
        for B in allsubs:
            assert rk_of[A|B] + rk_of[A&B] <= rk_of[A] + rk_of[B]
    circuits = []
    for S in allsubs:
        if S and rk_of[S] < len(S) and all(
                rk_of[frozenset(x for x in S if x != e)] == len(S)-1 for e in S):
            circuits.append(tuple(sorted(S)))
    assert sorted(circuits) == [(0,1,2),(0,1,3,4,5),(0,2,3,4,5),(1,2,3,4,5)], circuits
    def closure(S):
        S = set(S); r = rk_of[frozenset(S)]
        return frozenset(e for e in E if e in S or rk_of[frozenset(S|{e})] == r)
    flats = sorted({closure(S) for S in allsubs},
                   key=lambda f: (rk_of[f], len(f), sorted(f)))
    by_rank = {r: sum(1 for f in flats if rk_of[f] == r) for r in range(5)}
    assert by_rank == {0:1, 1:6, 2:13, 3:13, 4:1}, by_rank
    return E, flats, (lambda f: rk_of[f]), by_rank

# ---------- triple agreement: FY == Mobius == Bjorner ----------
def triple(name, E, flats, rank_of, expect_dims, expect_m):
    CR, alpha, beta, a_lift, b_lift, m_fy = analyse(
        name, E, flats, rank_of, expect_dims, expect_m)
    bj, _ex = bjorner_count(E, flats, rank_of, CR.d)
    assert bj == m_fy == expect_m, (name, bj, m_fy)
    print(f"  {name}: FY == Mobius == Bjorner == {tuple(expect_m)}: TRIPLE AGREEMENT")
    return CR, a_lift, b_lift

def ample_lift(CR, E):
    n = len(E)
    cfun = lambda S: len(S) * (n - len(S))
    allsubs = [frozenset(c) for s in range(n+1) for c in combinations(E, s)]
    cnt = 0
    for A in allsubs:
        for B in allsubs:
            if not (A <= B or B <= A):
                assert cfun(A|B) + cfun(A&B) < cfun(A) + cfun(B); cnt += 1
    print(f"  c(S)=|S|({n}-|S|) strictly submodular: {cnt} incomparable pairs")
    return {i: Fraction(cfun(f)) for f, i in CR.gidx.items()}

def sylvester_positive(G, label):
    """G: symmetric Fraction matrix.  Clear denominators by L>0 and
    assert every leading principal minor > 0 (exact Bareiss)."""
    L = 1
    for r in G:
        for v in r: L = L * v.denominator // gcd(L, v.denominator)
    Gi = [[int(v*L) for v in r] for r in G]
    assert all(Gi[i][j] == Gi[j][i] for i in range(len(Gi)) for j in range(len(Gi)))
    minors = [det_bareiss([r[:k] for r in Gi[:k]]) for k in range(1, len(Gi)+1)]
    assert all(v > 0 for v in minors), (label, minors)
    print(f"  {label}: all {len(minors)} Sylvester minors > 0 (min {min(minors)}, L={L}): OK")
    return minors, L

def kahler_rank3(CR, E, label):
    """d=2: deg(l^2)>0 and -deg(.,.) positive definite on P^1=ker(l)."""
    l_lift = ample_lift(CR, E)
    def top_of(a_coords, b_lift):
        return CR.deg(CR.mul(1, a_coords, b_lift))
    deg_l2 = top_of({0: Fraction(1)} if False else CR.mul(0, {0: Fraction(1)}, l_lift), l_lift)
    assert deg_l2 > 0; print(f"  {label}: deg(l^2) = {deg_l2} > 0: OK")
    row = [top_of({i: Fraction(1)}, l_lift) for i in range(len(CR.basis[1]))]
    piv = next(i for i, v in enumerate(row) if v != 0)
    prim = []
    for j in range(len(CR.basis[1])):
        if j == piv: continue
        v = {j: Fraction(1), piv: -row[j]/row[piv]}
        assert top_of(v, l_lift) == 0; prim.append(v)
    print(f"  {label}: dim P^1 = {len(prim)}")
    G = [[-CR.deg(CR.mul(1, prim[j], CR.lift1(prim[i])))
          for j in range(len(prim))] for i in range(len(prim))]
    return sylvester_positive(G, f"{label} Hodge-Riemann on P^1")

def kahler_rank4(CR, E, label, expect_deg_l3):
    """d=3: PD pairing, HL l:A^1->A^2 iso, deg(l^3)>0, HRR on P^1=ker(l^2)."""
    l_lift = ample_lift(CR, E)
    n1 = len(CR.basis[1])
    pair = [[Fraction(0)]*n1 for _ in range(n1)]
    for i in range(n1):
        lift_i = {CR.basis[1][i]: Fraction(1)}
        for j in range(n1):
            pair[i][j] = CR.deg(CR.mul(2, {j: Fraction(1)}, lift_i))
    assert dense_rank(pair) == n1; print(f"  {label}: Poincare pairing nondegenerate (rank {n1}): OK")
    assert dense_rank(CR.op_matrix(1, l_lift)) == n1
    print(f"  {label}: Hard Lefschetz l: A^1 -> A^2 is an isomorphism (rank {n1}): OK")
    cur = {0: Fraction(1)}; dd = 0
    for _ in range(3):
        cur = CR.mul(dd, cur, l_lift); dd += 1
    deg_l3 = CR.deg(cur)
    assert deg_l3 == expect_deg_l3 and deg_l3 > 0
    print(f"  {label}: deg(l^3) = {deg_l3} > 0: OK")
    def l2_func(a):
        return CR.deg(CR.mul(2, CR.mul(1, a, l_lift), l_lift))
    row = [l2_func({i: Fraction(1)}) for i in range(n1)]
    piv = next(i for i, v in enumerate(row) if v != 0)
    prim = []
    for j in range(n1):
        if j == piv: continue
        v = {j: Fraction(1), piv: -row[j]/row[piv]}
        assert l2_func(v) == 0; prim.append(v)
    print(f"  {label}: dim P^1 = {len(prim)}")
    G = [[-CR.deg(CR.mul(2, CR.mul(1, prim[j], l_lift), CR.lift1(prim[i])))
          for j in range(len(prim))] for i in range(len(prim))]
    return sylvester_positive(G, f"{label} Hodge-Riemann on P^1")

def main():
    summary = {}
    # Controls + rank-3 ledger precursor (Days 12/15 controls)
    E = [0,1,2,3]; flats, rk = uniform_flats(E, 3)
    triple("U(3,4) control", E, flats, rk, (1,7,1), [1,3,3])
    E = [0,1,2,3,4]; flats, rk = uniform_flats(E, 3)
    triple("U(3,5)", E, flats, rk, (1,11,1), [1,4,6])
    E = list(range(6)); flats, rk = uniform_flats(E, 3)
    _cr, _a, _b = triple("U(3,6)", E, flats, rk, (1,16,1), [1,5,10])
    E, flats, rk = m_line_rank3()
    CR, _a, _b = triple("M_line rank 3 (Day 12)", E, flats, rk, (1,14,1), [1,5,9])
    kahler_rank3(CR, E, "M_line rank 3")
    summary["rank3_cost"] = [u - m for u, m in zip([1,5,10], [1,5,9])]
    assert summary["rank3_cost"] == [0,0,1]
    print(f"  RANK-3 LEDGER: U(3,6) (1,5,10) -> M_line (1,5,9), cost (0,0,1)")

    # Rank 4: sparse ring, uniform then the defect (Days 14-15)
    E = [0,1,2,3,4]; flats, rk = uniform_flats(E, 4)
    CR, _a, _b = triple("U(4,5) (Day 14)", E, flats, rk, (1,21,21,1), [1,4,6,4])
    assert [len(CR.mons[k]) for k in range(4)] == [1,25,105,245]
    kahler_rank4(CR, E, "U(4,5)", 5280)
    E = list(range(6)); flats, rk = uniform_flats(E, 4)
    mb, _mu = mobius_chibar(E, flats, rk); bj, _ = bjorner_count(E, flats, rk, 3)
    assert mb == bj == [1,5,10,10], (mb, bj)
    print("  U(4,6): Mobius == Bjorner == (1,5,10,10): OK")
    E, flats, rk, by_rank = m_line_rank4()
    CR, _a, _b = triple("M_line^(4) (Day 15)", E, flats, rk, (1,27,27,1), [1,5,9,7])
    _mb, mu = mobius_chibar(E, flats, rk)
    assert mu[frozenset({0,1,2})] == 2 and mu[frozenset(E)] == 7, (
        mu[frozenset({0,1,2})], mu[frozenset(E)])
    print("  mu(line flat) = +2 (uniform +1), mu(E) = 7 (uniform 10): OK")
    kahler_rank4(CR, E, "M_line^(4)", 18144)
    cost = [u - m for u, m in zip([1,5,10,10], [1,5,9,7])]
    assert cost == [0,0,1,3], cost
    summary["rank4_cost"] = cost
    print("  HEADLINE LEDGER: one collinear triple costs m exactly (0,0,1,3)")
    verify_w3_bookkeeping()
    verify_circuit_hyperplanes()
    verify_try_it()
    verify_spectrum_exact()
    with open(os.path.join(HERE, "ch14_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print("\nALL CHAPTER-14 ASSERTIONS PASSED")



# ---------- Phase 2 verification groups (merged at ship) ----------

# Day 12's integral 13x13 Gram matrix of -deg(.,.) on P^1 for M_line,
# in Day 12's own integral basis (Gi in day12_nonuniform_chow.py).
# Embedded here so the Sec 14.2 spectrum claim has a runnable receipt
# inside this bundle (Duet referee note 3, 2026-11-21).
DAY12_GRAM = [[3472, 77, 77, 77, 77, 77, 77, 77, 77, 77, -300, -300, -300], [77, 3036, 509, 509, 870, -1296, -1296, 870, -1296, -1296, 117, 117, -225], [77, 509, 3036, 509, -1296, 870, -1296, -1296, 870, -1296, 117, -225, 117], [77, 509, 509, 3036, -1296, -1296, 870, -1296, -1296, 870, -225, 117, 117], [77, 870, -1296, -1296, 3036, 509, 509, 870, -1296, -1296, 117, 117, -225], [77, -1296, 870, -1296, 509, 3036, 509, -1296, 870, -1296, 117, -225, 117], [77, -1296, -1296, 870, 509, 509, 3036, -1296, -1296, 870, -225, 117, 117], [77, 870, -1296, -1296, 870, -1296, -1296, 3036, 509, 509, 117, 117, -225], [77, -1296, 870, -1296, -1296, 870, -1296, 509, 3036, 509, 117, -225, 117], [77, -1296, -1296, 870, -1296, -1296, 870, 509, 509, 3036, -225, 117, 117], [-300, 117, 117, -225, 117, 117, -225, 117, 117, -225, 78, 15, 15], [-300, 117, -225, 117, 117, -225, 117, 117, -225, 117, 15, 78, 15], [-300, -225, 117, 117, -225, 117, 117, -225, 117, 117, 15, 15, 78]]

def verify_spectrum_exact():
    import sympy as sp
    x = sp.symbols("x")
    M = sp.Matrix(DAY12_GRAM)
    assert sum(DAY12_GRAM[i][i] for i in range(13)) == 31030
    expect = ((x - 5776)**2 * (x - 361)**4
              * (x**2 - 6922*x + 81225)**2
              * (x**3 - 4190*x**2 + 2235172*x - 53686476))
    got = sp.expand(M.charpoly(x).as_expr())
    assert sp.expand(got - expect) == 0, "charpoly factorisation mismatch"
    print("  Sec 14.2 spectrum: charpoly of Day-12 Gram factors over Z as")
    print("    (x-5776)^2 (x-361)^4 (x^2-6922x+81225)^2")
    print("    (x^3-4190x^2+2235172x-53686476); 361=19^2 (x4),")
    print("    5776=76^2 (x2) are EXACT eigenvalues: OK")

def verify_w3_bookkeeping():
    """Duet fix 2: the 1 -> 3 jump, closed per-flat (Sec 14.7 Honest Box)."""
    E = list(range(6)); line = frozenset({0, 1, 2})
    fu, ru = uniform_flats(E, 4)
    _m, mu_u = mobius_chibar(E, fu, ru)
    Ed, fd, rd, _b = m_line_rank4()
    _m2, mu_d = mobius_chibar(Ed, fd, rd)
    w3u = sum(int(mu_u[f]) for f in fu if ru(f) == 3)
    w3d = sum(int(mu_d[f]) for f in fd if rd(f) == 3)
    assert (w3u, w3d) == (-20, -16), (w3u, w3d)
    # (+1) the line itself leaves rank 3
    assert ru(line) == 3 and mu_u[line] == -1 and rd(line) == 2
    # (+3) nine 2-line-point triples (mu -1 each in U) are not flats in
    # the defect matroid; three line+point planes carry mu -2 each
    nine = [f for f in fu if ru(f) == 3 and len(set(f) & set(line)) == 2]
    assert len(nine) == 9 and all(mu_u[f] == -1 for f in nine)
    assert all(f not in set(fd) for f in nine)
    planes = [f for f in fd if rd(f) == 3 and line <= set(f)]
    assert len(planes) == 3 and all(mu_d[f] == -2 for f in planes)
    rest_u = [f for f in fu if ru(f) == 3 and f != line and f not in nine]
    rest_d = [f for f in fd if rd(f) == 3 and f not in planes]
    assert len(rest_u) == len(rest_d) == 10
    assert all(mu_u[f] == -1 for f in rest_u) and all(mu_d[f] == -1 for f in rest_d)
    assert (w3d - w3u) == 1 + (-6 - (-9)) == 4
    print("  Sec 14.7 bookkeeping: w_3 -20 -> -16 = +1 (line leaves rank 3)")
    print("    + 3 (nine triples -9 merge into three planes -6);")
    print("    m_3 = |9-16| = 7 vs |10-20| = 10: OK")

def verify_circuit_hyperplanes():
    """Sec 14.7 sibling family: k circuit-hyperplanes cost (0,0,0,k)."""
    E = list(range(6))
    defects = [frozenset({0,1,2,3}), frozenset({0,1,4,5}),
               frozenset({2,3,4,5})]
    for k in (1, 2, 3):
        D = defects[:k]
        def rk(S, D=D):
            S = frozenset(S)
            return 3 if S in D else min(4, len(S))
        allsubs = [frozenset(c) for s in range(7) for c in combinations(E, s)]
        def closure(S, rk=rk):
            S = set(S); r = rk(S)
            return frozenset(e for e in E if e in S or rk(set(S) | {e}) == r)
        flats = sorted({closure(S) for S in allsubs},
                       key=lambda f: (rk(f), len(f), sorted(f)))
        mb, _mu = mobius_chibar(E, flats, rk)
        bj, _ex = bjorner_count(E, flats, rk, 3)
        assert mb == bj == [1, 5, 10, 10 - k], (k, mb, bj)
        print(f"  circuit-hyperplanes k={k}: Mobius == Bjorner =="
              f" (1,5,10,{10-k}), cost (0,0,0,{k}): OK")

def verify_try_it():
    assert 5*5 >= 1*9 and 9*9 >= 5*7          # Try It 1
    E, flats, rk, _b = m_line_rank4()          # Try It 3
    fs = set(flats)
    a, b = frozenset({0, 1, 2}), frozenset({3, 4})
    assert a in fs and b in fs and not (a <= b or b <= a)
    E5 = [0, 1, 2, 3, 4]                       # Try It 4
    f5, r5 = uniform_flats(E5, 4)
    _c, ex = bjorner_count(E5, f5, r5, 3)
    assert ((4,), (3, 4), (2, 3, 4)) in ex[3], ex[3]
    print("  Try It 1-4 answers asserted (log-concavity, ledger,")
    print("    incomparable flats, Bjorner flag membership): OK")


if __name__ == "__main__":
    main()
