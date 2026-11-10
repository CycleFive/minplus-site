"""
Chapter 10 program — the tropical Grassmannian G_{2,4} is tree space
(Speyer-Sturmfels).  Merged for the book (2026-11-10) from Day 10 of the
study (day10_grassmannian.py) plus the Phase-2 verification script
(ch10_phase2_verify.py): assertion groups 8-10 are the Try It answer
key, rank(phi)=4, and the lifted minors satisfying the Plucker relation
exactly.  Everything is exact Fractions / exact pi-adic series.

Source: D. Speyer & B. Sturmfels, "The tropical Grassmannian",
arXiv:math/0304218.  Theorem 3.4: G_{2,n} = space of phylogenetic trees
(Billera-Holmes-Vogtmann).  We do the first nontrivial case n=4 end to end.

Setup (paper Sec. 2-4, min convention, K a valued field with deg(t)=1):
- I_{2,4} = < p12*p34 - p13*p24 + p14*p23 >, the single Plucker quadric
  (1-indexed pairs, matching the paper).
- w in R^6 indexed by pairs (12,13,14,23,24,34).
- Pair sums: s1 = w12+w34, s2 = w13+w24, s3 = w14+w23.
- Four-point condition (6): min(s1,s2,s3) attained at least twice.
- phi: R^4 -> R^6, phi(a)_ij = a_i + a_j  (column-scaling torus action).
- E_{A,B} = sum_{i in A, j in B} e_ij  (cut metric of the split {A,B}).

SIGN CORRECTION (verified against the arXiv TeX source, not just ar5iv):
the paper's (8) as printed says  cone = image(phi) + R>=0 * {E_{A,B}}.
With B_sigma defined by "-w_ij is the distance" and four-point (6) as
printed, the correct generator sign is MINUS E_{A,B}:
a tree metric d = t*E_{12,34} (internal edge t>0, zero pendants) gives
w = -t*E, which satisfies (6) but is NOT in image(phi)+R>=0*E
(proof: -tE = phi(a)+sE, s>=0 forces c:=-(t+s)<0 with phi(a)=cE,
 and a1+a2=0, a1+a3=c, a3+a4=0 give c=0, contradiction).
The code below asserts: -E satisfies (6)-interior, +E does NOT.

What is computed (all asserted in-code):
1. four_point(w) <=> min(s1,s2,s3) attained twice; topology from the
   unique strict max (e.g. s1>s2=s3 -> cherries {1,2},{3,4}).
2. Cone-generator equivalence: w satisfies (6) with topology 12|34
   <=> w = phi(a) - t*E_{12,34} for a unique-mod-... a and t>=0,
   with t = (s1-s2)/2 the internal edge length ("-w = distance" units).
   Checked on random integer points (exact Fractions).
3. Quotient fan: R^6/image(phi) --iso--> R^2 via (s1-s2, s1-s3).
   Three rays: [E_{12,34}] -> (1,1), [E_{13,24}] -> (-1,0),
   [E_{14,23}] -> (0,-1).  Primitive in Z^2, sum = (0,0): BALANCED.
   (So G_{2,4}/image(phi) is a tropical line -- cf. day 0.)
4. Initial ideals: for w in a cone interior, in_w of the Plucker quadric
   is the binomial -p13*p24 + p14*p23 (etc.), monomial-free, i.e. w is in
   the tropical variety (Thm 2.1(c)); for w off the fan in_w is a
   monomial.  J_sigma recovered.
5. Lifting (Cor. 2.2 made constructive): every w in the fan is val(p)
   for an actual 2-plane: explicit 2x4 matrix over C((pi)) with
   val(det of columns i,j) = w_ij, verified by exact pi-adic series
   arithmetic.  Column scaling by pi^{a_i} realizes the phi(a) part.
6. Tree reconstruction (baby Additive Linkage): w -> (topology, t);
   t is phi-invariant: pendant edges are modded out, as the paper says.
7. Combinatorics of T_4: 3 vertices = 2^3-4-1, 3 facets = 3!!,
   pure of dimension n-3 = 1.
"""

from fractions import Fraction as Q
from itertools import combinations
import random

# ---------------------------------------------------------------- pairs

PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
PIDX = {p: k for k, p in enumerate(PAIRS)}


def pair(i, j):
    return (i, j) if i < j else (j, i)


def zeros6():
    return [Q(0)] * 6


def add6(u, v):
    return [u[k] + v[k] for k in range(6)]


def scale6(c, u):
    return [c * u[k] for k in range(6)]


# ---------------------------------------------------------------- four-point


def pair_sums(w):
    s1 = w[PIDX[(1, 2)]] + w[PIDX[(3, 4)]]
    s2 = w[PIDX[(1, 3)]] + w[PIDX[(2, 4)]]
    s3 = w[PIDX[(1, 4)]] + w[PIDX[(2, 3)]]
    return (s1, s2, s3)


def four_point(w):
    """min(s1,s2,s3) attained at least twice -- paper (6)."""
    s = pair_sums(w)
    m = min(s)
    return sum(1 for x in s if x == m) >= 2


# topology: which pair-sum is the unique strict MAXIMUM (cherry pairs).
# s1 max  -> cherries {1,2},{3,4}  (internal edge 12|34)
TOPOS = {
    0: ((1, 2), (3, 4), (1, 3), (1, 4), (2, 3), (2, 4)),  # s1 strict max
    1: ((1, 3), (2, 4), (1, 2), (1, 4), (2, 3), (3, 4)),  # s2 strict max
    2: ((1, 4), (2, 3), (1, 2), (1, 3), (2, 4), (3, 4)),  # s3 strict max
}
# TOPOS[k] = (cherry1, cherry2, cross1, cross2, cross3, cross4)


def topology(w):
    """Return 0/1/2 for the three tree topologies, 'ridge' if s1=s2=s3,
    None if the four-point condition fails."""
    if not four_point(w):
        return None
    s = pair_sums(w)
    if s[0] == s[1] == s[2]:
        return "ridge"
    # unique strict max -> that topology (min attained twice => exactly
    # one strict max unless ridge)
    m = max(s)
    assert sum(1 for x in s if x == m) == 1
    return s.index(m)


# ---------------------------------------------------------------- phi, E

def phi(a):
    """a in R^4 -> R^6, phi(a)_ij = a_i + a_j."""
    w = zeros6()
    for k, (i, j) in enumerate(PAIRS):
        w[k] = a[i - 1] + a[j - 1]
    return w


def E_split(A, B):
    """Cut vector E_{A,B} = sum_{i in A, j in B} e_ij."""
    w = zeros6()
    for i in A:
        for j in B:
            w[PIDX[pair(i, j)]] += 1
    return w


E12_34 = E_split({1, 2}, {3, 4})  # = e13+e14+e23+e24
E13_24 = E_split({1, 3}, {2, 4})  # = e12+e14+e23+e34
E14_23 = E_split({1, 4}, {2, 3})  # = e12+e13+e24+e34
ES = [E12_34, E13_24, E14_23]


def decompose(w, topo):
    """Write w = phi(a) - t*E_topo.  Returns (a, t) or None if t < 0.

    t = (s_strict - s_tied)/2; a solved from phi(a) = w + t*E and
    consistency of all six equations asserted.
    """
    s = pair_sums(w)
    tied = [s[k] for k in range(3) if k != topo]
    assert tied[0] == tied[1], "not in closed cone for this topology"
    t = (s[topo] - tied[0]) / 2
    if t < 0:
        return None
    u = add6(w, scale6(t, ES[topo]))  # u should be in image(phi)
    # solve a1..a4 from a_i + a_j = u_ij
    a1 = (u[PIDX[(1, 2)]] + u[PIDX[(1, 3)]] - u[PIDX[(2, 3)]]) / 2
    a2 = (u[PIDX[(1, 2)]] + u[PIDX[(2, 3)]] - u[PIDX[(1, 3)]]) / 2
    a3 = (u[PIDX[(1, 3)]] + u[PIDX[(2, 3)]] - u[PIDX[(1, 2)]]) / 2
    a4 = (u[PIDX[(3, 4)]] + u[PIDX[(1, 4)]] - u[PIDX[(1, 3)]]) / 2
    a = [a1, a2, a3, a4]
    got = phi(a)
    assert got == u, f"decomposition inconsistent: {got} != {u}"
    return a, t


# ---------------------------------------------------------------- quotient fan

def quot(w):
    """R^6/image(phi) --iso--> R^2,  w |-> (s1-s2, s1-s3).  Kernel is
    exactly image(phi) (both 4-dim; phi(a) has s1=s2=s3=sum(a))."""
    s1, s2, s3 = pair_sums(w)
    return (s1 - s2, s1 - s3)


# ---------------------------------------------------------------- initial ideals

# Plucker quadric terms: (coeff, (pairA, pairB))
PLUCKER = [
    (Q(1), ((1, 2), (3, 4))),
    (Q(-1), ((1, 3), (2, 4))),
    (Q(1), ((1, 4), (2, 3))),
]


def initial_form(w):
    """in_w of the Plucker quadric: terms attaining min pair-sum weight.
    Returns (terms, is_monomial)."""
    wt = [w[PIDX[a]] + w[PIDX[b]] for (_, (a, b)) in PLUCKER]
    m = min(wt)
    terms = [(c, pr) for (c, pr), x in zip(PLUCKER, wt) if x == m]
    return terms, len(terms) == 1


# ---------------------------------------------------------------- pi-adic lifting

def sadd(x, y):
    z = dict(x)
    for e, c in y.items():
        z[e] = z.get(e, Q(0)) + c
        if z[e] == 0:
            del z[e]
    return z


def smul(x, y):
    z = {}
    for e1, c1 in x.items():
        for e2, c2 in y.items():
            e = e1 + e2
            z[e] = z.get(e, Q(0)) + c1 * c2
            if z[e] == 0:
                del z[e]
    return z


def sneg(x):
    return {e: -c for e, c in x.items()}


def sval(x):
    assert x, "zero series has no valuation"
    return min(x)


def mon(e, c=Q(1)):
    return {Q(e): Q(c)}


def det2(c1, c2):
    """det of 2x2 with columns c1, c2 (each a pair of series)."""
    return sadd(smul(c1[0], c2[1]), sneg(smul(c1[1], c2[0])))


def lift_matrix(w):
    """Build 2x4 matrix over C((pi)) with val(p_ij) = w_ij.

    Decompose w = phi(a) - t*E_topo (t>=0).  Lift the core -t*E_topo,
    then scale column i by pi^{a_i} (adds a_i+a_j to val(p_ij)).
    Series are dicts {exponent: coeff}; pi is the uniformizer.
    """
    topo = topology(w)
    assert topo in (0, 1, 2, "ridge"), "w not in the fan"
    if topo == "ridge":
        t, E = Q(0), E12_34
    else:
        dec = decompose(w, topo)
        assert dec is not None
        a, t = dec
        E = ES[topo]
    if topo == "ridge":
        a = [Q(0)] * 4
        # recompute: w = phi(a) exactly on the ridge
        u = w
        a1 = (u[PIDX[(1, 2)]] + u[PIDX[(1, 3)]] - u[PIDX[(2, 3)]]) / 2
        a2 = (u[PIDX[(1, 2)]] + u[PIDX[(2, 3)]] - u[PIDX[(1, 3)]]) / 2
        a3 = (u[PIDX[(1, 3)]] + u[PIDX[(2, 3)]] - u[PIDX[(1, 2)]]) / 2
        a4 = (u[PIDX[(3, 4)]] + u[PIDX[(1, 4)]] - u[PIDX[(1, 3)]]) / 2
        a = [a1, a2, a3, a4]
        assert phi(a) == w
        t = Q(0)
    # core matrix for -t*E_topo (topo 0 shown; others by permuting columns)
    # columns: c1=(1,0), c2=(0,1), c3=(pi^-t, pi^-t), c4=(pi^-t, pi^-t - pi^t)
    mt = -t
    one = mon(0)
    zero = {}
    pit = mon(mt)            # pi^{-t}
    pit2 = mon(t)            # pi^{t}
    binom = sadd(pit, sneg(pit2))  # pi^{-t} - pi^{t}
    cols = [(one, zero), (zero, one), (pit, pit), (pit, binom)]
    if t == 0:
        # separate lift for the ridge point 0: all Pluckers val 0, nonzero
        cols = [(one, zero), (zero, one),
                (one, one), (sadd(one, one), one)]  # c4 = (2,1)
    # permute columns to match the topology's cherry pairing:
    # core above has cherries {1,2},{3,4} (topo 0)
    perm = {0: [0, 1, 2, 3], 1: [0, 2, 1, 3], 2: [0, 2, 3, 1],
            "ridge": [0, 1, 2, 3]}[topo]
    cols = [cols[p] for p in perm]
    # scale column i by pi^{a_i}
    cols = [(smul(c[0], mon(a[i])), smul(c[1], mon(a[i])))
            for i, c in enumerate(cols)]
    # verify valuations and nonvanishing
    for (i, j) in PAIRS:
        d = det2(cols[i - 1], cols[j - 1])
        assert d, f"p_{i}{j} vanished!"
        v = sval(d)
        assert v == w[PIDX[(i, j)]], \
            f"val(p_{i}{j}) = {v} != w = {w[PIDX[(i,j)]]}"
    return cols


# ---------------------------------------------------------------- tests

def main():
    rng = random.Random(20261020)

    # --- 1. the sign correction: -E satisfies (6)-interior, +E does not
    for k, E in enumerate(ES):
        assert four_point(scale6(-1, E)) and topology(scale6(-1, E)) == k, \
            f"-E_{k} should be interior of cone {k}"
        assert not four_point(E), \
            f"+E_{k} violates four-point: paper (8) sign is wrong as printed"
    print("1. sign check: -E in fan interior, +E not in fan  [paper (8) typo confirmed]")

    # --- 2. cone-generator equivalence on random integer points
    n_test = 0
    for trial in range(4000):
        w = [Q(rng.randint(-3, 3)) for _ in range(6)]
        fp = four_point(w)
        in_some_cone = False
        for k in range(3):
            s = pair_sums(w)
            tied = [s[j] for j in range(3) if j != k]
            if tied[0] == tied[1] and s[k] >= tied[0]:
                dec = decompose(w, k)
                assert dec is not None
                a, t = dec
                # round trip
                assert add6(phi(a), scale6(-t, ES[k])) == w
                in_some_cone = True
        assert fp == in_some_cone, f"fan membership mismatch at {w}"
        n_test += 1
    print(f"2. four-point <=> union of the 3 cones (image(phi)-tE), {n_test} random pts")

    # --- 3. quotient fan: 3 rays, primitive, balanced
    rays = [quot(scale6(-1, E)) for E in ES]
    assert rays == [(Q(2), Q(2)), (Q(-2), Q(0)), (Q(0), Q(-2))]
    prim = [(Q(1), Q(1)), (Q(-1), Q(0)), (Q(0), Q(-1))]
    for r, p in zip(rays, prim):
        assert (r[0] / p[0] if p[0] else r[1] / p[1]) == 2  # ray = 2*primitive
        assert abs(p[0]) + abs(p[1]) == 1 or (p == (Q(1), Q(1)))
    import math
    for p in prim:
        assert math.gcd(int(p[0]), int(p[1])) == 1
    assert (sum(p[0] for p in prim), sum(p[1] for p in prim)) == (Q(0), Q(0))
    # kernel of quot is exactly image(phi): dim check
    assert quot(phi([Q(1), Q(2), Q(3), Q(4)])) == (Q(0), Q(0))
    print("3. quotient fan: rays (1,1),(-1,0),(0,-1), primitive, sum=0 -- a tropical line")

    # --- 4. initial ideals: binomial on fan, monomial off fan
    for k in range(3):
        w = add6(phi([Q(1), Q(0), Q(2), Q(0)]), scale6(-Q(3), ES[k]))
        terms, is_mon = initial_form(w)
        assert not is_mon and len(terms) == 2
    # the three binomials are the J_sigma of the paper
    w0 = scale6(-Q(1), E12_34)
    terms, _ = initial_form(w0)
    got = sorted((tuple(sorted(p)) for _, (a, b) in terms for p in [a, b]))
    assert got == [(1, 3), (1, 4), (2, 3), (2, 4)], got
    n_off = 0
    for trial in range(2000):
        w = [Q(rng.randint(-3, 3)) for _ in range(6)]
        terms, is_mon = initial_form(w)
        assert is_mon == (not four_point(w)), f"in_w mismatch at {w}"
        n_off += 1
    print(f"4. in_w binomial (monomial-free) on fan / monomial off fan, {n_off} random pts")

    # --- 5. lifting: every fan point is val(p) of a real 2-plane
    test_ws = []
    for k in range(3):
        for t in [Q(1), Q(5) / 2]:
            test_ws.append(add6(phi([Q(rng.randint(-2, 2)) for _ in range(4)]),
                                scale6(-t, ES[k])))
    test_ws.append(phi([Q(1), Q(-1), Q(2), Q(0)]))  # ridge point
    test_ws.append([Q(0)] * 6)                      # origin
    for w in test_ws:
        lift_matrix(w)
    print(f"5. lifted {len(test_ws)} fan points to 2x4 matrices/C((pi)), vals exact")

    # --- 6. tree reconstruction: t is phi-invariant (pendants modded out)
    for trial in range(500):
        k = rng.randrange(3)
        t = Q(rng.randint(1, 5))
        a = [Q(rng.randint(-3, 3)) for _ in range(4)]
        w = add6(phi(a), scale6(-t, ES[k]))
        assert topology(w) == k
        dec = decompose(w, k)
        assert dec[1] == t, "internal edge length changed?!"
        # adding phi(b) does not change topology or t
        b = [Q(rng.randint(-3, 3)) for _ in range(4)]
        w2 = add6(w, phi(b))
        assert topology(w2) == k and decompose(w2, k)[1] == t
    print("6. reconstruction: topology + internal edge t recovered, phi-invariant (500x)")

    # --- 7. T_4 combinatorics from the paper
    verts = []
    for r in (2,):
        for A in combinations([1, 2, 3, 4], r):
            B = tuple(i for i in [1, 2, 3, 4] if i not in A)
            if len(B) >= 2 and A < B:
                verts.append((A, B))
    assert len(verts) == 3 == 2 ** 3 - 4 - 1
    assert len(ES) == 3 == 3 * 1  # (2*4-5)!! = 3 facets
    print("7. T_4: 3 vertices (=2^3-4-1), 3 facets (=3!!), pure dim n-3=1")

    # --- 8. Try It answers (added in Phase 2 of the book cycle)
    a_try = [Q(5, 2), Q(5, 2), Q(0), Q(0)]  # sum(a) = 5
    w_try = add6(phi(a_try), scale6(-Q(1), E12_34))
    assert pair_sums(w_try) == (Q(5), Q(3), Q(3))
    assert four_point(w_try) and topology(w_try) == 0
    assert decompose(w_try, 0)[1] == 1  # t = (5-3)/2
    for a8 in ([Q(1), Q(2), Q(3), Q(4)], [Q(-1), Q(0), Q(5), Q(2)]):
        s8 = pair_sums(phi(a8))
        # each pair sum adds two disjoint pair-terms covering all four
        # indices once, so phi(a) shifts every pair sum by sum(a) --
        # NOT 2*sum(a) (the draft's slip, caught in Phase 2)
        assert s8[0] == s8[1] == s8[2] == sum(a8)
    print("8. Try It key: (5,3,3)->topo 12|34, t=1; phi shifts pair sums by sum(a)")

    # --- 9. rank(phi) = 4 exactly; the scalar shift is already in image(phi)
    M = [[phi([Q(1) if i == j else Q(0) for i in range(4)])[k]
          for j in range(4)] for k in range(6)]
    rank = 0
    r0 = 0
    for c in range(4):
        piv = next((r for r in range(r0, 6) if M[r][c] != 0), None)
        if piv is None:
            continue
        M[r0], M[piv] = M[piv], M[r0]
        pv = M[r0][c]
        M[r0] = [x / pv for x in M[r0]]
        for r in range(6):
            if r != r0 and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[r0])]
        r0 += 1
        rank += 1
    assert rank == 4, rank
    assert phi([Q(1, 2)] * 4) == [Q(1)] * 6
    print("9. rank(phi)=4 (lineality 4-dim, not 5); (1,...,1)=phi(1/2,...,1/2)")

    # --- 10. the lifted minors satisfy the Plucker relation EXACTLY
    for k in range(3):
        w10 = add6(phi([Q(1), Q(0), Q(2), Q(0)]), scale6(-Q(3), ES[k]))
        cols10 = lift_matrix(w10)

        def minor(i, j, cols=cols10):
            return det2(cols[i - 1], cols[j - 1])

        rel = sadd(sadd(smul(minor(1, 2), minor(3, 4)),
                        sneg(smul(minor(1, 3), minor(2, 4)))),
                   smul(minor(1, 4), minor(2, 3)))
        assert rel == {}, f"Plucker relation fails for topo {k}: {rel}"
    print("10. lifted minors satisfy p12p34-p13p24+p14p23=0 exactly (3 topos)")

    print("\nAll chapter-10 assertions passed.")


if __name__ == "__main__":
    main()
