"""Chapter 12: The Sandpile Group of a Donut -- chip-firing on C_6.

Baker-Norine theory (arXiv:math/0608360) on the 6-cycle inside the
tropical cubic of Chapter 3:
  - divisors, Laplacian, principal divisors, Jac(G) = Div^0/Prin,
  - q-reduced divisors (G-parking functions) -- unique per class,
  - the Abel-Jacobi map S_q(v) = [(v)-(q)], bijective for the
    2-edge-connected genus-1 cycle,
  - the group law: addition is rotation of the hexagon,
  - tropical Riemann-Roch r(D) - r(K-D) = deg(D) + 1 - g, brute-forced
    over all divisors with entries in {-1,0,1},
  - plus the book cycle's Phase-2 verification groups: the Try It
    answers, the three violating sets of the Honest Box counterexample,
    Dhar's burning test <=> q-reduced, the surjectivity route to
    |Jac| = 6 (Thm 1.7, independent of Matrix-Tree), what the 729-box
    tests when K = 0, the lattice-length computation of Sec. 12.8
    (Euclidean total 4+2*sqrt(2), lattice total L = 6), and the
    min-closure of L(D) (BN Remark 1.13(iii)).

Run as a script: every claim in the chapter is asserted in code.
Source study: Day 8 (2026-10-18); book Phase 2 verified 2026-11-15.
"""

import itertools
import math
from fractions import Fraction
from math import gcd

import numpy as np
from sympy.matrices.normalforms import smith_normal_form
from sympy import Matrix

N = 6                      # the 6-cycle, abstractly (isomorphism type of day-2's cycle)
Q = 0                      # base point q = v_0
EDGES = [(i, (i + 1) % N) for i in range(N)]

# Laplacian Q = D - A  (Baker-Norine (1.3): [Delta(f)] = Q[f])
LAP = np.zeros((N, N), dtype=int)
for i, j in EDGES:
    LAP[i, i] += 1
    LAP[j, j] += 1
    LAP[i, j] -= 1
    LAP[j, i] -= 1

G = len(EDGES) - N + 1


def outdeg(A, v):
    A = set(A)
    return sum(1 for i, j in EDGES if (i == v and j not in A) or (j == v and i not in A))


def is_q_reduced(D):
    """P1: no debt away from q; P2: every nonempty A subset V\\{q} has a
    vertex v with D(v) < outdeg_A(v) (Dhar burning condition)."""
    if any(D[v] < 0 for v in range(N) if v != Q):
        return False
    others = [v for v in range(N) if v != Q]
    for r in range(1, len(others) + 1):
        for A in itertools.combinations(others, r):
            if not any(D[v] < outdeg(A, v) for v in A):
                return False
    return True


# BFS layers from q (BN Prop 3.1: order vertices so every v != q has a
# neighbor preceding it; used by the out-of-debt phase)
_ADJ = {v: set() for v in range(N)}
for _i, _j in EDGES:
    _ADJ[_i].add(_j)
    _ADJ[_j].add(_i)
_DIST = {Q: 0}
_dq = [Q]
for _u in _dq:
    for _w in _ADJ[_u]:
        if _w not in _DIST:
            _DIST[_w] = _DIST[_u] + 1
            _dq.append(_w)
_ORDER = sorted(range(N), key=lambda v: _DIST[v])  # q first


def _fire(D, f, v, n=1):
    return D - n * LAP[:, v], f + n * np.eye(N, dtype=int)[:, v]


def stabilize(D, subset_order="lex", _trace=False):
    """BN Prop 3.1 (informal proof, ar5iv L801-848), implemented literally.

    Phase 1 -- out of debt: process v != q farthest-first; a preceding
    neighbor w (dist(w) = dist(v)-1) fires -D[v] times, taking v to 0.
    Already-fixed vertices only gain chips, so each vertex is fixed once.
    Phase 2 -- fire violating sets: while some nonempty A subset V\\{q}
    has D[v] >= outdeg_A(v) for all v in A, fire all of A at once.
    Each firing strictly decreases sum_{v != q} D[v] >= 0, so it terminates.
    Tracks the firing script f with D - Dq = LAP @ f exactly.
    With _trace=True, also returns the list of fired sets (for figures).
    """
    D = np.array(D, dtype=int)
    D0 = D.copy()
    f = np.zeros(N, dtype=int)
    fired_sets = []
    # Phase 1
    for v in reversed(_ORDER):
        if v == Q:
            continue
        while D[v] < 0:
            preds = [w for w in _ADJ[v] if _DIST[w] == _DIST[v] - 1]
            assert preds, f"no preceding neighbor for {v}"
            w = preds[0]
            n = -int(D[v])
            D, f = _fire(D, f, w, n)
    assert all(D[v] >= 0 for v in range(N) if v != Q), "phase 1 failed"
    # Phase 2
    others = [v for v in range(N) if v != Q]
    subsets = [A for r in range(1, N)
               for A in itertools.combinations(others, r)]
    if subset_order == "rev":
        subsets = subsets[::-1]
    for _ in range(1000000):
        fired = False
        for A in subsets:
            if all(D[v] >= outdeg(A, v) for v in A):
                if _trace:
                    fired_sets.append((set(A), D.copy()))
                for v in A:
                    D, f = _fire(D, f, v)
                fired = True
                break
        if not fired:
            break
    else:
        raise AssertionError("phase 2 did not terminate")
    assert D.sum() == D0.sum(), "degree changed during stabilization"
    assert is_q_reduced(D), f"stabilized divisor {D} is not q-reduced"
    assert np.array_equal(D0 - D, LAP @ f), "D - Dq not principal!"
    return (D, fired_sets) if _trace else D


def reduced_class(D, subset_order="lex"):
    """q-reduced representative of the class [D]; asserts D - Dq in Prin."""
    return stabilize(np.array(D, dtype=int), subset_order)


def same_class(D1, D2):
    return np.array_equal(reduced_class(D1), reduced_class(D2))


def is_effective_class(D):
    """|D| != empty  <=>  q-reduced rep is effective (firing preserves
    effectivity, borrowing never triggers on an effective divisor)."""
    return bool(np.all(reduced_class(D) >= 0))


def effective_divisors(s):
    for combo in itertools.combinations_with_replacement(range(N), s):
        D = np.zeros(N, dtype=int)
        for v in combo:
            D[v] += 1
        yield D


_red_cache = {}


def reduced_cached(D):
    key = tuple(int(x) for x in D)
    if key not in _red_cache:
        _red_cache[key] = reduced_class(np.array(D))
    return _red_cache[key]


def rank(D):
    """BN: r(D) >= s iff |D - E| != empty for every effective E of degree s."""
    D = np.array(D, dtype=int)
    d = int(D.sum())
    if d < 0:
        return -1
    for s in range(d + 1):
        for E in effective_divisors(s):
            if not np.all(reduced_cached(D - E) >= 0):
                return s - 1
    return d


def aj(i):
    """Abel-Jacobi: S_q(v_i) = [(v_i) - (q)], q-reduced rep."""
    D = np.zeros(N, dtype=int)
    D[i] += 1
    D[Q] -= 1
    return reduced_class(D)


def add_classes(c1, c2):
    return reduced_class(np.array(c1) + np.array(c2))


def main():
    print(f"C_{N}: V={N}, E={len(EDGES)}, genus g = {G}")
    assert G == 1
    # 2-edge-connected: every cut has >= 2 edges (needed for AJ injectivity, Thm 1.8)
    for r in range(1, N):
        for S in itertools.combinations(range(N), r):
            S = set(S)
            cut = sum(1 for i, j in EDGES if (i in S) != (j in S))
            assert cut >= 2, f"cut {S} has size {cut}"
    print("  2-edge-connected: every cut has size >= 2 (Thm 1.8 applies)")

    # principal divisors have degree 0
    assert np.all(LAP.sum(axis=0) == 0) and np.all(LAP.sum(axis=1) == 0)
    print("  Laplacian rows/cols sum to 0: Prin(G) <= Div^0(G)")

    # ---- Kirchhoff: |Jac| = # spanning trees = N for a cycle ----
    RED = LAP[1:, 1:]
    ntrees = round(abs(np.linalg.det(RED)))
    print(f"  Matrix-Tree: det(reduced Laplacian) = {ntrees} spanning trees")
    assert ntrees == N, f"expected {N} spanning trees of C_{N}"
    SNF = smith_normal_form(Matrix(RED.tolist()))
    diag = [int(SNF[i, i]) for i in range(N - 1)]
    print(f"  Smith normal form of reduced Laplacian: diag{diag}")
    assert diag == [1, 1, 1, 1, N], f"SNF diag != [1,1,1,1,{N}]"
    print(f"  => Jac(C_{N}) ~= Z/{N}Z  (cyclic of order {N})")

    # confluence: stabilization independent of subset visit order
    # (uniqueness is BN Prop 3.1; this checks the implementation)
    rng = np.random.default_rng(8)
    for trial in range(30):
        D = rng.integers(-4, 5, size=N)
        D[Q] -= D.sum()  # force degree 0
        r1 = reduced_class(D, "lex")
        r2 = reduced_class(D, "rev")
        assert np.array_equal(r1, r2), f"stabilization not confluent on {D}"
    print("  stabilization confluent on 30 random degree-0 divisors (uniqueness, BN Prop 3.1)")

    # regression: stable-but-not-q-reduced divisor (broke the naive fire/borrow loop)
    D_cex = np.array([-3, 1, 1, 1, 0, 0])
    assert all(D_cex[v] >= 0 for v in range(N) if v != Q)  # no debt away from q
    assert not is_q_reduced(D_cex), "counterexample unexpectedly q-reduced"
    assert np.array_equal(reduced_class(D_cex), np.zeros(N, dtype=int)), \
        "counterexample should stabilize to 0"
    print("  regression: stable-but-not-reduced [-3,1,1,1,0,0] stabilizes to 0")

    # all q-reduced degree-0 divisors: P2 with A={v} forces D(v) < deg(v) = 2
    qreds = []
    for vals in itertools.product([0, 1], repeat=N - 1):
        D = np.zeros(N, dtype=int)
        D[1:] = vals
        D[Q] = -D.sum()
        if is_q_reduced(D):
            qreds.append(D)
    print(f"  q-reduced degree-0 divisors (G-parking functions): {len(qreds)}")
    assert len(qreds) == N, "expected #spanning trees = 6 q-reduced divisors"

    # ---- Abel-Jacobi: S_q(v_i) = [(v_i) - (q)] ----
    classes = [aj(i) for i in range(N)]
    print("  Abel-Jacobi classes S_q(v_i) = [(v_i)-(q)] (q-reduced reps):")
    for i, c in enumerate(classes):
        print(f"    v_{i} -> {list(c)}")
    seen = {tuple(c) for c in classes}
    assert len(seen) == N, "AJ map not injective!"
    print(f"  {N} distinct classes = |Jac| -> S^(1): V -> Jac is BIJECTIVE")
    print("  (BN Thm 1.7: S^(g) surjective; Thm 1.8: S^(1) injective iff 2-edge-connected)")

    # ---- the group law: it is Z/6Z ----
    order = 1
    acc = classes[1]
    while not np.array_equal(acc, classes[0]):
        acc = add_classes(acc, classes[1])
        order += 1
        assert order <= N
    print(f"  generator g = S_q(v_1) has order {order}")
    assert order == N
    for i in range(N):
        for j in range(N):
            s = add_classes(classes[i], classes[j])
            assert np.array_equal(s, classes[(i + j) % N]), \
                f"group law fails: {i}+{j}"
    print("  S_q(v_i) + S_q(v_j) = S_q(v_{i+j mod 6}): the cycle IS its Jacobian,")
    print("  addition = rotation. Tropical shadow of the elliptic-curve group law.")

    # ---- Riemann-Roch ----
    K = np.array([LAP[v, v] - 2 for v in range(N)])  # BN (1.11)
    print(f"  canonical divisor K = {[int(x) for x in K]} (deg {K.sum()} = 2g-2 = 0)")
    assert K.sum() == 2 * G - 2

    # hand checks first
    assert rank(np.zeros(N, dtype=int)) == 0, "r(0) should be 0"
    e0 = np.zeros(N, dtype=int)
    e0[0] = 1
    assert rank(e0) == 0, "r((v0)) should be 0"
    e0d = np.zeros(N, dtype=int)
    e0d[0] = 2
    assert rank(e0d) == 1, "r(2(v0)) should be 1 (RR: 1 - (-1) = 2)"
    assert rank(K) == G - 1 == 0, "r(K) = g-1 = 0 (BN Remark 2.4)"
    print("  hand checks: r(0)=0, r((v0))=0, r(2(v0))=1, r(K)=g-1=0")

    # brute force over {-1,0,1}^6: 729 divisors, K = 0 so RR is r(D) - r(-D) = deg D
    n_checked = 0
    for vals in itertools.product([-1, 0, 1], repeat=N):
        D = np.array(vals, dtype=int)
        lhs = rank(D) - rank(-D)
        rhs = int(D.sum()) + 1 - G
        assert lhs == rhs, f"RR fails for {vals}: {lhs} != {rhs}"
        n_checked += 1
    print(f"  Riemann-Roch r(D) - r(K-D) = deg(D) + 1 - g verified on {n_checked} divisors")
    print("  (K = 0 on the cycle, so this reads r(D) - r(-D) = deg D)")

    print("\nALL CHAPTER-12 CORE ASSERTIONS PASSED")


def phase2_checks():
    print("Phase-2 verification groups (Try It, referee checks, metric coda):")

    # ---------- A. Try It answers ----------
    classes = [aj(i) for i in range(N)]
    # Try It 1: identity is v0; S(v3)+S(v4) = S(v1)
    assert np.array_equal(classes[0], np.zeros(N, dtype=int))
    assert np.array_equal(add_classes(classes[3], classes[4]), classes[1])
    print("  A1 Try It 1: identity v0; S(v3)+S(v4)=S(v1)  OK")
    # Try It 2 is B below.
    # Try It 3: C_n has n spanning trees (det of reduced Laplacian), n=3..8;
    # group structure at n=6 comes from the SNF (asserted in the main program).
    for n in range(3, 9):
        L = np.zeros((n, n), dtype=int)
        for i in range(n):
            L[i, i] = 2
            L[i, (i + 1) % n] -= 1
            L[(i + 1) % n, i] -= 1
        det = round(abs(np.linalg.det(L[1:, 1:])))
        assert det == n, f"C_{n}: {det} trees != {n}"
    print("  A3 Try It 3: #trees(C_n) = n for n=3..8, so |Jac(C_n)| = n  OK")
    # Try It 4: D = (v0)+(v3); RR + r(-D) = -1 (negative degree) gives r(D) = 1,
    # and the program's rank agrees.
    D4 = np.zeros(N, dtype=int); D4[0] = 1; D4[3] = 1
    assert int(D4.sum()) == 2 and rank(-D4) == -1 and rank(D4) == 1
    assert rank(D4) - rank(-D4) == int(D4.sum())  # RR with K=0, g=1
    print("  A4 Try It 4: r((v0)+(v3)) = 1 (RR: 1-(-1) = deg 2)  OK")

    # ---------- B. violating sets + figure trace ----------
    Dcex = np.array([-3, 1, 1, 1, 0, 0])
    viol = []
    for r in range(1, N):
        for A in itertools.combinations(range(1, N), r):
            if all(Dcex[v] >= outdeg(A, v) for v in A):
                viol.append(tuple(A))
    assert viol == [(1, 2), (2, 3), (1, 2, 3)], viol
    # the draft's named set {1,2,3}: out-degrees 1,0,1, D = 1,1,1 there
    assert [outdeg((1, 2, 3), v) for v in (1, 2, 3)] == [1, 0, 1]
    _, fired = stabilize(Dcex, _trace=True)
    assert [sorted(s) for s, _ in fired] == [[1, 2], [3], [2, 3, 4], [1, 2, 3, 4, 5]]
    print("  B  violating sets exactly {1,2},{2,3},{1,2,3}; Fig 12.1B trace OK")

    # ---------- C. Dhar burning <=> q-reduced ----------
    def burns_all(D):
        if any(D[v] < 0 for v in range(N) if v != 0):
            return False
        burned = {0}
        changed = True
        while changed:
            changed = False
            for v in range(N):
                if v in burned:
                    continue
                e = sum(1 for i, j in EDGES
                        if (i == v and j in burned) or (j == v and i in burned))
                if D[v] < e:
                    burned.add(v)
                    changed = True
        return len(burned) == N

    count = 0
    for vals in itertools.product([-1, 0, 1, 2], repeat=N - 1):
        D = np.zeros(N, dtype=int)
        D[1:] = vals
        D[0] = -sum(vals)
        assert burns_all(D) == is_q_reduced(D), D
        count += 1
    assert not burns_all(Dcex)
    print(f"  C  Dhar burn <=> P1+P2 on all {count} test divisors  OK")

    # ---------- D. independent order count via surjectivity ----------
    qreds = set()
    for vals in itertools.product([0, 1], repeat=N - 1):
        D = np.zeros(N, dtype=int)
        D[1:] = vals
        D[0] = -D.sum()
        if is_q_reduced(D):
            qreds.add(tuple(D))
    aj_image = {tuple(c) for c in classes}
    assert aj_image == qreds and len(qreds) == 6
    print("  D  AJ image = all 6 q-reduced deg-0 divisors (S^(1) surjective,")
    print("     Thm 1.7 at k=g=1): |Jac| = 6 without Matrix-Tree  OK")

    # ---------- E. what the 729-box tests ----------
    n_pos = n_zero = 0
    for vals in itertools.product([-1, 0, 1], repeat=N):
        D = np.array(vals, dtype=int)
        d = int(D.sum())
        if d >= 1:
            assert rank(D) == d - 1, (vals, rank(D))
            n_pos += 1
        elif d == 0:
            assert (rank(D) == 0) == (rank(-D) == 0)
            assert rank(D) == rank(-D)
            n_zero += 1
    print(f"  E  box: r(D)=deg-1 on {n_pos} divisors with deg>=1; "
          f"D~-D symmetry on {n_zero} of degree 0  OK")

    # ---------- F. lattice length (referee point 1) ----------
    # Day 2 cycle vertices in cyclic order (face -> position, day02 log):
    cyc = [(-3, -2), (-2, -2), (-2, -3), (-3, -4), (-4, -4), (-4, -3)]
    eucl, latt, cum = [], [], [Fraction(0)]
    for i in range(6):
        x1, y1 = cyc[i]
        x2, y2 = cyc[(i + 1) % 6]
        dx, dy = x2 - x1, y2 - y1
        eucl.append(math.hypot(dx, dy))
        g = gcd(abs(dx), abs(dy))
        prim = (dx // g, dy // g)
        # displacement = g * primitive direction -> lattice length g
        assert (g * prim[0], g * prim[1]) == (dx, dy)
        latt.append(g)
        cum.append(cum[-1] + g)
    assert sorted(eucl) == sorted([1, 1, 1, 1, math.sqrt(2), math.sqrt(2)])
    assert abs(sum(eucl) - (4 + 2 * math.sqrt(2))) < 1e-12
    assert latt == [1] * 6 and sum(latt) == 6
    assert cum[:-1] == [Fraction(i) for i in range(6)]  # vertices at integers
    # Euclidean cumulative positions are NOT equally spaced:
    cum_e = np.cumsum([0] + eucl[:-1])
    assert not np.allclose(np.diff(cum_e), np.diff(cum_e)[0])
    print("  F  Euclidean total 4+2*sqrt(2) but lattice lengths all 1,")
    print("     L = 6, vertices at lattice arclengths 0..5  OK")

    # ---------- G. L(D) closed under min ----------
    def in_L(f, D):
        return bool(np.all(LAP @ np.array(f) >= -np.array(D)))

    for Dtest in ([1, 0, 0, 1, 0, 0], [0] * 6, [2, 0, 0, 0, 0, 0]):
        fs = [f for f in itertools.product(range(-2, 3), repeat=N)
              if in_L(f, Dtest)]
        assert fs
        for f in fs:
            for g_ in fs:
                assert in_L(np.minimum(f, g_), Dtest)
        print(f"  G  L(D) min-closed in box for D={Dtest} ({len(fs)} elements)  OK")

    print("\nALL CHAPTER-12 PHASE-2 ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
    phase2_checks()
    print("\nALL CHAPTER-12 ASSERTIONS PASSED")
