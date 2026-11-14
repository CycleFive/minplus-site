"""Day 8: chip-firing on the tropical cubic's 6-cycle.

The day-2 tropical cubic's bounded 6-cycle is a tropical elliptic curve.
Here it gets the Baker-Norine treatment (arXiv:math/0608360, verified
against the ar5iv HTML this session):
  - divisors, Laplacian, principal divisors, Jac(G) = Div^0/Prin,
  - q-reduced divisors (G-parking functions) -- unique per class,
  - the Abel-Jacobi map S_q(v) = [(v)-(q)], proved bijective for the
    2-edge-connected genus-1 cycle,
  - the group law on the cycle (tropical shadow of the elliptic-curve
    group law),
  - tropical Riemann-Roch r(D) - r(K-D) = deg(D) + 1 - g, brute-forced
    over all divisors with entries in {-1,0,1}.

Run as a script: every claim is asserted in code.
Run log: hidden_files/day08-run.log.
Importable: stabilize/reduced_class/rank/aj are reusable machinery.
"""

import itertools

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
    print(f"  canonical divisor K = {list(K)} (deg {K.sum()} = 2g-2 = 0)")
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

    print("\nALL DAY-8 ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
