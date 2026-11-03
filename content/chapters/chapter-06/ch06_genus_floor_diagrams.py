"""
Chapter 6: Genus Is Cycles -- floor diagrams of genus g > 0.

The full program behind Chapter 6 of "Learn You a Tropical Geometry
for Great Good!" Converted from Day 5 of the tropical geometry study
(2026-10-15). SELF-CONTAINED: the relative floor-diagram machinery of
Chapter 5 (itself converted from Day 4, definitions verified against
T. Dedieu, arXiv:2511.21123, Lecture IX) is inlined below, followed by
the genus extension. No imports from other study files.

What changes at genus g (Dedieu (1.1)-(1.2) are stated for general g):
  - a degree-d genus-g floor diagram has d floors, div = in - out = 1
    at every floor, no directed cycles, and first Betti number
        b1 = |E_fin| - d + 1 = g,   i.e. |E_fin| = d - 1 + g;
  - the underlying multigraph may have PARALLEL edges -- a doubled
    elevator is already a cycle -- so the enumeration runs over
    multisets of ordered pairs (combinations_with_replacement), not
    sets. (Reusing Chapter 5's set enumeration finds ZERO genus-1
    diagrams: on 3 floors the only 3-edge acyclic set is the
    transitive tournament, whose sink always violates div = 1.)
  - the automorphism count must preserve per-pair weight MULTISETS
    and include parallel-edge swap factors prod m_{i,j,w}! ; the
    action on labellings stays free, so #markings = LE / |Aut|
    exactly, as before.
  - multiplicity is unchanged: mu^C = I^beta * prod_{fin} w^2
    (Dedieu (1.1)); the absolute counts below are the relative type
    (0, (1^d)), so I^beta = 1 and alpha is empty.

Chapter claims asserted in __main__ (curves of degree d and genus g
through 3d - 1 + g points):
  N_{3,1} = 1    (1 shape: doubled 2->0 plus 0->1, mu=1, 1 marking;
                  classically the unique cubic through 9 points is
                  smooth, genus 1)
  N_{4,1} = 225  (11 shapes; per-shape contributions asserted;
                  classical witness: genus-1 quartics are 2-nodal
                  quartics, node polynomial
                  N2(d) = (3/2)(d-1)(d-2)(3d^2-3d-11)
                  (Vainsencher; Kleiman-Piene for higher delta)
                  gives 225 at d=4 and cross-checks at d=3 with
                  21 = C(7,2), the reducible line+conic count)
  N_{4,2} = 27   (5 shapes, contributions 3+4+5+7+8;
                  classical witness: 1-nodal quartics,
                  N1(d) = 3(d-1)^2 (Steiner) = 27 at d=4)
  N_{4,3} = 1    (1 shape: doubled 0->1, single 1->2, tripled 3->0;
                  the unique smooth quartic through 14 points)
  b1 = g is asserted on every enumerated shape, and both d=4 totals
  at genus 1 and 2 plus N_{3,1} are re-run with the finite-edge
  weight bound raised d -> d+1 with identical results (the bound is
  exhaustive). N_{4,3} at bound d+1 is not re-run: its bound-d
  enumeration already takes ~1 minute and the exhaustive-bound
  argument is the same upward-closed one as Chapters 4-5; the three
  re-runs that are performed cover both degrees and three genera.
"""


from itertools import combinations, product
from collections import Counter
from functools import lru_cache
import math


# --------------------------------------------------------------------------
# Diagram representation
# --------------------------------------------------------------------------
# ainf[v]: tuple of weights of alpha-infinite edges at floor v (fixed)
# binf[v]: tuple of weights of beta-infinite edges at floor v (moving)
# fin: tuple of (i, j, w), i != j


def check_relative(d, alpha, beta, ainf, binf, fin, genus=0):
    I = lambda prof: sum(i * c for i, c in prof.items())
    if I(alpha) + I(beta) != d:
        return False
    # weight profiles match
    if Counter(w for v in range(d) for w in ainf[v]) != Counter(
            {i: c for i, c in alpha.items() for _ in range(c)}):
        return False
    if Counter(w for v in range(d) for w in binf[v]) != Counter(
            {i: c for i, c in beta.items() for _ in range(c)}):
        return False
    # divergence = 1
    for v in range(d):
        inw = (sum(ainf[v]) + sum(binf[v])
               + sum(w for (i, j, w) in fin if j == v))
        outw = sum(w for (i, j, w) in fin if i == v)
        if inw - outw != 1:
            return False
    if not is_acyclic(d, fin):
        return False
    if len(fin) != d - 1 + genus:
        return False
    if not is_connected(d, ainf, binf, fin):
        return False
    return True


def is_acyclic(d, fin):
    indeg = [0] * d
    adj = [[] for _ in range(d)]
    for (i, j, w) in fin:
        adj[i].append(j)
        indeg[j] += 1
    stack = [v for v in range(d) if indeg[v] == 0]
    seen = 0
    while stack:
        v = stack.pop()
        seen += 1
        for u in adj[v]:
            indeg[u] -= 1
            if indeg[u] == 0:
                stack.append(u)
    return seen == d


def is_connected(d, ainf, binf, fin):
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for v in range(d):
        parent[v] = v
    nxt = d
    for v in range(d):
        for _ in ainf[v]:
            parent[nxt] = nxt
            union(v, nxt)
            nxt += 1
        for _ in binf[v]:
            parent[nxt] = nxt
            union(v, nxt)
            nxt += 1
    for (i, j, w) in fin:
        union(i, j)
    return len({find(v) for v in range(nxt)}) == 1


def canonical_key(d, ainf, binf, fin):
    from itertools import permutations
    best = None
    for p in permutations(range(d)):
        # p: old floor v -> new label p[v]; inverse: new label i <- old inv[i]
        inv = [0] * d
        for v in range(d):
            inv[p[v]] = v
        a2 = tuple(tuple(sorted(ainf[inv[i]])) for i in range(d))
        b2 = tuple(tuple(sorted(binf[inv[i]])) for i in range(d))
        fin2 = tuple(sorted((p[i], p[j], w) for (i, j, w) in fin))
        key = (a2, b2, fin2)
        if best is None or key < best:
            best = key
    return best


def multisets_summing_to(n, kinds):
    """All multisets (as lists of (kind, weight)) of total weight n.

    kinds: list of (kind, weight) with distinct weights per kind irrelevant;
    enumeration by multiplicity per kind-index is unique.
    """
    out = []

    def rec(k, rem, cur):
        if k == len(kinds):
            if rem == 0:
                out.append(list(cur))
            return
        kind, w = kinds[k]
        for c in range(rem // w + 1):
            cur.extend([(kind, w)] * c)
            rec(k + 1, rem - c * w, cur)
            del cur[len(cur) - c:]
    rec(0, n, [])
    return out


def inf_distributions(d, alpha, beta, need):
    """Ways to realize per-floor infinite-weight needs with the (alpha,beta)
    profile. Yields (ainf, binf), each a tuple of d weight-tuples."""
    target = {}
    for i, c in alpha.items():
        target[('a', i)] = c
    for i, c in beta.items():
        target[('b', i)] = c
    kinds = sorted(target.keys(), key=lambda k: (k[0], k[1]))
    kw = [(k, k[1]) for k in kinds]  # weight of kind k is k[1]
    results = []

    def rec(v, remaining, ca, cb):
        if v == d:
            if all(x == 0 for x in remaining.values()):
                results.append((tuple(ca), tuple(cb)))
            return
        for ms in multisets_summing_to(need[v], kw):
            use = Counter(k for (k, w) in ms)
            if any(use[k] > remaining.get(k, 0) for k in use):
                continue
            new_rem = {k: remaining[k] - use.get(k, 0) for k in remaining}
            a = tuple(sorted(w for (k, w) in ms if k[0] == 'a'))
            b = tuple(sorted(w for (k, w) in ms if k[0] == 'b'))
            ca.append(a)
            cb.append(b)
            rec(v + 1, new_rem, ca, cb)
            ca.pop()
            cb.pop()
    rec(0, dict(target), [], [])
    return results


def relative_diagrams(d, alpha, beta, genus=0, wbound=None):
    """Enumerate relative floor diagrams of type (alpha,beta), degree d."""
    I = sum(i * c for i, c in alpha.items()) + sum(i * c for i, c in beta.items())
    assert I == d, (alpha, beta, d)
    F = d - 1 + genus
    if wbound is None:
        wbound = d
    pairs = [(i, j) for i in range(d) for j in range(d) if i != j]
    seen = set()
    for epairs in combinations(pairs, F):
        if not is_acyclic(d, [(i, j, 1) for (i, j) in epairs]):
            continue
        for weights in product(range(1, wbound + 1), repeat=F):
            fin = tuple(sorted((i, j, w) for (i, j), w in zip(epairs, weights)))
            need = []
            ok = True
            for v in range(d):
                outw = sum(w for (i, j, w) in fin if i == v)
                inw = sum(w for (i, j, w) in fin if j == v)
                c = 1 + outw - inw
                if c < 0:
                    ok = False
                    break
                need.append(c)
            if not ok or sum(need) != d:
                continue
            for ainf, binf in inf_distributions(d, alpha, beta, need):
                if not is_connected(d, ainf, binf, fin):
                    continue
                key = (ainf, binf, fin)
                if key in seen:
                    continue
                seen.add(key)
                yield ainf, binf, fin


def relative_shapes(d, alpha, beta, genus=0, **kw):
    seen = set()
    for ainf, binf, fin in relative_diagrams(d, alpha, beta, genus, **kw):
        key = canonical_key(d, ainf, binf, fin)
        if key in seen:
            continue
        seen.add(key)
        yield key


# --------------------------------------------------------------------------
# Markings: poset of floors + finite edges + beta-infinite edges (Dedieu 1.6)
# --------------------------------------------------------------------------

def poset_elements(d, ainf, binf, fin):
    elems = []
    index = {}
    for v in range(d):
        index[('floor', v)] = len(elems)
        elems.append(('floor', v))
    for t, (i, j, w) in enumerate(fin):
        index[('fin', t)] = len(elems)
        elems.append(('fin', t))
    for v in range(d):
        for k, w in enumerate(binf[v]):
            index[('binf', v, k)] = len(elems)
            elems.append(('binf', v, k, w))
    n = len(elems)
    le = [[False] * n for _ in range(n)]
    for t, (i, j, w) in enumerate(fin):
        le[index[('floor', i)]][index[('fin', t)]] = True
        le[index[('fin', t)]][index[('floor', j)]] = True
    for v in range(d):
        for k, w in enumerate(binf[v]):
            le[index[('binf', v, k)]][index[('floor', v)]] = True
    for m in range(n):
        for a in range(n):
            if le[a][m]:
                for b in range(n):
                    if le[m][b]:
                        le[a][b] = True
    pred = []
    for b in range(n):
        mask = 0
        for a in range(n):
            if le[a][b]:
                mask |= (1 << a)
        pred.append(mask)
    return n, pred


def count_linear_extensions(n, pred):
    N = 1 << n
    dp = [0] * N
    dp[0] = 1
    for mask in range(N):
        c = dp[mask]
        if not c:
            continue
        for j in range(n):
            if mask & (1 << j):
                continue
            if pred[j] & ~mask == 0:
                dp[mask | (1 << j)] += c
    return dp[N - 1]


def automorphism_count(d, ainf, binf, fin):
    from itertools import permutations
    wmat = [[0] * d for _ in range(d)]
    for (i, j, w) in fin:
        wmat[i][j] = w
    count = 0
    for p in permutations(range(d)):
        if any(tuple(sorted(ainf[v])) != tuple(sorted(ainf[p[v]])) for v in range(d)):
            continue
        if any(tuple(sorted(binf[v])) != tuple(sorted(binf[p[v]])) for v in range(d)):
            continue
        if any(wmat[i][j] != wmat[p[i]][p[j]] for i in range(d) for j in range(d)):
            continue
        count += 1
    for v in range(d):
        ca = Counter(ainf[v])
        cb = Counter(binf[v])
        for c in ca.values():
            count *= math.factorial(c)
        for c in cb.values():
            count *= math.factorial(c)
    return count


def markings(d, alpha, beta, ainf, binf, fin):
    """#markings(D) per Dedieu (1.6).

    Labels {1..n} (n = 2d+g-1+|beta|) mark floors, finite edges, beta-infinite
    edges (mobile); labels {-|alpha|+1..0} mark the alpha-infinite edges
    (fixed), the alpha_k labels of weight-class k permuting freely among the
    alpha_k weight-k alpha-edges: factor prod_k (alpha_k)!.
    Two markings equivalent under oriented-graph automorphisms preserving m.
    """
    n, pred = poset_elements(d, ainf, binf, fin)
    le = count_linear_extensions(n, pred)
    aut = automorphism_count(d, ainf, binf, fin)
    assert le % aut == 0, (ainf, binf, fin, le, aut)
    alpha_label_factor = 1
    for c in alpha.values():
        alpha_label_factor *= math.factorial(c)
    return le // aut * alpha_label_factor


# --------------------------------------------------------------------------
# Multiplicity: mu^C(D) = I^beta * prod_{fin} w^2 (Dedieu (1.1) Definition).
# The `mode` switch keeps the pre-verification alternatives for the record;
# "guess" is the verified formula.
# --------------------------------------------------------------------------

def multiplicity(fin, binf, mode="guess"):
    """mode 'guess': prod_{fin} w^2 * prod_{beta-inf} w."""
    m = 1
    for (_, _, w) in fin:
        m *= w * w
    if mode == "guess":
        for v in binf:
            for w in v:
                m *= w
    elif mode == "fin-only":
        pass
    elif mode == "fin-and-beta-sq":
        for v in binf:
            for w in v:
                m *= w * w
    else:
        raise ValueError(mode)
    return m


def relative_count(d, alpha, beta, genus=0, mode="guess", verbose=True, **kw):
    total = 0
    n_diag = 0
    for ainf, binf, fin in relative_shapes(d, alpha, beta, genus, **kw):
        mu = multiplicity(fin, binf, mode)
        mk = markings(d, alpha, beta, ainf, binf, fin)
        contrib = mu * mk
        total += contrib
        n_diag += 1
        if verbose:
            print(f"  ainf={ainf} binf={binf} fin={fin} mu={mu} "
                  f"markings={mk} contrib={contrib}")
    if verbose:
        print(f"  {n_diag} shapes, N = {total}")
    return total


def breakdown(d, alpha, beta, genus=0, **kw):
    """Sorted (mu, markings, contribution) per canonical shape."""
    out = []
    for ainf, binf, fin in relative_shapes(d, alpha, beta, genus, **kw):
        assert check_relative(d, alpha, beta, ainf, binf, fin, genus)
        mu = multiplicity(fin, binf)
        mk = markings(d, alpha, beta, ainf, binf, fin)
        out.append((mu, mk, mu * mk))
    return sorted(out)


# ------------------------------------------------------------------------
# Genus extension (Day 5): parallel edges, multiset enumeration,
# multiset-preserving automorphisms.
# ------------------------------------------------------------------------

from itertools import combinations_with_replacement, permutations


def automorphism_count_genus(d, ainf, binf, fin):
    """|Aut(D)| for possibly non-tree diagrams.

    Aut = (floor perms preserving the weighted directed multigraph and the
    per-floor inf profiles) x (parallel finite-edge swaps) x (inf-edge
    swaps).  The action on labelings is free (labels distinct), so
    #markings = LE / |Aut| exactly.
    """
    # per ordered pair: sorted tuple of finite-edge weights
    wmulti = [[() for _ in range(d)] for _ in range(d)]
    tmp = [[[] for _ in range(d)] for _ in range(d)]
    for (i, j, w) in fin:
        tmp[i][j].append(w)
    for i in range(d):
        for j in range(d):
            wmulti[i][j] = tuple(sorted(tmp[i][j]))
    n_floor_perm = 0
    for p in permutations(range(d)):
        if any(tuple(sorted(ainf[v])) != tuple(sorted(ainf[p[v]]))
               for v in range(d)):
            continue
        if any(tuple(sorted(binf[v])) != tuple(sorted(binf[p[v]]))
               for v in range(d)):
            continue
        if any(wmulti[i][j] != wmulti[p[i]][p[j]]
               for i in range(d) for j in range(d)):
            continue
        n_floor_perm += 1
    aut = n_floor_perm
    # parallel finite edges (i,j) of equal weight are indistinguishable
    for i in range(d):
        for j in range(d):
            for c in Counter(tmp[i][j]).values():
                aut *= math.factorial(c)
    # infinite edges: indistinguishable within (floor, weight, a/b)
    for v in range(d):
        for c in Counter(ainf[v]).values():
            aut *= math.factorial(c)
        for c in Counter(binf[v]).values():
            aut *= math.factorial(c)
    return aut


def markings_genus(d, alpha, beta, ainf, binf, fin):
    n, pred = poset_elements(d, ainf, binf, fin)
    le = count_linear_extensions(n, pred)
    aut = automorphism_count_genus(d, ainf, binf, fin)
    assert le % aut == 0, (ainf, binf, fin, le, aut)
    alpha_label_factor = 1
    for c in alpha.values():
        alpha_label_factor *= math.factorial(c)
    return le // aut * alpha_label_factor


def relative_diagrams_genus(d, alpha, beta, genus, wbound=None):
    I = (sum(i * c for i, c in alpha.items())
         + sum(i * c for i, c in beta.items()))
    assert I == d, (alpha, beta, d)
    F = d - 1 + genus
    if wbound is None:
        wbound = d
    pairs = [(i, j) for i in range(d) for j in range(d) if i != j]
    seen = set()
    n_pair_sets = 0
    for epairs in combinations_with_replacement(pairs, F):
        if not is_acyclic(d, [(i, j, 1) for (i, j) in epairs]):
            continue
        n_pair_sets += 1
        for weights in product(range(1, wbound + 1), repeat=F):
            fin = tuple(sorted((i, j, w)
                               for (i, j), w in zip(epairs, weights)))
            need = []
            ok = True
            for v in range(d):
                outw = sum(w for (i, j, w) in fin if i == v)
                inw = sum(w for (i, j, w) in fin if j == v)
                c = 1 + outw - inw
                if c < 0:
                    ok = False
                    break
                need.append(c)
            if not ok or sum(need) != d:
                continue
            for ainf, binf in inf_distributions(d, alpha, beta, need):
                if not is_connected(d, ainf, binf, fin):
                    continue
                key = (ainf, binf, fin)
                if key in seen:
                    continue
                seen.add(key)
                yield ainf, binf, fin
    # print(f"  [debug] {n_pair_sets} acyclic pair-multisets", file=sys.stderr)


def relative_shapes_genus(d, alpha, beta, genus, **kw):
    seen = set()
    for ainf, binf, fin in relative_diagrams_genus(d, alpha, beta, genus,
                                                   **kw):
        key = canonical_key(d, ainf, binf, fin)
        if key in seen:
            continue
        seen.add(key)
        yield key


def relative_count_genus(d, alpha, beta, genus, verbose=True, **kw):
    total = 0
    n_diag = 0
    shapes = []
    for ainf, binf, fin in relative_shapes_genus(d, alpha, beta, genus,
                                                 **kw):
        assert check_relative(d, alpha, beta, ainf, binf, fin, genus=genus)
        # b1 check: E - V + 1 on the finite part (+inf stubs are trees)
        assert len(fin) - d + 1 == genus
        mu = multiplicity(fin, binf, "guess")
        mk = markings_genus(d, alpha, beta, ainf, binf, fin)
        contrib = mu * mk
        total += contrib
        n_diag += 1
        shapes.append((ainf, binf, fin, mu, mk, contrib))
        if verbose:
            print(f"  ainf={ainf} binf={binf} fin={fin} mu={mu} "
                  f"markings={mk} contrib={contrib}")
    if verbose:
        print(f"  {n_diag} shapes, N = {total}")
    return total, shapes


def breakdown_genus(d, genus, **kw):
    """Sorted (mu, markings, contribution) per canonical shape."""
    out = []
    for ainf, binf, fin in relative_shapes_genus(d, {}, {1: d}, genus, **kw):
        assert check_relative(d, {}, {1: d}, ainf, binf, fin, genus)
        assert len(fin) - d + 1 == genus
        mu = multiplicity(fin, binf)
        mk = markings_genus(d, {}, {1: d}, ainf, binf, fin)
        out.append((mu, mk, mu * mk))
    return sorted(out)


def node_poly_N2(d):
    """Vainsencher's 2-node polynomial (3/2)(d-1)(d-2)(3d^2-3d-11)."""
    num = 3 * (d - 1) * (d - 2) * (3 * d * d - 3 * d - 11)
    assert num % 2 == 0
    return num // 2


if __name__ == "__main__":
    import time

    # classical witnesses, checked as formulas before any enumeration
    assert node_poly_N2(4) == 225
    assert node_poly_N2(3) == 21 == math.comb(7, 2)
    assert 3 * (4 - 1) ** 2 == 27
    print("classical witnesses OK: N2(4)=225, N2(3)=21=C(7,2), "
          "N1(4)=3(d-1)^2=27")

    t0 = time.time()
    b31 = breakdown_genus(3, 1)
    assert b31 == [(1, 1, 1)], b31
    print(f"N_{{3,1}} = 1 OK, per-shape: {b31} "
          f"({time.time()-t0:.1f}s)")

    t0 = time.time()
    b41 = breakdown_genus(4, 1)
    assert b41 == [(1, 9, 9), (1, 15, 15), (1, 21, 21), (1, 21, 21),
                   (1, 26, 26), (4, 4, 16), (4, 6, 24), (4, 6, 24),
                   (4, 7, 28), (9, 1, 9), (16, 2, 32)], b41
    assert sum(c for _, _, c in b41) == 225
    print(f"N_{{4,1}} = 225 OK, per-shape: {b41} "
          f"({time.time()-t0:.1f}s)")

    t0 = time.time()
    b42 = breakdown_genus(4, 2)
    assert b42 == [(1, 3, 3), (1, 5, 5), (1, 7, 7), (4, 1, 4),
                   (4, 2, 8)], b42
    assert sum(c for _, _, c in b42) == 27
    print(f"N_{{4,2}} = 27 OK, per-shape: {b42} "
          f"({time.time()-t0:.1f}s)")

    t0 = time.time()
    b43 = breakdown_genus(4, 3)
    assert b43 == [(1, 1, 1)], b43
    print(f"N_{{4,3}} = 1 OK, per-shape: {b43} "
          f"({time.time()-t0:.1f}s)")

    # robustness: weight bound d is exhaustive (re-run with d+1)
    assert breakdown_genus(3, 1, wbound=4) == b31
    assert breakdown_genus(4, 1, wbound=5) == b41
    assert breakdown_genus(4, 2, wbound=5) == b42
    print("robustness: wbound=d+1 changes nothing for "
          "(3,1), (4,1), (4,2) OK")

    print("Chapter 6: all assertions pass. "
          "N_{3,1}=1; N_{4,1}=225; N_{4,2}=27; N_{4,3}=1; "
          "b1 = g on every shape.")
