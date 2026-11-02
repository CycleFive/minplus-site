"""
Chapter 5: Tangency Is Just Floors Being Stubborn -- relative floor diagrams.

The full program behind Chapter 5 of "Learn You a Tropical Geometry
for Great Good!" Converted from Day 4 of the tropical geometry study
(2026-10-14), where the definitions were verified against the primary
source before any number was claimed: T. Dedieu, "Two lectures on the
enumeration of curves by means of floor diagrams" (arXiv:2511.21123),
Lecture IX, reporting Brugalle-Mikhalkin. Section citations below
refer to Lecture IX.

Definitions used (all confirmed against the source):
  - Diagram axioms (1.1)-(1.2): connected weighted oriented graph, no
    directed cycles, div(v) = incoming - outgoing = 1 at every floor,
    sum of infinite-edge weights = d (the paper prints "-d" in (1.2)(iii),
    a typo: weights are positive and the proof of Lemma (1.5) needs "=d";
    the same typo was already noted in Chapter 4).
  - Type (alpha,beta), I(alpha)+I(beta)=d: exactly alpha_i+beta_i infinite
    edges of weight i (1.4)-(1.5). Lemma (1.5)(i): still d floors.
  - Marking (1.6): labels {-|alpha|+1..0} on the alpha-infinite edges
    (fixed contacts; alpha_k labels of weight-class k permute freely),
    labels {1..n}, n = 2d+g-1+|beta|, on floors + finite edges +
    beta-infinite edges (mobile), strictly increasing for the incidence
    poset; counted up to oriented-graph automorphisms preserving m.
  - Multiplicity (1.1): mu^C(D) = I^beta * prod_{fin} w(e)^2 with
    I^beta = prod_i i^{beta_i} (beta-infinite edge of weight w -> w,
    alpha-infinite edge -> 1).
  - Correspondence (1.2) Thm: N_g^{d,irr}(alpha,beta) = sum_{(D,m)} mu^C(D)
    = # irreducible degree-d genus-g plane curves with alpha_i order-i
    contacts at fixed points of a line L, beta_i order-i contacts at moving
    points, through 2d+g-1+|beta| general points.

Chapter claims asserted in __main__ (Dedieu Ex. (1.3), per shape):
  sanity: alpha={}, beta={1:3} recovers the absolute count N = 12;
  N_0^{3,irr}(0,[1,1]) = 36 with (mu, markings, contribution) breakdown
    (2,6,12), (2,4,8), (8,2,16)  [alpha={}, beta={1:1,2:1}];
  N_0^{3,irr}([0,1],1) = 10 with breakdown
    (1,3,3), (1,3,3), (4,1,4)    [alpha={2:1}, beta={1:1}];
  and both totals are unchanged with the finite-edge weight bound
  raised from d to d+1 (the bound is exhaustive).

NOTE: Welschinger real counts are NOT in these lectures (full-text
check on Day 4); a different source is needed for W_3 = 8. This
program makes no real-count claim.
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


if __name__ == "__main__":
    # sanity: the relative machinery recovers Chapter 4's absolute count
    assert relative_count(3, {}, {1: 3}, verbose=False) == 12
    print("sanity: alpha={}, beta={1:3} -> N = 12 OK")

    b36 = breakdown(3, {}, {1: 1, 2: 1})
    assert b36 == [(2, 4, 8), (2, 6, 12), (8, 2, 16)], b36
    assert sum(c for _, _, c in b36) == 36
    print("N_0^3(0,[1,1]) = 36 OK, per-shape:", b36)

    b10 = breakdown(3, {2: 1}, {1: 1})
    assert b10 == [(1, 3, 3), (1, 3, 3), (4, 1, 4)], b10
    assert sum(c for _, _, c in b10) == 10
    print("N_0^3([0,1],1) = 10 OK, per-shape:", b10)

    # robustness: weight bound d is exhaustive (re-run with d+1)
    assert breakdown(3, {}, {1: 1, 2: 1}, wbound=4) == b36
    assert breakdown(3, {2: 1}, {1: 1}, wbound=4) == b10
    print("robustness: wbound=d+1 changes nothing OK")

    print("Chapter 5: all assertions pass. "
          "36 = 12 + 8 + 16; 10 = 3 + 3 + 4; absolute sanity N = 12.")
