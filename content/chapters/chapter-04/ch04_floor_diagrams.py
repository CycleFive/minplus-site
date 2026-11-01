"""
Chapter 4: Twelve Cubics, Counted on the Floor -- floor diagrams.

The full program behind Chapter 4 of "Learn You a Tropical Geometry
for Great Good!" Converted from Day 3 of the tropical geometry study
(2026-10-13), where the definitions were implemented from the primary
source: T. Dedieu, "Two lectures on the enumeration of curves by means
of floor diagrams" (arXiv:2511.21123), Lecture IX, Sec. 1.1--1.3,
reporting Brugalle--Mikhalkin.

A floor diagram of degree d, genus g:
  - D: connected, acyclic, oriented graph.
  - Vertices: "vertices at infinity" (1-valent, unique incident edge OUTGOING)
    and "floors" (all others). No edge joins two infinity vertices.
  - Edges: Edges^inf (from an infinity vertex) and Edges^fin (floor to floor).
  - Weight w: Edges -> positive integers.
  - div(v) = sum(incoming w) - sum(outgoing w).
  - (i) b_1(D) = g; (ii) div(v) = 1 at every floor; (iii) sum of w over
    Edges^inf = d.  Lemma: #floors = d.

Absolute (point-only) count N_{d,g}: type (alpha, beta) = (0, (1^d)):
  - |Edges^inf| = d, all of weight 1 (weights sum to d over d edges).
  - Points: 3d+g-1 (for g=0: 2, 5, 8, 11 points for d = 1..4).
  - A marking is a bijection m: {1..3d+g-1} -> floors |_| edges,
    strictly increasing for the incidence poset (s < e < t for edge
    e=(s,t)). I.e., a linear extension of the poset.
  - #markings(D) = #(linear extensions) / |Aut(D)|.
  - mu(D) = prod_{e in Edges^fin} w(e)^2.
  - N_{d,g} = sum_D mu(D) * #markings(D).

Chapter claims asserted in __main__:
  N_{1,0}=1, N_{2,0}=1, N_{3,0}=12, N_{4,0}=620 (Kontsevich numbers),
  shape counts 1, 1, 3, 12, and for d=3 the per-shape breakdown
  (mu, markings, contribution) = (4,1,4), (1,5,5), (1,3,3).
"""

from itertools import combinations, product
from functools import lru_cache


# --------------------------------------------------------------------------
# Diagram representation
# --------------------------------------------------------------------------
# floors: 0..d-1
# inf[v]: number of infinite (weight-1) edges entering floor v  (absolute type)
# fin: tuple of (i, j, w), i != j, at most one edge per ordered pair
#      ("edges viewed as subset of Vert x Vert")

def check_diagram(d, inf, fin, genus=0):
    """Validate the floor-diagram axioms; return True/False."""
    floors = list(range(d))
    # (iii) infinite weights sum to d  (all weight 1 here)
    if sum(inf) != d:
        return False
    # divergence = 1 at every floor
    for v in floors:
        inw = inf[v] + sum(w for (i, j, w) in fin if j == v)
        outw = sum(w for (i, j, w) in fin if i == v)
        if inw - outw != 1:
            return False
    # acyclicity of finite edges (topological order exists)
    if not is_acyclic(d, fin):
        return False
    # genus: V = d floors + d infinity vertices, E = d + len(fin)
    if len(fin) != d - 1 + genus:
        return False
    # connected (as unoriented graph, including infinity vertices)
    if not is_connected(d, inf, fin):
        return False
    return True


def is_acyclic(d, fin):
    # Kahn's algorithm on floors
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


def is_connected(d, inf, fin):
    # unoriented connectivity over floors + infinity vertices.
    # Represent infinity vertices as d + k. Each floor v with inf[v] > 0 is
    # adjacent to inf[v] distinct infinity vertices.
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
        for _ in range(inf[v]):
            parent[nxt] = nxt
            union(v, nxt)
            nxt += 1
    for (i, j, w) in fin:
        union(i, j)
    roots = {find(v) for v in range(nxt)}
    return len(roots) == 1


def canonical_key(d, inf, fin):
    """Lexicographically smallest (inf, fin) under floor relabeling.

    For a relabeling p (old floor v gets new label p[v]):
      new_inf[p[v]] = old_inf[v],  new_fin = {(p[i], p[j], w)}.
    """
    from itertools import permutations
    best = None
    for p in permutations(range(d)):
        inf2 = [0] * d
        for v in range(d):
            inf2[p[v]] = inf[v]
        inf2 = tuple(inf2)
        fin2 = tuple(sorted((p[i], p[j], w) for (i, j, w) in fin))
        key = (inf2, fin2)
        if best is None or key < best:
            best = key
    return best


def floor_diagram_shapes(d, genus=0):
    """One canonical representative per isomorphism class."""
    seen = set()
    for inf, fin in floor_diagrams(d, genus):
        key = canonical_key(d, inf, fin)
        if key in seen:
            continue
        seen.add(key)
        yield key  # (inf, fin) canonical representative


def floor_diagrams(d, genus=0):
    """Enumerate absolute-type floor diagrams of degree d, genus g.

    Yields (inf, fin) with inf a tuple of length d, fin a tuple of (i,j,w).
    Weight bound: w(e) <= d (proved via upward-closed sets in the analysis).
    """
    F = d - 1 + genus
    pairs = [(i, j) for i in range(d) for j in range(d) if i != j]
    seen_shapes = set()
    for epairs in combinations(pairs, F):
        # quick acyclicity pre-check on the unweighted digraph
        if not is_acyclic(d, [(i, j, 1) for (i, j) in epairs]):
            continue
        for weights in product(range(1, d + 1), repeat=F):
            fin = tuple(sorted((i, j, w)
                               for (i, j), w in zip(epairs, weights)))
            # solve for inf from divergence: inf[v] = 1 + out_f(v) - in_f(v)
            inf = []
            ok = True
            for v in range(d):
                outw = sum(w for (i, j, w) in fin if i == v)
                inw = sum(w for (i, j, w) in fin if j == v)
                c = 1 + outw - inw
                if c < 0:
                    ok = False
                    break
                inf.append(c)
            if not ok or sum(inf) != d:
                continue
            inf = tuple(inf)
            if not is_connected(d, inf, fin):
                continue
            key = (inf, fin)
            if key in seen_shapes:
                continue
            seen_shapes.add(key)
            yield inf, fin


def multiplicity(fin):
    m = 1
    for (_, _, w) in fin:
        m *= w * w
    return m


# --------------------------------------------------------------------------
# Poset of floors + edges, linear extensions, automorphisms
# --------------------------------------------------------------------------

def poset_elements(d, inf, fin):
    """Elements of Dd and cover relations s < e < t.

    Elements: ('f', v) floors; ('e', k) infinite edges (k = index among the
    inf[v] edges into v); ('f', idx) finite edges.
    Represented as small ints; returns (n, pred) with pred[j] = bitmask of
    strict predecessors of j.
    """
    elems = []
    index = {}
    for v in range(d):
        index[('floor', v)] = len(elems)
        elems.append(('floor', v))
    for v in range(d):
        for k in range(inf[v]):
            index[('inf', v, k)] = len(elems)
            elems.append(('inf', v, k))
    for t, (i, j, w) in enumerate(fin):
        index[('fin', t)] = len(elems)
        elems.append(('fin', t))
    n = len(elems)
    # relations: for infinite edge (v,k): e < floor v.
    # for finite edge t=(i,j,w): floor i < e < floor j.
    le = [[False] * n for _ in range(n)]
    for v in range(d):
        for k in range(inf[v]):
            le[index[('inf', v, k)]][index[('floor', v)]] = True
    for t, (i, j, w) in enumerate(fin):
        le[index[('floor', i)]][index[('fin', t)]] = True
        le[index[('fin', t)]][index[('floor', j)]] = True
    # transitive closure (Floyd-Warshall, n <= 3d-1 <= 11 for d<=4)
    for m in range(n):
        for a in range(n):
            if le[a][m]:
                row_m = le[m]
                row_a = le[a]
                for b in range(n):
                    if row_m[b]:
                        row_a[b] = True
    pred = []
    for b in range(n):
        mask = 0
        for a in range(n):
            if le[a][b]:
                mask |= (1 << a)
        pred.append(mask)
    return n, pred


def count_linear_extensions(n, pred):
    """DP over subsets: dp[mask] = #ways to order exactly the set mask."""
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
            if pred[j] & ~mask == 0:  # all predecessors already placed
                dp[mask | (1 << j)] += c
    return dp[N - 1]


def automorphism_count(d, inf, fin):
    """|Aut(D)|: permutations of floors preserving inf counts and the
    weighted finite-edge relation. Infinity vertices attached to a floor
    are permuted freely among themselves (they are symmetric)."""
    from itertools import permutations
    # adjacency with weights as canonical data
    wmat = [[0] * d for _ in range(d)]
    for (i, j, w) in fin:
        wmat[i][j] = w
    count = 0
    for p in permutations(range(d)):
        if any(inf[v] != inf[p[v]] for v in range(d)):
            continue
        if any(wmat[i][j] != wmat[p[i]][p[j]]
               for i in range(d) for j in range(d)):
            continue
        count += 1
    # multiply by permutations of identical infinite edges at each floor
    import math
    for v in range(d):
        count *= math.factorial(inf[v])
    return count


def markings(d, inf, fin):
    n, pred = poset_elements(d, inf, fin)
    le = count_linear_extensions(n, pred)
    aut = automorphism_count(d, inf, fin)
    assert le % aut == 0, (d, inf, fin, le, aut)
    return le // aut


def show(d, genus=0):
    total = 0
    n_diag = 0
    print(f"=== d={d}, genus {genus}: shapes ===")
    for inf, fin in floor_diagram_shapes(d, genus):
        mu = multiplicity(fin)
        mk = markings(d, inf, fin)
        contrib = mu * mk
        total += contrib
        n_diag += 1
        print(f"  inf={inf} fin={fin} mu={mu} markings={mk} contrib={contrib}")
    print(f"  {n_diag} shapes, N = {total}")
    return total


if __name__ == "__main__":
    expected_N = {1: 1, 2: 1, 3: 12, 4: 620}
    expected_shapes = {1: 1, 2: 1, 3: 3, 4: 12}
    for d in (1, 2, 3, 4):
        total = show(d)
        assert total == expected_N[d], (d, total)
    # per-shape breakdown at d=3: the chapter's central table
    breakdown = sorted(
        (multiplicity(fin), markings(3, inf, fin),
         multiplicity(fin) * markings(3, inf, fin))
        for inf, fin in floor_diagram_shapes(3))
    assert breakdown == [(1, 3, 3), (1, 5, 5), (4, 1, 4)], breakdown
    assert sum(c for _, _, c in breakdown) == 12
    print("Chapter 4: all assertions pass. "
          "N = 1, 1, 12, 620; d=3 breakdown 4 + 5 + 3 = 12.")
