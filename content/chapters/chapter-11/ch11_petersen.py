"""
Day 11 -- G'''_{2,5} is the Petersen graph (Speyer-Sturmfels, Ex. 3.3).

Source: D. Speyer & B. Sturmfels, "The tropical Grassmannian",
arXiv:math/0304218v3, verified against the arXiv HTML this session:
- Sec. 3 prime reductions: G_{d,n} fan in R^{C(n,d)} (max cone dim
  (n-d)d+1); G' = mod R*1 (dim d(n-d)); G'' = mod image(phi)
  (dim (d-1)(n-d-1)); G''' = G'' cap unit sphere (max face dim
  nd-n-d^2).  For (2,5): 7 in R^10, 6, 2, and 1 -- so G''' is a GRAPH.
- Ex. 3.3: G'''_{2,5} is the Petersen graph, 10 vertices, 15 edges.
- Thm. 3.4: G'''_{2,n} = tree space T_n: 2^{n-1}-n-1 vertices,
  (2n-5)!! facets, homotopy = bouquet of (n-2)! spheres of dim n-4.
- Sec. 4: Vert(T_n) = splits {A,B} (|A|,|B| >= 2, A u B = [n]);
  compatibility (4): A<=A' or A<=B' or B<=A' or B<=B'; T_n is the flag
  complex on that graph; facets = trivalent trees, |facet| = n-3.
- Four-point condition (6) + cone presentation (8), with day 10's SIGN
  CORRECTION: generators are -E_{A,B} (paper prints +E; with
  "-w_ij = distance" and (6) as printed, +E is not in the fan).

Conventions as day 10 (1-indexed leaves, min convention):
w in R^10 indexed by pairs; phi(a)_ij = a_i+a_j;
E_{A,B} = sum_{i in A, j in B} e_ij;  -w = tree distance.

What is computed (all asserted in-code, exact Fractions where linear):
1. Splits for n=5: all 2|3, so a split is its 2-side: 10 splits
   = 2^{4}-5-1.  Compatibility (4) <=> the 2-sides are DISJOINT
   (checked pair by pair).  Compatibility graph = Kneser KG(5,2).
2. The graph is Petersen: 10 vertices, 15 edges, 3-regular, girth 5,
   diameter 2, chi = 3, independence number 4, connected, cycle rank
   E-V+1 = 6 = (5-2)! -- Thm 3.4's bouquet of 6 circles (S^1's).
   Nonplanar by the girth-5 Euler bound E <= (5/3)(V-2) = 13.3 < 15.
   S_5 induces 120 distinct automorphisms, transitive on vertices
   and on edges (the classical Aut = S_5; we verify S_5 <= Aut).
3. Facets (max cliques) = the 15 edges themselves: each is a pair of
   disjoint 2-splits = a trivalent tree on 5 leaves (2 internal
   edges); count = (2*5-5)!! = 15; each split lies in exactly 3 trees.
4. Fan level: -E_A in G_{2,5} (four-point on all 5 quadruples), +E_A
   not.  For each facet {A,B}: w = phi(a) - s E_A - t E_B (s,t > 0)
   satisfies four-point on all quadruples, and (a,s,t) is recovered
   uniquely by exact rational linear algebra (cone is simplicial:
   rank [phi | E_A | E_B] = 7 = (n-d)d+1, Cor. 3.1; quotient cone
   dim 2 in G'').  Membership: four-point on all 5 quadruples <=>
   w in the union of the 15 closed cones -- 3000 random integer pts.
5. Initial forms of the 5 Plucker quadrics (one per 4-subset):
   interior of a maximal cone -> all 5 binomial; a ray point -E_A ->
   exactly 3 binomial + 2 trinomial (the quadruples meeting A in one
   leaf), none monomial; off the fan -> some quadric monomial.
   (The 5 quadrics generate I_{2,5}; Thm 3.4's proof promotes the
   quadric check to the full initial ideal via J_sigma = in_w(I).)
6. BALANCING at every ray of G'', in closed form: for a split A with
   complement A^c (|A^c| = 3), the three splits B disjoint from A are
   the 2-subsets of A^c, and
        sum_{B disjoint A} E_B  =  E_A + 2*phi( (1/2) * 1_{A^c} )
   exactly (proved entrywise in the log; asserted here for all 10).
   Hence sum of the neighbouring cone generators is 0 mod
   span(E_A, image phi): balanced.  All rays have the same lattice
   index by S_5-symmetry, so E-generators vs primitive generators
   differ by one uniform factor and balancing is unaffected.
7b. Phase 2 additions (2026-11-12, referee anchors): Petersen =
   complement of L(K5) checked pair by pair; rank(image phi) = 5 with
   (1,...,1) = phi((1/2)*1) (7 = 5 + 2 lineality bookkeeping);
   hypohamiltonian asserted by brute force (no Hamiltonian cycle,
   every one-vertex deletion Hamiltonian); Try It 1 and 3 answers
   asserted explicitly.
7. Corroboration at n=4 and n=6 with the same machinery:
   n=4: 3 splits, no two compatible -> G''' = 3 points = wedge of
   (4-2)! = 2 copies of S^0.  n=6 (Ex. 4.1): 25 splits (15 of type
   2|4, 10 of type 3|3), 105 compatibility edges, 105 triangles,
   every facet a triangle (pure dim 2), every edge in exactly 3
   triangles; 105 = 7!! facets, 25 = 2^5-6-1 vertices.
"""

from fractions import Fraction as Q
from itertools import combinations, permutations
import random

LEAVES = (1, 2, 3, 4, 5)
PAIRS = list(combinations(LEAVES, 2))          # 10 pairs, lexicographic
PIDX = {p: k for k, p in enumerate(PAIRS)}


def pair(i, j):
    return (i, j) if i < j else (j, i)


def phi(a):
    return [a[i - 1] + a[j - 1] for (i, j) in PAIRS]


def E_split(A):
    """E_{A, A^c} for a 2-subset A of [5] (the split's small side)."""
    A = set(A)
    w = [Q(0)] * 10
    for i in A:
        for j in LEAVES:
            if j not in A:
                w[PIDX[pair(i, j)]] += 1
    return w


def add(u, v):
    return [x + y for x, y in zip(u, v)]


def scale(c, u):
    return [c * x for x in u]


# ------------------------------------------------------------ linear algebra

def rref_solve(M, b):
    """Solve M x = b over Q.  M: list of rows (lists of Q), b: list of Q.
    Returns (solution, rank) with a unique solution, or (None, rank) if
    inconsistent or underdetermined."""
    nrows, ncols = len(M), len(M[0])
    A = [row[:] + [bb] for row, bb in zip(M, b)]
    piv_cols = []
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, nrows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(nrows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
        if r == nrows:
            break
    for i in range(nrows):
        if all(A[i][c] == 0 for c in range(ncols)) and A[i][ncols] != 0:
            return None, len(piv_cols)
    if len(piv_cols) < ncols:
        return None, len(piv_cols)
    x = [Q(0)] * ncols
    for i, c in enumerate(piv_cols):
        x[c] = A[i][ncols]
    return x, len(piv_cols)


def rank(M):
    """Rank over Q of a matrix given as list of rows."""
    A = [row[:] for row in M]
    nrows, ncols = len(A), len(A[0])
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, nrows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(nrows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == nrows:
            break
    return r


def columns_matrix(cols):
    """Matrix (rows = pair coords) whose columns are the given vectors."""
    return [[col[k] for col in cols] for k in range(10)]


PHI_COLS = [[Q(1) if leaf in p else Q(0) for p in PAIRS] for leaf in LEAVES]
# PHI_COLS[i] = phi(e_i) as a vector: (phi(e_i))_jk = [i in {j,k}]


# ------------------------------------------------------------ splits / graph

SPLITS = [frozenset(c) for c in combinations(LEAVES, 2)]   # the 2-side
ES = {A: E_split(A) for A in SPLITS}


def splits_compat(A, B):
    """Paper's (4): with B' = complement etc. A,B are the 2-sides;
    the split is {A, A^c}.  (4) on representatives (A,A^c),(B,B^c):
    A<=B or A<=B^c or A^c<=B or A^c<=B^c (any choice of reps)."""
    Ac, Bc = frozenset(LEAVES) - A, frozenset(LEAVES) - B
    return (A <= B or A <= Bc or Ac <= B or Ac <= Bc)


def build_graph(n, split_list, compat):
    adj = {v: set() for v in split_list}
    edges = []
    for A, B in combinations(split_list, 2):
        if compat(A, B):
            adj[A].add(B)
            adj[B].add(A)
            edges.append((A, B))
    return adj, edges


ADJ, EDGES = build_graph(5, SPLITS, splits_compat)


# ------------------------------------------------------------ four-point

QUADS = list(combinations(LEAVES, 4))          # 5 quadruples


def quad_sums(w, quad):
    i, j, k, l = quad
    return (w[PIDX[pair(i, j)]] + w[PIDX[pair(k, l)]],
            w[PIDX[pair(i, k)]] + w[PIDX[pair(j, l)]],
            w[PIDX[pair(i, l)]] + w[PIDX[pair(j, k)]])


def four_point_quad(w, quad):
    s = quad_sums(w, quad)
    m = min(s)
    return sum(1 for x in s if x == m) >= 2


def in_fan(w):
    return all(four_point_quad(w, q) for q in QUADS)


def quadric_initial(w, quad):
    """Initial form of the Plucker quadric of `quad`: # terms attaining
    the minimum pair-sum (1 = monomial, 2 = binomial, 3 = trinomial)."""
    s = quad_sums(w, quad)
    m = min(s)
    return sum(1 for x in s if x == m)


# ------------------------------------------------------------ cone membership

def decompose(w, A, B):
    """Solve w = phi(a) - s E_A - t E_B.  Returns (a, s, t) or None."""
    cols = PHI_COLS + [scale(Q(-1), ES[A]), scale(Q(-1), ES[B])]
    M = columns_matrix(cols)
    sol, rk = rref_solve(M, w[:])
    if sol is None:
        return None
    return sol[:5], sol[5], sol[6]


def in_cone_union(w):
    """w in some closed maximal cone (subfaces included via s,t >= 0)."""
    for (A, B) in EDGES:
        dec = decompose(w, A, B)
        if dec is not None and dec[1] >= 0 and dec[2] >= 0:
            return True
    return False


# ------------------------------------------------------------ graph algorithms

def bfs_dist(adj, src):
    dist = {src: 0}
    frontier = [src]
    while frontier:
        nxt = []
        for v in frontier:
            for u in adj[v]:
                if u not in dist:
                    dist[u] = dist[v] + 1
                    nxt.append(u)
        frontier = nxt
    return dist


def girth(adj):
    """Shortest cycle length via BFS from every vertex (small graph)."""
    best = None
    for s in adj:
        dist = {s: 0}
        parent = {s: None}
        frontier = [s]
        while frontier:
            nxt = []
            for v in frontier:
                for u in adj[v]:
                    if u not in dist:
                        dist[u] = dist[v] + 1
                        parent[u] = v
                        nxt.append(u)
                    elif parent[v] != u and parent[u] != v:
                        cyc = dist[v] + dist[u] + 1
                        if best is None or cyc < best:
                            best = cyc
            frontier = nxt
    return best


def chromatic(adj):
    """Exact chromatic number + a coloring, by backtracking (10 verts)."""
    verts = list(adj)
    for k in (1, 2, 3, 4):
        color = {}

        def try_color(idx):
            if idx == len(verts):
                return True
            v = verts[idx]
            used = {color[u] for u in adj[v] if u in color}
            for c in range(k):
                if c not in used:
                    color[v] = c
                    if try_color(idx + 1):
                        return True
                    del color[v]
            return False

        if try_color(0):
            return k, dict(color)
    raise AssertionError("no 4-coloring found?!")


def has_hamiltonian_cycle(adj, verts):
    """Brute-force DFS for a Hamiltonian cycle on `verts` in `adj`."""
    verts = list(verts)
    if len(verts) < 3:
        return False
    start = verts[0]
    path = [start]
    seen = {start}

    def dfs():
        if len(path) == len(verts):
            return start in adj[path[-1]]
        for u in adj[path[-1]]:
            if u in seen or u not in verts:
                continue
            seen.add(u)
            path.append(u)
            if dfs():
                return True
            path.pop()
            seen.remove(u)
        return False

    return dfs()


def independence_number(adj):
    verts = list(adj)
    best = 0
    for r in range(len(verts) + 1):
        for S in combinations(verts, r):
            if all(B not in adj[A] for A, B in combinations(S, 2)):
                best = max(best, r)
    return best


def max_cliques(adj):
    """All maximal cliques (Bron-Kerbosch, tiny graphs)."""
    out = []

    def bk(R, P, X):
        if not P and not X:
            out.append(frozenset(R))
            return
        for v in list(P):
            bk(R | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}

    bk(set(), set(adj), set())
    return out


# ------------------------------------------------------------ general T_n

def all_splits(n):
    """Splits {A, A^c} of [n], |A|,|A^c| >= 2, canonical frozenset pair."""
    leaves = tuple(range(1, n + 1))
    seen = set()
    out = []
    for r in range(2, n - 1):
        for A in combinations(leaves, r):
            A = frozenset(A)
            Ac = frozenset(leaves) - A
            key = frozenset((A, Ac))
            if key not in seen:
                seen.add(key)
                out.append(key)
    return out


def splits_compat_general(S, T):
    A, Ac = tuple(S)
    B, Bc = tuple(T)
    return A <= B or A <= Bc or Ac <= B or Ac <= Bc


def double_fact_odd(m):
    """(2k-5)!! style: product of odd numbers down from m."""
    out = 1
    while m > 1:
        out *= m
        m -= 2
    return out


# ------------------------------------------------------------ main

def main():
    rng = random.Random(20261021)

    # --- 1. splits, compatibility = disjointness, Kneser graph
    assert len(SPLITS) == 10 == 2 ** 4 - 5 - 1
    for A, B in combinations(SPLITS, 2):
        assert splits_compat(A, B) == (len(A & B) == 0), (A, B)
    assert len(EDGES) == 15
    assert all(len(ADJ[A]) == 3 for A in SPLITS)
    # Petersen as complement of L(K5) (Phase 2, referee anchor 1):
    # L(K5) has the same vertices (2-subsets), adjacent iff they SHARE
    # an element; complement adjacency is exactly disjointness.
    for A, B in combinations(SPLITS, 2):
        lk5_adjacent = len(A & B) == 1
        assert (B in ADJ[A]) == (not lk5_adjacent), (A, B)
    # Try It 1, explicitly: neighbours of 12 are exactly 34, 35, 45.
    assert ADJ[frozenset((1, 2))] == {frozenset((3, 4)), frozenset((3, 5)),
                                      frozenset((4, 5))}
    print("1. 10 splits (2|3); compatibility (4) <=> 2-sides disjoint; "
          "graph = Kneser KG(5,2) = complement of L(K5), 15 edges, 3-regular")

    # --- 2. Petersen invariants
    assert len(bfs_dist(ADJ, SPLITS[0])) == 10                      # connected
    diam = max(max(bfs_dist(ADJ, v).values()) for v in SPLITS)
    assert diam == 2
    g = girth(ADJ)
    assert g == 5
    chi, coloring = chromatic(ADJ)
    assert chi == 3 and all(coloring[A] != coloring[B] for A, B in EDGES)
    assert independence_number(ADJ) == 4
    cyc_rank = len(EDGES) - len(SPLITS) + 1
    assert cyc_rank == 6 == 3 * 2 * 1                                # (n-2)!
    # nonplanarity: planar + girth >= 5 forces E <= (5/3)(V-2) = 13.33
    assert len(EDGES) > Q(5, 3) * (len(SPLITS) - 2)
    # S_5 action: 120 distinct automorphisms, transitive on V and on E
    autos = set()
    for perm in permutations(LEAVES):
        img = {A: frozenset(perm[x - 1] for x in A) for A in SPLITS}
        assert all(img[B] in ADJ[img[A]] for A in SPLITS for B in ADJ[A])
        autos.add(tuple(img[A] for A in SPLITS))
    assert len(autos) == 120
    assert all(any(img[A] == B for img in
                   [{A0: frozenset(p[x - 1] for x in A0) for A0 in SPLITS}
                    for p in permutations(LEAVES)]) for A in SPLITS for B in SPLITS)
    edge_images = set()
    for perm in permutations(LEAVES):
        img = {A: frozenset(perm[x - 1] for x in A) for A in SPLITS}
        for (A, B) in EDGES:
            edge_images.add(frozenset((img[A], img[B])))
    assert edge_images == {frozenset(e) for e in EDGES}
    # Hypohamiltonian (Phase 2, referee closing warning): the intro's
    # word is now backed -- no Hamiltonian cycle on all 10 vertices,
    # but deleting ANY one vertex leaves a Hamiltonian graph.
    assert not has_hamiltonian_cycle(ADJ, set(SPLITS))
    for v in SPLITS:
        assert has_hamiltonian_cycle(ADJ, set(SPLITS) - {v}), v
    print("2b. hypohamiltonian: no Hamiltonian cycle; every "
          "one-vertex deletion is Hamiltonian (brute-force DFS, all 10)")
    print("2. Petersen: connected, diameter 2, girth 5, chi=3, alpha=4, "
          "cycle rank 6=(5-2)! (bouquet of 6 circles, Thm 3.4), "
          "nonplanar (Euler girth bound), S_5 -> 120 automorphisms, "
          "vertex- and edge-transitive")

    # --- 3. facets = edges = trivalent trees
    cliques = max_cliques(ADJ)
    assert {frozenset(e) for e in EDGES} == set(cliques)
    assert all(len(c) == 2 == 5 - 3 for c in cliques)                # |facet| = n-3
    assert len(cliques) == 15 == double_fact_odd(2 * 5 - 5)
    for A in SPLITS:
        assert sum(1 for c in cliques if A in c) == 3
    print("3. facets of T_5 = the 15 edges, each a trivalent tree "
          "(2 internal edges); (2n-5)!! = 15; each split in exactly 3 trees")

    # --- 4. fan level: signs, cones, membership equivalence
    # Lineality bookkeeping (Phase 2, referee anchor 2): image(phi)
    # is 5-dimensional (phi injective for n >= 3), so a maximal cone
    # splits as 7 = 5 (lineality) + 2 (G'' dimension); and
    # (1,...,1) = phi((1/2)*1) already lies in image(phi), so the
    # G' step (mod R*1) is redundant inside the G'' step.
    assert rank(columns_matrix(PHI_COLS)) == 5
    assert phi([Q(1, 2)] * 5) == [Q(1)] * 10
    for A in SPLITS:
        assert in_fan(scale(Q(-1), ES[A]))
        assert not in_fan(ES[A])
    for (A, B) in EDGES:
        cols = PHI_COLS + [ES[A], ES[B]]
        assert rank(columns_matrix(cols)) == 7 == (5 - 2) * 2 + 1    # Cor 3.1
    for trial in range(600):
        (A, B) = EDGES[rng.randrange(len(EDGES))]
        a = [Q(rng.randint(-3, 3)) for _ in range(5)]
        s, t = Q(rng.randint(1, 4)), Q(rng.randint(1, 4))
        w = add(phi(a), scale(-s, ES[A]))
        w = add(w, scale(-t, ES[B]))
        assert in_fan(w)
        dec = decompose(w, A, B)
        assert dec is not None and dec[0] == a and dec[1] == s and dec[2] == t
    n_test = 0
    for trial in range(3000):
        w = [Q(rng.randint(-3, 3)) for _ in range(10)]
        assert in_fan(w) == in_cone_union(w), w
        n_test += 1
    print(f"4. fan: -E in, +E out (sign correction persists); each facet "
          f"cone simplicial of dim 7 in R^10 (2 in G''); decomposition "
          f"recovered exactly (600x); four-point on 5 quadruples <=> "
          f"union of 15 cones ({n_test} random pts)")

    # --- 5. initial forms of the 5 Plucker quadrics
    (A0, B0) = EDGES[0]
    w_int = add(scale(Q(-1), ES[A0]), scale(Q(-1), ES[B0]))
    assert all(quadric_initial(w_int, q) == 2 for q in QUADS)
    for A in SPLITS:
        w_ray = scale(Q(-1), ES[A])
        kinds = sorted(quadric_initial(w_ray, q) for q in QUADS)
        assert kinds == [2, 2, 2, 3, 3], (A, kinds)
        tri = [q for q in QUADS if quadric_initial(w_ray, q) == 3]
        assert all(len(set(q) & A) == 1 for q in tri)
    # Try It 3, explicitly: at A = {1,2} the trinomial quadruples are
    # exactly 1345 and 2345 (Phase 2 assertion of the draft answer).
    w12 = scale(Q(-1), ES[frozenset((1, 2))])
    tri12 = {q for q in QUADS if quadric_initial(w12, q) == 3}
    assert tri12 == {(1, 3, 4, 5), (2, 3, 4, 5)}
    n_off = 0
    for trial in range(2000):
        w = [Q(rng.randint(-3, 3)) for _ in range(10)]
        has_monomial = any(quadric_initial(w, q) == 1 for q in QUADS)
        assert has_monomial == (not in_fan(w))
        n_off += 1
    print(f"5. in_w of the 5 Plucker quadrics: binomial x5 in cone "
          f"interiors; ray point: 3 binomial + 2 trinomial (quadruples "
          f"meeting the split in one leaf); monomial somewhere <=> off "
          f"fan ({n_off} random pts)")

    # --- 6. balancing, closed form
    for A in SPLITS:
        nbrs = sorted(ADJ[A])
        assert len(nbrs) == 3
        lhs = [Q(0)] * 10
        for B in nbrs:
            lhs = add(lhs, ES[B])
        a_vec = [Q(1, 2) if leaf not in A else Q(0) for leaf in LEAVES]
        rhs = add(ES[A], scale(Q(2), phi(a_vec)))
        assert lhs == rhs, A
        # hence sum of -E_B over neighbours = -E_A mod image(phi):
        # balanced at the ray rho_A in G''
        diff = add(scale(Q(-1), lhs), ES[A])          # = -2 phi(a_vec)
        sol, rk = rref_solve(columns_matrix(PHI_COLS), diff)
        assert sol is not None and phi(sol) == diff
    print("6. balancing at all 10 rays, closed form: "
          "sum_{B disjoint A} E_B = E_A + 2*phi((1/2)*1_{A^c}); "
          "neighbour generators sum to 0 mod span(E_A, image phi)")

    # --- 7. n=4 and n=6 corroboration (same machinery)
    for n, (nv, ne, nfac, tri_per_edge) in ((4, (3, 0, 3, None)),
                                            (6, (25, 105, 105, 3))):
        sp = all_splits(n)
        assert len(sp) == nv == 2 ** (n - 1) - n - 1
        adjn, edn = build_graph(n, sp, splits_compat_general)
        cliq = max_cliques(adjn)
        assert len(edn) == ne
        assert len(cliq) == nfac == double_fact_odd(2 * n - 5)
        assert all(len(c) == n - 3 for c in cliq)                    # pure
        if tri_per_edge is not None:
            for e in edn:
                assert sum(1 for c in cliq if set(e) <= c) == tri_per_edge
        print(f"   T_{n}: {nv} vertices, {ne} edges, {nfac} facets "
              f"(all of size {n - 3})"
              + (f", each edge in {tri_per_edge} facets (Ex. 4.1)" if tri_per_edge else
                 " -- 3 isolated points = G'''_{2,4}, wedge of 2 copies of S^0"))
    types6 = {}
    for S in all_splits(6):
        A, _ = tuple(S)
        types6[min(len(A), 6 - len(A))] = types6.get(min(len(A), 6 - len(A)), 0) + 1
    assert types6 == {2: 15, 3: 10}
    print("   T_6 split types: 15 of type 2|4 + 10 of type 3|3 (Ex. 4.1)")

    print("\nAll day-11 assertions passed.")


if __name__ == "__main__":
    main()
