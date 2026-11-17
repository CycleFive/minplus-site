"""Chapter 13 -- Trees Are Invisible: chip-firing on the full tropical cubic.

Assembled for the book (Phase 1 draft, 2026-11-17) from Day 9 of the study:
chip_firing.py (the parameterized Baker-Norine machinery, inlined below as
the ChipFiring class) + day09_full_cubic.py (the analysis). Self-contained:
run as a script, every claim is asserted in code.

Chapter 12 put Baker-Norine theory on the 6-cycle alone. The actual smooth
tropical cubic of Chapter 3 has 9 vertices / 9 bounded edges: that cycle
plus 3 pendant vertices grafted on via bridges. This program rebuilds the
cubic from its Newton subdivision (not hardcoded), then asserts:

 1. Jac(G) ~= Z/6Z still: Matrix-Tree = 6 spanning trees; SNF of the
    reduced Laplacian = diag(1 x7, 6). Trees add no torsion.
 2. Push-in lemma: every divisor in {-1,0,1}^9 (19,683 of them) is
    linearly equivalent to one supported on the 2-core (the cycle).
 3. Abel-Jacobi factors through the cycle: S_q(pendant) =
    S_q(attachment) for all three trees; AJ is NOT injective on G --
    Baker-Norine Thm 1.8 needs 2-edge-connected, and the bridges are
    exactly the failure. The cycle generator still has order 6, and
    [pendant] + [x] = [attachment] + [x] for all x.
 4. K = sum (deg(v)-2)(v) is principal (one firing per leaf kills it);
    Riemann-Roch r(D) - r(K-D) = deg(D) verified on 8,953 divisors.
 5. Metric coda (with Chapter 12's Sec. 12.8 correction applied):
    lattice length L = 6 (not the Euclidean 4+2*sqrt(2)); the six cycle
    vertices sit at integer lattice arclengths, and the tropical
    Abel-Jacobi map collapses each tree to its attachment point
    (data exported to ch13_metric.json for the figure).

Sources: M. Baker & S. Norine, arXiv:math/0608360 (stabilization =
Prop. 3.1); Mikhalkin-Zharkov, arXiv:math/0612267 (metric Jacobian).
"""

import itertools
import json
import os
import sys
import time
from math import gcd

import numpy as np
from scipy.spatial import ConvexHull
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form


# ======================================================================
# Baker-Norine machinery (inlined from chip_firing.py)
# ======================================================================

class ChipFiring:
    """BN chip-firing on a finite connected graph.

    n: number of vertices (labeled 0..n-1); edges: list of (i, j) pairs;
    q: base vertex for q-reduced divisors and the Abel-Jacobi map.
    """

    def __init__(self, n, edges, q):
        self.n = n
        self.q = q
        self.edges = [(int(i), int(j)) for i, j in edges]
        self.adj = {v: set() for v in range(n)}
        for i, j in self.edges:
            assert i != j, "no loops"
            self.adj[i].add(j)
            self.adj[j].add(i)
        self.lap = np.zeros((n, n), dtype=int)
        for i, j in self.edges:
            self.lap[i, i] += 1
            self.lap[j, j] += 1
            self.lap[i, j] -= 1
            self.lap[j, i] -= 1
        # BFS layers from q (BN Prop 3.1: every v != q has a neighbor
        # in a strictly earlier layer)
        dist = {q: 0}
        dq = [q]
        for u in dq:
            for w in self.adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    dq.append(w)
        assert len(dist) == n, "graph is not connected"
        self.dist = dist
        self.order = sorted(range(n), key=lambda v: dist[v])
        others = [v for v in range(n) if v != q]
        self.subsets = [A for r in range(1, n)
                        for A in itertools.combinations(others, r)]
        # outdeg(A, v) = #{edges from v to V\A}, precomputed
        self._od = {}
        for A in self.subsets:
            s = set(A)
            for v in A:
                self._od[(A, v)] = sum(1 for w in self.adj[v] if w not in s)
        self._ev = [np.eye(n, dtype=int)[:, v] for v in range(n)]
        self.genus = len(self.edges) - n + 1
        self._cache = {}

    # ---- q-reduced divisors (G-parking functions) ---------------------
    def is_q_reduced(self, D):
        """P1: no debt away from q; P2: Dhar burning condition."""
        if any(D[v] < 0 for v in range(self.n) if v != self.q):
            return False
        for A in self.subsets:
            if not any(D[v] < self._od[(A, v)] for v in A):
                return False
        return True

    def stabilize(self, D, subset_order="lex", _trace=False):
        """BN Prop 3.1. Phase 1: out of debt, farthest-first via a
        preceding neighbor. Phase 2: repeatedly fire violating sets A
        (all v in A with D(v) >= outdeg_A(v)). Postconditions asserted:
        q-reduced, degree preserved, D - Dq = LAP @ f exactly."""
        D = np.array(D, dtype=int)
        D0 = D.copy()
        n, q = self.n, self.q
        f = np.zeros(n, dtype=int)
        fired_sets = []
        for v in reversed(self.order):
            if v == q:
                continue
            while D[v] < 0:
                preds = [w for w in self.adj[v]
                         if self.dist[w] == self.dist[v] - 1]
                assert preds, f"no preceding neighbor for {v}"
                w = preds[0]
                nn = -int(D[v])
                D = D - nn * self.lap[:, w]
                f = f + nn * self._ev[w]
        assert all(D[v] >= 0 for v in range(n) if v != q), "phase 1 failed"
        subsets = self.subsets if subset_order == "lex" else self.subsets[::-1]
        for _ in range(1000000):
            fired = False
            for A in subsets:
                if all(D[v] >= self._od[(A, v)] for v in A):
                    if _trace:
                        fired_sets.append((set(A), D.copy()))
                    for v in A:
                        D = D - self.lap[:, v]
                        f = f + self._ev[v]
                    fired = True
                    break
            if not fired:
                break
        else:
            raise AssertionError("phase 2 did not terminate")
        assert D.sum() == D0.sum(), "degree changed during stabilization"
        assert self.is_q_reduced(D), f"stabilized divisor {D} not q-reduced"
        assert np.array_equal(D0 - D, self.lap @ f), "D - Dq not principal!"
        return (D, fired_sets) if _trace else D

    def reduced_class(self, D, subset_order="lex"):
        """q-reduced representative of [D]; asserts D - Dq is principal."""
        return self.stabilize(np.array(D, dtype=int), subset_order)

    def reduced_cached(self, D):
        key = tuple(int(x) for x in D)
        if key not in self._cache:
            self._cache[key] = self.reduced_class(np.array(D, dtype=int))
        return self._cache[key]

    def same_class(self, D1, D2):
        return np.array_equal(self.reduced_class(D1), self.reduced_class(D2))

    def is_effective_class(self, D):
        """|D| nonempty <=> q-reduced rep is effective."""
        return bool(np.all(self.reduced_class(D) >= 0))

    # ---- rank / Riemann-Roch ------------------------------------------
    def effective_divisors(self, s):
        for combo in itertools.combinations_with_replacement(range(self.n), s):
            D = np.zeros(self.n, dtype=int)
            for v in combo:
                D[v] += 1
            yield D

    def rank(self, D):
        """BN: r(D) >= s iff |D - E| != empty for every effective E of
        degree s."""
        D = np.array(D, dtype=int)
        d = int(D.sum())
        if d < 0:
            return -1
        for s in range(d + 1):
            for E in self.effective_divisors(s):
                if not np.all(self.reduced_cached(D - E) >= 0):
                    return s - 1
        return d

    # ---- Abel-Jacobi ----------------------------------------------------
    def aj(self, i):
        """S_q(v_i) = [(v_i) - (q)], q-reduced rep."""
        D = np.zeros(self.n, dtype=int)
        D[i] += 1
        D[self.q] -= 1
        return self.reduced_class(D)

    def add_classes(self, c1, c2):
        return self.reduced_class(np.array(c1) + np.array(c2))

    # ---- the group itself -------------------------------------------------
    def spanning_tree_count(self):
        """Kirchhoff: |Jac(G)| = det of any reduced Laplacian."""
        q = self.q
        keep = [v for v in range(self.n) if v != q]
        red = self.lap[np.ix_(keep, keep)]
        return round(abs(np.linalg.det(red)))

    def snf_diagonal(self):
        """SNF of the reduced Laplacian: Jac(G) ~= Z/d_1 x ... x Z/d_k."""
        q = self.q
        keep = [v for v in range(self.n) if v != q]
        red = self.lap[np.ix_(keep, keep)]
        snf = smith_normal_form(Matrix(red.tolist()))
        # invariant factors are defined up to sign; sympy may return -6
        return sorted(abs(int(snf[i, i])) for i in range(self.n - 1))


# ======================================================================
# The full cubic (from day09_full_cubic.py)
# ======================================================================

TOL = 1e-9

# ---------------------------------------------------------------- curve machinery
# (copied from day02_tropical_cubic.py, which is not import-safe)


def lower_subdivision(pts, heights):
    pts = np.asarray(pts, float)
    h = np.asarray(heights, float)
    if np.allclose(h, h[0]):
        return [tuple(ConvexHull(pts).vertices)]
    P = np.column_stack([pts, h])
    hull = ConvexHull(P)
    groups = {}
    for s, eq in zip(hull.simplices, hull.equations):
        if eq[2] < -1e-9:
            n = eq[:3] / np.linalg.norm(eq[:3])
            key = tuple(np.round(np.append(n, eq[3] / np.linalg.norm(eq[:3])), 6))
            groups.setdefault(key, []).append(s)
    faces = []
    for tris in groups.values():
        verts = sorted({v for t in tris for v in t})
        if len(verts) == 3:
            faces.append(tuple(verts))
        else:
            sub = ConvexHull(pts[verts])
            faces.append(tuple(verts[i] for i in sub.vertices))
    return faces


def face_vertex(terms, face, check_tol=1e-6):
    idx = list(face)
    base = idx[0]
    c0, (i0, j0) = terms[base]
    A, b, used = [], [], []
    for k in idx[1:]:
        ck, (ik, jk) = terms[k]
        A.append([i0 - ik, j0 - jk])
        b.append(ck - c0)
        used.append(k)
        if len(used) == 2 and abs(np.linalg.det(np.array(A))) > 1e-9:
            break
    A = np.array(A[:2])
    b = np.array(b[:2])
    assert abs(np.linalg.det(A)) > 1e-9, f"face {face}: no 2 independent eqs"
    x, y = np.linalg.solve(A, b)
    vals = np.array([c + i * x + j * y for c, (i, j) in terms])
    m = vals.min()
    tying = [k for k, v in enumerate(vals) if abs(v - m) < check_tol]
    assert all(k in tying for k in idx), \
        f"face {face}: terms do not all tie at ({x:.4f},{y:.4f})"
    return x, y, tying


def lattice_length(p, q):
    d = np.round(np.subtract(q, p)).astype(int)
    return gcd(int(abs(d[0])), int(abs(d[1])))


def primitive(v):
    v = np.asarray(v, float)
    assert np.linalg.norm(v) > 1e-12, "zero direction"
    for k in range(1, 51):
        w = np.round(v * k).astype(int)
        if (w != 0).any() and np.allclose(v * k, w, atol=1e-6):
            g = gcd(int(abs(w[0])), int(abs(w[1])))
            return (int(w[0] // g), int(w[1] // g))
    raise RuntimeError(f"no integer direction found for {v}")


def inward_normal(p1, p2, interior):
    d = np.subtract(p2, p1)
    for n in [(d[1], -d[0]), (-d[1], d[0])]:
        if np.dot(n, np.subtract(interior, (np.add(p1, p2) / 2))) > 0:
            return primitive(n)
    raise RuntimeError("no inward normal found")


def curve_from_subdivision(pts, terms, faces):
    pts = np.asarray(pts, float)
    interior = pts.mean(axis=0)
    verts = {}
    for fid, face in enumerate(faces):
        x, y, tying = face_vertex(terms, face)
        verts[fid] = (x, y)
    edge_faces = {}
    for fid, face in enumerate(faces):
        m = len(face)
        for k in range(m):
            e = frozenset((face[k], face[(k + 1) % m]))
            edge_faces.setdefault(e, []).append(fid)
    bounded, rays = [], []
    for e, flist in edge_faces.items():
        a, b = tuple(e)
        p1, p2 = pts[a], pts[b]
        w = lattice_length(p1, p2)
        if len(flist) == 2:
            f1, f2 = flist
            v1, v2 = verts[f1], verts[f2]
            d = primitive(np.subtract(v2, v1))
            assert abs(np.dot(d, np.subtract(p2, p1))) < 1e-9, "duality perp failed"
            bounded.append({"edge": e, "faces": (f1, f2), "v1": v1, "v2": v2,
                            "dir": d, "weight": w})
        else:
            f = flist[0]
            n = inward_normal(p1, p2, interior)
            assert abs(np.dot(n, np.subtract(p2, p1))) < 1e-9, "ray perp failed"
            rays.append({"edge": e, "face": f, "origin": verts[f],
                         "dir": n, "weight": w})
    for fid, (x, y) in verts.items():
        bal = np.zeros(2)
        for be in bounded:
            f1, f2 = be["faces"]
            if fid == f1:
                s = int(np.sign(np.dot(be["dir"], np.subtract(be["v2"], be["v1"]))))
                bal += be["weight"] * s * np.array(be["dir"])
            elif fid == f2:
                s = int(np.sign(np.dot(be["dir"], np.subtract(be["v1"], be["v2"]))))
                bal += be["weight"] * s * np.array(be["dir"])
        for r in rays:
            if r["face"] == fid:
                bal += r["weight"] * np.array(r["dir"])
        assert (bal == 0).all(), f"vertex {fid}: balancing failed: {bal}"
    return verts, bounded, rays


# ---------------------------------------------------------------- graph combinatorics

def two_core(n, edges):
    """Iterative leaf-stripping; returns the surviving (2-core) vertex set."""
    adj = {v: set() for v in range(n)}
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    alive = set(range(n))
    while True:
        leaves = [v for v in alive
                  if sum(1 for w in adj[v] if w in alive) <= 1]
        # a single surviving vertex with a self... no; stop when no leaves
        # (an isolated vertex would be a leaf too -- none here)
        if not leaves:
            break
        # never strip the last two vertices of an edge being peeled: plain
        # leaf removal is fine since the core here is a genuine cycle
        for v in leaves:
            alive.discard(v)
    return alive


def tree_roots(n, edges, core):
    """For each non-core vertex, the core vertex its tree attaches to.

    The non-core vertices induce a forest; each of its trees meets the
    core in exactly one vertex (two distinct core neighbors would close
    a cycle through non-core vertices, contradicting maximality of the
    2-core)."""
    adj = {v: set() for v in range(n)}
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    seen = set()
    root = {}
    for v in range(n):
        if v in core or v in seen:
            continue
        # BFS the non-core component of v
        comp = set()
        stack = [v]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            comp.add(u)
            stack.extend(w for w in adj[u] if w not in core and w not in seen)
        core_nbrs = {w for u in comp for w in adj[u] if w in core}
        assert len(core_nbrs) == 1, \
            f"non-core component {sorted(comp)} touches core {sorted(core_nbrs)}"
        r = next(iter(core_nbrs))
        for u in comp:
            root[u] = r
    assert set(root) == {v for v in range(n) if v not in core}
    return root


def push_to_core(cf, D, core):
    """Move every off-core chip along its tree to the core attachment:
    at a leaf p with neighbor u, D[u] += D[p]; D[p] = 0 is firing (resp.
    borrowing) p exactly D[p] times -- hence principal. Returns the
    pushed divisor and the peel order."""
    n = cf.n
    D = np.array(D, dtype=int)
    cur = {v: set(w for w in cf.adj[v]) for v in range(n)}
    peeled = []
    while True:
        # peel a single leaf at a time: its neighbor is always still
        # present (in a tree rooted at the core, a leaf's unique neighbor
        # is its parent toward the core, never already peeled)
        leaves = [v for v in cur if v not in core and len(cur[v]) == 1]
        if not leaves:
            break
        p = leaves[0]
        u = next(iter(cur[p]))
        assert u in cur, f"leaf {p}'s neighbor {u} already peeled"
        # fire/borrow p exactly D[p] times toward u: D[u] += D[p], D[p] = 0
        D[u] += D[p]
        D[p] = 0
        del cur[p]
        for w in cur.values():
            w.discard(p)
        peeled.append((p, u))
    assert all(D[v] == 0 for v in range(n) if v not in core), \
        "push-in did not clear the trees"
    assert all(D[v] == 0 or v in core for v in range(n))
    return D, peeled


# ---------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=" * 70)
    print("CHAPTER 13: Trees Are Invisible -- the full tropical cubic")
    print("=" * 70)

    # ---- rebuild the day-2 smooth tropical cubic ----
    PTS = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1),
           (0, 2), (1, 2), (0, 3)]
    H = [i * i + j * j + i * j for (i, j) in PTS]
    terms = [(h, p) for h, p in zip(H, PTS)]
    faces = lower_subdivision(PTS, H)
    assert len(faces) == 9 and all(len(f) == 3 for f in faces)
    verts, bounded, rays = curve_from_subdivision(PTS, terms, faces)
    assert len(bounded) == 9 and len(rays) == 9
    print(f"  cubic rebuilt: 9 vertices, 9 bounded edges, 9 rays "
          f"({time.time()-t0:.1f}s)")

    # ---- the finite graph G: vertices = faces 0..8, edges = bounded ----
    n = 9
    gedges = []
    seen_pairs = set()
    for be in bounded:
        f1, f2 = be["faces"]
        assert f1 != f2
        key = tuple(sorted((f1, f2)))
        assert key not in seen_pairs, "parallel edge in curve graph?!"
        seen_pairs.add(key)
        gedges.append((f1, f2))
    assert len(gedges) == 9
    g = len(gedges) - n + 1
    print(f"  G: V={n}, E={len(gedges)}, genus g={g}")
    assert g == 1

    core = two_core(n, gedges)
    print(f"  2-core (minimal skeleton): {sorted(core)}")
    assert len(core) == 6
    # the core is 2-regular: the 6-cycle
    adj = {v: set() for v in range(n)}
    for i, j in gedges:
        adj[i].add(j)
        adj[j].add(i)
    assert all(sum(1 for w in adj[v] if w in core) == 2 for v in core), \
        "core is not 2-regular"
    trees = [v for v in range(n) if v not in core]
    print(f"  trees: {len(trees)} off-core vertices {sorted(trees)}")
    root = tree_roots(n, gedges, core)
    for v in sorted(trees):
        print(f"    vertex {v}: tree attaches at core vertex {root[v]}")
    bridges = [(i, j) for i, j in gedges
               if (i in core) != (j in core)]
    print(f"  bridges: {bridges} ({len(bridges)} of them)")

    q = sorted(core)[0]
    cf = ChipFiring(n, gedges, q)
    print(f"  base point q = {q} (on the cycle)")

    # ---- 1. the group: Matrix-Tree + SNF ----
    ntrees = cf.spanning_tree_count()
    print(f"  1. Matrix-Tree: {ntrees} spanning trees")
    assert ntrees == 6, "trees changed the tree count?!"
    diag = cf.snf_diagonal()
    print(f"     SNF(reduced Laplacian) = diag{tuple(diag)}")
    assert diag == [1] * 7 + [6], f"unexpected SNF {diag}"
    print("     => Jac(G) ~= Z/6Z: the grafted trees add no torsion")

    # ---- 2. push-in lemma: every class has a core-supported rep ----
    t1 = time.time()
    count = 0
    for tup in itertools.product((-1, 0, 1), repeat=n):
        D = np.array(tup, dtype=int)
        pushed, _ = push_to_core(cf, D, core)
        assert cf.same_class(D, pushed), f"push-in changed class of {tup}"
        count += 1
    print(f"  2. push-in lemma: all {count} divisors in {{-1,0,1}}^9 "
          f"are linearly equivalent to a core-supported divisor "
          f"({time.time()-t1:.1f}s)")

    # ---- 3. Abel-Jacobi factors through the cycle ----
    aj_classes = {}
    for v in range(n):
        aj_classes[v] = tuple(int(x) for x in cf.aj(v))
    image = sorted(set(aj_classes.values()))
    print(f"  3. AJ image: {len(image)} distinct classes (6 expected)")
    assert len(image) == 6
    for v in sorted(trees):
        assert aj_classes[v] == aj_classes[root[v]], \
            f"S_q({v}) != S_q({root[v]})"
        print(f"     S_q(pendant {v}) = S_q(attachment {root[v]}): "
            f"the tree is invisible to Abel-Jacobi")
    # AJ is not injective on G -- exactly the tree collapse
    assert len(set(aj_classes.values())) < n
    print("     AJ is not injective on G (BN Thm 1.8 needs 2-edge-connected;"
          " the bridges are exactly the failure)")
    # group law consistency: pendant acts like its attachment
    gen = aj_classes[sorted(core)[1]]  # a cycle neighbor of q
    assert not np.array_equal(np.array(gen), np.zeros(n, dtype=int)), \
        "generator is 0?!"
    multiples = set()
    c = np.zeros(n, dtype=int)
    for _ in range(6):
        multiples.add(tuple(int(x) for x in c))
        c = cf.add_classes(c, gen)
    assert np.array_equal(c, np.zeros(n, dtype=int)), "6*gen != 0"
    assert len(multiples) == 6, "generator does not have order 6"
    print("     S_q(cycle neighbor of q) has order 6: still generates Jac(G)")
    for v in sorted(trees):
        for x in range(n):
            assert np.array_equal(
                cf.add_classes(aj_classes[v], aj_classes[x]),
                cf.add_classes(aj_classes[root[v]], aj_classes[x]))
    print("     group law: [pendant] + [x] = [attachment] + [x] for all x")

    # ---- 4. K is principal; Riemann-Roch on G ----
    K = np.array([len(adj[v]) - 2 for v in range(n)], dtype=int)
    print(f"  4. canonical divisor K = {list(K)} "
          f"(+1 at attachments, -1 at leaves)")
    assert int(K.sum()) == 0 == 2 * g - 2
    assert cf.same_class(K, np.zeros(n, dtype=int)), "K not principal?!"
    print("     K is principal: K = sum over leaves of (laplacian column),"
          " i.e. one firing per leaf")
    t2 = time.time()
    n_rr = 0
    for tup in itertools.product((-1, 0, 1), repeat=n):
        D = np.array(tup, dtype=int)
        d = int(D.sum())
        if d not in (-1, 0, 1):
            continue
        rD = cf.rank(D)
        rKD = cf.rank(K - D)
        assert rD - rKD == d - g + 1 == d, \
            f"RR fails for {tup}: r(D)={rD}, r(K-D)={rKD}"
        n_rr += 1
    print(f"     Riemann-Roch r(D)-r(K-D)=deg(D) verified on {n_rr} "
          f"divisors (entries in {{-1,0,1}}, deg in {{-1,0,1}}) "
          f"({time.time()-t2:.1f}s)")

    # ---- 5. metric coda: export data for the figure ----
    # walk the cycle in order, cumulative arclength
    start = sorted(core)[0]
    order = [start]
    prev, cur = None, start
    while True:
        nbrs = [w for w in adj[cur] if w in core and w != prev]
        if prev is None:
            nxt = nbrs[0]  # first step: either neighbor
        else:
            assert len(nbrs) == 1, f"cycle walk stuck at {cur}"
            nxt = nbrs[0]
        if nxt == start:
            break
        prev, cur = cur, nxt
        order.append(cur)
    assert len(order) == 6 and set(order) == core, \
        f"cycle walk gave {order}"
    pos = {v: verts[v] for v in range(n)}
    # Lattice length of a cycle edge: gcd(|dx|,|dy|) of its displacement
    # (Ch. 12 Sec. 12.8 correction: the tropical metric is lattice length,
    # not Euclidean length -- each diagonal edge (1,1) has lattice length 1).
    def lat_len(a, b):
        dx = int(round(pos[b][0]-pos[a][0])); dy = int(round(pos[b][1]-pos[a][1]))
        return gcd(abs(dx), abs(dy))
    cum = {start: 0}
    L_lat = 0
    L_euc = 0.0
    s_lat = 0
    for a, b in zip(order, order[1:] + [start]):
        L_lat += lat_len(a, b)
        L_euc += float(np.linalg.norm(np.subtract(pos[b], pos[a])))
        s_lat += lat_len(a, b)
        if b != start:
            cum[b] = s_lat
    assert L_lat == 6, f"lattice cycle length {L_lat} != 6"
    assert abs(L_euc - (4 + 2 * np.sqrt(2))) < 1e-6
    L = float(L_lat)
    print(f"  5. metric cycle length: lattice L = {L_lat} "
          f"(Euclidean {L_euc:.6f} = 4+2*sqrt(2)); vertices at integer "
          f"lattice arclengths {cum}")
    metric = {
        "cycle_order": order,
        "arclength_frac": {str(v): cum[v] / L for v in order},
        "L": L,  # lattice length (= 6); Euclidean total is 4+2*sqrt(2)
        "pendant_root": {str(v): root[v] for v in trees},
        "pos": {str(v): [float(pos[v][0]), float(pos[v][1])] for v in range(n)},
        "bridges": bridges,
        "aj_class_index": {str(v): image.index(aj_classes[v]) for v in range(n)},
        "q": q,
    }
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch13_metric.json"), "w") as fh:
        json.dump(metric, fh, indent=1)
    print("     wrote ch13_metric.json")

    print(f"ALL DAY-9 ASSERTIONS PASSED ({time.time()-t0:.1f}s total)")


if __name__ == "__main__":
    main()
