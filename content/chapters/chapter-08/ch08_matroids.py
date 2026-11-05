#!/usr/bin/env python3
"""
Chapter 8: A Matroid Walks Into a Fan -- the Bergman fan.

The full program behind Chapter 8 of "Learn You a Tropical Geometry
for Great Good!" Converted from Day 6 of the tropical geometry study
(2026-10-16). SELF-CONTAINED apart from numpy/scipy.

A matroid M of rank r on E = {0, ..., n-1} IS a tropical linear
space: its Bergman fan B(M) in R^n / R*1. Three descriptions, from
three sources, are implemented and proved equal BY COMPUTATION:

  1. Initial matroids (Ardila-Klivans; Fry thesis arXiv:1910.00627):
     B(M) = {w : M^w is loop-free}, M^w = matroid of w-maximal bases.
  2. Chains of flats (arXiv:1010.3967, max convention): a chain
     empty < F_1 < ... < F_k < E of flats gives the cone generated
     by the indicator vectors e_{F_i}; B(M) is their union mod R*1.
  3. Circuit tropical hyperplanes (arXiv:1901.06893): w in B(M) iff
     on every circuit C the minimum of w over C is attained at
     least twice.

Sign convention (the chapter's Honest Box): the literature splits
between min-weight/negative-ray and max-weight/positive-ray; they
are negatives of each other and each description pairs with exactly
one circuit condition. This program uses max-weight/positive-ray
throughout, and the pairing it asserts is circuit MIN-twice -- the
first write-up paired max with max-twice and the assertions killed
it; the code determines the pairing, not the other way round.

Chapter claims asserted below (the table of section 8.2):
  matroid          rank  flats  circuits  rays  dim  2-cones
  U(2,4)             2      6        4       4    1      0
  parallel pair      2      5        3       3    1      0
  U(3,4)             3     12        1      10    2     12
  U(3,5)             3     17        5      15    2     20
For every matroid: the fan is pure of dimension r-1; the three
descriptions agree on 300 random integer weights (150 for U(3,5)),
and loop-free == circuit-min-twice on 2000 more; the ray count
equals the number of proper nonempty flats; the rank-2 fans are
balanced at the origin (ray generators sum to a multiple of 1);
U(3,4) and U(3,5) are balanced at every ray. The parallel-pair fan
has the merged ray e_0+e_1 and only 3 rays: the fan remembers the
matroid. Try It answers asserted too: B(U(2,5)) has 5 rays (its 5
proper nonempty flats are the singletons; triple agreement checked
on 100 samples), and dim B(U(3,4)) = r-1 = 2.

Running this file writes ch08_summary.json (the data the chapter's
figures are drawn from) into the current directory.
"""

import itertools
import numpy as np
from scipy.optimize import linprog


# ---------------------------------------------------------------- matroid

class Matroid:
    def __init__(self, n, bases):
        self.n = n
        self.bases = [frozenset(b) for b in bases]
        self.r = len(self.bases[0])
        assert all(len(b) == self.r for b in self.bases), "bases must be equinumerous"
        assert self.bases, "need at least one basis"

    def rank(self, S):
        S = set(S)
        return max(len(B & S) for B in self.bases)

    def closure(self, S):
        S = set(S)
        rS = self.rank(S)
        return frozenset(e for e in range(self.n) if self.rank(S | {e}) == rS)

    def flats(self):
        F = set()
        for r in range(self.n + 1):
            for S in itertools.combinations(range(self.n), r):
                F.add(self.closure(S))
        return sorted(F, key=lambda f: (len(f), sorted(f)))

    def circuits(self):
        deps = [frozenset(S)
                for r in range(1, self.n + 1)
                for S in itertools.combinations(range(self.n), r)
                if self.rank(S) < len(S)]
        return [D for D in deps if not any(C < D for C in deps)]

    @staticmethod
    def uniform(n, r):
        return Matroid(n, list(itertools.combinations(range(n), r)))

    @staticmethod
    def parallel_pair():
        """Rank 2 on {0,1,2,3} with 0||1: every 2-set a basis except {0,1}."""
        bases = [b for b in itertools.combinations(range(4), 2) if set(b) != {0, 1}]
        return Matroid(4, bases)


# ------------------------------------------------------- Bergman fan (1)

def bergman_cones(M):
    """All chains ∅=F_0 ⊂ F_1 ⊂ ... ⊂ F_k ⊂ F_{k+1}=E of flats (proper
    inclusions); cone = positive span of {e_{F_1},..,e_{F_k}} mod R*1.
    Returns list of chains (tuples of flats); chain [] = lineality only."""
    E = frozenset(range(M.n))
    inner = [F for F in M.flats() if F and F != E]
    chains = []

    def extend(chain):
        chains.append(tuple(chain))
        top = chain[-1] if chain else frozenset()
        for F in inner:
            if top < F:
                extend(chain + [F])

    extend([])
    return chains


def cone_contains(w, chain, n):
    """w in cone(e_F : F in chain) + R*1 ?  (LP feasibility)"""
    w = np.asarray(w, dtype=float)
    gens = []
    for F in chain:
        g = np.zeros(n)
        g[list(F)] = 1.0
        gens.append(g)
    if not gens:
        return bool(np.allclose(w - w[0] * np.ones(n), 0, atol=1e-9))
    G = np.column_stack(gens + [np.ones(n)])
    res = linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=w,
                  bounds=[(0, None)] * len(gens) + [(None, None)],
                  method="highs")
    return bool(res.success)


def in_bergman_fan(w, M, chains):
    return any(cone_contains(w, ch, M.n) for ch in chains)


# ------------------------------------------- Bergman fan (2): M^w loop-free

def max_weight_bases(M, w):
    wts = [sum(w[i] for i in B) for B in M.bases]
    best = max(wts)
    return [B for B, t in zip(M.bases, wts) if abs(t - best) < 1e-9]


def initial_loopless(M, w):
    covered = set().union(*max_weight_bases(M, w))
    return covered == set(range(M.n))


# --------------------------------- Bergman fan (3): circuit min-twice

def circuit_min_twice(M, w):
    for C in M.circuits():
        vals = sorted(w[i] for i in C)
        if vals[1] - vals[0] > 1e-9:  # minimum attained only once
            return False
    return True


# ------------------------------------------------------------- balancing

def ray_generators(M, chains):
    """Primitive integer generators of the 1-dim cones (rays), mod lineality."""
    rays = []
    for ch in chains:
        if len(ch) == 1:
            (F,) = ch
            g = np.zeros(M.n, dtype=int)
            g[list(F)] = 1
            rays.append(g)
    # dedupe mod positive scale (they are 0/1 vectors, so exact dedupe)
    uniq = []
    for g in rays:
        if not any(np.array_equal(g, h) for h in uniq):
            uniq.append(g)
    return uniq


def check_balanced_origin(rays, n):
    """dim-1 fan: sum of primitive ray generators == 0 mod Z*1."""
    s = np.sum(rays, axis=0)
    assert np.all(s == s[0]), f"unbalanced at origin: sum={s}"
    return s


def check_balanced_2d(M, chains, rays):
    """Pure 2-dim fan: at each ray rho (gen u), sum of the other generator
    of each 2-cone above rho is 0 mod span(u, 1). Returns (ray, alpha, beta)."""
    E = frozenset(range(M.n))
    twocones = [ch for ch in chains if len(ch) == 2]
    out = []
    for u in rays:
        total = np.zeros(M.n)
        for ch in twocones:
            F1, F2 = ch
            g1 = np.zeros(M.n); g1[list(F1)] = 1
            if np.array_equal(g1, u):
                g2 = np.zeros(M.n); g2[list(F2)] = 1
                total = total + g2
            else:
                g2 = np.zeros(M.n); g2[list(F2)] = 1
                if np.array_equal(g2, u):
                    total = total + g1
        # solve total = alpha*u + beta*1
        A = np.column_stack([u.astype(float), np.ones(M.n)])
        coef, res, *_ = np.linalg.lstsq(A, total, rcond=None)
        assert np.allclose(A @ coef, total, atol=1e-9), \
            f"unbalanced at ray {u}: sum={total}"
        out.append((u, coef))
    return out


# ------------------------------------------------------------------ main

def study(name, M, n_grid=300, n_pair=2000, seed=0):
    rng = np.random.default_rng(seed)
    chains = bergman_cones(M)
    flats = M.flats()
    print(f"== {name}: n={M.n} r={M.r}, flats={len(flats)}, "
          f"circuits={len(M.circuits())}, chains={len(chains)}")
    maxlen = max(len(ch) for ch in chains)
    assert maxlen == M.r - 1, f"not pure of dim r-1: {maxlen} vs {M.r - 1}"
    print(f"   pure of dimension r-1 = {M.r - 1}: OK")

    # triple agreement on random integer weights
    Ws = rng.integers(-2, 3, size=(n_grid, M.n)).astype(float)
    for i, w in enumerate(Ws):
        a = in_bergman_fan(w, M, chains)
        b = initial_loopless(M, w)
        c = circuit_min_twice(M, w)
        assert a == b == c, f"descriptions disagree at {w}: fan={a} loopfree={b} circ={c}"
    print(f"   fan == M^w-loop-free == circuit-min-twice on {n_grid} samples: OK")

    # pairwise loopfree vs circuit on more samples (cheap)
    Ws2 = rng.integers(-3, 4, size=(n_pair, M.n)).astype(float)
    for w in Ws2:
        assert initial_loopless(M, w) == circuit_min_twice(M, w), \
            f"loopfree/circuit disagree at {w}"
    print(f"   loopfree == circuit-min-twice on {n_pair} more samples: OK")

    rays = ray_generators(M, chains)
    print(f"   rays ({len(rays)}): {[list(map(int, r)) for r in rays]}")
    if M.r == 2:
        s = check_balanced_origin(rays, M.n)
        print(f"   balanced at origin: sum = {list(map(int, s))} = {int(s[0])}*1: OK")
    elif M.r == 3:
        res = check_balanced_2d(M, chains, rays)
        for u, coef in res:
            print(f"   balanced at ray {list(map(int, u))}: "
                  f"sum = {coef[0]:+.0f}*u {coef[1]:+.0f}*1: OK")
    return M, chains, rays


def fan_stats(M, chains, rays):
    """Counts the chapter table is made of, read off the fan itself."""
    E = frozenset(range(M.n))
    proper = [F for F in M.flats() if F and F != E]
    twocones = [ch for ch in chains if len(ch) == 2]
    # rays are exactly the proper nonempty flats, one e_F each
    assert len(rays) == len(proper), (len(rays), len(proper))
    return {
        "n": M.n, "rank": M.r, "flats": len(M.flats()),
        "circuits": len(M.circuits()), "rays": len(rays),
        "dim": M.r - 1, "chains": len(chains),
        "twocones": len(twocones), "proper_flats": len(proper),
        "ray_vectors": [[int(x) for x in r] for r in rays],
    }


if __name__ == "__main__":
    import json

    results = {}
    results["U(2,4)"] = study("U(2,4)", Matroid.uniform(4, 2))
    results["parallel-pair"] = study("parallel-pair (0||1)",
                                     Matroid.parallel_pair())
    results["U(3,4)"] = study("U(3,4)", Matroid.uniform(4, 3))
    results["U(3,5)"] = study("U(3,5)", Matroid.uniform(5, 3),
                              n_grid=150, seed=1)

    EXPECTED = {
        # name: (flats, circuits, rays, dim, twocones)
        "U(2,4)": (6, 4, 4, 1, 0),
        "parallel-pair": (5, 3, 3, 1, 0),
        "U(3,4)": (12, 1, 10, 2, 12),
        "U(3,5)": (17, 5, 15, 2, 20),
    }
    summary = {}
    for name, (M, chains, rays) in results.items():
        st = fan_stats(M, chains, rays)
        exp = EXPECTED[name]
        got = (st["flats"], st["circuits"], st["rays"], st["dim"],
               st["twocones"])
        assert got == exp, f"{name}: table mismatch {got} != {exp}"
        summary[name] = st
        print(f"[ok] {name}: table row (flats, circuits, rays, dim, "
              f"2-cones) = {got}")

    # The fan remembers the matroid: the parallel pair's merged ray.
    pp_rays = summary["parallel-pair"]["ray_vectors"]
    assert [1, 1, 0, 0] in pp_rays and len(pp_rays) == 3
    assert [1, 0, 0, 0] not in pp_rays and [0, 1, 0, 0] not in pp_rays
    print("[ok] parallel pair: rays are e_0+e_1, e_2, e_3 -- the two "
          "parallel elements share one ray")

    # U(3,4): the 12 triangular 2-cones are the fine subdivision of
    # the 6 coarse cones cone(e_i, e_j): every coarse pair (i, j)
    # contains exactly 2 fine cones (via e_i+e_j from each side).
    M34, chains34, _ = results["U(3,4)"]
    fine = [ch for ch in chains34 if len(ch) == 2]
    assert len(fine) == 12
    coarse_pairs = set()
    for F1, F2 in fine:
        coarse_pairs.add(frozenset(set(F1) | set(F2)))
    assert len(coarse_pairs) == 6, coarse_pairs
    print("[ok] U(3,4): 12 fine 2-cones subdivide exactly 6 coarse "
          "cones cone(e_i, e_j)")

    # Try It 1, asserted: B(U(2,5)) has 5 rays.
    M25, chains25, rays25 = study("U(2,5) (Try It)", Matroid.uniform(5, 2),
                                  n_grid=100, n_pair=500, seed=2)
    st25 = fan_stats(M25, chains25, rays25)
    assert st25["flats"] == 7 and st25["rays"] == 5
    summary["U(2,5)"] = st25
    print("[ok] Try It 1: U(2,5) has 7 flats, 5 of them proper "
          "nonempty, so B(U(2,5)) has 5 rays")

    with open("ch08_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print("\nwrote ch08_summary.json")
    print("Chapter 8: ALL ASSERTIONS PASSED -- three descriptions of "
          "the Bergman fan agree on every sampled weight for U(2,4), "
          "the parallel pair, U(3,4), U(3,5) (and U(2,5)); every fan "
          "is pure of dimension r-1 and balanced; the table, the "
          "parallel-pair merger and the Try It answers are as printed "
          "in the chapter.")
