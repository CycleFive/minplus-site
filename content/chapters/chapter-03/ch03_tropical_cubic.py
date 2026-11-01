"""Day 2: the tropical cubic — a tropical elliptic curve, genus 1.

f = min over (c_ij + i*x + j*y), (i,j) in the degree-3 lattice triangle,
with heights c(i,j) = i^2 + j^2 + i*j (a regular unimodular triangulation:
9 unit triangles; the star of the interior point (1,1) is the full 6-face star).

Predictions (from combinatorics, before running):
  - 9 faces -> 9 vertices; 18 subdivision edges, 9 boundary (lattice length 1)
    -> 9 weight-1 rays: 3 x (0,1), 3 x (1,0), 3 x (-1,-1) (degree 3 preserved);
  - 9 bounded edges; genus g = E - V + 1 = 9 - 9 + 1 = 1:
    exactly one cycle, the 6-cycle dual to the star of (1,1)
    (the other 3 bounded edges are bridges = trees grafted on the cycle);
  - hand-computed vertex checks: face {(0,0),(1,0),(0,1)} -> (-1,-1);
    face {(1,0),(2,0),(1,1)} -> (-3,-2); face {(1,0),(0,1),(1,1)} -> (-2,-2).
  - balancing holds at all 9 vertices (asserted in curve_from_subdivision).

Machinery (tropical_mask, lower_subdivision, face_vertex, curve_from_subdivision,
check_on_mask, infer_ray_directions) is copied unchanged from day01_tropical_conic.py
(not imported: that module runs its warm-up at import time).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from math import gcd
from scipy.spatial import ConvexHull
from collections import defaultdict

FIGDIR = "/home/hatch/workspace/goals/tropical-geometry-study/files/figures"
TOL = 5e-3

# ---------------------------------------------------------------- core (from day01)

def tropical_mask(terms, X, Y, tol=TOL):
    vals = np.array([c + i * X + j * Y for c, (i, j) in terms])
    m = vals.min(axis=0)
    return (np.abs(vals - m) < tol).sum(axis=0) >= 2


def make_grid(xmin, xmax, n):
    xs = np.linspace(xmin, xmax, n)
    X, Y = np.meshgrid(xs, xs)
    return xs, X, Y


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
    A = np.array(A[:2]); b = np.array(b[:2])
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


def check_on_mask(mask, xs, pt, label):
    ix = np.argmin(np.abs(xs - pt[0])); iy = np.argmin(np.abs(xs - pt[1]))
    ok = bool(mask[iy, ix])
    print(f"    mask check {label} at ({pt[0]:.3f},{pt[1]:.3f}): "
          f"{'ON CURVE' if ok else 'MISSING!'}")
    return ok


def infer_ray_directions(mask, X, Y, r_min, r_max, gap=0.35):
    R = np.sqrt(X ** 2 + Y ** 2)
    sel = mask & (R > r_min) & (R < r_max)
    angs = np.arctan2(Y[sel], X[sel])
    pts = np.column_stack([X[sel], Y[sel]])
    order = np.argsort(angs)
    angs, pts = angs[order], pts[order]
    diffs = np.diff(angs)
    big = [c for c in np.argsort(diffs)[::-1] if diffs[c] > gap]
    bounds = sorted(big)
    clusters = []
    start = 0
    for b in bounds:
        clusters.append((start, b + 1)); start = b + 1
    clusters.append((start, len(angs)))
    results = []
    for a, b in clusters:
        if b - a < 8:
            continue
        C = pts[a:b] - pts[a:b].mean(axis=0)
        _, _, Vt = np.linalg.svd(C, full_matrices=False)
        d = Vt[0]
        best, best_ang = None, 1e9
        for i in range(-3, 4):
            for j in range(-3, 4):
                if (i, j) == (0, 0) or gcd(abs(i), abs(j)) != 1:
                    continue
                v = np.array([i, j], float); v /= np.linalg.norm(v)
                ang = abs(np.arccos(np.clip(abs(np.dot(d, v)), -1, 1)))
                if ang < best_ang:
                    best, best_ang = (i, j), ang
        centroid = pts[a:b].mean(axis=0)
        if np.dot(best, centroid) < 0:
            best = (-best[0], -best[1])
        results.append({"dir": best, "snap_err_deg": float(np.degrees(best_ang)),
                        "n_pts": b - a,
                        "mean_angle_deg": float(np.degrees(angs[a:b].mean()))})
    return results


# ---------------------------------------------------------------- the cubic
PTS = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (0, 3)]
H = [i * i + j * j + i * j for (i, j) in PTS]
INTERIOR = PTS.index((1, 1))  # 5

print("=" * 70)
print("TROPICAL CUBIC: heights c(i,j) = i^2 + j^2 + i*j")
print("=" * 70)
terms = [(h, p) for h, p in zip(H, PTS)]
faces = lower_subdivision(PTS, H)
print(f"  faces: {len(faces)} (predicted 9)")
assert len(faces) == 9, "not a full triangulation"
P = np.array(PTS)
for f in faces:
    A = P[list(f)]
    assert len(f) == 3, f"non-triangular face {f}"
    area2 = abs(np.cross(np.append(A[1] - A[0], 0), np.append(A[2] - A[0], 0))[2])
    assert area2 == 1, f"face {f} not a unit triangle"
print("  all 9 faces are unit triangles -> smooth tropical cubic")

# hand-computed vertex checks (tie equations solved analytically)
face_id = {frozenset(f): fid for fid, f in enumerate(faces)}
hand = {frozenset((0, 1, 4)): (-1.0, -1.0),
        frozenset((1, 2, 5)): (-3.0, -2.0),
        frozenset((1, 4, 5)): (-2.0, -2.0)}
for fs, (ex, ey) in hand.items():
    fid = face_id[fs]
    x, y, _ = face_vertex(terms, faces[fid])
    assert abs(x - ex) < 1e-9 and abs(y - ey) < 1e-9, \
        f"face {tuple(fs)}: got ({x},{y}), hand-computed ({ex},{ey})"
    print(f"  hand check: face {sorted(fs)} -> ({x:.1f},{y:.1f}) MATCHES analytic")

verts, bounded, rays = curve_from_subdivision(PTS, terms, faces)
print(f"  vertices: {len(verts)}, bounded edges: {len(bounded)} (predicted 9), "
      f"rays: {len(rays)} (predicted 9)")
assert len(bounded) == 9 and len(rays) == 9
print("  balancing holds at all 9 vertices (asserted in curve_from_subdivision)")

ray_dirs = sorted(r["dir"] for r in rays)
from collections import Counter
print(f"  ray directions: {dict(Counter(ray_dirs))} (predicted 3x each of "
      f"(0,1),(1,0),(-1,-1), all weight 1)")
assert Counter(ray_dirs) == Counter({(0, 1): 3, (1, 0): 3, (-1, -1): 3})
assert all(r["weight"] == 1 for r in rays)

# ---- the cycle: leaf-strip the bounded graph; what remains is the cycle
adj = defaultdict(set)
for be in bounded:
    f1, f2 = be["faces"]
    adj[f1].add(f2); adj[f2].add(f1)
remaining = {v: set(n) for v, n in adj.items()}
while True:
    leaves = [v for v, n in remaining.items() if len(n) == 1]
    if not leaves:
        break
    for v in leaves:
        (u,) = remaining[v]
        remaining[u].discard(v)
        del remaining[v]
cycle_v = set(remaining)
print(f"  leaf-stripping leaves {len(cycle_v)} vertices: {sorted(cycle_v)}")
assert len(cycle_v) == 6, f"expected 6-cycle, got {len(cycle_v)}"
assert all(len(n) == 2 for n in remaining.values()), "remainder is not a 2-regular cycle"
cycle_e = [be for be in bounded
           if be["faces"][0] in cycle_v and be["faces"][1] in cycle_v]
assert len(cycle_e) == 6
# duality with the star of the interior lattice point
star = {fid for fid, f in enumerate(faces) if INTERIOR in f}
print(f"  star of (1,1): faces {sorted(star)} (predicted 6)")
assert star == cycle_v, "cycle vertices != faces of the star of (1,1)"
for be in cycle_e:
    assert INTERIOR in be["edge"], \
        f"cycle edge {be['edge']} not dual to a star edge of (1,1)"
print("  CYCLE = dual of the star of (1,1): 6 vertices, 6 edges, each dual edge "
      "incident to (1,1) -- MATCHES prediction")
# every non-cycle bounded edge is a bridge
nontree = [be for be in bounded if be not in cycle_e]
assert len(nontree) == 3
for be in nontree:
    f1, f2 = be["faces"]
    seen, stack = {f1}, [f1]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if (u == f1 and w == f2) or (u == f2 and w == f1):
                continue
            if w not in seen:
                seen.add(w); stack.append(w)
    assert f2 not in seen, f"edge {be['faces']} is not a bridge"
print(f"  the other {len(nontree)} bounded edges are all bridges "
      f"(trees grafted on the cycle)")
# genus + cycle length
g = len(bounded) - len(verts) + 1
cyc_len = sum(np.linalg.norm(np.subtract(be["v2"], be["v1"])) for be in cycle_e)
print(f"  genus g = E - V + 1 = {len(bounded)} - {len(verts)} + 1 = {g} "
      f"(predicted 1 = # interior lattice points)")
assert g == 1
print(f"  cycle length (Euclidean): {cyc_len:.6f}")
for be in cycle_e:
    assert be["weight"] == 1

print("\nvertex positions (face -> (x, y), * = on the cycle):")
for fid in sorted(verts):
    mark = "*" if fid in cycle_v else " "
    print(f"   {mark} face {fid} {faces[fid]}: ({verts[fid][0]:.4f}, {verts[fid][1]:.4f})")

# ---- on-mask checks
xs, X, Y = make_grid(-8, 4, 1500)
mask = tropical_mask(terms, X, Y)
ok = True
for fid, v in verts.items():
    ok &= check_on_mask(mask, xs, v, f"vertex {fid}")
for be in bounded:
    mid = ((be["v1"][0] + be["v2"][0]) / 2, (be["v1"][1] + be["v2"][1]) / 2)
    ok &= check_on_mask(mask, xs, mid, f"bounded edge {be['faces']}")
for r in rays:
    tip = (r["origin"][0] + 1.2 * r["dir"][0], r["origin"][1] + 1.2 * r["dir"][1])
    ok &= check_on_mask(mask, xs, tip, f"ray {r['dir']} w={r['weight']}")
assert ok, "some predicted curve points missing from mask!"
print("  all vertices/edges/rays present in the sampled mask")

# ---- ray inference from the mask alone (centered on the mask centroid,
# since the curve sits far from the grid origin: angular separation of
# parallel rays is measured about the curve's own center)
print("  ray directions INFERRED from far-field mask (mask-centroid-centered):")
cy, cx = Y[mask].mean(), X[mask].mean()
ang = np.arctan2(Y - cy, X - cx)
rr = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
inferred = None
for r_rng, gap in [((3.5, 5.9), 0.2), ((4.0, 5.9), 0.2), ((3.5, 5.9), 0.15)]:
    r0, r1 = r_rng
    sel = mask & (rr > r0) & (rr < r1)
    a = np.sort(ang[sel])
    pts = np.column_stack([X[sel], Y[sel]])
    order = np.argsort(ang[sel])
    a, pts = ang[sel][order], pts[order]
    diffs = np.diff(a)
    big = sorted(c for c in np.argsort(diffs)[::-1] if diffs[c] > gap)
    bounds = [(0, big[0] + 1)] + \
        [(big[k] + 1, big[k + 1] + 1) for k in range(len(big) - 1)] + \
        ([(big[-1] + 1, len(a))] if big else [(0, len(a))])
    inf = []
    for s0, s1 in bounds:
        if s1 - s0 < 8:
            continue
        C = pts[s0:s1] - pts[s0:s1].mean(axis=0)
        _, _, Vt = np.linalg.svd(C, full_matrices=False)
        d = Vt[0]
        best, best_ang = None, 1e9
        for i in range(-3, 4):
            for j in range(-3, 4):
                if (i, j) == (0, 0) or gcd(abs(i), abs(j)) != 1:
                    continue
                v = np.array([i, j], float); v /= np.linalg.norm(v)
                ang_err = abs(np.arccos(np.clip(abs(np.dot(d, v)), -1, 1)))
                if ang_err < best_ang:
                    best, best_ang = (i, j), ang_err
        centroid = pts[s0:s1].mean(axis=0)
        if np.dot(best, centroid - [cx, cy]) < 0:
            best = (-best[0], -best[1])
        inf.append({"dir": best, "snap_err_deg": float(np.degrees(best_ang)),
                    "n_pts": s1 - s0})
    print(f"    r in {r_rng} (about mask centroid), gap {gap}: {len(inf)} clusters -> "
          f"{sorted(s['dir'] for s in inf)}")
    if len(inf) == 9:
        inferred = inf
        break
assert inferred is not None, "could not separate all 9 rays by angular clustering"
bal = np.zeros(2, dtype=int)
for s in inferred:
    bal += np.array(s["dir"])
cnt = Counter(s["dir"] for s in inferred)
print(f"  inferred multiset: {dict(cnt)}, balancing sum {tuple(bal)}, "
      f"max snap err {max(s['snap_err_deg'] for s in inferred):.2f} deg")
assert cnt == Counter({(0, 1): 3, (1, 0): 3, (-1, -1): 3})
assert (bal == 0).all()
assert max(s["snap_err_deg"] for s in inferred) < 1.0
print("  all 9 rays recovered from the mask alone; inferred directions balance")

print("\nALL CHECKS PASSED")
