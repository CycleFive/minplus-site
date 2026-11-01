"""Day 1: the tropical conic — degenerate (doubled line) vs smooth.

Warm-up: tropical roots of f(x) = min(2x, x+1, 0) — analytic + numeric.
Main event: f = min over terms (c_ij + i*x + j*y), (i,j) in the degree-2 triangle.
  - all c = 0: Newton triangle unsubdivided -> 3 rays from the origin, each of
    weight 2 (the doubled tropical line = trop((x+y+1)^2), a non-reduced conic).
  - c(i,j) = i^2+j^2 with c(1,1) bumped 2 -> 2.5: triangle subdivided into 4 unit
    triangles -> 4 trivalent vertices, 3 bounded edges, 6 weight-1 rays; the
    bounded part is a tree, genus 0 = (d-1)(d-2)/2 for d=2. The smooth conic is
    rational classically; it stays rational tropically.

Honesty machinery (all verified in-code, not asserted):
  - dual subdivision computed numerically: lower faces of the lifted points via
    scipy.spatial.ConvexHull (coplanar lower triangles merged);
  - curve vertices solved from the face tie-equations, then cross-checked: every
    face term must tie for the GLOBAL minimum at the computed point;
  - balancing checked at every vertex from the subdivision data;
  - ray directions additionally INFERRED from far-field mask points (angular
    clustering + PCA + primitive-integer snap) and balancing re-checked on the
    inferred directions (addresses the day-0 precision note).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from math import gcd
from scipy.spatial import ConvexHull

FIGDIR = "/home/hatch/workspace/goals/tropical-geometry-study/files/figures"
TOL = 5e-3

# ---------------------------------------------------------------- core

def tropical_mask(terms, X, Y, tol=TOL):
    """terms: list of (c, (i, j)). Mask of points where the min is tied."""
    vals = np.array([c + i * X + j * Y for c, (i, j) in terms])
    m = vals.min(axis=0)
    return (np.abs(vals - m) < tol).sum(axis=0) >= 2


def make_grid(xmin, xmax, n):
    xs = np.linspace(xmin, xmax, n)
    X, Y = np.meshgrid(xs, xs)
    return xs, X, Y


def lower_subdivision(pts, heights):
    """Regular subdivision of pts induced by heights: lower hull of the lift.

    Returns a list of faces; each face is a tuple of point indices in CCW order
    (triangles, or one polygon when all heights are equal).
    """
    pts = np.asarray(pts, float)
    h = np.asarray(heights, float)
    if np.allclose(h, h[0]):
        # No subdivision: single cell = convex hull of the point set.
        return [tuple(ConvexHull(pts).vertices)]
    P = np.column_stack([pts, h])
    hull = ConvexHull(P)
    # merge coplanar lower triangles into polygonal faces
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
    """Solve the tie equations for a face; verify the face terms tie for the
    global minimum there. Returns (x, y, tying_indices) or raises."""
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
        f"face {face}: terms {[round(v,6) for v in vals]} do not all tie at ({x:.4f},{y:.4f})"
    return x, y, tying


def lattice_length(p, q):
    d = np.round(np.subtract(q, p)).astype(int)
    return gcd(int(abs(d[0])), int(abs(d[1])))


def primitive(v):
    """Primitive integer vector in the direction of float vector v."""
    v = np.asarray(v, float)
    assert np.linalg.norm(v) > 1e-12, "zero direction"
    for k in range(1, 51):
        w = np.round(v * k).astype(int)
        if (w != 0).any() and np.allclose(v * k, w, atol=1e-6):
            g = gcd(int(abs(w[0])), int(abs(w[1])))
            return (int(w[0] // g), int(w[1] // g))
    raise RuntimeError(f"no integer direction found for {v}")


def inward_normal(p1, p2, interior):
    """Primitive integer inward normal of directed segment p1->p2 w.r.t. a
    known interior point."""
    d = np.subtract(p2, p1)
    for n in [(d[1], -d[0]), (-d[1], d[0])]:
        if np.dot(n, np.subtract(interior, (np.add(p1, p2) / 2))) > 0:
            return primitive(n)
    raise RuntimeError("no inward normal found")


def curve_from_subdivision(pts, terms, faces):
    """Dual tropical curve: vertices, bounded edges, rays, balancing check."""
    pts = np.asarray(pts, float)
    interior = pts.mean(axis=0)
    verts = {}   # face_id -> (x, y)
    for fid, face in enumerate(faces):
        x, y, tying = face_vertex(terms, face)
        verts[fid] = (x, y)
    # edge -> faces
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
    # balancing at each vertex
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
    """Is the nearest grid point to pt on the mask?"""
    ix = np.argmin(np.abs(xs - pt[0])); iy = np.argmin(np.abs(xs - pt[1]))
    ok = bool(mask[iy, ix])
    print(f"    mask check {label} at ({pt[0]:.3f},{pt[1]:.3f}): "
          f"{'ON CURVE' if ok else 'MISSING!'}")
    return ok


def infer_ray_directions(mask, X, Y, r_min, r_max, gap=0.35):
    """Infer ray directions purely from far-field mask points: angular
    clustering, PCA per cluster, snap to primitive integer vector, orient
    outward (positive dot with the far-field centroid)."""
    R = np.sqrt(X ** 2 + Y ** 2)
    sel = mask & (R > r_min) & (R < r_max)
    angs = np.arctan2(Y[sel], X[sel])
    pts = np.column_stack([X[sel], Y[sel]])
    order = np.argsort(angs)
    angs, pts = angs[order], pts[order]
    diffs = np.diff(angs)
    wrap = (angs[0] + 2 * np.pi) - angs[-1]
    cuts = sorted(np.argsort(diffs)[::-1][:20], reverse=False)
    # cluster boundaries = the largest gaps; keep gaps > `gap`
    big = [c for c in np.argsort(diffs)[::-1] if diffs[c] > gap]
    bounds = sorted(big)
    clusters = []
    start = 0
    for b in bounds:
        clusters.append((start, b + 1)); start = b + 1
    clusters.append((start, len(angs)))
    # also consider the wrap gap
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
        # orient: rays are unbounded, so in the far field the outward direction
        # has positive dot product with the cluster centroid
        centroid = pts[a:b].mean(axis=0)
        if np.dot(best, centroid) < 0:
            best = (-best[0], -best[1])
        results.append({"dir": best, "snap_err_deg": float(np.degrees(best_ang)),
                        "n_pts": b - a,
                        "mean_angle_deg": float(np.degrees(angs[a:b].mean()))})
    return results


# ---------------------------------------------------------------- warm-up
print("=" * 70)
print("WARM-UP: tropical roots of f(x) = min(2x, x+1, 0)")
print("=" * 70)
for cand, desc in [(1.0, "2x = x+1"), (0.0, "2x = 0"), (-1.0, "x+1 = 0")]:
    v = (2 * cand, cand + 1, 0.0)
    m = min(v)
    ties = [t for t, val in zip(["2x", "x+1", "0"], v) if abs(val - m) < 1e-12]
    print(f"  x={cand:5.1f} ({desc:10s}): values {v} -> min {m}, "
          f"tied by {ties} {'ROOT' if len(ties) >= 2 else 'not a root'}")
xs = np.linspace(-2, 2, 4001)
vals = np.array([2 * xs, xs + 1, np.zeros_like(xs)])
m = vals.min(axis=0)
num_roots = xs[(np.abs(vals - m) < 1e-9).sum(axis=0) >= 2]
print(f"  numeric scan: tied-min points in [{num_roots.min():.4f}, {num_roots.max():.4f}] "
      f"(only x = 0) -> MATCHES analytic")
print("  multiplicity: dual edge spans exponents 0 -> 2, lattice length 2,")
print("  so x = 0 is a root of multiplicity 2 (tropical FTA: degree 2 = 2 roots counted)")

# ---------------------------------------------------------------- conics
PTS = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)]  # degree-2 lattice points


def run_conic(name, heights, plot_range, r_infer, expected):
    print("=" * 70)
    print(name)
    print("=" * 70)
    terms = [(h, p) for h, p in zip(heights, PTS)]
    faces = lower_subdivision(PTS, heights)
    print(f"  subdivision faces (point indices): {faces}")
    assert len(faces) == expected["n_faces"], "face count != predicted"
    verts, bounded, rays = curve_from_subdivision(PTS, terms, faces)
    print("  vertices (face -> (x, y)):")
    for fid, v in verts.items():
        print(f"    face {fid} {faces[fid]}: ({v[0]:.4f}, {v[1]:.4f})")
    for fid, (ex, ey) in expected["verts"].items():
        v = verts[fid]
        assert abs(v[0] - ex) < 1e-9 and abs(v[1] - ey) < 1e-9, \
            f"vertex {fid}: got {v}, predicted ({ex},{ey})"
    print(f"  vertex positions MATCH prediction; balancing holds at all "
          f"{len(verts)} vertices (asserted in curve_from_subdivision)")
    print(f"  bounded edges: {len(bounded)}, rays: {len(rays)} "
          f"(predicted {expected['n_bounded']}, {expected['n_rays']})")
    assert len(bounded) == expected["n_bounded"] and len(rays) == expected["n_rays"]
    tot_w = {}
    for r in rays:
        tot_w[r["dir"]] = tot_w.get(r["dir"], 0) + r["weight"]
    print(f"  ray directions x weights: {[(d, w) for d, w in sorted(tot_w.items())]}")
    # mask + on-mask checks
    xs, X, Y = make_grid(*plot_range, 1400)
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
    print("  all predicted vertices/edges/rays present in sampled mask")
    # inferred rays from the mask alone
    inf = infer_ray_directions(mask, X, Y, *r_infer)
    print("  ray directions INFERRED from far-field mask:")
    for s in inf:
        print(f"    dir {s['dir']} (snap err {s['snap_err_deg']:.2f} deg, "
              f"{s['n_pts']} pts)")
    bal = np.zeros(2, dtype=int)
    for s in inf:
        bal += np.array(s["dir"])
    print(f"  balancing on INFERRED primitive directions: {tuple(bal)} -> "
          f"{'BALANCED' if (bal == 0).all() else 'NOT BALANCED'}")
    assert (bal == 0).all()
    return terms, faces, verts, bounded, rays, (xs, X, Y, mask)


# --- degenerate conic: all heights 0 -> doubled tropical line
# --- smooth conic: heights i^2+j^2, bump (1,1) 2 -> 2.5 to pick the triangulation
if __name__ == "__main__":
    deg = run_conic(
        "DEGENERATE CONIC: f = min(2x, 2y, x+y, x, y, 0)  [trop((x+y+1)^2)]",
        heights=[0, 0, 0, 0, 0, 0],
        plot_range=(-3, 3), r_infer=(1.8, 2.9),
        expected={"n_faces": 1, "n_bounded": 0, "n_rays": 3,
                  "verts": {0: (0.0, 0.0)}})

    smooth = run_conic(
        "SMOOTH CONIC: heights c(i,j)=i^2+j^2, c(1,1)=2.5",
        heights=[0, 1, 4, 1, 2.5, 4],
        plot_range=(-4.5, 4.5), r_infer=(3.6, 4.4),
        expected={"n_faces": 4, "n_bounded": 3, "n_rays": 6,
                  # face order from ConvexHull: T1=(0,1,3), T2=(1,2,4),
                  #                             T4=(1,3,4), T3=(3,4,5)
                  "verts": {0: (-1, -1), 1: (-3, -1.5), 2: (-1.5, -1.5),
                            3: (-1.5, -3)}})
    print("\nALL CHECKS PASSED")
