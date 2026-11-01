"""Day 0: the tropical line.

Tropical arithmetic: replace (+) with min, (x) with (+).
A tropical polynomial f = min over terms of (c + i*x + j*y).
Its tropical curve V(f) = points where the minimum is attained >= twice.

Here: f = min(x, y, 0) -- the tropical line. Three rays from the origin,
dual to the three edges of its Newton polygon (convex hull of exponents).
Includes a numerical check of the balancing condition: primitive integer
direction vectors of the rays sum to zero.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Terms: (constant c, exponent (i, j)) for value c + i*x + j*y
TERMS = [(0.0, (1, 0)),   # x
         (0.0, (0, 1)),   # y
         (0.0, (0, 0))]   # 0

def eval_terms(X, Y):
    return np.array([c + i * X + j * Y for c, (i, j) in TERMS])

def tropical_curve_mask(X, Y, tol=1e-9):
    vals = eval_terms(X, Y)
    m = vals.min(axis=0)
    ties = (np.abs(vals - m) < tol).sum(axis=0)
    return ties >= 2

# --- sample the curve on a grid ---
n = 1200
xs = np.linspace(-3, 3, n)
ys = np.linspace(-3, 3, n)
X, Y = np.meshgrid(xs, ys)
mask = tropical_curve_mask(X, Y, tol=4e-3)

# --- balancing check: the three rays' primitive directions sum to 0 ---
rays = [((0, 1), "north", "x=0,y>=0 (terms 0,x tie)"),
        ((1, 0), "east", "y=0,x>=0 (terms 0,y tie)"),
        ((-1, -1), "southwest", "x=y<=0 (terms x,y tie)")]
bal = (rays[0][0][0] + rays[1][0][0] + rays[2][0][0],
       rays[0][0][1] + rays[1][0][1] + rays[2][0][1])
print("ray directions:", [r[0] for r in rays])
print("balancing sum:", bal, "->", "BALANCED" if bal == (0, 0) else "NOT BALANCED")

# --- figure: Newton polygon (dual) | tropical line ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.2))

# Newton polygon: convex hull of exponents (0,0), (1,0), (0,1)
tri = [(0, 0), (1, 0), (0, 1)]
ax1.add_patch(Polygon(tri, closed=True, fill=False, edgecolor="black", lw=2))
# color-code edges <-> dual rays
edge_colors = {"AB": "tab:blue", "AC": "tab:red", "BC": "tab:green"}
ax1.plot([0, 1], [0, 0], color=edge_colors["AB"], lw=5, solid_capstyle="round")   # bottom
ax1.plot([0, 0], [0, 1], color=edge_colors["AC"], lw=5, solid_capstyle="round")   # left
ax1.plot([1, 0], [0, 1], color=edge_colors["BC"], lw=5, solid_capstyle="round")   # hypotenuse
for (px, py), lab in [((0, 0), "0"), ((1, 0), "x"), ((0, 1), "y")]:
    ax1.text(px, py, f"  {lab}", fontsize=13)
ax1.set_xlim(-0.6, 1.6); ax1.set_ylim(-0.6, 1.6); ax1.set_aspect("equal")
ax1.set_title("Newton polygon of min(x, y, 0)\n(exponents of the terms)", fontsize=12)
ax1.set_xlabel("i"); ax1.set_ylabel("j")

# Tropical line
ax2.contour(X, Y, mask.astype(int), levels=[0.5], colors="black", linewidths=0)
ax2.imshow(mask, extent=(-3, 3, -3, 3), origin="lower", cmap="Greys", vmin=0, vmax=2)
for (dx, dy), name, where, col in [
        ((0, 1), "north", "x=0, y\u22650", edge_colors["AB"]),
        ((1, 0), "east", "y=0, x\u22650", edge_colors["AC"]),
        ((-1, -1), "southwest", "x=y\u22640", edge_colors["BC"])]:
    ax2.annotate("", xy=(2.2 * dx, 2.2 * dy), xytext=(0.15 * dx, 0.15 * dy),
                 arrowprops=dict(arrowstyle="->", color=col, lw=2.5))
    mx, my = 2.45 * dx, 2.45 * dy
    ax2.text(mx, my, f"{name}\n({dx},{dy})", color=col, fontsize=10,
             ha="center", va="center")
ax2.plot(0, 0, "ko", ms=7)
ax2.set_xlim(-3, 3); ax2.set_ylim(-3, 3); ax2.set_aspect("equal")
ax2.set_title('Tropical line V(min(x, y, 0))\n"where the minimum is attained twice"',
              fontsize=12)
ax2.set_xlabel("x"); ax2.set_ylabel("y")

fig.suptitle("Day 0 — tropical geometry: min-plus arithmetic turns a line into 3 rays;\n"
             "each ray is dual (perpendicular) to an edge of the Newton polygon, "
             "and the directions balance: (0,1)+(1,0)+(-1,-1) = 0",
             fontsize=11, y=1.02)
fig.tight_layout()
out = "/home/hatch/workspace/goals/tropical-geometry-study/files/figures/day00-tropical-line.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)
