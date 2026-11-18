"""Chapter 13 Phase 2 verification (2026-11-18).

Independent assertions NOT covered by ch13_chapter_program_draft.py's main():
 - the signed K equality exactly as it must read in the draft (Duet point 1)
 - firing vs borrowing directions, computed from the transformations
 - the full Abel-Jacobi table, row by row, as printed in draft section 13.5
 - raw pendant divisor (p)-(q) is NOT q-reduced; its reduced rep is the table row
 - the section 13.4 worked push-in instance
 - Try It 1-4 answers
 - bounded-graph degrees vs degrees counting rays (Duet point 3)
 - brute-force spanning-tree count (Try It 2 / section 13.2)
Imports the chapter program as a module (import-safe: main() guarded).
"""
import itertools, sys
import numpy as np
sys.path.insert(0, ".")
import ch13_chapter_program_draft as ch13

# Rebuild graph exactly as main() does
PTS = [(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(0,2),(1,2),(0,3)]
H = [i*i+j*j+i*j for (i,j) in PTS]
terms = [(h,p) for h,p in zip(H,PTS)]
faces = ch13.lower_subdivision(PTS,H)
verts, bounded, rays = ch13.curve_from_subdivision(PTS,terms,faces)
gedges = []
seen=set()
for be in bounded:
    f1,f2 = be["faces"]; key=tuple(sorted((f1,f2))); assert key not in seen; seen.add(key); gedges.append((f1,f2))
n=9; q=3
cf = ch13.ChipFiring(n, gedges, q)
Q = cf.lap
adj = cf.adj

# --- Duet point 3: degrees ---
bdeg = {v: len(adj[v]) for v in range(n)}
assert [bdeg[v] for v in range(9)] == [1,1,1,2,3,3,2,3,2], [bdeg[v] for v in range(9)]
ray_count = {v:0 for v in range(n)}
for r in rays: ray_count[r["face"]] += 1
assert [ray_count[v] for v in range(n)]==[2,2,2,1,0,0,1,0,1], ray_count
assert all(bdeg[v]+ray_count[v]==3 for v in range(n)), "every vertex trivalent counting its ray"
print("degrees: bounded", [bdeg[v] for v in range(9)], "; rays per vertex", [ray_count[v] for v in range(9)], "; bounded+rays = 3 everywhere OK")

# --- Signed K (Duet point 1) ---
K = np.array([len(adj[v])-2 for v in range(n)], dtype=int)
assert list(K)==[-1,-1,-1,0,1,1,0,1,0]
leaf_cols_sum = Q[:,0]+Q[:,1]+Q[:,2]
assert np.array_equal(leaf_cols_sum, -K), leaf_cols_sum
assert np.array_equal(-leaf_cols_sum, K)
f = np.array([-1,-1,-1,0,0,0,0,0,0], dtype=int)
assert np.array_equal(Q @ f, K), Q@f
# firing transformation D -> D - Q e_p, applied once per leaf to 0:
D_fire = np.zeros(n,dtype=int)
for p in (0,1,2): D_fire = D_fire - Q[:,p]
assert np.array_equal(D_fire, K), D_fire
# borrowing transformation D -> D + Q e_p, applied once per leaf to K:
D_borr = K.copy()
for p in (0,1,2): D_borr = D_borr + Q[:,p]
assert np.array_equal(D_borr, np.zeros(n,dtype=int)), D_borr
# and borrowing applied to 0 gives -K (so draft's old wording was inverted):
D_b0 = np.zeros(n,dtype=int)
for p in (0,1,2): D_b0 = D_b0 + Q[:,p]
assert np.array_equal(D_b0, -K)
print("K signed equality OK: K = -(Q[:,0]+Q[:,1]+Q[:,2]) = Q@(-1,-1,-1,0,...); firing leaves 0->K; borrowing leaves K->0; borrowing leaves 0->-K")

# --- AJ table exactly as draft 13.5 ---
expected = {
 0:[0,0,0,-1,1,0,0,0,0], 1:[0,0,0,-1,0,1,0,0,0], 2:[0,0,0,-1,0,0,0,1,0],
 3:[0,0,0,0,0,0,0,0,0], 4:[0,0,0,-1,1,0,0,0,0], 5:[0,0,0,-1,0,1,0,0,0],
 6:[0,0,0,-1,0,0,1,0,0], 7:[0,0,0,-1,0,0,0,1,0], 8:[0,0,0,-1,0,0,0,0,1]}
for v in range(9):
    row = [int(x) for x in cf.aj(v)]
    assert row==expected[v], (v,row)
    assert cf.is_q_reduced(np.array(row)), v
print("AJ table: all 9 rows match draft and are q-reduced OK")
# raw pendant divisor not q-reduced
for p in (0,1,2):
    raw=np.zeros(n,dtype=int); raw[p]+=1; raw[q]-=1
    assert not cf.is_q_reduced(raw), p
    assert np.array_equal(cf.reduced_class(raw), np.array(expected[p]))
print("raw (p)-(q) for pendants is NOT q-reduced; reduced rep = table row OK")

# --- 13.4 worked instance ---
D=np.array([2,-1,-1,-1,1,0,0,0,0],dtype=int)
core=ch13.two_core(n,gedges)
pushed,_=ch13.push_to_core(cf,D,core)
assert list(pushed)==[0,0,0,-1,3,-1,0,-1,0], list(pushed)
assert cf.same_class(D,pushed)
print("13.4 worked instance pushes to [0,0,0,-1,3,-1,0,-1,0], same class OK")

# --- Try It 2: brute force spanning trees ---
bridges={(0,4),(1,5),(2,7)}
def is_tree(edge_subset):
    parent=list(range(n))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for a,b in edge_subset:
        ra,rb=find(a),find(b)
        if ra==rb: return False
        parent[ra]=rb
    return len({find(v) for v in range(n)})==1
trees=[es for es in itertools.combinations(gedges,8) if is_tree(es)]
assert len(trees)==6, len(trees)
assert all(bridges.issubset({tuple(sorted(e)) for e in es}) for es in trees)
print("TryIt2: brute-force spanning trees = 6, every one contains all 3 bridges, 5 cycle edges OK")

# --- Try It 1: positions ---
# position = index in Z/6 via generator S_q(4) (position 1)
gen = tuple(int(x) for x in cf.aj(4))
pos={}
c=np.zeros(n,dtype=int)
for k in range(6):
    pos[tuple(int(x) for x in c)]=k
    c=cf.add_classes(c, np.array(gen))
for v,exp_pos in [(3,0),(4,1),(8,2),(7,3),(6,4),(5,5),(0,1),(1,5),(2,3)]:
    assert pos[tuple(int(x) for x in cf.aj(v))]==exp_pos, v
s = cf.add_classes(cf.aj(0), cf.aj(6))
assert pos[tuple(int(x) for x in s)]==5
assert np.array_equal(s, cf.aj(5)) and np.array_equal(s, cf.aj(1))
print("TryIt1: S_q(0)+S_q(6) = position 5 = class of vertex 5 and pendant 1 OK")

# --- Try It 3 coordinates ---
assert bdeg[0]==1 and K[0]==-1 and bdeg[4]==3 and K[4]==1
assert (Q@f)[0]==-1 and (Q@f)[4]==1
# firing leaf 0 alone: 0 -> -Q e0 has -1 at 0, +1 at 4
assert (-Q[:,0])[0]==-1 and (-Q[:,0])[4]==1
print("TryIt3: K(0)=-1, K(4)=+1, (Qf)(0)=-1, (Qf)(4)=+1 OK")

# --- Try It 4 ---
assert len({tuple(int(x) for x in cf.aj(v)) for v in range(9)})==6  # not injective
assert cf.spanning_tree_count()==6
print("TryIt4: Jac order still 6 but AJ not injective (6 classes for 9 vertices) OK")
print("ALL PHASE 2 EXTRA ASSERTIONS PASSED")
