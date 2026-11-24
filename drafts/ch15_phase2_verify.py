#!/usr/bin/env python3
"""Ch.15 Phase 2 verification (2026-11-24).

Independent claim-by-claim checks for the Chapter 15 draft
("The Matroid With No Matrices", Vamos V8), beyond what the draft
program itself asserts:

  * Try It 1-4 answers, asserted individually;
  * the full Whitney/chi/chi-bar polynomials for V8 and U(4,8)
    (draft sec 15.1/15.3 state them; the draft program asserts
    only m and two mu values);
  * log-concavity with the correct products m_{k-1}*m_{k+1}
    (catches the draft's shorthand "21^2 = 441 >= 30", which is
    not the log-concavity inequality -- the true RHS is 7*30=210);
  * Duet fix 2 as a PROOF: per circuit-hyperplane, the Mobius
    accounting (4 absorbed triple flats mu=-1 each, replaced by
    one flat mu=-3) moves w_3 by exactly +1; five defects move
    w_3 from -56 to -51, hence m_3 from 35 to 30;
  * sec 15.6 counts: 77 proper NONEMPTY flats, quotient dim 7,
    fan dim r-1 = 3;
  * sec 15.5's 52,670 incomparable pairs are ORDERED pairs
    (26,335 unordered) -- Duet note 4;
  * exact deg(l^3) = 150,720 is checked against the draft run's
    summary JSON by the caller (the draft program asserts only
    deg > 0).
"""
import importlib.util
import itertools
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "day14", os.path.join(HERE, "day14_rank4_chow.py"))
day14 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(day14)
mobius_chibar = day14.mobius_chibar

E = list(range(8))
BLOCKS = [frozenset({0, 1}), frozenset({2, 3}),
          frozenset({4, 5}), frozenset({6, 7})]
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]  # all but (2,3)=C,D
H = [BLOCKS[i] | BLOCKS[j] for i, j in PAIRS]
Hset = set(H)


def rho(S):
    S = frozenset(S)
    return 3 if S in Hset else min(len(S), 4)


allsubs = [frozenset(c) for s in range(9)
           for c in itertools.combinations(E, s)]
rk_of = {S: rho(S) for S in allsubs}


def closure(S):
    S = set(S)
    r = rk_of[frozenset(S)]
    return frozenset(e for e in E
                     if e in S or rk_of[frozenset(S | {e})] == r)


flats = sorted({closure(S) for S in allsubs},
               key=lambda f: (rk_of[f], len(f), sorted(f)))
flatset = set(flats)

# ---- Try It 1: one rank by hand ----
A, B, C, D = BLOCKS
assert A | C in Hset and rho(A | C) == 3
assert C | D not in Hset and rho(C | D) == 4
print("TryIt1: rho(A+C)=3 (in H), rho(C+D)=4 (basis): OK")

# ---- Try It 2: Ingleton RHS terms ----
rhs_sets = [A | B, A | C, A | D, B | C, B | D]
assert rhs_sets == H and all(rho(s) == 3 for s in rhs_sets)
assert sum(rho(s) for s in rhs_sets) == 15
lhs = rho(A) + rho(B) + rho(A | B | C) + rho(A | B | D) + rho(C | D)
assert lhs == 16 and lhs > 15
print("TryIt2: RHS five terms are exactly H, each rank 3, sum 15; LHS 16: OK")

# ---- Whitney numbers, chi, chi-bar (V8 and uniform) ----
def whitney(flats_, rk_):
    m_, mu_ = mobius_chibar(E, flats_, rk_)
    r = max(rk_(f) for f in flats_)
    w = [0] * (r + 1)
    for f in flats_:
        w[rk_(f)] += int(mu_[f])
    return w, m_, mu_

w_v, m_v, mu_v = whitney(flats, lambda f: rk_of[f])
assert w_v == [1, -8, 28, -51, 30], w_v          # chi coeffs (signed)
assert m_v == [1, 7, 21, 30], m_v                # chi-bar unsigned coeffs
# chi(t) = t^4 - 8t^3 + 28t^2 - 51t + 30; chi-bar = t^3 - 7t^2 + 21t - 30
cum, s = [], 0
for x in w_v[:-1]:
    s += x
    cum.append(abs(s))
assert cum == m_v  # m_k = |w_0+...+w_k| (synthetic division by t-1)
print(f"V8 Whitney w={w_v} -> chi t^4-8t^3+28t^2-51t+30, "
      f"chi-bar coeffs {m_v}: OK")

uflats = [frozenset(c) for s in range(4)
          for c in itertools.combinations(E, s)] + [frozenset(E)]
urk = lambda f: min(len(f), 4)
w_u, m_u, mu_u = whitney(uflats, urk)
assert w_u == [1, -8, 28, -56, 35], w_u
assert m_u == [1, 7, 21, 35], m_u
print(f"U(4,8) Whitney w={w_u}, m={m_u}: OK")

# ---- Try It 3: the ledger ----
cost = [u - v for u, v in zip(m_u, m_v)]
assert cost == [0, 0, 0, 5], cost
# log-concavity, correct products (draft shorthand fixed in Phase 2)
assert 7 ** 2 >= 1 * 21 and 21 ** 2 >= 7 * 30
assert 21 ** 2 == 441 and 7 * 30 == 210
print("TryIt3: cost (0,0,0,5); log-concavity 49>=1*21, 441>=7*30=210: OK")

# ---- Duet fix 2: per-defect price is a proof, not an induction ----
# In U(4,8) every 3-subset is a flat with mu = -1.
# In V8, for each circuit-hyperplane Hf (rank r-1 = 3, |Hf| = r = 4):
#   * its r = 4 subsets of size r-1 = 3 are NOT flats (closure = Hf);
#   * the flats strictly below Hf are exactly: empty, 4 singletons,
#     6 pairs -- so mu(Hf) = -(1 - 4 + 6) = -3 = (-1)^(r-1) * (r-1);
#   * w_3 therefore moves by -3 - 4*(-1) = +1 per defect.
absorbed = []
for Hf in H:
    subs3 = [frozenset(c) for c in itertools.combinations(sorted(Hf), 3)]
    assert len(subs3) == 4
    for t in subs3:
        assert t in set(uflats) and int(mu_u[t]) == -1
        assert t not in flatset and closure(t) == Hf
    below = [f for f in flats if f < Hf]
    assert Counter((rk_of[f], len(f)) for f in below) == Counter(
        {(0, 0): 1, (1, 1): 4, (2, 2): 6}), below
    assert int(mu_v[Hf]) == -3 == (-1) ** 3 * (4 - 1)
    absorbed.extend(subs3)
assert len(absorbed) == 20 and len(set(absorbed)) == 20  # distinct
n_rk3_u = sum(1 for f in uflats if urk(f) == 3)
n_rk3_v = sum(1 for f in flats if rk_of[f] == 3)
assert n_rk3_u == 56 and n_rk3_v == 36 + 5 == 41 == 56 - 20 + 5
assert w_v[3] == 36 * (-1) + 5 * (-3) == -51
assert w_v[3] - w_u[3] == 5  # exactly +1 per defect, five defects
assert abs(w_v[0] + w_v[1] + w_v[2] + w_v[3]) == 30
assert abs(w_u[0] + w_u[1] + w_u[2] + w_u[3]) == 35
print("Fix2 proof: each H absorbs 4 triple flats (mu -1) into one flat "
      "(mu -3); w_3: -56 -> -51 (+1/defect x5); m_3: 35 -> 30: OK")

# ---- Try It 4: a monomial that dies (CH2) ----
X, Y = A | B, A | C
assert X in flatset and Y in flatset
assert not (X <= Y or Y <= X)
print("TryIt4: A+B and A+C are incomparable flats -> x_X * x_Y = 0 by CH2: OK")

# ---- sec 15.6: Bergman-fan counts, Ch. 8 wording ----
proper_nonempty = [f for f in flats if f and f != frozenset(E)]
assert len(proper_nonempty) == 77 and len(flats) == 79
assert 8 - 1 == 7 and 4 - 1 == 3  # quotient dim n-1, fan dim r-1
print("Sec15.6: 79 flats, 77 proper nonempty; R^8/R.1 dim 7, fan dim 3: OK")

# ---- sec 15.5: incomparable pairs are ordered (Duet note 4) ----
cnt_ord = 0
for X_ in allsubs:
    for Y_ in allsubs:
        if not (X_ <= Y_ or Y_ <= X_):
            cnt_ord += 1
assert cnt_ord == 52670 and cnt_ord % 2 == 0 and cnt_ord // 2 == 26335
print("Sec15.5: 52,670 incomparable pairs = ordered (26,335 unordered): OK")

print("\nALL PHASE-2 INDEPENDENT CHECKS PASSED")
