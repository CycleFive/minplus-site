#!/usr/bin/env python3
"""Ch.14 Phase 2 (pre-run, proactivity 2026-11-20 evening):
assert the draft's four Try It answers individually."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch14_chapter_program_draft import uniform_flats, m_line_rank4, bjorner_count

# Try It 1: log-concavity of (1,5,9,7) by hand
assert 5**2 >= 1*9 and 9**2 >= 5*7
print("TryIt1: 25>=9 and 81>=35: OK")

# Try It 2: ledger subtraction
cost = [u-m for u, m in zip((1,5,10,10), (1,5,9,7))]
assert cost == [0,0,1,3]
print("TryIt2: cost (0,0,1,3): OK")

# Try It 3: {0,1,2} and {3,4} are incomparable flats of M_line^(4) -> CH2 kills product
E, flats, rk, _ = m_line_rank4()
A, B = frozenset({0,1,2}), frozenset({3,4})
assert A in flats and B in flats and not (A <= B or B <= A)
print("TryIt3: both flats, incomparable, CH2 product = 0: OK")

# Try It 4: {4} < {3,4} < {2,3,4} in U(4,5) is initial descending, counted by Bjorner m3=4
E5 = [0,1,2,3,4]; flats5, rk5 = uniform_flats(E5, 4)
chain = (frozenset({4}), frozenset({3,4}), frozenset({2,3,4}))
assert all(c in flats5 for c in chain) and chain[0] < chain[1] < chain[2]
assert [rk5(c) for c in chain] == [1,2,3]
mins = [min(c) for c in chain]
assert mins == [4,3,2] and all(mins[i] > mins[i+1] for i in range(2)) and mins[-1] > 0
counts, examples = bjorner_count(E5, flats5, rk5, 3)
assert counts == [1,4,6,4], counts
assert tuple(tuple(sorted(c)) for c in chain) in examples[3], examples[3]
print(f"TryIt4: flag initial descending, among Bjorner examples for m3=4: OK")

print("\nALL FOUR TRY IT ANSWERS VERIFIED")
