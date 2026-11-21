import io, contextlib
src = open("code/day12_nonuniform_chow.py").read()
ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, "day12", "exec"), ns)
Gi = ns["Gi"]
import sympy as sp
x = sp.symbols("x")
M = sp.Matrix(Gi)
cp = sp.factor(M.charpoly(x).as_expr())
print("size:", len(Gi))
print("charpoly factored over ZZ:", cp)
for val, mult_expected in [(361,4),(5776,2)]:
    # multiplicity of root
    p = M.charpoly(x).as_expr()
    m = 0
    while sp.rem(p, x-val, x) == 0:
        p = sp.quo(p, x-val, x); m += 1
    print(f"eigenvalue {val}: exact algebraic multiplicity {m} (draft claims x{mult_expected})")
print("trace check:", sum(Gi[i][i] for i in range(len(Gi))))
