"""The v-degree of the closure defect, split into the two halves of (T).

Writing v = t*v0 and expanding the defect
    E = F_{2n}(gamma)(TX) - T(F_n(gamma)(X))
in powers of t separates the obstructions by v-degree. Two things come out:

  * Theta(F(gamma)) = 2^{loops} F(gamma) is FALSE in the exotic class:
    |DX|_F^2 X carries NO loop yet has base-block multiplier 2.
  * In the BASE block, pure loops obstruct at v-degree 0 and pure stolons only
    at v-degree 2 -- disjoint, so they cannot cancel against one another.

Differentiation here is by exact polynomial interpolation, so the differentials
can be written uniformly in the dimension with no nested dual numbers.
"""
from fractions import Fraction as F


def tcoef(fn, u, w, K):
    """Coefficients of s^0..s^K in fn(u + s*w); exact for polynomial fn."""
    pts = [F(k) for k in range(K + 1)]
    vals = [fn([u[i] + s * w[i] for i in range(len(u))]) for s in pts]
    out = []
    for comp in range(len(vals[0])):
        A = [[pts[k] ** j for j in range(K + 1)] + [vals[k][comp]] for k in range(K + 1)]
        for c in range(K + 1):                       # exact Gauss-Jordan
            p = next(r for r in range(c, K + 1) if A[r][c] != 0)
            A[c], A[p] = A[p], A[c]
            pv = A[c][c]; A[c] = [x / pv for x in A[c]]
            for r in range(K + 1):
                if r != c and A[r][c] != 0:
                    f = A[r][c]; A[r] = [A[r][j] - f * A[c][j] for j in range(K + 2)]
        out.append([A[j][K + 1] for j in range(K + 1)])
    return out


def d1(fn, u, w, K=10): return [c[1] for c in tcoef(fn, u, w, K)]
def d2(fn, u, w, K=10): return [2 * c[2] for c in tcoef(fn, u, w, K)]
def basis(m, j): return [F(1) if i == j else F(0) for i in range(m)]
def jac(fn, u, m):
    cols = [d1(fn, u, basis(m, j)) for j in range(m)]
    return [[cols[j][i] for j in range(m)] for i in range(m)]


def g_loop(f, u, m):                      # div(X) X        one loop
    J = jac(f, u, m); dv = sum(J[i][i] for i in range(m)); fx = f(u)
    return [dv * fx[i] for i in range(m)]
def g_liana(f, u, m):                     # Delta X         liana only
    return [sum(d2(f, u, basis(m, j))[i] for j in range(m)) for i in range(m)]
def g_stolon(f, u, m):                    # <X,X> X         stolon, no loop
    fx = f(u); s = sum(fx[i] * fx[i] for i in range(m))
    return [s * fx[i] for i in range(m)]
def g_ls(f, u, m):                        # |DX|_F^2 X      liana+stolon, NO loop
    J = jac(f, u, m); fx = f(u)
    s = sum(J[i][j] * J[i][j] for i in range(m) for j in range(m))
    return [s * fx[i] for i in range(m)]


def Xc(u):
    x1, x2 = u
    return [x1 * x1 * x1 + x1 * x2 * x2, x2 * x2 * x2 - x1 * x1 * x2]
def TXc(u):
    x1, x2, v1, v2 = u
    a = 3 * x1 * x1 + x2 * x2; b = 2 * x1 * x2
    c = -2 * x1 * x2; d = 3 * x2 * x2 - x1 * x1
    return [x1 * x1 * x1 + x1 * x2 * x2, x2 * x2 * x2 - x1 * x1 * x2,
            a * v1 + b * v2, c * v1 + d * v2]


X0, V0, K = [F(2, 3), F(-5, 4)], [F(3, 5), F(7, 2)], 8


def degrees(G):
    """Non-zero v-degrees of the defect, as (base-block, fibre-block)."""
    A = tcoef(lambda p: G(TXc, p, 4), X0 + [F(0), F(0)], [F(0), F(0)] + V0, K)
    base, lin = G(Xc, X0, 2), d1(lambda p: G(Xc, p, 2), X0, V0)
    B = [[base[i]] + [F(0)] * K for i in range(2)] \
        + [[F(0), lin[i]] + [F(0)] * (K - 1) for i in range(2)]
    E = [[A[i][k] - B[i][k] for k in range(K + 1)] for i in range(4)]
    pick = lambda idx: sorted({k for i in idx for k in range(K + 1) if E[i][k] != 0})
    return pick((0, 1)), pick((2, 3))


def main():
    ok = True
    print(f"  {'differential':22s} {'structure':23s} {'BASE':10s} FIBRE")
    exp = {"div(X) X": ([0], [1]), "Delta X": ([], []),
           "<X,X> X": ([2], [1, 3]), "|DX|_F^2 X": ([0, 2], [1, 3])}
    for name, G, note in [("div(X) X", g_loop, "one loop"),
                          ("Delta X", g_liana, "liana only"),
                          ("<X,X> X", g_stolon, "stolon, no loop"),
                          ("|DX|_F^2 X", g_ls, "liana+stolon, NO loop")]:
        b, f = degrees(G)
        good = ([b, f] == [exp[name][0], exp[name][1]]); ok &= good
        print(f"  {name:22s} {note:23s} {str(b or 'zero'):10s} {f or 'zero'}"
              + ("" if good else "   <-- UNEXPECTED"))
    bl, _ = degrees(g_loop); bs, _ = degrees(g_stolon); bls, _ = degrees(g_ls)
    sep = (0 in bl and 0 not in bs)
    refute = (0 in bls)                     # no loop, yet base v-degree 0
    print(f"\n  loops (base deg {bl}) and stolons (base deg {bs}) are disjoint at"
          f" v-degree 0: {sep}")
    print(f"  |DX|_F^2 X has NO loop yet obstructs at base v-degree 0: {refute}")
    print(f"    -> Theta acts by N_0(gamma), a count of zero-cost level assignments,"
          f" NOT by 2^{{loops}}")
    return ok and sep and refute


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
