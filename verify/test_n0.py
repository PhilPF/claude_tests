"""Stress-test of the N_0 law on forests with several interacting decorations.

N_0(gamma) counts families of vertex-disjoint GENERALISED CYCLES -- closed walks
that traverse arrows head-to-tail and reverse direction at a liana (two lower
ends) or a stolon (two upper ends). The base multiplier on the zero section
should equal N_0.

The decisive pair: the theta graph and the minimal residual forest carry the
SAME decorations (two lianas, one stolon) yet have N_0 = 3 and N_0 = 1. If both
hold, N_0 is structural rather than a tally of decorations -- which is what
distinguishes it from the 2^aromas and 2^loops laws that preceded it and failed.
"""
from fractions import Fraction as F
from test_vdegree import tcoef, d1, basis, Xc, TXc


def jac(f, u, m):
    cols = [d1(f, u, basis(m, j)) for j in range(m)]
    return [[cols[j][i] for j in range(m)] for i in range(m)]          # J[i][j] = d_j f^i


def hess(f, u, m):
    """H[i][j][k] = d_j d_k f^i, exactly, by polarisation."""
    q2 = lambda w: [2 * c[2] for c in tcoef(f, u, w, 10)]
    H = [[[None] * m for _ in range(m)] for _ in range(m)]
    for j in range(m):
        for k in range(m):
            e = [basis(m, j)[q] + basis(m, k)[q] for q in range(m)]
            a, b, c = q2(e), q2(basis(m, j)), q2(basis(m, k))
            for i in range(m):
                H[i][j][k] = (a[i] - b[i] - c[i]) / 2
    return H


def g_theta(f, u, m):          # |D^2X|_F^2 X          two lianas + one stolon
    H = hess(f, u, m); fx = f(u)
    s = sum(H[a][j][k] * H[a][j][k] for a in range(m) for j in range(m) for k in range(m))
    return [s * fx[i] for i in range(m)]


def g_residual(f, u, m):       # d_jk X^i d_j X^a d_k X^a   two lianas + one stolon
    H, J = hess(f, u, m), jac(f, u, m)
    return [sum(H[i][j][k] * J[a][j] * J[a][k]
                for j in range(m) for k in range(m) for a in range(m)) for i in range(m)]


def g_2stolon(f, u, m):        # <X,X>^2 X             two stolons, no liana
    fx = f(u); s = sum(fx[i] * fx[i] for i in range(m))
    return [s * s * fx[i] for i in range(m)]


def g_lianastolon(f, u, m):    # |DX|_F^2 X            one liana + one stolon
    J = jac(f, u, m); fx = f(u)
    s = sum(J[i][j] * J[i][j] for i in range(m) for j in range(m))
    return [s * fx[i] for i in range(m)]


def g_mixed(f, u, m):          # div(X) |DX|_F^2 X     loop + liana + stolon, disjoint
    J = jac(f, u, m); fx = f(u)
    dv = sum(J[i][i] for i in range(m))
    s = sum(J[i][j] * J[i][j] for i in range(m) for j in range(m))
    return [dv * s * fx[i] for i in range(m)]


CASES = [("|D^2X|_F^2 X  (theta)",    g_theta,       "2 lianas + 1 stolon",      3),
         ("d_jk X^i d_j X^a d_k X^a", g_residual,    "2 lianas + 1 stolon",      1),
         ("<X,X>^2 X",                g_2stolon,     "2 stolons, no liana",      1),
         ("|DX|_F^2 X",               g_lianastolon, "1 liana + 1 stolon",       2),
         ("div(X) |DX|_F^2 X",        g_mixed,       "loop + liana + stolon",    4)]


def main():
    x0 = [F(2, 3), F(-5, 4)]
    z0 = x0 + [F(0), F(0)]                                   # zero section, v = 0
    ok = True
    print(f"  {'differential':28s} {'structure':24s} {'measured':>9s}  N_0")
    for name, G, struct, claim in CASES:
        lo, hi = G(Xc, x0, 2), G(TXc, z0, 4)
        nz = [i for i in range(2) if lo[i] != 0]
        rs = {F(hi[i], lo[i]) for i in nz}
        r = rs.pop() if len(rs) == 1 else rs
        good = (r == claim); ok &= good
        print(f"  {name:28s} {struct:24s} {str(r):>9s}  {claim}" + ("" if good else "  <-- MISMATCH"))
    print("\n  the decisive pair: identical decorations (2 lianas + 1 stolon), multipliers 3 and 1")
    print("  -> N_0 is structural, not a tally of decorations")
    return ok


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
