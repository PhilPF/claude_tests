"""Claims that survived the final adjudication round, checked here directly.

1. The lift multiplier counts LOOPS, not aromas.
2. The liana constant c_A beyond the dual numbers.
3. A polynomial closed method with no equivariance at all (so closure alone
   cannot be a definability condition).
4. Gradient fields are not stable under the tangent lift.
"""
from fractions import Fraction as F
from dual import *


def main():
    ok = True

    # --- 1. multiplier counts loops, not aromas -------------------------------
    # evaluate on the zero section (v = 0), where the base multiplier is visible
    z = list(PT) + [F(0), F(0)]
    JT, Jb = JTX(z), JX(PT)
    frob = lambda M: sum(M[i][j] * M[i][j] for i in range(len(M)) for j in range(len(M)))
    loop = F(div(JT), div(Jb))                                   # div X : loop
    liana_stolon = F(frob(JT), frob(Jb))                         # |DX|_F^2 : liana+stolon, NO loop
    fx, tfx = X(PT), TX(z)
    aroma = F(sum(tfx[i] * tfx[i] for i in range(4)), sum(fx[i] * fx[i] for i in range(2)))
    print("  multipliers on the zero section:")
    print(f"    div X            (one loop)                 x{loop}")
    print(f"    |DX|_F^2         (liana+stolon, NO loop)    x{liana_stolon}")
    print(f"    <X,X>            (aroma, NO loop)           x{aroma}")
    good = (loop == 2 and liana_stolon == 2 and aroma == 1)
    ok &= good
    print(f"    -> multiplier is not a function of aroma count: {good}")

    # --- 2. the liana constant c_A ------------------------------------------
    print("\n  liana defect over A = R[e]/(e^r):  predicted (c_A - 1) * T^A(Delta X)")
    for r in (2, 3, 4):
        n, N = 2, 2 * r
        u = [F(2, 3), F(-5, 4), F(3, 5), F(7, 2), F(-1, 6), F(4, 9), F(5, 3), F(-2, 7)][:r * n]
        TA = lambda w: unpackA(Xc(packA(w, r, n)), r, n)
        got = lap(TA, u, r * n)
        pk = packA(u, r, n)
        TAlap = unpackA([Trunc(r, 8) * pk[0], Trunc(r, 4) * pk[1]], r, n)   # Delta Xc = (8x1, 4x2)
        cA = Trunc(r)
        for k in range(r):
            e2k = [0] * r
            if 2 * k < r: e2k[2 * k] = 1; cA = cA + Trunc(r, *e2k)
        pred = unpackA([(cA - 1) * Trunc(r, 8) * pk[0], (cA - 1) * Trunc(r, 4) * pk[1]], r, n)
        defect = [got[i] - TAlap[i] for i in range(N)]
        m = defect == pred; ok &= m
        print(f"    r={r}: c_A = {cA}   defect matches (c_A-1)*T^A(Delta X): {m}"
              + ("   (c_A = 1: liana is free)" if r == 2 else ""))

    # --- 3. a polynomial closed method with NO equivariance -------------------
    # Psi(u) = u + h f + h^2 (f^1)^2 f, with the scalar squared in the RING.
    def m_ring(f, u, h, n):
        fx = f(u); s = fx[0] * fx[0]
        return [u[i] + h * fx[i] + h * h * s * fx[i] for i in range(n)]

    def m_Dstruct(f, u, h, n):        # on R^4 = D^2: square the first D-coordinate
        fx = f(u)
        s = Dual(fx[0], fx[2]) * Dual(fx[0], fx[2])          # (f^1 + e f^3)^2 in D
        out = []
        for j in range(2):                                    # D-scalar times D-vector
            p = s * Dual(fx[j], fx[j + 2])
            out.append((j, p))
        base = [u[j] + h * fx[j] + h * h * out[j][1].a for j in range(2)]
        fib = [u[j + 2] + h * fx[j + 2] + h * h * out[j][1].b for j in range(2)]
        return base + fib

    o = m_ring(X, [Dual(PT[0], TG[0]), Dual(PT[1], TG[1])], D(H), 2)
    L = ([x.a for x in o], [x.b for x in o])
    r4 = m_Dstruct(TX, list(PT) + list(TG), H, 4)
    r4n = m_ring(TX, list(PT) + list(TG), H, 4)
    print("\n  method  u + hX + h^2 (X^1)^2 X  (uses a preferred coordinate):")
    ok &= report("D-structured scalar (ring-valued)", L, (r4[:2], r4[2:]), True)
    ok &= report("naive real scalar on R^4", L, (r4n[:2], r4n[2:]), False)
    A = [[F(1), F(1)], [F(0), F(1)]]; Ai = [[F(1), F(-1)], [F(0), F(1)]]
    def Xp(u): return matvec(A, X(matvec(Ai, u)))
    eq = matvec(A, m_ring(X, PT, H, 2)) == m_ring(Xp, matvec(A, PT), H, 2)
    print(f"    affine equivariant? {eq}   -> polynomial, closed, NO equivariance: {not eq}")
    ok &= not eq

    # --- 4. gradient fields are not stable under the lift ---------------------
    # X = grad V with V = x1^3 x2 + x2^2  =>  X = (3x1^2 x2, x1^3 + 2x2)
    def Xg(u): x1, x2 = u; return [3 * x1 * x1 * x2, x1 * x1 * x1 + 2 * x2]
    def TXg(u):
        x1, x2, v1, v2 = u
        return [3 * x1 * x1 * x2, x1 * x1 * x1 + 2 * x2,
                6 * x1 * x2 * v1 + 3 * x1 * x1 * v2, 3 * x1 * x1 * v1 + 2 * v2]
    def jac(f, u, n):
        R = [[None] * n for _ in range(n)]
        for j in range(n):
            out = f([Dual(u[k], 1 if k == j else 0) for k in range(n)])
            for i in range(n): R[i][j] = out[i].b
        return R
    Jg, Jt = jac(Xg, PT, 2), jac(TXg, list(PT) + list(TG), 4)
    sym = lambda M: all(M[i][j] == M[j][i] for i in range(len(M)) for j in range(len(M)))
    print(f"\n  X = grad V: Jacobian symmetric? {sym(Jg)}    T X: symmetric? {sym(Jt)}")
    print(f"    -> the gradient class is NOT stable under T, so closure is not"
          f" well posed on it: {sym(Jg) and not sym(Jt)}")
    ok &= sym(Jg) and not sym(Jt)
    return ok


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
