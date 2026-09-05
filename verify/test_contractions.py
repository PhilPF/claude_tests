"""The trichotomy: loops die, LIANAS (Laplacian) survive over D only, STOLONS die.

Terminology follows Laurent & Munthe-Kaas (arXiv:2305.10993, Def. 2.7): a liana
identifies two ARROWS (double derivation = Laplacian); a stolon identifies two
NODES (double evaluation = scalar product).
"""
from fractions import Fraction as F
from dual import *

def m_liana(f, u, h, n):   # x + hX + h^2 Delta X
    fx = f(u); lp = lap(f, u, n)
    return [u[i] + h * fx[i] + h * h * lp[i] for i in range(n)]

def m_stolon(f, u, h, n):  # x + hX + h^2 |X|^2 X
    fx = f(u); s = sum((fx[i] * fx[i] for i in range(n)), 0)
    return [u[i] + h * fx[i] + h * h * s * fx[i] for i in range(n)]

def pack(u):   return [Jet(u[0], u[2], u[4]), Jet(u[1], u[3], u[5])]
def unpack(w): return [w[0].c[0], w[1].c[0], w[0].c[1], w[1].c[1], w[0].c[2], w[1].c[2]]
def TAXc(u):   return unpack(Xc(pack(u)))

def main():
    ok = True
    lapTX = lap(TXc, list(PT) + list(TG), 4)
    o = lap(Xc, [Dual(PT[0], TG[0]), Dual(PT[1], TG[1])], 2)
    T = ([x.a for x in o], [x.b for x in o])
    nat = (lapTX[:2] == T[0] and lapTX[2:] == T[1]); ok &= nat
    print(f"  Laplacian natural:  Delta_R4(TX) == T(Delta X)?      {nat}")
    lhs, rhs = closure(m_liana, Xc, TXc)
    ok &= report("LIANA method  x+hX+h^2 Delta X", lhs, rhs, True)
    A = [[F(1), F(1)], [F(0), F(1)]]; Ai = [[F(1), F(-1)], [F(0), F(1)]]
    def Xp(u): return matvec(A, Xc(matvec(Ai, u)))
    eq = matvec(A, m_liana(Xc, PT, H, 2)) == m_liana(Xp, matvec(A, PT), H, 2)
    print(f"  ...affine equivariant?  {eq}   -> (T)-natural yet NOT a B-series: {not eq}")
    ok &= not eq
    lhs, rhs = closure(m_stolon, Xc, TXc)
    ok &= report("STOLON method x+hX+h^2 |X|^2 X (aroma-free)", lhs, rhs, False)

    # the liana survives only because eps^2 = 0: over R[e]/(e^3) it fails by exactly Delta X
    u6 = [F(2, 3), F(-5, 4), F(3, 5), F(7, 2), F(-1, 6), F(4, 9)]
    lapTA = lap(TAXc, u6, 6)
    lapX_A = unpack([Jet(8, 0, 0) * pack(u6)[0], Jet(4, 0, 0) * pack(u6)[1]])  # Delta Xc = (8x1,4x2)
    d = [lapTA[i] - lapX_A[i] for i in range(6)]
    good = (d[0:2] == [0, 0] and d[2:4] == [0, 0] and d[4:6] == [8 * u6[0], 4 * u6[1]])
    ok &= good
    print(f"  liana over R[e]/(e^3): defect in eps^2 block = Delta X?  {good}"
          f"   -> survives (T) only by the accident eps^2=0")
    return ok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
