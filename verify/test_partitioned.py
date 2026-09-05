"""PRK satisfies (T) for the D-block lift only, and is not affine equivariant."""
from fractions import Fraction as F
from dual import *

def prk(Fld, u, h, blocks):
    w = list(u)
    for blk in blocks:
        Fv = Fld(w)
        for i in blk:
            w[i] = u[i] + h * Fv[i]        # original u[i], field at CURRENT w
    return w

def main():
    o = prk(X, [Dual(PT[0], TG[0]), Dual(PT[1], TG[1])], D(H), [[0], [1]])
    L = ([x.a for x in o], [x.b for x in o])
    ok = True
    print("same PRK on R^4, three lifts of the partition:")
    for label, blocks, expect in [("D-submodule (q,qd)|(p,pd)", [[0, 2], [1, 3]], True),
                                  ("real base|fibre (q,p)|(qd,pd)", [[0, 1], [2, 3]], False),
                                  ("crossed (q,pd)|(p,qd)", [[0, 3], [1, 2]], False)]:
        r = prk(TX, list(PT) + list(TG), H, blocks)
        ok &= report(label, L, (r[:2], r[2:]), expect)
    A = [[F(1), F(1)], [F(0), F(1)]]; Ai = [[F(1), F(-1)], [F(0), F(1)]]
    def Xpush(u): return matvec(A, X(matvec(Ai, u)))
    eq = matvec(A, prk(X, PT, H, [[0], [1]])) == prk(Xpush, matvec(A, PT), H, [[0], [1]])
    print(f"  affine equivariant under a block-mixing A?         {eq}"
          f"   -> (T)-natural yet NOT a B-series: {not eq}")
    return ok and not eq

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
