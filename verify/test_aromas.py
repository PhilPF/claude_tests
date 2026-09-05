"""Aromas scale by 2^alpha; genuine trees are natural."""
from fractions import Fraction as F
from dual import *

def d2XXX(f, u):
    """D^2 f(u)(f(u), f(u)); exact for quadratic f via f(u+w)+f(u-w)-2f(u)."""
    n = len(u); fx = f(u)
    f1 = f([u[i] + fx[i] for i in range(n)])
    fm = f([u[i] - fx[i] for i in range(n)])
    return [f1[i] + fm[i] - 2 * fx[i] for i in range(n)]

def main():
    p = list(PT) + list(TG)
    print("aroma scaling under the lift: (div X)^a")
    ok = True
    for a in range(1, 4):
        lo, hi = div(JX(p[:2])) ** a, div(JTX(p)) ** a
        r = F(hi, lo); ok &= (r == 2 ** a)
        print(f"  a={a}: ratio = {r}   predicted 2^{a} = {2 ** a}   {'OK' if r == 2**a else 'FAIL'}")
    hi, lo = d2XXX(TX, p), d2XXX(X, p[:2])
    du = d2XXX(X, [Dual(p[0], p[2]), Dual(p[1], p[3])])
    T = ([x.a for x in du], [x.b for x in du])
    nat = (hi[:2] == T[0] and hi[2:] == T[1]); ok &= nat
    print(f"\n  tree D^2X(X,X) natural: F_TX(tau) == T(F_X(tau))?  {nat}")
    return ok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
