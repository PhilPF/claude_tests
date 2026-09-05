"""A = R[e]/(e^3): aromas scale by dim_R A = 3; RK4 is T^A-natural."""
from fractions import Fraction as F
from dual import *

def pack(u):   return [Jet(u[0], u[2], u[4]), Jet(u[1], u[3], u[5])]
def unpack(w): return [w[0].c[0], w[1].c[0], w[0].c[1], w[1].c[1], w[0].c[2], w[1].c[2]]
def TAX(u):    return unpack(X(pack(u)))

def jac(f, u, n):
    rows = [[None] * n for _ in range(n)]
    for j in range(n):
        out = f([Dual(u[k], 1 if k == j else 0) for k in range(n)])
        for i in range(n): rows[i][j] = out[i].b
    return rows

def main():
    u6 = [F(2, 3), F(-5, 4), F(3, 5), F(7, 2), F(-1, 6), F(4, 9)]
    dA, d1 = div(jac(TAX, u6, 6)), div(JX(u6[:2]))
    r = F(dA, d1); ok = (r == 3)
    print(f"  div(T^A X)/div(X) = {r}   predicted dim_R A = 3   {'OK' if ok else 'FAIL'}")
    lhs = unpack(rk(X, pack(u6), J(H), *RK4))
    rhs = rk(TAX, u6, H, *RK4)
    ok &= (lhs == rhs)
    print(f"  RK4 is T^A-natural for dim_R A = 3?  {lhs == rhs}")
    return ok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
