"""(T) for standard methods; failure for aromatic and basis-dependent ones."""
from fractions import Fraction as F
from dual import *

def m_rk4(f, u, h, n):  return rk(f, u, h, *RK4)
def m_heun(f, u, h, n): return rk(f, u, h, *HEUN)

def m_taylor2(f, u, h, n):                     # x + hX + (h^2/2) X'X : B-series, not RK-defined
    Jac = JX if n == 2 else JTX
    fx = f(u); jf = matvec(Jac(u), fx)
    return [u[d] + h * fx[d] + F(1, 2) * h * h * jf[d] for d in range(n)]

def m_aromatic(f, u, h, n):                    # x + hX + h^2 div(X) X : AROMATIC
    Jac = JX if n == 2 else JTX
    fx = f(u); dv = div(Jac(u))
    return [u[d] + h * fx[d] + h * h * dv * fx[d] for d in range(n)]

def m_basis(f, u, h, n):                       # x + hX + h^2 <e1,X> e1 : basis-dependent
    fx = f(u)
    return [u[d] + h * fx[d] + (h * h * fx[0] if d == 0 else 0) for d in range(n)]

def main():
    print("(T): T(Psi^X) == Psi^(TX), exact over Q")
    ok = True
    for name, m, expect in [("RK4", m_rk4, True), ("Heun (RK2)", m_heun, True),
                            ("Taylor-2 (B-series, not RK)", m_taylor2, True),
                            ("aromatic  x+hX+h^2 div(X)X", m_aromatic, False),
                            ("basis-dependent <e1,X>e1", m_basis, False)]:
        lhs, rhs = closure(m, X, TX)
        ok &= report(name, lhs, rhs, expect)
    lhs, rhs = closure(m_basis, X, TX)
    print(f"  ...basis-dependent splits (T): base {'OK' if lhs[0]==rhs[0] else 'FAIL'}, "
          f"fibre {'OK' if lhs[1]==rhs[1] else 'FAIL'}  -> (T-base) and (T-fibre) independent")
    return ok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
