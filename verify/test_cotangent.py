"""Cotangent closure on linear fields is R(z)R(-z) = 1."""
from fractions import Fraction as F

def pmul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b): c[i + j] += ai * bj
    return c
def refl(a): return [c * (-1) ** i for i, c in enumerate(a)]

def main():
    ok = True
    print("  R(z)R(-z) - 1   (explicit methods; 0 iff cotangent closure holds)")
    for name, R in [("explicit Euler", [F(1), F(1)]), ("Heun", [F(1), F(1), F(1, 2)]),
                    ("RK3", [F(1), F(1), F(1, 2), F(1, 6)]),
                    ("RK4", [F(1), F(1), F(1, 2), F(1, 6), F(1, 24)])]:
        p = pmul(R, refl(R)); s = len(R) - 1
        tail = {i: c for i, c in enumerate(p) if (c != 0 and i > 0) or (i == 0 and c != 1)}
        top, pred = p[2 * s], (-1) ** s * R[s] * R[s]
        ok &= (top == pred != 0)
        print(f"    {name:16s} " + " + ".join(f"{c}z^{i}" for i, c in tail.items())
              + f"      top z^{2*s} = {top} = (-1)^s c_s^2  (nonzero)")
    num, den = [F(1), F(1, 2)], [F(1), F(-1, 2)]
    mid = pmul(num, refl(num)) == pmul(den, refl(den)); ok &= mid
    print(f"    implicit midpoint  0 exactly -> holds: {mid}")
    A = [[F(0), F(0)], [F(1, 2), F(1, 2)]]; b = [F(1, 2), F(1, 2)]
    sym = all(b[i] * A[i][j] + b[j] * A[j][i] == b[i] * b[j] for i in range(2) for j in range(2))
    print(f"    trapezoidal rule passes the linear test but is symplectic? {sym}"
          f"  -> linear test strictly weaker: {not sym}")
    return ok and not sym

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
