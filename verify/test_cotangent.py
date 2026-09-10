"""Cotangent closure: the linear criterion R(z)R(-z) = 1, and (T*) on nonlinear fields.

In AD terms (T) is forward mode -- evaluation over D -- and (T*) is reverse mode: the
lam-component the method returns on the cotangent-lifted field is the discrete adjoint
of its own step, so (T*) says "differentiate-then-discretise = discretise-then-
differentiate" for the reverse sweep.

Two things are checked beyond the linear criterion.  Leapfrog satisfies (T*) on a
nonlinear separable field -- but only when the SPLITTING is transported by the crossed
rule {q, lam_p} | {p, lam_q}, not by the obvious {q, lam_q} | {p, lam_p}.  The reason is
that the lifted Hamiltonian <lam, X(x)> = lam_q f(p) + lam_p g(q) is separable exactly
for the crossed grouping, which is also the grouping that keeps the method explicit.  So
the shape principle of section 7 extends to T*, but the transport is dictated by the
lift and differs from the one T uses: T carries a splitting blockwise, T* crosses it.

And (T*) carries a prerequisite (T) does not.  T* is functorial only on the groupoid of
diffeomorphisms, so (T*) is statable only where Psi^X_h is invertible.  Over A that costs
nothing extra: a matrix over A is invertible iff its real part is.
"""
from fractions import Fraction as F
from poly import P
import algebra as al

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
    return ok and not sym and leapfrog_cotangent() and unit_matrices()


# ---- (T*) on a nonlinear separable field, and how the shape transports ----------

N_, M_ = 2, 4                                   # base (q,p);  cotangent (q,p,lq,lp)
MS = M_ + 1                                     # last slot is h


def _sep_field():
    """X(q,p) = (f(p), g(q)) with cubic f, g: leapfrog is explicit on it."""
    q, p = P.var(MS, 0), P.var(MS, 1)
    return [p*p*p + P.const(MS, 3)*p - P.const(MS, 2),
            q*q*q - P.const(MS, 2)*q*q + P.const(MS, 5)*q]


def _cotangent_lift(X):
    """X^{T*}(x,lam) = (X(x), -DX(x)^T lam)."""
    lam = [P.var(MS, 2), P.var(MS, 3)]
    return list(X) + [-sum((X[j].diff(i)*lam[j] for j in range(N_)), P(MS))
                      for i in range(N_)]


def _leapfrog(X, dim, part):
    """half kick on I2, drift on I1, half kick on I2; the field is only evaluated."""
    I1, I2 = part
    h = P.var(MS, M_)
    u = [P.var(MS, i) for i in range(dim)]
    def at(pt):
        im = list(pt) + [P.var(MS, k) for k in range(len(pt), M_)] + [h]
        return [al.subs(q, im, MS) for q in X]
    w = list(u)
    for i in I2: w[i] = u[i] + h*at(w)[i]*F(1, 2)
    Xw = at(w)
    for i in I1: w[i] = u[i] + h*Xw[i]
    Xw2 = at(w)
    for i in I2: w[i] = w[i] + h*Xw2[i]*F(1, 2)
    return w


def _euler(X, dim):
    h = P.var(MS, M_)
    return [P.var(MS, i) + h*X[i] for i in range(dim)]


def _isTstar(base_out, cot_out):
    """(T*) iff the x-part agrees and (D Psi^X)^T . Lambda == lam, which avoids
    inverting a polynomial matrix."""
    lam = [P.var(MS, 2), P.var(MS, 3)]
    Lam = cot_out[N_:]
    lhs = [sum((base_out[j].diff(i)*Lam[j] for j in range(N_)), P(MS)) for i in range(N_)]
    return (cot_out[:N_] == base_out) and lhs == lam


def leapfrog_cotangent():
    ok = True
    X, XT = _sep_field(), None
    XT = _cotangent_lift(X)
    print("\n  (T*) on a NONLINEAR separable cubic field, base R^2 -> cotangent R^4")
    got = _isTstar(_euler(X, N_), _euler(XT, M_))
    ok &= not got
    print(f"    explicit Euler                                     (T*): {got}  expected False")
    base = _leapfrog(X, N_, ([0], [1]))
    for lbl, part, exp in [("crossed          {q,lam_p}|{p,lam_q}", ([0, 3], [1, 2]), True),
                           ("cotangent-lifted {q,lam_q}|{p,lam_p}", ([0, 2], [1, 3]), False),
                           ("base|fibre       {q,p}|{lam_q,lam_p}", ([0, 1], [2, 3]), False)]:
        got = _isTstar(base, _leapfrog(XT, M_, part))
        ok &= (got == exp)
        print(f"    leapfrog, {lbl}   (T*): {str(got):5s} expected {exp}"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    print("    -> leapfrog IS cotangent-closed, but the splitting must be transported by")
    print("       the crossed rule: <lam, X> = lam_q f(p) + lam_p g(q) is separable exactly")
    print("       for {q,lam_p}|{p,lam_q}, which is also where the method stays explicit.")
    print("       T carries a splitting blockwise (test_jettransport); T* crosses it.")
    return ok


def unit_matrices():
    """(T*) needs Psi^X_h invertible, since T* is functorial only on diffeomorphisms.
    Over A that costs nothing extra: invertible iff the real part is."""
    ok = True
    print("\n  the (T*) prerequisite: invertibility over A reduces to the real part")
    for A in (al.D2, al.J3, al.R2s):
        for lbl, blk in (("real part invertible", [[F(2), F(1)], [F(0), F(3)]]),
                         ("real part singular  ", [[F(1), F(1)], [F(1), F(1)]])):
            # 2x2 matrix over A: entry (i,j) is blk[i][j]*1_A + (nilpotent junk)
            det_real = blk[0][0]*blk[1][1] - blk[0][1]*blk[1][0]
            a = [det_real * c for c in A.one]                 # det, as an element of A
            inv = _solvable([[sum(a[x]*A.c[x][y][g] for x in range(A.N))
                              for y in range(A.N)] for g in range(A.N)],
                            [A.one[g] for g in range(A.N)], A.N)
            exp = (det_real != 0)
            ok &= (inv == exp)
            print(f"    {A.name:26s} {lbl}  det a unit in A: {str(inv):5s} expected {exp}")
    print("    -> so once (T*) is statable over R it is statable over every A.")
    return ok


def _solvable(M, b, w):
    R = [M[i][:] + [b[i]] for i in range(len(M))]
    r = 0
    for c in range(w):
        p_ = next((i for i in range(r, len(R)) if R[i][c]), None)
        if p_ is None: continue
        R[r], R[p_] = R[p_], R[r]; pv = R[r][c]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c]/pv; R[i] = [R[i][j] - f*R[r][j] for j in range(w + 1)]
        r += 1
    return not any(all(R[i][c] == 0 for c in range(w)) and R[i][w] != 0 for i in range(len(R)))

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
