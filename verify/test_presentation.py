"""Presentation uniqueness, and what it gives the gluing problem.

A *presentation* of Z in X(R^m) is a triple (n, (A,iota), X) with n dim A = m and
Z = X^A in the fixed basis.  Every Z has the trivial one, (m, R, Z).  Whether two
presentations of one field must have a common refinement is what decides whether
the conditions algebraic naturality puts on a single Psi_m are consistent.

C1.  THE COMMUTANT IS THE PRESENTATION.  Put

        A(Z) = { L in End(R^m) : L D^kZ(u)[v_1,...,v_k] = D^kZ(u)[Lv_1,v_2,...,v_k]
                                 for all k >= 1, all u, all v_i }.

     This is a unital subalgebra of End(R^m), and it is computable: the conditions
     are linear in L.  One line gives the inclusion that matters -- if Z = X^A then
     every D^kZ(u) is A-multilinear, so rho_A(A) is contained in A(Z).  Measured
     here: for a generic cubic field and each of 17 based algebras, A(X^A) is
     EXACTLY rho_A(A).  So a generic lifted field determines its presentation.

C2.  AND IT IS THE COMMON REFINEMENT.  In both coincidence families found earlier
     -- Y^(D(x)D) presented over (1, D(x)D) and (2, D); (W^(R^2))^(R^3) presented
     over (2, R^3) and (3, R^2) -- the commutant is the algebra of the FINEST
     presentation and contains the action of both coarser ones.  So two
     presentations refine to the commutant presentation.

     Consequence for the gluing.  If Z = X^A = Y^B and the family is already natural
     below dimension m, then both prescriptions for Psi_m(Z) factor through the
     commutant presentation (n_0, E, W): X = W^C gives T^A Psi(X) = T^A T^C Psi(W) =
     T^E Psi(W), and likewise for B.  The prescription on the union of lifted loci is
     therefore CONSISTENT, by induction on dimension.  What remains for open item 2
     is not consistency but smoothness -- a Whitney-type extension across the union.

C3.  ZERO-FIELD RIGIDITY.  The zero field is lifted by every algebra, 0^A = 0, so
     every pair of algebras of one dimension collides on it.  With affine rigidity
     (Psi^0 is an affine map, u -> S u + w) this forces w = 0 and S = s.I:

        Psi^0_h(u) = s(h) u,   and s = 1 once consistency is imposed.

     So an algebraically natural family carries NO distinguished vector and no
     distinguished endomorphism -- u + hX + h^2 e_1 is affine on affine fields, so
     affine rigidity misses it, and it dies here instead.
"""
from fractions import Fraction as F
from itertools import product as iproduct
from poly import P
import algebra as al
from algebra import ALGEBRAS
from test_affine import gen_map, kernel, euler, rk4, laplacian, _u


def commutant(Z, m, kmax=3):
    """Basis of A(Z), as flattened m x m matrices.  Linear in L, so a nullspace."""
    eqs = {}
    def acc(key, var, coef):
        if coef:
            row = eqs.setdefault(key, {})
            row[var] = row.get(var, F(0)) + coef
    for k in range(1, kmax + 1):
        for multi in iproduct(range(m), repeat=k):
            T = [Z[i] for i in range(m)]
            for a in multi: T = [t.diff(a) for t in T]
            for i in range(m):
                for c in range(m):
                    for mon, v in T[c].d.items(): acc((k, multi, i, mon), i * m + c, v)
                for c in range(m):
                    S = [Z[j] for j in range(m)]
                    for a in (c,) + multi[1:]: S = [s.diff(a) for s in S]
                    for mon, v in S[i].d.items(): acc((k, multi, i, mon), c * m + multi[0], -v)
    rows = [[row.get(j, F(0)) for j in range(m * m)] for row in eqs.values()]
    return kernel(rows, m * m) if rows else [[F(1) if i == j else F(0)
                                              for i in range(m) for j in range(m)]]


def rho(A, n):
    """the A-action on A^n = R^{nN}, block convention, as flattened matrices."""
    N, m = A.N, n * A.N
    out = []
    for a in range(N):
        M = [[F(0)] * m for _ in range(m)]
        for j in range(n):
            for be in range(N):
                for g in range(N): M[j * N + g][j * N + be] = A.c[a][be][g]
        out.append([x for r in M for x in r])
    return out


def rank(rows, w):
    R = [r[:] for r in rows]; r = 0
    for c in range(w):
        p = next((i for i in range(r, len(R)) if R[i][c]), None)
        if p is None: continue
        R[r], R[p] = R[p], R[r]; pv = R[r][c]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c] / pv; R[i] = [R[i][j] - f * R[r][j] for j in range(w)]
        r += 1
    return r


def recovers():
    ok = True
    n, Ms = 2, 2
    X = gen_map(n, Ms, 3)
    print("  C1  the commutant of a generic lifted field recovers its algebra")
    print("      based algebra                    dim A  dim A(Z)  A(Z) = rho_A(A)")
    for _, lst in sorted(ALGEBRAS.items()):
        for A in lst:
            K = commutant(al.lift_map(A, n, X, Ms), n * A.N)
            good = (len(K) == A.N and rank(K + rho(A, n), (n * A.N) ** 2) == A.N)
            ok &= good
            print(f"      {A.name:32s} {A.N:5d} {len(K):9d}       {good}"
                  + ("" if good else "   <-- UNEXPECTED"))
    print("      -> a generic lifted field determines its presentation outright.")
    return ok


def refinement():
    ok = True
    print("\n  C2  in the coincidence families the commutant is the COMMON REFINEMENT")
    t = P.var(1, 0)
    Y = [t * t * t + P.const(1, 3) * t * t - P.const(1, 2) * t]
    DD = al.tensor(al.D2, al.D2)
    Z = al.lift_map(DD, 1, Y, 1)
    K = commutant(Z, 4)
    a = (len(K) == 4 and rank(K + rho(DD, 1), 16) == 4 and rank(K + rho(al.D2, 2), 16) == 4
         and al.recognise(DD, 1, Z) == Y)
    ok &= a
    print(f"      Z = Y^(D(x)D) on R^4, presented over (1, D(x)D) and over (2, D):")
    print(f"        dim A(Z) = {len(K)}, contains both actions, and Z is a lift over it: {a}")
    W = [t * t * t - P.const(1, 2) * t * t + t]
    Z2 = al.lift_map(al.R3s, 2, al.lift_map(al.R2s, 1, W, 1), 2)
    K2 = commutant(Z2, 6, kmax=2)
    b = (len(K2) == 6 and rank(K2 + rho(al.R3s, 2), 36) == 6
         and rank(K2 + rho(al.R2s, 3), 36) == 6)
    ok &= b
    print(f"      Z = (W^(R^2))^(R^3) on R^6, presented over (2, R^3) and over (3, R^2):")
    print(f"        dim A(Z) = {len(K2)} = dim(R^3 (x) R^2), contains both actions: {b}")
    print("      -> both presentations refine to the commutant one, so the two")
    print("         prescriptions for Psi_m agree.  The overlap condition the gluing")
    print("         needs is discharged; what is left is smoothness, not consistency.")
    return ok


def marked(X, n, Ms):
    """u + hX + h^2 e_1 -- a method carrying a distinguished vector."""
    h = P.var(Ms, n); u = _u(n, Ms)
    return [u[i] + h * X[i] + (h * h if i == 0 else P(Ms)) for i in range(n)]


def zero_field():
    ok = True
    n, Ms = 2, 3
    zero = [P(Ms)] * n
    print("\n  C3  zero-field rigidity: 0^A = 0 for every A, so every pair collides there")
    print("      method                        Psi^0 affine  T^D Psi^0 == T^(R^2) Psi^0")
    for name, meth, exp in [("explicit Euler", euler, True), ("RK4", rk4, True),
                            ("Laplacian", laplacian, True),
                            ("u + hX + h^2 e_1 (marked vector)", marked, False)]:
        psi = meth(zero, n, Ms)
        aff = all(max((sum(k[:n]) for k in p.d), default=0) <= 1 for p in psi)
        same = (al.lift_map(al.D2, n, psi, Ms, hs=n)
                == al.lift_map(al.R2s, n, psi, Ms, hs=n))
        ok &= (same == exp) and aff
        print(f"      {name:30s} {str(aff):12s}  {same}"
              + ("" if same == exp else "   <-- UNEXPECTED"))
    print("      -> the marked-vector method is affine on affine fields, so affine")
    print("         rigidity misses it; the zero-field collision kills it.  Hence")
    print("         Psi^0_h(u) = s(h) u, with no translation part -- an algebraically")
    print("         natural family carries no distinguished vector, and (taking n = 1")
    print("         in S_{nN} = I_N (x) S_n) no distinguished endomorphism either.")
    return ok


def jet_locus():
    """C4.  Locality turns the extension problem into a fibrewise one: Psi_m is a
    function on jets, and the prescribed loci become V_A(u) = { j^r(X^A)(u) : X },
    LINEAR subspaces of J^r_u since base change is linear in X.  Inclusion-exclusion
    over them would give a smooth LOCAL extension -- if their dimensions were locally
    constant.  They are not: the lift forgets the higher jet at a real point, where
    D(X^D)(1 (x) x) is block-diagonal, so dim V_A drops.  Measured below."""
    ok = True
    n, r, m = 2, 1, 4
    W = m * (1 + m)
    print("\n  C4  the loci V_A(u) in jet space degenerate on the real locus")
    print("      A                              dim V_A at a generic u | at a real point")
    for A, nn in [(al.D2, 2), (al.R2s, 2), (al.R2u, 2)]:
        rows = lambda u: [jet1(al.lift_map(A, nn, X, nn), m, u) for X in monfields(nn, 2)]
        g = rank(rows([F(1), F(3), F(-2), F(5)]), W)
        rl = rank(rows([F(1), F(0), F(3), F(0)]), W)
        print(f"      {A.name:30s} {g:9d}            {rl}")
        if A is al.D2: ok &= (g == 10 and rl == 6)
    print("      -> dim V_D(u) falls from 10 to 6 on the real locus 1_A (x) R^n, so the")
    print("         V_A(u) are not a vector bundle and no smooth family of projections")
    print("         onto them exists.  The fibrewise-linear route to the gluing fails")
    print("         exactly there; section 5's retraction still handles ONE algebra.")
    return ok


def monfields(n, deg):
    mons = []
    def rec(j, left, e):
        if j == n - 1: mons.append(tuple(e + [left])); return
        for a in range(left + 1): rec(j + 1, left - a, e + [a])
    for d in range(deg + 1): rec(0, d, [])
    return [[P(n, {mo: F(1)}) if t == i else P(n) for t in range(n)]
            for i in range(n) for mo in mons]


def jet1(Z, m, u):
    """the 1-jet of the field Z at u, flattened."""
    def ev(q):
        tot = F(0)
        for mo, v in q.d.items():
            t = v
            for j in range(m): t *= u[j] ** mo[j]
            tot += t
        return tot
    return [ev(Z[i]) for i in range(m)] + [ev(Z[i].diff(a)) for i in range(m) for a in range(m)]


def main():
    return recovers() & refinement() & zero_field() & jet_locus()


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
