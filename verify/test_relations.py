"""How the closure spectrum behaves under the relations between algebras.

The Weil algebras form a category, and one would like to check naturality on a
generating family and deduce the rest.  That fails, and the failures are sharp.

D1.  NOT CLOSED UNDER QUOTIENTS.  tower(D (x) D) is T^{D(x)D}-natural but not
     T^D-natural, and D is a quotient of D (x) D (send x, y -> eps).  Since EVERY
     Weil algebra is a quotient of a jet algebra R[x_1..x_k]/m^{r+1}, this kills
     the obvious reduction: one cannot verify naturality on the jet algebras and
     deduce it for the rest.

D2.  NOT CLOSED UNDER PRODUCTS.  tower(D) is T^D-natural but not T^{DxD}-natural.

     With the earlier failure under subalgebras, the spectrum is closed under (x)
     and under nothing else.  The reason is structural: (x) is the only relation
     among Weil algebras that is a relation between the FUNCTORS, T^{A(x)B} =
     T^A o T^B, so one condition can be applied twice.  Sub, quotient and product
     give natural transformations BETWEEN functors, and naturality is a condition
     on each functor separately -- Psi_m and Psi_{m'} in different dimensions are
     a priori unrelated, so nothing transfers.

D3.  WHAT THE MORPHISMS DO GIVE.  For phi: A -> B the map phi_n = phi (x) id is
     real-linear, intertwines X^A with X^B, and satisfies T^B F . phi_n = phi_n .
     T^A F.  So an algebraically natural Psi obeys

         phi_n . Psi^{X^A} = Psi^{X^B} . phi_n     on lifted fields,

     which for A = B is Aut(A)-equivariance of Psi_m on the A-lifted locus.  This
     is not vacuous: Aut(R[x_1..x_k]/m^2) = GL(k), so T^{J^1_k}-naturality forces
     honest GL(k)-equivariance -- the first genuine GL in this programme.  It acts
     in the ALGEBRA directions, not the R^n directions, which is exactly the gap
     section 11 records against "algebraically natural => affine equivariant".
"""
from fractions import Fraction as F
from itertools import combinations
from poly import P
import algebra as al
from test_affine import gen_map, _u, euler
from test_spectrum import tower, power, natural


def is_hom(pi, A, B):
    """pi: a list of B-coordinate vectors, one per basis element of A."""
    def bmul(x, y):
        o = [F(0)] * B.N
        for a in range(B.N):
            for b in range(B.N):
                if x[a] and y[b]:
                    for g in range(B.N): o[g] += x[a] * y[b] * B.c[a][b][g]
        return o
    if [sum(A.one[a] * pi[a][g] for a in range(A.N)) for g in range(B.N)] != B.one:
        return False
    for a in range(A.N):
        for b in range(A.N):
            lhs = [sum(A.c[a][b][c] * pi[c][g] for c in range(A.N)) for g in range(B.N)]
            if lhs != bmul(pi[a], pi[b]): return False
    return True


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


def closure_ops():
    ok = True
    DD = al.tensor(al.D2, al.D2)
    pi = [[F(1), F(0)], [F(0), F(1)], [F(0), F(1)], [F(0), F(0)]]   # x, y -> eps
    h, sur = is_hom(pi, DD, al.D2), rank(pi, 2) == 2
    ok &= h and sur
    print("  D1  quotients.  pi: D(x)D -> D sending x, y -> eps is a surjective")
    print(f"      algebra homomorphism: hom {h}, surjective {sur}")
    a, b = natural(DD, tower(DD), 3), natural(al.D2, tower(DD), 3)
    ok &= (a and not b)
    print(f"      tower(D(x)D) is T^(D(x)D)-natural {a}, T^D-natural {b}")
    print("      -> the spectrum contains the parent and not the quotient.  Since every")
    print("         Weil algebra is a quotient of a jet algebra, naturality cannot be")
    print("         reduced to the jet algebras.")
    c, d = natural(al.D2, tower(al.D2), 3), natural(al.DxD, tower(al.D2), 3)
    ok &= (c and not d)
    print(f"\n  D2  products.  tower(D) is T^D-natural {c}, T^(DxD)-natural {d}")
    print("      -> with the earlier failure under subalgebras: the spectrum is closed")
    print("         under (x) and under nothing else.  (x) is the only relation that is")
    print("         one between the functors, T^{A(x)B} = T^A o T^B.")
    return ok


def aut_equivariance():
    ok = True
    W2, n, Ms = al.W2, 2, 3
    m = n * W2.N
    Lm = lambda g: [[F(1), F(0), F(0)], [F(0), g[0][0], g[0][1]], [F(0), g[1][0], g[1][1]]]

    def post(vecs, g, M):                       # apply phi_n to values: v -> L v
        L = Lm(g); out = []
        for j in range(n):
            blk = vecs[j * 3:(j + 1) * 3]
            out += [sum((blk[b] * L[a][b] for b in range(3)), P(M)) for a in range(3)]
        return out

    def pre(polys, g, M):                       # precompose with phi_n
        L = Lm(g)
        im = [sum((P.var(M, j * 3 + b) * L[a][b] for b in range(3)), P(M))
              for j in range(n) for a in range(3)]
        if M > m: im.append(P.var(M, m))
        return [al.subs(p, im, M) for p in polys]

    def marked_block(Z, mm, MM):                # singles out the x-block coordinate
        hh = P.var(MM, mm); u = _u(mm, MM); s = Z[1] * Z[1]
        return [u[i] + hh * Z[i] + hh * hh * s * Z[i] for i in range(mm)]

    g = [[F(2), F(-1)], [F(3), F(1, 2)]]
    XA = al.lift_map(W2, n, gen_map(n, Ms, 3), Ms, hs=n)
    print("\n  D3  Aut(A)-equivariance on the lifted locus, for A = R[x,y]/m^2 (Aut = GL(2))")
    for label, obj, exp in [("the lifted field X^A", XA, True),
                            ("Euler on the locus", euler(XA, m, m + 1), True),
                            ("a method singling out the x-block", marked_block(XA, m, m + 1), False)]:
        got = (post(obj, g, m + 1) == pre(obj, g, m + 1))
        ok &= (got == exp)
        print(f"      commutes with phi_n:  {label:34s} {got}"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    print("      -> T^{J^1_k}-naturality forces honest GL(k)-equivariance of Psi on the")
    print("         lifted locus, in the ALGEBRA directions, not the R^n directions.")
    return ok


def wedge_even(k):
    """The even part of Lambda R^k as an ordinary based algebra: basis the even wedge
    monomials, which commute because even elements are central in a supercommutative
    algebra.  An ordinary manifold sees only this part -- the A-points of R^n are
    (A_even)^n, there being no odd parameters to pair with odd basis vectors."""
    basis = []
    for d in range(0, k + 1, 2):
        basis += [frozenset(c) for c in combinations(range(k), d)]
    basis.sort(key=lambda t: (len(t), sorted(t)))
    idx = {b: i for i, b in enumerate(basis)}
    N = len(basis)

    def sign(I, J):
        arr, s = sorted(I) + sorted(J), 1
        for _ in range(len(arr)):
            for b in range(len(arr) - 1):
                if arr[b] > arr[b + 1]:
                    arr[b], arr[b + 1] = arr[b + 1], arr[b]; s = -s
        return s

    table = {}
    for i, I in enumerate(basis):
        for j, J in enumerate(basis):
            if i > j: continue
            v = [F(0)] * N
            if not (I & J): v[idx[I | J]] = F(sign(I, J))
            table[(i, j)] = v
    return al.Alg(f"even part of Lambda R^{k}", N, table, [F(1)] + [F(0)] * (N - 1))


def graded():
    """D4.  The OTHER notion named after Weil -- Cartan's W(g) = Lambda g^* (x) S g^*,
    the model of EG behind Chern-Weil theory -- lives in the graded world, where
    R[theta] with theta ODD gives Pi T = T[1] and C^oo(T[1]M) = Omega(M).  It is not a
    stronger test for us: on an ordinary manifold a super Weil algebra acts through its
    EVEN part, and that is an ordinary Weil algebra already of the kind in our family."""
    ok = True
    print("\n  D4  graded Weil algebras act through their even parts")
    print("      k   dim   unital/assoc/comm   m^2 = 0   identification")
    ids = {2: ("D = R[e]/(e^2)", al.D2), 3: ("J^1_3 = R[x,y,z]/m^2", al.W3)}
    for k in (2, 3, 4, 5):
        A = wedge_even(k)
        m2 = all(all(A.c[a][b][g] == 0 for g in range(A.N))
                 for a in range(1, A.N) for b in range(1, A.N))
        note = ""
        if k in ids:
            name, ref = ids[k]
            same = (A.c == ref.c and A.one == ref.one)
            ok &= same
            note = f"IS {name}: {same}"
        else:
            note = "an ordinary Weil algebra of dim " + str(A.N)
        ok &= A.check()
        print(f"      {k}   {A.N:3d}   {str(A.check()):17s}   {str(m2):7s}   {note}")
    print("      -> even(Lambda R^3) has exactly the structure constants of W_3, which is")
    print("         already in the family (it is the algebra whose Aut = GL(3) generated")
    print("         gl(m-1) in test_aut).  So the graded notion adds no test algebra.")
    print("      -> it would add strength only by enlarging the category the METHOD acts")
    print("         on -- methods on supermanifolds, whose input field has odd components.")
    return ok


def main():
    return closure_ops() & aut_equivariance() & graded()


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
