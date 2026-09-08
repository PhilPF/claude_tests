"""The closure spectrum: which Weil algebras a method closes under.

For a family Psi put  W(Psi) = { based algebras (A, iota) : Psi is T^A-natural }.
This is the object a Galois-correspondence-style account would classify.  What is
established here:

B1.  W(Psi) is a MONOID under (x).  Proof in one line -- Psi^{X^{A(x)B}} =
     Psi^{(X^B)^A} = T^A Psi^{X^B} = T^A T^B Psi^X = T^{A(x)B} Psi^X -- and
     structurally because all three invariants of section 4 are multiplicative:
     dim(A(x)B) = dim A dim B, c_{A(x)B} = c_A c_B, and L_{a(x)b} = L_a (x) L_b.
     The invariants depend on (A, g) alone, g the inner product declaring the basis
     orthonormal; the SAME algebra in two bases can have different spectra.

B2.  For a contraction method the spectrum is cut out by those invariants, and the
     law is coarser than section 4's table suggests: a loop OR A STOLON survives
     only for A = R, while a liana survives iff c_A = 1.  Balancedness is necessary
     for the real form to factor through the A-valued one but is NOT sufficient for
     a stolon to survive: R^k in its standard basis is balanced and has c = 1, yet
     kills stolons.  The operative reason unifies the loop and stolon rows -- both
     deliver a REAL scalar where base change demands an A-scalar, which is scalar
     rigidity.  The defect is exactly (r.1_A - q).X^A, verified below.

B3.  W(Psi) is NOT closed under subalgebras, even with the induced metric.

B4.  The realized spectra form a LATTICE, NOT A CHAIN.  The tower construction of
     section 5 -- take the scalar RING-valued -- realizes the principal monoid
     <A>_(x) for each based A, and <D>_(x), <R^2>_(x) are incomparable.

     A CORRECTION falls out of this.  The write-up called u + hX + h^2 (X^1)^2 X
     "the (T)-natural non-equivariant method"; with the scalar left REAL-valued it
     is not (T)-natural at all (measured: the defect is the fibre term
     2 X^1 (DX^1 v) X).  Section 5's own prescription -- make the scalar ring-valued
     -- is what is needed, and the tower method below is the genuine witness.
"""
from fractions import Fraction as F
from poly import P
import algebra as al
from algebra import ALGEBRAS
from test_affine import gen_map, _u, euler, taylor2, laplacian


def stolon(X, n, Ms):
    """u + hX + h^2 <X,X> X -- a stolon (two upper indices)."""
    h = P.var(Ms, n); u = _u(n, Ms)
    q = sum((X[i] * X[i] for i in range(n)), P(Ms))
    return [u[i] + h * X[i] + h * h * q * X[i] for i in range(n)]


def aroma(X, n, Ms):
    """u + hX + h^2 div(X) X -- a loop."""
    h = P.var(Ms, n); u = _u(n, Ms)
    d = sum((X[i].diff(i) for i in range(n)), P(Ms))
    return [u[i] + h * X[i] + h * h * d * X[i] for i in range(n)]


def liana_stolon(X, n, Ms):
    """u + hX + h^2 |DX|_F^2 X -- one liana feeding one stolon, no loop."""
    h = P.var(Ms, n); u = _u(n, Ms)
    q = sum((X[i].diff(j) * X[i].diff(j) for i in range(n) for j in range(n)), P(Ms))
    return [u[i] + h * X[i] + h * h * q * X[i] for i in range(n)]


def ring_real(X, n, Ms):
    """u + hX + h^2 (X^1)^2 X, scalar left REAL-valued -- NOT (T)-natural."""
    h = P.var(Ms, n); u = _u(n, Ms); s = X[0] * X[0]
    return [u[i] + h * X[i] + h * h * s * X[i] for i in range(n)]


def power(A, k):
    B = al.R1
    for _ in range(k): B = al.tensor(A, B, f"{A.name}^(x){k}")
    return B


def tower(base):
    """The same method with the scalar taken RING-valued: on R^m = A^q with
    m = 2^k q, q odd, and A = base^(x)k, the square and the product happen in A.
    Every dimension is hit exactly once, so this is a family indexed by plain
    dimension, and it is (T)-natural for base = D."""
    def meth(X, m, Ms):
        k, q = 0, m
        while q % 2 == 0: q //= 2; k += 1
        A = power(base, k); N = A.N
        h = P.var(Ms, m); u = _u(m, Ms)
        sig = al.amul(A, X[:N], X[:N], Ms)
        out = []
        for j in range(q):
            blk = al.amul(A, sig, X[j * N:(j + 1) * N], Ms)
            out += [u[j * N + a] + h * X[j * N + a] + h * h * blk[a] for a in range(N)]
        return out
    return meth


def natural(A, meth, n=2):
    """Is meth T^A-natural?  Exact, on a cubic test field with symbolic h."""
    Ms, N = n + 1, A.N
    X = gen_map(n, Ms, 3)
    psi = meth(X, n, Ms)
    return meth(al.lift_map(A, n, X, Ms, hs=n), n * N, n * N + 1) == \
        al.lift_map(A, n, psi, Ms, hs=n)


LOOP, LIANA, STOLON = 1, 2, 4
METHODS = [("explicit Euler", euler, 0), ("Taylor-2  u+hX+h^2/2 DX.X", taylor2, 0),
           ("Laplacian  u+hX+h^2 Delta X", laplacian, LIANA),
           ("stolon  u+hX+h^2 <X,X> X", stolon, STOLON),
           ("aroma  u+hX+h^2 div(X) X", aroma, LOOP),
           ("liana+stolon  u+hX+h^2 |DX|_F^2 X", liana_stolon, LIANA | STOLON)]


def predict(A, deco):
    """A loop or a stolon needs dim A = 1; a liana needs c_A = 1."""
    N, c, c1, bal = al.spectrum_invariants(A)
    return not ((deco & (LOOP | STOLON) and N != 1) or (deco & LIANA and not c1))


def spectrum():
    ok = True
    print("  B2  measured closure spectrum (Y = T^A-natural), against the invariants")
    for i, (name, _, _) in enumerate(METHODS): print(f"        [{i}] {name}")
    print("      based algebra                    dim c=1 bal | " +
          " ".join(str(i) for i in range(len(METHODS))))
    for A in [A for _, lst in sorted(ALGEBRAS.items()) for A in lst]:
        N, c, c1, bal = al.spectrum_invariants(A)
        row, good = [], True
        for _, meth, deco in METHODS:
            m = natural(A, meth)
            row.append("Y" if m else "n")
            good &= (m == predict(A, deco))
        ok &= good
        print(f"      {A.name:32s} {N:3d} {'Y' if c1 else 'n':3s} {'Y' if bal else 'n':3s} | "
              + " ".join(row) + ("" if good else "   <-- UNEXPECTED"))
    print("      -> loops AND stolons survive only for A = R; lianas iff c_A = 1.")
    print("      -> R^k standard is balanced with c = 1 and still kills stolons, so")
    print("         balancedness is necessary for the real form to factor through the")
    print("         A-valued one, not sufficient for the stolon to survive.")
    print("      -> R^2 standard vs unit-first: one algebra, two bases, two spectra.")

    # the mechanism, as an exact identity: the stolon defect is (r.1_A - q).X^A
    n, Ms, A = 2, 3, al.R2s
    X = gen_map(n, Ms, 3)
    XA = al.lift_map(A, n, X, Ms, hs=n)
    Mt, N = n * A.N + 1, A.N
    defect = [a - b for a, b in zip(stolon(XA, n * N, Mt),
                                    al.lift_map(A, n, stolon(X, n, Ms), Ms, hs=n))]
    q = sum((X[i] * X[i] for i in range(n)), P(Ms))                      # <X,X> on R^n
    qA = al.basechange(A, n, [q], al.generic_point(A, n, Mt), Ms, Mt, n, n * N)[0]
    r = sum((p * p for p in XA), P(Mt))                                  # the real form
    h2 = P.var(Mt, n * N) * P.var(Mt, n * N)
    pred = []
    for j in range(n):
        blk = al.amul(A, [r * A.one[a] - qA[a] for a in range(N)], XA[j * N:(j + 1) * N], Mt)
        pred += [h2 * blk[a] for a in range(N)]
    good = (defect == pred)
    ok &= good
    print(f"      mechanism: stolon defect == h^2 (r.1_A - q).X^A over R^2: {good}"
          + ("" if good else "  <-- UNEXPECTED"))
    return ok


def monoidal():
    ok = True
    print("\n  B1  the spectrum is a monoid under (x), because the invariants are")
    print("      A                     B                     dim  c mult  bal mult  T^AxB = T^A T^B")
    for A, B in [(al.D2, al.D2), (al.D2, al.R2s), (al.R2s, al.R2s), (al.D2, al.W2),
                 (al.R2u, al.D2)]:
        T = al.tensor(A, B)
        cA, cB, cT = al.c_A(A), al.c_A(B), al.c_A(T)
        prod = [F(0)] * T.N
        for a in range(A.N):
            for b in range(B.N):
                if cA[a] and cB[b]: prod[b * A.N + a] += cA[a] * cB[b]
        cm = (prod == cT)
        bm = (al.balanced(T) == (al.balanced(A) and al.balanced(B)))
        dm = (T.N == A.N * B.N)
        n, Ms = 2, 2
        X = gen_map(n, Ms, 3)
        fun = (al.lift_map(T, n, X, Ms) ==
               al.lift_map(A, n * B.N, al.lift_map(B, n, X, Ms), n * B.N))
        ok &= cm and bm and dm and fun
        print(f"      {A.name:21s} {B.name:21s} {str(dm):5s} {str(cm):7s} {str(bm):9s} {fun}"
              + ("" if (cm and bm and dm and fun) else "   <-- UNEXPECTED"))
    a, b = natural(al.D2, laplacian), natural(al.tensor(al.D2, al.D2), laplacian)
    ok &= (a and b)
    print(f"      Laplacian closes under D: {a}; hence under D (x) D: {b} (measured)")
    return ok


def subalgebras():
    ok = True
    print("\n  B3  the spectrum is NOT closed under subalgebras")
    S = al.Alg("R[s]/(s^3)", 3, {(1, 1): [0, 0, 2], (1, 2): [0] * 3, (2, 2): [0] * 3}, [1, 0, 0])
    g = [[F(1), F(0), F(0)], [F(0), F(2), F(0)], [F(0), F(0), F(1)]]
    amb, sub, c_ind = natural(al.tensor(al.D2, al.D2), laplacian), natural(al.J3, laplacian), al.c_A(S, g)
    ok &= (amb and not sub and c_ind != S.one)
    print(f"      R[s]/(s^3) is the S_2-invariant subalgebra of D (x) D.  With the metric")
    print(f"      induced from the ambient orthonormal basis, diag(1,2,1), its c is")
    print(f"      1 + s^2/2 = 1 + xy, not 1.  Laplacian: T^(D(x)D)-natural {amb}, T^(R[e]/(e^3))-natural {sub}")
    print("      -> the obstruction is not inherited by subalgebras, because c is not.")
    return ok


def realized():
    ok = True
    print("\n  B4  which spectra are realized -- a lattice, not a chain")
    bad = natural(al.D2, ring_real)
    ok &= not bad
    print(f"      u+hX+h^2 (X^1)^2 X with the scalar REAL-valued is (T)-natural: {bad}"
          + ("   <-- UNEXPECTED" if bad else "   (so the write-up's label was wrong)"))
    tD, tR = tower(al.D2), tower(al.R2s)
    print("      the same method with the scalar RING-valued, base dimension 3 (odd):")
    print("      A                              tower(D)  tower(R^2)")
    exp = {"D = R[e]/(e^2)": (True, False), "R^2 (standard basis)": (False, True),
           "R[x,y]/m^2": (False, False)}
    for A in [al.D2, al.R2s, al.tensor(al.D2, al.D2), al.tensor(al.R2s, al.R2s), al.W2]:
        d, r = natural(A, tD, 3), natural(A, tR, 3)
        if A.name in exp: ok &= ((d, r) == exp[A.name])
        print(f"      {A.name:30s} {str(d):8s}  {r}")
    print("      -> tower(D) realizes <D>_(x) and tower(R^2) realizes <R^2>_(x); each")
    print("         contains an algebra the other excludes, so the realized spectra are")
    print("         incomparable.  With {R}, {c_A = 1} (Laplacian) and everything (RK),")
    print("         the correspondence method-class -> (x)-closed class is not a chain.")
    return ok


def main():
    return spectrum() & monoidal() & subalgebras() & realized()


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
