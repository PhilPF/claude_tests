"""DECOUPLING: the axiom that connects dimensions ADDITIVELY, and what it is.

Closure connects dimensions MULTIPLICATIVELY -- it relates Psi_n to Psi_{n dim A} and
to nothing else, which is exactly why the extension theorem can choose Psi_q freely
for every odd q.  The missing dimension-uniformity axiom has to connect n and m to
n + m.  The obvious candidate is

    (D)   Psi_{n+m}^{X (+) Y}  =  Psi_n^X (+) Psi_m^Y        for INDEPENDENT X, Y.

and it has an exact description in this project's own terms.  Base change along R^k
sends ONE field to the DIAGONAL X (+) ... (+) X -- verified below -- so closure over
R^k is precisely (D) restricted to equal fields.  (D) is its multi-field completion,
and that completion is something base change structurally cannot express: a base
change acts on one field at a time.

  E1  X^{R^k} realizes as the diagonal product: R^k-closure IS diagonal decoupling
  E2  (D) and closure are logically independent, with witnesses in both directions
  E3  leapfrog with the plain shape 'first half | second half' is closed over the
      WHOLE class -- C and the non-split family included -- and fails (D).  So (D) is
      not implied by closure over FCA, and it localizes how leapfrog escapes: closed
      for each single field, not for independent fields in separate blocks
  E4  what (D) costs the extension theorem: Psi_q for odd q is no longer free, it is
      determined on every product locus, and the free data shrinks to Psi_d on
      IRREDUCIBLE fields -- neither a product nor a lift
  E5  the honest negative: at n = 1 the overlap of the two prescriptions is exactly
      the affine fields, and equating them there yields only LINEAR RIGIDITY, which
      closure already gives.  (D) buys nothing at n = 1.

Conjecture, not proved: algebraically natural over FCA + (D) implies affine
equivariance, hence B-series.  It survives every witness here, and it holds for
contraction-type methods by section 4 (full closure forces N_0 = 1, i.e. trees, and
trees are GL-equivariant).  The attack is whether the extension-theorem construction
can be redone respecting (D) on the irreducible locus.
"""
from fractions import Fraction as F
from poly import P
import algebra as al
from test_affine import gen_map, euler, heun, taylor2, laplacian, ring
from test_jettransport import leapfrog
from test_spectrum import natural, tower

CLASS = [al.D2, al.R2s, al.R2u, al.CC, al.Alam(F(-1, 2)), al.Alam(F(2)), al.J3,
         al.W2, al.DD]


def oplus(X, n, Y, m):
    """X (+) Y on R^{n+m}: X's coordinates first, then Y's, h last."""
    d = n + m; Mt = d + 1
    sx = [P.var(Mt, j) for j in range(n)] + [P.var(Mt, d)]
    sy = [P.var(Mt, n + j) for j in range(m)] + [P.var(Mt, d)]
    return [al.subs(p, sx, Mt) for p in X] + [al.subs(p, sy, Mt) for p in Y]


def decouples(meth, n=2, m=2, deg=3, same=False):
    X = gen_map(n, n + 1, deg)
    Y = X if same and n == m else [p * F(3, 2) for p in gen_map(m, m + 1, deg)]
    return meth(oplus(X, n, Y, m), n + m, n + m + 1) == \
        oplus(meth(X, n, n + 1), n, meth(Y, m, m + 1), m)


def halves(X, n, Ms):
    """leapfrog for 'first half | second half' -- a PLAIN dimension-indexed family,
    since the block realization carries the first half of n to the first half of nN."""
    return leapfrog(X, n, Ms, (list(range(n // 2)), list(range(n // 2, n))))


def diagonal(X, n, Ms):
    """u + hX + h^2 (d_i X^i) X^i, no sum: coordinate-wise, so it decouples."""
    h = P.var(Ms, n)
    return [P.var(Ms, i) + h * X[i] + h * h * X[i].diff(i) * X[i] for i in range(n)]


def mark0(X, n, Ms):
    """u + hX + h^2 (d_0 X^0) X: a marked coordinate in the scalar."""
    h = P.var(Ms, n); s = X[0].diff(0)
    return [P.var(Ms, i) + h * X[i] + h * h * s * X[i] for i in range(n)]


def diagonal_identification():
    """E1: base change along R^k gives the DIAGONAL, so R^k-closure is (D) for X = Y."""
    ok = True
    n, Ms = 2, 3
    X = gen_map(n, Ms, 3)
    print("  E1  base change along R^k produces the diagonal product, not a general one")
    for A, k in [(al.R2s, 2), (al.R3s, 3)]:
        want = X
        for c in range(1, k): want = oplus(want, n * c, X, n)
        got = (al.lift_map(A, n, X, Ms, hs=n, shuffle=True) == want)
        ok &= got
        print(f"      X^({A.name}) in the shuffled realization = X (+) ... (+) X "
              f"({k} copies): {got}" + ("" if got else "   <-- UNEXPECTED"))
    print("      -> so closure over R^k is exactly (D) with the SAME field in every")
    print("         factor.  A base change acts on one field at a time and cannot")
    print("         express independent fields; (D) is that completion.")
    return ok


def independence():
    """E2: neither axiom implies the other."""
    ok = True
    print("\n  E2  (D) against closure -- independent in both directions")
    print(f"      {'method':32s} {'closed/D':>9s} {'diag (D)':>9s} {'full (D)':>9s}")
    rows = [("explicit Euler", euler, (2, 2), (True, True, True)),
            ("Heun", heun, (2, 2), (True, True, True)),
            ("Taylor-2", taylor2, (2, 2), (True, True, True)),
            ("Laplacian u+hX+h^2 Delta X", laplacian, (2, 2), (True, True, True)),
            ("diagonal (d_i X^i) X^i", diagonal, (2, 2), (False, True, True)),
            ("mark0  (d_0 X^0) X", mark0, (2, 2), (False, False, False)),
            ("ring   (X^0)^2 X", ring, (2, 2), (False, False, False)),
            ("tower(D), dims 1+1", tower(al.D2), (1, 1), (True, False, False)),
            ("leapfrog, halves", halves, (2, 2), (True, False, False))]
    for lbl, meth, (n, m), exp in rows:
        got = (natural(al.D2, meth), decouples(meth, n, n, same=True), decouples(meth, n, m))
        ok &= (got == exp)
        print(f"      {lbl:32s} {str(got[0]):>9s} {str(got[1]):>9s} {str(got[2]):>9s}"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    print("      -> the diagonal method decouples and is not closed; leapfrog and the")
    print("         extension-theorem witness tower(D) are closed and do not decouple.")
    print("         And diagonal (D) is strictly weaker than full (D): the tower fails")
    print("         both, but the gap is real -- see E1.")
    return ok


def leapfrog_across():
    """E3: leapfrog's plain shape is closed over the whole class and still fails (D)."""
    ok = True
    print("\n  E3  leapfrog, shape 'first half | second half', over the whole class")
    for A in CLASS:
        c = natural(A, halves)
        ok &= c
        print(f"      {A.name:30s} closed: {c}" + ("" if c else "   <-- UNEXPECTED"))
    d = decouples(halves, 2, 2)
    ok &= not d
    print(f"      and (D): {d}   (at dimensions 2 + 2: the shape of R^4 is {{0,1}}|{{2,3}},")
    print("      while the two factors contribute {0}|{1} and {2}|{3}, i.e. {0,2}|{1,3})")
    print("      -> closure over ALL of FCA does not imply (D).  Stormer-Verlet is closed")
    print("         for each single field and not for independent fields in two blocks.")
    return ok


def extension_cost():
    """E4: (D) connects dimensions additively, which is what the tower exploited."""
    ok = True
    print("\n  E4  what (D) costs the extension theorem")
    t = tower(al.D2)
    for n, m in [(1, 1), (2, 2), (1, 2)]:
        got = decouples(t, n, m)
        ok &= not got
        print(f"      tower(D) decouples at dims {n} + {m}: {got}"
              + ("" if not got else "   <-- UNEXPECTED"))
    print("      -> closure relates Psi_n only to Psi_{n dim A}, so Psi_q was free for")
    print("         every odd q.  (D) relates Psi_{n+m} to Psi_n and Psi_m, so Psi_q is")
    print("         now determined on every product locus by lower dimensions.  The free")
    print("         data shrinks to Psi_d on IRREDUCIBLE fields -- neither a product nor")
    print("         a lift.  That is a smaller residue, not an empty one.")
    return ok


def honest_negative():
    """E5: the n = 1 overlap gives only linear rigidity, which closure already has."""
    ok = True
    n, Ms = 1, 2
    print("\n  E5  where the two prescriptions overlap at n = 1, and why it buys nothing")
    for deg, lbl, exp in [(1, "affine  ax+b", True), (2, "quadratic", False),
                          (3, "cubic", False)]:
        X = gen_map(n, Ms, deg)
        XD = al.lift_map(al.D2, n, X, Ms)
        prod = XD[0].diff(1).is_zero() and XD[1].diff(0).is_zero()
        ok &= (prod == exp)
        print(f"      X^D is a product field for a {lbl:14s} X: {str(prod):5s} (expected {exp})"
              + ("" if prod == exp else "   <-- UNEXPECTED"))
    print("      -> the overlap is exactly the AFFINE fields.  There closure gives")
    print("         Psi_2^{X^D} = T Psi_1^X and (D) gives Psi_1^{ax+b} (+) Psi_1^{ax};")
    print("         equating the fibres forces Psi_1^{ax}(v) linear with no constant.")
    print("         That is LINEAR RIGIDITY, which closure alone already gives, so (D)")
    print("         adds nothing at n = 1.  Recorded so the dead end is not re-walked.")
    return ok


def main():
    return (diagonal_identification() & independence() & leapfrog_across()
            & extension_cost() & honest_negative())


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
