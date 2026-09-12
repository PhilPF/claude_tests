"""THE PROPERTY ON FINITE-DIMENSIONAL COMMUTATIVE ALGEBRAS, AND WHAT IT NECESSARILY
FORCES.  Sufficient conditions (be a program) do not characterise; these are the
conditions a method CANNOT avoid.

Definition.  Let FCA be the based finite-dimensional commutative associative unital
R-algebras (A, iota), N = dim_R A, iota : R^N -> A part of the data, with the blockwise
realization iota_n : A^n -> R^{nN}.  Fix a data class D (polynomial, or analytic on a
domain).  A method is a family Psi = {Psi_d}, local, smooth in h, consistent.  Then

    Psi is ALGEBRAICALLY NATURAL over D  iff  for all n, all (A, iota) in FCA,
    all X in D-fields on R^n and all small real h,

        Psi_{nN}[ iota_n* X^A ]_h  o  iota_n  =  iota_n  o  ( Psi_n[X]_h )^A .        (AN)

Six necessary conditions, in increasing strength of what they see.

P1  MODULE RIGIDITY.  The right side of (AN) is iota_n of an A-map, so Psi^{X^A}_h must
    be A-differentiable: D Psi(u) lies in the commutant of rho_A(A).  For A = C this is
    literally "holomorphic field in, holomorphic map out" -- generalized Cauchy-Riemann.
    Intrinsic form, since the commutant recovers the algebra: Psi must be equivariant
    for A(Z), the commutant of the field it is handed.  STRICTLY WEAKER than closure:
    the Laplacian method satisfies it over C and still fails (AN) there.
P2  WHICH PROJECTIONS EXIST.  The unit eta : R -> A always exists, so the real locus is
    always there and Psi must restrict to Psi_n on it.  An augmentation pi : A -> R
    exists iff A has a real point, and #Hom(A,R) is the SIGNATURE p - q of the trace
    form.  C has none.  So the (T-base)/(T-fibre) split of section 1 is a feature of
    SPLIT algebras, not of the property.
P3  THE LOCUS SATISFIES DIFFERENTIAL IDENTITIES.  A-linearity of DZ is a linear PDE
    system on the lifted locus.  Over C it is the Cauchy-Riemann system, whose
    consequence is that every component is HARMONIC.  So the real Laplacian annihilates
    the C-lifted locus -- which is exactly why c_C = 0, and why the Laplacian method
    dies over C while living over D.
P4  THE DEFORMATION.  A_lam = R[x]/(x^2 - lam) based by (1,x) is ONE algebraic family
    of based algebras, all of dimension 2, sweeping R^2 (lam > 0), D (lam = 0) and
    C (lam < 0), with c_A = 1 + lam taking EVERY real value.  For a method analytic in
    the jet the closure defect is analytic (here polynomial) in lam, so its vanishing
    set is all of the family or a discrete set.  (T) is the single point lam = 0, and
    the Laplacian method's defect has lam-degree 1 vanishing exactly there: closure over
    dimension 2 is STRICTLY STRONGER than (T), at the same pair of dimensions.
P5  THE SPECTRUM IS ANALYTICALLY CLOSED.  Consequence of P4: no method analytic in the
    jet can have closure spectrum exactly the split (Weil) algebras, because split
    meets the family in [0, infinity), which is neither the whole line nor discrete.
    The Weil/non-Weil distinction is not the closure spectrum of ANY such method.
P6  AND THE CEILING DOES NOT MOVE.  Collisions across the family are still affine only,
    so the enlarged class prunes harder without fattening the lifted locus: the
    extension-theorem junk survives, and affine rigidity remains the ceiling.
"""
from fractions import Fraction as F
from poly import P
import algebra as al
from test_affine import gen_map, euler, heun, laplacian
from test_spectrum import METHODS

GRID = [F(-2), F(-1), F(-1, 2), F(0), F(1, 2), F(1), F(2)]


def _rho(A, n, a):
    N, m = A.N, n * A.N
    L = [[sum(a[x] * A.c[x][b][g] for x in range(N)) for b in range(N)] for g in range(N)]
    M = [[F(0)] * m for _ in range(m)]
    for j in range(n):
        for g in range(N):
            for b in range(N): M[j * N + g][j * N + b] = L[g][b]
    return M


def _A_differentiable(A, n, Fm, Mt):
    m = n * A.N
    J = [[Fm[i].diff(k) for k in range(m)] for i in range(m)]
    for x in range(A.N):
        R = _rho(A, n, [F(1) if b == x else F(0) for b in range(A.N)])
        for i in range(m):
            for k in range(m):
                if sum((J[i][t] * R[t][k] for t in range(m)), P(Mt)) != \
                   sum((R[i][t] * J[t][k] for t in range(m)), P(Mt)): return False
    return True


def coord1(X, n, Ms):
    """u + hX + h^2 X^{(0)} e_0: a marked coordinate.  Violates P1 outright."""
    h = P.var(Ms, n)
    out = [P.var(Ms, i) + h * X[i] for i in range(n)]
    out[0] = out[0] + h * h * X[0]
    return out


def module_rigidity():
    """P1: A-differentiability of the output is necessary, and strictly weaker."""
    ok = True
    n, Ms = 2, 3
    X = gen_map(n, Ms, 3)
    Alist = [al.CC, al.D2, al.Alam(1)]
    print("  P1  module rigidity: is Psi^{X^A} A-differentiable?   (necessary for (AN))")
    print(f"      {'method':12s} " + "  ".join(f"{'A = ' + A.name.split(' ')[0]:24s}" for A in Alist))
    exp = {"euler": [(1, 1), (1, 1), (1, 1)], "heun": [(1, 1), (1, 1), (1, 1)],
           "laplacian": [(1, 0), (1, 1), (1, 0)], "coord1": [(0, 0), (0, 0), (0, 0)]}
    for lbl, meth in [("euler", euler), ("heun", heun), ("laplacian", laplacian),
                      ("coord1", coord1)]:
        cells, got = [], []
        for A in Alist:
            m, Mt = n * A.N, n * A.N + 1
            XA = al.lift_map(A, n, X, Ms, hs=n)
            out = meth(XA, m, Mt)
            p1 = _A_differentiable(A, n, out, Mt)
            cl = (out == al.lift_map(A, n, meth(X, n, Ms), Ms, hs=n))
            got.append((int(p1), int(cl)))
            cells.append(f"P1 {str(p1):5s} (AN) {str(cl):5s}")
        ok &= (got == exp[lbl])
        print(f"      {lbl:12s} " + "  ".join(f"{c:24s}" for c in cells)
              + ("" if got == exp[lbl] else "   <-- UNEXPECTED"))
    print("      -> the Laplacian row is the point: P1 holds over C and (AN) fails, so P1")
    print("         is necessary and NOT sufficient.  coord1 shows P1 is not vacuous.")
    return ok


def projections():
    """P2: eta is universal, pi is not; #Hom(A,R) is the trace-form signature."""
    ok = True
    print("\n  P2  which projections exist: #Hom(A,R) = signature p - q of the trace form")
    for lam, exp in [(F(-4), 0), (F(-1), 0), (F(0), 1), (F(1), 2), (F(4), 2), (F(1, 4), 2)]:
        A = al.Alam(lam)
        got = al.real_points(A)
        ok &= (got == exp)
        print(f"      lam = {str(lam):>5s}   real points {got}   (roots of t^2 = {lam}: {exp})"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    agree = all(al.real_points(A) == al.num_real_points(A)
                for d in al.ALGEBRAS for A in al.ALGEBRAS[d])
    ok &= agree
    print(f"      on the 17 split based algebras p-q = dim A/Nil: {agree}; off them they")
    print("      differ, e.g. C: 0 real points but dim A/Nil = 2.")
    print("      -> no augmentation over C, so (T-base)/(T-fibre) is not part of (AN).")
    print("         The unit eta : R -> A is what is universal, and the real locus with it.")
    return ok


def harmonicity():
    """P3: over C the lifted locus is harmonic, which IS the vanishing of c_A."""
    ok = True
    n, Ms = 2, 2
    X = gen_map(n, Ms, 3)
    print("\n  P3  the lifted locus satisfies differential identities")
    for A, lbl, exp in [(al.CC, "C", True), (al.D2, "D", False), (al.Alam(1), "R^2 (1,x)", False)]:
        m = n * A.N
        XA = al.lift_map(A, n, X, Ms)
        got = all(p.is_zero() for p in
                  [sum((XA[i].diff(k).diff(k) for k in range(m)), P(m)) for i in range(m)])
        ok &= (got == exp)
        print(f"      {lbl:12s} the REAL Laplacian kills X^A: {str(got):6s} (expected {exp})"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    lap = [sum((X[i].diff(j).diff(j) for j in range(n)), P(Ms)) for i in range(n)]
    ok &= not all(p.is_zero() for p in lap)
    print(f"      while (Delta X)^A, the transported one, is nonzero: {not all(p.is_zero() for p in lap)}")
    print("      and both are instances of one EXACT OPERATOR IDENTITY on the locus:")
    print("         Delta_{R^{nN}} (X^A)  =  rho_A(c_A) . (Delta_A X)^A")
    for A in [al.D2, al.CC, al.Alam(F(1)), al.Alam(F(-1, 2)), al.J3, al.W2, al.R2s,
              al.DD, al.CxR]:
        N, m = A.N, n * A.N
        XA = al.lift_map(A, n, X, Ms)
        lhs = [sum((XA[i].diff(k).diff(k) for k in range(m)), P(m)) for i in range(m)]
        lapA = [sum((_dA(A, _dA(A, XA[i], j, m), j, m) for j in range(n)), P(m))
                for i in range(m)]
        c = al.c_A(A)
        rhs = []
        for j in range(n):
            rhs += al.amul(A, [P.const(m, v) for v in c],
                           [lapA[j * N + a] for a in range(N)], m)
        got = (lhs == rhs)
        ok &= got
        print(f"         {A.name:32s} c_A = {str([str(v) for v in c]):24s} {got}"
              + ("" if got else "   <-- UNEXPECTED"))
    print("      -> Cauchy-Riemann forces each component harmonic, so the real Laplacian")
    print("         annihilates the C-locus: that is c_C = 1 + i^2 = 0, read as a PDE.  And")
    print("         section 4's liana law is now a consequence rather than a table -- a method")
    print("         using Delta at dimension n dim A meets rho_A(c_A) where transport wants")
    print("         the identity, so it closes iff c_A = 1_A.")
    return ok


def _dA(A, p, j, Mt):
    """the A-derivative in slot j: the directional derivative along 1_A."""
    out = P(Mt)
    for a in range(A.N):
        if A.one[a]: out = out + p.diff(j * A.N + a) * A.one[a]
    return out


def _defect(A, meth, n=2):
    Ms, m = n + 1, n * A.N
    X = gen_map(n, Ms, 3)
    got = meth(al.lift_map(A, n, X, Ms, hs=n), m, m + 1)
    want = al.lift_map(A, n, meth(X, n, Ms), Ms, hs=n)
    return [got[i] - want[i] for i in range(m)]


def _lam_degree(meth, grid=range(0, 8)):
    seqs = {}
    for lam in grid:
        for i, p in enumerate(_defect(al.Alam(F(lam)), meth)):
            for k, v in p.d.items(): seqs.setdefault((i, k), {})[lam] = v
    best = -1
    for d in seqs.values():
        col = [d.get(l, F(0)) for l in grid]
        deg = len(col) - 1
        while deg >= 0:
            c = col[:]
            for _ in range(deg): c = [c[i + 1] - c[i] for i in range(len(c) - 1)]
            if any(x != 0 for x in c): break
            deg -= 1
        best = max(best, deg)
    return best


def deformation():
    """P4: one family, all of dimension 2, and the defect is polynomial in lam."""
    ok = True
    print("\n  P4  A_lam = R[x]/(x^2 - lam) based by (1,x): one family, dimension 2 throughout")
    print(f"      {'lam':>6s} {'c_A':>6s} {'sig':>11s} {'split':>6s} {'real pts':>9s} {'dim Der':>8s}")
    for lam in GRID:
        A = al.Alam(lam)
        p, q, z = al.trace_signature(A)
        ok &= (al.c_A(A)[0] == 1 + lam) & (al.is_split(A) == (lam >= 0)) \
            & (len(al.derivations(A)) == (1 if lam == 0 else 0))
        print(f"      {str(lam):>6s} {str(al.c_A(A)[0]):>6s} {f'(+{p},-{q},0^{z})':>11s} "
              f"{str(al.is_split(A)):>6s} {p - q:>9d} {len(al.derivations(A)):>8d}")
    print("      -> c_A = 1 + lam takes every real value INSIDE dimension 2, and lam = 0 is")
    print("         the unique member with a derivation: T is the singular point.")
    print(f"\n      {'method':34s} {'lam-degree':>10s}   (AN) holds at lam =")
    exp = {"explicit Euler": (-1, GRID), "Taylor-2  u+hX+h^2/2 DX.X": (-1, GRID),
           "Laplacian  u+hX+h^2 Delta X": (1, [F(0)])}
    for lbl, meth, _ in METHODS:
        deg = _lam_degree(meth)
        good = [l for l in GRID if all(p.is_zero() for p in _defect(al.Alam(l), meth))]
        if lbl in exp:
            e = exp[lbl]; ok &= (deg == e[0]) and (good == e[1])
        print(f"      {lbl:34s} {deg:>10d}   {[str(l) for l in good] or 'nowhere'}"
              + ("" if lbl not in exp or (deg == exp[lbl][0] and good == exp[lbl][1])
                 else "   <-- UNEXPECTED"))
    print("      -> degree -1 means the defect vanishes identically.  The Laplacian defect")
    print("         is a degree-1 polynomial in lam with its only root at lam = 0: (T) alone")
    print("         is STRICTLY WEAKER than closure over dimension 2.  Killing it used to")
    print("         need dimension 3; one family inside dimension 2 now does it.")
    return ok


def spectrum_closed():
    """P5: the spectrum cannot be the split algebras."""
    ok = True
    print("\n  P5  the spectrum, met against the family, is all of it or discrete")
    split = [l for l in GRID if al.is_split(al.Alam(l))]
    ok &= (split == [l for l in GRID if l >= 0])
    print(f"      split members of the family: lam = {[str(l) for l in split]}  (i.e. lam >= 0)")
    print("      For Psi analytic in the jet the defect is analytic in lam, so its zero set")
    print("      is the whole line or locally finite.  [0, infinity) is neither.")
    print("      -> NO method analytic in the jet has closure spectrum exactly the split")
    print("         (Weil) algebras.  The Weil/non-Weil distinction is not the spectrum of")
    print("         any such method -- the sharpest form of 'the property is not about Weil")
    print("         algebras'.  The measured spectra in the family are the whole line (Euler)")
    print("         or the single point lam = 0 (Laplacian), never the half-line.")
    return ok


def ceiling():
    """P6: collisions across the family are still only affine -- the locus is no fatter."""
    ok = True
    n, Ms = 2, 2
    print("\n  P6  collisions across the family: the ceiling of the two-lift mechanism")
    print("      A_lam vs A_mu                 lin  aff  quad   predicted")
    for lam, mu in [(F(0), F(1)), (F(0), F(-1)), (F(-1), F(1)), (F(1), F(2)), (F(-1), F(-2))]:
        A, B = al.Alam(lam), al.Alam(mu)
        res = [al.lift_map(A, n, gen_map(n, Ms, d, c), Ms) ==
               al.lift_map(B, n, gen_map(n, Ms, d, c), Ms)
               for d, c in ((1, False), (1, True), (2, True))]
        e = [True, True, False]                      # same unit (1,0) -> affine, never more
        ok &= (res == e)
        print(f"      lam = {str(lam):>4s}, mu = {str(mu):>4s}            "
              + "  ".join("Y" if r else "n" for r in res) + "    affine"
              + ("" if res == e else "   <-- UNEXPECTED"))
    print("      -> affine and no further, for every pair in the family.  So the union of")
    print("         lifted loci over lam meets itself only in affine fields: the enlarged")
    print("         class PRUNES harder but does not FATTEN the locus.  The extension")
    print("         theorem survives, and affine rigidity is still the ceiling.")
    return ok


def main():
    return (module_rigidity() & projections() & harmonicity() & deformation()
            & spectrum_closed() & ceiling())


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
