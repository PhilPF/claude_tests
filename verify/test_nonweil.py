"""IS THE PROPERTY ABOUT WEIL ALGEBRAS AT ALL?  No -- it is about finite-dimensional
commutative R-algebras, and the Weil restriction comes from the REGULARITY of the
data, not from the property.

Kolar-Michor-Slovak classify the product-preserving functors on Mf as the Weil
functors, and this project took that as the legitimate class.  But Mf is the SMOOTH
category, and that is where the restriction comes from.  A finite-dimensional
commutative R-algebra is a product of local ones, and the only finite field
extensions of R are R and C, so each residue field is R or C:

    A is a product of WEIL algebras   <=>   every residue field is R
                                      <=>   the trace form on A/Nil(A) is positive
                                            definite                          [W1]
                                      <=>   SMOOTH base change is defined,

the last because f^A(u) = sum_a d^a f(pi u)(u - pi u)^a / a! terminates only when
u - pi(u) is nilpotent, and a C residue field leaves a direction in which no real
expansion point makes it so.                                                   [W4]

Drop smoothness and the class grows.  For POLYNOMIAL (or, locally, analytic) data
base change is just substitution, it works over any commutative R-algebra, and the
closure property is perfectly well posed.  Verified: Euler, Heun, the midpoint rule
and RK4 are T^C-natural, C is not a Weil algebra, and section 4's trichotomy
predicts every row of the table over the non-Weil algebras too.                [W2]

C is not a formality: it is a STRICTLY NEW test.  Over a product of Weil algebras
c_A is a unit in every basis (the image of sum_a e_a^2 in A/Nil = R^k has j-th
component sum_a (e_a)_j^2 > 0, the e_a spanning), so the liana factor can never
vanish.  Over C in the basis (1, i) it is 0.  So the Laplacian method DIES over C
while it SURVIVES over D -- a separation no Weil algebra can produce.           [W3]

And the second half of the question -- relating the lift to other vector fields.
The fields pi_A-related to X are far more than the lifts: every derivation of A
gives a vertical linear field delta_A on A^n, and X^A + delta_A is pi_A-related to X
without being a lift.  Verified: [delta_A, X^A] = 0 always (X^A is Aut(A)-invariant,
being natural in A), so the exact flows split as
phi^{X^A+delta}_t = phi^{X^A}_t o T^{exp(t delta)}.  A method keeps the BASE half of
that but does NOT split the commuting pair, so the Der(A) part must be exponentiated
exactly -- which is what Aut(A)-equivariance supplies, and it is why nothing new has
to be assumed.                                                                 [W6]

  W1  the boundary is computable: the trace-form signature
  W2  closure holds over non-Weil algebras, and the trichotomy predicts it
  W3  c_A is a unit over every split algebra; over C it vanishes
  W4  where Weil IS needed: smooth base change, and nowhere else
  W5  the collision trichotomy extends across the boundary
  W6  pi_A-related fields beyond the lifts: the Der(A) twists
"""
from fractions import Fraction as F
from poly import P
import algebra as al
from test_affine import gen_map, euler, heun, midpoint, rk4, laplacian
from test_spectrum import natural, METHODS, predict

CD = al.tensor(al.CC, al.D2, "C (x) D")
CC2 = al.tensor(al.CC, al.CC, "C (x) C")
NEW = [al.CC, al.CxR, CD, CC2]


def boundary():
    """W1: split (= a product of Weil algebras) iff the trace form has no negative
    square.  A real local factor contributes +, a complex one the pair diag(+,-)."""
    ok = True
    print("  W1  the trace form B(x,y) = tr(L_xy): signature (+p, -q, 0^z), z = dim Nil")
    print("      A                                sig            split  is a Weil product")
    rows = [(A, True) for d in sorted(al.ALGEBRAS) for A in al.ALGEBRAS[d]]
    rows += [(al.Rx2, True)] + [(A, False) for A in NEW]
    for A, exp in rows:
        p, q, z = al.trace_signature(A)
        got = al.is_split(A)
        ok &= (got == exp) & A.check()
        print(f"      {A.name:32s} (+{p},-{q},0^{z})     {str(got):6s} {exp}"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    print("      -> R[x]/(x^2-1) is the control: it LOOKS exotic and is just R^2, and the")
    print("         test says so.  C, C x R, C(x)D, C(x)C are genuinely outside.")
    return ok


def closure_outside():
    """W2: the closure property over algebras that are not Weil products."""
    ok = True
    print("\n  W2  closure over non-Weil algebras (cubic field; rk4 on a quadratic one,")
    print("      where degree 3^4 in 2n dim A slots is otherwise the cost ceiling)")

    def nat2(A, meth, n=2, deg=2):
        Ms, N = n + 1, A.N
        X = gen_map(n, Ms, deg)
        return meth(al.lift_map(A, n, X, Ms, hs=n), n * N, n * N + 1) == \
            al.lift_map(A, n, meth(X, n, Ms), Ms, hs=n)

    hdr = ["euler", "heun", "midpt"] + [l.split()[0] for l, _, _ in METHODS[2:]]
    print(f"      {'algebra':18s} " + " ".join(f"{h:9s}" for h in hdr) + " rk4")
    for A in [al.D2] + NEW[:3]:
        ms = [euler, heun, midpoint] + [m for _, m, _ in METHODS[2:]]
        got = [natural(A, m) for m in ms]
        exp = [True, True, True] + [predict(A, d) for _, _, d in METHODS[2:]]
        r4 = nat2(A, rk4)
        ok &= (got == exp) & r4
        print(f"      {A.name:18s} " + " ".join(f"{str(v):9s}" for v in got) + f" {r4}"
              + ("" if got == exp and r4 else "   <-- UNEXPECTED"))
    print("      -> every row matches the section-4 prediction (a loop or stolon needs")
    print("         dim A = 1, a liana needs c_A = 1_A).  Nothing in the trichotomy, and")
    print("         nothing in the machinery, ever used nilpotency.")
    return ok


def liana_factor():
    """W3: over a split algebra c_A is a unit in EVERY basis; over C it is 0."""
    ok = True
    print("\n  W3  the liana factor c_A = sum_a e_a^2, and what only C can do")
    for A in [A for d in sorted(al.ALGEBRAS) for A in al.ALGEBRAS[d]]:
        ok &= al.is_unit(A, al.c_A(A))
    print(f"      c_A is a unit on all 17 based (split) algebras: {ok}")
    S = {2: [[F(1), F(1)], [F(1), F(-1)]], 3: [[F(1), F(0), F(1)], [F(0), F(1), F(1)],
                                               [F(1), F(1), F(0)]]}
    for A in [al.D2, al.R2s, al.R2u, al.J3, al.W2, al.DxR]:
        B = al.rebase(A, S[A.N])
        u = al.is_unit(B, al.c_A(B))
        ok &= u & B.check() & al.is_split(B)
        print(f"      {A.name:26s} rebased: split {al.is_split(B)}, c_A a unit: {u}")
    for A, lbl in [(al.CC, "C, basis (1,i) "), (al.rebase(al.CC, [[F(1), F(0)], [F(0), F(2)]]),
                                                "C, basis (1,2i)")]:
        c = al.c_A(A)
        print(f"      {lbl:26s} c_A = {[str(v) for v in c]}   a unit: {al.is_unit(A, c)}")
    ok &= (al.c_A(al.CC) == [F(0)] * 2) and not al.is_unit(al.CC, al.c_A(al.CC))
    print("      -> c_A = 0 needs a complex residue field; in the split world every basis")
    print("         gives a unit, since the components of sum_a e_a^2 in A/Nil = R^k are")
    print("         sums of squares of spanning vectors.  A basis change can still make a")
    print("         NON-split c_A a unit (-3 above), so the implication runs one way only.")
    return ok


def regularity():
    """W4: where Weil is genuinely needed -- and it is only smoothness."""
    ok = True
    print("\n  W4  what the Weil condition actually buys: base change of SMOOTH data")
    for A, lbl, exp in [(al.D2, "D", True), (al.J3, "R[e]/(e^3)", True),
                        (al.CC, "C", False), (al.CxR, "C x R", False)]:
        m = [[F(1) if b == a else F(0) for b in range(A.N)] for a in range(1, A.N)]
        nil = [al.is_nilpotent(A, v) for v in m]
        ok &= (all(nil) == exp)
        print(f"      {lbl:12s} e_1..e_{A.N-1} nilpotent: {nil}"
              + ("" if all(nil) == exp else "   <-- UNEXPECTED"))
    bad = [al.is_nilpotent(al.CC, [-c, F(1)]) for c in (F(0), F(1), F(2), F(-1))]
    ok &= not any(bad)
    print(f"      in C at u = i, u - c*1 nilpotent for c = 0,1,2,-1: {bad}")
    print("      (and never: (i-c)^2 = -1-c^2-2ci = 0 needs c^2 = -1.)  So (*) does not")
    print("      terminate and f^C is undefined for general smooth f.")
    n, Ms, deg = 2, 2, 3
    f = gen_map(n, Ms, deg); g = [f[1], f[0] + f[1]]
    gof = [al.subs(p, f, Ms) for p in g]
    for A in NEW:
        Mt = n * A.N
        L, Lg = al.lift_map(A, n, f, Ms), al.lift_map(A, n, g, Ms)
        fun = ([al.subs(p, L, Mt) for p in Lg] == al.lift_map(A, n, gof, Ms))
        ok &= fun
        print(f"      {A.name:18s} POLYNOMIAL base change functorial: {fun}"
              + ("" if fun else "   <-- UNEXPECTED"))
    print("      -> so the split condition is exactly the smooth/analytic divide, not a")
    print("         condition of the closure property.  Complexification is a functor on")
    print("         real-analytic manifolds and not on smooth ones; that, and only that,")
    print("         is what Kolar-Michor-Slovak's classification is recording.")
    return ok


def collisions():
    """W5: the collision trichotomy does not notice the boundary either."""
    ok = True
    print("\n  W5  C(A,B) = {F : F^A = F^B} with one side outside the Weil world")
    n, Ms = 2, 2
    for A, B in [(al.CC, al.D2), (al.CC, al.R2u), (al.CC, al.R2s), (al.CC, al.Rx2)]:
        res = []
        for deg, const in ((1, False), (1, True), (2, True)):
            Fm = gen_map(n, Ms, deg, const)
            res.append(al.lift_map(A, n, Fm, Ms) == al.lift_map(B, n, Fm, Ms))
        same = (A.c == B.c and A.one == B.one)
        exp = [True, same or A.one == B.one, same]
        ok &= (res == exp)
        lab = "all" if same else ("affine" if A.one == B.one else "linear")
        print(f"      {A.name:12s} vs {B.name:26s} " + "  ".join("Y" if r else "n" for r in res)
              + f"   {lab}" + ("" if res == exp else "   <-- UNEXPECTED"))
    print("      -> same trichotomy, same rule (equal units -> affine, else linear).  So")
    print("         affine rigidity now has a witness pair inside dimension 2 that is not")
    print("         two Weil algebras.")
    return ok


def _twist(A, n, d, Mt):
    return [sum((P.var(Mt, j * A.N + a) * d[g * A.N + a] for a in range(A.N)), P(Mt))
            for j in range(n) for g in range(A.N)]


def twists():
    """W6: the fields pi_A-related to X are far more than the lifts."""
    ok = True
    n, Ms = 2, 3
    X = gen_map(n, Ms, 3)
    psi = heun(X, n, Ms)
    print("\n  W6  beyond the lifts: X^A + delta_A for delta in Der(A) is pi_A-related to X")
    print(f"      {'A':26s} {'dim Der':8s} {'[d_A,X^A]=0':12s} {'pi o Psi = Psi o pi':20s} splits")
    for A in [al.D2, al.J3, al.W2, al.DD, al.CC]:
        m, Mt = n * A.N, n * A.N + 1
        XA = al.lift_map(A, n, X, Ms, hs=n); hA = P.var(Mt, m)
        Xr = al.lift_map(A, n, X, Ms)
        ders = [d for d in al.derivations(A)
                if not all(p.is_zero() for p in _twist(A, n, d, m))]
        com, base, spl = True, True, True
        for d in ders:
            dv = _twist(A, n, d, m)
            com &= all(p.is_zero() for p in
                       [sum((Xr[i].diff(k) * dv[k] - dv[i].diff(k) * Xr[k]
                             for k in range(m)), P(m)) for i in range(m)])
            dA = _twist(A, n, d, Mt)
            Z = [XA[i] + dA[i] for i in range(m)]
            pz = heun(Z, m, Mt)
            base &= ([pz[j * A.N] for j in range(n)] ==
                     [al.subs(p, [P.var(Mt, j * A.N) for j in range(n)] + [hA], Mt)
                      for p in psi])
            pxa, pdl = heun(XA, m, Mt), heun(dA, m, Mt)
            spl &= (pz == [al.subs(p, pdl + [hA], Mt) for p in pxa])
        exp_spl = (len(ders) == 0)
        ok &= com & base & (spl == exp_spl)
        print(f"      {A.name:26s} {len(ders):7d}  {str(com):12s} {str(base):20s} {spl}"
              + ("" if (com and base and spl == exp_spl) else "   <-- UNEXPECTED"))
    print("      -> Der(C) = 0: the twists are a nilpotent phenomenon, living exactly where")
    print("         jets do.  [delta_A, X^A] = 0 because X^A is Aut(A)-invariant, so the")
    print("         exact flows split as phi^{X^A}_t o T^{exp(t delta)}.  A method keeps the")
    print("         base half and does NOT split the pair -- exp(t delta) is transcendental")
    print("         (for D it is e^{ct} on the fibre), so the Der(A) part has to be taken")
    print("         exactly.  Aut(A)-equivariance is what supplies it; no new axiom.")
    return ok


def main():
    return boundary() & closure_outside() & liana_factor() & regularity() & collisions() & twists()


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
