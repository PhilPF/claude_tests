"""WHAT JET TRANSPORT IS -- and what it is NOT.

A first version of this module identified jet transport with the closure property
itself ("a method jet-transports iff it is a term in the C^infty-theory generated
by field evaluation and the slot structure of a base-changing shape").  That is
WRONG, and T6 below refutes it on two concrete, unremarkable methods.  The error was
to fuse two different statements.  Kept separate they read:

  TRANSPORT is an OPERATION on programs, not a property of methods.  Given an
  algorithm whose every step is a C^infty-ring operation, it has an interpretation
  over any Weil algebra A -- the state ranging over A^n, n unchanged.  That is what
  jet arithmetic runs.  Its correctness needs no hypothesis beyond being defined:

        P<A>  =  T^A (P<R>)      always,     by functoriality.               [T6]

  So there is no interesting class of "methods that jet-transport": EVERY
  C^infty-typed program does, and being C^infty-typed is a property of the code.

  CLOSURE is a COHERENCE condition on a FAMILY {Psi_d}, comparing two different ways
  to grow: raising the ring at fixed dimension (transport) against raising the
  dimension over R (the family's own member at d = n dim A).  With the realization
  iota : R^{n dim A} -> A^n,

        Psi_{n dim A}(iota^* X^A)  =  iota^*( Psi_n<A> ).

  Jet transport supplies the RIGHT-hand side and makes the equation meaningful.
  Closure is the equation.  Confusing the two makes the property look automatic.

That relocates J1-J4.  J1-J3 are prerequisites for TRANSPORT: conditions on the code,
so that P<A> exists at all.  J4 is not a transport condition -- a metric base-changes
perfectly well -- it is the UNIFORMITY of the family in the dimension direction: the
structure the method uses at dimension n dim A must be the realization of the one it
uses at n over A.  That is where c_A, balancedness, the trace and the splitting live,
and it is section 4's trichotomy, not a prerequisite for jet arithmetic to run.

A Weil algebra IS a C^infty-ring -- a set with an operation f_R : R^n -> R^k for every
smooth f, compatible with composition, identities and products -- and base change IS
its C^infty-structure:

    f^A(u) = sum_alpha  d^alpha f(pi u) / alpha!  (u - pi u)^alpha,          (*)

a FINITE sum because u - pi(u) is nilpotent; and (*) is exactly T^A f.  Everything
transport needs is in that sentence.

  T1  base change is the C^infty-structure: the Taylor formula (*), and functoriality
  T2  the expansion point is forced -- one real point per local factor
  T3  division by units is a C^infty operation, agreeing with the algebraic inverse
  T4  implicit definitions transport: no hypothesis on implicitness is needed
  T5  the coherence side: the shape is an action of the (x)-monoid
  T6  transport vs closure, on the two methods the old definition got wrong
"""
from fractions import Fraction as F
from itertools import product
from poly import P, deriv
import algebra as al
from test_affine import gen_map
from test_jettransport import leapfrog


def _factorial(a):
    f = 1
    for e in a:
        for i in range(2, e + 1): f *= i
    return F(f)


def _multi(n, d):
    """all multi-indices in n slots of total degree <= d"""
    for k in product(range(d + 1), repeat=n):
        if sum(k) <= d: yield k


def _unit_first_local(A):
    """A is local with 1_A = e_0, so m = span(e_1..e_{N-1}) is the nilradical."""
    return (al.num_real_points(A) == 1
            and A.one == [F(1)] + [F(0)] * (A.N - 1))


def taylor():
    """T1: f^A(u) = sum_a d^a f(pi u) (u - pi u)^a / a!, and (g o f)^A = g^A o f^A."""
    ok = True
    n, deg = 2, 3
    Ms = n
    f = gen_map(n, Ms, deg)
    g = [f[1], f[0] + f[1]]
    gof = [al.subs(p, f, Ms) for p in g]
    print("  T1  base change IS the C^infty-ring structure of A")
    print("      Taylor formula about the real part vs. substitution;  functoriality")
    print("      A                                 (*) holds   (g o f)^A = g^A o f^A")
    for d in sorted(al.ALGEBRAS):
        for A in al.ALGEBRAS[d]:
            Mt = n * A.N
            u = al.generic_point(A, n, Mt)
            sub = al.basechange(A, n, f, u, Ms, Mt)          # substitution = T^A f
            tay = None
            if _unit_first_local(A):
                a = [P.var(Mt, j * A.N) for j in range(n)]    # pi(u_j), the real part
                nu = [[P(Mt)] + u[j][1:] for j in range(n)]   # u_j - pi(u_j), nilpotent
                tay = [[P(Mt)] * A.N for _ in range(n)]
                for k in _multi(n, deg):
                    e = al.aone(A, Mt)
                    for j in range(n):
                        for _ in range(k[j]): e = al.amul(A, e, nu[j], Mt)
                    for i in range(n):
                        c = al.subs(deriv(f[i], sum(([j] * k[j] for j in range(n)), [])),
                                    a, Mt) * (1 / _factorial(k))
                        tay[i] = [tay[i][gm] + e[gm] * c for gm in range(A.N)]
                ok &= (tay == sub)
            L = al.lift_map(A, n, f, Ms)
            Lg = al.lift_map(A, n, g, Ms)
            fun = ([al.subs(p, L, Mt) for p in Lg] == al.lift_map(A, n, gof, Ms))
            ok &= fun
            print(f"      {A.name:32s} {('yes ' if tay == sub else 'NO  ') if tay is not None else 'n/a ':10s}  {fun}"
                  + ("" if (fun and (tay is None or tay == sub)) else "   <-- UNEXPECTED"))
    print("      -> (*) is well defined and functorial: A is a C^infty-ring and T^A|_CartSp")
    print("         is its structure.  'A-arithmetic' is the polynomial part of it.")
    return ok


def expansion_point():
    """T2: (*) terminates only at the real part, and there is one PER LOCAL FACTOR."""
    ok = True
    print("\n  T2  the expansion point in (*) is forced, and counted by A/Nil(A)")
    print("      A                                dim  dim Nil  real pts  span(e_1..) nil")
    for d in sorted(al.ALGEBRAS):
        for A in al.ALGEBRAS[d]:
            nil = al.nilradical(A)
            ok &= all(al.is_nilpotent(A, v) for v in nil)      # self-check of rad(tr)
            k = al.num_real_points(A)
            comp = all(al.is_nilpotent(A, [F(1) if b == a else F(0) for b in range(A.N)])
                       for a in range(1, A.N))
            ok &= (comp == _unit_first_local(A))
            print(f"      {A.name:32s} {A.N:3d}  {len(nil):7d}  {k:8d}  {comp}")
    print("      -> for POLYNOMIAL f any expansion point works (Taylor is exact), so the")
    print("         datum is invisible; for smooth f the series terminates only when")
    print("         u - a is nilpotent, i.e. a = pi(u) in each factor.  A product algebra")
    print("         therefore carries SEVERAL real points and no single one will do:")
    a = [F(0), F(1)]                                            # (0,1) in R^2
    bad = [al.is_nilpotent(al.R2s, [a[0] - c, a[1] - c]) for c in (F(0), F(1), F(2))]
    ok &= not any(bad)
    print(f"         in R^2 at u = (0,1), u - c*1 nilpotent for c = 0,1,2: {bad}")
    print("         (each c fixes one factor and breaks the other -- the recorded trap)")
    return ok


def reciprocal():
    """T3: 1/x is a C^infty operation on units, and its lift is the algebraic inverse."""
    ok = True
    print("\n  T3  division by a unit is a C^infty operation")
    print("      the Taylor lift of 1/x at pi(a) equals the algebraic inverse a^{-1}")
    M = 1
    for d in sorted(al.ALGEBRAS):
        for A in al.ALGEBRAS[d]:
            if not _unit_first_local(A): continue
            a0 = F(2)
            nu = [F(0)] + [F(1, 3 + i) for i in range(A.N - 1)]
            a = [a0] + nu[1:]
            nuP = [P.const(M, c) for c in nu]
            r, pw = [P(M)] * A.N, al.aone(A, M)
            for k in range(A.N):                                # sum (-1)^k nu^k/a0^{k+1}
                s = F(-1) ** k / a0 ** (k + 1)
                r = [r[g] + pw[g] * P.const(M, s) for g in range(A.N)]
                pw = al.amul(A, pw, nuP, M)
            prod = al.amul(A, [P.const(M, c) for c in a], r, M)
            got = (prod == [P.const(M, c) for c in A.one])
            ok &= got
            print(f"      {A.name:32s} a = 2 + nu   a * (lift of 1/x) = 1: {got}"
                  + ("" if got else "   <-- UNEXPECTED"))
    print("      -> so J2's 'division only by units' is not an extra rule: 1/x is smooth")
    print("         exactly on the units, and (*) computes the inverse there.")
    return ok


def implicit():
    """T4: the implicit stage equation transports -- J3 is a theorem, not an axiom."""
    ok = True
    n, Ms, hs, K = 2, 3, 2, 2
    X = gen_map(n, Ms, 3)
    h = P.var(Ms, hs)
    u = [P.var(Ms, j) for j in range(n)]

    def picard(X, m, Ms, hs, K):
        """K steps of z <- u + (h/2) X(z) from z = u; each step is a program."""
        h = P.var(Ms, hs); z = [P.var(Ms, j) for j in range(m)]
        for _ in range(K):
            Xz = [al.subs(p, z + [h], Ms) for p in X]
            z = [P.var(Ms, j) + h * Xz[j] * F(1, 2) for j in range(m)]
        return z

    zr = picard(X, n, Ms, hs, K)
    psi = [u[j] + h * al.subs(p, zr + [h], Ms) for j, p in enumerate(X)]
    print("\n  T4  implicit definitions transport (implicit midpoint, K = 2 Picard steps)")
    print("      the STAGE equation z = u + (h/2) X(z), solved over A:")
    heavy = {"D = R[e]/(e^2)", "R^2 (standard basis)", "D (x) D"}
    for A in [al.D2, al.J3, al.W2, al.R2s, al.DxR, al.tensor(al.D2, al.D2, "D (x) D")]:
        m, Mt = n * A.N, n * A.N + 1
        XA = al.lift_map(A, n, X, Ms, hs=hs)
        z1 = picard(XA, m, Mt, m, K)
        st = (z1 == al.lift_map(A, n, zr, Ms, hs=hs))
        good = st
        line = f"      {A.name:26s} z^A = z_A: {st}"
        if A.name in heavy:                      # the full step, where it is affordable
            hA = P.var(Mt, m)
            p1 = [P.var(Mt, i) + hA * al.subs(p, z1 + [hA], Mt) for i, p in enumerate(XA)]
            fu = (p1 == al.lift_map(A, n, psi, Ms, hs=hs))
            good &= fu
            line += f"      Psi^A = Psi_A: {fu}"
        ok &= good
        print(line + ("" if good else "   <-- UNEXPECTED"))
    print("      -> each Picard step is a term, so the iteration transports; nilpotence")
    print("         makes it terminate and uniqueness makes the limit the base change.")
    print("         Nothing about implicitness had to be assumed.  (The full step is")
    print("         checked on three algebras; degree 27 in 7 slots is the cost ceiling.)")
    return ok


def shape_action():
    """T5: sigma^{A(x)B} = (sigma^B)^A -- the block convention is an ACTION, not just
    a per-algebra choice, and leapfrog is closed for the composite."""
    ok = True
    n, Ms = 2, 3
    A = B = al.D2
    C = al.tensor(A, B)
    print("\n  T5  the shape is an action of the (x)-monoid, not a per-algebra choice")

    f = gen_map(n, Ms, 3)
    once = al.lift_map(C, n, f, Ms, hs=n)
    twiceB = al.lift_map(B, n, f, Ms, hs=n)
    twice = al.lift_map(A, n * B.N, twiceB, n * B.N + 1, hs=n * B.N)
    nose = (once == twice)
    ok &= nose
    print(f"      T^(A(x)B) = T^A o T^B on the nose in the block realization: {nose}"
          + ("" if nose else "   <-- UNEXPECTED"))

    idx = all((j * C.N + (be * A.N + a)) == ((j * B.N + be) * A.N + a)
              for j in range(n) for a in range(A.N) for be in range(B.N))
    ok &= idx
    print(f"      slot cocycle  j*N_C + (be*N_A + a) = (j*N_B + be)*N_A + a: {idx}")

    m = n * C.N
    X = gen_map(n, Ms, 3)
    XC = al.lift_map(C, n, X, Ms, hs=n)
    want = al.lift_map(C, n, leapfrog(X, n, Ms, ([0], [1])), Ms, hs=n)
    blk = lambda I: [j * C.N + g for j in I for g in range(C.N)]
    cases = [("lifted twice   {0..3}|{4..7}", (blk([0]), blk([1])), True),
             ("A-blocks swapped               ", ([0, 1, 4, 5], [2, 3, 6, 7]), False),
             ("real | rest                    ", ([0, 4], [1, 2, 3, 5, 6, 7]), False)]
    print("      leapfrog on R^2, splitting {0}|{1}, lifted by D (x) D to R^8:")
    for lbl, part, exp in cases:
        got = (leapfrog(XC, m, m + 1, part) == want)
        ok &= (got == exp)
        print(f"        {lbl:34s} closed: {str(got):5s}  expected {exp}"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    print("      -> lifting the splitting twice by D gives exactly the D(x)D splitting, so")
    print("         the shape is a genuine action and 'closed over A and over B' composes.")
    return ok


def _dA(A, p, j, Mt):
    """the A-derivative in slot j of a realized polynomial: the directional
    derivative along 1_A, since d_{A,j}F . 1_A = D_{(j,1_A)}F."""
    out = P(Mt)
    for a in range(A.N):
        if A.one[a]: out = out + p.diff(j * A.N + a) * A.one[a]
    return out


def transport_vs_closure():
    """T6: transport is unconditional; closure is the coherence equation.  Run on the
    two methods that refute the earlier definition -- both use auxiliary structure
    that does NOT base-change uniformly in dimension (a metric, a trace), so the
    earlier definition denied them jet transport.  They transport perfectly."""
    ok = True
    n, Ms, hs = 2, 3, 2
    X = gen_map(n, Ms, 3)
    h = P.var(Ms, hs); u = [P.var(Ms, j) for j in range(n)]
    lap = [sum((X[i].diff(j).diff(j) for j in range(n)), P(Ms)) for i in range(n)]
    div = sum((X[j].diff(j) for j in range(n)), P(Ms))
    psiL = [u[i] + h * X[i] + h * h * lap[i] for i in range(n)]          # Laplacian
    psiA = [u[i] + h * X[i] + h * h * div * X[i] for i in range(n)]      # aromatic
    print("\n  T6  TRANSPORT is unconditional; CLOSURE is the coherence equation")
    print("      the same program, run two ways: at dimension n over A (transport),")
    print("      and at the realized dimension n*dim A over R (the family's member)")
    print(f"      {'algebra':30s} {'c_A = 1':8s} {'Laplacian':19s} {'aromatic'}")
    print(f"      {'':30s} {'':8s} {'transp':9s} {'dim nN':9s} {'transp':9s} {'dim nN'}")
    for d in sorted(al.ALGEBRAS):
        for A in al.ALGEBRAS[d]:
            N = A.N; m = n * N; Mt = m + 1
            XA = al.lift_map(A, n, X, Ms, hs=hs)
            hA = P.var(Mt, m); uu = [P.var(Mt, i) for i in range(m)]
            # the same program at dimension n, with the state ranging over A
            lapA = [sum((_dA(A, _dA(A, XA[i], j, Mt), j, Mt) for j in range(n)), P(Mt))
                    for i in range(m)]
            tL = [uu[i] + hA * XA[i] + hA * hA * lapA[i] for i in range(m)]
            dvA = [sum((_dA(A, XA[j * N + a], j, Mt) for j in range(n)), P(Mt))
                   for a in range(N)]
            pr = [al.amul(A, dvA, [XA[j * N + a] for a in range(N)], Mt) for j in range(n)]
            tA = [uu[j * N + a] + hA * XA[j * N + a] + hA * hA * pr[j][a]
                  for j in range(n) for a in range(N)]
            # the family's own member at the realized dimension, over R
            lapR = [sum((XA[i].diff(k).diff(k) for k in range(m)), P(Mt)) for i in range(m)]
            rL = [uu[i] + hA * XA[i] + hA * hA * lapR[i] for i in range(m)]
            dvR = sum((XA[k].diff(k) for k in range(m)), P(Mt))
            rA = [uu[i] + hA * XA[i] + hA * hA * dvR * XA[i] for i in range(m)]
            wL = al.lift_map(A, n, psiL, Ms, hs=hs)
            wA = al.lift_map(A, n, psiA, Ms, hs=hs)
            c1 = (al.c_A(A) == A.one)
            got = (tL == wL, rL == wL, tA == wA, rA == wA)
            exp = (True, c1, True, A.N == 1)          # c_A = 1 for lianas; dim 1 for loops
            ok &= (got == exp)
            print(f"      {A.name:30s} {str(c1):8s} {str(got[0]):9s} {str(got[1]):9s} "
                  f"{str(got[2]):9s} {got[3]}" + ("" if got == exp else "   <-- UNEXPECTED"))
    print("      -> transport: True on all 17, for BOTH methods.  A metric and a trace are")
    print("         perfectly transportable; jet arithmetic runs them and returns T^A Psi_n.")
    print("      -> closure: the Laplacian agrees iff c_A = 1_A, the aromatic one only for")
    print("         A = R -- exactly section 4's liana and loop laws, here DERIVED as the")
    print("         failure of the family to be uniform in the dimension direction.")
    print("      -> so 'uses only base-changing structure' is not a condition for jet")
    print("         transport.  The earlier definition in this file denied transport to")
    print("         two methods that transport; it was a definition of the wrong thing.")
    return ok


def main():
    return (taylor() & expansion_point() & reciprocal() & implicit()
            & shape_action() & transport_vs_closure())


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
