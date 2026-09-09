"""What jet transport requires of a GENERAL method.

The condition Psi_{n dim A}(X^A) = T^A(Psi_n(X)) has an asymmetry that the
equivariance framing hides.  On the left, Psi receives X^A as a BARE field on
R^{n dim A}: it is not told the field is a lift, nor from which (n, A).  On the
right, T^A is applied to the output using the full structure.  So the property
says: a method blind to the structure must nonetheless return the structured
answer.

For a method that is an ALGORITHM OVER THE BASE RING this is automatic, and is
what jet transport does in practice -- run the same program with R-arithmetic
replaced by A-arithmetic.  The prerequisites for that are a short list, and this
module checks the two that are not obvious.

J1  Field access by evaluation only.  Evaluation at an A-point IS X^A, by
    definition of base change, so it costs nothing.
J2  Only R-algebra operations on values: +, -, real scalars, products.  No norms,
    no order, no branching on them, and DIVISION ONLY BY UNITS -- checked in J3
    below, since over A a general element is not invertible.
J3  Implicit definitions are fine when the defining equation is natural in the
    base and has a unique solution; nilpotence of m gives this whenever the
    R-linearization is invertible.  (This is why implicit RK jet-transports.)
J4  Any auxiliary structure must be NATURAL -- it must base-change functorially.
    A SPLITTING does.  A METRIC does not: that is exactly the c_A / balancedness /
    trace trichotomy of section 4.

The consequence, checked here on leapfrog: the property is not a condition on a
family indexed by plain dimension.  It is a condition relative to a SHAPE that
base-changes.  Leapfrog's shape is "dimension with a splitting", and the splitting
must be lifted the way base change lifts it -- one of three matchings, not the
other two.
"""
from fractions import Fraction as F
from poly import P
import algebra as al
from test_affine import gen_map


def leapfrog(X, m, Ms, part):
    """Stormer-Verlet / leapfrog for the splitting `part = (I1, I2)` of the
    coordinates: half-kick on I2, drift on I1, half-kick on I2.  The field is
    touched only by evaluation, and the only extra datum is the splitting."""
    I1, I2 = part
    h = P.var(Ms, m); u = [P.var(Ms, i) for i in range(m)]
    ev = lambda pt: [al.subs(p, list(pt) + [h], Ms) for p in X]
    Xu = ev(u)
    w = list(u)
    for i in I2: w[i] = u[i] + h * Xu[i] * F(1, 2)      # half kick
    Xw = ev(w)
    for i in I1: w[i] = u[i] + h * Xw[i]                # drift
    Xw2 = ev(w)
    for i in I2: w[i] = w[i] + h * Xw2[i] * F(1, 2)     # half kick
    return w


def shape():
    """J1/J4: leapfrog is closed for the LIFTED splitting and for no other matching."""
    ok = True
    n, Ms, A = 2, 3, al.D2
    N = A.N; m = n * N
    X = gen_map(n, Ms, 3)
    base = leapfrog(X, n, Ms, ([0], [1]))
    want = al.lift_map(A, n, base, Ms, hs=n)
    XA = al.lift_map(A, n, X, Ms, hs=n)
    print("  J1/J4  leapfrog on R^2 with splitting {0}|{1}, lifted by D to R^4")
    print("         (block convention: coordinate j becomes the A-block {jN..jN+N-1})")
    cases = [("lifted splitting   {0,1}|{2,3}", ([0, 1], [2, 3]), True),
             ("real|fibre         {0,2}|{1,3}", ([0, 2], [1, 3]), False),
             ("crossed            {0,3}|{1,2}", ([0, 3], [1, 2]), False)]
    for label, part, exp in cases:
        got = (leapfrog(XA, m, m + 1, part) == want)
        ok &= (got == exp)
        print(f"         {label}   closed: {str(got):5s}  expected {exp}"
              + ("" if got == exp else "   <-- UNEXPECTED"))
    print("         -> the splitting must be lifted as base change lifts it.  The")
    print("            property is relative to a SHAPE, not to plain dimension.")
    return ok


def not_equivariant():
    """Leapfrog is closed and NOT affine equivariant -- the splitting is not GL-invariant."""
    n, Ms = 2, 3
    X = gen_map(n, Ms, 3)
    g  = [[F(1), F(1)], [F(0), F(1)]]                    # a block-mixing shear
    gi = [[F(1), F(-1)], [F(0), F(1)]]
    push = [sum((al.subs(X[j], [sum((P.var(Ms, b) * gi[j2][b] for b in range(n)), P(Ms))
                                for j2 in range(n)] + [P.var(Ms, n)], Ms) * g[i][j]
                 for j in range(n)), P(Ms)) for i in range(n)]
    lhs = [sum((leapfrog(X, n, Ms, ([0], [1]))[j] * g[i][j] for j in range(n)), P(Ms))
           for i in range(n)]
    rhs = [al.subs(p, [sum((P.var(Ms, b) * gi[j][b] for b in range(n)), P(Ms))
                       for j in range(n)] + [P.var(Ms, n)], Ms)
           for p in leapfrog(push, n, Ms, ([0], [1]))]
    rhs = [sum((rhs[j] * g[i][j] for j in range(n)), P(Ms)) for i in range(n)]
    eq = (lhs == rhs)
    print(f"\n  J2     leapfrog affine equivariant under a block-mixing shear: {eq}")
    print("         -> closed and not equivariant: the recorded factorisation")
    print("            [equivariance] n [closure] is not the only way to characterise.")
    return not eq


def units():
    """J2/J3: over A, division is available only by UNITS -- a is invertible iff its
    image in each local factor's residue field is nonzero.  Checked by solving ax = 1."""
    ok = True
    print("\n  J3     invertibility in A (division is legal only by units)")
    tests = [(al.D2, [F(1), F(5)], True, "1 + 5e"), (al.D2, [F(0), F(1)], False, "e"),
             (al.R2s, [F(2), F(3)], True, "(2,3)"), (al.R2s, [F(1), F(0)], False, "(1,0)"),
             (al.J3, [F(0), F(1), F(2)], False, "e + 2e^2"),
             (al.J3, [F(3), F(1), F(0)], True, "3 + e")]
    for A, a, exp, lbl in tests:
        rows = [[sum(a[al_] * A.c[al_][be][g] for al_ in range(A.N)) for be in range(A.N)]
                for g in range(A.N)]
        solvable = _solvable(rows, [A.one[g] for g in range(A.N)], A.N)   # solve a.x = 1
        ok &= (solvable == exp)
        print(f"         {A.name:26s} a = {lbl:12s} invertible: {str(solvable):5s} "
              f"expected {exp}" + ("" if solvable == exp else "   <-- UNEXPECTED"))
    print("         -> a method that divides by a field-dependent quantity jet-transports")
    print("            only where that quantity stays a unit; RK never divides, which is")
    print("            why the prerequisite is invisible for it.")
    return ok


def _solvable(M, b, w):
    R = [M[i][:] + [b[i]] for i in range(len(M))]
    r = 0
    for c in range(w):
        p = next((i for i in range(r, len(R)) if R[i][c]), None)
        if p is None: continue
        R[r], R[p] = R[p], R[r]; pv = R[r][c]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c] / pv; R[i] = [R[i][j] - f * R[r][j] for j in range(w + 1)]
        r += 1
    return not any(all(R[i][c] == 0 for c in range(w)) and R[i][w] != 0 for i in range(len(R)))


def main():
    return shape() & not_equivariant() & units()


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
