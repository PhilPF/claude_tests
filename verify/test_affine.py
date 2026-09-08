"""Affine rigidity: the per-dimension consequence of algebraic naturality, and why
it is the last one obtainable by comparing two lifts of one field.

A1.  TRICHOTOMY.  For based algebras (A,iota), (B,kappa) of the same dimension N,
     put C(A,B) = { F : F^A = F^B }.  Then

         C(A,B) = all maps        if the based algebras coincide,
                = the AFFINE maps if they differ but 1_A = 1_B,
                = the LINEAR maps if 1_A != 1_B.

     The same set of conditions applies to a vector field and to a map, so a
     coincidence of lifts transfers verbatim from X to Psi^X.

A2.  AFFINE RIGIDITY.  Hence, for algebraically natural Psi: if X is an affine
     field then Psi^X_h is an affine map.  Cleanest witness pair, in dimension 2:
     the D-lift of X realizes as (X(u0), DX(u0)u1) -- a tangent -- and the
     R^2-lift in the unit-first basis (1,p) as (X(u0), X(u0+u1)-X(u0)) -- a
     secant.  Secant = tangent exactly on affine maps.  In dimension 3 the pair
     (R[e]/(e^3), R[x,y]/m^2) does the same with both algebras in their standard
     bases, so no basis convention is needed.

A3.  SHARPNESS.  A2 is strictly stronger than the known linear rigidity (a local
     method satisfying the latter and failing the former is exhibited), and no
     comparison of two lifts of a single field can reach beyond affine: distinct
     based algebras already differ on some product of two basis vectors.

A4.  CROSS-FACTORIZATION.  Comparisons across DIFFERENT source dimensions
     (X^A = Y^B with X on R^{n1}, Y on R^{n2}, n1 != n2) add nothing: over the
     enumerated range every such coincidence is either vacuous (B = A (x) C with
     Y^B = (Y^C)^A, so both conditions are the same condition) or confined to the
     affine maps again.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial
from poly import P
import algebra as al
from algebra import ALGEBRAS, D2, R2s, R2u, J3, W2, DxR, R3u, R3s


# ---- test maps: generic in each degree, with a cubic term (quadratic fields
#      degenerate -- standing trap) --------------------------------------------

def gen_map(n, Ms, deg, const=True):
    u = [P.var(Ms, j) for j in range(n)]
    c = [F(2, 3), F(-5, 7), F(4, 5)][:n]
    M = [[F(1, 2), F(-2), F(3, 5)][:n], [F(3), F(5, 4), F(-1, 3)][:n], [F(7, 2), F(-4, 9), F(2)][:n]][:n]
    out = [sum((u[j] * M[i][j] for j in range(n)), P(Ms)) for i in range(n)]
    if const: out = [out[i] + P.const(Ms, c[i]) for i in range(n)]
    if deg >= 2:
        out = [out[i] + u[i % n] * u[(i + 1) % n] * F(3, 2) + u[0] * u[0] * F(-1, 5) for i in range(n)]
    if deg >= 3:
        out = [out[i] + u[0] * u[(i + 1) % n] * u[n - 1] * F(2, 7) for i in range(n)]
    return out


def deg_u(p, n):
    return max((sum(k[:n]) for k in p.d), default=-1)


# ---- A1: the trichotomy -------------------------------------------------------

def trichotomy():
    ok = True
    print("  A1  C(A,B) = {F : F^A = F^B}, over pairs of based algebras of equal dimension")
    print("      N  A                              B                              lin  aff  quad predicted")
    for N, lst in ALGEBRAS.items():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                A, B = lst[i], lst[j]
                n, Ms = 2, 2
                res = []
                for deg, const in ((1, False), (1, True), (2, True)):
                    Fm = gen_map(n, Ms, deg, const)
                    res.append(al.lift_map(A, n, Fm, Ms) == al.lift_map(B, n, Fm, Ms))
                same = (A.c == B.c and A.one == B.one)
                pred = [True, same or A.one == B.one, same]
                good = res == pred
                ok &= good
                lab = "all" if same else ("affine" if A.one == B.one else "linear")
                print(f"      {N}  {A.name:30s} {B.name:30s} "
                      + "  ".join("Y" if r else "n" for r in res)
                      + f"    {lab}" + ("" if good else "   <-- UNEXPECTED"))
    print("      -> C(A,B) is: all maps / the affine maps / the linear maps, according as")
    print("         the based algebras agree / differ with equal units / have different units")
    return ok


# ---- methods, as polynomial maps in u and h ----------------------------------

def ev(X, args, Ms, n):
    im = list(args) + [P.var(Ms, n)]
    return [al.subs(p, im, Ms) for p in X]


def _u(n, Ms): return [P.var(Ms, j) for j in range(n)]


def euler(X, n, Ms):
    h = P.var(Ms, n); u = _u(n, Ms)
    return [u[i] + h * X[i] for i in range(n)]


def heun(X, n, Ms):
    h = P.var(Ms, n); u = _u(n, Ms)
    k2 = ev(X, [u[i] + h * X[i] for i in range(n)], Ms, n)
    return [u[i] + (X[i] + k2[i]) * F(1, 2) * h for i in range(n)]


def midpoint(X, n, Ms):
    h = P.var(Ms, n); u = _u(n, Ms)
    k2 = ev(X, [u[i] + h * X[i] * F(1, 2) for i in range(n)], Ms, n)
    return [u[i] + h * k2[i] for i in range(n)]


def rk4(X, n, Ms):
    h = P.var(Ms, n); u = _u(n, Ms)
    k1 = X
    k2 = ev(X, [u[i] + h * k1[i] * F(1, 2) for i in range(n)], Ms, n)
    k3 = ev(X, [u[i] + h * k2[i] * F(1, 2) for i in range(n)], Ms, n)
    k4 = ev(X, [u[i] + h * k3[i] for i in range(n)], Ms, n)
    return [u[i] + h * (k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) * F(1, 6) for i in range(n)]


def taylor2(X, n, Ms):
    h = P.var(Ms, n); u = _u(n, Ms)
    dxx = [sum((X[i].diff(j) * X[j] for j in range(n)), P(Ms)) for i in range(n)]
    return [u[i] + h * X[i] + h * h * dxx[i] * F(1, 2) for i in range(n)]


def laplacian(X, n, Ms):
    """u + hX + h^2 (Delta X): (T)-natural, not a B-series method."""
    h = P.var(Ms, n); u = _u(n, Ms)
    lap = [sum((X[i].diff(j).diff(j) for j in range(n)), P(Ms)) for i in range(n)]
    return [u[i] + h * X[i] + h * h * lap[i] for i in range(n)]


def ring(X, n, Ms):
    """u + hX + h^2 (X^1)^2 X: the round-3 (T)-natural non-equivariant method."""
    h = P.var(Ms, n); u = _u(n, Ms)
    s = X[0] * X[0]
    return [u[i] + h * X[i] + h * h * s * X[i] for i in range(n)]


def secant(X, n, Ms):
    """u + hX + h^2 <X - DX.u, u> X.  The bracket vanishes identically on LINEAR
    fields (Euler's identity), so this satisfies linear rigidity; on an affine
    field it is <c,u>, and the method is quadratic in u."""
    h = P.var(Ms, n); u = _u(n, Ms)
    c = [X[i] - sum((X[i].diff(j) * u[j] for j in range(n)), P(Ms)) for i in range(n)]
    q = sum((c[i] * u[i] for i in range(n)), P(Ms))
    return [u[i] + h * X[i] + h * h * q * X[i] for i in range(n)]


METHODS = [("explicit Euler", euler, True, True), ("Heun (RK2)", heun, True, True),
           ("midpoint", midpoint, True, True), ("RK4", rk4, True, True),
           ("Taylor-2", taylor2, True, True),
           ("Laplacian  u+hX+h^2 Delta X", laplacian, True, True),
           ("ring  u+hX+h^2 (X^1)^2 X", ring, False, False),
           ("secant  u+hX+h^2 <X-DX.u,u> X", secant, True, False)]


# ---- A2/A3: rigidity of the methods ------------------------------------------

def rigidity():
    ok = True
    n, Ms = 2, 3                                    # variables 0,1 = u; 2 = h
    lin = gen_map(n, Ms, 1, const=False)
    aff = gen_map(n, Ms, 1, const=True)
    cub = gen_map(n, Ms, 3, const=True)
    a = [F(3, 4), F(-2, 5)]
    sh = [P.var(Ms, 0) + P.const(Ms, a[0]), P.var(Ms, 1) + P.const(Ms, a[1]), P.var(Ms, n)]
    print("\n  A2  linear rigidity (Psi^M linear) vs affine rigidity (Psi^X affine)")
    print("      method                          lin-rigid  aff-rigid  transl-eqv  predicted")
    for name, meth, plin, paff in METHODS:
        rl = all(deg_u(p, n) == 1 for p in meth(lin, n, Ms))
        ra = all(deg_u(p, n) <= 1 for p in meth(aff, n, Ms))
        # translation equivariance: Psi^{X(.+a)}(u) = Psi^X(u+a) - a.  Checked on the
        # single-evaluation methods, which is where the claim carries information; for
        # a multi-stage RK method it is evident from the stage equations, and the
        # composed polynomials are of degree 3^4 on a cubic field.
        if meth in (heun, midpoint, rk4):
            rt = None
        else:
            Xa = [al.subs(q, sh, Ms) for q in cub]
            rt = (meth(Xa, n, Ms)
                  == [al.subs(q, sh, Ms) - P.const(Ms, a[i])
                      for i, q in enumerate(meth(cub, n, Ms))])
        good = (rl, ra) == (plin, paff) and not (rt and rl and not ra)
        ok &= good
        print(f"      {name:32s}  {str(rl):5s}      {str(ra):5s}      "
              f"{('-' if rt is None else str(rt)):5s}      ({plin}, {paff})" + ("" if good else "  <-- UNEXPECTED"))
    print("      -> 'secant' satisfies linear rigidity and fails affine rigidity:")
    print("         affine rigidity is STRICTLY stronger than linear rigidity as a condition.")
    print("      -> but it is not translation equivariant, and cannot be: for a local method")
    print("         u + Phi(j^r X(u), h) the affine field Mu+c has the SAME jet slots as the")
    print("         linear field, only a shifted value, so linear rigidity already forces")
    print("         Phi(.,M,0;h) linear and hence Psi^{Mu+c} affine.  The gap between the two")
    print("         rigidities lives exactly where translation equivariance is not assumed.")

    print("\n      the derivation, run on the two witness pairs:")
    for A, B, tag in ((D2, R2u, "dim 2: tangent vs secant, R^2 in the unit-first basis"),
                      (J3, W2, "dim 3: both algebras in their standard bases")):
        same_field = al.lift_map(A, n, aff, Ms, hs=n) == al.lift_map(B, n, aff, Ms, hs=n)
        print(f"      [{A.name} | {B.name}]  {tag}")
        print(f"        X^A == X^B on the affine field: {same_field}")
        ok &= same_field
        for name, meth, plin, paff in METHODS:
            psi = meth(aff, n, Ms)
            eq = al.lift_map(A, n, psi, Ms, hs=n) == al.lift_map(B, n, psi, Ms, hs=n)
            ok &= (eq == paff)
            if not paff:
                print(f"        T^A Psi^X == T^B Psi^X for '{name}': {eq}"
                      + ("  -> not algebraically natural" if not eq else "  <-- UNEXPECTED"))
    print("        (a single Psi_{3n} cannot serve both algebras, since X^A = X^B is one field)")
    return ok


# ---- A4: cross-factorization --------------------------------------------------

def rank(rows):
    rows = [r[:] for r in rows if any(r)]
    w = len(rows[0]) if rows else 0
    r = 0
    for c in range(w):
        piv = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][c]
        for i in range(r + 1, len(rows)):
            if rows[i][c]:
                f = rows[i][c] / pv
                rows[i] = [rows[i][j] - f * rows[r][j] for j in range(w)]
        r += 1
    return r


def homog_basis(n, k):
    """monomial basis of the degree-k homogeneous polynomial maps R^n -> R^n."""
    mons = []
    def rec(j, left, e):
        if j == n - 1: mons.append(tuple(e + [left])); return
        for a in range(left + 1): rec(j + 1, left - a, e + [a])
    rec(0, k, [])
    return [(i, m) for i in range(n) for m in mons]


def cross_kernel(B, n2, A, n1, k):
    """Basis of { Y homogeneous of degree k on R^{n2} : Y^B is the A-lift of a field
    on R^{n1} }.  Both Y -> Y^B and the A-lift residual are linear, so this is a
    kernel over Q; it is returned as an explicit basis, not just a dimension."""
    basis = homog_basis(n2, k)
    Ys = [[P(n2, {mon: F(1)}) if t == i else P(n2) for t in range(n2)] for (i, mon) in basis]
    cols, keys = [], {}
    for Y in Ys:
        col = {}
        for t, p in enumerate(al.lift_residual(A, n1, al.lift_map(B, n2, Y, n2))):
            for e, v in p.d.items():
                key = (t, e); keys.setdefault(key, len(keys)); col[key] = v
        cols.append(col)
    rows = [[cols[j].get(key, F(0)) for j in range(len(basis))] for key in keys]
    ker = kernel(rows, len(basis))
    return [[sum((Ys[j][t] * c for j, c in enumerate(v) if c), P(n2)) for t in range(n2)]
            for v in ker], len(basis)


def kernel(rows, w):
    """basis of the null space of the matrix rows (w columns), over Q."""
    R = [r[:] for r in rows]
    piv, r = [], 0
    for c in range(w):
        i = next((i for i in range(r, len(R)) if R[i][c]), None)
        if i is None: continue
        R[r], R[i] = R[i], R[r]
        R[r] = [x / R[r][c] for x in R[r]]
        for i2 in range(len(R)):
            if i2 != r and R[i2][c]:
                f = R[i2][c]; R[i2] = [R[i2][j] - f * R[r][j] for j in range(w)]
        piv.append(c); r += 1
    out = []
    for c in range(w):
        if c in piv: continue
        v = [F(0)] * w; v[c] = F(1)
        for i, pc in enumerate(piv): v[pc] = -R[i][c]
        out.append(v)
    return out


def certificate(A, n1, B, n2, Y):
    """Is the coincidence X^A = Y^B forced by functoriality?  Two ways: X is itself
    a C-lift of Y (so B = A (x) C and both conditions are T^A applied to the same
    one), or X and Y are lifts of one common field W (so both conditions are the
    condition at W).  Either way the relation holds for every algebraically natural
    Psi and constrains nothing."""
    X = al.recognise(A, n1, al.lift_map(B, n2, Y, n2))
    if X is None: return None
    if n2 and n1 % n2 == 0:
        for C in ALGEBRAS.get(n1 // n2, []):
            if C.N > 1 and al.recognise(C, n2, X) == Y: return "X = Y^C"
    for C in ALGEBRAS.get(n2, []):
        W = al.recognise(C, 1, Y)
        if W is None: continue
        for Cp in ALGEBRAS.get(n1, []):
            if al.recognise(Cp, 1, X) == W: return "X, Y lift one W"
    return None


def cross():
    ok = True
    print("\n  A4  cross-factorization: X^A = Y^B with n1 != n2 -- can it reach a non-affine X?")
    print("      m  (n2, B)                          (n1, A)                    ker(k=2) ker(k=3)  verdict")
    for m, n2, Bs, n1, As in [(4, 1, ALGEBRAS[4], 2, ALGEBRAS[2]),
                              (6, 2, ALGEBRAS[3], 3, ALGEBRAS[2])]:
        for B in Bs:
            for A in As:
                dims, verdict = [], "no coincidence beyond affine"
                for k in (2, 3):
                    ker, full = cross_kernel(B, n2, A, n1, k)
                    dims.append(len(ker))
                    if not ker: continue
                    cert = {certificate(A, n1, B, n2, Y) for Y in ker}
                    # a certificate for each basis vector, and the certified
                    # subspace (all of it, or the one-variable lifts) exhausts the kernel
                    if None in cert:
                        verdict = "UNCERTIFIED"; ok = False
                    else:
                        c = cert.pop()
                        exp = full if c == "X = Y^C" else 1
                        verdict = f"functorial ({c})"
                        if len(ker) != exp: verdict = "UNCERTIFIED"; ok = False
                print(f"      {m}  ({n2}, {B.name:28s}) ({n1}, {A.name:22s}) "
                      f"{dims[0]:5d}    {dims[1]:5d}     {verdict}")
    print("      -> every cross coincidence carrying a non-affine field is functorial: the two")
    print("         conditions are the same condition, so they constrain nothing.  At m = 6 no")
    print("         tensor factorization exists (dim C = 3/2) and the only coincidences are")
    print("         those where X and Y are both lifts of one field on R.")
    return ok


def compositions(k, j):
    if j == 1:
        yield (k,); return
    for a in range(1, k - j + 2):
        for rest in compositions(k - a, j - 1): yield (a,) + rest


def variational():
    """A0.  The A_r-base-changed field IS the tower of higher variational equations.

    For A_r = R[e]/(e^{r+1}) and u = sum_k e^k u_k, expanding X(u) in A_r gives, by
    Faa di Bruno,
        d u_k / dt = sum_j 1/j! sum_{i_1+...+i_j = k, i_l >= 1} D^j X(u_0)[u_{i_1},...,u_{i_j}],
    which is VE_1, ..., VE_r along the solution u_0(t).  So the algebras algebraic
    naturality quantifies over are exactly the ones that produce the variational
    equations, and closure says the method commutes with passing to them."""
    n, Ms, A = 2, 2, al.J4                                    # A = R[e]/(e^4), r = 3
    X = gen_map(n, Ms, 3)
    Mt = n * A.N
    Z = al.lift_map(A, n, X, Ms)
    emb = lambda q: al.subs(q, [P.var(Mt, A.N * j) for j in range(n)], Mt)   # q(u_0)
    blk = lambda k: [P.var(Mt, A.N * j + k) for j in range(n)]

    def term(idxs):                          # D^j X(u_0)[u_{i_1}, ..., u_{i_j}]
        out = []
        for i in range(n):
            acc = P(Mt)
            for multi in product(range(n), repeat=len(idxs)):
                q = X[i]
                for a in multi: q = q.diff(a)
                t = emb(q)
                for slot, a in zip(idxs, multi): t = t * blk(slot)[a]
                acc = acc + t
            out.append(acc)
        return out

    ok = True
    for k in range(A.N):
        ve = [emb(q) for q in X] if k == 0 else [P(Mt)] * n
        for j in range(1, k + 1):
            for c in compositions(k, j):
                t = term(c)
                ve = [ve[i] + t[i] * F(1, factorial(j)) for i in range(n)]
        good = [Z[A.N * jj + k] for jj in range(n)] == ve
        ok &= good
        print(f"  A0  block {k} of X^(R[e]/(e^4)) is "
              + ("X itself" if k == 0 else f"VE_{k}") + f": {good}"
              + ("" if good else "  <-- UNEXPECTED"))
    print("      -> the Weil algebras quantified over are the higher variational equations;")
    print("         closure says the method commutes with passing to them.")
    return ok


def main():
    return variational() & trichotomy() & rigidity() & cross()


if __name__ == "__main__": raise SystemExit(0 if main() else 1)
