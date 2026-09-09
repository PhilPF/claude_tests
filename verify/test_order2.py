"""Order-2 classification: what (*) B_{n dim A}(X^A) = (B_n X)^A forces on a local family.

Everything here is *pointwise* (jet-level), because locality makes B_m(Z)(u) a function of
j^r Z(u) alone.  The germ-level commutant of test_presentation.py is replaced by the pointwise
commutant A_r(u, j^r Z(u)); the two differ exactly where the loci degenerate, which is the
point at issue.
"""
from fractions import Fraction as F
from itertools import product as iproduct
from poly import P
import algebra as al
from algebra import ALGEBRAS
from test_affine import gen_map, kernel, _u
from test_presentation import rho, rank


# ---- pointwise jets and the pointwise commutant --------------------------------

def evalp(p, pt):
    """evaluate an exact polynomial at a rational point (extra slots read as 0)."""
    tot = F(0)
    for k, v in p.d.items():
        t = v
        for j, e in enumerate(k):
            if e:
                x = pt[j] if j < len(pt) else F(0)
                if x == 0: t = F(0); break
                t *= x ** e
        tot += t
    return tot


def pjet(Z, m, u, kmax):
    """{k: T_k} with T_k[i][(j_1..j_k)] = d^k Z^i / du_{j_1}..du_{j_k} at u."""
    out = {}
    cur = [[(tuple(), Z[i])] for i in range(m)]          # per component: (multi, poly)
    for k in range(1, kmax + 1):
        nxt, T = [], []
        for i in range(m):
            lst, d = [], {}
            for multi, p in cur[i]:
                for j in range(m):
                    q = p.diff(j)
                    if not q.is_zero():
                        lst.append((multi + (j,), q))
                        d[multi + (j,)] = evalp(q, u)
            nxt.append(lst); T.append(d)
        out[k] = T; cur = nxt
    return out


def pcommutant(T, m, kmax):
    """basis of {L : L T_k[v1..vk] = T_k[L v1, v2..vk], 1<=k<=kmax}, flattened m x m."""
    rows = []
    for k in range(1, kmax + 1):
        Tk = T.get(k)
        if Tk is None: continue
        for i in range(m):
            for multi in iproduct(range(m), repeat=k):
                row = [F(0)] * (m * m)
                for a in range(m):
                    row[i * m + a] += Tk[a].get(multi, F(0))
                rest = multi[1:]
                for b in range(m):
                    row[b * m + multi[0]] -= Tk[i].get((b,) + rest, F(0))
                if any(row): rows.append(row)
    if not rows:
        return [[F(1) if i == j else F(0) for i in range(m) for j in range(m)]]
    return kernel(rows, m * m)


def pcom_dim(Z, m, u, kmax=3):
    return len(pcommutant(pjet(Z, m, u, kmax), m, kmax))


# ---- test points ---------------------------------------------------------------

def pt(m, seed=0):
    """a rational point of R^m."""
    return [F((7 * (i + 1) + 3 * seed) % 11 - 5, (i % 3) + 2) for i in range(m)]


def real_pt(A, n, seed=0):
    """the honest real locus point 1_A (x) v of A^n = R^m.

    NOTE.  Zeroing the nilpotent coordinates (u[j*N+a] = 0 for a > 0) is the real locus
    only when 1_A = e_0.  For R^2 and R^3 in their *standard* bases the unit is
    (1,..,1), so that shortcut tests the wrong point; use this instead.
    """
    v = pt(n, seed)
    return [v[j] * A.one[a] for j in range(n) for a in range(A.N)]



# ---- jet vectors and lifted loci ------------------------------------------------

def multis(m, r):
    """all multi-indices (as sorted tuples) of length 0..r on m variables."""
    out = [()]
    cur = [()]
    for _ in range(r):
        nxt = []
        for t in cur:
            for j in range(t[-1] if t else 0, m):
                nxt.append(t + (j,))
        out += nxt; cur = nxt
    return out


def jetvec(Z, m, u, r, MS=None):
    """j^r Z(u) flattened: component i, multi-index alpha."""
    MS = MS or multis(m, r)
    out = []
    for i in range(m):
        for a in MS:
            p = Z[i]
            for j in a: p = p.diff(j)
            out.append(evalp(p, u))
    return out


def monfields(n, d):
    """spanning set of polynomial fields on R^n of degree <= d: e_i * x^alpha."""
    out = []
    for i in range(n):
        for a in multis(n, d):
            e = [0] * n
            for j in a: e[j] += 1
            X = [P(n) for _ in range(n)]
            X[i] = P(n, {tuple(e): F(1)})
            out.append((X, sum(e)))
    return out


def locus_rows(A, n, u, r, d):
    """rows spanning V_A(u) = {j^r (X^A)(u)}, plus the degree of each generator."""
    m = n * A.N
    MS = multis(m, r)
    rows, degs = [], []
    for X, dg in monfields(n, d):
        Z = al.lift_map(A, n, X, n)
        rows.append(jetvec(Z, m, u, r, MS)); degs.append(dg)
    return rows, degs


# ---- the intrinsic recover-and-lift formula ------------------------------------
# On a unit-first based algebra (1_A = e_0, m = span(e_1..e_{N-1})) the recovery of X
# from Z = X^A is the SAME linear map for every algebra of that dimension: read the
# e_0-coordinates.  And the re-lift can be written with no structure constants at all,
# using only the commutant subalgebra a = A_r(J) ⊆ End(R^m) and the generators
# g_j = e_{jN}.  So one formula per (n, N) covers every A of dimension N at once.

def eps_jet(Z, m, n, N, u, kmax):
    """recover (x, D^k X(x)) from j Z(u) by the e_0-coordinate map.  u-independent."""
    x = [u[j * N] for j in range(n)]
    T = pjet(Z, m, u, kmax)
    out = {0: [evalp(Z[i * N], u) for i in range(n)]}
    for k in range(1, kmax + 1):
        d = {}
        for multi in iproduct(range(n), repeat=k):
            d[multi] = [T[k][i * N].get(tuple(j * N for j in multi), F(0)) for i in range(n)]
        out[k] = d
    return x, out


def matmul(A, B, m):
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(m)] for i in range(m)]


def unflat(v, m):
    return [[v[i * m + j] for j in range(m)] for i in range(m)]


def intrinsic_lift(Kbasis, m, n, N, u, jetY, kmax):
    """Y^a(u) computed from the commutant a = span(Kbasis) alone: no structure
    constants, no basis of A.  jetY = {k: {multi: D^k Y(x)}} with x the e_0-part of u."""
    mats = [unflat(v, m) for v in Kbasis]
    t_ = len(mats)
    g = [[F(1) if i == j * N else F(0) for i in range(m)] for j in range(n)]   # generators
    x = [u[j * N] for j in range(n)]
    nu = [u[i] - sum(x[j] * g[j][i] for j in range(n)) for i in range(m)]
    nus = []                                   # nu_j in a, with nu_j . g_j = nu|block j
    for j in range(n):
        rows = [[mats[t][i][j * N] for t in range(t_)] + [nu[i]] for i in range(m)]
        c = solve(rows, t_)
        nus.append([[sum(c[t] * mats[t][a][b] for t in range(t_)) for b in range(m)]
                    for a in range(m)])
    out = [F(0)] * m
    for i in range(n):
        for a in range(m): out[a] += jetY[0][i] * g[i][a]
    fac = F(1)
    for k in range(1, kmax + 1):
        fac *= k
        for multi in iproduct(range(n), repeat=k):
            co = jetY[k].get(multi)
            if co is None or all(c == 0 for c in co): continue
            prod = [[F(1) if a == b else F(0) for b in range(m)] for a in range(m)]
            for j in multi: prod = matmul(prod, nus[j], m)
            for i in range(n):
                if co[i] == 0: continue
                for a in range(m): out[a] += co[i] * prod[a][i * N] / fac
    return out


def solve(rows, w):
    """least-structure exact solve of an overdetermined consistent system (w unknowns)."""
    R = [r[:] for r in rows]; piv = []; r = 0
    for c in range(w):
        p = next((i for i in range(r, len(R)) if R[i][c]), None)
        if p is None: continue
        R[r], R[p] = R[p], R[r]; pv = R[r][c]
        R[r] = [x / pv for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c]:
                f = R[i][c]; R[i] = [R[i][j] - f * R[r][j] for j in range(w + 1)]
        piv.append(c); r += 1
    sol = [F(0)] * w
    for i, c in enumerate(piv): sol[c] = R[i][w]
    return sol


def a_coords(Kbasis, m, n, N, u, gens=None):
    """write u in R^m = (+)_j a.g_j as (u_1..u_n), u_j in a, returned as matrices.
    gens defaults to g_j = e_{jN}, correct for every *unit-first* based algebra."""
    mats = [unflat(v, m) for v in Kbasis]; t_ = len(mats)
    out = []
    for j in range(n):
        g = gens[j] if gens else [F(1) if i == j * N else F(0) for i in range(m)]
        rows = [[sum(mats[t][i][b] * g[b] for b in range(m)) for t in range(t_)] + [u[i]]
                for i in range(m)]
        c = solve(rows, t_)
        out.append([[sum(c[t] * mats[t][a][b] for t in range(t_)) for b in range(m)]
                    for a in range(m)])
    return out


def intrinsic_subst(Kbasis, m, n, N, u, Y, gens=None):
    """(Y)^a(u) for a POLYNOMIAL field Y on R^n: substitute the a-point u into Y with
    products taken in a = span(Kbasis).  Uses only a ⊆ End(R^m) and g_j = e_{jN}:
    no basis of A, no structure constants, no augmentation, no real point."""
    U = a_coords(Kbasis, m, n, N, u, gens)
    G = gens or [[F(1) if i == j * N else F(0) for i in range(m)] for j in range(n)]
    I = [[F(1) if a == b else F(0) for b in range(m)] for a in range(m)]
    out = [F(0)] * m
    for i in range(n):
        acc = [[F(0)] * m for _ in range(m)]
        for k, v in Y[i].d.items():
            term = [row[:] for row in I]
            for j in range(n):
                for _ in range(k[j]): term = matmul(term, U[j], m)
            acc = [[acc[a][b] + v * term[a][b] for b in range(m)] for a in range(m)]
        for a in range(m): out[a] += sum(acc[a][b] * G[i][b] for b in range(m))
    return out


# ---- the checks -----------------------------------------------------------------

def gen_field(n, seed=0):
    """a generic cubic field on R^n for any n (gen_map of test_affine caps at n = 3)."""
    st = [seed * 37 + 11]
    def nxt():
        st[0] = (st[0] * 1103515245 + 12345) % 2147483647
        return F(st[0] % 17 - 8, st[0] % 5 + 2)
    out = []
    for i in range(n):
        p_ = P.const(n, nxt())
        for a in multis(n, 3):
            if not a: continue
            e = [0] * n
            for j in a: e[j] += 1
            p_ = p_ + P(n, {tuple(e): nxt()})
        out.append(p_)
    return out


def treeseed(X, n):
    """B_n(X) = D^2 X[X,X]: a tree elementary differential, zero on affine fields."""
    out = []
    for i in range(n):
        s = P(n)
        for a in range(n):
            for b in range(n): s = s + X[i].diff(a).diff(b) * X[a] * X[b]
        out.append(s)
    return out


def detector():
    """O1: the POINTWISE commutant of a lifted field is exactly rho_A(A) at every point --
    including the real locus, where dim V_A(u) collapses.  O2: room for junk."""
    ok = True
    print("  O1  pointwise commutant A_3(u, j^3 (X^A)(u)) of a lifted cubic field")
    print("      based algebra                    n  dimA   generic  real-locus  origin")
    for N, lst in sorted(ALGEBRAS.items()):
        if N < 2: continue
        for A in lst:
            for n in (1, 2):
                m = n * A.N
                if m > 6: continue
                Z = al.lift_map(A, n, gen_map(n, n, 3), n)
                cells = []
                for u in (pt(m, 1), real_pt(A, n, 1), [F(0)] * m):
                    K = pcommutant(pjet(Z, m, u, 3), m, 3)
                    good = (len(K) == A.N and rank(K + rho(A, n), m * m) == A.N)
                    ok &= good; cells.append(f"{len(K)}" + ("" if good else "!"))
                print(f"      {A.name:32s} {n:2d} {A.N:5d}   {cells[0]:>7s}  {cells[1]:>10s}  {cells[2]:>6s}")
    print("      -> = dim A in every cell (! would mark a failure of A_r = rho_A(A));")
    print("         the real-locus column is the one that matters: no collapse.")
    print("  O2  room for junk: dim A_3 of a generic cubic field, and of an affine field")
    for m in (2, 3, 4):
        g = min(pcom_dim(gen_field(m, s_), m, pt(m, 1)) for s_ in range(3))
        aff = [P(m, {tuple(1 if j == ((i + 1) % m) else 0 for j in range(m)): F(i + 2)})
               + P.const(m, F(i + 1)) for i in range(m)]
        a = pcom_dim(aff, m, pt(m, 1))
        good = (g == 1 and a >= m)
        ok &= good
        print(f"      m={m}:  generic cubic {g}   affine {a}  (>= m = {m})"
              + ("" if good else "   <-- UNEXPECTED"))
    print("      -> Omega_m = {dim A_r = 1} is open, nonempty, disjoint from every lifted")
    print("         locus, and contains no affine field: junk there is affine-rigid for free.")
    return ok


def multiplicative():
    """O3: the detector reads the FINEST presentation: dim A(X^A) = dim A * dim A(X)."""
    ok = True
    print("\n  O3  multiplicativity of the pointwise commutant under a second lift")
    print("      X = W^C          then lift by A       dim A(X)  dim A(X^A)  expected")
    for C in (ALGEBRAS[2][0], ALGEBRAS[2][1]):
        for A in (ALGEBRAS[2][0], ALGEBRAS[2][1], ALGEBRAS[3][0]):
            X = al.lift_map(C, 1, gen_map(1, 1, 3), 1)
            nX = C.N
            dX = pcom_dim(X, nX, pt(nX, 1), 3)
            m = nX * A.N
            dZ = pcom_dim(al.lift_map(A, nX, X, nX), m, pt(m, 1), 3)
            good = (dX == C.N and dZ == A.N * C.N); ok &= good
            print(f"      C = {C.name[:14]:16s} A = {A.name[:14]:16s} {dX:6d} {dZ:11d} {A.N*C.N:9d}"
                  + ("" if good else "   <-- UNEXPECTED"))
    return ok


def separates():
    """O4: the detector labels the factorization -- the rho_A(A) are pairwise distinct."""
    ok = True
    print("\n  O4  do the commutants separate the factorizations of m?")
    for m in (4, 6):
        facs = [(m // N, A) for N, lst in sorted(ALGEBRAS.items()) if N >= 2 and m % N == 0
                for A in lst]
        bad = []
        for i in range(len(facs)):
            for j in range(i + 1, len(facs)):
                n1, A1 = facs[i]; n2, A2 = facs[j]
                R1, R2 = rho(A1, n1), rho(A2, n2)
                if len(R1) == len(R2) and rank(R1 + R2, m * m) == len(R1):
                    bad.append((A1.name, A2.name))
        ok &= not bad
        print(f"      m = {m}: {len(facs)} factorizations (n,A), dim A >= 2;  coinciding "
              f"rho_A(A): {bad if bad else 'none'}")
    print("      -> distinct dimensions separate by dim A_r; equal dimensions by WHICH")
    print("         subalgebra the commutant is.  So the cutoffs chi are honest functions of J.")
    return ok


def retraction():
    """O5: the prescription (B_n X)^A(u) is computable from the commutant ALONE --
    one formula per (n, N), for every based algebra of that dimension at once."""
    ok = True
    print("\n  O5  intrinsic retraction: substitute the a-point u into B_n X, a = A_r(J)")
    print("      based algebra                    n  generic  real  origin   (unit-first?)")
    for N, lst in sorted(ALGEBRAS.items()):
        if N < 2: continue
        for A in lst:
            uf = (A.one == [F(1)] + [F(0)] * (A.N - 1))
            for n in (1, 2):
                m = n * A.N
                if m > 4: continue
                X = gen_map(n, n, 3); Z = al.lift_map(A, n, X, n)
                Y = treeseed(X, n); ZY = al.lift_map(A, n, Y, n)
                gens = None if uf else [[A.one[i % A.N] if i // A.N == j else F(0)
                                         for i in range(m)] for j in range(n)]
                cells = []
                for u in (pt(m, 1), real_pt(A, n, 1), [F(0)] * m):
                    K = pcommutant(pjet(Z, m, u, 3), m, 3)
                    got = intrinsic_subst(K, m, n, A.N, u, Y, gens)
                    good = (got == [evalp(ZY[i], u) for i in range(m)])
                    ok &= good; cells.append("ok" if good else "FAIL")
                print(f"      {A.name:32s} {n:2d} {cells[0]:>8s} {cells[1]:>5s} {cells[2]:>7s}"
                      f"   {'yes' if uf else 'no, generators supplied'}")
    print("      -> with the unit-first normal form the generators are g_j = e_{jN} for EVERY")
    print("         algebra, so the sum in Phi_m runs over divisors N | m, not over algebras:")
    print("         the moduli of commutative algebras in dim >= 7 never enter.")
    print("      Trap: expanding around the e_0-real part instead is correct for LOCAL algebras")
    print("         and wrong for every algebra with two or more local factors.")
    return ok


def main():
    return detector() & multiplicative() & separates() & retraction()


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
