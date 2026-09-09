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

def pt(m, seed=0, real_of=None):
    """a rational point of R^m; real_of=(n,N) zeroes the nilpotent slots (real locus)."""
    vals = [F((7 * (i + 1) + 3 * seed) % 11 - 5, (i % 3) + 2) for i in range(m)]
    if real_of:
        n, N = real_of
        vals = [vals[j * N + a] if a == 0 else F(0) for j in range(n) for a in range(N)]
    return vals


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
