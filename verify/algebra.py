"""Based finite-dimensional commutative R-algebras, base change, and lift recognition.

A *based* algebra is a pair (A, iota) with iota an R-linear isomorphism R^N -> A.
The basis is part of the data, not a convenience: what an algebraically natural
method sees is the realization A^n = R^{nN}, and that depends on iota.  Algebras
are recorded here by their structure constants in the chosen basis.

Conventions.  An element of A^n is a list of n algebra elements, each a list of N
entries; the realization orders coordinates blockwise, (j, alpha) -> j*N + alpha.
Entries are exact polynomials (poly.P) so every identity below is checked over Q.
A trailing variable is reserved for the step size h and is inert under base change
(h is a real scalar, not an A-point).
"""
from fractions import Fraction as F
from poly import P


class Alg:
    """table maps (a,b) with a<=b to the coordinate vector of e_a*e_b; one is 1_A."""

    def __init__(s, name, N, table, one):
        s.name, s.N, s.one = name, N, [F(x) for x in one]
        s.c = [[[F(0)] * N for _ in range(N)] for _ in range(N)]
        if s.one == [F(1)] + [F(0)] * (N - 1):        # unit-first basis: e_0 * e_b = e_b
            for b in range(N): s.c[0][b][b] = s.c[b][0][b] = F(1)
        for (a, b), vec in table.items():
            for g, v in enumerate(vec):
                s.c[a][b][g] = F(v)
                s.c[b][a][g] = F(v)

    def prod(s, a, b):
        return s.c[a][b]

    def check(s):
        """unital, commutative, associative -- guards against typos in the tables."""
        N = s.N
        ok = all(s.c[a][b] == s.c[b][a] for a in range(N) for b in range(N))
        for a in range(N):                                     # 1 * e_a = e_a
            v = [sum(s.one[b] * s.c[b][a][g] for b in range(N)) for g in range(N)]
            ok &= all(v[g] == (F(1) if g == a else F(0)) for g in range(N))
        for a in range(N):                                     # associativity
            for b in range(N):
                for d in range(N):
                    l = [sum(s.c[a][b][e] * s.c[e][d][g] for e in range(N)) for g in range(N)]
                    r = [sum(s.c[b][d][e] * s.c[a][e][g] for e in range(N)) for g in range(N)]
                    ok &= (l == r)
        return ok


# ---- arithmetic in A tensor (ring of exact polynomials) -----------------------

def amul(A, x, y, M):
    out = [P(M)] * A.N
    for a in range(A.N):
        if x[a].is_zero(): continue
        for b in range(A.N):
            if y[b].is_zero(): continue
            pr = x[a] * y[b]
            for g in range(A.N):
                cf = A.c[a][b][g]
                if cf: out[g] = out[g] + pr * cf
    return out


def aone(A, M):
    return [P.const(M, v) for v in A.one]


def apow(A, x, k, M):
    r = aone(A, M)
    for _ in range(k): r = amul(A, r, x, M)
    return r


def generic_point(A, n, M, shuffle=False):
    """the generic point of A^n in the realization R^{nN} (variables 0..nN-1)."""
    idx = (lambda j, al: al * n + j) if shuffle else (lambda j, al: j * A.N + al)
    return [[P.var(M, idx(j, al)) for al in range(A.N)] for j in range(n)]


def realize(A, n, elts, M, shuffle=False):
    idx = (lambda j, al: al * n + j) if shuffle else (lambda j, al: j * A.N + al)
    out = [None] * (n * A.N)
    for j in range(n):
        for al in range(A.N):
            out[idx(j, al)] = elts[j][al]
    return out


def basechange(A, n, polys, args, Ms, Mt, hs=None, ht=None):
    """polys: n polynomials in Ms variables (0..n-1 are the R^n coordinates, hs the
    step-size slot).  args: n algebra elements over the target ring (Mt variables).
    Returns the A-point F^A(args); base change acts on the R^n slots only."""
    out = [[P(Mt)] * A.N for _ in range(n)]
    for i, p in enumerate(polys):
        for k, v in p.d.items():
            elt = aone(A, Mt)
            for j in range(n):
                if k[j]: elt = amul(A, elt, apow(A, args[j], k[j], Mt), Mt)
            sc = P.const(Mt, v)
            if hs is not None and k[hs]:
                hp = P.var(Mt, ht)
                for _ in range(k[hs]): sc = sc * hp
            out[i] = [out[i][g] + elt[g] * sc for g in range(A.N)]
    return out


def lift_map(A, n, polys, Ms, hs=None, shuffle=False):
    """the realization of F^A as a polynomial map R^{nN} -> R^{nN}."""
    Mt = n * A.N + (0 if hs is None else 1)
    ht = None if hs is None else n * A.N
    args = generic_point(A, n, Mt, shuffle)
    return realize(A, n, basechange(A, n, polys, args, Ms, Mt, hs, ht), Mt, shuffle)


# ---- substitution and recognition of lifts ------------------------------------

def subs(p, images, Mt):
    out = P(Mt)
    for k, v in p.d.items():
        t = P.const(Mt, v)
        for j, e in enumerate(k):
            for _ in range(e): t = t * images[j]
        out = out + t
    return out


def remap(p, mapping, Mnew):
    d = {}
    for k, v in p.d.items():
        e = [0] * Mnew
        for j, a in enumerate(k):
            if a:
                if mapping[j] is None: return None      # variable that must not occur
                e[mapping[j]] += a
        d[tuple(e)] = d.get(tuple(e), F(0)) + v
    return P(Mnew, {k: v for k, v in d.items() if v})


def recognise(A, n, Z, hs=None, shuffle=False):
    """Is the polynomial map Z: R^{nN} -> R^{nN} the realization of X^A for some X?
    Returns X (n polynomials) or None.  The candidate X is read off the real locus
    1_A (x) R^n, which is forced, so the test is complete."""
    N, m = A.N, n * A.N
    Mt = m + (0 if hs is None else 1)
    ht = None if hs is None else m
    idx = (lambda j, al: al * n + j) if shuffle else (lambda j, al: j * A.N + al)
    a0 = next((a for a in range(N) if A.one[a]), None)
    if a0 is None: return None
    images = [P(Mt)] * Mt
    for j in range(n):
        for al in range(N):
            images[idx(j, al)] = P.const(Mt, A.one[al]) * P.var(Mt, idx(j, a0))
    if ht is not None: images[ht] = P.var(Mt, ht)
    Zr = [subs(z, images, Mt) for z in Z]
    for j in range(n):                                  # Z(1 (x) x) in 1 (x) R^n ?
        for al in range(N):
            if Zr[idx(j, al)] != Zr[idx(j, a0)] * A.one[al]: return None
    Ms = n + (0 if hs is None else 1)
    mp = [None] * Mt
    for j in range(n): mp[idx(j, a0)] = j
    if ht is not None: mp[ht] = hs
    X = []
    for j in range(n):
        q = remap(Zr[idx(j, a0)] * (1 / A.one[a0]), mp, Ms)
        if q is None: return None
        X.append(q)
    return X if lift_map(A, n, X, Ms, hs, shuffle) == Z else None


def lift_residual(A, n, Z, shuffle=False):
    """The obstruction to Z being an A-lift, as a list of polynomials depending
    LINEARLY on Z: Z is an A-lift iff every one of them vanishes.  (The candidate
    X is read off the real locus, which is a linear operation, and X -> X^A is
    linear, so the whole residual is.)  Used to size coincidence loci by rank."""
    N, m = A.N, n * A.N
    idx = (lambda j, al: al * n + j) if shuffle else (lambda j, al: j * A.N + al)
    a0 = next(a for a in range(N) if A.one[a])
    images = [P(m)] * m
    for j in range(n):
        for al in range(N):
            images[idx(j, al)] = P.const(m, A.one[al]) * P.var(m, idx(j, a0))
    Zr = [subs(z, images, m) for z in Z]
    res = [Zr[idx(j, al)] - Zr[idx(j, a0)] * A.one[al]
           for j in range(n) for al in range(N) if al != a0]
    mp = [None] * m
    for j in range(n): mp[idx(j, a0)] = j
    X = []
    for j in range(n):
        q = remap(Zr[idx(j, a0)] * (1 / A.one[a0]), mp, n)
        if q is None: return res + list(Z)          # real locus not even in 1 (x) R^n
        X.append(q)
    return res + [Z[t] - q for t, q in enumerate(lift_map(A, n, X, n, None, shuffle))]


# ---- the based algebras used in the tests -------------------------------------

D2   = Alg("D = R[e]/(e^2)", 2, {(1, 1): [0, 0]}, [1, 0])
R2s  = Alg("R^2 (standard basis)", 2, {(0, 0): [1, 0], (1, 1): [0, 1], (0, 1): [0, 0]}, [1, 1])
R2u  = Alg("R^2 (unit-first basis 1,p)", 2, {(1, 1): [0, 1]}, [1, 0])

J3   = Alg("R[e]/(e^3)", 3, {(1, 1): [0, 0, 1], (1, 2): [0, 0, 0], (2, 2): [0, 0, 0]}, [1, 0, 0])
W2   = Alg("R[x,y]/m^2", 3, {(1, 1): [0, 0, 0], (1, 2): [0, 0, 0], (2, 2): [0, 0, 0]}, [1, 0, 0])
DxR  = Alg("D x R (unit-first)", 3, {(1, 1): [0, 0, 0], (1, 2): [0, 0, 0], (2, 2): [0, 0, 1]}, [1, 0, 0])
R3u  = Alg("R^3 (unit-first)", 3, {(1, 1): [0, 1, 0], (1, 2): [0, 0, 0], (2, 2): [0, 0, 1]}, [1, 0, 0])
R3s  = Alg("R^3 (standard basis)", 3,
           {(0, 0): [1, 0, 0], (1, 1): [0, 1, 0], (2, 2): [0, 0, 1],
            (0, 1): [0, 0, 0], (0, 2): [0, 0, 0], (1, 2): [0, 0, 0]}, [1, 1, 1])

J4   = Alg("R[e]/(e^4)", 4, {(1, 1): [0, 0, 1, 0], (1, 2): [0, 0, 0, 1], (1, 3): [0] * 4,
                             (2, 2): [0] * 4, (2, 3): [0] * 4, (3, 3): [0] * 4}, [1, 0, 0, 0])
DD   = Alg("D (x) D = R[x,y]/(x^2,y^2)", 4,
           {(1, 1): [0] * 4, (1, 2): [0, 0, 0, 1], (1, 3): [0] * 4,
            (2, 2): [0] * 4, (2, 3): [0] * 4, (3, 3): [0] * 4}, [1, 0, 0, 0])
W3   = Alg("R[x,y,z]/m^2", 4, {(a, b): [0] * 4 for a in range(1, 4) for b in range(a, 4)},
           [1, 0, 0, 0])
Axy  = Alg("R[x,y]/(xy, x^2-y^2)", 4,
           {(1, 1): [0, 0, 0, 1], (1, 2): [0] * 4, (1, 3): [0] * 4,
            (2, 2): [0, 0, 0, 1], (2, 3): [0] * 4, (3, 3): [0] * 4}, [1, 0, 0, 0])
DxD  = Alg("D x D (unit-first 1,e,p,eta)", 4,
           {(1, 1): [0] * 4, (1, 2): [0] * 4, (1, 3): [0] * 4,
            (2, 2): [0, 0, 1, 0], (2, 3): [0, 0, 0, 1], (3, 3): [0] * 4}, [1, 0, 0, 0])
J3xR = Alg("R[e]/(e^3) x R (unit-first)", 4,
           {(1, 1): [0, 0, 1, 0], (1, 2): [0] * 4, (1, 3): [0] * 4,
            (2, 2): [0] * 4, (2, 3): [0] * 4, (3, 3): [0, 0, 0, 1]}, [1, 0, 0, 0])
W2xR = Alg("R[x,y]/m^2 x R (unit-first)", 4,
           {(1, 1): [0] * 4, (1, 2): [0] * 4, (1, 3): [0] * 4,
            (2, 2): [0] * 4, (2, 3): [0] * 4, (3, 3): [0, 0, 0, 1]}, [1, 0, 0, 0])
R4u  = Alg("R^4 (unit-first)", 4,
           {(1, 1): [0, 1, 0, 0], (2, 2): [0, 0, 1, 0], (3, 3): [0, 0, 0, 1],
            (1, 2): [0] * 4, (1, 3): [0] * 4, (2, 3): [0] * 4}, [1, 0, 0, 0])

ALGEBRAS = {2: [D2, R2s, R2u], 3: [J3, W2, DxR, R3u, R3s],
            4: [J4, DD, W3, Axy, DxD, J3xR, W2xR, R4u]}
