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


def tensor(A, B, name=None):
    """A (x) B in the tensor basis, ordered so that the realization of (A(x)B)^n is
    the realization of A^{(B^n)} -- i.e. e_alpha (x) f_beta sits at beta*N_A + alpha.
    With that ordering T^{A(x)B} = T^A o T^B holds on the nose, not just up to a
    shuffle; test_spectrum.py checks it."""
    NA, NB = A.N, B.N
    N = NA * NB
    idx = lambda al, be: be * NA + al
    table = {}
    for a1 in range(NA):
        for b1 in range(NB):
            for a2 in range(NA):
                for b2 in range(NB):
                    i, j = idx(a1, b1), idx(a2, b2)
                    if i > j: continue
                    v = [F(0)] * N
                    for a3 in range(NA):
                        if not A.c[a1][a2][a3]: continue
                        for b3 in range(NB):
                            if not B.c[b1][b2][b3]: continue
                            v[idx(a3, b3)] += A.c[a1][a2][a3] * B.c[b1][b2][b3]
                    table[(i, j)] = v
    one = [F(0)] * N
    for a in range(NA):
        for b in range(NB):
            one[idx(a, b)] = A.one[a] * B.one[b]
    return Alg(name or f"({A.name}) (x) ({B.name})", N, table, one)


# ---- the three invariants of section 4 -------------------------------------
# A contraction survives base change to (A, iota) according to these, and to
# nothing else: a loop scales by dim_R A, a liana by c_A, a stolon needs the
# declared-orthonormal metric to be A-balanced.

def c_A(A, g=None):
    """c_A = mu(g^{-1}) = sum_{alpha,beta} (g^{-1})^{alpha beta} e_alpha e_beta, the
    liana factor of section 4.  g is the inner product for which the chosen basis is
    declared orthonormal (default: the identity, i.e. the basis itself).  Lianas
    survive base change iff c_A = 1_A -- note 1_A, not the coordinate vector e_0:
    for a product algebra in its standard basis the unit is (1,...,1)."""
    gi = _inverse(g) if g else None
    out = [F(0)] * A.N
    for al in range(A.N):
        for be in range(A.N):
            w = gi[al][be] if gi else (F(1) if al == be else F(0))
            if w:
                for gm in range(A.N): out[gm] += w * A.c[al][be][gm]
    return out


def balanced(A, g=None):
    """Is multiplication by every a in A self-adjoint for g?  With g = I this reads
    c^beta_{gamma alpha} = c^alpha_{gamma beta}.  Stolons survive iff it holds."""
    N = A.N
    for gm in range(N):
        L = [[A.c[gm][al][be] for al in range(N)] for be in range(N)]   # (L_gm)_{be,al}
        for i in range(N):
            for j in range(N):
                lhs = sum((g[i][k] if g else (F(1) if i == k else F(0))) * L[k][j] for k in range(N))
                rhs = sum(L[k][i] * (g[k][j] if g else (F(1) if k == j else F(0))) for k in range(N))
                if lhs != rhs: return False
    return True


def _inverse(g):
    n = len(g)
    M = [list(g[i]) + [F(1) if i == j else F(0) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c])
        M[c], M[p] = M[p], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]; M[i] = [M[i][j] - f * M[c][j] for j in range(2 * n)]
    return [r[n:] for r in M]


def spectrum_invariants(A, g=None):
    """(dim, c_A, c_A == 1_A, balanced) -- loops scale by dim, lianas by c_A,
    stolons need balancedness.  Everything depends on (A, g) only, not on the
    basis: two bases with the same induced inner product give the same triple."""
    c = c_A(A, g)
    return A.N, c, c == A.one, balanced(A, g)


def derivations(A):
    """Basis of Der(A) = { d : d(ab) = d(a)b + a d(b) }, flattened N x N matrices with
    d(e_alpha) = sum_gamma d[gamma][alpha] e_gamma.  Der(A) is the Lie algebra of
    Aut(A); every derivation kills 1, since d(1) = d(1.1) = 2d(1)."""
    from test_affine import kernel
    N = A.N; rows = []
    for a in range(N):
        for b in range(N):
            for de in range(N):
                r = [F(0)] * (N * N)
                for g in range(N):
                    if A.c[a][b][g]: r[de * N + g] += A.c[a][b][g]
                for g in range(N):
                    if A.c[g][b][de]: r[g * N + a] -= A.c[g][b][de]
                    if A.c[a][g][de]: r[g * N + b] -= A.c[a][g][de]
                if any(r): rows.append(r)
    return kernel(rows, N * N) if rows else []


# ---- the based algebras used in the tests -------------------------------------

def Wk(k):
    """the jet algebra J^1_k = R[x_1..x_k]/m^2, basis (1, x_1, ..., x_k).  Its
    automorphism group is GL(k) acting on m = span(x_i) -- the largest Aut in its
    dimension, and the reason the Aut-generated algebra is all of gl(m-1)."""
    N = k + 1
    return Alg(f"R[x_1..x_{k}]/m^2", N,
               {(a, b): [F(0)] * N for a in range(1, N) for b in range(a, N)},
               [1] + [0] * k)


def Jr(r):
    """R[e]/(e^{r+1}), basis (1, e, ..., e^r)."""
    N = r + 1
    tab = {}
    for a in range(1, N):
        for b in range(a, N):
            v = [F(0)] * N
            if a + b < N: v[a + b] = F(1)
            tab[(a, b)] = v
    return Alg(f"R[e]/(e^{N})", N, tab, [1] + [0] * r)


R1   = Alg("R", 1, {}, [1])
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

ALGEBRAS = {1: [R1], 2: [D2, R2s, R2u], 3: [J3, W2, DxR, R3u, R3s],
            4: [J4, DD, W3, Axy, DxD, J3xR, W2xR, R4u]}


# ---- the semisimple quotient: how many real points A has ---------------------
# For a SMOOTH (non-polynomial) f, base change f^A(u) = sum_a d^a f(pi u) (u-pi u)^a/a!
# terminates only because u - pi(u) is nilpotent.  So the expansion point is not a
# convention: it is forced to be the real part, and there is one real part PER LOCAL
# FACTOR.  Nil(A) is the radical of the trace form (char 0, commutative), and
# dim_R A/Nil(A) counts those points.

def _trace_form(A):
    """B(x,y) = tr(L_{xy}) in the chosen basis."""
    t = [sum(A.c[a][b][b] for b in range(A.N)) for a in range(A.N)]      # tr(L_{e_a})
    return [[sum(A.c[a][b][g] * t[g] for g in range(A.N)) for b in range(A.N)]
            for a in range(A.N)]


def nilradical(A):
    """basis of Nil(A) = rad(trace form).  Each returned vector is checked nilpotent."""
    M = [row[:] for row in _trace_form(A)]
    N, piv, r = A.N, [], 0
    for c in range(N):
        p = next((i for i in range(r, N) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        M[r] = [x / M[r][c] for x in M[r]]
        for i in range(N):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [M[i][j] - f * M[r][j] for j in range(N)]
        piv.append(c); r += 1
    free = [c for c in range(N) if c not in piv]
    out = []
    for c in free:
        v = [F(0)] * N; v[c] = F(1)
        for i, p in enumerate(piv): v[p] = -M[i][c]
        out.append(v)
    return out


def is_nilpotent(A, x):
    y = list(x)
    for _ in range(A.N):
        if all(v == 0 for v in y): return True
        y = [sum(y[a] * x[b] * A.c[a][b][g] for a in range(A.N) for b in range(A.N))
             for g in range(A.N)]
    return all(v == 0 for v in y)


def num_real_points(A):
    """dim_R A/Nil(A).  For a product of Weil algebras this is the number of local
    factors, i.e. the number of expansion points a smooth base change needs."""
    return A.N - len(nilradical(A))
