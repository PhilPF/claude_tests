"""Shared exact-arithmetic machinery.

Forward-mode AD *is* evaluation over the dual numbers D = R[e]/(e^2), so running a
method over Dual gives exactly T(Psi^X); running it on the lifted field TX over R^{2n}
gives Psi^{TX}. Everything is over Fraction, so agreement is exact, not numerical.
"""
from fractions import Fraction as F


class Dual:
    """a + b*eps,  eps^2 = 0."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0): self.a, self.b = a, b
    def __add__(s, o): o = D(o); return Dual(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o): o = D(o); return Dual(s.a - o.a, s.b - o.b)
    def __rsub__(s, o): o = D(o); return Dual(o.a - s.a, o.b - s.b)
    def __neg__(s): return Dual(-s.a, -s.b)
    def __mul__(s, o): o = D(o); return Dual(s.a * o.a, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    def __repr__(s): return f"({s.a}+{s.b}e)"


def D(o): return o if isinstance(o, Dual) else Dual(o, 0)


class Jet:
    """a0 + a1*e + a2*e^2,  e^3 = 0.  Entries may live in any commutative ring."""
    __slots__ = ('c',)

    def __init__(s, *c): s.c = list(c)
    def __add__(s, o): o = J(o); return Jet(*[s.c[i] + o.c[i] for i in range(3)])
    __radd__ = __add__
    def __sub__(s, o): o = J(o); return Jet(*[s.c[i] - o.c[i] for i in range(3)])
    def __rsub__(s, o): o = J(o); return Jet(*[o.c[i] - s.c[i] for i in range(3)])
    def __neg__(s): return Jet(*[-a for a in s.c])
    def __mul__(s, o):
        o = J(o); a, b = s.c, o.c
        return Jet(a[0] * b[0], a[0] * b[1] + a[1] * b[0],
                   a[0] * b[2] + a[1] * b[1] + a[2] * b[0])
    __rmul__ = __mul__
    def __repr__(s): return f"({s.c[0]}+{s.c[1]}e+{s.c[2]}e2)"


def J(o): return o if isinstance(o, Jet) else Jet(o, 0, 0)


# ---- quadratic test field on R^2, valid over any commutative ring ----
def X(u):
    x1, x2 = u
    return [x1 * x2 + x2 * x2, x1 * x1 - 3 * x2]


def JX(u):
    x1, x2 = u
    return [[x2, x1 + 2 * x2], [2 * x1, F(-3)]]


def TX(u):
    """tangent lift of X, written out explicitly on R^4 (no duals used)."""
    x1, x2, v1, v2 = u
    return [x1 * x2 + x2 * x2,
            x1 * x1 - 3 * x2,
            x2 * v1 + (x1 + 2 * x2) * v2,
            2 * x1 * v1 - 3 * v2]


def JTX(u):
    x1, x2, v1, v2 = u
    return [[x2, x1 + 2 * x2, F(0), F(0)],
            [2 * x1, F(-3), F(0), F(0)],
            [v2, v1 + 2 * v2, x2, x1 + 2 * x2],
            [2 * v1, F(0), 2 * x1, F(-3)]]


# ---- cubic test field, so that Delta X is non-constant ----
def Xc(u):
    x1, x2 = u
    return [x1 * x1 * x1 + x1 * x2 * x2, x2 * x2 * x2 - x1 * x1 * x2]


def TXc(u):
    x1, x2, v1, v2 = u
    a = 3 * x1 * x1 + x2 * x2; b = 2 * x1 * x2
    c = -2 * x1 * x2; d = 3 * x2 * x2 - x1 * x1
    return [x1 * x1 * x1 + x1 * x2 * x2, x2 * x2 * x2 - x1 * x1 * x2,
            a * v1 + b * v2, c * v1 + d * v2]


def div(Jac): return sum((Jac[i][i] for i in range(len(Jac))), 0)


def matvec(M, w):
    return [sum((M[i][j] * w[j] for j in range(len(w))), 0) for i in range(len(M))]


def lap(f, u, n):
    """Exact Laplacian of f at u in R^n, via second-order jets:
       f(u + e*e_j) = f(u) + e*d_j f + (e^2/2) d^2_j f."""
    tot = None
    for j in range(n):
        out = f([Jet(u[k], 1 if k == j else 0, 0) for k in range(n)])
        col = [2 * o.c[2] for o in out]
        tot = col if tot is None else [tot[i] + col[i] for i in range(len(col))]
    return tot


def rk(f, u, h, A, b):
    """Explicit Runge-Kutta step."""
    s, k = len(b), []
    for i in range(s):
        arg = [u[d] + h * sum((A[i][j] * k[j][d] for j in range(i)), 0)
               for d in range(len(u))]
        k.append(f(arg))
    return [u[d] + h * sum((b[i] * k[i][d] for i in range(s)), 0) for d in range(len(u))]


RK4 = ([[0, 0, 0, 0], [F(1, 2), 0, 0, 0], [0, F(1, 2), 0, 0], [0, 0, 1, 0]],
       [F(1, 6), F(1, 3), F(1, 3), F(1, 6)])
HEUN = ([[0, 0], [1, 0]], [F(1, 2), F(1, 2)])

# standard test data
PT, TG, H = [F(2, 3), F(-5, 4)], [F(3, 5), F(7, 2)], F(1, 7)


def closure(method, f2, f4, pt=PT, tg=TG, h=H, n=2):
    """Return (T(Psi^f), Psi^{Tf}) so the caller can compare them.

    method(field, u, h, n) -> list.  f2 lives on R^n, f4 = its tangent lift on R^{2n}.
    """
    o = method(f2, [Dual(pt[i], tg[i]) for i in range(n)], D(h), n)
    lhs = ([x.a for x in o], [x.b for x in o])
    r = method(f4, list(pt) + list(tg), h, 2 * n)
    return lhs, (r[:n], r[n:])


def report(name, lhs, rhs, expect=True):
    ok = lhs == rhs
    mark = "OK  " if ok else "FAIL"
    flag = "" if ok == expect else "   <-- UNEXPECTED"
    print(f"  {name:46s} {mark}{flag}")
    if not ok and expect:
        print(f"      T(Psi^X) = {lhs}\n      Psi^(TX) = {rhs}")
    return ok == expect


class Trunc:
    """a_0 + a_1 e + ... + a_{r-1} e^{r-1},  e^r = 0.  Generalises Jet to any order."""
    __slots__ = ('r', 'c')

    def __init__(s, r, *c):
        s.r = r; s.c = list(c) + [0] * (r - len(c))
    def _lift(s, o): return o if isinstance(o, Trunc) else Trunc(s.r, o)
    def __add__(s, o):
        o = s._lift(o); return Trunc(s.r, *[s.c[i] + o.c[i] for i in range(s.r)])
    __radd__ = __add__
    def __sub__(s, o):
        o = s._lift(o); return Trunc(s.r, *[s.c[i] - o.c[i] for i in range(s.r)])
    def __rsub__(s, o):
        o = s._lift(o); return Trunc(s.r, *[o.c[i] - s.c[i] for i in range(s.r)])
    def __neg__(s): return Trunc(s.r, *[-a for a in s.c])
    def __mul__(s, o):
        o = s._lift(o); out = [0] * s.r
        for i in range(s.r):
            for j in range(s.r - i):
                out[i + j] = out[i + j] + s.c[i] * o.c[j]
        return Trunc(s.r, *out)
    __rmul__ = __mul__
    def __eq__(s, o): o = s._lift(o); return s.c == o.c
    def __repr__(s): return "+".join(f"{a}e{i}" for i, a in enumerate(s.c))


def packA(u, r, n=2):
    """R^{rn} -> A^n with A = R[e]/(e^r); level-major layout u[k*n + j] = x^(k)_j."""
    return [Trunc(r, *[u[k * n + j] for k in range(r)]) for j in range(n)]


def unpackA(w, r, n=2):
    return [w[j].c[k] for k in range(r) for j in range(n)]
