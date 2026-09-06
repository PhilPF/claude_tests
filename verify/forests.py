"""Exotic aromatic forests with ARROWS, LOOPS, LIANAS and STOLONS, their elementary
differentials, and the level-flow count N_0.

Half-edge model.  Nodes 0..p-1.  Node nu has q[nu] LOWER half-edges (derivative slots)
and exactly one UPPER half-edge (its output index).  Half-edges:
    ('U', nu)          upper of nu
    ('L', nu, s)       s-th lower of nu
A forest is a perfect matching of all half-edges except one distinguished free upper
half-edge (the root).  Pair types:
    U-L  arrow  (source nu -> target nu'; nu==nu' is a LOOP)
    U-U  stolon
    L-L  liana
"""
from fractions import Fraction as F
from itertools import product
from poly import P, deriv


class Forest:
    def __init__(self, q, pairs, root):
        self.q = list(q); self.p = len(q); self.pairs = [tuple(x) for x in pairs]
        self.root = root                      # ('U', nu)
        HE = [('U', nu) for nu in range(self.p)] + \
             [('L', nu, s) for nu in range(self.p) for s in range(self.q[nu])]
        seen = []
        for a, b in self.pairs: seen += [a, b]
        assert sorted(map(str, seen + [root])) == sorted(map(str, HE)), "not a perfect matching"

    def kind(self, e):
        (a, b) = e
        if a[0] == 'U' and b[0] == 'U': return 'stolon'
        if a[0] == 'L' and b[0] == 'L': return 'liana'
        return 'loop' if a[1] == b[1] else 'arrow'

    def stats(self):
        from collections import Counter
        return Counter(self.kind(e) for e in self.pairs)

    # ---------- elementary differential -------------------------------------
    def F(self, Xc, n):
        """Xc[i] : list of n polynomials in n vars.  Returns list of n polynomials."""
        edges = self.pairs
        m = len(edges)
        out = [P(n) for _ in range(n)]
        # map half-edge -> index variable number (edge id), root -> free
        loc = {}
        for k, (a, b) in enumerate(edges): loc[a] = k; loc[b] = k
        for i in range(n):
            loc_root = i
            for assign in product(range(n), repeat=m):
                idx = lambda he: (loc_root if he == self.root else assign[loc[he]])
                term = P.const(n, 1)
                for nu in range(self.p):
                    up = idx(('U', nu))
                    low = [idx(('L', nu, s)) for s in range(self.q[nu])]
                    term = term * deriv(Xc[up], low)
                    if term.is_zero(): break
                out[i] = out[i] + term
        return out

    # ---------- level-flow count N_0 ----------------------------------------
    def N0(self):
        """number of level assignments (subsets E1 of edges, root at level 0) with
        every node satisfying K_nu = k_nu in {0,1}."""
        cnt = 0
        m = len(self.pairs)
        for sub in product((0, 1), repeat=m):
            K = [0] * self.p; k = [0] * self.p
            ok = True
            for e, lv in zip(self.pairs, sub):
                if not lv: continue
                for he in e:
                    if he[0] == 'U': k[he[1]] += 1
                    else: K[he[1]] += 1
            for nu in range(self.p):
                if K[nu] != k[nu] or K[nu] > 1: ok = False; break
            if ok: cnt += 1
        return cnt

    def vdeg_spectrum(self, root_level):
        """multiset of v-degrees of surviving level assignments with given root level."""
        res = {}
        m = len(self.pairs)
        for sub in product((0, 1), repeat=m):
            K = [0] * self.p; k = [0] * self.p
            k[self.root[1]] += root_level
            for e, lv in zip(self.pairs, sub):
                if not lv: continue
                for he in e:
                    if he[0] == 'U': k[he[1]] += 1
                    else: K[he[1]] += 1
            d = 0; ok = True
            for nu in range(self.p):
                if K[nu] >= 2 or k[nu] >= 2: ok = False; break
                if K[nu] == 1 and k[nu] == 0: ok = False; break
                if K[nu] == 0 and k[nu] == 1: d += 1
            if ok: res[d] = res.get(d, 0) + 1
        return res

    def balance_formula(self, root_level):
        """predicted v-degrees 2(s1-l1)+[root=1] over surviving assignments"""
        res = {}
        m = len(self.pairs)
        for sub in product((0, 1), repeat=m):
            K = [0] * self.p; k = [0] * self.p
            k[self.root[1]] += root_level
            s1 = l1 = 0
            for e, lv in zip(self.pairs, sub):
                if not lv: continue
                kd = self.kind(e)
                if kd == 'stolon': s1 += 1
                if kd == 'liana': l1 += 1
                for he in e:
                    if he[0] == 'U': k[he[1]] += 1
                    else: K[he[1]] += 1
            ok = all(not (K[nu] >= 2 or k[nu] >= 2 or (K[nu] == 1 and k[nu] == 0))
                     for nu in range(self.p))
            if ok:
                d = 2 * (s1 - l1) + root_level
                res[d] = res.get(d, 0) + 1
        return res


# ---------- the tangent lift as an explicit field on R^{2n} -----------------
def lift(Xc, n):
    """TX on R^{2n}: components 0..n-1 = X(x); n..2n-1 = DX(x) v."""
    Xe = [c.embed(2 * n) for c in Xc]
    out = list(Xe)
    for i in range(n):
        s = P(2 * n)
        for j in range(n):
            s = s + Xe[i].diff(j) * P.var(2 * n, n + j)
        out.append(s)
    return out
