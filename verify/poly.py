"""Exact multivariate polynomials over Q, and exotic-aromatic elementary differentials."""
from fractions import Fraction as F
from itertools import product

class P:
    """dict: exponent tuple (len m) -> Fraction"""
    __slots__=('m','d')
    def __init__(s,m,d=None): s.m=m; s.d=dict(d or {})
    @staticmethod
    def const(m,c): return P(m,{(0,)*m:F(c)}) if c else P(m)
    @staticmethod
    def var(m,j):
        e=[0]*m; e[j]=1; return P(m,{tuple(e):F(1)})
    def __add__(s,o):
        d=dict(s.d)
        for k,v in o.d.items():
            d[k]=d.get(k,F(0))+v
            if d[k]==0: del d[k]
        return P(s.m,d)
    def __neg__(s): return P(s.m,{k:-v for k,v in s.d.items()})
    def __sub__(s,o): return s+(-o)
    def __mul__(s,o):
        if isinstance(o,(int,F)):
            return P(s.m,{k:v*o for k,v in s.d.items() if v*o!=0})
        d={}
        for k1,v1 in s.d.items():
            for k2,v2 in o.d.items():
                k=tuple(a+b for a,b in zip(k1,k2))
                d[k]=d.get(k,F(0))+v1*v2
        return P(s.m,{k:v for k,v in d.items() if v!=0})
    __rmul__=__mul__
    def diff(s,j):
        d={}
        for k,v in s.d.items():
            if k[j]:
                e=list(k); c=e[j]; e[j]-=1
                d[tuple(e)]=d.get(tuple(e),F(0))+v*c
        return P(s.m,{k:v for k,v in d.items() if v!=0})
    def subs_zero(s,idxs):
        """set variables in idxs to 0"""
        d={}
        for k,v in s.d.items():
            if all(k[j]==0 for j in idxs):
                d[k]=d.get(k,F(0))+v
        return P(s.m,{k:v for k,v in d.items() if v!=0})
    def embed(s,m2):
        """same poly, ambient dimension m2>=m (pad exponents with zeros)"""
        return P(m2,{k+(0,)*(m2-s.m):v for k,v in s.d.items()})
    def is_zero(s): return not s.d
    def __eq__(s,o): return s.d==o.d
    def __repr__(s):
        if not s.d: return "0"
        return "+".join(f"{v}*x^{k}" for k,v in sorted(s.d.items()))

def deriv(p, multi):
    for j in multi: p=p.diff(j)
    return p
