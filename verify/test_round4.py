"""Round-4 results, verified exactly over Q.

R1.  The base-block v-degree-0 multiplier of an exotic aromatic forest gamma is
     N_0(gamma) = # level assignments with K_nu = k_nu in {0,1} for every node
                = # families of vertex-disjoint GENERALISED CYCLES
     (arrows traversed forward, direction reversed at a liana or a stolon).
     N_0 is NOT 2^(any graph statistic): the theta-graph |D^2X|_F^2 X has N_0 = 3.
     On STOLON-FREE forests N_0 = 2^(#aromas), which is why the old law looked right.

R2.  The v-degree balance  v-deg = 2(s_1 - l_1) + [root level = 1]  is CORRECT,
     but the node rule printed in section B of the brief is not: the case giving
     v-degree +1 is  K_nu = 0, k_nu = 1  (not K_nu = 1, k_nu = 0), and the case
     K_nu = k_nu = 0 (contributing d^q X at v-degree 0) is missing.

R3.  The (T)-defect of a single elementary differential vanishes identically
     IFF gamma has no stolon and no directed arrow-cycle -- i.e. iff gamma is an
     exotic tree in the sense of Laurent-Munthe-Kaas (arXiv:2305.10993).

R4.  Coherence across factorizations (C3).  For ANY finite-dimensional commutative
     R-algebra A and ANY R-basis, the A-lift of a LINEAR field M is I_N (x) M.
     Hence the D-lifted and R^2-lifted loci meet in the linear fields and
     algebraic naturality forces D Psi^M(x) v = Psi^M(v), i.e. Psi^M LINEAR.
     This kills the round-3 counterexample u + hX + h^2 (X^1)^2 X.
"""
from fractions import Fraction as F
from poly import P
from forests import Forest, lift
from dual import Dual, D, RK4, rk, matvec

U = lambda v: ('U', v)
L = lambda v, s: ('L', v, s)

CASES = {
 "X":                         (Forest([0], [], U(0)), 1, True),
 "DX.X  (Butcher tree)":      (Forest([1,0], [(U(1),L(0,0))], U(0)), 1, True),
 "div(X) X  (1 aroma)":       (Forest([1,0], [(U(0),L(0,0))], U(1)), 2, False),
 "div(X)^2 X  (2 aromas)":    (Forest([1,1,0], [(U(0),L(0,0)),(U(1),L(1,0))], U(2)), 4, False),
 "Delta X  (liana)":          (Forest([2], [(L(0,0),L(0,1))], U(0)), 1, True),
 "d_jkl X^i X^l (liana+arrow)":(Forest([3,0], [(L(0,0),L(0,1)),(U(1),L(0,2))], U(0)), 1, True),
 "<X,X> X  (stolon)":         (Forest([0,0,0], [(U(0),U(1))], U(2)), 1, False),
 "|DX|_F^2 X (liana+stolon)": (Forest([1,1,0], [(L(0,0),L(1,0)),(U(0),U(1))], U(2)), 2, False),
 "|D2X|_F^2 X (theta graph)": (Forest([2,2,0], [(L(0,0),L(1,0)),(L(0,1),L(1,1)),(U(0),U(1))], U(2)), 3, False),
 "2 lianas + stolon at root":  (Forest([2,1,1], [(L(0,0),L(1,0)),(L(0,1),L(2,0)),(U(1),U(2))], U(0)), 1, False),
}

def Xcubic(n=2):
    x1, x2 = P.var(n,0), P.var(n,1)
    return [x1*x1*x1 + P.const(n,3)*x1*x2*x2 + x2*x2*x2*x2 + P.const(n,5)*x1*x1,
            x2*x2*x2*x1 - P.const(n,2)*x1*x1*x1 + P.const(n,7)*x1*x2 + x2*x2]

def main():
    ok = True
    n = 2
    Xc = Xcubic(n); TXc = lift(Xc, n)
    print("  R1/R2/R3  forest                       N0  mult=N0  balance  defect=0  predicted")
    for name,(g,N0exp,zero_exp) in CASES.items():
        N0 = g.N0()
        big = g.F(TXc, 2*n); ref = g.F(Xc, n)
        b0 = [c.subs_zero([n+j for j in range(n)]) for c in big[:n]]
        mult = all(b0[i] == ref[i].embed(2*n)*N0 for i in range(n))
        bal = (g.vdeg_spectrum(0)==g.balance_formula(0) and
               g.vdeg_spectrum(1)==g.balance_formula(1))
        # full defect: F_2n(gamma)(TX) - T(F_n(gamma)(X)) == 0 ?
        Tref = [r.embed(2*n) for r in ref]
        for i in range(n):
            s = P(2*n)
            for j in range(n): s = s + ref[i].embed(2*n).diff(j)*P.var(2*n, n+j)
            Tref.append(s)
        zero = all(big[i]==Tref[i] for i in range(2*n))
        st = g.stats()
        pred = (st.get('stolon',0)==0 and st.get('loop',0)==0)   # no stolon, no cycle
        good = (N0==N0exp) and mult and bal and (zero==zero_exp) and (zero==pred)
        ok &= good
        print(f"            {name:32s} {N0:>2d}   {str(mult):5s}   {str(bal):5s}"
              f"   {str(zero):5s}     {str(pred):5s}" + ("" if good else "  <-- UNEXPECTED"))
    print(f"    -> N_0 is not a power of 2 in general (theta graph: N_0 = 3):"
          f" {CASES['|D2X|_F^2 X (theta graph)'][0].N0()==3}")

    # ---- R4: coherence across factorizations --------------------------------
    print("\n  R4  algebraic naturality forces Psi^M to be LINEAR on linear fields M")
    M = [[F(1,2),F(-2)],[F(3),F(5,4)]]; H = F(1,7)
    x = [F(2,3),F(-5,4)]; v = [F(3,5),F(7,2)]
    def euler(f,u,h,n): fx=f(u); return [u[i]+h*fx[i] for i in range(n)]
    def rk4(f,u,h,n): return rk(f,u,h,*RK4)
    def m_ring(f,u,h,n):
        fx=f(u); s=fx[0]*fx[0]
        return [u[i]+h*fx[i]+h*h*s*fx[i] for i in range(n)]
    # the D-lift and the R^2-lift of a linear field agree (both are I_2 (x) M)
    lin = lambda u: matvec(M,u)
    Dlift = lambda u: matvec(M,u[:2])+matvec(M,u[2:])
    IM = [[(M[i%2][j%2] if i//2 == j//2 else F(0)) for j in range(4)] for i in range(4)]
    agree = all(Dlift(w) == matvec(IM, w)
                for w in ([F(1),F(2),F(3),F(4)], [F(-1,3),F(5,7),F(2),F(0)]))
    print(f"      M^D == M^(R^2) == I_2 (x) M : {agree}")
    ok &= agree
    for name,meth,exp in [("explicit Euler",euler,True),("RK4",rk4,True),
                          ("u+hX+h^2(X^1)^2 X (ring)",m_ring,False)]:
        o = meth(lin,[Dual(x[i],v[i]) for i in range(2)],D(H),2)
        good = ([p.b for p in o] == meth(lin,v,H,2))
        print(f"      {name:28s} D Psi^M(x)v == Psi^M(v): {good}"
              + ("" if good==exp else "  <-- UNEXPECTED"))
        ok &= (good==exp)
    print("      -> the round-3 (T)-natural non-equivariant method is NOT algebraically"
          " natural,\n         because a single Psi_4 cannot serve both (2,D) and (2,R^2).")
    return ok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
