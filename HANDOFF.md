# Handoff — read this first

State of play for continuing in a fresh session. Kept current; if it disagrees with a
conversation, trust this file and `flow-closure-under-differentiation.md`.

## The question

For a vector field $X$ with flow $\varphi^X_t$, "the derivative of the flow is the flow of the
derivative": $T(\varphi^X_t)=\varphi^{TX}_t$. Some integrators share it. **What does the property
force on an approximation of the flow?**

## The answer, in one paragraph

$T$ is the Weil functor of the dual numbers $D=\mathbb R[\varepsilon]/(\varepsilon^2)$, and $TX$ is not a new
vector field but $X$ **base-changed** $\mathbb R\to D$. So closure says: *the method commutes with base
change*. Anything that is the unique solution of an equation natural in the algebra base is closed —
which is why the exact flow (Picard) and Runge–Kutta (stage equations) are closed **for the same
reason**, nothing to do with the Butcher tableau. That criterion is sufficient, **not necessary**.

## What closure does NOT do — the framing that matters

Closure is a **per-tower** condition, not a per-dimension one. It relates $\Psi_n$ to $\Psi_{2n}$ on
the thin locus of lifted fields and, by the extension theorem, imposes **no condition whatever on any
individual $\Psi_d$**: every local smooth $\Psi_q$ ($q$ odd) extends to a closed sequence, and the
closed class has the cardinality of the continuum. Consequently:

* **"Closure $\Rightarrow$ equivariance" is false.** Partitioned RK and the Laplacian method are
  closed and not affine equivariant; there is a *polynomial* closed method with no equivariance.
* Every published-style characterisation here factors as
  $[\text{equivariance}]\cap[\text{closure}]=[\text{pruned series}]$. Equivariance + locality is what
  produces a *series* at all; **closure only prunes it**, killing exactly the contraction schemes with
  $N_0\neq1$. Do not present a relative result as if closure alone gave it.

## Current target: closure under ALL Weil algebras

Strengthen from one algebra to all of them, **coherently across factorizations** $m=n\dim A$. Call it
**algebraic naturality**: for each such $A$, $\Psi_{n\dim A}$ restricted to $A$-lifted fields is the
$\mathbb R$-realization of the $A$-base change of $\Psi_n$.

**Corrected: the class is finite-dimensional commutative $\mathbb R$-algebras, not Weil algebras.**
Kolář–Michor–Slovák's Weil algebras are exactly the product-preserving functors on $\mathbf{Mf}$, and
this file used to conclude "$\mathbb C$ is not one of them" — true, and the wrong conclusion.
Every residue field of a finite-dimensional commutative $\mathbb R$-algebra is $\mathbb R$ or
$\mathbb C$; the Weil products are exactly those with all residue fields $\mathbb R$ (equivalently
the trace form on $A/\mathrm{Nil}$ has no negative square — computable, `al.is_split`), and that is
exactly the condition for base change of **smooth** data, since $(*)$ terminates only on nilpotents.
For polynomial or (locally) analytic data any finite-dimensional commutative algebra works, and the
closure property is well posed there: **Euler, Heun, midpoint and RK4 are $T^{\mathbb C}$-natural**,
and §4's trichotomy predicts every non-Weil row (`verify/test_nonweil.py`). $\mathbb C$ is strictly
stronger, not a formality: $c_A$ is a **unit in every basis** over a split algebra, while
$c_{\mathbb C}=1+i^2=0$, so the **Laplacian method dies over $\mathbb C$ and survives over $D$** — a
separation no Weil algebra can produce. Complexification is a functor on real-analytic manifolds and
not on smooth ones; that, and only that, is what the KMS restriction records. Since jet transport is
run on Taylor methods with analytic fields, the larger class is the one the application lives in.

Already established (verified, **no equivariance assumed**):

* **Linear rigidity.** For any $A$ and any basis a linear field realises as $I_N\otimes M$, so the
  $D$-lift and the $\mathbb R^2$-lift coincide *exactly* on linear fields; demanding both forces
  $D\Psi^M(x)v=\Psi^M(v)$ — the method must be **linear on linear fields**. Euler ✓, RK4 ✓, the
  coordinate-1 method ✗.
* **Affine rigidity, and the ceiling above it.** The same comparison in the **unit-first** basis
  $(1,p)$ of $\mathbb R^2$ realises the $\mathbb R^2$-lift as a *secant*
  $\bigl(X(u^0),X(u^0{+}u^1)-X(u^0)\bigr)$ against the $D$-lift's *tangent*
  $\bigl(X(u^0),DX(u^0)u^1\bigr)$; secant $=$ tangent is affineness, so **affine fields go to affine
  maps**. Linear rigidity is the unit-misaligned shadow of this. And it is the **last** consequence
  of the mechanism: two distinct based algebras always differ at a product of two basis vectors, so
  a collision $X^A=X^B$ never reaches a non-affine field (collision trichotomy, §7). The gap over
  linear rigidity is real but lives only where translation equivariance fails.
* **Scalar rigidity.** A real-valued coefficient $q_m(j^rX)$ must be constant; base change demands an
  $A$-scalar and no real-valued contraction supplies one.
* **The closure spectrum.** $\mathcal W(\Psi)=\{(A,\iota):\Psi\text{ is }T^A\text{-natural}\}$ is a
  **monoid under $\otimes$** — one line, and structurally because $\dim$, $c_A$ and balancedness are
  all multiplicative — and is **not closed under subalgebras**: the Laplacian closes under
  $D\otimes D$ but not under its $S_2$-invariant $\mathbb R[s]/(s^3)$, even with the induced metric
  $\operatorname{diag}(1,2,1)$, where $c=1+xy$. The invariants depend on $(A,g)$ and not on $A$, so
  one algebra in two bases can have two spectra ($\mathbb R^2$ standard vs unit-first). For contraction
  methods the survival law is **coarser than §4's table**: a loop *or a stolon* survives only for
  $A=\mathbb R$, a liana iff $c_A=1$; balancedness is necessary but **not** sufficient, and the
  unifying reason is scalar rigidity, the defect being exactly $h^2(r\cdot1_A-q)\cdot X^A$. Realized
  spectra form a **lattice, not a chain**: the ring-valued tower realizes $\langle A\rangle_\otimes$
  for every based $A$, and $\langle D\rangle_\otimes$, $\langle\mathbb R^2\rangle_\otimes$ are
  incomparable. (§7, `verify/test_spectrum.py`.)
* **Presentations, via the commutant.** For $Z$ on $\mathbb R^m$, the commutant
  $\mathcal A(Z)=\{L:\ L\,D^kZ(u)[v_1,\dots]=D^kZ(u)[Lv_1,v_2,\dots]\ \forall k,u,v\}$ is a unital
  subalgebra, computable (the conditions are linear in $L$), and contains $\rho_A(A)$ whenever
  $Z=X^A$. Measured: $\mathcal A(X^A)=\rho_A(A)$ **exactly** for a generic field over all 17 based
  algebras, and in both coincidence families it is the algebra of the *finest* presentation and
  contains both coarser actions. So two presentations always refine to the commutant one, and the
  prescription on the union of lifted loci is **consistent** by induction on dimension. (§7,
  `verify/test_presentation.py`.)
* **Zero-field rigidity** *(corrected)*. $0^A=0$ for every $A$, so every pair of algebras of one
  dimension collides on the zero field, giving $\Psi^0_h(u)=s(h)u$ with $s\equiv1$ under
  consistency. But the zero field is **linear**, so the no-translation half is just linear rigidity
  at $M=0$; an earlier note claimed $u+hX+h^2e_1$ "escapes affine rigidity and dies to the zero-field
  collision", which is misleading — linear rigidity kills it outright. The new content is only
  $S_m=s(h)I_m$, so **no distinguished endomorphism** either.
* **What $\operatorname{Aut}$-equivariance imposes — much less than it sounds.** Every automorphism
  is unital, so $\phi_n$ fixes the real locus $1_A\otimes\mathbb R^n$ **pointwise**, and
  $\Psi_{n\dim A}$ restricted there is $\Psi_n$: the condition says **nothing about the base
  method**, only about the nilpotent directions. Collected over all factorizations of $m$ the
  generated Lie algebra is exactly $\mathfrak{gl}(m-1)$ (measured $9,25,49$ at $m=4,6,8$ against
  $\mathfrak{gl}(m)=16,36,64$), already realised by $J^1_{m-1}$ alone. **Correction:** the recorded
  claim that "the group generated is large enough ($\mathfrak{gl}(nN)$), the locus too thin" is
  false for the automorphism group — it is the *wrong group*, since it fixes the base pointwise, so
  it can never give $GL(n)$ there. Positively, for $A=J^1_k$ the lift is the $k$-fold Whitney sum
  and $\operatorname{Aut}=GL(k)$ mixes the copies, so the whole imposition is: no preferred basis of
  $\mathfrak m$. (§7, `verify/test_aut.py`.)
* **Only $\otimes$ transfers between algebras.** The spectrum is closed under $\otimes$ and under
  **nothing else**: not subalgebras (Laplacian; $\mathbb R[s]/(s^3)\subset D\otimes D$), not
  quotients ($\operatorname{tower}(D\otimes D)$; $D\otimes D\twoheadrightarrow D$), not products
  ($\operatorname{tower}(D)$; $D\times D$). $\otimes$ is the only relation that is one between the
  *functors*, so only it lets a condition be applied twice. Consequence: **there is no generating
  family of test algebras** — in particular one cannot verify naturality on the jet algebras
  $\mathbb R[x]/\mathfrak m^{r+1}$ and deduce it for their quotients, though every Weil algebra is
  one. This is a concrete reason (i)$\Rightarrow$(ii) resists.
* **Aut$(A)$-equivariance, and the first genuine $GL$.** For $\phi:A\to B$, $\phi_n$ intertwines the
  lifts, so $\phi_n\Psi^{X^A}=\Psi^{X^B}\phi_n$ on lifted fields; for $A=B$ this is
  $\operatorname{Aut}(A)$-equivariance of $\Psi_m$ on the $A$-lifted locus. Not vacuous:
  $\operatorname{Aut}(\mathbb R[x_1..x_k]/\mathfrak m^2)=GL(k)$, so that one algebra forces honest
  $GL(k)$-equivariance — in the **algebra** directions, not the $\mathbb R^n$ ones, which is exactly
  the recorded gap.
* **Correction — the non-equivariant witness.** $u+hX+h^2(X^1)^2X$ with the scalar left
  *real*-valued is **not** (T)-natural (defect $2X^1(DX^1v)X$, measured); earlier notes labelled it
  so. §5's own prescription — make the scalar *ring*-valued — gives the genuine witness: for
  $m=2^kq$ with $q$ odd, let $\Psi_m$ be the $D^{\otimes k}$-base change of that formula on
  $\mathbb R^q$. It is polynomial, local, non-equivariant, $T^{D^{\otimes k}}$-natural for all $k$, and
  still fails linear rigidity — so every conclusion drawn from the old witness survives.
* No non-trivial partitioned method is algebraically natural; dimension-dependent coefficients die.
* Separation of duties: $A=\mathbb R^k$ forces (diagonal) decoupling, killing multi-aromas; $A=D$ kills
  loops and stolons; $A=\mathbb R[\varepsilon]/(\varepsilon^3)$ kills lianas — **but that last is
  basis-dependent** (see the $c_A$ trap below), so state it relative to a fixed identification.

Two traps specific to this target:

* **Do not quantify over all bases $\iota$.** With $A=\mathbb R$, $\iota\in GL(n)$ is arbitrary and the
  condition *states* $GL(n)$-equivariance by fiat, trivialising the whole question. Fix one $\iota$
  per algebra.
* **Coherence does not recover $GL(n)$-equivariance**, and not for the obvious reason: the group
  generated is large enough ($\mathfrak{gl}(nN)$), but every base-change map **preserves the lifted
  locus**, so naturality never compares a lifted field against a generic one.

**The other Weil algebra — checked, and it is not stronger** (§7). Cartan's
$W(\mathfrak g)=\Lambda\mathfrak g^*\otimes S\mathfrak g^*$, the $EG$-model behind Chern–Weil and
equivariant cohomology, is a different object from Weil's algebras of infinitely near points; no
theorem identifies them. The common home is the graded world ($\mathbb R[\theta]$ odd gives
$\Pi T=T[1]$, $C^\infty(T[1]M)=\Omega(M)$). But it supplies **no new test algebra**: an ordinary
manifold's $A$-points are $(A_{\text{even}})^n$, so a super Weil algebra acts through its even part,
and $\mathrm{even}(\Lambda\mathbb R^3)$ is *exactly* $J^1_3$, already in the family. It would add
strength only by enlarging the category the **method** acts on (supermanifolds/graded bundles). The
topological notion's real point of contact is §8's connection-dependence, not a closure axiom.

**AD dictionary, and the shape under $T^*$** (§9). (T) is forward mode, (T*) is reverse mode — the
$\lambda$-part the method returns on the cotangent-lifted field is the discrete adjoint of its own
step. For RK, (T*) $\iff$ symplecticity (already recorded). **New and verified:** the shape
transports *differently* under the two lifts. Leapfrog is (T*)-closed on a nonlinear separable field,
but only with the **crossed** splitting $\{q,\lambda_p\}\mid\{p,\lambda_q\}$ — the obvious
cotangent-lifted one fails — because $\langle\lambda,X\rangle=\lambda_qf(p)+\lambda_pg(q)$ is
separable exactly there. $T$ carries a splitting blockwise; $T^*$ crosses it. So **leapfrog is
(T)-closed, (T*)-closed and not affine equivariant.** (T*) also carries a prerequisite (T) lacks —
$T^*$ is functorial only on diffeomorphisms, so $\Psi^X_h$ must be invertible; over $A$ that is free,
since a matrix over $A$ is invertible iff its real part is.

**Reframing (new, and it changes the target).** The property has an asymmetry the equivariance
framing hides: $\Psi$ receives $X^A$ as a *bare* field and is not told it is a lift, yet must return
$T^A$ of the unlifted answer. For a method that is an **algorithm over the base ring** — which is
what jet transport runs in practice — this is automatic, and the prerequisites are a short list
(§7): evaluation-only access to the field; $\mathbb R$-algebra operations only, with **division by
units alone**; implicit definitions when the equation is natural and uniquely solvable; and any
auxiliary structure must base-change functorially. A **splitting** does, a **metric** does not.
Consequently the property is stated relative to a **shape functor**, not to plain dimension —
leapfrog's shape is dimension-with-a-splitting, and of the three ways to match that splitting against
the $D$-lift exactly one works (verified). "PRK is natural only for partitioned indexing" is then an
instance, not a caveat. J1–J4 are sufficient and **not** necessary: programmable $\Rightarrow$
closed, and the gap is exactly the extension-theorem junk.

**Transport is not closure (new, and it retracts a definition given one round earlier).** A first
attempt defined jet transport as "a term in the $C^\infty$-theory generated by field evaluation and
the slot structure of a base-changing shape". **Retracted**: the Laplacian and aromatic methods use a
metric and a trace, so that definition denied them jet transport, yet both transport correctly over
all 17 based algebras (`verify/test_transport.py`, T6). Keep the two apart:

* **Transport** is an *operation on programs*: interpret the algorithm over $A$, the state ranging
  over $A^n$ with **$n$ unchanged**. If every step is a $C^\infty$-ring operation and the field is
  touched only by evaluation, then $P\langle A\rangle=T^A(P\langle\mathbb R\rangle)$ — always, by
  functoriality. A Weil algebra **is** a $C^\infty$-ring and base change **is** its structure, the
  finite Taylor sum $f^A(u)=\sum_\alpha\partial^\alpha f(\pi u)(u-\pi u)^\alpha/\alpha!$ (verified).
  So there is no interesting class of "methods that jet-transport": every $C^\infty$-typed program
  does, and that is a property of the *code*.
* **Closure** is the *coherence equation* between two ways of growing,
  $\Psi_{n\dim A}(\iota^*X^A)=\iota^*(\Psi_n\langle A\rangle)$ — raise the dimension over $\mathbb R$,
  versus raise the ring at fixed dimension. Transport supplies the right-hand side and makes the
  equation meaningful; closure is the equation, and it constrains the *family*.

That relocates the J-list: **J1–J3 are prerequisites for transport** (conditions on the code — field
access by evaluation; only $C^\infty$-operations, so branches are out and the real part $\pi$ is what
implementations branch on; division on units because $1/x$ is smooth exactly there; implicit
definitions by the implicit function theorem, all verified). **J4 is not a transport condition at
all** — a metric base-changes perfectly well — it is *uniformity of the family in the dimension
direction*, and §4's liana and loop laws are exactly its failure: verified, the Laplacian closes iff
$c_A=1_A$ and the aromatic method iff $\dim_{\mathbb R}A=1$, on all 17. The **shape** (a set with a
realization, $(\sigma^B)^A=\sigma^{A\otimes B}$, the cocycle
$\iota_{\sigma,A\otimes B}=T^A(\iota_{\sigma,B})\circ\iota_{\sigma^B,A}$, and an admissible class stable
under lift) is the bookkeeping of that uniformity, hence belongs to closure. Its payoffs: the
$\otimes$-monoid **is** the cocycle (verified over $D\otimes D$), $\iota$ is the standing
identification trap made a datum, and §9's crossed $T^*$ transport **is** stability of the admissible
class.

**The property, formally, and the NECESSARY conditions** (§7, `verify/test_necessary.py`). On
$\mathbf{FCA}$ = based finite-dimensional commutative $\mathbb R$-algebras $(A,\iota)$ with data class
$\mathcal D$: $\Psi$ is *algebraically natural* iff
$\Psi_{nN}[\iota_{n*}X^A]_h\circ\iota_n=\iota_n\circ(\Psi_n[X]_h)^A$ for all $n,(A,\iota),X,h$.
Basis, data class and the reality of $h$ are all **data**. Sufficient conditions do not characterise,
so here are six things a method cannot avoid:

* **P1 module rigidity.** $D\Psi^{X^A}_h(u)$ must lie in $\rho_A(A)'$ — the output must be
  $A$-differentiable; over $\mathbb C$, *holomorphic in, holomorphic out*. Intrinsic form via the
  commutant theorem: $\Psi$ must be equivariant for $\mathcal A(Z)$. Necessary and **strictly weaker**:
  the Laplacian satisfies it over $\mathbb C$ and still fails closure; the marked-coordinate method
  violates it.
* **P2 projections.** $\eta:\mathbb R\to A$ always exists (so the real locus does, and $\Psi$ restricts
  to $\Psi_n$ there); an augmentation $\pi:A\to\mathbb R$ exists iff $A$ has a real point, and
  $\#\operatorname{Hom}(A,\mathbb R)=p-q$, the **signature** of the trace form ($A/\mathrm{Nil}\cong
  \mathbb R^r\times\mathbb C^s$ gives signature $(r+s,s)$). $\mathbb C$ has none — so
  **(T-base)/(T-fibre) is a feature of split algebras, not of the property**. Note `num_real_points`
  returns $p+q=\dim A/\mathrm{Nil}$; use `real_points` for the honest count.
* **P3 the locus satisfies PDEs.** $A$-linearity of $DZ$ is Cauchy–Riemann over $\mathbb C$, hence
  harmonicity, hence the real Laplacian annihilates the $\mathbb C$-locus. Exact identity, verified on
  nine algebras: $\Delta_{\mathbb R^{nN}}(X^A)=\rho_A(c_A)\,(\Delta_AX)^A$. **§4's liana law is now a
  consequence, not a table**, and $c_{\mathbb C}=0$ *is* harmonicity.
* **P4 one family inside one dimension.** $A_\lambda=\mathbb R[x]/(x^2-\lambda)$ based by $(1,x)$:
  dimension 2 throughout, $\lambda>0$ gives $\mathbb R^2$, $\lambda=0$ gives $D$, $\lambda<0$ gives
  $\mathbb C$, $c_{A_\lambda}=1+\lambda$ takes **every** real value, and $\lambda=0$ is the unique member
  with a derivation ($T$ is the singular point). For $\Psi$ analytic in the jet the defect is analytic
  in $\lambda$; verified, the Laplacian's is a **degree-1 polynomial vanishing only at $\lambda=0$**.
  So **(T) alone is strictly weaker than closure over dimension 2**, at the same pair of dimensions —
  killing the Laplacian used to need dimension 3.
* **P5 the spectrum is analytically closed, hence never the Weil world.** $\mathcal W(\Psi)$ meets an
  irreducible algebraic family in an analytic subset; on a line that is everything or locally finite;
  $\{\lambda:A_\lambda$ split$\}=[0,\infty)$ is neither. **No method analytic in the jet has closure
  spectrum exactly the split algebras.** Sharpest form of "not about Weil algebras": it cannot be.
* **P6 the ceiling does not move.** Collisions across the family are affine and no more, so the union
  of lifted loci meets itself only in affine fields: the enlarged class **prunes harder without
  fattening the locus**. The gluing of §5 survives, affine rigidity is still the ceiling, and the
  missing axiom still has the shape of P1 extended off the lifted locus.

**Beyond the lifts: the $\operatorname{Der}(A)$ twists** (§7, `verify/test_nonweil.py`). The fields
$\pi_A$-related to $X$ are far more than the lifts: every $\delta\in\operatorname{Der}(A)$ gives a
vertical linear $\delta_A$ on $A^n$ with $X^A+\delta_A$ still $\pi_A$-related to $X$. Verified:
$[\delta_A,X^A]=0$ always (X^A is $\operatorname{Aut}(A)$-invariant), so the exact flows split as
$\varphi^{X^A}_t\circ T^{\exp(t\delta)}$; a method keeps the **base** half of that and does **not**
split the commuting pair. That is the right division of labour rather than a defect —
$\exp(t\delta)$ is transcendental ($e^{ct}$ on $D$'s fibre), so the algebra part must be taken
exactly, which is what $\operatorname{Aut}(A)$-equivariance supplies. So the invariant-manifold use of
jet transport (the twist is the internal dynamics; the invariance equation seeks a zero of
$X^A-\delta_A$) needs no new axiom. $\operatorname{Der}(\mathbb C)=0$: the twists are a nilpotent
phenomenon. More generally every property here is naturality on a **subcategory** of
$\mathbf{Vect}$ (objects $(M,X)$, morphisms $p$ with $Tp\circ Y=X\circ p$): all of it gives the flow
(§8), the affine morphisms give B-series, base change gives closure. Which subcategory is the dial.

The open questions, in order of interest:

0. **Which dimension-indexed structures are uniform, and must a closed family come from a
   dimension-indexed rule at all?** (The old form, "is every closed method a program?", is the wrong
   question: programs transport, and transport is not closure.)
   This replaces the old headline: leapfrog refutes "closure $\Rightarrow$ affine equivariance"
   outright, so that question was mis-posed. The reframed one splits in two — *which shapes
   base-change* (§4's trichotomy answers it for tensors) and *which programs over a fixed shape are
   closed* — and leapfrog is an instance of it rather than a counterexample.
1. Is "algebraically natural $\Rightarrow$ affine equivariant $\Rightarrow$ B-series" true for
   $f$-analytic families? The first implication is the open half; conjectured **false** for smooth
   methods (1-jet loci are nowhere dense, $\dim n^2N<M^2$), blocked on gluing junk up the tower.
2. Construct the flat, coherent, non-equivariant algebraically natural method, or prove none
   exists. **A candidate construction now exists at second order** (§7, `notes/order2-classification.md`,
   `verify/test_order2.py`), built on a *verified* detector: the **pointwise** commutant reads
   $\rho_A(A)$ exactly and — unlike $\dim V_A(u)$ — does **not** collapse on the real locus; it is
   multiplicative, separates the factorizations, and the recover-and-lift prescription needs no
   per-algebra data on polynomial fields. The proposal sums over divisors $N\mid m$ rather than over
   algebras, so the dimension-$\ge7$ moduli never enter. **It is a lead, not a result**: three gaps
   remain — mutually exclusive smooth cutoffs; commutativity/freeness of $\mathcal A_r$ (item 3's
   residual point, now load-bearing); and that the intrinsic formula is verified for *polynomial*
   $Y$ while the witness needs a compactly supported *smooth* bump, whose base change expands about
   the real point the formula claims not to use.
3. *(Closed modulo one point — §7.)* Cross-dimensional collisions are functorial because the
   commutant is the common refinement. What is left is to prove that $\mathcal A(Z)$ is always
   commutative with $\mathbb R^m$ free over it, so the finest presentation exists; verified in all
   19 cases here, not proved.
4. **Which $\otimes$-closed classes of based algebras are realized as closure spectra?** The tower
   construction realizes every *principal* one $\langle A\rangle_\otimes$; $\{c_A=1\}$ and
   everything are realized by the Laplacian and by RK. Is that the whole image? This is the
   surjectivity half of the Galois-shaped correspondence, and it is the natural next target.

## Also open

* The ambient $GL(n_1)\times GL(n_2)$ classification (aromatic P-series). Explicitly future work in
  Laurent–Munthe-Kaas §5; the closure half is already proved.

## Settled, do not re-derive

The write-up is the single source of truth. Headlines: closure $\iff$ base change; the obstruction is
**not** the trace but a trichotomy (loops die by $\operatorname{tr}_{A/\mathbb R}$, stolons by failure of
$A$-balancedness, **lianas survive** — and only over $D$, by the accident $c_D=1$); the multiplier is
$N_0$, counting families of vertex-disjoint generalised cycles; closure is **not** a jet condition
(flat counterexample); diffeomorphism-naturality collapses to $\varphi^X_{c(h)}$ with $c$ an arbitrary
reparametrisation; reverse mode costs symplecticity, for duality rather than functoriality reasons.
The per-dimension consequences of algebraic naturality stop at **affine $\Rightarrow$ affine**, by the
collision trichotomy (not to be confused with the contraction trichotomy above). The $O(n)$
characterisation is complete — Lemma L proved (`notes/lemma-L-proof.md`) and independently audited
(`notes/lemma-L-audit.md`).

## Provenance

Corroborated by two or more independent passes plus checks here: the trichotomy, $N_0$, the parity
law, closure $\not\Rightarrow$ B-series, the non-jet property. **Single source + one audit:** the
Lemma L proof. **Single source:** linear rigidity and scalar rigidity (both spot-checked here);
affine rigidity and the collision trichotomy (one pass, but the trichotomy is proved outright and
machine-checked on all 41 pairs from 16 based algebras — the *cross-dimensional* half is the part
resting on a bounded enumeration). The closure-spectrum results are also single-source, but
each is either a one-line proof or machine-measured on 17 based algebras.
Treat single-source items as leads when building on them.

## Keeping this current

Update this file whenever the target changes, a claim is retracted, or an open item closes — before
ending a session, not after. It is the only thing a fresh session reads to know where the work is.
