# Flow closure under differentiation

**What structure does "the derivative of the flow is the flow of the derivative" force on an
approximation of the flow?**

For a smooth vector field $X$ on a manifold $M$ with flow $\varphi^X_t$ it is classical that
$T(\varphi^X_t)=\varphi^{TX}_t$, where $T$ is the tangent functor and $TX$ the complete
(tangent) lift. Certain numerical integrators share the property. This note determines what
the property *forces*, rather than exhibiting families that happen to enjoy it.

Every computational claim below is verified in exact rational arithmetic; see
[`verify/`](verify/) and the [Verified computations](#12-verified-computations) index.

---

## 1. Formalisation

**The lift.** For $X\in\mathfrak X(M)$ the complete lift is $X^{C}=\kappa_M\circ T(X)\in\mathfrak X(TM)$,
with $\kappa_M:T^2M\to T^2M$ the canonical involution. In induced coordinates on
$T\mathbb R^n\cong\mathbb R^n\times\mathbb R^n$,
$$TX(x,v)=\bigl(X(x),\,DX(x)v\bigr),$$
whose second component is the variational equation. $TX$ is the unique field on $TM$ that is
$\pi_M$-related to $X$ and linear over $X$.

**Methods are sequences, not single maps.** The closure property relates the method on
$\mathbb R^n$ to the method on $\mathbb R^{2n}$, so it cannot be stated in fixed dimension. A
*method* is a sequence $\Psi=\{\Psi_n\}_{n\ge1}$, $\Psi_n^X{}_{,h}:\mathbb R^n\to\mathbb R^n$, one
per dimension — the setting of McLachlan–Modin–Munthe-Kaas–Verdier (MMMV). Standing
hypotheses: **locality** ($\Psi^X_h(x)$ depends only on the germ of $X$ at $x$; by Peetre's
theorem this forces locally finite jet order), smoothness, and **consistency**.

**The closure property.**
$$\textbf{(T)}\qquad T\bigl(\Psi^X_h\bigr)=\Psi^{TX}_h \qquad\text{for all }n,\;X,\;h\ \text{small.}$$
Equivalently: forward-mode automatic differentiation of the integrator equals the integrator
applied to the variational equation. It splits into two logically independent halves:

* **(T-base)** the base component of $\Psi^{TX}$ is $v$-independent and equals $\Psi^X$;
* **(T-fibre)** the fibre component equals $D\Psi^X(x)v$.

These are genuinely independent: $x+hX+h^2\langle e_1,X\rangle e_1$ satisfies (T-base) and
fails (T-fibre).

---

## 2. The reduction: (T) is base change

Let $D=\mathbb R[\varepsilon]/(\varepsilon^2)$ be the dual numbers. Then $T$ is the Weil functor of $D$:
$$T\mathbb R^n=D\otimes_\mathbb R\mathbb R^n,$$
a free $D$-module of rank $n$. Crucially, **$TX$ is not a new vector field**: it is $X$ with its
coefficients extended from $\mathbb R$ to $D$,
$$X(x_0+\varepsilon x_1)=X(x_0)+\varepsilon\,DX(x_0)x_1,$$
which is literally the $\varepsilon^2=0$ truncation of the Taylor series. Hence:

> **Reduction.** (T) holds iff $\Psi_{2n}$, *restricted to the locus of lifted fields*, agrees
> with the $D$-base-change of $\Psi_n$.

Two consequences that are easy to overstate and must not be:

1. (T) constrains only the **doubling tower** $n,2n,4n,\dots$, and only on the thin locus of
   lifted fields. It says nothing about $\Psi_{2n}$ elsewhere, and **nothing whatsoever
   relating dimension $n$ to $n+1$**.
2. The identification $\mathbb R^{2n}\cong D^n$ is part of the data. Which real $2n$-frame is used
   matters: a partitioned Runge–Kutta method satisfies (T) for the $D$-submodule block lift
   $(q,\dot q)\,|\,(p,\dot p)$ and **fails** for the real splittings base$|$fibre and the crossed
   one. So (T) is a *coherence condition on the tower*: the structure used in dimension $2n$
   must be the base change of the structure used in dimension $n$.

---

## 3. Sufficiency — and why the flow and the integrators agree for the same reason

> **Theorem (sufficiency).** Any operator defined as the **unique** solution of an equation
> built from evaluation of $X$ together with operations that are defined over an arbitrary
> commutative $\mathbb R$-algebra base, and natural in that base, satisfies (T) — indeed
> $T^A$-naturality for **every** Weil algebra $A$. The reason is that *uniqueness of a solution
> is preserved under base change*.

*Exact flow.* The Picard operator $\mathcal P(\gamma)(t)=x_0+\int_0^tX(\gamma(s))\,ds$ is built
from evaluation of $X$ and $\mathbb R$-linear operations. Since $\dim_\mathbb R D<\infty$ we have
$C([0,T],E\otimes D)=C([0,T],E)\otimes D$ with the Bochner integral commuting with base
change. The base-changed fixed point is a fixed point of the base-changed operator;
uniqueness forces equality.

*Runge–Kutta.* The stage equations $k_i=X\bigl(x+h\sum_ja_{ij}k_j\bigr)$, $\Psi^X_h=x+h\sum_ib_ik_i$
lie in the same language. Applying $T$ and using $T(u+w)=Tu+Tw$, $T(\lambda u)=\lambda Tu$ for
$\lambda\in\mathbb R$, the stages $K_i:=Tk_i$ satisfy the *same* equations for $TX$. Over $D$ the
stage system is **block-triangular** — solve the base component, then the $\varepsilon$-component
solves a linear system with the same invertible operator $I-h(a_{ij}\otimes DX)$ — so
solvability and uniqueness hold on the same $h$-range, implicit methods included.

This is the answer to *"what makes the integrators satisfy it?"*: **nothing about the Butcher
tableau.** The flow and the RK/PRK methods are (T)-natural for one and the same reason.

**Sufficiency is not necessary.** A tempting but *incorrect* rendering of the criterion is
"evaluation of $X$ plus $\mathbb R$-affine operations". That fails to separate the good from the
bad: $\operatorname{div}X$ is $\mathbb R$-linear in $X$ and a limit of affine combinations of evaluations,
exactly like $\Delta X$, yet $\operatorname{div}$ breaks (T) while $\Delta$ does not; and the Picard
argument already requires limits, so the language cannot exclude them. The correct statement
is the functor-of-points one above. Moreover it is **strictly stronger than (T)**: the
Laplacian method $x+hX+h^2\Delta X$ satisfies (T) but fails $T^A$-naturality for
$A=\mathbb R[\varepsilon]/(\varepsilon^3)$ (§7), hence is not definable in this sense. So

$$\text{definable}\;\Longrightarrow\;(T),\qquad\text{but}\quad (T)\;\not\Longrightarrow\;\text{definable}.$$

---

## 4. What (T) forbids: a trichotomy of contractions, not "the trace"

The natural first guess — *the trace is the unique obstruction* — is **false**. The correct
statement is:

> (T) annihilates exactly those invariant contraction schemes whose $\mathbb R$-realization on
> $\mathbb R^{2n}$ differs from their $D$-realization on $D^n$.

Write $A=\mathbb R[\varepsilon]/(\varepsilon^r)$, so $\mathbb R^{rn}=A^n$. On base-changed objects a derivative
index at level $k$ acts as $\varepsilon^k\partial_z$, and an evaluation index at level $k$ reads off the
$\varepsilon^k$-component. Per contracted pair:

| contraction | graph name | example | real realization on lifted objects | (T)? |
|---|---|---|---|---|
| one upper, one lower, **different** factors | **arrow** (tree edge) | $DX\,X$ | $\sum_\beta e^\beta(y)e_\beta=y$ — exact, metric-free, no defect | ✓ always |
| upper $\leftrightarrow$ lower, **same** factor | **loop / aroma** | $\operatorname{div}X$ | identity replaced by $\operatorname{tr}_{A/\mathbb R}(y)=N\cdot\operatorname{aug}(y)$ — real-valued: fibre killed, scaled by $N$ | ✗ |
| two **lower** indices, $\delta^{jk}$ | **liana** | $\Delta X$ | global factor $c_A$ | ✓ iff $c_A=1$ |
| two **upper** indices, $\delta_{jk}$ | **stolon** | $\langle X,X\rangle$ | $A$-multiplication $\mu$ replaced by the $\mathbb R$-valued form $g$ — not $A$-bilinear | ✗ |

Here $N=\dim_\mathbb R A$ and, for an $\mathbb R$-basis $\{e_\alpha\}$ of $A$ declared orthonormal (for $A=D$ the
basis $\{1,\varepsilon\}$, giving the Sasaki metric),
$$c_A:=\sum_\alpha e_\alpha^2=\mu(g^{-1})\in A .$$
This single constant governs the liana row and explains everything about it:
$$c_D=1^2+\varepsilon^2=1,\qquad c_{D^{\otimes m}}=1,\qquad c_{\mathbb R[\varepsilon]/(\varepsilon^{r+1})}=\sum_k\varepsilon^{2k}=1+\varepsilon^2+\cdots\ (r\ge2).$$
**$c_A$ is basis-dependent, and this matters.** Diagonal rescalings of the standard basis never
give $c_A=1$: with $g=\operatorname{diag}(1,w_1,\dots)$ one gets $c_A=1+w_1^{-1}\varepsilon^2+\cdots$. But a
non-diagonal basis can. For $A=\mathbb R[\varepsilon]/(\varepsilon^3)$ take
$$e_0=1-\tfrac{c^2}{2}\varepsilon^2,\qquad e_1=c\varepsilon+d\varepsilon^2,\qquad e_2=f\varepsilon^2\qquad(c,f\neq0),$$
a genuine basis (change-of-basis determinant $cf$), for which $c_A=\sum_\alpha e_\alpha^2=1$ **exactly**
— verified over $\mathbb Q$. The phenomenon is not confined to that algebra: it bites already for the
**dual numbers**, where the sheared basis $\{1+\varepsilon,\ \varepsilon\}$ (determinant $1$) gives
$c_D=(1+\varepsilon)^2+\varepsilon^2=1+2\varepsilon\neq1$. So even "lianas are free under $T$" is relative to the
*canonical* tangent identification $(x,v)\mapsto x+\varepsilon v$ — which is what makes the standard basis
the right one, but which must be stated rather than assumed. The liana obstruction is a property of
the pair (algebra, chosen orthonormal basis), not of the algebra alone, and every statement below about lianas failing over
$\mathbb R[\varepsilon]/(\varepsilon^3)$ is relative to the standard identification $\{1,\varepsilon,\varepsilon^2\}$. Any
basis-free formulation must quantify over all $\mathbb R$-bases of $A$. This is one instance of a
pattern worth stating plainly: **be suspicious of any claim here that silently fixes an
identification $\mathbb R^{n\dim A}\cong A^n$** — the partitioned-RK lift (§2) was the first instance,
this is the second.

*(Terminology is that of Laurent–Munthe-Kaas, Def. 2.7: a **liana** identifies two **arrows** —
a double derivation, i.e. the Laplacian; a **stolon** identifies two **nodes** — a double
evaluation, i.e. the scalar product.)*

**Loops.** For $M\in\operatorname{End}_D(D^n)$ one has $\operatorname{tr}_\mathbb R=\operatorname{tr}_{D/\mathbb R}\circ\operatorname{tr}_D$ with
$\operatorname{tr}_{D/\mathbb R}(a+b\varepsilon)=2a$, because multiplication by $a+b\varepsilon$ has matrix
$\left(\begin{smallmatrix}a&0\\b&a\end{smallmatrix}\right)$. In general $\operatorname{tr}_{A/\mathbb R}(u)=(\dim_\mathbb R A)\operatorname{aug}(u)$,
since multiplication by a nilpotent has trace $0$. Hence the $\varepsilon$-part is **discarded** — which is why a loop evaluated on a lifted field is
fibre-free — and each loop scales by $\dim_\mathbb R A$. Verified: ratios $2,4,8$ for $A=D$ and ratio
$3$ for $\dim_\mathbb R A=3$.

**No power-of-two law survives at all.** Sharper than the loops-vs-aromas correction below: the
theta graph $\lvert D^2X\rvert_F^2\,X$ — two nodes joined by two lianas and one stolon — has base
multiplier exactly $\mathbf 3$ (verified). The true multiplier is
$$N_0(\gamma)=\#\{\text{families of vertex-disjoint \emph{generalised cycles}}\},$$
a generalised cycle being a closed walk that traverses arrows head-to-tail and *reverses direction*
at a liana (two lower ends) or a stolon (two upper ends). Half-edge counting gives
$\#\text{arrows}_1+2\,\#\text{stolons}_1=\#\text{arrows}_1+2\,\#\text{lianas}_1$, so a non-empty
circulation needs either a directed arrow-cycle **or** at least one liana *and* one stolon. On
stolon-free forests this collapses to $N_0=2^{\#\text{aromas}}$ — which is exactly why the naive law
looked right for so long, and exactly why it broke once stolons entered.

**The multiplier counts loops, not aromas.** In the $GL$/aromatic world the two coincide: since
$\sigma$ is a fixed-point-free involution the root's upper index is consumed by the ghost arrow,
so every non-root component carries exactly one cycle, and a forest with $\alpha$ aromas scales by
$(\dim_\mathbb R A)^{\alpha}$. In the exotic ($O(n)$) world they **diverge**, and both halves of that are verified on the zero
section: $\langle X,X\rangle$ is a loopless aroma with multiplier $1$, while
$\lvert DX\rvert_F^2\,X$ — a liana feeding a stolon, with *no loop at all* — has multiplier
exactly $2$. So the grading automorphism $\delta_N$ and the slogan
"$\operatorname{Fix}(\delta_2)=$ aroma-free" are valid **only within the aromatic class**; in the exotic
class the multiplier is a count of admissible $\varepsilon$-level assignments, not a power of $N$.

**Stolons.** The real obstruction is not "two functionals" but **$D$-balancedness**. Any pairing
of the form $\lambda(\langle a,b\rangle_D)$ satisfies $\langle\varepsilon a,b\rangle=\langle a,\varepsilon b\rangle$. The
Euclidean (Sasaki) form on $\mathbb R^{2n}$ has $\langle\varepsilon a,b\rangle_\mathbb R=\sum a_0b_1\neq\sum a_1b_0=\langle a,\varepsilon b\rangle_\mathbb R$
— **not balanced**, hence not of that form for any $\lambda$. Metric contractions of evaluation
slots therefore cannot survive base change.

**Lianas survive by an accident.** The Laplacian passes *not* for a functorial reason but
because $c_D=1$, i.e. exactly because $\varepsilon^2=0$. This is precisely why the Laplacian method is
(T)-natural yet not $T^A$-natural once $c_A\neq1$ (§7).

---

## 5. The characterisation is relative, and (T) alone is weak

(T) does **not** imply "B-series". Three independent witnesses are local, consistent and
(T)-natural without being B-series methods: **partitioned RK**, the **Laplacian method**, and
tower-incoherent constructions. Indeed, since (T) constrains only the doubling tower, one may
choose for each odd $m$ an arbitrary partition of $\{1,\dots,m\}$ with an arbitrary PRK per
block and propagate along each tower by the $D$-lift: every such sequence is local, consistent
and (T)-natural, so the (T)-natural class has **the cardinality of the continuum**.

Far more strongly:

> **Extension theorem.** (T) is *vacuous on any single dimension*. Every local smooth $\Psi_q$ in
> odd dimension $q$ extends to a (T)-natural sequence: $T^A$ of a merely smooth map exists by the
> finite Taylor formula (as $\mathfrak m$ is nilpotent), and setting
> $$\Psi_{2q}(Y)(x_0,v_0):=T\bigl(\Psi_q(Y_1(\cdot,v_0))\bigr)(x_0,v_0)$$
> is local — note that the naive $Y_1(\cdot,0)$ would *not* be — smooth, and reduces to
> $T(\Psi_q(X))$ when $Y=TX$, since then $Y_1(\cdot,v_0)=X$. Iterate up the tower.

So the map "(T)-natural sequences $\to$ arbitrary sequences in odd dimensions" is **surjective**.
There is also an explicit *polynomial* (T)-natural method with no equivariance at all, obtained
by making the offending scalar ring-valued rather than real-valued.

The useful statements are therefore *relative* to an equivariance class, and — because (T) is
cross-dimensional — that class must be an affine group **scheme**, so that $H(D)$ exists and
the group acting in dimension $2n$ is $H(D)$, not $H(\mathbb R)$ in dimension $2n$.

| equivariance imposed (+ locality, trivial decoupling) | ambient class | (T) selects |
|---|---|---|
| $GL(n)$ | aromatic B-series (Munthe-Kaas–Verdier) | **B-series** (aromas killed by the factor $2$) |
| $O(n)$ | exotic aromatic B-series (Laurent–Munthe-Kaas Thm 2.12) | **exotic B-series** — trees + lianas; loops doubled away, stolons destroyed |
| $GL(n_1)\times GL(n_2)$ | aromatic P-series *(classification open in print)* | **P-series** *(conjectural)* |

The partitioned row carries a further caveat: "blockwise" is a hypothesis, not a theorem. PRK is
a (T)-natural family only when indexed by **partitioned** dimensions $(n_1,n_2)\mapsto(2n_1,2n_2)$;
a family indexed by plain dimension with a fixed first-half partition is not stable under base
change, since the induced partition of $D^{n_1}\times D^{n_2}$ interleaves. Note also that
splitting/IMEX methods are not a counterexample to the row at all: they are maps of a
**pair** of fields on one space (additive / bicoloured B-series), a different functor, rather than
equivariant maps of a single field on a partitioned space. The group must also be taken
**affine**, $(GL(V_1)\ltimes V_1)\times(GL(V_2)\ltimes V_2)$: with the linear group alone, explicit
$x$-dependence survives the classification.

Two remarks worth recording.

* **A limitation not previously flagged.** In the application where exotic series actually
  arise (Laurent–Munthe-Kaas §4.3; Laurent–Vilmart), the fields are gradients $X=\nabla V$, and
  there the classification collapses. But **the gradient class is not stable under $T$**: if
  $X=\nabla V$ then $TX(x,v)=(\nabla V(x),\nabla^2V(x)v)$, whereas
  $\nabla_{(x,v)}\langle\nabla V(x),v\rangle=(\nabla^2V(x)v,\nabla V(x))$. So (T) is not even well posed
  on $\mathfrak X^{\nabla}$, and their Thm 4.6 cannot be combined with the row above.
* **A bridge to their categorical properties.** Laurent–Munthe-Kaas Prop. 4.3 states that
  *exotic B-series are right-orthogonal-equivariant* (equivariance under affine $A$ with
  $AA^\top=I$). So within the local, orthogonal-equivariant, trivially decoupling class, (T) is
  equivalent — at the level of Taylor expansions — to **right-orthogonal-equivariance**. This
  identification does not appear in the literature.
* **The $GL(n_1)\times GL(n_2)$ row is genuinely open.** Their Conclusion lists partitioned
  B-series among the extensions *not yet* given an equivariance characterisation. The ambient
  classification must be established before the P-series row can be asserted; the mathematics
  is routine-looking (colour-preserving upper–lower contractions only, so no lianas or stolons
  arise and aromas are the sole obstruction), but it is not a citable theorem.

---

## 6. (T) is not a jet condition

Every classification theorem in this area (MMMV; Munthe-Kaas–Verdier; Laurent–Munthe-Kaas)
classifies the **Taylor expansion at the zero vector field**. (T) is an *exact identity at each
fixed $h$*. The two are not interchangeable. Let $\theta(t)=e^{-1/t^2}$, $\theta(0)=0$, and put
$$\Psi^X_h(x)=x+hX(x)+\theta\bigl(h\operatorname{div}X(x)\bigr)\,hX(x).$$
This is local, smooth, consistent, $GL$-equivariant, and **every** Taylor coefficient equals
explicit Euler's — so it is a B-series map in the sense of MMMV Def. 2.3. But
$\operatorname{div}(TX)=2\operatorname{div}X$ and $\theta(2t)\neq\theta(t)$, so **(T) fails at every $h\neq0$**.

Consequently:

* "aroma-free $\Rightarrow$ (T)" holds **for jets**, and on the nose only for methods *defined*
  in the language of §3.
* Peetre/extension arguments do not close the gap: they convert locality into finite jet order
  *in $x$*, after the transfer step has already replaced the method by its Taylor terms.
* Since $X\mapsto TX$ is $\mathbb R$-linear, (T) *does* descend to every homogeneous Taylor term.
  The failure is entirely in the flat remainder.
* The same example shows the **converse** directions of MMMV Thm 2.4 and Laurent–Munthe-Kaas
  Thm 2.13 ("B-series map $\Rightarrow$ affine / semi-orthogonal equivariant") are false as
  literally stated with those papers' own definitions: it is a B-series map by MMMV Def. 2.3
  yet is not decoupling. Those implications must be read as statements about the Taylor
  expansion, with the map identified with its series. The substantive direction
  (equivariance $\Rightarrow$ series of the stated type) is untouched.

Two distinct pathologies must be kept apart: **flatness in the field**, a $C^\infty$ artefact that
disappears in the analytic category; and **dimension-incoherence** (§5), which analyticity does
*not* cure, since the PRK-tower family is polynomial in $X$. Note also that flatness *in $h$* is
harmless: $\Psi^X_h=\varphi^X_{h+e^{-1/h^2}}$ is (T)-natural and perfectly definable.

---

## 7. Beyond the dual numbers

$T^{A\otimes B}=T^A\circ T^B$, so (T) gives $T^{D^{\otimes m}}=T^m$-naturality for free — all
iterated tangent bundles. It does **not** give $T^A$ for general Weil $A$.

> **Witness.** The Laplacian method satisfies (T) and every $T^m$, yet fails $T^A$-naturality
> for $A=\mathbb R[\varepsilon]/(\varepsilon^3)$ **with the standard identification**: the defect is concentrated
> in the $\varepsilon^2$ block and equals exactly $\Delta X$. (Verified.) The qualifier is not cosmetic —
> by §4 there is a basis of that same algebra in which $c_A=1$ and the method passes.

The mechanism is the liana constant of §4: $c_{D^{\otimes m}}=1$, so lianas survive every $T^m$,
whereas $c_{\mathbb R[\varepsilon]/(\varepsilon^3)}=1+\varepsilon^2$, so the liana realizes $(1+\varepsilon^2)\Delta_A$ instead
of $\Delta_A$ and the defect is $\varepsilon^2\,T^A(\Delta X)$ — exactly the verified $\varepsilon^2$-block
discrepancy. This *derives* the failure rather than asserting it, and it is the operative
reason; the observation that $\mathbb R[\varepsilon]/(\varepsilon^{r+1})\hookrightarrow D^{\otimes r}$ is a subalgebra
rather than a quotient is a corollary, not the cause.

*Positive counterpart.* (T) together with equivariance under linear permutation matrices gives
$T^{(r)}$-naturality for all $r$ — every higher variational equation — because
$\mathbb R[s]/(s^{r+1})\subset D^{\otimes r}$ is the $S_r$-invariant subalgebra.

*What the tower of algebras is.* For $A_r=\mathbb R[\varepsilon]/(\varepsilon^{r+1})$ and
$u=\sum_k\varepsilon^ku_k$, expanding $X(u)$ over $A_r$ gives, by Faà di Bruno,
$$\dot u_k=\sum_{j\ge1}\frac1{j!}\sum_{\substack{i_1+\cdots+i_j=k\\ i_l\ge1}}D^jX(u_0)[u_{i_1},\dots,u_{i_j}],$$
so the realization of $X^{A_r}$ is exactly $(X,VE_1,\dots,VE_r)$ — the higher variational equations
along the solution, the objects of Morales–Ramis–Simó and of jet-transport integration (verified to
$r=3$). Algebraic naturality therefore reads: *the method commutes with passing to variational
equations*; that the exact flow does so is the fact that $VE_k$ is solved by $\partial^k\varphi$.
Two cautions before importing anything from differential Galois theory. These algebras are **not**
differential rings — base change is a change of scalars with $\partial_t$ acting trivially, so every
element of $A$ is a constant and $\operatorname{Aut}(A)$ is not the Galois group of anything here;
and the Galois group of a *discretised* system is a difference-Galois group, so a "discrete
Morales–Ramis" transfer would compare two theories rather than assert an identity. What genuinely is
shared is the functor-of-points move: $\operatorname{Lie}(G)=\ker\bigl(G(C[\varepsilon])\to G(C)\bigr)$
for the Picard–Vessiot group scheme is the same use of $D$ that §5 needs in order for $H(D)$ to
exist.

### Affine rigidity, and the ceiling of the two-lift mechanism

Every clause in the definition of algebraic naturality (§11) ties $\Psi_{nN}$ to $\Psi_n$; taken one
at a time, none of them constrains a single $\Psi_n$ at all. A **per-dimension** consequence can
therefore only come from a *collision*: two based algebras $(A,\iota)$, $(B,\kappa)$ of the same
dimension $N$ whose lifts of one and the same field agree. Then $X^A=X^B$ is a single field, the map
$\Psi_{nN}$ assigns to it is a single map, and both prescriptions must return it:
$$X^A=X^B\quad\Longrightarrow\quad T^A(\Psi^X_n)=T^B(\Psi^X_n).$$
Base change is the same Taylor formula for a vector field as for a map, so the collision locus
$\mathcal C(A,B)=\{F:F^A=F^B\}$ is cut out by the same conditions on both sides: *algebraic
naturality preserves every collision class.*

> **Theorem (collision trichotomy).** For based algebras $(A,\iota),(B,\kappa)$ of equal dimension
> $N$, with $1_A,1_B\in\mathbb R^N$ the coordinate vectors of the two units, and $F$ polynomial,
> $$\mathcal C(A,B)=\begin{cases}\text{everything}&\text{if the based algebras coincide,}\\
> \text{the affine maps}&\text{if they differ and }1_A=1_B,\\
> \text{the linear maps}&\text{if }1_A\neq1_B.\end{cases}$$

*Proof.* Base change of a polynomial map is that polynomial evaluated over $A$:
$F^A(u)=\sum_k\frac1{k!}\partial^kF(0)[u^{\otimes k}]$, powers taken in $A$. Its $k$-th term is
homogeneous of degree $k$ in $u\in\mathbb R^{nN}$, so $F^A=F^B$ holds degree by degree. Degree $0$ reads
$F(0)\otimes 1_A=F(0)\otimes 1_B$; degree $1$ is $\operatorname{id}\otimes\partial F(0)$ on both sides,
always equal. For $k\ge2$ it reads
$\sum\delta^{\gamma}_{\alpha_1\cdots\alpha_k}\partial^kF(0)[u^{\alpha_1},\dots,u^{\alpha_k}]=0$, with
$\delta$ the difference of the $k$-fold structure constants. If the two multiplications differ,
choose $\gamma,\alpha,\beta$ with $\delta^\gamma_{\alpha\beta}\neq0$ and set $u^\alpha=w_1$,
$u^\beta=w_2$, all other blocks $0$: the $k=2$ condition becomes $\partial^2F(0)[w_1,w_2]=0$ for all
$w_1,w_2$, so $\partial^2F(0)=0$; the same at every base point, so $F$ is affine — whereupon every
higher condition holds automatically. If the multiplications agree, the based algebras agree, units
included. Degree $0$ then forces $F(0)=0$ exactly when $1_A\neq1_B$. $\square$

> **Corollary (affine rigidity).** If $\Psi$ is algebraically natural then for every $n$, every
> **affine** field $X$ on $\mathbb R^n$ and every $h$, the map $\Psi^X_h$ is **affine**.

The cheapest witness sits in dimension **two**, and writing it in coordinates explains the earlier
result. Realise $D$ in the basis $(1,\varepsilon)$ and $\mathbb R^2$ in the unit-first basis $(1,p)$ with
$p=(0,1)$:
$$X^{D}(u^0,u^1)=\bigl(X(u^0),\;DX(u^0)u^1\bigr),\qquad
X^{\mathbb R^2}(u^0,u^1)=\bigl(X(u^0),\;X(u^0+u^1)-X(u^0)\bigr).$$
One is a tangent, the other a **secant**, and they agree exactly on affine fields. The same identity
applied to $\Psi^X$ reads $\Psi^X(a+b)-\Psi^X(a)=D\Psi^X(a)\,b$, so $\Psi^X$ is affine — no jet
argument, no regularity beyond differentiability. In the **standard** basis of $\mathbb R^2$ the unit is
$(1,1)\neq(1,0)$, the trichotomy drops to its third case, and what survives is exactly the
previously recorded **linear rigidity**: that result was the unit-misaligned shadow of this one. No
convention is needed at all in dimension **three**, where $\mathbb R[\varepsilon]/(\varepsilon^3)$ and
$\mathbb R[x,y]/\mathfrak m^2$ both carry $1$ as the first vector of their standard bases and collide
precisely on the affine maps.

> **Corollary (ceiling).** No collision of two lifts of one field constrains $\Psi$ on any field
> that is not affine.

— because distinct based algebras already differ on a product of *two* basis vectors, so the
obstruction is always second order. "Affine $\Rightarrow$ affine" is therefore the last consequence
available from this mechanism, not merely the best one found. It is corroborated from the other
side: every per-dimension consequence must hold for the exact flow *and* for every Runge–Kutta
method, and RK methods raise polynomial degree ($u+hX(u+hX(u))$ is quartic on a quadratic field), so
no statement of the form "degree $\le d$ is preserved" can survive for $d\ge2$.

**Calibration — where the new content lives.** Affine rigidity is strictly stronger than linear
rigidity *as a condition on families*: $\Psi^X_h(u)=u+hX(u)+h^2\langle X(u)-DX(u)u,\,u\rangle X(u)$ is
local, satisfies linear rigidity — the bracket vanishes identically on linear fields by Euler's
identity — and fails affine rigidity. That witness is not translation equivariant, and no witness
can be: for a local method $\Psi^X_h(u)=u+\Phi(j^rX(u),h)$ the affine field $Mu+c$ occupies the same
jet slots as the linear field $Mu$ and differs only in the value, so linear rigidity already forces
$\Phi(\cdot,M,0;h)$ to be linear in its first argument and hence $\Psi^{Mu+c}$ affine. The gap
between the two rigidities is exactly as wide as the failure of translation equivariance — which is
the right size, since the point of a per-dimension consequence is that it assumes no equivariance
whatever.

**Cross-factorization adds nothing** *(evidence, not proof)*. A collision may also join *different*
source dimensions: $X^A=Y^B$ with $X$ on $\mathbb R^{n_1}$, $Y$ on $\mathbb R^{n_2}$,
$n_1\dim A=n_2\dim B$ and $n_1\neq n_2$, forcing $T^A\Psi^X_{n_1}=T^B\Psi^Y_{n_2}$. Over the
enumerated range — $m=4$ pairing $(n_2,\dim B)=(1,4)$ against $(n_1,\dim A)=(2,2)$, and $m=6$ pairing
$(2,3)$ against $(3,2)$; 16 based algebras; all homogeneous $Y$ of degree $\le3$, by exact kernel
computation — every collision carrying a non-affine field is **functorial**, in one of exactly two
ways: $B\cong A\otimes C$ with $X=Y^C$, so both sides are $T^A$ applied to one condition; or $X$ and
$Y$ are both lifts of a common field $W$, so both sides are the condition at $W$. Either way the
relation holds for every algebraically natural $\Psi$ and constrains nothing. The $m=6$ column is
the sharper test: there $\dim C=3/2$, so no tensor factorization is available at all, and only the
second kind occurs. What is missing for a proof is that a lifted field determines its presentation
up to common refinement.

### The closure spectrum

For a family $\Psi$ put
$$\mathcal W(\Psi):=\{(A,\iota)\ :\ \Psi\text{ is }T^A\text{-natural}\},$$
its **closure spectrum**. Algebraic naturality is $\mathcal W=$ everything; (T) alone is
$D\in\mathcal W$. Asking which classes occur is the Galois-correspondence-shaped question in this
subject, and the spectrum has real structure.

> **Proposition (monoidality).** $\mathcal W(\Psi)$ is a monoid under $\otimes$.

*Proof.* $\Psi^{X^{A\otimes B}}=\Psi^{(X^B)^A}=T^A\Psi^{X^B}=T^AT^B\Psi^X=T^{A\otimes B}\Psi^X$.
$\square$

Structurally this is because the three invariants of §4 are multiplicative:
$\dim(A\otimes B)=\dim A\dim B$, $c_{A\otimes B}=c_Ac_B$, and $L_{a\otimes b}=L_a\otimes L_b$ is
self-adjoint when both factors are. All three are verified, along with $T^{A\otimes B}=T^A\circ T^B$
holding on the nose in the tensor basis ordered $e_\alpha\otimes f_\beta\mapsto\beta N_A+\alpha$.

**The invariants are metric data, not basis data.** $c_A=\mu(g^{-1})$ and balancedness depend only
on the pair $(A,g)$, $g$ being the inner product for which the chosen basis is declared orthonormal;
they are therefore invariant under $O(N)$ changes of basis and nothing more. §4's warning is exactly
this, and it is now systematic rather than anecdotal: $\mathbb R^2$ in its standard basis has $c=1$ and
is balanced, in the unit-first basis $(1,p)$ it is neither, and the two have different spectra.

**The survival law is coarser than §4's table suggests** *(correction)*. Measured on 17 based
algebras against six methods: a **loop or a stolon** survives only for $A=\mathbb R$, while a **liana**
survives iff $c_A=1$. Balancedness is necessary for the real form to factor as
$\lambda\circ\langle\cdot,\cdot\rangle_A$ — which is what §4 uses to explain the failure over $D$ —
but it is **not sufficient** for a stolon to survive: $\mathbb R^k$ in its standard basis is balanced
*and* has $c=1$, and still kills stolons. The reason unifies the loop and stolon rows and is scalar
rigidity: both contractions deliver a **real** scalar $r$ where base change demands the $A$-scalar
$q$, and the defect is exactly
$$h^2\,(r\cdot 1_A-q)\cdot X^A$$
— verified as an identity over $\mathbb R^2$. A loop differs only in that $r=\dim_\mathbb R A\cdot\operatorname{aug}(q)$,
which is why its defect presented as a clean multiplier.

> **Corollary.** If the $h^2$-coefficient of $\Psi$ is a combination of exotic-aromatic elementary
> differentials, $\mathcal W(\Psi)$ is everything (no decorations), $\{A:c_A=1\}$ (lianas only), or
> $\{\mathbb R\}$ (some loop or stolon).

**Not closed under subalgebras.** $\mathbb R[s]/(s^3)$ is the $S_2$-invariant subalgebra of $D\otimes D$,
and the Laplacian method closes under $D\otimes D$ but not under it. The failure survives the obvious
repair: with the metric *induced* from the ambient orthonormal basis, $g=\operatorname{diag}(1,2,1)$,
one gets $c=1+\tfrac12s^2=1+xy\neq1$. So $\mathcal W$ is a $\otimes$-monoid and nothing finer — $c$
is not inherited by subalgebras.

**The realized spectra form a lattice, not a chain** — and a correction is due here. §5 obtains a
polynomial (T)-natural method with no equivariance "by making the offending scalar ring-valued rather
than real-valued". The round-4 notes then attached that label to $u+hX+h^2(X^1)^2X$ with the scalar
left **real**-valued, which is **not** (T)-natural at all: its defect is the fibre term
$2X^1(DX^1v)X$ (measured). The genuine witness follows §5 literally. Every $m$ is uniquely $2^kq$
with $q$ odd, so set
$$\Psi_m:=\text{the }D^{\otimes k}\text{-base change of }u+hX+h^2(X^1)^2X\text{ on }\mathbb R^{q}.$$
This is polynomial, local, defined in every dimension, manifestly non-equivariant, and
$T^{D^{\otimes k}}$-natural for every $k$ (verified) — and it fails linear rigidity, so the round-4
conclusion that algebraic naturality kills it stands, now with a witness that is actually
(T)-natural. Replacing $D$ by any based $A$ gives $\mathcal W\supseteq\langle A\rangle_\otimes$, and
$\langle D\rangle_\otimes$, $\langle\mathbb R^2\rangle_\otimes$ are **incomparable**: each contains an
algebra the other omits. With $\{\mathbb R\}$, $\{c_A=1\}$ and everything, the image of
$\Psi\mapsto\mathcal W(\Psi)$ is therefore a lattice rather than a chain. Which $\otimes$-closed
classes are realized is open; the tower construction realizes every principal one.

### Presentations, and where the gluing stands

A *presentation* of $Z\in\mathfrak X(\mathbb R^m)$ is a triple $(n,(A,\iota),X)$ with $n\dim A=m$ and
$Z=X^A$. Whether two presentations of one field must have a common refinement is what decides
whether the conditions algebraic naturality places on a single $\Psi_m$ are mutually consistent.

> **The commutant is the presentation.** Put
> $$\mathcal A(Z):=\bigl\{L\in\operatorname{End}(\mathbb R^m)\ :\ L\,D^kZ(u)[v_1,\dots,v_k]=D^kZ(u)[Lv_1,v_2,\dots,v_k]\ \ \forall k\ge1,u,v_i\bigr\}.$$
> It is a unital subalgebra of $\operatorname{End}(\mathbb R^m)$, and it is computable: the conditions are
> linear in $L$. If $Z=X^A$ then $\rho_A(A)\subseteq\mathcal A(Z)$, because base change makes every
> $D^kZ(u)$ an $A$-multilinear $A$-valued map.

Measured: for a generic cubic field and each of the 17 based algebras, $\mathcal A(X^A)$ is
**exactly** $\rho_A(A)$ — a generic lifted field determines its presentation outright. And in both
coincidence families above, $\mathcal A(Z)$ is the algebra of the *finest* presentation and contains
the actions of both coarser ones: $Z=Y^{D\otimes D}$, presented over $(1,D\otimes D)$ and $(2,D)$,
has $\dim\mathcal A=4$; $Z=(W^{\mathbb R^2})^{\mathbb R^3}$, presented over $(2,\mathbb R^3)$ and
$(3,\mathbb R^2)$, has $\dim\mathcal A=6=\dim(\mathbb R^3\otimes\mathbb R^2)$.

> **Consequence (the overlap condition).** If $\Psi$ is natural below dimension $m$ and
> $Z=X^A=Y^B$, both prescriptions for $\Psi_m(Z)$ factor through the commutant presentation
> $(n_0,E,W)$: from $X=W^C$, $T^A\Psi(X)=T^AT^C\Psi(W)=T^E\Psi(W)$, and likewise for $B$. So the
> prescription on the union of lifted loci is **consistent**, by induction on dimension.

That is the input the gluing needed, and it closes open item 3 up to one point: that $\mathcal A(Z)$
is always commutative with $\mathbb R^m$ free over it, so that the finest presentation exists.
Verified in all 19 cases here; not proved.

**Zero-field rigidity — with a correction.** The zero field is lifted by *every* algebra, $0^A=0$,
so every pair of based algebras of one dimension collides there. Base change of an affine map gives
$S_{nN}=I_N\otimes S_n$ and $w_{nN}=1_A\otimes w_n$; reading the second for $A=D$ against
$\mathbb R^2$ in its standard basis ($1_D=(1,0)$ but $1_{\mathbb R^2}=(1,1)$) forces $w=0$, and the
first at $n=1$ gives $S_m=S_1I_m$.

> **Corollary.** $\Psi^0_h(u)=s(h)\,u$, with $s\equiv1$ once consistency is imposed.

*What is new here, and what is not.* The zero field is **linear**, so $w=0$ is nothing but linear
rigidity at $M=0$. An earlier draft claimed that $u+hX+h^2e_1$ "escapes affine rigidity and is killed
by the zero-field collision"; that is misleading — linear rigidity kills it directly, since
$\Psi^0_h(u)=u+h^2e_1$ is not a linear map. The claim survived a first pass only because the check
tested $\max_u\deg=1$, which admits a constant term; it now requires every monomial to have
$u$-degree exactly one. The genuinely new content is $S_m=s(h)I_m$ — the linear part on the zero
field is a **scalar**, forced by the $n=1$ instance of the tensor compatibility — so the family
carries no distinguished endomorphism either.

**Where the gluing stands.** Building a non-equivariant algebraically natural family (open item 2)
has the shape: $\Psi_1$ is free up to affine rigidity, and for $m>1$, $\Psi_m$ is prescribed on
$L_m=\bigcup_{(n,A),\,\dim A>1}\{X^A\}$ and free off it. Consistency is now settled, so what remains
is a **smooth and local** extension across $L_m$.

Locality makes that fibrewise: $\Psi_m$ is a function of the jet, and $V_A(u)=\{j^r(X^A)(u):X\}$ is a
*linear* subspace of $J^r_u$, base change being linear in $X$. Inclusion–exclusion over finitely many
subspaces with compatible projections would then give a smooth local extension — and one can say
exactly why that route fails: $\dim V_A(u)$ is **not** locally constant. At a real point
$u=1_A\otimes x$ the Jacobian of $X^D$ is block diagonal and the higher jet of $X$ is forgotten, so
for $A=D$, $n=2$, $r=1$ the dimension drops from $10$ to $6$ (measured); each algebra degenerates on
its own real locus, which is precisely where the loci meet. Two further obstacles: the projections
would have to be compatible with locality, and once $\dim A\ge7$ there are *moduli* of commutative
algebras, so $L_m$ is a union of a continuum of loci rather than finitely many.

What survives is §5's retraction, which handles one algebra at a time:
$\Psi_{nN}(Z)(u):=T^A\bigl(\Psi_n(X_{Z,u})\bigr)(u)$, freezing the nilpotent part of $u$, is smooth,
local, and correct on $\{X^A\}$. So the open problem is now sharp: **is there a simultaneous
retraction** — one smooth local formula restricting correctly on every $\{X^A\}$ at once? The
plausible route is to support the junk near a field far from every lifted locus *and* from every
retraction image, then propagate; that is an inductive genericity argument, and the moduli in high
dimension are what make it delicate.

**Only the tensor relation transfers.** The Weil algebras form a category, and one would like to
test naturality on a generating family: every Weil algebra is a quotient of a jet algebra
$\mathbb R[x_1,\dots,x_k]/\mathfrak m^{r+1}$, so quotient-closure would reduce everything to the jet
algebras together with finite products. It fails — and so does every relation except $\otimes$.

| relation | transfers? | witness |
|---|---|---|
| $A,B\in\mathcal W\Rightarrow A\otimes B\in\mathcal W$ | **yes** | one line, from $T^{A\otimes B}=T^AT^B$ |
| $A\in\mathcal W$ and $B\subseteq A$ | no | Laplacian; $\mathbb R[s]/(s^3)\subset D\otimes D$ |
| $A\in\mathcal W$ and $A\twoheadrightarrow B$ | no | $\operatorname{tower}(D\otimes D)$; $D\otimes D\twoheadrightarrow D$ by $x,y\mapsto\varepsilon$ |
| $A,B\in\mathcal W\Rightarrow A\times B\in\mathcal W$ | no | $\operatorname{tower}(D)$; $D\times D$ |

The reason is uniform. $\otimes$ is the only one of these that is a relation between the *functors*,
$T^{A\otimes B}=T^A\circ T^B$, so one condition can be applied twice. Subalgebra, quotient and
product give natural transformations *between* functors, and naturality is a condition on each
functor separately — $\Psi_m$ and $\Psi_{m'}$ in different dimensions are a priori unrelated, so
nothing passes along a transformation. **Algebraic naturality therefore has no generating family of
test algebras**, which is one concrete reason the implication (i)$\Rightarrow$(ii) of §11 resists:
there is no finite set of algebras on which to run an argument.

**What the morphisms do give: a genuine $GL$.** For $\phi:A\to B$ the map $\phi_n=\phi\otimes\operatorname{id}$
is real-linear, intertwines $X^A$ with $X^B$, and satisfies $T^BF\circ\phi_n=\phi_n\circ T^AF$. So an
algebraically natural $\Psi$ obeys
$$\phi_n\circ\Psi^{X^A}_{n\dim A}=\Psi^{X^B}_{n\dim B}\circ\phi_n\qquad\text{on lifted fields,}$$
a consequence of naturality rather than an extra condition. For $A=B$ it reads: $\Psi_m$ is
**$\operatorname{Aut}(A)$-equivariant on the $A$-lifted locus** — and that is not vacuous. Since
$\operatorname{Aut}(\mathbb R[x_1,\dots,x_k]/\mathfrak m^2)=GL(k)$, naturality under that one algebra
forces honest $GL(k)$-equivariance there (verified: the lifted field is $GL(2)$-invariant, Euler
commutes with the action, a method singling out one algebra block does not). This is the first
genuine $GL$ to appear in this programme, and it acts in the **algebra** directions rather than the
$\mathbb R^n$ directions — precisely the gap §11 records against "algebraically natural
$\Rightarrow$ affine equivariant".

### What $\operatorname{Aut}$-equivariance imposes

The condition above reads like equivariance. It is much weaker, and one can say exactly how weak.

> **It is trivial on the base.** Every automorphism is unital, so $\phi_n$ fixes the real locus
> $1_A\otimes\mathbb R^n$ **pointwise**. Since $\Psi_{n\dim A}$ restricted to that locus is $\Psi_n$,
> the condition imposes nothing whatever on the base method: its entire content lies in the
> nilpotent directions.

Infinitesimally the same fact is $d(1)=d(1\cdot1)=2d(1)$, hence $d(1)=0$, together with
$d(\mathfrak m)\subseteq\mathfrak m$: a derivation has zero $e_0$-row and column. Computed:
$\operatorname{Der}(J^1_k)=\mathfrak{gl}(k)$ of dimension $k^2$;
$\operatorname{Der}(\mathbb R[\varepsilon]/(\varepsilon^{r+1}))$ of dimension $r$; and
$\operatorname{Der}(\mathbb R^k)=0$ — a product algebra has no infinitesimal automorphisms at all,
only the finite $S_k$.

> **What it generates is $\mathfrak{gl}(m-1)$, not $\mathfrak{gl}(m)$.** Over *all* factorizations
> $m=n\dim A$ the generated Lie algebra is exactly $\mathfrak{gl}(m-1)$, acting on the nilpotent
> directions and fixing the unit direction: measured $9,\,25,\,49$ at $m=4,6,8$ against
> $\mathfrak{gl}(m)=16,\,36,\,64$, and against the stabiliser of $e_0$, larger still at $m^2-m$. It
> is realised already by $J^1_{m-1}$ alone, so combining factorizations gains nothing.

*This corrects half of a recorded reason.* §11 locates the failure to recover $GL(n)$-equivariance
in "the group generated is large enough ($\mathfrak{gl}(nN)$) but the lifted locus too thin". For the
automorphism-generated group that is **false**: the group is not large enough, and thinness is not
the operative obstruction. It fixes the base pointwise, so no amount of it can produce equivariance
in the $\mathbb R^n$ directions — the wrong group, not a large group on a small set.

**What it does impose.** For $A=J^1_k$ the lift is the $k$-fold Whitney sum,
$$X^{J^1_k}(u^0,u^1,\dots,u^k)=\bigl(X(u^0),\ DX(u^0)u^1,\ \dots,\ DX(u^0)u^k\bigr),$$
and $\operatorname{Aut}(J^1_k)=GL(k)$ mixes the $k$ tangent copies. So the content is exactly:
**$\Psi$ may not distinguish directions inside the nilpotent part** — no preferred basis of
$\mathfrak m$ — and nothing else.

### The pointwise detector, and a proposed simultaneous retraction

At second order — write $\Psi^X_{m,h}(u)=u+hX(u)+h^2B_m(X)(u)+O(h^3)$, so naturality is exactly
$B_{n\dim A}(X^A)=(B_nX)^A$ — the retraction problem changes character, because locality makes the
relevant object the **pointwise** commutant
$$\mathcal A_r\bigl(u,j^rZ(u)\bigr)=\bigl\{L:\ L\,D^kZ(u)[v_1,\dots,v_k]=D^kZ(u)[Lv_1,v_2,\dots,v_k],\ 1\le k\le r\bigr\}.$$

> **Verified.** $\dim\mathcal A_3(u,j^3(X^A)(u))=\dim A$, and the commutant *is* $\rho_A(A)$, at a
> generic point, at the real locus $1_A\otimes\mathbb R^n$, and at the origin — all 17 based
> algebras, $n=1,2$. So the pointwise commutant does **not** collapse on the real locus, whereas
> $\dim V_A(u)$ does — and that collapse is what killed the fibrewise-linear route above. It is
> multiplicative, $\dim\mathcal A(X^A)=\dim A\cdot\dim\mathcal A(X)$, and it *separates* the
> factorizations: the $\rho_A(A)$ are pairwise distinct (11 factorizations at $m=4$, 8 at $m=6$).

> **Verified.** The prescription needs no per-algebra data. With $\mathfrak a=\mathcal A_r$ and
> generators $g_j=1_A\otimes e_j$, a **polynomial** field $Y$ on $\mathbb R^n$ satisfies
> $Y^A(u)=\sum_iY_i(u_1,\dots,u_n)\cdot g_i$, products taken in $\mathfrak a$ — no basis of $A$, no
> structure constants, no augmentation. Checked against `lift_map` for two different $Y$, over all 17
> algebras, at all three kinds of point.

Together these suggest summing over the **divisors** $N\mid m$ rather than over algebras,
$$\Phi_m:=\kappa_m+\sum_{N\mid m,\ N\ge2}\chi_N(J)\,E_{m/N,N}(u,J),$$
with $E$ the intrinsic recover-and-lift, $\chi_N$ a smooth cutoff keyed to $\dim\mathcal A_r$, and
$\kappa_m$ free junk supported in $\Omega_m=\{\dim\mathcal A_r=1\}$ — nonempty, disjoint from every
lifted locus, and containing no affine field, so junk there is affine-rigid for free. If it works,
both recorded obstacles dissolve, and a witness ($\kappa_1=0$, $\kappa_2$ a bump in $\Omega_2$) would
refute "algebraically natural $\Rightarrow$ affine equivariant" for smooth families.

**This is a lead, not a result, and three gaps separate it from one.**

1. *The cutoffs.* $\{\dim\mathcal A_r=N\}$ is locally closed, not open, and on a lifted locus has
   empty interior. Keying $\chi_N$ to a spectral gap is the standard repair, but it must also be made
   mutually exclusive — first $N$ singular values below $\varepsilon$ **and** the $(N{+}1)$-st above
   $\delta$ — or a jet with $N$ small values could satisfy a gap condition at some $N'>N$ and pick up
   a spurious term.
2. *Commutativity and freeness.* The construction needs $\mathcal A_r(J)$ commutative with
   $\mathbb R^m$ free over it wherever the prescription is nonzero — precisely the residual point of
   the commutant theorem above, verified in cases and never proved, and now **load-bearing**.
3. *Polynomial versus smooth.* The intrinsic formula is verified for **polynomial** $Y$, where
   substitution needs no expansion point. The witness requires a compactly supported **smooth**
   $\kappa_2$, whose base change is defined by Taylor expansion about the real point — the datum the
   formula is advertised as not using. For an algebra with several local factors there is one real
   point *per factor*, the recorded trap that already invalidated two earlier forms of this
   prescription. The idempotents of $\mathfrak a$ are intrinsic, so this looks closable; it is not
   closed, and the verification does not reach it.

Until those close, the headline reading — that the open half of §11's conjecture is false for smooth
families — is **not** established. What *is* established is the detector: pointwise, non-degenerate,
multiplicative, separating, and per-algebra-free on polynomial data.

### Jet transport for a general method, and the shape it is relative to

The condition $\Psi_{n\dim A}(X^A)=T^A(\Psi_n(X))$ carries an asymmetry the equivariance framing
hides. On the left, $\Psi$ receives $X^A$ as a **bare** field on $\mathbb R^{n\dim A}$: it is not
told the field is a lift, nor from which $(n,A)$. On the right, $T^A$ is applied to the output using
the full structure. So the property says: *a method blind to the structure must nonetheless return
the structured answer.*

For a method that is an **algorithm over the base ring** this is automatic, and is exactly what jet
transport does in practice — run the same program with $\mathbb R$-arithmetic replaced by
$A$-arithmetic. The prerequisites are a short list.

| | requirement | why, and who violates it |
|---|---|---|
| **J1** | the field is touched only by **evaluation** of $X$ and finitely many derivatives, at computed points | evaluation at an $A$-point *is* $X^A$, by definition of base change — so it costs nothing |
| **J2** | only $\mathbb R$-algebra operations on values: $+$, $-$, real scalars, products. No norms, no order, no branching on them; **division only by units** | over $A$ a general element is not invertible: $a$ is a unit iff its image in each local factor's residue field is nonzero (verified) |
| **J3** | implicit definitions are admissible when the defining equation is natural in the base and has a **unique** solution | nilpotence of $\mathfrak m$ supplies this whenever the $\mathbb R$-linearization is invertible: solve over $\mathbb R$, lift by finitely many Newton steps. This is why implicit RK jet-transports |
| **J4** | any **auxiliary structure** must base-change functorially | a *splitting* does; a *metric* does not — precisely the $c_A$/balancedness/trace trichotomy of §4 |

Every recorded example falls out of the list: the exact flow, RK and Taylor use J1–J3 with no
auxiliary structure, so they are closed for every $A$, equivariant, B-series; **leapfrog** — and PRK
generally — uses J1–J3 *plus a splitting*, so it is closed and **not** affine equivariant, a
splitting being no $GL$-invariant; the **Laplacian method** uses the metric, violating J4, and
survives over $D$ only by the accident $c_D=1$; **aromatic** methods use the trace, violating J4, and
die by the factor $\dim_\mathbb R A$; and **adaptive step control** on
$\lVert\mathrm{err}\rVert<\mathrm{tol}$ violates J2 twice, by a norm and by a branch, so it does not
jet-transport at all.

**The missing ingredient is the shape.** J4 forces a restatement of the property itself. Leapfrog's
auxiliary datum is a splitting of the coordinates, and the condition only means anything once one
says *how the splitting is carried along the lift*. Verified: of the three ways to match a splitting
of $\mathbb R^2$ against the $D$-lift to $\mathbb R^4$, exactly one — the one base change dictates,
coordinate $j\mapsto$ the $A$-block $\{jN,\dots,jN+N-1\}$ — makes leapfrog closed; real-versus-fibre
and crossed both fail. So the honest statement of the property is not "a family indexed by dimension
commuting with base change" but

> a **shape functor** $\mathcal S$ — dimension, or dimension-with-a-splitting, or whatever structure
> the method needs — itself carrying a base change $A\mapsto\mathcal S^A$, together with a family
> $\Psi$ over $\mathcal S$ commuting with it.

Plain dimension is the trivial shape; partitioned dimension is leapfrog's. The recorded caveat that
"PRK is a (T)-natural family only when indexed by **partitioned** dimensions" is then not a caveat
but an instance: the shape is part of the data, and had been silently fixed.

**Where the strength is.** J1–J4 are *sufficient*. They are **not** necessary: the extension theorem
produces closed families by gluing rather than computing, and those are programs in no interface. So

$$\text{programmable}\ \Longrightarrow\ \text{closed},\qquad\text{and the gap is exactly the extension-theorem junk.}$$

This reframes the open problem in a way leapfrog does not refute. Instead of "does closure imply
affine equivariance" — false, leapfrog — the question becomes: **is every closed method given by an
algorithm one of J1–J4, relative to some shape?** Leapfrog is then an instance rather than a
counterexample, and the classification splits into two independent halves: *which shapes
base-change* (a question about natural structures, where §4's trichotomy answers it for tensors), and
*which programs over a fixed shape are closed*.

### The other Weil algebra, and why it is not a stronger test

Two different objects carry the name. Ours is **Weil's algebra of infinitely near points** — a
finite-dimensional commutative $\mathbb R$-algebra $\mathbb R\oplus\mathfrak m$ with $\mathfrak m$
nilpotent — and by Kolář–Michor–Slovák these classify the product-preserving endofunctors of
$\mathbf{Mf}$. The other is **Cartan's Weil algebra of a Lie algebra**,
$W(\mathfrak g)=\Lambda\mathfrak g^*\otimes S\mathfrak g^*$: an acyclic differential graded algebra
modelling $EG$, the object behind Chern–Weil theory and the Cartan model of equivariant cohomology.
They share a name and a namesake; I know of no theorem identifying them.

There *is* a genuine common home, and it is the graded one. Taking the generator of
$D=\mathbb R[\varepsilon]/(\varepsilon^2)$ to be **odd** rather than even replaces $TM$ by
$\Pi TM=T[1]M=\operatorname{Map}(\mathbb R^{0|1},M)$, whose function algebra is
$C^\infty(T[1]M)=\Omega^\bullet(M)$ with the de Rham differential as a homological vector field; and
$W(\mathfrak g)$ is the function algebra of $\mathfrak g[1]\oplus\mathfrak g[2]$. So both notions sit
inside "algebras of functions on graded or infinitesimally thickened points", which is presumably the
kinship being felt.

**It is nevertheless not a stronger condition here, and the reason is one line.** An ordinary
manifold's $A$-points are $(A_{\text{even}})^n$ — there are no odd parameters to pair with odd basis
vectors — so a super Weil algebra acts through its **even part**, which is an ordinary Weil algebra.
Verified: $\mathrm{even}(\Lambda\mathbb R^2)=D$, and $\mathrm{even}(\Lambda\mathbb R^3)$ has *exactly*
the structure constants of $J^1_3=\mathbb R[x,y,z]/\mathfrak m^2$ — an algebra already in the family,
indeed the one whose $\operatorname{Aut}=GL(3)$ generated the $\mathfrak{gl}(m-1)$ of the previous
subsection. $\mathrm{even}(\Lambda\mathbb R^4)$ and $\mathrm{even}(\Lambda\mathbb R^5)$ are ordinary
Weil algebras of dimensions 8 and 16. So the graded notion contributes no new test algebra.

Where it *would* add strength is not by supplying an algebra but by enlarging the category the
**method** acts on: methods defined on supermanifolds or graded bundles, whose input field itself has
odd components. That is a different — and larger — undertaking than adding a test object, and the
J-list of the previous subsection is what it would have to be re-run against, since a program in
$+,\times$ and evaluations extends to a supercommutative ring only once its signs are fixed.

And the topological notion does touch this project, but in **§8**, not here. §8 records that
$\tfrac12X'X$ is not a natural vector field — it needs a connection — and that by
Kolář–Michor–Slovák the natural operators $\mathfrak X\rightsquigarrow\mathfrak X$ are only
$X\mapsto cX$, so that the affine category is exactly what makes non-exact one-step methods possible
at all. Chern–Weil theory is the theory of what a connection buys; that is the honest point of
contact, and it is about §8's connection-dependence rather than about a stronger closure axiom.

**What the Galois analogy gives, and what it does not.** The shape is right — a correspondence
between classes of methods and $\otimes$-closed classes of algebras, with the monoid structure in the
role of the subgroup lattice — and it is what organises the results above. Two things it does not
give. There is no group: $\mathcal W$ is a monoid of *objects*, not of automorphisms, and
$\operatorname{Aut}(A)$ cannot play that role because base change preserves the lifted locus (§11's
second trap), so it acts trivially on everything in sight. And the algebras are not differential
rings, so no Picard–Vessiot statement is available; the genuine import from that theory is the
functor-of-points move recorded in §7 above.

---

## 8. Naturality under all diffeomorphisms collapses to the flow

> **Theorem.** If $\Psi$ is local and natural under all local diffeomorphisms, then
> $\Psi^X_h=\varphi^X_{c(h)}$ for a smooth $c$ with $c(0)=0$ (consistency adds $c'(0)=1$).

*Proof sketch.* Taking $\phi=\varphi^X_s$ shows $\Psi^X_h$ commutes with the flow of $X$. On
$\{X\neq0\}$ straighten to $X=\partial_1$ by the flow-box theorem. Naturality under
$x_1$-translations makes the displacement $x_1$-independent; naturality under all diffeomorphisms
of the transverse slice forces the transverse component to be a natural vector field, hence $0$,
and the $\partial_1$-component to be a natural scalar, hence constant. So $\Psi^{\partial_1}_h$ is
translation by $c(h)e_1$. Locality and density of $\{X\neq0\}$ extend this to all $X$. $\square$

**The conjecture "$c(h)=ch$" is false**: $\Psi^X_h=\varphi^X_{\sin h}$ is diffeomorphism-natural and
consistent. Linearity follows only if one additionally imposes the one-parameter-group law
$\Psi_{h_1}\circ\Psi_{h_2}=\Psi_{h_1+h_2}$. The reason no genuine one-step method survives is that
$\tfrac12X'X$ is **not a natural vector field** (it needs a connection); by Kolář–Michor–Slovák the
natural operators $\mathfrak X\rightsquigarrow\mathfrak X$ are only $X\mapsto cX$. Hence: *the affine
category is exactly what makes non-exact one-step methods possible at all.*

---

## 9. The cotangent asymmetry: why forward mode is free and reverse mode is not

The cotangent lift of a field is the Hamiltonian field of $H(x,\lambda)=\langle\lambda,X(x)\rangle$,
i.e. $\dot x=X(x)$, $\dot\lambda=-DX(x)^{\!\top}\lambda$; of a map,
$T^*\psi(x,\lambda)=(\psi(x),(D\psi(x))^{-\top}\lambda)$. The analogue of (T) is
$\Psi^{X^{T^*}}_h=T^*(\Psi^X_h)$.

**On linear fields this is exactly $R(z)R(-z)=1$** for the stability function $R$. Verified:

| method | $R(z)R(-z)-1$ | (T*) |
|---|---|---|
| explicit Euler | $-z^2$ | ✗ |
| Heun | $z^4/4$ | ✗ |
| RK3 | $-z^6/36$ | ✗ |
| RK4 | $z^6/72+z^8/576$ | ✗ |
| implicit midpoint | $0$ | ✓ |

These are instances of a theorem, not a coincidence. For polynomial $R$, $R(z)R(-z)=1$ forces
$\deg R=0$ by a degree count — equivalently the top coefficient of $R(z)R(-z)$ is
$(-1)^sc_s^2\neq0$ — so **no consistent explicit Runge–Kutta method can ever satisfy cotangent
closure.**

In general the discrete adjoint of RK $(A,b)$ with $b_i\neq0$ is RK with
$\hat a_{ij}=b_j-b_ja_{ji}/b_i$, $\hat b=b$ (Hager; Sanz-Serna), and closure for a *single*
method requires
$$b_ia_{ij}+b_ja_{ji}=b_ib_j,$$
i.e. **symplecticity**. (The linear test is strictly weaker: the trapezoidal rule passes
$R(z)R(-z)=1$ but is not symplectic and fails (T*) on nonlinear fields.)

**In AD terms.** (T) is **forward mode** — evaluation over $D$, which is what §2 says — and (T*) is
**reverse mode**: the $\lambda$-component the method returns on the cotangent-lifted field *is* the
discrete adjoint of its own step, so (T*) is the statement that differentiate-then-discretise and
discretise-then-differentiate agree on the reverse sweep. Read that way the recorded result is the
familiar one: the adjoint of an RK method is again RK (Hager; Sanz-Serna), and it is the *same*
method exactly when the method is symplectic.

**The shape transports differently for $T$ and for $T^*$** *(verified)*. §7 records that the property
is stated relative to a shape which itself carries a base change; leapfrog's shape is a splitting of
the coordinates. Under $T$ the splitting travels blockwise — coordinate $j$ becomes the $A$-block
$\{jN,\dots\}$, and that matching alone makes leapfrog closed. Under $T^*$ it does **not**. On a
nonlinear separable field $X(q,p)=(f(p),g(q))$, leapfrog satisfies (T*) with the **crossed** matching
$\{q,\lambda_p\}\mid\{p,\lambda_q\}$, and fails with the obvious cotangent-lifted one
$\{q,\lambda_q\}\mid\{p,\lambda_p\}$ and with base-versus-fibre. The reason is visible in the lifted
Hamiltonian:
$$\langle\lambda,X(x)\rangle=\lambda_q f(p)+\lambda_p g(q),$$
whose two terms involve $(p,\lambda_q)$ and $(q,\lambda_p)$ — so it is separable exactly for the
crossed grouping, which is also the grouping that keeps the method explicit. Leapfrog splits a
Hamiltonian into Lagrangian halves, and the adjoint of a position behaves like a momentum. So the
shape is genuinely part of the data and its transport is *dictated by the lift*, not chosen: $T$
carries a splitting blockwise, $T^*$ crosses it.

Two consequences worth recording. **Leapfrog is (T)-closed, (T*)-closed and not affine equivariant** —
both AD modes coexist with non-equivariance, which is the sharpest form of the point that the
equivariance framing is not the only lens. And (T*) has a **prerequisite (T) does not**: since $T^*$
is functorial only on the groupoid of diffeomorphisms, (T*) is statable only where $\Psi^X_h$ is
invertible. Over $A$ that costs nothing extra — by the unit criterion of §7 a matrix over $A$ is
invertible iff its real part is (verified) — so once (T*) is statable over $\mathbb R$ it is statable
over every $A$.

**The reason** is *not* "$T^*$ is not a functor" — $T^*$ is a functor on the groupoid of
diffeomorphisms, and the cotangent lift of a vector field is perfectly natural. The operative
point is that $T=-\otimes_\mathbb R D$ is a **covariant base change**, so the Picard and stage
equations base-change and uniqueness transfers; whereas $T^*$ uses the **dual module** and the
contragredient $(D\psi)^{-\top}$. Inversion and transposition are not base-change operations, and
the discrete transpose of an RK method is a *different* method. In one line:
**(T) is a tensor property and is free; (T\*) is a duality property and costs a condition.**

---

## 10. Summary

1. **(T) $\iff$ commutes with base change $\mathbb R\rightsquigarrow D$**, on the lifted locus, along the doubling tower only.
2. **Sufficiency:** anything that is the unique solution of an equation natural in the algebra base is (T)-natural — the flow and RK/PRK for the *same* reason. Necessity **fails**.
3. **The obstruction is not the trace.** Loops die (factor $\dim_\mathbb R A$), stolons die (non-balancedness), **lianas survive** — and only over $D$.
4. **(T) alone is weak**: continuum-many (T)-natural methods; it does not relate $n$ to $n+1$. Useful characterisations are relative to an equivariance group scheme.
5. **(T) is not a jet condition**, unlike everything in the classification literature.
6. **Beyond $D$:** $T^m$ free, general Weil $A$ not.
7. **Diffeomorphism-naturality** collapses to $\varphi^X_{c(h)}$.
8. **Reverse mode costs symplecticity**, for duality rather than functoriality reasons.

## 11. What remains open

* **Definability is impossible — settled negatively.** One might hope to characterise (T)-natural methods intrinsically, presupposing neither an equivariance class nor a series. The extension theorem of §5 rules this out in **every** category — smooth, analytic *or* algebraic — since the counterexamples are already polynomial. What is missing is not regularity but a **dimension-uniformity axiom**; the equivariance hypothesis is not a removable technical convenience, it is exactly the missing content.   The constructive replacement is to stop privileging $D$. Call $\Psi$ **algebraically natural**
  if for *every* finite-dimensional commutative $\mathbb R$-algebra $A$, $\Psi_{n\dim A}$ restricted to
  $A$-lifted fields is the $\mathbb R$-realization of the $A$-base change of $\Psi_n$. By
  Kolář–Michor–Slovák this is exactly naturality with respect to **all product-preserving
  endofunctors of $\mathbf{Mf}$** — intrinsic, series-free and equivariance-free — and the three
  algebra types kill the three obstructions separately:

  | algebra | condition it imposes | what it kills |
  |---|---|---|
  | $A=\mathbb R^k$ | decoupling $\varphi(f_1\oplus f_2)=\varphi(f_1)\oplus\varphi(f_2)$ | multi-aromas (disconnectedness) |
  | $A=D$ | (T) | loops and stolons |
  | $A=\mathbb R[\varepsilon]/(\varepsilon^3)$ | $c_A=1+\varepsilon^2\neq1$ | lianas |

  **What it forces on a single dimension — answered, and the answer is small (§7).** A collision
  between two based algebras of equal dimension forces $\Psi^X_h$ to be affine whenever $X$ is, and
  the collision trichotomy shows *nothing more is available by that route*: the obstruction between
  two distinct based algebras is always second order, so no collision ever reaches a non-affine
  field. Together with the extension theorem of §5 this is the honest answer to "what does closure
  provide unaided": rigidity on the affine locus, and off it nothing that this mechanism can see.

  **Conjecture.** For local, $f$-analytic families: algebraically natural $\iff$ affine
  equivariant $\iff$ B-series. Here (ii)$\iff$(iii) is MMMV and (iii)$\Rightarrow$(i) is immediate
  (trees are contractions over any base); (i)$\Rightarrow$(iii) is the open half. The honest doubt:
  base changes produce only *block* maps $\phi\otimes\operatorname{id}_{\mathbb R^n}$, never mixing the
  $\mathbb R^n$ directions, so they may not generate $GL(n)$-equivariance on their own — and the
  partitioned case shows the boundary is delicate.
* **The $O(n)$ row.** "(T) + orthogonal equivariance + locality + trivial decoupling $\Rightarrow$ exotic B-series" holds one-directionally at jet level. Two corrections to the obvious attack, both verified:

  *The $\Theta$-operator does not act by $2^{\ell}$.* One is tempted to say that
  $\Theta=\pi\circ(\cdot)\circ\zeta$ satisfies $\Theta(F(\gamma))=2^{\ell(\gamma)}F(\gamma)$ and conclude
  $b(\gamma)=0$ for every $\gamma$ carrying a loop. This is **false in the exotic class**:
  $\lvert DX\rvert_F^2X$ has *no loop at all* yet obstructs at base $v$-degree $0$ with multiplier $2$.
  $\Theta$ acts by $N_0(\gamma)$, the number of zero-cost $\varepsilon$-level assignments, which is not a
  power of two attached to any graph statistic. The $v$-degree-$0$ equation is therefore
  $\sum_\gamma b(\gamma)\bigl(N_0(\gamma)-1\bigr)F(\gamma)=0$, killing every $\gamma$ with $N_0\neq1$ — strictly
  more than the loop-carrying ones. That still drives the theorem (exotic trees have $N_0=1$;
  $\Delta X$ has identically zero defect), but the surviving class must be characterised as
  "multiplier $1$", never as "loop-free".

  *The $v$-degree separation is real, but only in the base block.* Expanding the defect in
  $v=t v_0$:

  | differential | structure | base block | fibre block |
  |---|---|---|---|
  | $\operatorname{div}(X)X$ | one loop | $[0]$ | $[1]$ |
  | $\Delta X$ | liana only | zero | zero |
  | $\langle X,X\rangle X$ | stolon, no loop | $[2]$ | $[1,3]$ |
  | $\lvert DX\rvert_F^2X$ | liana+stolon, no loop | $[0,2]$ | $[1,3]$ |

  Pure loops obstruct at base $v$-degree $0$ and pure stolons only at base $v$-degree $2$, so those
  two cannot cancel and independence is needed only *within* a fixed $v$-degree — weaker than a
  full independence lemma. The claim that stolons have "minimal $v$-degree exactly $2$" is true of
  the base block only: the fibre block carries a degree-$1$ defect, $-2\langle X,DXv\rangle X$, from
  $D(\lvert X\rvert^2X)v=2\langle X,DXv\rangle X+\lvert X\rvert^2DXv$.

  The remaining gap was the *stolon* half: since (T) lives on the thin locus $\{TX\}$, the obstruction
  of a stolon-carrying forest with $N_0=1$ appears only at base $v$-degree 2, so an independence
  statement there — **Lemma L** — is unavoidable.

  **Lemma L is now proved** *(single source; audit in progress)*, and not by the route anyone
  proposed. No dual vector fields and no $\theta$-parameters are needed. Set $x=0$ and pair with
  $u$; the jet tensors $c_k=\partial^kX(0)$ are freely prescribable by a polynomial field, so the
  relation becomes a polynomial identity and multihomogeneity separates valence profile from
  $v$-degree. One profile then gives a slot set $H$ of size $2N$,
  $N=1+d+\#\text{arrows}+\#\text{lianas}+\#\text{stolons}$ (upper slots $m=1+a+2s$, lower slots
  $\sum_tk_t=a+2l+d$), and each elementary differential becomes a perfect matching $\Phi_M$ of $H$.
  The hinge is that
  $$\tilde G=\Bigl(\textstyle\prod_k\bigl(S_{m_k}\ltimes(S_k)^{m_k}\bigr)\Bigr)\times S_d$$
  has orbits *exactly* the isomorphism classes of $v$-decorated forests — neither too small (which
  would split one forest across orbits and manufacture a false relation) nor too large (which would
  merge two). The jet data span $(V^{\otimes H})^{\tilde G}$ by polarisation, $\Phi_M$ restricted
  there is the $\tilde G$-orbit sum, distinct forests give disjoint orbits and hence disjoint
  supports, and matchings are linearly independent once $n\ge N$. So every coefficient vanishes.

  The one external ingredient is independence of Brauer matchings for $n\ge\lvert H\rvert/2$ — the
  **easy** half of $O(n)$ invariant theory; the argument never uses that contractions *span*. The
  result is more general than needed: any $v$-degree, any number of colours, no restriction to the
  residual class.

  **Audited independently; it stands, with corrections.** The bound is on $N=m-s+l+d$, *not* on the
  node count — the theta graph has $m=3,N=4$ while $\langle X,X\rangle X$ has $m=3,N=2$; $n\ge N$ is
  probably not sharp ($n\ge\lvert V\rvert+1$ is conjectured). The audit also found the first rank
  drop in this project: the profile $(2,2,2)$ at $d=0$ has 16 iso-classes but rank only 15 in
  dimension 3, becoming full at $n=4$ — so the dimension hypothesis genuinely bites rather than
  being a formality. Theoretical and empirical ranks agree even where deficient, sharp evidence that
  the polarisation step loses nothing. Step 1's correspondence was checked on 1690 matchings and
  Step 3's invariant subspace against Burnside in 12 cases.

  One methodological correction, now a standing trap in `CLAUDE.md`: the natural test of whether
  $\tilde G$ is *too large* — comparing its orbits against equality classes of the computed
  differentials — is **blind**, since $\tilde G$ fixes the jet tensor and so orbit-equality implies
  value-equality for any subgroup at all. The hinge was instead confirmed by enumerating iso-classes
  independently of $\tilde G$, on 11 valence profiles.

  This also sharpens the standing trap about the shortcut $v:=X(x)$: it enlarges $S_{m_0}\times S_d$
  to $S_{m_0+d}$ and fuses orbits **iff** the forest already has a derivative-free node, so it is
  invisible on profiles with $m_0=0$ — which is why it is easy to miss.
* **The $GL(n_1)\times GL(n_2)$ row.** The ambient aromatic-P-series classification is open in print.
* **Substitution.** $\delta_N$ is an automorphism of the aromatic Butcher *composition* group; its compatibility with the *substitution* law is unverified.

### Round-4 status (two independent passes, converging)

Two independent agents attacked the open items; they agree on every substantive point below, and
the starred items were re-verified here. Items still resting on a single pass are marked.

* **The multiplier law, exactly.** $N_0(\gamma)$ counts families of vertex-disjoint generalised
  cycles (§4). $N_0=2^{\#\text{aromas}}$ on stolon-free forests; $\mathbf{N_0=3}$ for the theta graph
  $\lvert D^2X\rvert_F^2X$ ★. The $v$-degree-$0$ equation $\sum_\gamma b(\gamma)(N_0(\gamma)-1)F(\gamma)=0$
  therefore kills $\gamma$ iff it carries an aroma **or** a liana–stolon alternating cycle.
* **Defect characterisation.** For a single elementary differential the defect vanishes
  identically **iff** $\gamma$ has no stolon and no directed arrow-cycle — i.e. iff $\gamma$ is an
  exotic tree. Verified exhaustively on all 84 iso-classes with $\le4$ nodes. This matches
  Laurent–Munthe-Kaas Thm 2.13's right-orthogonal class exactly, so the $O(n)$ target is precisely
  *(T) $\Rightarrow$ right-orthogonal equivariance*.
* **Parity, with a clean proof.** Base block carries only even $v$-degrees, fibre only odd —
  because $\varphi=(\operatorname{id},-\operatorname{id})\in O(2n)$ fixes $TX$ pointwise, so $O(2n)$-equivariance
  forces the parity. (The level-counting proof of Appendix A gives the same.) ★ consistent with the
  measured table.
* **Route (b) does not close; route (a) is required — but is triangular.** Trivial decoupling kills
  disconnected forests, and $v$-degree $0$ kills $N_0\neq1$. What survives is non-empty: *connected
  exotic aromatic forests carrying a stolon with $N_0=1$*, minimal member (3 nodes)
  $\partial_{jk}X^i\,\partial_jX^a\,\partial_kX^a$. For these the base defect starts at $v$-degree 2 and the
  fibre $v$-degree-1 defect vanishes identically, so no equation among undecorated differentials
  sees them. Independence at $v$-degree 2 is unavoidable. It is however *triangular*: erasing the
  $v$-leaves recovers $\gamma$, so distinct $\gamma$ contribute disjoint families and no cross-$\gamma$
  cancellation is possible.
* **Lemma L, reduced and tested.** Attaching a fresh $Y$-leaf at each $v$-slot is injective onto
  bicoloured forests whose $Y$-nodes are sourceless derivative-free leaves, so Lemma L follows from
  the two-colour version of Laurent–Munthe-Kaas Prop. 4.1 with $Y$ constant. The proof should be
  theirs with one extra colour and the $\theta$-parameters retained — their Prop. 4.1 explicitly
  notes the $\theta$-free dual field is *insufficient*. Conclusion verified on the whole residual
  class up to 4 nodes (full rank). *Superseded: Lemma L is now proved outright — see §11.*
  The two enumerations quoted in this project (62 forests at $\le3$ nodes; 84 iso-classes at $\le4$
  nodes and $\le4$ slots) have been reconciled: they are different truncations and **neither
  contains the other** (46 shared, 16 and 38 exclusive). Contrary to an earlier guess here, no
  connectivity filter distinguishes them.
* **C2: the (T) half is proved.** With no metric there are no lianas or stolons, so $s_1=0$ forces
  base $v$-degree $\equiv0$: a single equation, killed by bicoloured aromatic independence, leaving
  exactly the bicoloured trees. The fibre $v$-degree-1 defect vanishes because in an aroma-free
  forest each node has a unique path to the root. **No Lemma L needed.** Only the ambient step
  (transcribing LMK §3 with the $GL(V_1)\times GL(V_2)$ FFT) remains open.
* **C3, the coherence subtlety — settled, but by a different mechanism than proposed** ★. For *any*
  finite-dimensional commutative $A$ and *any* basis, a **linear** field $M$ realises as
  $M^A=I_N\otimes M$. Hence the $D$-lifted and $\mathbb R^2$-lifted loci coincide exactly on linear
  fields, and demanding both prescriptions forces
  $$D\Psi^M(x)\,v=\Psi^M(v)\quad\text{for all }x,v,$$
  i.e. $\Psi^M$ is **linear** for every linear field. Verified: Euler ✓, RK4 ✓, and the
  method $u+hX+h^2(X^1)^2X$ ✗. *(Correction: that formula, with the scalar left real-valued, is
  **not** (T)-natural — see §7. The R4 conclusion is unaffected and now has a genuine witness, the
  ring-valued tower method, which is (T)-natural and also fails linear rigidity.)* So the killing
  factorisation is the pair
  $(2,D)$ vs $(2,\mathbb R^2)$ — a one-line identity — not the $n=1$ argument I had guessed.
* **Scalar rigidity.** If $\Psi_m=u+hX+h^2q_m(j^rX(u))X(u)+O(h^3)$ with $q_m$ **real-valued** and
  $\Psi$ is (T)-natural as a family indexed by *real dimension*, then every $q_m$ is constant
  ($\operatorname{div}X$, $\lvert X\rvert^2$, $(X^1)^2$ all fail; $q\equiv1$ passes). A real-valued
  contraction multiplies by a real scalar, whereas base change demands an $A$-scalar.
* **A definitional horn dilemma.** Quantifying algebraic naturality over *all* $\mathbb R$-bases $\iota$
  is **too strong**: taking $A=\mathbb R$ makes $\iota\in GL(n)$ arbitrary and the condition *states*
  $GL(n)$-equivariance by fiat, trivialising the conjecture. Fixing one $\iota$ per $A$ is the right
  definition; the resulting basis-relativity is then $GL$-conjugation, hence invisible for
  equivariant $\Psi$ and visible only where the conjecture lives. *(This corrects the instruction I
  gave — "quantify over all bases" — which would have been vacuous.)*
* **Two scope corrections.** The legitimate class is Weil algebras **and their finite products**,
  not all finite-dimensional commutative $\mathbb R$-algebras ($\mathbb C$ gives no product-preserving
  functor on real manifolds). And $A=\mathbb R^k$ yields only *diagonal* decoupling, not
  $\Psi^{X\oplus Y}=\Psi^X\oplus\Psi^Y$ — though that still kills every forest with a non-root
  component.
* **Coherence does not recover $GL(n)$-equivariance.** Both passes say no. They locate it
  differently — one finds the group large enough ($\mathfrak{gl}(nN)$ is generated) but the lifted
  locus too thin — *this half is corrected: the automorphism-generated Lie algebra is exactly
  $\mathfrak{gl}(m-1)$ and fixes the base pointwise, so it is the wrong group rather than a large
  group on a small set; see §7*; the other observes that coherence is not a group action at all, being a system of
  restriction conditions on overlapping thin loci. The disagreement traces to the horn dilemma
  above: the "all bases" reading supplies the extra generators but trivialises the question.
* **A trap worth recording.** Testing on a *quadratic* field makes $\Delta X$ constant and produces
  a spurious exception; and taking $v:=X(x)$ as a shortcut in Lemma L destroys injectivity and
  manufactures exactly the false relations one is trying to rule out. Forests must be canonicalised
  up to renaming of internal indices before independence is tested.

## 12. Verified computations

All exact over $\mathbb Q$; run `python3 verify/run_all.py`.

| file | establishes |
|---|---|
| `verify/dual.py` | dual numbers $D$ and jets $\mathbb R[\varepsilon]/(\varepsilon^3)$ over `Fraction` |
| `verify/test_basic.py` | RK4, Heun, Taylor-2 satisfy (T); aromatic and basis-dependent methods fail |
| `verify/test_aromas.py` | aromas scale by $2^{\alpha}$ (ratios 2, 4, 8); trees are natural |
| `verify/test_partitioned.py` | PRK satisfies (T) for the $D$-block lift, fails for the real lifts; not affine equivariant |
| `verify/test_weil.py` | $\operatorname{div}(T^AX)=3\operatorname{div}X$ for $\dim_\mathbb R A=3$; RK4 is $T^A$-natural |
| `verify/algebra.py` | based algebras (Weil and products), base change over any of them, recognition of lifts |
| `verify/test_jettransport.py` | leapfrog is closed for the lifted splitting and for neither other matching; it is not affine equivariant; invertibility in $A$ is exactly unit-hood |
| `verify/test_aut.py` | $\operatorname{Der}(A)$ for the based algebras; every derivation kills the unit; the $\operatorname{Aut}$-generated Lie algebra is $\mathfrak{gl}(m-1)$ at $m=4,6$; the $J^1_k$-lift is the $k$-fold Whitney sum |
| `verify/test_order2.py` | the pointwise commutant reads $\rho_A(A)$ without collapsing on the real locus; multiplicativity; separation of the factorizations; the per-algebra-free lift formula on polynomial data |
| `verify/test_relations.py` | graded Weil algebras act through their even parts, which are ordinary ones; failure of quotient- and product-closure with explicit witnesses; $\operatorname{Aut}(A)$-equivariance on the lifted locus and the $GL(k)$ it produces |
| `verify/test_presentation.py` | the commutant recovers the algebra on 17 based algebras and is the common refinement in both coincidence families; zero-field rigidity; degeneration of the jet loci on the real locus |
| `verify/test_spectrum.py` | the closure spectrum: monoidality under $\otimes$, the corrected survival law on 17 based algebras, the stolon-defect identity, failure of subalgebra-closure, and the incomparable tower spectra |
| `verify/test_affine.py` | the collision trichotomy on all 41 pairs from 16 based algebras; affine rigidity and its two witness pairs; its separation from linear rigidity; functoriality of the cross-dimensional collisions |
| `verify/test_contractions.py` | liana (Laplacian) is (T)-natural and not affine equivariant; stolon fails; liana defect over $\varepsilon^3$ equals $\Delta X$ |
| `verify/test_cotangent.py` | $R(z)R(-z)$ table; $(-1)^sc_s^2\neq0$; trapezoidal rule passes the linear test but is not symplectic |
| `verify/test_vdegree.py` | $v$-degree of the defect, split base/fibre; refutes $\Theta=2^{\ell}$; loop/stolon separation at base degree $0$ vs $2$ |
| `verify/test_exotic.py` | multiplier counts loops not aromas ($\lvert DX\rvert_F^2$: no loop, $\times2$; $\langle X,X\rangle$: loopless aroma, $\times1$); the $c_A$ law for $r=2,3,4$; a **polynomial** closed method with no equivariance; gradient fields not stable under $T$ |

## Appendix A. The index calculus — intuition, and what it cannot see

The multipliers of §4 can be *derived* by index bookkeeping. This is worth recording because it
reproduces every verified number from one rule; it is worth fencing off because it is a
coordinate shadow of the algebra, not the reason.

**The rule.** On $T^A\mathbb R^n=A^n$ write a real index as $I=(i,\alpha)$ with $\alpha$ the
$\varepsilon$-level, and let $\{e_\alpha\}$ be an $\mathbb R$-basis of $A$ declared orthonormal, $\{e^\alpha\}$
its dual. For a base-changed field the node tensor is
$$f^{(i,\alpha)}_{(j_1\beta_1)\cdots(j_k\beta_k)}
   =e^{\alpha}\!\left(\partial_{j_1\cdots j_k}\tilde X^{i}\cdot e_{\beta_1}\cdots e_{\beta_k}\right),$$
because $\partial/\partial x^j=\partial_{z^j}$ while $\partial/\partial v^j=\varepsilon\,\partial_{z^j}$: a lower
index at level $\beta$ carries $e_\beta$, an upper index at level $\alpha$ reads off $e^\alpha$.
Summing each contracted pair over its level gives four outcomes, and only four:

| pattern | sum over levels | value | over $D$ |
|---|---|---|---|
| **arrow** (upper·lower, distinct legs) | $\sum_\alpha e^\alpha(y)\,e_\alpha$ | $y$ — identity on $A$ | exact |
| **loop** (upper·lower, same cycle) | $\sum_\alpha e^\alpha(e_\alpha y)$ | $\operatorname{tr}_{A/\mathbb R}(y)=N\operatorname{aug}(y)$ | $2\operatorname{aug}(y)$ |
| **liana** $\delta^{jk}$ (lower·lower) | $\sum_\alpha e_\alpha^2$ | $c_A$ | $1+\varepsilon^2=1$ |
| **stolon** $\delta_{jk}$ (upper·upper) | $\sum_\alpha e^\alpha\otimes e^\alpha$ | $g$ | Euclidean form |

The loop entry also explains why a loop evaluated on a lifted field is *fibre-free*: $e^\alpha(e_\alpha y)$
destroys the $\varepsilon$-part rather than merely rescaling it.

**Level bookkeeping and a parity law.** Write $k_v$ for the upper level of a node and
$K_v=\sum(\text{incoming arrow levels})$. Since $\partial^kX(z)=\partial^kX(x)+\varepsilon\,\partial_m\partial^kX(x)v^m$,
the node contributes $0$ when $K_v\ge2$; $\partial^qX$ when $K_v=k_v$; and $\partial_m\partial^qX\,v^m$
(v-degree $+1$) exactly in the state $(K_v,k_v)=(0,1)$. Counting arrow sources twice,
$$\textstyle\sum_vK_v=\sum_{\rm arrows}k_{\rm src}+2l_1,\qquad
  \sum_vk_v=\sum_{\rm arrows}k_{\rm src}+2s_1+k_{\rm out},$$
so $\sum_v(k_v-K_v)=2(s_1-l_1)+k_{\rm out}$ with $l_1,s_1$ the numbers of level-1 lianas and
stolons. On non-vanishing terms $k_v-K_v\in\{0,1\}$, giving
$$\textbf{v-degree}=2(s_1-l_1)+[k_{\rm out}=1].$$
**Corollary: the base block carries only even v-degrees and the fibre block only odd ones** — which
is exactly the pattern of the measured table in §11 ($[0]$ vs $[1]$; $[2]$ vs $[1,3]$; $[0,2]$ vs
$[1,3]$).

**Three things the calculus cannot see.**

1. *It presupposes the classification it appears to explain.* Which contraction patterns exist at
   all is the invariant-tensor theorem plus Peetre — i.e. Munthe-Kaas–Verdier and
   Laurent–Munthe-Kaas. Fix $H=GL$ and only upper–lower pairings exist, so "loops are the only
   obstruction" looks like a fact about differentiation; fix $H=O(n)$ and $\delta^{jk},\delta_{jk}$
   appear, and the same rule yields a different answer. The calculus computes *inside* an ambient
   class; it never produces one. This is precisely why the characterisation is relative (§5).
2. *It is jet-level.* Every index expression is a statement about a Taylor coefficient, whereas
   (T) is an exact identity at each fixed $h$; the flat method of §6 has no index expression
   distinguishing it from Euler.
3. *It is per-dimension*, and that is where (T) is actually weak. By the extension theorem of §5,
   (T) constrains no individual $\Psi_d$ at all. No monomial calculus can detect this, since
   indices live in a fixed dimension while the phenomenon is a relation between $\Psi_n$ and
   $\Psi_{2n}$.

**What the table really says.** $A^n$ viewed over $\mathbb R$ carries strictly more structure than it
carries over $A$: $\operatorname{End}_A(A^n)\subsetneq\operatorname{End}_\mathbb R(A^n)$, and an $\mathbb R$-bilinear form on
$A^n$ need not be $A$-balanced. The four entries are four invariants of the pair (algebra,
metric) — the identity, the trace form $\operatorname{tr}_{A/\mathbb R}$, the element $c_A=\mu(g^{-1})$, and $g$
itself. The graph vocabulary merely presents *which $\mathbb R$-multilinear invariants happen to be
base changes of $A$-multilinear ones*. The calculus reports that self-traces and
metric-on-evaluations are forbidden; the reason they are forbidden is that base change is a
functor and those operations are not defined over the base ring.

## References

* R. I. McLachlan, K. Modin, H. Munthe-Kaas, O. Verdier, *B-series methods are exactly the affine equivariant methods*, Numer. Math. **133** (2016) 599–622. (Def. 2.1–2.3, Thm 2.4, Rmk 2.5.)
* H. Munthe-Kaas, O. Verdier, *Aromatic Butcher series*, Found. Comput. Math. **16** (2016) 183–215. (Thm 2.4; §7–8 for general $H$.)
* A. Laurent, H. Munthe-Kaas, *The universal equivariance properties of exotic aromatic B-series*, Found. Comput. Math. (2024); arXiv:2305.10993. (Def. 2.7 lianas/stolons; Thms 2.12–2.13; Prop. 4.1, 4.3; Table 2.)
* A. Bogfjellmo, *Algebraic structure of aromatic B-series*, arXiv:1505.01973.
* I. Kolář, P. W. Michor, J. Slovák, *Natural Operations in Differential Geometry*, Springer 1993. (Ch. VIII, Weil functors = product-preserving bundle functors, Thm 35.13.)
* K. Yano, S. Ishihara, *Tangent and Cotangent Bundles*, Marcel Dekker 1973. (Complete lift; $\varphi^{X^C}_t=T\varphi^X_t$.)
* W. W. Hager, *Runge–Kutta methods in optimal control and the transformed adjoint system*, Numer. Math. **87** (2000) 247–282.
* J. M. Sanz-Serna, *Symplectic Runge–Kutta schemes for adjoint equations, automatic differentiation, optimal control, and more*, SIAM Review **58** (2016) 3–33.
* E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd ed., Springer 2006.
