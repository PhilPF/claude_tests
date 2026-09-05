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
Nor can $c_A=1$ be restored by rescaling the metric: with $g=\operatorname{diag}(1,w_1,\dots)$ one gets
$c_A=1+w_1^{-1}\varepsilon^2+\cdots$, never $1$.

*(Terminology is that of Laurent–Munthe-Kaas, Def. 2.7: a **liana** identifies two **arrows** —
a double derivation, i.e. the Laplacian; a **stolon** identifies two **nodes** — a double
evaluation, i.e. the scalar product.)*

**Loops.** For $M\in\operatorname{End}_D(D^n)$ one has $\operatorname{tr}_\mathbb R=\operatorname{tr}_{D/\mathbb R}\circ\operatorname{tr}_D$ with
$\operatorname{tr}_{D/\mathbb R}(a+b\varepsilon)=2a$, because multiplication by $a+b\varepsilon$ has matrix
$\left(\begin{smallmatrix}a&0\\b&a\end{smallmatrix}\right)$. In general $\operatorname{tr}_{A/\mathbb R}(u)=(\dim_\mathbb R A)\operatorname{aug}(u)$,
since multiplication by a nilpotent has trace $0$. Hence the $\varepsilon$-part is **discarded** — which is why a loop evaluated on a lifted field is
fibre-free — and each loop scales by $\dim_\mathbb R A$. Verified: ratios $2,4,8$ for $A=D$ and ratio
$3$ for $\dim_\mathbb R A=3$.

**The multiplier counts loops, not aromas.** In the $GL$/aromatic world the two coincide: since
$\sigma$ is a fixed-point-free involution the root's upper index is consumed by the ghost arrow,
so every non-root component carries exactly one cycle, and a forest with $\alpha$ aromas scales by
$(\dim_\mathbb R A)^{\alpha}$. In the exotic ($O(n)$) world they **diverge**: $\langle X,X\rangle$ is a
loopless aroma, and $\lvert DX\rvert_F^2\,X$ — a liana feeding a stolon, with *no loop at all* —
still has multiplier exactly $2$. So the grading automorphism $\delta_N$ and the slogan
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
splitting/IMEX methods are (T)-natural and are *not* P-series — nor are they
$GL(n_1)\times GL(n_2)$-equivariant, so there is no contradiction.

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
> for $A=\mathbb R[\varepsilon]/(\varepsilon^3)$: the defect is concentrated in the $\varepsilon^2$ block and equals
> exactly $\Delta X$. (Verified.)

The mechanism is the liana constant of §4: $c_{D^{\otimes m}}=1$, so lianas survive every $T^m$,
whereas $c_{\mathbb R[\varepsilon]/(\varepsilon^3)}=1+\varepsilon^2$, so the liana realizes $(1+\varepsilon^2)\Delta_A$ instead
of $\Delta_A$ and the defect is $\varepsilon^2\,T^A(\Delta X)$ — exactly the verified $\varepsilon^2$-block
discrepancy. This *derives* the failure rather than asserting it, and it is the operative
reason; the observation that $\mathbb R[\varepsilon]/(\varepsilon^{r+1})\hookrightarrow D^{\otimes r}$ is a subalgebra
rather than a quotient is a corollary, not the cause.

*Positive counterpart.* (T) together with equivariance under linear permutation matrices gives
$T^{(r)}$-naturality for all $r$ — every higher variational equation — because
$\mathbb R[s]/(s^{r+1})\subset D^{\otimes r}$ is the $S_r$-invariant subalgebra.

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

  **Conjecture.** For local, $f$-analytic families: algebraically natural $\iff$ affine
  equivariant $\iff$ B-series. Here (ii)$\iff$(iii) is MMMV and (iii)$\Rightarrow$(i) is immediate
  (trees are contractions over any base); (i)$\Rightarrow$(iii) is the open half. The honest doubt:
  base changes produce only *block* maps $\phi\otimes\operatorname{id}_{\mathbb R^n}$, never mixing the
  $\mathbb R^n$ directions, so they may not generate $GL(n)$-equivariance on their own — and the
  partitioned case shows the boundary is delicate.
* **The $O(n)$ row.** "(T) + orthogonal equivariance + locality + trivial decoupling $\Rightarrow$ exotic B-series" holds one-directionally at jet level. The *loop* half is gap-free: the operator $\Theta=\pi\circ(\cdot)\circ\zeta$ gives $\Theta(F(\gamma))=2^{\ell(\gamma)}F(\gamma)$ for all $\gamma$, and independence of exotic aromatic elementary differentials (Laurent–Munthe-Kaas Prop. 4.1) forces $b(\gamma)=0$ whenever $\gamma$ has a loop. The *stolon* half has one isolated gap: since (T) lives on the thin locus $\{TX\}$, one restricts along $s_W(x)=(x,W(x))$ and needs the **two-coloured** analogue of Prop. 4.1. Their $\theta$-parametrised dual vector fields should carry over, but this is not in the literature.
* **The $GL(n_1)\times GL(n_2)$ row.** The ambient aromatic-P-series classification is open in print.
* **Substitution.** $\delta_N$ is an automorphism of the aromatic Butcher *composition* group; its compatibility with the *substitution* law is unverified.

## 12. Verified computations

All exact over $\mathbb Q$; run `python3 verify/run_all.py`.

| file | establishes |
|---|---|
| `verify/dual.py` | dual numbers $D$ and jets $\mathbb R[\varepsilon]/(\varepsilon^3)$ over `Fraction` |
| `verify/test_basic.py` | RK4, Heun, Taylor-2 satisfy (T); aromatic and basis-dependent methods fail |
| `verify/test_aromas.py` | aromas scale by $2^{\alpha}$ (ratios 2, 4, 8); trees are natural |
| `verify/test_partitioned.py` | PRK satisfies (T) for the $D$-block lift, fails for the real lifts; not affine equivariant |
| `verify/test_weil.py` | $\operatorname{div}(T^AX)=3\operatorname{div}X$ for $\dim_\mathbb R A=3$; RK4 is $T^A$-natural |
| `verify/test_contractions.py` | liana (Laplacian) is (T)-natural and not affine equivariant; stolon fails; liana defect over $\varepsilon^3$ equals $\Delta X$ |
| `verify/test_cotangent.py` | $R(z)R(-z)$ table; $(-1)^sc_s^2\neq0$; trapezoidal rule passes the linear test but is not symplectic |

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
