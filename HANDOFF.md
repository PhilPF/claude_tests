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
$\mathbb R$-realization of the $A$-base change of $\Psi_n$. By Kolář–Michor–Slovák the legitimate class
is **Weil algebras and their finite products** — these are exactly the product-preserving functors on
manifolds; $\mathbb C$ is not one of them.

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

The open questions, in order of interest:

1. Is "algebraically natural $\Rightarrow$ affine equivariant $\Rightarrow$ B-series" true for
   $f$-analytic families? The first implication is the open half; conjectured **false** for smooth
   methods (1-jet loci are nowhere dense, $\dim n^2N<M^2$), blocked on gluing junk up the tower.
2. Construct the flat, coherent, non-equivariant algebraically natural method, or prove none exists.
3. *(Residue of the per-dimension question, now answered — §7.)* Cross-dimensional collisions
   $X^A=Y^B$ with $n_1\neq n_2$ were found to be functorial over a bounded range only. A proof needs:
   a lifted field determines its presentation up to common refinement.
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
