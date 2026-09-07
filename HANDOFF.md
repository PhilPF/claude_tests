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
* **Scalar rigidity.** A real-valued coefficient $q_m(j^rX)$ must be constant; base change demands an
  $A$-scalar and no real-valued contraction supplies one.
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

1. **Does algebraic naturality have per-dimension consequences beyond linearity on linear fields?**
   This is the real form of "what structure does closure provide, unaided".
2. Is "algebraically natural $\Rightarrow$ affine equivariant $\Rightarrow$ B-series" true for
   $f$-analytic families? The first implication is the open half; conjectured **false** for smooth
   methods (1-jet loci are nowhere dense, $\dim n^2N<M^2$), blocked on gluing junk up the tower.
3. Construct the flat, coherent, non-equivariant algebraically natural method, or prove none exists.

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
The $O(n)$ characterisation is complete — Lemma L proved (`notes/lemma-L-proof.md`) and independently
audited (`notes/lemma-L-audit.md`).

## Provenance

Corroborated by two or more independent passes plus checks here: the trichotomy, $N_0$, the parity
law, closure $\not\Rightarrow$ B-series, the non-jet property. **Single source + one audit:** the
Lemma L proof. **Single source:** linear rigidity and scalar rigidity (both spot-checked here).
Treat single-source items as leads when building on them.

## Keeping this current

Update this file whenever the target changes, a claim is retracted, or an open item closes — before
ending a session, not after. It is the only thing a fresh session reads to know where the work is.
