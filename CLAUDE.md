# Working notes for this repository

Research project on **flow closure under differentiation**: what the property
$T(\Psi^X_h)=\Psi^{TX}_h$ forces on an approximation of the flow.

> **Start by reading `HANDOFF.md`.** It carries the current target, what is settled, what is open,
> and the provenance of each claim, so that a fresh session can continue without re-deriving
> anything. **Keep it current**: update it whenever the target changes, a claim is retracted, or an
> open item closes — before ending a session, not after.

`flow-closure-under-differentiation.md` is the state of knowledge and the single source of truth; it
is kept *correct* rather than append-only — claims in it have been retracted several times when
verification contradicted them.

**Framing discipline.** Closure is a *per-tower* condition and imposes nothing on any individual
$\Psi_d$. Every characterisation in this area factors as
$[\text{equivariance}]\cap[\text{closure}]$, where equivariance produces the series and closure only
prunes it. Do not report such a result as though closure alone delivered it — that framing has
already drifted once.

## Standing practices

Cheap precautions. Following them costs almost nothing; their *scope* may still be wrong, so
widen or narrow them freely when a case demands it.

**Verify, don't trust.** Every computational claim in the write-up is checked in exact rational
arithmetic. `python3 verify/run_all.py` runs all of it (no dependencies; `sympy`/`numpy` are
unavailable and pip has no index). Lines printed `FAIL` are *expected* failures — methods that
must **not** satisfy closure; the suite passes iff every outcome matches its prediction.

**"Max degree $\le1$" is not "linear".** Testing linear rigidity by `max u-degree == 1` admits a
constant term — exactly what a distinguished vector contributes — and it passed a marked-vector
method that linear rigidity in fact kills, producing a wrong novelty claim. Require every monomial
to have $u$-degree exactly one.

**Use cubic or higher test fields.** With a quadratic field $\Delta X$ is constant and low-order
differentials degenerate. This produced spurious results twice.

**Canonicalise graphs before testing independence.** Two presentations of one forest under
renaming of internal indices manufacture a false linear relation. This happened twice.

**Distrust any statement that silently fixes an identification $\mathbb R^{n\dim A}\cong A^n$.**
Source of three separate errors: the partitioned-RK lift, $c_A$ over $\mathbb R[\varepsilon]/(\varepsilon^3)$,
and $c_A$ over the dual numbers themselves.

**A symmetry group cannot be validated against the values it acts on.** To check that a group $G$
used to classify contraction diagrams is not *too large*, comparing $G$-orbits with equality classes
of the computed differentials is **blind**: $G$ fixes the jet tensor $J$, so
$\Phi_{gM}(J)=\Phi_M(g^{-1}J)=\Phi_M(J)$ for *any* subgroup, and orbit-equality forces
value-equality automatically. Only an isomorphism-class count computed **independently of $G$** —
e.g. a slot-forgetting encoding canonicalised over all relabellings — decides it. A rank test is
equally blind in that direction.

## Findings that could be wrong — and how to overturn them

These are *empirical*, not axioms. Each carries its evidence and its falsifier. **Challenging them
is welcome and cheap**; a challenge just has to be evidenced — an exact computation, a source
quote, or an explicit counterexample. An unevidenced objection is not a challenge, and re-deriving
a settled result without new evidence is the waste this file exists to prevent.

| finding | evidence | what would overturn it |
|---|---|---|
| The multiplier is $N_0$, counting families of vertex-disjoint generalised cycles | two independent passes; five-case stress test here (`verify/test_n0.py`), incl. the decisive pair and a multiplicative mixed case | a forest whose measured base multiplier differs from its $N_0$ |
| Defect $\equiv0$ iff the forest is an exotic tree | exhaustive over 84 iso-classes, $\le4$ nodes | an exotic tree with non-zero defect, or a non-exotic forest with zero defect |
| Base block even $v$-degrees, fibre odd | two proofs (level counting; $(\mathrm{id},-\mathrm{id})\in O(2n)$ fixes $TX$) plus measurement | any measured defect of the wrong parity |
| Closure does not imply B-series | PRK and the Laplacian method, both verified | — (a counterexample stands) |
| Closure is not a jet condition | the flat method $\theta(h\operatorname{div}X)hX$ | — |
| Affine fields go to affine maps | two independent witness pairs — secant vs tangent in dim 2, and $(\mathbb R[\varepsilon]/\varepsilon^3,\mathbb R[x,y]/\mathfrak m^2)$ in dim 3 (`verify/test_affine.py`) | an algebraically natural method and an affine $X$ whose $\Psi^X_h$ is not affine |
| Collisions of two lifts cap at affine | proof (distinct based algebras differ at a product of two basis vectors) + all 41 pairs from 16 based algebras | based algebras of equal dimension whose lifts agree on a non-affine map |
| Loops *and stolons* survive only for $A=\mathbb R$; lianas iff $c_A=1$ | 17 based algebras $\times$ 6 methods, plus the exact defect identity $h^2(r\cdot1_A-q)X^A$ (`verify/test_spectrum.py`) | a based algebra of dimension $>1$ over which a loop- or stolon-carrying differential is $T^A$-natural |
| The closure spectrum is a $\otimes$-monoid, not subalgebra-closed | one-line proof plus multiplicativity of $\dim$, $c_A$, balancedness; Laplacian over $D\otimes D$ vs its $S_2$-invariant | a method natural over $A$ and over $B$ but not over $A\otimes B$ |
| The commutant recovers the presentation | $\mathcal A(X^A)=\rho_A(A)$ on 17 based algebras, and is the common refinement in both coincidence families (`verify/test_presentation.py`) | a lifted field whose commutant is strictly larger than the algebra of its finest presentation |
| $\Psi^0_h(u)=s(h)u$ — no distinguished endomorphism | zero-field collision at $n=1$; the no-translation half is just linear rigidity at $M=0$, not new | an algebraically natural family whose $\Psi^0_h$ is not a scalar |
| Only $\otimes$ transfers between algebras: not sub, quotient or product | three explicit witnesses (`verify/test_relations.py`) | a method natural over $A$ and not over some $A\otimes B$, or a proof that some other relation transfers |
| $\operatorname{Aut}(A)$-equivariance is trivial on the base; it generates $\mathfrak{gl}(m-1)$, not $\mathfrak{gl}(m)$ | $d(1)=0$ for every derivation; measured $9,25,49$ at $m=4,6,8$ (`verify/test_aut.py`) | an automorphism-induced map moving the real locus, or a generated algebra exceeding $\mathfrak{gl}(m-1)$ |
| The *pointwise* commutant reads $\rho_A(A)$ and does not collapse on the real locus | 17 based algebras $\times$ $n=1,2$ at generic/real-locus/origin points; multiplicative; separating (`verify/test_order2.py`) | a lifted jet whose pointwise commutant differs from $\rho_A(A)$, at any point |

**$N_0$ has now survived the test that killed its predecessors.** It is the newest item, and it is
the kind of clean combinatorial law that failed three times here: $2^{\text{aromas}}$, then
$2^{\text{loops}}$, then any power of two — each held in the $GL$ world where the examples live and
broke once metric contractions entered. So it was stress-tested on forests with several interacting
decorations (`verify/test_n0.py`). The decisive pair: the theta graph $\lvert D^2X\rvert_F^2X$ and the
minimal residual forest $\partial_{jk}X^i\partial_jX^a\partial_kX^a$ carry **identical** decorations —
two lianas and one stolon — yet measure $3$ and $1$. So $N_0$ is structural, not a tally. It is also
multiplicative over disjoint cycles ($\operatorname{div}(X)\lvert DX\rvert_F^2X$ measures $4=2\times2$).
The remaining gap has since been closed too: with three or more cycles sharing vertices,
$\lvert D^3X\rvert_F^2X$ measures $4$ (verified here on a quartic field, and matching the count
$1+3$: one stolon can raise exactly one of three lianas), $\lvert D^4X\rvert_F^2X$ measures $5$, and
a four-node necklace measures $5$ — a genuine independent-set count rather than a product. $N_0$ is
now the best-corroborated of the combinatorial claims.

## Cost discipline

The binding constraint is the account's rolling limit; **subagents dominate it** — seven completed
research agents cost ~1.09M tokens, individually up to 240k. Before spawning, ask whether a short
exact computation settles the question instead: historically, direct verification has caught as
many errors as agent cross-checking, far more cheaply.

**Before spawning, check the budget.** `get_session` (claude-code-remote) returns
`external_metadata.rate_limit_info.status`. On `allowed_warning`, spawn at most one agent; if the
window is nearly spent, do the work directly and say so.

**Spawn sequentially, never concurrently.** Three agents launched together race the same budget and
die together — that is exactly how eight agents were lost, one round losing all three at their
first tool call. Running them one at a time preserves independence (they still cannot see each
other) and guarantees that a budget exhaustion costs the last agent, not all of them.

**Require incremental checkpointing.** The eight crashed agents produced *nothing* salvageable,
because everything was written at the end. Instruct: *append each finding to `<file>` as soon as it
is established, before starting the next one.* Then a crash costs the remainder, not the whole run.

Also: cap the final message at ~400 words with detail routed to a file (a subagent's final message
is pasted verbatim into the caller's context, so long reports are paid for twice); point at
`references/literature.md` and forbid re-fetching papers; give the settled state as given and ask
exactly one open question; have them run `verify/run_all.py` rather than rebuild machinery —
`verify/dual.py` provides dual numbers, truncated jets $\mathbb R[\varepsilon]/(\varepsilon^r)$, Runge–Kutta
and helpers, and `verify/test_vdegree.py` differentiates by exact polynomial interpolation, which
avoids nested dual numbers.

## Layout

| path | contents |
|---|---|
| `HANDOFF.md` | **read first**: current target, open items, provenance |
| `flow-closure-under-differentiation.md` | the result; open items in §11, index calculus in Appendix A |
| `notes/` | full text of the Lemma L proof and its audit (rescued from ephemeral scratch) |
| `references/literature.md` | verbatim source quotes, marked verified vs agent-reported |
| `verify/` | exact-arithmetic checks, one module per claim family |
