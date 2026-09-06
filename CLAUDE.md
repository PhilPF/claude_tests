# Working notes for this repository

Research project on **flow closure under differentiation**: what the property
$T(\Psi^X_h)=\Psi^{TX}_h$ forces on an approximation of the flow. Read
`flow-closure-under-differentiation.md` for the state of knowledge; it is the single source of
truth and is kept correct rather than append-only — claims there have been *retracted* several
times when verification contradicted them.

## Ground rules

**Verify, don't trust.** Every computational claim in the write-up is checked in exact rational
arithmetic. `python3 verify/run_all.py` runs all of it (no dependencies; `sympy`/`numpy` are
unavailable and pip has no index). Lines printed `FAIL` are *expected* failures — methods that
must **not** satisfy closure; the suite passes iff every outcome matches its prediction.

**Test fields must be cubic or higher.** With a quadratic field $\Delta X$ is constant and
low-order differentials degenerate, which has produced spurious results twice.

**Canonicalise graphs before testing independence.** Two presentations of the same forest under
renaming of internal indices will manufacture a false linear relation.

**Distrust any statement that silently fixes an identification $\mathbb R^{n\dim A}\cong A^n$.**
This has been the source of three separate errors (the partitioned-RK lift; $c_A$ over
$\mathbb R[\varepsilon]/(\varepsilon^3)$; $c_A$ over the dual numbers themselves).

**Clean combinatorial laws in this problem are usually wrong.** $2^{\text{aromas}}$, then
$2^{\text{loops}}$, then any power of two — each looked right because it holds in the $GL$ world
where the examples live, and each broke once metric contractions entered. The correct multiplier
is $N_0$, a count of generalised cycles.

## Cost discipline

The binding constraint is the account's rolling token limit, and **subagents dominate it** — a
single research agent has cost up to 240k tokens. Before spawning one, check whether the question
can be settled by a short exact computation instead; historically, direct verification has caught
as many errors as agent cross-checking, far more cheaply.

When a subagent is warranted:

* **Two, not three.** Independent passes have converged closely; the third mostly replicates.
* **Cap the report.** Require: *final message ≤ 400 words; write full detail to a named file.*
  A subagent's final message is pasted verbatim into the caller's context, so long reports are
  paid for twice.
* **Point at `references/literature.md`** and forbid re-fetching papers. Several agents each
  downloaded the same paper and re-read the same definitions.
* **Forbid re-deriving settled material.** Give the settled state as given; ask exactly one open
  question.
* **Have them run `verify/run_all.py` rather than rebuild machinery.** `verify/dual.py` already
  provides dual numbers, truncated jets $\mathbb R[\varepsilon]/(\varepsilon^r)$, Runge–Kutta, and
  helpers; `verify/test_vdegree.py` differentiates by exact polynomial interpolation, which avoids
  nested dual numbers.

## Layout

| path | contents |
|---|---|
| `flow-closure-under-differentiation.md` | the result, with open items in §11 and Appendix A on the index calculus |
| `references/literature.md` | verbatim source quotes, marked verified vs agent-reported — quote from here |
| `verify/` | exact-arithmetic checks, one module per claim family |
