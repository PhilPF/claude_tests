# Flow closure under differentiation

What structure does the closure property
$$T(\Psi^X_h) = \Psi^{TX}_h$$
force on an approximation of the flow — "the derivative of the flow is the flow of the
derivative" — with the exact flow itself as one member?

**Read [`flow-closure-under-differentiation.md`](flow-closure-under-differentiation.md).**

## Headline

* $T$ is the Weil functor of the dual numbers $D=\mathbb{R}[\varepsilon]/(\varepsilon^2)$, and $TX$ is not a new
  vector field but $X$ **base-changed** $\mathbb{R}\to D$. So the closure property says exactly:
  *the method commutes with base change.*
* **The flow and the integrators satisfy it for the same reason**: both are the *unique*
  solution of an equation natural in the algebra base, and uniqueness is preserved under base
  change. Nothing about the Butcher tableau is involved. This criterion is sufficient but
  **not necessary**.
* **The obstruction is not the trace.** Of the invariant contractions: loops (aromas) die by a
  factor $\dim_\mathbb{R} A$; stolons (scalar products) die by failure of $D$-balancedness;
  **lianas (Laplacians) survive** — and only over $D$, by the accident $\varepsilon^2=0$.
* **Closure does not imply B-series.** Partitioned RK and the Laplacian method are
  (T)-natural and not B-series; the (T)-natural class has the cardinality of the continuum.
  Useful characterisations are relative to an equivariance group scheme.
* **Closure is not a jet condition**, unlike every classification theorem in the area.

## Verification

Every computational claim is checked in exact rational arithmetic — no floating point, no
dependencies (`sympy`/`numpy` are not required):

```
python3 verify/run_all.py
```

Forward-mode AD *is* evaluation over the dual numbers, so running a method over `Dual` computes
$T(\Psi^X)$ directly, and it is compared against $\Psi^{TX}$ computed on the lifted field.
Lines marked `FAIL` are **expected** failures — methods that must *not* satisfy the closure
property; the suite passes iff every outcome matches its prediction.

## Status

Produced by an adversarial multi-agent process: an initial synthesis, then three rounds in
which independent agents attacked it, with every load-bearing computation re-verified
independently. Several central claims of the first synthesis were **refuted** along the way
(trace-uniqueness; "closure $\Rightarrow$ B-series"; the sufficiency criterion being necessary).
Open questions are listed in §11 of the main document.
