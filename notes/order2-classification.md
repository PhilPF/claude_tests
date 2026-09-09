# Order-2 classification of algebraically natural families

Working file for the question: write an algebraically natural family as
`Psi^X_{m,h}(u) = u + h X(u) + h^2 B_m(X)(u) + O(h^3)`; algebraic naturality is exactly

  (*)   B_{n dim A}(X^A) = (B_n X)^A   for every based Weil algebra or finite product (A, iota),
        every n, every X.

**Deliverable:** classify the local families `B = (B_m)` satisfying (*) — index calculus (a
constant-coefficient combination of tree elementary differentials) or not.

Findings are appended as soon as they are established.  `[C]` = checked in exact rational
arithmetic (file named); `[A]` = argument only.

---

## F1. Bookkeeping: what (*) is, and what it is not  `[A]`

(*) is a condition indexed by **factorizations** `m = n·N` together with a based algebra `(A,iota)`
of dimension `N`.  For fixed `m` it pins `B_m` on the union of lifted loci
`L_m = union over (n,A), N = dim A >= 2, of { X^A : X in Xf(R^n) }`,
and says nothing off `L_m`.  Three immediate bookkeeping consequences.

1. **`N = 1` is vacuous.**  `A = R` gives `B_n(X) = B_n(X)`.  So `B_1` carries **no** condition
   from (*) at all — dimension 1 has no nontrivial factorization.  Every constraint on `B_1` must
   come from a *collision* in a higher dimension, i.e. from affine rigidity, which is all the
   collision trichotomy allows (§7).  Concretely `B_1` is constrained only by:
   `B_1(X)` is affine whenever `X` is affine.

2. **The support set is an up-set for divisibility.**  Put `S = { m : B_m != 0 }`.  If `B_n != 0`
   then for every `N >= 2` and every `A` of dimension `N`, `B_{nN}` restricted to `{X^A}` equals
   `(B_n X)^A`, which is nonzero because base change is injective (it restricts to the identity on
   the real locus).  So `n in S  ==>  nN in S` for all `N`.  Dually `B_m = 0` forces `B_n = 0` for
   every proper divisor `n | m`.  **Junk cannot be confined to one dimension**: it always
   propagates upwards by base change.  This is why "put a bump where no lifted field lives" fails
   as stated — see F5.

3. **Flatness in `h` is unavailable.**  §6's counterexample `x + hX + theta(h div X) hX` is flat in
   `h`, so its `h^2` coefficient is **zero**: it is invisible at any fixed order.  The order-2
   question is therefore entirely about flatness *in the field* (in jet space), not about the
   `C^infty`-in-`h` artefact.  This narrows the problem genuinely, and it means an order-2
   counterexample must be a *jet-space* construction.

So the order-2 problem has exactly the architecture recorded for the full problem, with no
simplification from truncating in `h`: seed free in dimension 1 up to affine rigidity, prescribed
on `L_m` for `m > 1`, free off `L_m`, and the whole question is the **simultaneous smooth local
retraction**.  What follows is an attack on that retraction, using the recorded constraints
together.

---

## F2. The pointwise commutant is a non-degenerate detector  `[C]` (`verify/test_order2.py`)

Locality makes `B_m(Z)(u)` a function of `j^r Z(u)`, so the object to work with is the **pointwise**
commutant of a jet, not the germ commutant of §7:

```
A_r(u, j^r Z(u)) = { L in End(R^m) : L D^k Z(u)[v_1..v_k] = D^k Z(u)[L v_1, v_2..v_k], 1<=k<=r }.
```

**(a) Containment, always** `[A]`.  If `Z = X^A` then each `D^k Z(u)` is `A`-multilinear over `A`,
so `rho_A(A) ⊆ A_r(u, j^r Z(u))` **at every point `u`**, for every `X`, every `r >= 1`.  Hence
`dim A_r >= dim A`.

**(b) Equality, and no degeneration on the real locus** `[C]`.  Measured with `r = 3`, for a generic
cubic `X`, over all 17 based algebras of `verify/algebra.py` and `n = 1, 2` (m <= 6), at three kinds
of point — generic `u`, a point of the **real locus** `1_A (x) R^n`, and `u = 0`:

    dim A_3(u, j^3(X^A)(u)) = dim A = N  and equals rho_A(A), in all 3 x 24 cases.

This is the decisive contrast with the recorded obstruction.  `dim V_A(u)` (the dimension of the
space of lifted jets) **does** collapse on the real locus — that is what killed the fibrewise-linear
route in §7 — whereas the pointwise commutant does **not**: it reads `N` uniformly, including at
`u = 0`.  So the *presentation* is detectable pointwise even where the *locus* degenerates.

**(c) Consequences** `[A]`.
* `Omega_m := {(u,J) : dim A_r(J) = 1}` is **open** (kernel dimension is upper semicontinuous),
  nonempty for every `m >= 2` (a generic cubic field, `[C]` above), and **disjoint from every
  lifted locus** with `dim A >= 2`, by (a).
* For an **affine** field `Z(u) = Mu + c` all `D^k Z = 0`, `k >= 2`, so `A_r` is the centralizer of
  the single matrix `M`, of dimension `>= m`.  Hence **no affine field lies in `Omega_m` for
  `m >= 2`**: junk supported in `Omega_m` is automatically affine-rigid, for free.

So per dimension there is plenty of room: for every `m >= 2` there is a nonzero local smooth
operator vanishing identically on every lifted locus *and* on every affine field.  By F1(2) this is
**not** by itself a family satisfying (*) — junk propagates — but it fixes where junk may sit.
