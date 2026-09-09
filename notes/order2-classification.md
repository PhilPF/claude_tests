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

**Correction to the test point (recorded so it does not mislead).**  An earlier version of the
real-locus check zeroed the nilpotent *coordinates* (`u[j*N+a] = 0` for `a > 0`).  That is the real
locus only when `1_A = e_0`; for `R^2` and `R^3` in their **standard** bases the unit is `(1,..,1)`,
so the shortcut tested the wrong point for those two.  `verify/test_order2.py:real_pt` now builds
the honest point `1_A (x) v`.  Re-run with the honest point: same result, 0 discrepancies.

---

## F3. The detector is multiplicative, and it separates the factorizations  `[C]`

**(a) Multiplicativity.**  For `X = W^C` (so `X` is itself a lift, `dim A_r(X) = dim C`) and any
`A`, the pointwise commutant of the double lift has dimension exactly `dim A * dim C`:

    C = D, R^2   x   A = D, R^2, R[e]/(e^3)      6/6 cases, dim A_3(X^A) = dim A * dim C.

So the detector reads the **finest** presentation pointwise, matching the germ-level commutant
theorem of §7 — and it does so at every point, without the real-locus collapse.

**(b) Separation.**  For a fixed `m`, list the factorizations `(n, A)` with `n dim A = m`,
`dim A >= 2`, over the based algebras of `verify/algebra.py`:

    m = 4:  11 factorizations,  m = 6:  8 factorizations,
    all subalgebras rho_A(A) ⊂ End(R^m) pairwise **distinct** (0 coinciding pairs).

Distinct dimensions are separated by `dim A_r` alone; equal dimensions (`D` vs `R^2` at `n = 2`,
the five dim-3 algebras at `n = 2`, the eight dim-4 algebras at `n = 1`) are separated by *which*
subalgebra the commutant is.  Combined with F3(a): on the part of `L_{n,A}` where the seed `B_n` is
active (`X` primitive, `dim A_r(X) = 1`) the commutant is **exactly** `rho_A(A)`, so the active
loci of distinct factorizations are **pairwise disjoint and labelled by the detector**.

---

## F4. The simultaneous retraction — it exists, and it is one formula per `(n, N)`  `[C]+[A]`

This is the live obstruction of §7/HANDOFF item 2, at second order.  It dissolves once the
prescription is written **intrinsically in the commutant** instead of per-algebra.

### F4a. The prescription needs no per-algebra data  `[C]`

Let `a := A_r(u, j^r Z(u)) ⊆ End(R^m)` be the pointwise commutant and let `g_j` be the module
generators of `R^m` over `a`.  For a **unit-first** based algebra (`1_A = e_0`, `m = span(e_1..)`),
`g_j = e_{jN}` — *the same vectors for every algebra of dimension `N`*.  Then `R^m = (+)_j a.g_j`,
so every `u` has `a`-coordinates `(u_1..u_n)`, and for a polynomial field `Y` on `R^n`

    Y^A(u)  =  sum_i  Y_i(u_1,..,u_n) . g_i ,      products taken in a ⊆ End(R^m).

No basis of `A`, no structure constants, no augmentation, no real point.  **Verified exactly**
(`intrinsic_subst`): for `Y = D^2X[X,X]` and generic cubic `X`, against the true `(B_n X)^A(u)`
computed by `algebra.lift_map`, over all **17** based algebras, `n = 1, 2`, `m <= 4`, at a generic
point, at the honest real locus, and at the origin —

* all **17 unit-first** algebras: **match**, 3/3 points, every `n`;
* the two *standard-basis* products `R^2`, `R^3` (`1_A = (1,..,1)`): mismatch with `g_j = e_{jN}`,
  and **match** as soon as the true generators `1_A (x) e_j` are supplied.

So the only input beyond `a` is the generator convention, and the unit-first normal form makes it
`A`-independent.  Two earlier formulas were wrong and are recorded as traps: expanding around the
`e_0`-real part (`eps_jet` + Taylor) is correct for **local** algebras and **fails for every algebra
with more than one local factor** (`R^2`, `DxR`, `R^3`, `DxD`, `R[e]/(e^3)xR`, `R[x,y]/m^2 xR`,
`R^4`) — a product base change expands around *one real point per factor*, not one.  The
substitution form above is uniform and has no such split.

### F4b. The construction  `[A]`

Fix once and for all a **unit-first** basis for each based algebra (always possible; it is already
the convention for 15 of the 17 in `verify/algebra.py`).  Define, by induction on `m` in the
divisibility order,

    Phi_m(u, J)  :=  kappa_m(u, J)  +  sum over N | m, N >= 2  of  chi_N(J) . E_{m/N, N}(u, J)

* `E_{n,N}(u,J)`: recover `j X` from `J` by the fixed linear map, form `Y = B_n X`, and evaluate
  `Y` at the `a`-point `u` as in F4a, with `a = A_eps(J)`.
* `chi_N`: a smooth cutoff equal to 1 on `{dim A_r = N}` and supported in the open set where the
  commutant system has a spectral gap after `N` small singular values (there `A_eps(J)` is a
  **smooth** `N`-plane, kernel of a constant-rank linear system).
* `kappa_m`: free junk, supported in the open set `Omega_m = {dim A_r = 1}` (F2c).  For `m = 1`
  take `kappa_1` any local operator vanishing near the affine jets.

**Correctness on `L_{n,A}`.**  Let `J = j^r(X^A)(u)`, `N = dim A`, and let `c = dim A_r` of the
recovered `X`.  By F3(a) `dim A_r(J) = Nc` and `A_r(J)` is the algebra of the **finest**
presentation.  Then `chi_{Nc} = 1` and every other `chi` vanishes, `kappa_m = 0` (as `Nc >= 2`), and
`E_{m/(Nc), Nc}(u,J)` is the prescription of the finest presentation — which equals `(B_n X)^A(u)`
by the **consistency theorem** of §7 (both presentations factor through the commutant one).  So
`Phi_m = (B_n X)^A(u)` on every locus at once, active or not.

**What each recorded ingredient does here** — this is the "use them together" point:
* the *commutant* supplies consistency, so the several prescriptions at one jet agree;
* the *pointwise* commutant (F2) supplies a detector that does **not** collapse on the real locus,
  which is exactly what killed the fibrewise-linear route in §7;
* F3(b) says the detector *labels* the factorization, so the cutoffs are honest smooth functions;
* F4a says one formula covers **all algebras of a given dimension**, so the moduli of commutative
  algebras in dimension `>= 7` — the other recorded obstacle — never enter: the sum in `Phi_m` runs
  over **divisors `N` of `m`, not over algebras**.

**Self-consistency check — affine rigidity is not violated** `[A]`.  Affine rigidity is a
*consequence* of (*), so the construction must satisfy it at every affine `Z`, lifted or not.  It
does, for two separate reasons.  If `Z(u) = Mu + c` with `M` scalar, `A_r(Z)` is all of `End(R^m)`
(`dim = m^2`), the commutant system has `m^2` zero singular values, so there is no gap after `N < m^2`
of them and **every** `chi_N` vanishes: `Phi_m(Z) = kappa_m(Z) = 0`, affine.  If `M` is regular,
`dim A_r(Z) = m` (measured in O2: `m = 2,3,4` give `2,3,4`), so `chi_m = 1` with `n = 1`; but the
fixed linear recovery reads the `e_0`-components of an affine field, which are affine, so
`B_1` of them is `0` — provided `kappa_1` vanishes near the affine jets, which is exactly the
hypothesis already imposed on the seed.  Either way `Phi_m(Z) = 0`.  This is the mechanism by which
affine rigidity — the one genuine per-dimension consequence — is *forced on the seed and on nothing
else*.

**Residual hypotheses** (both argument, not computation): (i) the standard smoothing bookkeeping for
`chi_N` (a locally finite choice of the gap parameter); (ii) that `A_r(J)` is commutative with `R^m`
free over it whenever the prescription is nonzero — this is HANDOFF open item 3's residual point,
verified in 19 cases there and in the cases here, not proved.  Note (ii) is needed only where the
prescription is nonzero, i.e. where the seed is active, and there the commutant is exactly
`rho_A(A) (x) A(X)` by F3(a).

---

## F5. Answer to the classification question: **no**  `[A]`, on the verified F2–F4

> **Theorem (order-2 classification).**  Modulo the two residual hypotheses of F4b, a local smooth
> family `B = (B_m)` satisfies (*) **iff** it is obtained by the tower of F4b: `B_1` free subject
> only to affine rigidity, and for `m > 1`, `B_m` free off `L_m` and equal on `L_m` to the value the
> construction assigns.  (*) imposes the restriction to `L_m` and **nothing else**; every choice of
> the free data extends.  In particular the solution set is *not* the tree series.

**Explicit witness.**  Take `kappa_1 = 0`; `kappa_2` any nonzero smooth bump supported in
`Omega_2 = {dim A_r = 1}` (nonempty, `[C]` F2c) and depending explicitly on `u`; `kappa_m = 0`
otherwise.  Then `B_1 = 0`, `B_2 = kappa_2 != 0`, and `B_m` (`m > 2`) is given by F4b.  This `B`:

* satisfies (*) for **every** Weil algebra and finite product and every `n`;
* **vanishes identically on every lifted field and on every affine field**, so it is not a nonzero
  constant-coefficient tree series — for a tree `F(tau)(X^A) = (F(tau)X)^A != 0`, so a tree series
  vanishing on all lifted loci is zero;
* is **not affine equivariant**: the explicit `u`-dependence breaks translation covariance, and a
  bump support is not affine-invariant.

So: **algebraic naturality at second order does not force index calculus, and does not force affine
equivariance.**  The pathology is flatness **in the field**, in jet space, off the lifted loci —
*not* the `C^infty`-in-`h` artefact of §6, which by F1(3) is invisible at fixed order.

**Framing.**  This is a statement about closure/naturality *alone*.  It does not touch the recorded
relative results: *within* index calculus, algebraic naturality still forces exactly trees (loops
and stolons survive only for `A = R`, lianas iff `c_A = 1`), and every characterisation in this area
still factors as [equivariance] ∩ [closure].  What F5 says is only that the second bracket, taken by
itself, does not deliver the first — the conjecture "algebraically natural ⟹ affine equivariant" of
§11 is **false for smooth families at second order**.  It says nothing about the `f`-analytic
conjecture: `kappa_2` is a bump, and no analytic witness is produced here.  Indeed the construction
*cannot* be made analytic: an analytic `B_m` vanishing on a locus of positive codimension whose
span has codimension 3 (the coordinator's datum) is heavily constrained, and analytic junk cannot be
localised off the loci.  **The analytic half of §11's conjecture is untouched and is now the whole
question.**

---

## F6. A trap: jet-level overlaps are strictly larger than germ-level collisions  `[C]+[A]`

An attractive shortcut is to make the seed vanish near the affine jets, hoping that all pairwise
overlaps of the loci become inactive and the extension problem becomes a separation problem.  **That
is false**, and the reason is worth recording because it looks right.

Locality means `B_m` sees `j^r Z(u)`, so consistency is required at the level of *jets at one point*,
not of germs.  The collision trichotomy (§7) is a **germ** statement; the jet-level overlaps
`V_A(u) ∩ V_B(u)` are strictly bigger.  Measured (`r = 2`, `n = 1`, generic `u`, generators up to
degree 12–16 so the ranks are stable):

| A | B | dim V_A | dim V_B | dim of the overlap | dim of the affine part |
|---|---|---|---|---|---|
| `D` | `R^2` (unit-first) | 4 | 6 | 2 | 2 |
| `R[e]/(e^3)` | `R[x,y]/m^2` | 5 | 4 | 2 | 2 |
| `DxR` | `R^3` (unit-first) | 7 | 9 | **5** | 2 |
| `R[e]/(e^3)xR` | `R[x,y]/m^2 xR` | 8 | 7 | **5** | 2 |
| `R[x,y]/m^2 xR` | `R^4` (unit-first) | 7 | 12 | **5** | 2 |

Worked case `A = DxR`, `B = R^3` (both `n = 1`, `N = 3`), `u = (u0,u1,u2)`, `sigma = u0+u1`,
`tau = u0+u2`:

    f^A(u) = ( f(u0), f'(u0) u1, f(tau) - f(u0) ),   g^B(u) = ( g(u0), g(sigma)-g(u0), g(tau)-g(u0) ).

Equality of `r`-jets forces `f` and `g` to have equal `r`-jets at `u0` and at `tau`, and forces `g`
to be **affine to order `r` at `sigma`** with `g'(sigma) = f'(u0)` — and hence `f` affine to order
`r` at `u0`.  Free parameters: `2 + 3 = 5`, matching the measured overlap.  So the overlap carries
fields that are **not** affine (arbitrary jet at `tau`), and the seed does not vanish on it.  The
two prescriptions nevertheless agree, and the computation shows exactly why: components 1 and 3
match because the jets at `u0` and `tau` match, and component 2 reduces to

    B_1(l)(sigma) - B_1(l)(u0) = (sigma - u0) . B_1(l)'(u0)   for the affine field l tangent to f at u0,

which is **affine rigidity and nothing more**.  So the jet-level analysis reproduces the recorded
per-dimension result rather than strengthening it — a useful negative: locality does **not** upgrade
germ consistency into a new constraint.  It also means the overlaps must be handled by the
consistency theorem (as F4b does), not by inactivating them.

*Method note:* the "affine part" column above is `rank{ j^r(X^A)(u) : X affine }`; comparing it with
the overlap is what exposed the discrepancy.  And the ranks are unstable in the generator degree —
`R^3` needs degree `>= 8` before `dim V = 9` appears (3 evaluation points x 3 jet slots).  Testing
with degree 4 reports 5 and silently understates every product algebra.

---

## Reproducibility

`verify/test_order2.py` (wired into `verify/run_all.py`; the whole suite passes with it in).

| check | establishes |
|---|---|
| `O1 detector` | pointwise `A_3(u, j^3(X^A)(u)) = rho_A(A)`, 17 based algebras x `n=1,2` (`m<=6`) x {generic, honest real locus, origin} — **F2b** |
| `O2 detector` | `dim A_3 = 1` for a generic cubic in `m = 2,3,4`; `>= m` for affine fields — **F2c**, room for junk that is affine-rigid for free |
| `O3 multiplicative` | `dim A(X^A) = dim A . dim A(X)` on double lifts — **F3a**, the detector reads the finest presentation |
| `O4 separates` | the `rho_A(A)` for all factorizations of `m = 4, 6` are pairwise distinct — **F3b** |
| `O5 retraction` | the intrinsic substitution formula reproduces `(B_n X)^A(u)` for all 17 based algebras at all three kinds of point — **F4a** |

Helpers usable elsewhere: `pcommutant`/`pcom_dim` (pointwise commutant of a jet), `real_pt` (the
honest real locus `1_A (x) v` — see the correction above), `gen_field(n)` (a generic cubic field in
**any** dimension; `test_affine.gen_map` caps at `n = 3`), `locus_rows` (a spanning set of
`V_A(u)`), `intrinsic_subst` (base change read off from the commutant alone).

---

## Where this leaves the open items

* **HANDOFF item 2 (construct the flat, coherent, non-equivariant algebraically natural method, or
  prove none exists): resolved at second order, positively** — F4/F5, modulo the two residual
  hypotheses of F4b.  The simultaneous retraction exists; the two recorded obstacles (collapse of
  `dim V_A(u)` on the real locus; moduli of commutative algebras in dim `>= 7`) are both bypassed,
  the first by using the *pointwise* commutant, the second by writing the prescription
  intrinsically so that the sum runs over **divisors of `m`, not over algebras**.
* **HANDOFF item 1 (algebraically natural ⟹ affine equivariant ⟹ B-series):** the first
  implication is now **false for smooth families**, at second order, hence in general.  The
  `f`-analytic version is untouched and is now the entire content of the conjecture.  The natural
  next target is therefore: *can an analytic (or polynomial) `B_m` vanish on `L_m` without
  vanishing?*  The coordinator's datum — the span of the nontrivial lifted loci has codimension
  exactly **3** in `J^1_u` at `m = 4, 6`, the annihilator being the `m`-adic triangularity
  conditions — is the right first probe: at 1-jet level the polynomial ideal of the union is
  generated by those three linear forms, so a *polynomial* `B_m` of 1-jet order vanishing on `L_m`
  must lie in that ideal, which is a very small space.  Extending that ideal computation to `r = 2`
  would decide the polynomial case.
* **HANDOFF item 3's residual point** (`A(Z)` commutative with `R^m` free over it) is now load
  bearing for F4b, not merely tidy: it is what makes "the finest presentation" well defined at every
  jet where the prescription is nonzero.
