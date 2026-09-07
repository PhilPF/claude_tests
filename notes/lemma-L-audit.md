# Audit of the Lemma L proof — full text

Verdict: the proof stands, with corrections. Rescued from the ephemeral scratchpad.

# Audit of the Lemma L proof (`r5/lemmaL-findings.md`)

Adversarial audit. Status markers: **[OK]** checked and correct, **[ERR]** error found,
**[FIX]** defect that is repairable and how, **[T]** taken on trust, **[C]** my own conjecture.
Audit code lives in `r5/audit_work/` (written from scratch, not reusing the author's `r5/work/`
except where explicitly noted).

---

## 0. Log begins

## A. Reproduction of the author's own computations — [OK]

`python3 r5/work/run_lemmaL_checks.py` runs and prints `ALL LEMMA-L CHECKS PASSED`. Every number
quoted in `lemmaL-findings.md` §1, §3, §4b, §5 reproduces exactly (Brauer ranks 1/3/15/105 with
thresholds n>=k; bridge 45/45; orbit counts 2, 16, 19, 6; ranks; collapse 19->19, 15->13, 12->10,
7->6). Nothing is fabricated and no file referenced in §7 is missing. The mod-p argument
(rank_{F_p} <= rank_Q, so a *full* rank mod p is rigorous) is correctly applied: only full-rank
entries are used.

Code read line by line; `orbits.py` is a correct computation of
dim span{elementary differentials} for a fixed profile:

* rows = Gtil-orbit sums, columns = set partitions pi of H with <= n blocks, entry
  #{M' in orbit : M' refines pi}. This is (Phi_M o s)(e_I) up to the row-constant
  |Stab|/|Gtil|, so the rank is right.
* dropping columns that coarsen no admissible matching is legitimate (Gtil preserves
  admissibility, so those columns are identically zero).
* rank{Phi_M|_W} = rank{Phi_M o s} and the jet tensors span W, so this rank really is the
  dimension of the span of the differentials. **This part of the verification is sound and
  non-circular.**

## B. [ERR, methodological] §4b's "adversarial" test of Step 2 is vacuous in the direction it
claims to test

§4b asserts: "A Gtil that is too large fuses matchings with different values". **This cannot
happen, for any subgroup whatsoever that fixes the jet tensors.** Concretely: Phi_M o g =
Phi_{g^{-1}M}, and the jet tensor J = (x)_k c_k^{ox m_k} (x) v^{ox d} (x) u is by construction
Gtil-invariant; hence for g in Gtil,

        Phi_{gM}(J) = Phi_M(g^{-1}J) = Phi_M(J).

So *every* matching in a Gtil-orbit has the same value on every jet tensor, automatically,
whether or not Gtil is the right group. `step2.py` therefore verifies only the converse
inclusion (equal value => same orbit), which is logically the same content as the rank test of
§3, not an independent check. The claim in §4b that the test detects an over-large group is
false, and the same objection applies to the framing in §5 ("Trap 1").

Consequence: the brief's proposed remedy ("compare group orbits against equality classes of
actual computed differentials") is itself blind to an over-large group. An over-large Gtil would
make two non-isomorphic forests have *identical* differentials; the value test would then still
report AGREE, and the rank test would still report rank = #orbits. Only a count of iso-classes
computed **independently of Gtil** can decide it. That is done in section C below.

## C. Step 2 (the hinge), tested against an independent encoding — [OK]

`r5/audit_work/indep.py` (written from scratch: my own matching enumerator, my own
admissibility test, and a **slot-forgetting** encoding of a decorated forest as
`(valences, root, multiset of arrows src->tgt, multiset of lianas, multiset of stolons,
v-multiplicity vector)`, canonicalised by brute-force minimisation over all m! node
relabellings). This never mentions Gtil, so comparing the two partitions of the admissible
matchings is a real test of "Gtil-orbits = iso-classes", including the too-large direction.

Result: the two partitions coincide on every profile tested, including six the author never
tried:

| k | d | \|H\| | adm. M | Gtil-orbits | independent iso-classes |
|---|---|---|---|---|---|
| (1,1) | 1 | 6 | 4 | 2 | 2 |
| (2,1,1) | 0 | 8 | 45 | 16 | 16 |
| (2,1,1) | 2 | 10 | 108 | 19 | 19 |
| (3,1) | 1 | 8 | 24 | 6 | 6 |
| (2,1,0) | 1 | 8 | 27 | 15 | 15 |
| (2,2,0) | 2 | 10 | 108 | 12 | 12 |
| (1,1,0,0) | 1 | 8 | 24 | 7 | 7 |
| (2) | 2 | 6 | 2 | 1 | 1 |
| **(2,2,2)** | 0 | 10 | 315 | **16** | **16** |
| **(3,2,1)** | 0 | 10 | 315 | **48** | **48** |
| **(4,2)** | 3 | 12 | 720 | **10** | **10** |

So **Gtil is neither too small nor too large**, relative to the ordinary notion of isomorphism
of a decorated exotic aromatic forest. The `(4,2), d=3` row is the strongest: every one of its
10 forests has a node carrying two or three v's, and several carry a liana and v's on the same
node simultaneously; the `(S_k)^{m_k}` factor is legitimate there.

Independently of the computation, the too-large direction is in fact a triviality: for
g in Gtil the node permutation underlying g *is* an isomorphism gamma_M -> gamma_{gM}, so
orbits can never be coarser than iso-classes. The genuine content of Step 2 is the injectivity
half, and that is also correct.

### C1. The "decoration set" wording in Step 1 — sloppy but **not** an error

Step 1 defines the decoration of gamma_M as the *set* `{t : some l_{t,s} is matched to some w_j}`,
which discards the multiplicity, while admissible matchings (and the code, and the application at
v-degree 2, e.g. d^2X[v,v]) allow several v's on one node. I expected this to break injectivity.
It does not: in gamma_M node t keeps all k_t in-slots, so

        p_t  =  k_t - #(arrow heads at t) - #(liana ends at t)

is *recoverable* from gamma_M, and the set is redundant data. Checked: replacing the multiplicity
vector by its support changes no class count on any profile above, including (2,1,1) d=2 (6 of the
19 forests are multi-v) and (4,2) d=3 (all 10 are). **[OK]**

The residual defect is only in the *statement* in §0, "S subseteq V ... d = |S|", which literally
excludes the multi-v forests that the d=2 application needs. The parenthetical there says the
general case works verbatim, and it does; the statement should simply say "a multiplicity function
p : V -> Z_{>=0} with sum p = d".

## D. A direct certificate of Lemma L that bypasses Steps 1-3 — [OK], and a sharper picture of
the dimension

`r5/audit_work/direct.py` + `xcheck.py`. For each profile I enumerate iso-classes with the
independent encoding of section C, evaluate `<u, F(hat gamma)(X)(0)[v^{ox d}]>` with an evaluator
written straight from the definition (sum over one index variable per arrow / liana / stolon /
v-slot, plus the root index), and take the exact rank over Q of the
(iso-classes) x (random rational fields, v, u) value matrix. Full rank is a *rigorous certificate*
of independence in that dimension which uses **no** claim about W, about polarisation, or about
Brauer matchings — it therefore checks Steps 1, 3 and 4 as a black box.

My evaluator was first validated against the project's own `verify/forests.py::Forest.F`: exact
agreement over Q on **all 45 + 315 + 315 = 675** admissible matchings of profiles (2,1,1), (2,2,2),
(3,2,1) with degree-3/degree-4 rational fields.

| k | d | n | N | iso-classes | rank over Q |
|---|---|---|---|---|---|
| (1,1) | 1 | 3 | 3 | 2 | 2 |
| (2,1,1) | 0 | 2 / 3 / 4 | 4 | 16 | 11 / **16** / 16 |
| (2,1,1) | 2 | 3 | 5 | 19 | **19** |
| (3,1) | 1 | 3 | 4 | 6 | 6 |
| (2,1,0) | 1 | 3 | 4 | 15 | 15 |
| (2,2,0) | 2 | 3 | 5 | 12 | 12 |
| (1,1,0,0) | 1 | 3 | 4 | 7 | 7 |
| (2,2,2) | 0 | 3 / 4 / 5 | 5 | 16 | **15** / 16 / 16 |
| (3,2,1) | 0 | 3 | 5 | 48 | 48 |
| (4,2) | 3 | 3 | 6 | 10 | 10 |

Lemma L holds on every profile tested, including three the author never tried and the
`(4,2), d=3` case where every forest carries multiple v's on one node.

**New data point on the dimension.** Profile (2,2,2), d=0 is the first case found in this project
where the differentials are **dependent**: at n=3 the 16 iso-classes span only 15 dimensions, and
full rank is reached at n=4. This is not a counterexample — n_0 = N = 5 there, so the lemma claims
nothing at n = 3 — but it is a genuine rank drop and it is consistent with, and mild evidence for,
the author's **[C]** conjecture n >= |V| + 1 (here |V| = 3, deficient at n = 3, full at n = 4).
It also shows the conjectured bound is not necessary: (2,1,1) with d=0 and with d=2 both have
|V| = 3 and are already full at n = 3.

## E. Steps 1 and 3 — [OK]

`r5/audit_work/step13.py`.

* **Step 1 bookkeeping.** On every one of the 1690 admissible matchings across ten profiles,
  `m = 1 + a + 2s`, `sum_t k_t = a + 2l + d`, `|H| = 2N` with `N = 1 + a + l + s + d`, and
  `sum_t p_t = d`. So |H| is always even and n_0 = |H|/2 = #arrows + #lianas + #stolons + d + 1
  exactly as the lemma states. The relation (*) itself I re-derived by hand (both sides are the
  same sum over index assignments, with no combinatorial factor) and checked numerically in
  section D via `Forest.F`.
  One convention has to be read correctly and the write-up does state it: k_nu counts **all**
  lower slots of node nu **including** its v-slots, so decorating a node changes its valence and
  hence the profile. Step 0's multihomogeneity separation is consistent with that.

* **Step 3, the identification of W.** Computed dim (V^{ox H})^{Gtil} by Burnside (orbits of Gtil
  on maps H -> [n]) and compared with
  `n * C(n+d-1,d) * prod_k C(D_k + m_k - 1, m_k)`, `D_k = n*C(n+k-1,k)`, i.e. with
  `(x)_k Sym^{m_k}(V (x) Sym^k V) (x) Sym^d V (x) V`. Exact agreement in all 12 cases tested
  (profiles (1,1),(2,1,1),(2,2,0),(2,2,2),(3,2,1); d = 0,1,2; n = 2,3). The symmetry assumed of
  c_k -- symmetric in its k *lower* indices only, not involving the upper index -- is the true
  one; a stolon or a liana constrains the *matching*, not the tensor, so it cannot interfere.

* **The polarisation/span step loses nothing** -- sharp evidence: the *theoretical* rank computed
  by the author's `orbits.py` from Phi_M o s (i.e. assuming jet tensors span W) and my
  *empirical* rank from actual random polynomial fields agree even where both are **deficient**:
  profile (2,1,1), d=0 gives 11 at n=2 and 16 at n=3 by both routes; (2,1,1), d=2 gives 19 at n=3
  by both. Agreement on a deficient value is much stronger evidence than agreement on a full one.

## F. Smaller points

* **The x=0 reduction — [OK], the worry in the brief does not arise.** We are *proving*
  independence, so we may choose the base point: a relation valid for all (X,x) is in particular
  valid at x=0, and every jet at 0 is realised by a polynomial field, so no differential is lost.
  A differential vanishing at 0 but not elsewhere is impossible here for the same reason (the jet
  at 0 is free). The reduction *would* fail if the relation were only known for X in a proper
  subclass -- the write-up says this itself in §6, and it is the right caveat.
* **The Brauer threshold is applied correctly — [OK].** |H| = 2N, so Proposition B is invoked with
  k = N and needs n >= N = n_0. Only the direction "n >= k => independent" is used, which is the
  easy half. Proposition B itself I take from the literature **[T]** (Brauer/Weyl; Brauer algebra
  B_k(n) -> End_{O(n)}(V^{ox k}) injective iff n >= k) and it is confirmed numerically for
  k <= 4 by the author's `sft.py`, which I re-ran.
* **[ERR, provenance] §5 of the findings misdescribes the cache.** It states "The cached quotes do
  not contain the statement of Prop. 4.1". They do: `references/literature.md` carries a **[V]**
  verbatim quote of its conclusion, "In particular, the elementary differential map F is injective
  on Span(Gamma)", together with LMK's warning about the theta-free dual field. What the cache
  lacks is only the *construction* of the theta-parametrised fields. In a section whose whole
  purpose is provenance this should be corrected.
* **Trap 1 ("the thetas are mandatory") — the write-up's handling is right, but for a reason worth
  stating.** The bolded line in `literature.md`, "Any proof of the two-coloured analogue must carry
  the thetas", is a project **gloss**, not a quote: what LMK actually assert is that *their*
  theta-free **witness field** fails to separate one pair. The new proof constructs no witness
  field; it tests against the whole jet space and buys its separation from the dimension
  hypothesis n >= N instead. LMK work in dimension |gamma|+1, far below N, which is exactly why
  they need extra parameters. The two accounts are consistent, and the gloss is an over-reading.
  I could not check LMK's own dimension proviso "2d >= |kappa| + |kappa'| + 1" against n >= N
  because the cache does not record what d, kappa, kappa' denote there; the proof does not use it.

## G. Secondary: the N_0 law where three or more generalised cycles share vertices — [OK]

`r5/audit_work/n0.py`. Test used is the project's own (base block of F_{2n}(gamma)(TX) at v=0
must equal N_0(gamma) * F_n(gamma)(X), exactly over Q), on a degree-5 rational field in n=2.
The brief was right that `verify/test_n0.py` never exercises overlapping cycles: its five cases
have at most two cycles. New cases:

| forest | families (=N_0) | single cycles | conflicting pairs | mult = N_0 |
|---|---|---|---|---|
| \|D^2X\|_F^2 X (theta, control) | 3 | 2 | 1 | yes |
| **\|D^3X\|_F^2 X** | **4** | 3 | **3** | yes |
| **\|D^4X\|_F^2 X** | **5** | 4 | **6** | yes |
| **necklace: 4 nodes, 2 stolons, 4 lianas** | **5** | 3 | **2** | yes |
| div(X)^3 X (3 disjoint cycles, control) | 8 | 3 | 0 | yes |

The necklace is the sharpest: three generalised cycles of which one pair is vertex-disjoint and
two pairs are not, so N_0 = 1 + 3 + 1 = 5 is a genuine independent-set count -- neither 2^c nor
1 + c. |D^kX|_F^2 X gives N_0 = k + 1 for every k, from k pairwise-conflicting cycles.

Also swept **exhaustively**: every iso-class of profiles (2,2,0), (2,2,2), (3,3,0), (3,2,1),
(2,1,1) at d=0 -- 106 forests in total -- with mult = N_0 in every single case, 0 failures.
The disjointness condition therefore really is the right one, and it bites.

## H. Secondary: the 62 vs 84 enumerations — reconciled, no contradiction [OK]

Both numbers reproduce, and I re-derived both independently with the `indep.py` enumerator:

* **62** = `scratchpad/run_enum.py`: all rooted exotic aromatic forests with **<= 3 nodes** and
  **<= 2 derivative slots per node**, with *no cap on the total number of slots* (so up to 6).
* **84** = `scratchpad/r4/xforest.py::enumerate_forests(pmax=4, totq=4)`: **<= 4 nodes**, **<= 2
  slots per node**, **and total slots <= 4**.

They are two different truncations of the same set and **neither contains the other**:

        |A| = 62,  |B| = 84,  |A n B| = 46,  |A \ B| = 16,  |B \ A| = 38,  |A u B| = 100

`A \ B` is exactly the 16 forests with 3 nodes and 6 total slots (killed by `totq=4`);
`B \ A` is the 38 forests on 4 nodes (killed by `pmax=3`). So there is no disagreement to
resolve, and, contrary to the brief's guess, **neither enumerator applies a connectivity
filter** -- I read both; the difference is purely the node/slot caps.

Two labelling defects fall out:
* `scratchpad/r4/test_defect.py` prints "<=3 nodes, <=3 derivative slots" while calling
  `enumerate_forests(pmax=4, totq=4)`, i.e. <=4 nodes and <=4 slots. The printed scope is wrong.
* `flow-closure-under-differentiation.md` line 472 says the defect characterisation was "Verified
  exhaustively on all 84 iso-classes with <=4 nodes". That is under-specified: "all iso-classes
  with <=4 nodes" is an infinite set. It should read "<= 4 nodes, <= 2 derivative slots per node,
  <= 4 slots in total". As it stands the sentence overstates the verified scope.

## I. Two further defects in the *statement*, and one thing the proof gets right by luck

* **[FIX] The "equivalently" clause of the Theorem in §4 is wrong as written.** It reads
  "*equivalently, F_n is injective on the span of the v-decorated exotic aromatic forests of total
  size <= n*". Substituting m = 1 + a + 2s into N = 1 + a + l + s + d gives

        N  =  m - s + l + d          (m = #nodes, s = #stolons, l = #lianas, d = v-degree)

  so N is **not** the node count, and it errs in *both* directions: the theta graph has m = 3 and
  N = 4 (the clause would claim injectivity one dimension too early), while `<X,X>X` has m = 3 and
  N = 2 (the clause would claim less than is proved). Read with "size" := N the clause is correct;
  as it stands it is not, and it is exactly the quantity the **[C]** conjecture n >= |V| + 1 is
  about, so the two must not be conflated. Verified the identity N = m - s + l + d on all 1690
  admissible matchings of section E.
* **[OK] The bound "for n large" is harmless for the target theorem** -- worth stating because the
  brief asks whether n >= N is even correct. §1 of the main write-up defines a method as a
  *sequence* {Psi_n}_{n>=1}, so the coefficients b(gamma) are dimension-independent and killing
  them in one sufficiently large dimension kills them in all. The gap between n >= N and the
  conjectured n >= |V|+1 therefore costs the O(n) theorem nothing, as §5 of the findings says.
  I hunted for a counterexample to the conjecture on profiles where |V|+1 < N -- (3,2), (4,3),
  (3,3,0) at d=0, (2,2), (3,3), (4,2) at d=1 -- and found none; in every one independence already
  holds at n = |V|, so the conjecture survives but is itself not sharp. **[C] remains conjecture.**
* **[T] Out of scope, but load-bearing for the Corollary.** Lemma L gives independence of the
  F(hat gamma). The Corollary's step "the v-degree-2 equation forces each b(gamma) separately to
  vanish" additionally needs that, for each residual gamma, *some* decorated hat gamma occurs with a
  nonzero coefficient in gamma's v-degree-2 defect. That is asserted here (and in §11 of the main
  write-up), not proved in this document. I did not check it; the brief says everything outside
  Lemma L is already established.
* **[C] The multi-colour extension in §6** ("several fields X^(1..c), colours arbitrary; let Gtil
  permute only nodes of equal (colour,valence); Step 2 goes through verbatim") is plausible and I
  see no obstruction -- the jet tensors of distinct fields are independently free, so W just
  acquires more factors -- but I did not verify it computationally. Flag it as unverified.

---

# VERDICT

**The proof STANDS, WITH CORRECTIONS.** I could not break it. The chain
Step 0 -> Step 1 -> Step 2 -> Step 3 -> Step 4 is valid, and every link was checked both by
re-derivation and by exact computation with code written independently of the author's.

**Verified here (exact over Q):**
Step 1's slot bookkeeping on 1690 matchings; Step 2's bijection against a slot-forgetting
encoding on 11 profiles (three new), including profiles where every forest carries multiple v's
on one node; Step 3's identification of W by Burnside on 12 cases; the whole lemma by a *direct*
rank certificate on 10 profiles at several dimensions, using an evaluator validated against
`verify/forests.py::Forest.F` on 675 forests and against hand-written formulas at d = 1, 2;
the author's own five scripts, re-run and reproduced exactly.

**Taken on trust:** Proposition B (Brauer/Weyl; verified numerically only to k = 4); everything in
the O(n) characterisation outside Lemma L, in particular the triangularity input the Corollary
needs.

**Corrections required, none fatal:**
1. **[ERR]** §4b's claim that its test detects an over-large Gtil is false (section B). The test is
   sound but proves only what §3 already proved; the "adversarial check of the hinge" label should
   go, or the test should be replaced by the iso-class comparison of section C, which does have
   that content. The same over-claim appears in §5's discussion of Trap 1.
2. **[ERR]** §5's provenance claim that `references/literature.md` does not contain the statement of
   LMK Prop. 4.1 is false; it contains a [V] verbatim quote of its conclusion (section F).
3. **[FIX]** The "total size <= n" clause of the Theorem (section I).
4. **[FIX]** §0's "S subseteq V ... d = |S|" should be a multiplicity function; Step 1's "decoration
   set" should be the multiplicity (harmless as it stands -- section C1 -- but only by an accident
   worth stating explicitly, namely that p_t is recoverable from gamma_M).
5. **[FIX]** The main write-up's "all 84 iso-classes with <=4 nodes" is under-specified (section H).

**Not a defect but worth adding**, since it is the only real risk in the argument: the too-large
direction of Step 2 cannot be tested against computed differentials at all, for the structural
reason in section B, and the reader should be told that it is instead a triviality (the node
permutation underlying g in Gtil *is* an isomorphism gamma_M -> gamma_{gM}).
