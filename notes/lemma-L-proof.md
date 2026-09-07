# Proof of Lemma L — full text

Rescued from the ephemeral scratchpad; this is the only surviving copy of the argument
summarised in `flow-closure-under-differentiation.md` §11. Audited in `lemma-L-audit.md`,
which found the bound should be on $N=m-s+l+d$ and corrected a methodological claim about
validating the symmetry group. Provenance: one research pass, one independent audit.

# Lemma L — findings

Running log; appended as each item is established. Status markers:
**[P]** proved here, **[V]** verified by exact computation here, **[L]** taken from the
literature, **[C]** conjecture.

---

## 0. Formalisation used throughout (matches `verify/forests.py`)

Half-edge model. An exotic aromatic forest $\gamma$ is:

* a finite node set $V$; each $\nu\in V$ carries **one upper half-edge** $o_\nu$ (its output
  index) and a finite set $\mathrm{In}_\nu$ of **lower half-edges** (derivative slots),
  $k_\nu=|\mathrm{In}_\nu|$;
* a **perfect matching** $M$ of $H:=\{o_\nu\}\sqcup\bigsqcup_\nu \mathrm{In}_\nu$ minus one
  distinguished free upper half-edge $o_r$ (the root).

Matched pair types: upper–lower = **arrow** (same node: loop), lower–lower = **liana**,
upper–upper = **stolon**. Isomorphism = bijection of $V$ together with bijections
$\mathrm{In}_\nu\to\mathrm{In}_{\alpha(\nu)}$ carrying $M$ to $M'$ and root to root; $\gamma$
denotes the iso-class.

Because in Euclidean coordinates every contraction (arrow, liana, stolon) is "give the two
half-edges the same index and sum", the elementary differential is uniform:
$$F_n(\gamma)(X)(x)^{i}=\sum_{\phi:M\to[n]}\ \prod_{\nu\in V}\ \partial_{\phi(\mathrm{In}_\nu)}X^{\phi(o_\nu)}(x),
\qquad \phi(o_r):=i .$$
This is the formula implemented in `Forest.F`.

**$v$-decoration.** $\hat\gamma=(\gamma,S)$, $S\subseteq V$: each $\nu\in S$ gets one extra
lower half-edge, matched to one of $d=|S|$ tensor slots all carrying the same fixed vector
$v$. (Nothing below uses $|S\cap\{\nu\}|\le1$; several $v$-slots per node work verbatim.)

**Lemma L (to prove).** For $n$ large the functions
$X\mapsto F_n(\hat\gamma)(X)(x)[v^{\otimes d}]$, over pairwise distinct $v$-decorated exotic
aromatic forests, are linearly independent.

---

## 1. [L]+[V] The one external input: independence of Brauer matchings

**Proposition B.** Let $V=\mathbb R^n$ with the standard inner product, $H$ a finite set with
$|H|=2k$, and for a perfect matching $M$ of $H$ let
$$\Phi_M\in (V^{\otimes H})^*,\qquad \Phi_M\Big(\bigotimes_{h}z_h\Big)=\prod_{\{a,b\}\in M}\langle z_a,z_b\rangle .$$
Then $\{\Phi_M\}_{M}$ is linearly independent **iff $n\ge k$**, and $\#\{M\}=(2k-1)!!$.

*Source.* Classical: this is the second fundamental theorem for $O(n)$ / the faithfulness of
the Brauer algebra, $B_k(n)\to\operatorname{End}_{O(n)}(V^{\otimes k})$ injective iff $n\ge k$
(Brauer 1937; Weyl, *The Classical Groups*; Goodman–Wallach §10.1). It is **not** taken from
Laurent–Munthe-Kaas, and it is not stated in `references/literature.md`; it is the standard
statement whose shadow is the LMK proviso "[$\delta$] is a bijection if and only if
$2d\ge|\kappa|+|\kappa'|+1$" quoted there. Only the direction "$n\ge k\Rightarrow$
independent" is used below.

*Verified here* (`work/sft.py`): $\Phi_M$ written in the monomial basis has columns indexed by
set partitions $\pi$ of $[2k]$ with $\le n$ blocks and entry $[\,M\preceq\pi\,]$; rank computed
by elimination.

| $k$ | rank at $n=1,2,\dots$ | $(2k-1)!!$ |
|---|---|---|
| 1 | 1, 1, 1 | 1 |
| 2 | 1, **3**, 3, 3 | 3 |
| 3 | 1, 10, **15**, 15, 15 | 15 |
| 4 | 1, 35, 91, **105**, 105, 105 | 105 |

Exact over $\mathbb Q$ for $k\le3$; mod $p=2^{61}-1$ for $k=4$ (so the *full-rank* entries there
are rigorous, since $\operatorname{rank}_{\mathbb Q}\ge\operatorname{rank}_{\mathbb F_p}$ for an
integer matrix; the rank-deficient entries at $k=4$ are only upper bounds and are not used).
Threshold matches $n\ge k$ in every case.

---

## 2. [V] Formalisation cross-validated against the project's own code

`work/bridge.py`: for every one of the 45 admissible matchings of profile $(k)=(2,1,1)$,
$d=0$, in $n=3$, with a random **cubic** rational field $X$,
$$\langle u,\ F(\gamma)(X)(0)\rangle\ \ \text{computed by }\texttt{verify/forests.py::Forest.F}
\ \ =\ \ \Phi_M\big(\textstyle\bigotimes_\nu c_{k_\nu}\otimes u\big),\qquad c_k=\partial^kX(0),$$
exactly over $\mathbb Q$, in all 45 cases. So the "slot/matching functional" picture used in
the proof is the same object as the project's elementary differential. (Cubic field, per the
CLAUDE.md rule; a quadratic field would degenerate the valence-3 slots.)

## 3. [V] Lemma L verified end-to-end on four profiles, including a $d=2$ one

`work/orbits.py` computes, for a fixed valence profile $(k_1,\dots,k_m)$ and $v$-degree $d$:
(i) the admissible matchings of the slot set $H$; (ii) their $\tilde G$-orbits, $\tilde G$ =
(in-slot permutations) $\rtimes$ (permutations of equal-valence nodes) $\times\ S_d$; (iii) the
rank, per dimension $n$, of the restricted functionals $\Phi_M|_W$ — i.e. of the elementary
differentials as functions of the free jet data.

| profile | $d$ | $\lvert H\rvert$ | admissible $M$ | orbits = #forests | rank at $n=1,2,3,\dots$ |
|---|---|---|---|---|---|
| $(1,1)$ | 1 | 6 | 4 | **2** | 1, **2**, 2 |
| $(2,1,1)$ | 0 | 8 | 45 | **16** | 1, 11, **16**, 16 |
| $(2,1,1)$ | 2 | 10 | 108 | **19** | 1, 16, **19**, 19, 19 |
| $(3,1)$ | 1 | 8 | 24 | **6** | 1, **6**, 6, 6 |

Exact over $\mathbb Q$ for $|H|\le8$, mod $p=2^{61}-1$ otherwise (full-rank entries rigorous
either way). Rank $=$ number of iso-classes in every case once $n$ is large enough, and
**strictly smaller for small $n$** — so the dimension hypothesis in Lemma L is not decorative.
Note the third row is a genuine $v$-degree-2 test: 19 pairwise distinct $v$-decorated forests,
independent from $n=3$ on.

---

## 4. [P] **Lemma L, proved.**

> **Theorem.** Let $\mathcal F$ be a finite set of pairwise non-isomorphic $v$-decorated exotic
> aromatic forests, and put
> $$n_0:=\max_{\hat\gamma\in\mathcal F}\ \tfrac12\lvert H(\hat\gamma)\rvert
> =\max_{\hat\gamma}\ \big(\#\mathrm{arrows}+\#\mathrm{lianas}+\#\mathrm{stolons}+d(\hat\gamma)+1\big).$$
> Then for every $n\ge n_0$ the maps
> $\mathfrak X(\mathbb R^n)\times\mathbb R^n\times\mathbb R^n\to\mathbb R^n$,
> $(X,x,v)\mapsto F_n(\hat\gamma)(X)(x)[v^{\otimes d(\hat\gamma)}]$, $\hat\gamma\in\mathcal F$,
> are linearly independent over $\mathbb R$. Equivalently, $F_n$ is injective on the span of the
> $v$-decorated exotic aromatic forests of total size $\le n$. No restriction on $d$; $d=2$ is
> the case needed.

The only external ingredient is Proposition B of §1. Everything else is proved here.

### Step 0 — reduction to free jet data at one point, one profile at a time

Assume $\sum_{\hat\gamma}b(\hat\gamma)F_n(\hat\gamma)(X)(x)[v^{\otimes d}]=0$ for all $X,x,v$.
$F_n(\hat\gamma)(X)(x)$ depends on $X$ only through its jet at $x$; translating, take $x=0$;
pairing with $u\in\mathbb R^n$ makes the identity scalar. Every tuple of tensors
$c_k\in V\otimes\operatorname{Sym}^k V^*$ ($k\le K$, $V=\mathbb R^n$) is realised as
$\big(\partial^kX(0)\big)_k$ by the polynomial field $X(y)=\sum_k\frac1{k!}c_k[y^{\otimes k}]$,
so the identity is a polynomial identity in the **free** variables $(c_k)_k,v,u$ over the
infinite field $\mathbb R$, hence an identity of coefficients. Now
$F_n(\hat\gamma)(X)(0)[v^{\otimes d}]$ is multihomogeneous: degree
$m_k(\hat\gamma)=\#\{\nu\in V(\hat\gamma):k_\nu=k\}$ in $c_k$, degree $d(\hat\gamma)$ in $v$,
degree 1 in $u$. Distinct multidegrees separate, so **we may assume all $\hat\gamma\in\mathcal F$
share one valence profile $\mathbf m=(m_k)_k$ (so one node count $m=\sum_km_k$) and one $d$.**

### Step 1 — the slot set

Fix node labels $V=\{1,\dots,m\}$ and valences $k_1,\dots,k_m$ realising $\mathbf m$, and put
$$H:=\{u\}\sqcup\{o_t\}_{t\in V}\sqcup\{\ell_{t,s}\}_{t\in V,\,1\le s\le k_t}\sqcup\{w_1,\dots,w_d\},
\qquad \lvert H\rvert=1+m+\textstyle\sum_tk_t+d=2N,$$
$N=1+d+\#\mathrm{arrows}+\#\mathrm{lianas}+\#\mathrm{stolons}$ (count upper slots
$m=1+a+2s$ and lower slots $\sum k_t=a+2l+d$).

Call a perfect matching $M$ of $H$ **admissible** if $u$ is matched to some $o_t$ and each $w_j$
to some $\ell_{t,s}$. To an admissible $M$ attach the $v$-decorated forest $\gamma_M$: nodes $V$,
node $t$ with in-slots $\ell_{t,\ast}$, root the node matched to $u$, decoration set
$\{t:\ \exists s,j,\ \{\ell_{t,s},w_j\}\in M\}$, and remaining pairs read as arrows ($o$–$\ell$),
lianas ($\ell$–$\ell$), stolons ($o$–$o$). Unwinding the definitions (and cross-checked against
`Forest.F` in §2):
$$\big\langle u,\,F_n(\gamma_M)(X)(0)[v^{\otimes d}]\big\rangle
=\Phi_M\Big(\bigotimes_{t\in V}c_{k_t}\otimes v^{\otimes d}\otimes u\Big)\tag{$*$}$$
with $o_t$ carrying the upper index of $c_{k_t}$, $\ell_{t,s}$ its $s$-th lower index
($\operatorname{Sym}^kV^*\cong\operatorname{Sym}^kV$ via $\delta$), $w_j$ carrying $v$, $u$ carrying $u$.

### Step 2 — the symmetry group is *exactly* forest isomorphism

Let $\ \tilde G:=\Big(\prod_k\big(S_{m_k}\ltimes(S_k)^{m_k}\big)\Big)\times S_d\ \le\operatorname{Sym}(H)$
act by permuting equal-valence node blocks, permuting the in-slots inside a node, and permuting
the $w_j$, fixing $u$.

**Claim.** $M\mapsto\gamma_M$ induces a *bijection*
$$\{\text{admissible matchings of }H\}/\tilde G\ \xrightarrow{\ \sim\ }\
\{\text{iso-classes of }v\text{-decorated exotic aromatic forests of profile }\mathbf m,\ v\text{-degree }d\}.$$

*Well defined and surjective:* $\tilde G$ preserves admissibility and $\gamma_{gM}\cong\gamma_M$;
every such forest is $\gamma_M$ for some $M$ by choosing a labelling of its nodes compatible with
valences. *Injective:* an isomorphism $\gamma_M\to\gamma_{M'}$ is a bijection $\alpha$ of $V$ with
$k_{\alpha(t)}=k_t$ — i.e. an element of $\prod_kS_{m_k}$ — together with bijections
$\operatorname{In}_t\to\operatorname{In}_{\alpha(t)}$ — an element of $(S_k)^{m_k}$ — carrying
matched pairs to matched pairs, root to root and decorations to decorations; the induced bijection
on $v$-slots is an element of $S_d$. The resulting $g\in\tilde G$ satisfies $gM=M'$. $\square$

This is the hinge, and it is where the CLAUDE.md rule *"canonicalise before testing
independence"* is discharged: $\tilde G$ is neither too small (which would split one forest into
several orbits and manufacture a false relation) nor too large (which would merge two forests).

### Step 3 — restriction to the jet subspace is an orbit sum

Let $s=\frac1{\lvert\tilde G\rvert}\sum_{g\in\tilde G}g$ act on $V^{\otimes H}$ by permuting
tensor factors. $s$ is idempotent, $V^{\otimes H}=\operatorname{im}s\oplus\ker s$, and
$$W:=\operatorname{im}s=(V^{\otimes H})^{\tilde G}
=\Big(\bigotimes_k\operatorname{Sym}^{m_k}\!\big(V\otimes\operatorname{Sym}^kV\big)\Big)
\otimes\operatorname{Sym}^dV\otimes V .$$
The tensors on the right-hand side of $(*)$ are $\bigotimes_kc_k^{\otimes m_k}\otimes v^{\otimes d}\otimes u$;
as $(c_k),v,u$ vary freely these **span** $W$, because $z^{\otimes r}$ spans $\operatorname{Sym}^r$
in characteristic $0$ (applied to $v\in V$ and to $c_k\in V\otimes\operatorname{Sym}^kV$).

For any $g\in\operatorname{Sym}(H)$, $\ \Phi_M\circ g=\Phi_{g^{-1}M}$ (immediate from the
definition of $\Phi$), whence
$$\Phi_M\circ s=\frac1{\lvert\tilde G\rvert}\sum_{g\in\tilde G}\Phi_{gM}
=\frac{\lvert\operatorname{Stab}_{\tilde G}(M)\rvert}{\lvert\tilde G\rvert}\sum_{M'\in\tilde G\cdot M}\Phi_{M'} .$$
Since $s|_W=\operatorname{id}$ and $s$ kills a complement, a combination $\sum_Mc_M\Phi_M$
vanishes on $W$ **iff** $\sum_Mc_M(\Phi_M\circ s)=0$ on all of $V^{\otimes H}$.

### Step 4 — conclusion

Pick for each $\hat\gamma\in\mathcal F$ an admissible $M_{\hat\gamma}$ with
$\gamma_{M_{\hat\gamma}}\cong\hat\gamma$. By Step 0 and $(*)$, $\sum_{\hat\gamma}b(\hat\gamma)\Phi_{M_{\hat\gamma}}$
vanishes on a spanning set of $W$, hence on $W$; by Step 3,
$$\sum_{\hat\gamma\in\mathcal F}b(\hat\gamma)\,
\frac{\lvert\operatorname{Stab}(M_{\hat\gamma})\rvert}{\lvert\tilde G\rvert}
\sum_{M'\in\tilde G\cdot M_{\hat\gamma}}\Phi_{M'}=0\quad\text{in }(V^{\otimes H})^*.$$
By Step 2 the orbits $\tilde G\cdot M_{\hat\gamma}$ are pairwise **disjoint** (distinct forests,
distinct orbits), so this is a combination of *distinct* matchings with the coefficient of each
$M'\in\tilde G\cdot M_{\hat\gamma}$ equal to $b(\hat\gamma)\lvert\operatorname{Stab}\rvert/\lvert\tilde G\rvert$.
Since $n\ge n_0\ge\lvert H\rvert/2$, Proposition B says the $\Phi_{M'}$ are linearly independent;
therefore $b(\hat\gamma)=0$ for every $\hat\gamma$. $\blacksquare$

**Corollary (what the target theorem needs).** At base $v$-degree $2$, the residual class —
connected exotic aromatic forests carrying a stolon with $N_0=1$ — contributes $v$-decorated
differentials that are linearly independent, so the $v$-degree-2 equation forces each coefficient
$b(\gamma)$ separately to vanish. Combined with the triangularity already established (erasing
$v$-leaves recovers $\gamma$, so distinct $\gamma$ contribute disjoint decorated families), the
residual class is killed and the $O(n)$ characterisation closes.


## 4b. [V] Adversarial check of Step 2 of the proof, its hinge

`work/step2.py`. The rank test of §3 alone cannot detect a *too large* symmetry group (a group
that wrongly fuses two forests still yields independent orbit sums). So Step 2 was checked
directly: partition the admissible matchings (A) by $\tilde G$-orbit, and (B) by the actual value
of $\langle u,F(\gamma_M)(X)(0)[v^{\otimes d}]\rangle$ on three random **cubic** rational fields.
A $\tilde G$ that is too large fuses matchings with different values; one that is too small splits
matchings with equal values. The two partitions must coincide — and do:

| profile | $d$ | matchings | $\tilde G$-orbits | distinct value-tuples | evaluator |
|---|---|---|---|---|---|
| $(2,1,1)$ | 0 | 45 | 16 | 16 | `verify/forests.py::Forest.F` |
| $(3,1)$ | 1 | 24 | 6 | 6 | $\Phi_M$ (validated in §2) |
| $(2,1,1)$ | 2 | 108 | 19 | 19 | $\Phi_M$ |
| $(2,1,0)$ | 1 | 27 | 15 | 15 | $\Phi_M$ |

All exact over $\mathbb Q$, $n=4$. The first row runs entirely through the project's own,
independently written elementary-differential code.

---

## 5. Provenance, and how the two recorded traps are avoided

**What is taken from the literature.** Only Proposition B (§1) — a classical fact about
$O(n)$/Brauer diagrams, *not* an LMK result, and verified here. In particular the proof uses
only the **easy half** of invariant theory (contractions are independent in high dimension);
it never needs the first fundamental theorem (that contractions *span* the invariants), and it
never needs $O(n)$-equivariance of anything. Nothing is quoted from LMK beyond the vocabulary.

**This is not LMK's proof.** The cached quotes do not contain the statement of Prop. 4.1 or the
construction of the $\theta$-parametrised dual fields, and the papers were not re-fetched. The
argument above is independent, and yields the two-coloured/decorated statement directly rather
than by transcription. It implies the cited consequence "$F$ is injective on
$\operatorname{Span}(\Gamma)$" for $n$ large.

**Trap 1 — "the $\theta$'s are mandatory."** Respected, and side-stepped rather than repaired.
As described in the brief (this is a reconstruction of the $\theta$-route, not a quotation — the
cached literature does not carry LMK's construction), that route builds for each $\gamma$ a single
*witness* field $f_\gamma$ in dimension $\lvert\gamma\rvert+1$ meant to separate $\gamma$ from every
other forest, and LMK's recorded warning is that the $\theta$-free witness fails to separate some
pair. The proof above needs **no witness**: it tests against the whole jet space $W$, and all the
separating work is done by Step 2 — that $\tilde G$-orbits are *exactly* iso-classes. In this
language the failure mode is an implicit enlargement of the acting group: a group properly
containing $\tilde G$ may fuse two orbits, and any such fusion would "prove" a false independence
statement (Trap 2 below is exactly such a fusion, exhibited). So the two accounts agree on where
the danger lies; carrying the $\theta$-parameters is one way of keeping nodes distinguishable, and
pinning $\tilde G$ down exactly is another.
*Equivalent variant.* Substituting $X=\sum_{t=1}^m\lambda_tZ_t$ and extracting the coefficient of
$\lambda_1\cdots\lambda_m$ (full polarisation) replaces $\tilde G$ by the smaller
$G=\prod_tS_{k_t}\times S_d$, at the cost of an extra bookkeeping step
($\sum_\beta$ over node labellings, each forest occurring $\lvert\operatorname{Aut}\rvert$ times).
This is the precise analogue of the $\theta$'s; the proof above avoids needing it.

**Trap 2 — "$v:=X(x)$ fails."** Confirmed, and located exactly (`work/collapse.py`). In the slot
picture the substitution turns each $v$-slot $w_j$ into the upper slot of an extra valence-$0$
node, so the group jumps from $S_{m_0}\times S_d$ to $S_{m_0+d}$ ($m_0$ = number of
derivative-free nodes already present). Orbits therefore merge **iff $m_0\ge1$**:

| profile | $d$ | decorated forests | distinct images after $v:=X(x)$ |
|---|---|---|---|
| $(2,1,1)$ | 2 | 19 | 19 — no collapse ($m_0=0$) |
| $(2,1,0)$ | 1 | 15 | **13** |
| $(2,2,0)$ | 2 | 12 | **10** |
| $(1,1,0,0)$ | 1 | 7 | **6** |

The first row explains why the trap is easy to miss: on a forest with no $X$-leaf the shortcut
looks sound. The general residual class has $X$-leaves, so it is not.

**Dimension.** The proof gives $n\ge n_0=$ (number of contracted pairs). Not sharp: §3 finds full
rank already at $n=3$ for profile $(2,1,1)$, $d=2$, where $n_0=5$. **[C] Conjecture** (not needed,
not proved): $n\ge\lvert V\rvert+1$ suffices, which is the dimension in which a $\theta$-style dual
field would live. Lemma L is stated "for $n$ large", so the gap is immaterial; but a
*fixed-dimension* version of the $O(n)$ theorem would need it.

## 6. Two remarks on scope

* **The $Y$-leaf reduction is not needed.** The brief's reduction (attach a fresh constant-$Y$
  leaf at each $v$-slot, then apply a two-colour Prop. 4.1) is correct but can be skipped: the
  proof handles the fixed vector $v$ directly, via the factor $\operatorname{Sym}^dV$ in $W$. The
  same proof also gives the **general multi-coloured** statement — several distinct fields
  $X^{(1)},\dots,X^{(c)}$, colours arbitrary (not confined to leaves) — by letting $\tilde G$
  permute only nodes of equal *(colour, valence)*; Step 2 goes through verbatim.
* **Interface with the application (worth checking when the lemma is used).** The theorem
  hypothesises a relation valid for *all* $X\in\mathfrak X(\mathbb R^n)$ and all $v$. This is what
  the $O(n)$ argument supplies: although (T) lives on the thin locus $\{TX\}\subset\mathfrak
  X(\mathbb R^{2n})$, expanding the defect at $(x,v)$ re-expresses it in elementary differentials
  of an *arbitrary* base field $X$ on $\mathbb R^n$, decorated with $v$'s. So the hypothesis is met.
  It would **not** be met by a relation known only for $X$ in a proper subclass (e.g. gradient
  fields), and the proof would then fail at Step 0.

## 7. Reproducing the computations

`work/run_lemmaL_checks.py` runs all five (`sft`, `bridge`, `orbits`, `step2`,
`collapse`); no
dependencies; prints `ALL LEMMA-L CHECKS PASSED`. `bridge.py` and `step2.py` read
`/home/user/claude_tests/verify/{poly,forests}.py` and does not modify the repository. The four
modules are in drop-in shape for a `verify/test_lemmaL.py` should §11 of the write-up be updated.
