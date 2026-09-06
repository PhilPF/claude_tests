# Cached literature

**Purpose: stop re-fetching.** Several agents independently downloaded the same papers and
re-read the same definitions, at roughly 10–15k tokens per fetch. Everything below was read
from the source in this project; quote from here instead of fetching, and fetch only to settle
a *disputed* citation or to look up something genuinely not recorded here.

Provenance is marked. **[V]** = verified verbatim against the source in this project.
**[R]** = reported by an agent, not independently confirmed — treat as a lead, not a citation.

---

## Laurent & Munthe-Kaas, *The universal equivariance properties of exotic aromatic B-series*
arXiv:2305.10993v2 (= Found. Comput. Math., 2024).

**[V] Def. 2.7 vocabulary — the terminology that was inverted twice in this project.**
> "If $\sigma(a_1)=a_2$, we say that the unordered pair $(a_1,a_2)$ is a **liana** and we
> represent it with a dashed edge between the two nodes $\tau(a_1)$ and $\tau(a_2)$, that can be
> identical. If $\sigma(v_1)=v_2$, we call the unordered pair $(v_1,v_2)$ a **stolon** and we
> draw it with a double edge between $v_1$ and $v_2$."

The $a_i$ are **arrows**, the $v_i$ are **nodes**. So liana = two arrows = double derivation;
stolon = two nodes = double evaluation.

**[V] Which is which analytically** (§2.2, preceding Def. 2.7):
> "The Laplacian is then represented with a double arrow (i.e., a double derivation), while the
> scalar product is given by two nodes identified through the source map $\sigma$."

So **liana = Laplacian** and **stolon = scalar product**.

**[V] Exotic vs stolonic trees:**
> "We define an **exotic tree** (respectively a **stolonic tree**) as an exotic aromatic tree
> that reduces to a standard Butcher tree when removing all the **lianas** (respectively by
> removing all the **stolons**)."

Hence "trees + Laplacians" = exotic trees = **exotic B-series**.

**[V] Prop. 4.3:**
> "Connected exotic aromatic B-series are decoupling. stolonic B-series are
> left-orthogonal-equivariant and exotic B-series are right-orthogonal-equivariant."

**[V] Thm 2.13 classification (Table 2).** orthogonal $\leftrightarrow$ exotic aromatic B-series;
GL $\leftrightarrow$ aromatic; left-orthogonal $\leftrightarrow$ stolonic; right-orthogonal
$\leftrightarrow$ exotic; affine $\equiv$ semi-orthogonal $\leftrightarrow$ B-series. Right-orthogonal
means affine $A$ with $AA^\top=I$; left-orthogonal means $A^\top A=I$.

**[V] End of the Thm 2.13 proof:**
> "If $\phi$ is semi-orthogonal-equivariant, then $\gamma$ is a linear combination of connected
> exotic aromatic trees without lianas, loops, and stolons, that is, a combination of standard
> Butcher trees."

**[V] Prop. 4.1 — independence, and why the $\theta$'s are mandatory:**
> "In particular, the elementary differential map $F$ is injective on $\mathrm{Span}(\Gamma)$."

The paper explicitly shows the $\theta$-free dual vector field is **insufficient**
("the dual vector field without the $\theta$ parameter fails to identify the difference between
$\gamma$ and $\hat\gamma$"), and notes this "reveals a typographical error in [27, Rk. 4.8] … and a
minor error in [37, Sec. 4.2]". **Any proof of the two-coloured analogue must carry the $\theta$'s.**

**[V] Dimension proviso** (proof of Thm 3.3): $\delta$ "is a bijection if and only if
$2d \ge |\kappa| + |\kappa'| + 1$" — so injectivity arguments must be run in high enough dimension.

**[V] Gradient fields** (§4.3, Prop. 4.5, Thm 4.6): on $\mathfrak X^\nabla$ the classification
collapses via edge–liana and edge–stolon inversions. *Note (this project): the gradient class is
not stable under $T$, so closure is not well posed there and Thm 4.6 cannot be combined with it.*

**[V] §5 Conclusion — the partitioned case is explicitly open:**
> "There exists a handful extensions of B-series used in numerical analysis such as partitioned
> B-series, exponential B-series or Lie-Butcher series, and a variety of equivariance properties
> … It would be interesting to draw geometric links between the different equivariance properties
> and the different types of B-series."

Table 2 has **no partitioned row**.

---

## McLachlan, Modin, Munthe-Kaas & Verdier, *B-series methods are exactly the affine equivariant methods*
Numer. Math. **133** (2016) 599–622; arXiv:1409.1019.

**[R] Def. 2.1** integrator map (sequence of smooth maps on compactly supported fields).
**[R] Def. 2.3** "B-series map" is a condition on the derivatives $D^k\phi_n(0)$ **at the zero
vector field only** — this is what makes the classification a *jet* theorem and admits the
flat-remainder counterexample.
**[R] Thm 2.4** local + affine equivariant $\iff$ B-series. **[R] Rmk 2.5** only surjective affine
maps plus the trivial injections are needed.

## Munthe-Kaas & Verdier, *Aromatic Butcher series*
Found. Comput. Math. **16** (2016) 183–215; arXiv:1308.5824. **[R]** Thm 2.4: local + affine
equivariant in fixed dimension $\Rightarrow$ aromatic B-series.

## Others
**[R]** Hager, *Numer. Math.* **87** (2000) 247–282 — discrete adjoint of an RK method.
**[R]** Sanz-Serna, *SIAM Rev.* **58** (2016) 3–33 — symplectic RK for adjoint equations.
**[R]** Kolář–Michor–Slovák, *Natural Operations in Differential Geometry*, 1993 — Weil functors
are exactly the product-preserving bundle functors; natural operators $\mathfrak X\to\mathfrak X$
are only $X\mapsto cX$.
