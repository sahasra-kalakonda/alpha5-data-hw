# Report 2 — The Math & Algorithms

Every metric, derived from scratch, with *why we chose it* and *how to defend it*. The
rubric explicitly warns you'll be asked "why this metric" — this report is your answer sheet.

Notation: a sample is a vector of counts `a = (a₁, a₂, …)` over OTUs; `N = Σ aᵢ` is its total
number of reads. Two samples are `a` and `b`.

---

## 1. Richness `S`
**Definition:** number of OTUs with `aᵢ > 0`. Code: `helperFunctions.richness`.

**Reading:** "how many kinds of bacteria did we see." Ignores how common each is.

**Why / caveats:** simplest alpha metric, but **biased by sequencing depth** — sequence more
reads and you discover more rare OTUs. In our data richness vs depth is `r ≈ 0.56` (plot 7).
With exact-match OTUs it's also inflated by sequencing errors. *Defense:* we report richness
but do **not** hang conclusions on it; we use evenness and composition, which are far less
depth-sensitive.

---

## 2. Simpson's index — dominance `D` and Gini-Simpson `1 − D`
**Definition:** `D = Σ (aᵢ/N)²`. Code: `helperFunctions.simpsons_index`.

**Probabilistic meaning:** `D` is the probability that two reads drawn at random are the
**same** OTU. So **high `D` = low diversity** (a few OTUs dominate). The complement
`1 − D` (Gini-Simpson) is the probability they're **different** OTUs — an **evenness**
measure where **higher = more even**. Code: `simpsonsMap.simpsons_map` returns `1 − D`.

**Why we report `1 − D`:** it reads intuitively ("higher = more diverse") and is bounded in
`[0, 1)`. **Why Simpson over Shannon:** Simpson weights common taxa and is robust to the long
tail of rare, error-prone singletons that exact-match OTUs produce — exactly our situation.

**Watch the wording.** The function `simpsons_index` returns `D`; the table/plots use `1 − D`.
Our `alpha_diversity.csv` writes **both** columns (`simpson_D`, `gini_simpson`) so there's no
confusion about which "Simpson's index" we mean. (This was a real ambiguity we fixed — see
`REPORT_4`.)

---

## 3. Bray-Curtis distance `BC`
Let `C = Σ min(aᵢ, bᵢ)` = shared abundance. Then
```
BC = 1 − 2C / (N₁ + N₂)
```
Code: `helperFunctions.bray_curtis_distance` (with `average(N₁,N₂) = (N₁+N₂)/2`, so
`C / average = 2C/(N₁+N₂)` ✓).

**Reading:** fraction of reads **not** shared. `0` = identical communities, `1` = no overlap.
Abundance-sensitive: dominant OTUs matter.

**Why it's our main metric:** it's the standard for microbiome composition, uses abundances
(not just presence/absence), and bounded in `[0,1]`. Our whole site-vs-season argument lives
in this matrix.

---

## 4. Jaccard distance — the abundance-weighted (Ružička) form we implemented
With `C = Σ min` and `Σmax = Σ max(aᵢ, bᵢ)` over the OTU **union**:
```
J = 1 − C / Σmax
```
Code: `helperFunctions.jaccard_distance`.

**Important:** this is **not** the textbook presence/absence Jaccard `1 − |A∩B|/|A∪B|`. Using
counts makes it the **abundance-weighted Jaccard = Ružička distance**. It collapses to the
classic set version only if every count is 1. *Be ready to say this out loud* — it's the
single most likely "gotcha" on metric choice.

**Relationship to Bray-Curtis (why they agree at r ≈ 0.99):** algebraically
```
J = 2·BC / (1 + BC)
```
a fixed, increasing curve. So `J ≥ BC` for every pair and the two are near-perfectly
correlated **by construction** — visible as the tight curve hugging above `y = x` in plot 6.
*Interpretation:* our conclusion is metric-robust, but plot 6 is **not** independent evidence;
it's a consistency check. (A genuinely independent contrast would need *presence/absence*
Jaccard.)

---

## 5. The beta-diversity matrix algorithm
`betaDiversityMatrix.beta_diversity_matrix`:
- sorts sample names for a reproducible, fixed ordering;
- selects the distance function by name (raises on anything unknown — no silent default);
- fills an `n×n` matrix computing only the **upper triangle** (`j ≥ i`) and mirroring, since
  `d(A,B) = d(B,A)` and `d(A,A) = 0`. This halves ~1,200 distance computations to ~600.

Complexity: `O(n² · U)` where `U` is OTUs per pair. For `n = 49` it's instant.

---

## 6. PCoA / classical MDS (`analyze.pcoa`)
Goal: place the 49 samples as points in 2-D so that plotted distances approximate the
Bray-Curtis distances — then we can *see* the structure (plot 5).

Algorithm (classical multidimensional scaling), done by hand with numpy:
1. `D2 = D²` element-wise.
2. Double-centre: `B = −½ · J · D2 · J`, with `J = I − (1/n)·11ᵀ`.
3. Eigen-decompose symmetric `B` (`numpy.linalg.eigh`).
4. Coordinates = top-2 eigenvectors × `√eigenvalue`; variance explained per axis =
   `eigenvalue / Σ(positive eigenvalues)`. Our axes capture **36%** and **13%**.

**Caveat to own:** Bray-Curtis is not perfectly Euclidean, so a few eigenvalues are slightly
negative; we use the positive ones for the variance fractions (standard practice). 36% on
axis 1 is normal for noisy microbiome data and is enough to reveal the seasonal clusters.

---

## 7. Quantifying "site vs season" honestly
Raw means: site comparison `0.54`, season comparison `0.88`. Naively that's only `1.6×`. We
instead measure **above the replicate baseline `0.42`** (the floor set by measurement noise),
which isolates the *added* difference each factor causes:
```
site effect   = 0.54 − 0.42 = 0.11
season effect = 0.88 − 0.42 = 0.45
ratio         = 0.45 / 0.11 ≈ 4×
```
This descriptive "~4×" is the *plain-English* effect size. The **formal** version is §8.

---

## 8. PERMANOVA — the proper statistical test
Comparing mean pairwise distances has a flaw: the pairs aren't independent (each sample is in
many pairs — *pseudoreplication*), so those `n=` counts are **not** real sample sizes. The
correct test treats the **44 samples** as the units. We implemented one-way **PERMANOVA**
(Anderson 2001) by hand (`analyze.permanova`):
```
SS_total  = (1/N) Σ_{i<j} d²ᵢⱼ
SS_within = Σ_groups (1/n_g) Σ_{i<j in g} d²ᵢⱼ
pseudo-F  = (SS_among/(a−1)) / (SS_within/(N−a)),   R² = SS_among / SS_total
```
with a p-value from **999 label permutations**. Results on the river samples (Bray-Curtis):

| Factor | R² (variation explained) | pseudo-F | p (999 perms) |
|---|---|---|---|
| **Season** | **0.59** | 19.2 | **0.001** |
| **Site** | 0.08 | 0.8 | 0.675 |

Season explains **~59%** of compositional variation and is **highly significant**; site explains
**~8%** and is **statistically indistinguishable from random** (p ≈ 0.68, pseudo-F < 1). This is
the rigorous form of "season ≫ site." *Caveat:* these are one-way tests per factor; because the
design is unbalanced and season dominates, the one-way **site** test understates site (season
variance lands in its "within" term). A two-way/nested PERMANOVA would partition cleaner — but
that can only *help* site, which is already the loser.

## 9. Robustness: three independent stress tests
1. **Leave-one-sample-out** — drop each of 44 samples; season beats site **44/44** (ratio
   3.7×–4.3×). *Weak alone* (dropping 1 of ~600 pairs barely moves a mean), so:
2. **Leave-one-group-out** — drop a **whole season** (holds 4/4) and a **whole site** (holds
   5/5) and the ordering survives. This is the real structural test.
3. **Presence/absence Jaccard** — classic binary Jaccard (`helperFunctions.jaccard_binary_distance`)
   shares **no numerator** with Bray-Curtis (unlike Ružička). Under it, season still beats site
   by **3.6×** — so the result is **genuinely metric-independent**, not a tautology of the shared
   min-overlap term.
