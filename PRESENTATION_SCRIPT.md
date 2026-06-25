# Presentation Script — Three Rivers Metagenomics (4 people, ~7 minutes)

Structure follows the rubric: **the question → what we did → what we found → why we believe
it → what we're unsure about.** Each person owns one block. Bold = the number/idea to land.
Times are targets; total ≈ 7 min, leaving room for Q&A (the part that's actually graded).

> Slide order matches the plots in `Output/`. Speak to evidence, not to slides.

---

## 👤 Person 1 — The question, the data, the approach  (~1:30)
*(Slides: title + one data slide)*

> "Our question is simple to state and hard to answer: **how do bacterial communities differ
> across Pittsburgh's three rivers?** We have about **49 water samples** from 2019 — five river
> sites plus a distilled-water control, across four seasons. Each sample is a file of DNA
> reads from the bacterial 16S barcode gene; **identical reads mean the same organism**, which
> we call an OTU, so each sample becomes a frequency map — which bacteria, and how many.
>
> From those maps we measure two things. **Alpha diversity** — how diverse *one* sample is:
> its richness and its evenness. And **beta diversity** — how *different two* samples are,
> using Bray-Curtis distance. The plan: measure all 49 samples, then ask whether **site** or
> **season** explains what we see — and use the control to check ourselves."

**Hand off:** "Sarah will start with what one sample looks like — alpha diversity."

---

## 👤 Person 2 — Alpha diversity + the control sanity check  (~1:45)
*(Slides: plot 1, then plot 2)*

> "First, alpha diversity. **Richness** is how many distinct OTUs a sample has; **evenness** —
> the Gini-Simpson index, one minus Simpson's D — is whether it's a balanced mix or dominated
> by a few. Here's the key disappointment: **the sites overlap**. The within-site spread is
> bigger than the between-site difference, and richness and evenness barely correlate
> (**r ≈ −0.01**) — so a rich sample isn't necessarily an even one. **Alpha diversity alone
> can't tell our rivers apart.**
>
> Now the control. Distilled water should be nearly empty. But by **richness** it doesn't look
> special — Fall controls were among our *richest* samples, because with exact-match OTUs,
> sequencing noise invents fake species. So does our pipeline fail? No — switch to
> **composition**: the control's mean Bray-Curtis distance to river samples is **0.92**, versus
> **0.78** among river samples. **Beta diversity catches the control that alpha missed.** Our
> sanity check passes, and the lesson is: always cross-check the two."

**Hand off:** "So composition is where the signal is. Marcus takes it from here."

---

## 👤 Person 3 — Beta diversity: site vs season  (~1:45)  ← the centerpiece
*(Slides: plot 3, then plot 5 PCoA, then plot 4 heatmap)*

> "We compared every pair of river samples and sorted them into four groups. Two replicates —
> same site, same season — differ by **0.42**; that's our noise floor. Now change **only the
> site**: it rises to **0.54** — a small bump of **+0.11**. Change **only the season**: it
> jumps to **0.88** — a bump of **+0.45**. Above the replicate baseline, **season moves
> communities about four times more than site does.**
>
> You can *see* it without any of our bucketing. This is a PCoA — a map where nearby points
> are similar communities. **Same points, two colorings.** By **season**, four clean clusters —
> Spring, Summer, Fall, Winter — and the controls fly off to the side. By **site**, the exact
> same points are a random mix; no site clusters. And the heatmap agrees: the dark, similar
> blocks line up **by season**. Three independent views, one answer: **season, not site,
> organizes these communities.**"

**Hand off:** "But a good result has to survive scrutiny — Priya on why we believe it."

---

## 👤 Person 4 — Why we believe it + what we're unsure about  (~1:30)
*(Slides: plot 8 PERMANOVA, then plot 7, then a limitations slide)*

> "Is this real or just a suggestive average? We ran a proper test — a **PERMANOVA**, which
> treats the 44 samples as the units and asks how much variation each factor explains.
> **Season explains 59% of the variation, p = 0.001. Site explains 8%, and is no better than
> random — p = 0.68.** So statistically, season is the driver and site is undetectable.
>
> And it's robust three ways: dropping each sample, season wins **44/44**; dropping a whole
> season or a whole site, it still wins; and under a completely **different distance metric**
> — presence/absence Jaccard — season still wins by 3.6×.
>
> **What we're *not* sure about** — because honesty is the point. Richness partly just tracks
> sequencing depth (**r ≈ 0.56**), so we lean on composition, not raw counts — and we haven't
> rarefied to fully rule depth out. Our OTUs are exact-match, which inflates the control's
> apparent richness. And 'season' bundles temperature, runoff, and nutrients — we show the
> *pattern*, not the *mechanism*.
>
> **Bottom line: across the Three Rivers in 2019, season — not site — was the dominant driver
> of bacterial community composition, our control confirms the method is sound, and the result
> survives a formal test, a metric swap, and dropping any group.** Thank you — questions?"

---

## 30-second team-wide summary (everyone should be able to say this)
> "Reads → frequency maps → alpha diversity (within a sample) and beta diversity (between
> samples). Alpha can't separate the sites and doesn't flag the control; **beta does both**.
> **PERMANOVA: season R² = 0.59, p = 0.001; site R² = 0.08, p = 0.68** — season is the driver,
> site is undetectable. It's visible as **seasonal clusters in PCoA**, and it **survives
> leave-one-group-out and an independent presence/absence metric**. What we're still unsure
> about: depth-confounded richness (not yet rarefied), exact-match OTUs, and that 'season' is
> a proxy for temperature/runoff/nutrients (and possibly sequencing batch)."

---

## Optional: "The lost cities" analogy (a *going-farther* prompt)
The rubric mentions a "lost cities" analogy from your bootcamp. The natural mapping: you can
**compare places you can never fully visit by sampling the fragments they leave behind** — for
lost cities, artifacts; for our rivers, **DNA reads**. A diversity metric is how you compare
two such places *quantitatively* without a full census: richness = how many kinds of artifact,
beta diversity = how different two sites' artifact-piles are. **⚠️ Confirm the exact framing
against your bootcamp notes before using it — this is a reconstruction, not a quote.**

---

## Q&A DEFENSE BANK

Built from an independent verification (which reproduced every number and found **no code
bugs**) plus a four-persona examiner panel. Ordered by how likely you are to get asked. Each
has a full answer and a **one-line version to memorize**. Golden rule: **concede the real
limit, then show why the headline still stands.**

> Before anything: the headline rests on **PERMANOVA (season R²=0.59 p=0.001 vs site R²=0.08
> p=0.68)** and the full distance matrix — *not* on the PCoA picture or the raw "4×". If a
> question attacks a soft spot, retreat to PERMANOVA + leave-one-group-out.

### 🔴 Highest-probability questions

**1. (Stats) "Your 4× comes from 946 pairwise distances, but you only have 44 independent
samples — each is in 43 pairs. Isn't that pseudoreplication?"**
Yes, the pairs aren't independent, which is exactly why we never put a p-value on those 946
numbers — we used them only as a descriptive effect size. The real test treats the 44 samples
as units: **PERMANOVA, which gives season R²=0.59, p=0.001 and site R²=0.08, p=0.675.** That
respects the dependence structure and still says season dominates.
> *Memorize:* "Pairs aren't independent — so we ran PERMANOVA on the 44 samples: season R²=0.59 p=0.001, site p=0.68."

**2. (Robustness) "Depth is collinear with season — Spring ~3300 reads, Winter ~1560. Isn't
the season effect partly a sequencing-depth artifact?"**
This is our most honest open weakness. Bray-Curtis is more depth-robust than richness, which
is why we lead with composition — but with exact-match OTUs, two identical communities at
different depths still share a smaller fraction of exact sequences, and depth correlates with
season. The clean fix is **rarefying to a common depth and re-running; we haven't done that
yet**, so we say the effect is consistent with season but not fully de-confounded from depth.
> *Memorize:* "Honest gap — depth tracks season; only rarefaction fully separates them, and we haven't rarefied yet."

**3. (Metric) "You call it 'metric-independent' because BC and Jaccard correlate at r=0.995 —
but your Jaccard is abundance-weighted Ružička, built from the same min-overlap as BC. Isn't
that tautological?"**
You're right that r=0.995 is near-tautological — Ružička and BC share the `sum_of_minima`
numerator (it's even provable that Ružička ≥ BC always). So we *don't* lean on that. The
genuinely independent test is **classic presence/absence Jaccard, which ignores abundance
entirely — and season still beats site 3.6× under it.** That's the real metric-independence.
> *Memorize:* "BC vs Ružička r=0.995 is tautological; the real test is presence/absence Jaccard — season still wins 3.6×."

**4. (Biology) "Walk me through what one sequencing error does under your exact byte-for-byte
OTU definition."**
Each error turns one true read into a brand-new distinct string, so it **invents a fake
species** — that's why richness tracks read depth (r=0.56) and why the control looks "rich."
For Bray-Curtis, an error sequence is absent in the other sample, so it never enters the
min-overlap and only inflates distance. That inflation is roughly a constant offset that
raises our 0.42 replicate floor; since the headline is a *difference above that floor*, the
shared offset largely cancels — and we'd confirm with 97%/ASV clustering.
> *Memorize:* "Each error mints a fake OTU — inflates richness and the BC floor, but the headline is a difference above the floor, so it largely cancels."

**5. (Stats) "Your design is unbalanced — Fall 19, Spring 7 — and your buckets are unweighted
means. How do you know the 4× isn't just Fall over-representation?"**
The bucket means are vulnerable to that, true. Two defenses: the gap is huge (season +0.45 vs
site +0.11 above baseline), so re-weighting can't plausibly flip it, and **PERMANOVA already
accounts for the group structure** and still gives season R²=0.59 vs site 0.08. Leave-one-
season-out (4/4) confirms no single season — including Fall — carries the result.
> *Memorize:* "PERMANOVA handles the imbalance and still says R²=0.59 vs 0.08; leave-one-season-out is 4/4."

**6. (Metric) "Your PCoA reports 36% and 13% by dividing by positive eigenvalues only. How
much is in the negatives, and is 49% in two axes enough to claim clusters?"**
Reporting against positive eigenvalues is the standard classical-MDS convention; in our data
there's only one eigenvalue at ≈ −1e-16, so Bray-Curtis is essentially Euclidean here and the
distortion is negligible. But you're right that an eyeball cluster on 49% isn't proof — which
is why **the claim rests on the PERMANOVA and the full-matrix category means, not the 2-D
picture.** The PCoA is illustration; the statistics are the evidence.
> *Memorize:* "Only one ≈−1e-16 eigenvalue, so it's near-Euclidean; but the proof is PERMANOVA + the full matrix, not the 2-D plot."

**7. (Biology) "Your control has richness 104–136 — richer than most river samples. Distilled
water should be near-sterile. Does that prove contamination is *absent*, or the opposite?"**
The opposite, and we're careful with the word "flag." A water blank with 100+ sequences is a
loud **contamination/error floor** — reagent 'kitome' bacteria plus sequencing noise — running
through every sample. Beta diversity shows the control is *compositionally distinct* (0.92 vs
0.78), i.e. it has its own contaminant signature rather than being empty. The right next step
is a `decontam`-style blocklist subtraction, which we haven't done. The lesson stands: a
high-richness blank is dangerous, and **alpha diversity can't see it — beta can.**
> *Memorize:* "A 100+ richness blank is a contamination signal, not a clean negative; beta sees it's distinct, but we should decontam-subtract it."

### 🟠 Likely code / metric questions

**8. (Code) "In `betaDiversityMatrix` you loop `for j in range(i, n)` and write both `[i][j]`
and `[j][i]`. Why start at `i`, and how do you *know* mirroring is safe?"**
Both metrics are built from `sum_of_minima`/`sum_of_maxima`, which are symmetric in the two
samples, so d(A,B)=d(B,A) — mirroring is a correctness fact, not just a speed trick (it halves
~1,200 distance calls to ~600). Starting `j` at `i` (not `i+1`) includes the diagonal, which
computes d(sample, sample)=0. The one precondition is no empty sample, since d(A,A) on an
empty map would raise.
> *Memorize:* "Both metrics are provably symmetric, so we mirror the upper triangle; j starts at i to fill the zero diagonal."

**9. (Metric) "`simpsons_index` returns D but `simpsons_map` returns 1−D, and your CSV has
both. If `simpson_D`=0.8 for one sample and `gini_simpson`=0.8 for another, which is more
diverse?"**
The `gini_simpson`=0.8 sample. `D`=0.8 means 80% chance two random reads are the *same*
species — high dominance, *low* diversity. `gini_simpson`=1−D=0.8 means 80% chance they're
*different* — high evenness, *high* diversity. Plot 1's right panel plots the 1−D values, and
we verified `simpson_D + gini_simpson = 1` for all 49 rows. (Caveat: Gini-Simpson rises with
both richness and evenness, so "evenness" is a loose label.)
> *Memorize:* "D is dominance (high=less diverse); Gini=1−D is what we plot (high=more diverse); they sum to 1 for all 49 rows."

**10. (Stats) "Why compute the 'different site, different season' bucket at all? It's 0.88 —
identical to the season-only bucket."**
That equality is the point: once you've already changed the season, *additionally* changing
the site adds essentially nothing — **season saturates the distance**, which *supports*
season-dominance. Site has a small standalone effect (+0.11) but ~0 on top of a season change.
(One caveat: equal 0.88s could partly reflect a Bray-Curtis ceiling near 1.)
> *Memorize:* "Season-only ≈ both ≈ 0.88 means site adds nothing once season changes — season saturates the signal."

**11. (Robustness) "Dropping 1 of 44 samples barely moves a bucket of hundreds of pairs.
Why does 44/44 convince me?"**
On its own it doesn't — sample-level leave-one-out is near-trivial here, it only rules out a
single freak outlier. That's why we also did **leave-one-GROUP-out: drop a whole season (4/4
hold) and a whole site (5/5 hold)**, plus PERMANOVA. Those are the real structural tests; the
44/44 is just the rubric's literal "one sample" question answered.
> *Memorize:* "Sample-LOO is weak — so we also drop whole seasons (4/4) and whole sites (5/5), and ran PERMANOVA."

### 🟡 Limits & framing

**12. (Limits) "Site and season are partially confounded — not every site is sampled every
season. Can you actually separate them, or are they aliased?"**
Partially aliased, yes — an incomplete season×site grid can't be fully orthogonalized. The
mitigation: our 'different-site/same-season' bucket holds season fixed while varying site (and
vice-versa), and PERMANOVA estimates each factor's R² given the structure. What we can't fully
rule out is that the specific site-pairs co-occurring within a season happen to be the similar
ones — a two-way/nested PERMANOVA with a design-rank check is the proper fix.
> *Memorize:* "Partly aliased; each effect is estimated with the other held fixed in-bucket, but full de-aliasing needs a two-way model."

**13. (Limits) "What is 'season' actually a proxy for, and can you rule out a batch effect?"**
Season bundles water temperature, rainfall/runoff and flow, and nutrient pulses — we have no
environmental covariates to separate them, so our claim is correlational: composition tracks
*sampling time*. The confound we most want to name is **batch**: if all Fall samples were
extracted and sequenced together, a run effect is perfectly aliased with season. We'd need
extraction/run metadata (with controls in every run) to rule it out.
> *Memorize:* "Season = temperature+flow+nutrients, and possibly sequencing batch; with no metadata our claim is correlational, not mechanistic."

**14. (Limits) "Your 0.42 replicate baseline is mostly Fall/Winter pairs. Aren't you
subtracting a baseline measured in the wrong seasons?"**
The replicate bucket is unbalanced — Spring/Summer have many n=1 cells contributing zero
replicate pairs — so 0.42 is a rough, Fall/Winter-weighted noise floor, not a per-season
value. That mainly weakens the *site* comparison, not the headline: even a higher true
baseline leaves season (+0.45) dwarfing site (+0.11). And PERMANOVA doesn't use a baseline at
all, yet agrees.
> *Memorize:* "0.42 is Fall/Winter-heavy and mainly biases the site comparison; PERMANOVA uses no baseline and still says season wins."

---

### If you only memorize three rebuttals
- **Pseudoreplication / "is it real?"** → "PERMANOVA: season R²=0.59 p=0.001, site p=0.68."
- **"Metric-independent?"** → "Presence/absence Jaccard (shares nothing with BC) — season still 3.6×."
- **"Throw out a sample?"** → "44/44 samples, 4/4 seasons, 5/5 sites — nothing drives it."
