# Report 4 — Bugs Fixed, Improvements, and Limitations

What was wrong in the original code, what we changed, and what we'd still do with more time.
Being able to say "here's what was broken and here's how we know it's fixed" is itself a
strong defense.

---

## A. Bugs we found and fixed

All numbers below were verified by recomputing from the data.

| # | Bug | Severity | Fix |
|---|---|---|---|
| 1 | **Printed "conclusions" were hardcoded** and contradicted the data — e.g. it claimed Allegheny had the highest richness (~95) when Allegheny is actually the **lowest** river site (66); the ~95 was the **control**. Also wrong: "Neville most even" (it's Sharpsburg), "Simpson range 0.79–0.91" (really 0.87–0.97), "biggest gap Allegheny–Braddock" (it's Braddock–Sharpsburg), "most similar Neville–Sharpsburg" (it's Allegheny–Neville). | 🔴 High | Conclusions are now **computed from the data** with f-strings; nothing is hardcoded. Backed by `effect_summary.csv`. |
| 2 | **"Season explains 4× more variance"** — the raw numbers (0.54 vs 0.88) are only 1.6×, and "variance" implied a statistic never computed. | 🟠 Med | Reworded to **"~4× above the replicate baseline"** (0.45 vs 0.11) and we avoid the word "variance." |
| 3 | **A value labeled "Simpson's index" actually held `1 − D`** while the function of the same name returned `D`. | 🟠 Med | Table now writes **both** `simpson_D` and `gini_simpson`; plots/labels say "Gini-Simpson (1−D)"; docstrings clarified. |
| 4 | **Misleading error message** in `sum_of_maxima` ("no species common…") — it only triggers when **both** maps are empty. | 🟡 Low | Corrected message + comment. |
| 5 | **Unseeded RNG** made plot 1's jitter change every run. | 🟡 Low | `np.random.seed(0)`. |
| 6 | **Dead entry point** — `main.py` just printed "Metagenomics!"; `io_python.main`/`helperFunctions.main` were empty stubs. | 🟡 Low | `main.py` now calls `analyze.main()`. |
| 7 | **`.gitignore` listed `Data/`** while the data was committed — new samples would be silently ignored. | 🟡 Low | Removed `Data/` from ignore with an explanatory note. |

**Not a bug, but clarified:** `jaccard_distance` is the **abundance-weighted (Ružička)**
Jaccard, not classic presence/absence Jaccard. Documented loudly so we don't get caught on it.

---

## B. Improvements we added (beyond bug-fixing)

1. **PERMANOVA — a real statistical test** (plot 8). One-way permutational ANOVA on the
   distance matrix, treating samples (not pairs) as units. **Season R² = 0.59, p = 0.001;
   Site R² = 0.08, p = 0.675.** This turns "season ≫ site" from a descriptive average into a
   significant, quantified result and defeats the pseudoreplication objection head-on.
2. **Leave-one-out + leave-one-GROUP-out robustness.** Season beats site in **44/44**
   single-sample drops, **4/4** whole-season drops, and **5/5** whole-site drops. Directly
   answers the rubric's "what if one sample were thrown out?" — and the stronger group version.
3. **Genuinely independent metric check.** Added classic **presence/absence Jaccard** (shares
   no numerator with Bray-Curtis); season still beats site **3.6×**, so metric-independence is
   real, not a tautology.
4. **PCoA ordination plot** (plot 5) — the clearest single view of season-driven clustering,
   computed by hand (classical MDS) so we can explain every step.
5. **Bray-Curtis vs Jaccard agreement** (plot 6) — and we're honest that abundance-weighted
   Jaccard is mathematically linked to BC (the binary Jaccard in #3 is the independent one).
6. **Richness-vs-depth plot** (plot 7) — proactively surfaces the depth confound.
7. **`effect_summary.csv`** — an auditable record so every claim is traceable to a number.
8. **Better figures** — colorblind-safe palette, site means, season-blocked heatmap with
   labels, fixed legends/titles.
9. **Full documentation** — `DOCUMENTATION.md` + four reports + a presentation script.

---

## C. Honest limitations (say these before you're asked)

1. **Exact-match OTUs.** Two reads are the "same species" only if byte-for-byte identical, so
   one sequencing error = one fake new OTU. This **inflates richness** and is *why the
   distilled-water control looks rich*. Real pipelines cluster at ~97% or denoise into ASVs
   (DADA2). Our beta-diversity conclusions are more robust to this than richness is.
2. **No rarefaction / depth normalization.** Read depth ranges ~520–5,600; richness tracks it
   (`r ≈ 0.56`). We'd subsample to equal depth before comparing richness.
3. **PERMANOVA is one-way per factor.** We *did* run PERMANOVA (season R²=0.59 p=0.001; site
   R²=0.08 p=0.675), which fixes the pseudoreplication problem of comparing raw pairwise means.
   But it's one-way for each factor separately; a **two-way/nested** model would partition site
   and season jointly. (That can only *raise* site's share, and site is already non-significant.)
4. **Unbalanced, sparse design.** Fall has 19 samples, Spring 7; some site×season cells have
   a single sample, so per-site rankings are preliminary, not definitive.
5. **Season is a proxy.** "Season" bundles temperature, rainfall/runoff, and nutrient load —
   we show season *correlates* with composition, not which mechanism causes it.
6. **One year, one region, one gene.** 2019 only, Pittsburgh only, 16S only (bacteria, not
   eukaryotes/viruses); conclusions don't automatically generalize.

---

## D. If we had another week (future work)
*(PERMANOVA, presence/absence Jaccard, and leave-one-group-out are already done — section B.)*
- **Rarefy** all samples to equal sequencing depth and recompute richness and Bray-Curtis, to
  rule out the depth–season confound (our biggest remaining open weakness).
- Run a **two-way / nested PERMANOVA** to partition site and season jointly.
- Try **97%-similarity OTU clustering** (or ASV denoising) and check whether the season signal
  strengthens once sequencing noise is removed.
- Use the high-richness control to build a **`decontam`-style contaminant blocklist** and
  subtract it from every sample.
- Track a **single site across all four seasons** to show the within-site annual cycle directly.
