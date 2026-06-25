# Report 3 — The Visualizations

Eight figures, each built to make **one point**. For every plot: what it shows, how to read
it, and the claim it supports. (Rubric: "every figure should be there to make a point.")

---

## Plot 1 — `plot1_alpha_by_site.png` — Alpha diversity by site
**Shows:** two panels — richness (left) and Gini-Simpson evenness (right) — one dot per
sample, grouped by site, control in gray. The colored dash is each site's mean.

**How to read:** look for vertical separation between sites. There isn't much: the site
clouds overlap heavily and the within-site spread is large.

**Claim it backs:** *Alpha diversity alone cannot tell the sites apart, and richness and
evenness disagree* (across samples `r ≈ −0.01`). Note the control (gray) is **not** unusual
here — setting up plot 2.

---

## Plot 2 — `plot2_control_check.png` — The control sanity check
**Shows:** boxplots of Bray-Curtis distance for river-vs-river pairs vs control-vs-river
pairs (white dot = mean).

**How to read:** the control-vs-river box sits higher (mean **0.92** vs **0.78**).

**Claim it backs:** *Beta diversity flags the distilled-water control even though alpha
diversity didn't.* Our sanity check passes, and the lesson is "cross-check alpha with beta."

---

## Plot 3 — `plot3_site_vs_season.png` — **Centerpiece: site vs season**
**Shows:** mean Bray-Curtis for the four kinds of river pair, with pair counts (`n=`) and the
gray replicate-baseline line at 0.42.

**How to read:** the "different site, same season" bar barely clears the baseline (0.54); the
"same site, different season" bar towers (0.88), level with "everything different."

**Claim it backs:** *Season moves communities ~4× more than site, measured above the
replicate baseline.* This is the single most important slide.

---

## Plot 4 — `plot4_heatmap_by_season.png` — Pairwise distance heatmap
**Shows:** the full 49×49 Bray-Curtis matrix, samples ordered by season then site; white
lines + colored labels separate the four season blocks; dark = similar, yellow = different.

**How to read:** dark squares sit **on the diagonal within each season block**; off the
season blocks it's mostly yellow.

**Claim it backs:** *Samples cluster by season* — a visual, all-pairs confirmation of plot 3
that doesn't depend on our four-bucket averaging.

---

## Plot 5 — `plot5_pcoa.png` — **PCoA ordination (the money shot)**
**Shows:** every sample as a point in a 2-D map of community composition (axes capture 36% +
13% of variation). **Same coordinates, two colorings:** by season (left), by site (right);
controls are black ✕.

**How to read:** left panel — four tidy clusters (Spring, Summer, Fall, Winter) with controls
flung off to the side. Right panel — identical points, but coloring by site is a random mix:
no site clusters.

**Claim it backs:** *Composition organizes by season, not site* — the whole thesis in one
picture, plus the control separation, all at once.

---

## Plot 6 — `plot6_bc_vs_jaccard.png` — Do the two metrics agree?
**Shows:** every sample pair as a point: Bray-Curtis (x) vs abundance-weighted Jaccard (y),
with the `y = x` line; `r = 0.995`.

**How to read:** points form a tight curve sitting just above `y = x`.

**Claim it backs:** *Our conclusion is not a metric artifact* — the two distances rank pairs
almost identically. **Honesty note:** they agree because they're mathematically linked
(`J = 2·BC/(1+BC)`), so this is a *consistency check*, not independent evidence. A truly
independent contrast would use presence/absence Jaccard.

---

## Plot 7 — `plot7_richness_vs_depth.png` — Honesty about limits
**Shows:** richness (y) vs sequencing depth = total reads (x) for river samples, colored by
site, with a trend line; `r = 0.56`.

**How to read:** clear upward slope — deeper-sequenced samples look "richer."

**Claim it backs:** *Richness is partly an artifact of sequencing effort*, which is exactly
why we lean on evenness and composition. Putting this in front of the room pre-empts the
"isn't your richness just read depth?" question by answering it first.

---

## Plot 8 — `plot8_permanova.png` — The formal test (PERMANOVA)
**Shows:** the fraction of compositional variation (R²) explained by season vs by site, from a
one-way PERMANOVA on the river samples, with permutation p-values on each bar.

**How to read:** season bar towers at **R² = 0.59, p = 0.001**; site bar is tiny at
**R² = 0.08, p = 0.675**.

**Claim it backs:** *the rigorous version of "season ≫ site."* Season explains ~59% of the
variation and is highly significant; site explains ~8% and is **not statistically distinguishable
from random**. This is the slide that answers "did you do real statistics, or just compare
averages?" — pair it with plot 3 (the intuition) for the one-two punch.

---

## Design choices you can defend
- **Colorblind-safe palette** (Paul Tol "bright"); seasons get fixed, intuitive colors.
- **Seeded RNG** (`np.random.seed(0)`) so the jittered scatter in plot 1 is reproducible.
- **Pair counts shown** on plot 3 so the audience sees the unbalanced `n`'s honestly.
- **Every figure's number is in `effect_summary.csv`**, so any value on a slide is traceable.
