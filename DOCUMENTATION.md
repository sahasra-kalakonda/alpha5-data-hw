# Technical Documentation — Three Rivers Metagenomics

This is the full reference for the codebase: the data, every module and function, the
math, and the outputs. If you can explain this document, you can defend the project.

---

## 1. The data

`Data/2019_Samples/` holds ~49 plain-text files, one per water sample, collected in 2019.

- **File name** encodes the metadata: `Season_Site_Replicate.txt`, e.g. `Fall_Sharpsburg_4.txt`
  → season `Fall`, site `Sharpsburg`, replicate `4`.
- **File contents**: one line = one **DNA read** (a ~240-base fragment of the bacterial
  **16S rRNA gene**, the standard "barcode" gene for identifying bacteria).
- **Sites**: `Allegheny`, `Braddock`, `Monongahela`, `Neville`, `Sharpsburg`, and
  `Control` (distilled water — a negative control).
- **Seasons**: Spring, Summer, Fall, Winter.

The design is **unbalanced**: Fall has 19 samples, Winter 15, Summer 8, Spring 7, and not
every site was sampled in every season (5 of those samples are controls: 3 Fall, 1 Spring,
1 Winter). There are 44 river samples and 5 controls.

### What is an OTU here?
An **OTU** (Operational Taxonomic Unit) is our stand-in for "a species." In this code,
**two reads belong to the same OTU only if their sequences are byte-for-byte identical.**
This is the simplest possible definition (exact-match), and it has a consequence worth
owning: a single sequencing error produces a brand-new OTU. (Real pipelines cluster reads
at ~97% similarity or denoise into ASVs to avoid this — see `REPORT_4`.)

### The frequency map
Each sample is reduced to a **frequency map**: a dictionary `{sequence: count}` giving how
many reads of each distinct sequence appeared. Everything downstream operates on these maps.

```
Fall_Sharpsburg_4.txt          ->   {"TGGG...AA": 86, "TAAG...GT": 12, ...}
(one DNA read per line)              (distinct sequence -> # of reads)
```

---

## 2. Module-by-module reference

### `io_python.py` — reading samples
- `read_freq_map_from_file(filename) -> dict[str,int]`
  Reads one file, strips each line, and counts occurrences → a frequency map.
- `read_samples_from_directory(directory) -> dict[str, dict[str,int]]`
  Runs the above over every `.txt` in the directory, keyed by sample name (filename minus
  `.txt`). Returns the **"map of frequency maps"** the whole analysis consumes.

### `helperFunctions.py` — the math (see §3 for derivations)
- `richness(sample) -> int` — number of OTUs with count > 0. Raises on negative counts.
- `simpsons_index(sample) -> float` — Simpson **dominance** `D = Σ(n/N)²`.
- `frequency_map(patterns) -> dict[str,int]` — build a frequency map from a list.
- `sample_total(freq_map) -> int` — total reads `N` in a sample.
- `sum_of_minima(map1, map2) -> int` — `Σ min(aᵢ, bᵢ)` over shared OTUs (`C`).
- `sum_of_maxima(map1, map2) -> int` — `Σ max(aᵢ, bᵢ)` over the union of OTUs.
- `bray_curtis_distance(map1, map2) -> float` — `1 − 2C/(N₁+N₂)`.
- `jaccard_distance(map1, map2) -> float` — `1 − C/Σmax` (**abundance-weighted Ružička**).
- `jaccard_binary_distance(map1, map2) -> float` — classic **presence/absence** Jaccard
  `1 − |A∩B|/|A∪B|` over the OTU *sets* (our independent metric check).
- `max_2`, `min_2`, `average` — tiny arithmetic helpers (written out for the course).

### `richnessMap.py`
- `richness_map(all_maps) -> dict[str,int]` — applies `richness` to every sample.

### `simpsonsMap.py`
- `simpsons_map(all_maps) -> dict[str,float]` — returns **Gini-Simpson `1 − D`** per sample
  (an *evenness/diversity* measure; higher = more even). Note this is *not* `D` itself.

### `betaDiversityMatrix.py`
- `beta_diversity_matrix(all_maps, metric) -> (names, matrix)`
  Sorts sample names (reproducible order), selects the metric by name (`"Bray-Curtis"`,
  `"Jaccard"`, or `"Jaccard-binary"`), and fills an `n×n` distance matrix. Because the metrics
  are **symmetric** (`d(A,B)=d(B,A)`) and `d(A,A)=0`, it computes only the **upper triangle**
  (`j ≥ i`) and mirrors each value — roughly halving the work. Unknown names raise `ValueError`.

### `analyze.py` — the pipeline (see §4)
Loads data, computes alpha + beta, writes 4 CSVs, renders 8 figures, runs PERMANOVA plus
leave-one-out / leave-one-group-out robustness checks, and prints conclusions whose numbers
are all computed (never hardcoded).
Key internal helpers: `parse_sample_name`, `is_control`, `category_means` (sorts river pairs
into the four site/season categories), and `pcoa` (classical MDS by hand).

### `main.py`
Thin entry point: `python main.py` simply calls `analyze.main()`.

---

## 3. The math, precisely

Notation: a sample is counts `a = (a₁,…)` over OTUs; `N = Σ aᵢ` is its total reads.

### Alpha diversity (one sample at a time)
- **Richness** `S` = number of OTUs with `aᵢ > 0`. "How many kinds?" Ignores abundances.
- **Simpson dominance** `D = Σ (aᵢ/N)²` = probability two random reads are the **same** OTU.
  Ranges `1/S … 1`; **higher = less diverse** (one taxon dominates).
- **Gini-Simpson** `1 − D` = probability two random reads are **different** OTUs.
  **Higher = more even.** This is what we report and plot (it reads as "more diverse = higher").

> Why both forms exist in our output: the function `simpsons_index` returns `D`, but the
> table column and plots use `1 − D`. We deliberately write **both** `simpson_D` and
> `gini_simpson` to `alpha_diversity.csv` so there is zero ambiguity about which is which.

### Beta diversity (between two samples)
Let `C = Σ min(aᵢ, bᵢ)` (shared abundance) and `Σmax = Σ max(aᵢ, bᵢ)` over the OTU union.

- **Bray-Curtis** `BC = 1 − 2C / (N₁ + N₂)`.
  `0` = identical communities, `1` = no shared reads. Abundance-sensitive.
- **Jaccard (as implemented)** `J = 1 − C / Σmax`. This is the **abundance-weighted Jaccard
  = Ružička distance**. It reduces to the classic presence/absence Jaccard
  (`1 − |A∩B|/|A∪B|`) only if every count is 1.

**Useful identity** (why our two metrics agree so well): writing `s = Σmin`, `t = Σmax`,
and noting `N₁+N₂ = Σmin + Σmax` for the union form, one gets
`J = 2·BC / (1 + BC)` — a fixed monotonic curve. So `J ≥ BC` always and `r ≈ 0.99` is
**expected**, not a coincidence (see `plot6`).

### PCoA (Principal Coordinates Analysis) — `analyze.pcoa`
Turns the `n×n` distance matrix into 2-D points you can scatter-plot. Classical MDS:
1. `D2 = D²` (element-wise).
2. Double-centre: `B = −½ · J · D2 · J`, where `J = I − (1/n)·11ᵀ`.
3. Eigen-decompose the symmetric `B` (`numpy.linalg.eigh`).
4. The top-2 eigenvectors scaled by `√eigenvalue` are the 2-D coordinates; each
   eigenvalue ÷ (sum of positive eigenvalues) is the fraction of variation that axis explains.

Because Bray-Curtis isn't perfectly Euclidean, a few eigenvalues can be slightly negative;
we report variance fractions against the sum of the **positive** ones (a standard convention).

---

## 4. The analysis logic (`analyze.py`)

1. **Load** all samples → map of frequency maps.
2. **Alpha**: `richness_map` + `simpsons_map` → `alpha_diversity.csv`.
3. **Beta**: build Bray-Curtis and Jaccard matrices → two CSVs.
4. **Categorise river pairs** (`category_means`): every river-vs-river pair is one of
   - same site, same season → **replicates** (our noise floor),
   - different site, same season → the **site** effect,
   - same site, different season → the **season** effect,
   - different site, different season.
5. **Quantify the effect above baseline**:
   `site_effect = mean(diff-site) − mean(replicate)`,
   `season_effect = mean(diff-season) − mean(replicate)`,
   `ratio = season_effect / site_effect`.
6. **PERMANOVA** (`permanova`): one-way permutation test (999 perms) for season and for site,
   reporting R², pseudo-F, and p — the statistically proper "variance explained."
7. **Robustness**: leave-one-sample-out (44/44), leave-one-season-out (4/4), leave-one-site-out
   (5/5), and a presence/absence-Jaccard re-check (season still wins 3.6×).
8. **Render 8 plots**, **write `effect_summary.csv`**, and **print** conclusions built from
   the computed numbers.

---

## 5. Outputs

CSVs and PNGs are listed in `README.md`. The crucial one for defense is
`effect_summary.csv`: every number in the printed conclusions and in the slides comes from
there, so any claim can be traced back to a computed value.

---

## 6. Known limitations (own these — see `REPORT_4` for detail)
- **Exact-match OTUs** inflate richness with sequencing noise (this is *why* the control
  looks "rich").
- **No rarefaction**: richness tracks sequencing depth (`r ≈ 0.56`), so we lean on evenness
  and composition, not raw richness. (This is our biggest *remaining* open weakness — season
  and depth are partly collinear.)
- **PERMANOVA is one-way per factor**: we *do* run PERMANOVA (season R²=0.59 p=0.001; site
  R²=0.08 p=0.675), which fixes the pseudoreplication of raw pairwise means; a two-way/nested
  model would partition site and season jointly.
- **Unbalanced design** and small per-cell `n` (some site×season cells have one sample).
- **Season is a proxy** for temperature/runoff/nutrients — and possibly sequencing batch.
