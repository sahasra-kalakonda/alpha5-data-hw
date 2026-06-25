# Report 1 — What the Program Does

A plain-English walkthrough of the pipeline, in the order the code runs it. Goal: anyone in
the group can narrate the program end to end without looking at the source.

## The one-sentence version
The program reads ~49 water-sample read files, turns each into a count of bacterial OTUs,
measures how diverse each sample is and how different every pair of samples is, and uses
those numbers to argue that **season — not site — is what mostly shapes the river's
bacterial communities**, while double-checking itself with the distilled-water control.

## Step 0 — Input
`Data/2019_Samples/` contains one text file per sample, named `Season_Site_Replicate.txt`.
Each line in a file is one DNA read (a 16S "barcode" fragment). Identical reads = same OTU.

## Step 1 — Load (`io_python.read_samples_from_directory`)
Every `.txt` is read into a **frequency map** `{sequence: count}`. The whole dataset becomes
a dictionary of these maps, keyed by sample name. This is the "map of frequency maps" the
project is built around. Output: `Loaded 49 samples.`

## Step 2 — Alpha diversity (diversity *within* a sample)
For each sample we compute two numbers:
- **Richness** — how many distinct OTUs are present (`richnessMap.richness_map`).
- **Gini-Simpson index `1 − D`** — how *even* the community is, i.e. is it a balanced mix or
  dominated by a few OTUs (`simpsonsMap.simpsons_map`).

These are written to `Output/alpha_diversity.csv` (one row per sample, with both `simpson_D`
and `gini_simpson` to avoid ambiguity).

## Step 3 — Beta diversity (difference *between* two samples)
We build two 49×49 **distance matrices** — one **Bray-Curtis**, one **Jaccard**
(`betaDiversityMatrix.beta_diversity_matrix`). Entry `(i, j)` is how different samples `i`
and `j` are: 0 = identical communities, 1 = nothing in common. Saved as two CSVs. Because
distance is symmetric, the code computes half the matrix and mirrors it.

## Step 4 — Turn distances into an argument
This is where measurement becomes science. Looking only at **river** samples, every pair is
sorted into one of four buckets and we average the Bray-Curtis distance in each:

| Pair type | What it isolates | Mean BC |
|---|---|---|
| same site, same season (replicates) | measurement noise floor | **0.42** |
| different site, same season | the **SITE** effect | **0.54** |
| same site, different season | the **SEASON** effect | **0.88** |
| different site, different season | both | 0.88 |

Changing the site nudges communities a little (`+0.11` above the replicate floor); changing
the season moves them a lot (`+0.45`) — about **4× more**. That is the headline.

## Step 5 — Test it properly (PERMANOVA)
Comparing average distances is intuitive but not a real statistical test (the pairs aren't
independent). So the program runs a **PERMANOVA** — a permutation test that treats the 44
samples as the units and asks how much variation each factor explains:

| Factor | R² (variation explained) | p-value |
|---|---|---|
| **Season** | **0.59** | **0.001** |
| Site | 0.08 | 0.675 |

**Season explains ~59% of the variation and is highly significant; site explains ~8% and is
no better than random.** This is the rigorous version of the headline.

## Step 6 — Stress-test it (three ways)
- **Drop each sample** → season beats site **44/44** times.
- **Drop a whole season** (4/4) and **drop a whole site** (5/5) → still holds. No group drives it.
- **Switch to presence/absence Jaccard** (a metric that shares nothing with Bray-Curtis) →
  season still wins by **3.6×**, so it's not an artifact of how we measure distance.

(This directly answers the rubric's "what if one sample were thrown out?")

## Step 7 — The control (sanity check)
Distilled water should contain essentially no river community. By **richness** it doesn't
look special (Fall controls were among the *richest* samples — sequencing noise). But by
**composition** it's an outlier: its mean Bray-Curtis distance to river samples is **0.92**
vs **0.78** among river samples. Beta diversity catches what alpha diversity misses — so our
pipeline passes its own sanity check.

## Step 8 — Report the limits honestly
- **Richness vs read depth** correlate (`r ≈ 0.56`), which is *why* we trust evenness and
  composition over raw richness rather than ranking sites by richness.

## Step 9 — Outputs
Four CSV tables (alpha, two beta matrices, an auditable `effect_summary.csv`) and **eight**
figures land in `Output/`, and the program prints conclusions whose every number was
computed in the steps above — nothing is hardcoded.
