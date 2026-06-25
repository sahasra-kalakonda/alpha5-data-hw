# Three Rivers Metagenomics Analysis

CMU Pre-College Computational Biology — Module 1 project.

We take bacterial 16S sequencing reads from water samples around Pittsburgh's three
rivers and turn them into a scientific argument about **how bacterial communities
differ across the Three Rivers**. The short answer our analysis supports:

> **Season explains far more of the variation in community composition than site does**,
> the distilled-water **control** is correctly flagged by beta diversity (our sanity check
> passes), and the result is **robust** (it survives dropping any single sample and holds
> under a second distance metric).

## How to run

You need Python 3.10+ with `numpy` and `matplotlib`.

**Windows (PowerShell):**
```powershell
python -m venv venv
venv\Scripts\Activate.ps1
pip install numpy matplotlib
python analyze.py
```

**macOS / Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy matplotlib
python3 analyze.py
```

`python analyze.py` (or `python main.py`) reads everything in `Data/2019_Samples/`,
writes tables + figures into `Output/`, and prints data-derived conclusions.

## What you get in `Output/`

| File | What it is |
|---|---|
| `alpha_diversity.csv` | per-sample richness, Simpson D, and Gini-Simpson (1−D) |
| `beta_diversity_bray_curtis.csv` | 49×49 Bray-Curtis distance matrix |
| `beta_diversity_jaccard.csv` | 49×49 abundance-weighted Jaccard (Ružička) matrix |
| `effect_summary.csv` | the key numbers behind every printed conclusion (auditable) |
| `plot1_alpha_by_site.png` | richness + evenness by site — alpha can't separate sites |
| `plot2_control_check.png` | beta diversity flags the control |
| `plot3_site_vs_season.png` | **centerpiece**: season ≫ site |
| `plot4_heatmap_by_season.png` | pairwise BC heatmap, season-blocked |
| `plot5_pcoa.png` | **money shot**: PCoA, clusters by season not site |
| `plot6_bc_vs_jaccard.png` | the two distance metrics agree (r ≈ 0.99) |
| `plot7_richness_vs_depth.png` | honesty: richness is confounded by read depth |
| `plot8_permanova.png` | **the formal test**: PERMANOVA R² — season 0.59 (p=0.001) vs site 0.08 (p=0.68) |

## Repo layout

```
io_python.py           read .txt files -> {sample: {sequence: count}}
helperFunctions.py     the math: richness, Simpson, Bray-Curtis, Jaccard
richnessMap.py         richness across all samples
simpsonsMap.py         Gini-Simpson (1-D) across all samples
betaDiversityMatrix.py pairwise distance matrix for a chosen metric
analyze.py             the pipeline: load -> measure -> plot -> conclude
main.py                thin entry point (calls analyze.main)
```

## Documentation

- `DOCUMENTATION.md` — full technical reference (data format, every function, the math).
- `REPORT_1_what_the_program_does.md` — plain-English pipeline walkthrough.
- `REPORT_2_math_and_algorithms.md` — the diversity math, derived and justified.
- `REPORT_3_visualizations.md` — every figure: how to read it and the claim it backs.
- `REPORT_4_bugs_fixed_and_improvements.md` — what was wrong, what we fixed, what's next.
- `PRESENTATION_SCRIPT.md` — 7-minute 4-person script + Q&A defense bank.
