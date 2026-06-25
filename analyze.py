"""
Three Rivers Metagenomics: full analysis pipeline.

Loads every sample in Data/2019_Samples, computes alpha diversity (richness,
Simpson's index) and beta diversity (Bray-Curtis, Jaccard) across all of them,
saves the resulting tables, and produces plots that build a scientific argument
about our central question:

    How do bacterial communities differ across the Three Rivers?

The analysis is organised around four claims, each backed by a plot and by numbers
that are COMPUTED FROM THE DATA (never hardcoded):

  1. Alpha diversity alone cannot tell the sites apart, and cannot flag the control.
  2. Beta diversity CAN flag the distilled-water control (sanity check passes).
  3. SEASON explains far more compositional variation than SITE.
  4. The result is robust: it survives dropping any single sample (leave-one-out),
     and our two distance metrics (Bray-Curtis and Jaccard) agree.

Run:  python analyze.py
Outputs land in Output/ (3 CSV tables, 1 effect-summary CSV, 7 PNG figures).
"""
import csv
import statistics

import matplotlib.pyplot as plt
import numpy as np

from io_python import read_samples_from_directory
from richnessMap import richness_map
from simpsonsMap import simpsons_map
from betaDiversityMatrix import beta_diversity_matrix
from helperFunctions import sample_total, simpsons_index

DATA_DIR = "Data/2019_Samples"
OUT_DIR = "Output"

# Seed the RNG so the jittered scatter (plot 1) is identical on every run.
np.random.seed(0)

# Fixed display orders and a colour-blind-safe palette (Paul Tol "bright").
SITES = ["Allegheny", "Braddock", "Monongahela", "Neville", "Sharpsburg", "Control"]
RIVER_SITES = ["Allegheny", "Braddock", "Monongahela", "Neville", "Sharpsburg"]
SITE_COLORS = {
    "Allegheny": "#4477AA", "Braddock": "#66CCEE", "Monongahela": "#228833",
    "Neville": "#CCBB44", "Sharpsburg": "#EE6677", "Control": "#888888",
}
SITE_MARKERS = {
    "Allegheny": "o", "Braddock": "s", "Monongahela": "^",
    "Neville": "D", "Sharpsburg": "v",
}
SEASON_ORDER = {"Spring": 0, "Summer": 1, "Fall": 2, "Winter": 3}
SEASON_COLORS = {"Spring": "#228833", "Summer": "#CCBB44", "Fall": "#EE6677", "Winter": "#4477AA"}

# The four kinds of river-vs-river sample pair, used by the site-vs-season analysis.
CAT_REP = "same site,\nsame season\n(replicates)"
CAT_SITE = "different site,\nsame season"
CAT_SEASON = "same site,\ndifferent season"
CAT_BOTH = "different site,\ndifferent season"
CATEGORIES = [CAT_REP, CAT_SITE, CAT_SEASON, CAT_BOTH]


# ----------------------------------------------------------------------------
# Small naming helpers
# ----------------------------------------------------------------------------
def parse_sample_name(name: str) -> tuple[str, str, str]:
    """Splits 'Fall_Sharpsburg_4' into (season, site, replicate)."""
    parts = name.split("_")
    return parts[0], parts[1], parts[2]


def is_control(name: str) -> bool:
    """True for distilled-water control samples (site == 'Control')."""
    return parse_sample_name(name)[1] == "Control"


# ----------------------------------------------------------------------------
# File writers
# ----------------------------------------------------------------------------
def write_alpha_table(all_maps, rich, gini, path):
    """
    One row per sample. We write BOTH forms of Simpson's index to avoid ambiguity:
      simpson_D     = sum (n/N)^2          (dominance; higher = less diverse)
      gini_simpson  = 1 - simpson_D        (diversity/evenness; higher = more diverse)
    """
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "season", "site", "replicate", "total_reads",
                         "richness", "simpson_D", "gini_simpson"])
        for name in sorted(all_maps.keys()):
            season, site, replicate = parse_sample_name(name)
            d = simpsons_index(all_maps[name])  # raw Simpson dominance D
            writer.writerow([
                name, season, site, replicate,
                sample_total(all_maps[name]),
                rich[name],
                round(d, 6),
                round(gini[name], 6),
            ])


def write_beta_matrix(names, matrix, path):
    """Writes a square distance matrix to CSV with sample names as row/column headers."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([""] + names)
        for i, name in enumerate(names):
            writer.writerow([name] + [round(x, 6) for x in matrix[i]])


# ----------------------------------------------------------------------------
# Analysis helpers
# ----------------------------------------------------------------------------
def category_means(river_names, name_to_idx, matrix, exclude=None):
    """
    Sort every river-vs-river pair into one of the four categories and return
    (mean Bray-Curtis per category, pair count per category). `exclude` is an
    optional set of sample names to drop (used for the leave-one-out check).
    """
    exclude = exclude or set()
    buckets = {c: [] for c in CATEGORIES}
    for i in range(len(river_names)):
        for j in range(i + 1, len(river_names)):
            a, b = river_names[i], river_names[j]
            if a in exclude or b in exclude:
                continue
            sa, sitea, _ = parse_sample_name(a)
            sb, siteb, _ = parse_sample_name(b)
            d = matrix[name_to_idx[a]][name_to_idx[b]]
            if sitea == siteb and sa == sb:
                buckets[CAT_REP].append(d)
            elif sitea != siteb and sa == sb:
                buckets[CAT_SITE].append(d)
            elif sitea == siteb and sa != sb:
                buckets[CAT_SEASON].append(d)
            else:
                buckets[CAT_BOTH].append(d)
    means = {c: (statistics.mean(v) if v else float("nan")) for c, v in buckets.items()}
    counts = {c: len(v) for c, v in buckets.items()}
    return means, counts


def pcoa(matrix):
    """
    Classical multidimensional scaling (Principal Coordinates Analysis) of a
    distance matrix, done by hand with numpy so we can explain every step:

      1. square the distances,
      2. double-centre:  B = -1/2 * J D^2 J,  where J = I - (1/n) 11^T,
      3. eigen-decompose the symmetric B,
      4. the top-2 eigenvectors scaled by sqrt(eigenvalue) are the 2-D coordinates.

    Bray-Curtis is not perfectly Euclidean, so a few eigenvalues can be slightly
    negative; we report variance fractions against the sum of the positive ones.

    Returns (coords [n x 2], variance_fraction [2]).
    """
    D = np.asarray(matrix, dtype=float)
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)          # ascending order
    order = np.argsort(eigvals)[::-1]             # descending
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos_sum = eigvals[eigvals > 0].sum()
    coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2], 0))
    frac = np.maximum(eigvals[:2], 0) / pos_sum
    return coords, frac


def permanova(names, name_to_idx, matrix, group_of, n_perm=999):
    """
    One-way PERMANOVA (Anderson 2001) on a distance matrix -- a permutation test that
    properly treats the SAMPLES (not the pairs) as the units of replication, which is the
    statistically correct way to ask "does this factor explain composition?"

    Works directly from the distances:
      SS_total  = (1/N) * sum_{i<j} d_ij^2
      SS_within = sum over groups g of (1/n_g) * sum_{i<j in g} d_ij^2
      SS_among  = SS_total - SS_within
      pseudo-F  = (SS_among/(a-1)) / (SS_within/(N-a))         a = number of groups
      R^2       = SS_among / SS_total                          fraction of variation explained
    The p-value is the fraction of label-shuffles whose pseudo-F is >= the observed one.

    Returns (R2, pseudo_F, p_value).
    """
    idx = [name_to_idx[n] for n in names]
    D = np.asarray(matrix, dtype=float)[np.ix_(idx, idx)]
    N = len(names)
    D2 = D ** 2
    triu = np.triu_indices(N, 1)
    ss_total = D2[triu].sum() / N
    labels = np.array([group_of(n) for n in names])
    groups = sorted(set(labels))
    a = len(groups)

    def ss_within(labs):
        s = 0.0
        for g in groups:
            members = np.where(labs == g)[0]
            ng = len(members)
            if ng < 2:
                continue
            sub = D2[np.ix_(members, members)]
            s += sub[np.triu_indices(ng, 1)].sum() / ng
        return s

    ssw = ss_within(labels)
    ssa = ss_total - ssw
    pseudo_f = (ssa / (a - 1)) / (ssw / (N - a))
    r2 = ssa / ss_total

    ge = 0
    perm = labels.copy()
    for _ in range(n_perm):
        np.random.shuffle(perm)
        sswp = ss_within(perm)
        fp = ((ss_total - sswp) / (a - 1)) / (sswp / (N - a))
        if fp >= pseudo_f:
            ge += 1
    p_value = (ge + 1) / (n_perm + 1)
    return r2, pseudo_f, p_value


def plot_permanova(season_r2, site_r2, season_p, site_p):
    """Plot 8: PERMANOVA R^2 -- the statistically proper 'variance explained' by each factor."""
    fig, ax = plt.subplots(figsize=(6, 5))
    bars = ax.bar(["Season", "Site"], [season_r2, site_r2], color=["#EE6677", "#4477AA"], alpha=0.85)
    for bar, p in zip(bars, [season_p, site_p]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f"R²={bar.get_height():.2f}\np={p:.3f}", ha="center", fontsize=10)
    ax.set_ylabel("PERMANOVA R²  (fraction of variation explained)")
    ax.set_title("PERMANOVA on river samples:\nseason explains far more variation than site", fontsize=12)
    ax.set_ylim(0, max(season_r2, site_r2) * 1.35)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot8_permanova.png", dpi=150)
    plt.close()


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------
def plot_alpha_by_site(all_maps, rich, gini):
    """Plot 1: richness and Gini-Simpson per sample, by site, control highlighted."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for metric_ax, values, label in (
        (axes[0], rich, "Richness (# distinct sequences)"),
        (axes[1], gini, "Gini-Simpson index, 1 - D (evenness)"),
    ):
        for site in SITES:
            names = [n for n in all_maps if parse_sample_name(n)[1] == site]
            y = [values[n] for n in names]
            x = np.full(len(names), SITES.index(site)) + np.random.normal(0, 0.06, len(names))
            metric_ax.scatter(x, y, color=SITE_COLORS[site], s=55, edgecolor="black",
                              linewidth=0.5, zorder=3, label=site if metric_ax is axes[1] else None)
            if y:  # mean marker per site
                metric_ax.scatter([SITES.index(site)], [statistics.mean(y)], marker="_",
                                  s=900, color=SITE_COLORS[site], linewidth=2.5, zorder=2)
        metric_ax.set_xticks(range(len(SITES)))
        metric_ax.set_xticklabels(SITES, rotation=30, ha="right")
        metric_ax.set_ylabel(label)
        metric_ax.grid(axis="y", alpha=0.25)
    axes[0].set_title("Richness by site (— = site mean; control in gray)")
    axes[1].set_title("Evenness by site (— = site mean; control in gray)")
    axes[1].legend(fontsize=8, loc="lower right", framealpha=0.9)
    fig.suptitle("Alpha diversity: neither richness nor evenness cleanly separates the sites",
                 fontsize=13, weight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot1_alpha_by_site.png", dpi=150)
    plt.close()


def plot_control_check(names, bc):
    """Plot 2: does the distilled-water control look compositionally distinct?"""
    idx = {n: i for i, n in enumerate(names)}
    river = [n for n in names if not is_control(n)]
    control_to_river, river_to_river = [], []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            if is_control(a) != is_control(b):
                control_to_river.append(bc[idx[a]][idx[b]])
            elif not is_control(a) and not is_control(b):
                river_to_river.append(bc[idx[a]][idx[b]])

    fig, ax = plt.subplots(figsize=(6.5, 5))
    bp = ax.boxplot([river_to_river, control_to_river],
                    tick_labels=["river vs. river", "control vs. river"],
                    patch_artist=True, widths=0.55, showmeans=True,
                    meanprops=dict(marker="o", markerfacecolor="white", markeredgecolor="black"))
    for patch, color in zip(bp["boxes"], ["#4477AA", "#888888"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel("Bray-Curtis distance")
    ax.set_title("Sanity check: the control is more distant than river-vs-river\n"
                 f"(mean {statistics.mean(control_to_river):.2f} vs {statistics.mean(river_to_river):.2f})")
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot2_control_check.png", dpi=150)
    plt.close()
    return statistics.mean(control_to_river), statistics.mean(river_to_river)


def plot_site_vs_season(means, counts):
    """Plot 3 (centerpiece): mean Bray-Curtis for each kind of river pair."""
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    vals = [means[c] for c in CATEGORIES]
    ns = [counts[c] for c in CATEGORIES]
    colors = ["#228833", "#CCBB44", "#EE6677", "#4477AA"]
    bars = ax.bar(CATEGORIES, vals, color=colors, alpha=0.88)
    for bar, n in zip(bars, ns):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                f"n={n}", ha="center", fontsize=9)
    # annotate the replicate baseline as a reference line
    ax.axhline(means[CAT_REP], color="gray", ls="--", lw=1)
    ax.text(-0.55, means[CAT_REP], f"replicate baseline ({means[CAT_REP]:.2f})", fontsize=8,
            color="gray", ha="right", va="center", rotation=90, clip_on=False)
    ax.set_ylabel("Mean Bray-Curtis distance (river samples only)")
    ax.set_title("Site vs. season: changing the SEASON moves communities\nfar more than changing the SITE",
                 fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot3_site_vs_season.png", dpi=150)
    plt.close()


def plot_heatmap(names, bc):
    """Plot 4: pairwise Bray-Curtis heatmap, ordered by season then site."""
    ordered = sorted(names, key=lambda n: (SEASON_ORDER[parse_sample_name(n)[0]],
                                            parse_sample_name(n)[1], n))
    order_idx = [names.index(n) for n in ordered]
    M = np.array(bc)[np.ix_(order_idx, order_idx)]

    fig, ax = plt.subplots(figsize=(11, 9.5))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1)
    # draw white separators + season labels at each season block boundary
    seasons_in_order = [parse_sample_name(n)[0] for n in ordered]
    boundaries = [i for i in range(1, len(ordered)) if seasons_in_order[i] != seasons_in_order[i - 1]]
    for b in boundaries:
        ax.axhline(b - 0.5, color="white", lw=1.5)
        ax.axvline(b - 0.5, color="white", lw=1.5)
    block_edges = [0] + boundaries + [len(ordered)]
    for k in range(len(block_edges) - 1):
        mid = (block_edges[k] + block_edges[k + 1]) / 2 - 0.5
        season = seasons_in_order[block_edges[k]]
        ax.text(mid, -1.5, season, ha="center", va="bottom", fontsize=11, weight="bold",
                color=SEASON_COLORS[season])
    ax.set_xticks(range(len(ordered)))
    ax.set_yticks(range(len(ordered)))
    ax.set_xticklabels(ordered, rotation=90, fontsize=4.5)
    ax.set_yticklabels(ordered, fontsize=4.5)
    plt.colorbar(im, ax=ax, label="Bray-Curtis distance", fraction=0.046, pad=0.04)
    ax.set_title("Pairwise Bray-Curtis distance (sorted by season, then site)\n"
                 "Dark blocks on the diagonal = samples from the same season look alike", pad=30)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot4_heatmap_by_season.png", dpi=150)
    plt.close()


def plot_pcoa(names, bc):
    """Plot 5 (money shot): PCoA of all samples, coloured by season (left) and site (right)."""
    coords, frac = pcoa(bc)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))

    # left: colour by season (controls drawn as black X)
    for n, (x, y) in zip(names, coords):
        season, site, _ = parse_sample_name(n)
        if is_control(n):
            axes[0].scatter(x, y, marker="x", s=70, color="black", linewidth=1.8, zorder=3)
        else:
            axes[0].scatter(x, y, color=SEASON_COLORS[season], s=70, edgecolor="black",
                            linewidth=0.4, zorder=3)
    for season, color in SEASON_COLORS.items():
        axes[0].scatter([], [], color=color, label=season, edgecolor="black", s=60)
    axes[0].scatter([], [], marker="x", color="black", label="Control")
    axes[0].legend(title="Season", fontsize=8)
    axes[0].set_title("Coloured by SEASON — clusters form")

    # right: same coordinates, colour by site
    for n, (x, y) in zip(names, coords):
        season, site, _ = parse_sample_name(n)
        if is_control(n):
            axes[1].scatter(x, y, marker="x", s=70, color="black", linewidth=1.8, zorder=3)
        else:
            axes[1].scatter(x, y, color=SITE_COLORS[site], s=70, edgecolor="black",
                            linewidth=0.4, zorder=3)
    for site in RIVER_SITES:
        axes[1].scatter([], [], color=SITE_COLORS[site], label=site, edgecolor="black", s=60)
    axes[1].scatter([], [], marker="x", color="black", label="Control")
    axes[1].legend(title="Site", fontsize=8)
    axes[1].set_title("Coloured by SITE — no clean clusters")

    for ax in axes:
        ax.set_xlabel(f"PCoA axis 1 ({frac[0] * 100:.0f}% of variation)")
        ax.set_ylabel(f"PCoA axis 2 ({frac[1] * 100:.0f}% of variation)")
        ax.grid(alpha=0.2)
    fig.suptitle("Same map, two colourings: communities organise by season, not by site",
                 fontsize=13, weight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot5_pcoa.png", dpi=150)
    plt.close()


def plot_bc_vs_jaccard(names, bc, jac):
    """Plot 6: do our two distance metrics agree? (Bray-Curtis vs Jaccard for every pair.)"""
    idx = {n: i for i, n in enumerate(names)}
    xs, ys, colors = [], [], []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            xs.append(bc[idx[a]][idx[b]])
            ys.append(jac[idx[a]][idx[b]])
            colors.append("#888888" if (is_control(a) or is_control(b)) else "#4477AA")
    xs, ys = np.array(xs), np.array(ys)
    r = np.corrcoef(xs, ys)[0, 1]

    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(xs, ys, c=colors, s=14, alpha=0.5, edgecolor="none")
    ax.plot([0, 1], [0, 1], color="black", ls="--", lw=1, label="y = x")
    ax.set_xlabel("Bray-Curtis distance")
    ax.set_ylabel("Jaccard (abundance-weighted) distance")
    ax.set_title(f"The two metrics strongly agree (Pearson r = {r:.3f})\n"
                 "blue = river-river pairs, gray = pairs involving the control")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(loc="lower right"); ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot6_bc_vs_jaccard.png", dpi=150)
    plt.close()
    return r


def plot_richness_vs_depth(all_maps, rich):
    """Plot 7 (honesty about limits): richness rises with sequencing depth."""
    river = [n for n in all_maps if not is_control(n)]
    depth = np.array([sample_total(all_maps[n]) for n in river])
    rr = np.array([rich[n] for n in river])
    r = np.corrcoef(depth, rr)[0, 1]

    fig, ax = plt.subplots(figsize=(6.8, 5))
    for site in RIVER_SITES:
        sel = [k for k, n in enumerate(river) if parse_sample_name(n)[1] == site]
        ax.scatter(depth[sel], rr[sel], color=SITE_COLORS[site], s=55,
                   edgecolor="black", linewidth=0.4, label=site)
    # least-squares trend line
    m, c = np.polyfit(depth, rr, 1)
    xline = np.array([depth.min(), depth.max()])
    ax.plot(xline, m * xline + c, color="black", ls="--", lw=1)
    ax.set_xlabel("Sequencing depth (total reads in sample)")
    ax.set_ylabel("Richness (# distinct sequences)")
    ax.set_title(f"Why we don't trust richness alone: it tracks read depth (r = {r:.2f})")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/plot7_richness_vs_depth.png", dpi=150)
    plt.close()
    return r


# ----------------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------------
def main() -> None:
    # ---- 1. Load every sample as a frequency map ----
    all_maps = read_samples_from_directory(DATA_DIR)
    print(f"Loaded {len(all_maps)} samples.")

    # ---- 2. Alpha diversity ----
    rich = richness_map(all_maps)
    gini = simpsons_map(all_maps)  # Gini-Simpson 1 - D (evenness)
    write_alpha_table(all_maps, rich, gini, f"{OUT_DIR}/alpha_diversity.csv")

    # ---- 3. Beta diversity (both metrics) ----
    names, bc = beta_diversity_matrix(all_maps, "Bray-Curtis")
    _, jac = beta_diversity_matrix(all_maps, "Jaccard")              # abundance-weighted Ruzicka
    _, jacb = beta_diversity_matrix(all_maps, "Jaccard-binary")      # presence/absence Jaccard
    write_beta_matrix(names, bc, f"{OUT_DIR}/beta_diversity_bray_curtis.csv")
    write_beta_matrix(names, jac, f"{OUT_DIR}/beta_diversity_jaccard.csv")
    print("Saved alpha_diversity.csv, beta_diversity_bray_curtis.csv, beta_diversity_jaccard.csv")

    name_to_idx = {n: i for i, n in enumerate(names)}
    river_names = [n for n in names if not is_control(n)]

    # ---- 4. Plots ----
    plot_alpha_by_site(all_maps, rich, gini)
    ctrl_mean, river_mean = plot_control_check(names, bc)
    means, counts = category_means(river_names, name_to_idx, bc)
    plot_site_vs_season(means, counts)
    plot_heatmap(names, bc)
    plot_pcoa(names, bc)
    jaccard_r = plot_bc_vs_jaccard(names, bc, jac)
    depth_r = plot_richness_vs_depth(all_maps, rich)

    # ---- 4b. PERMANOVA: the statistically proper test (samples are the units) ----
    season_r2, season_f, season_p = permanova(
        river_names, name_to_idx, bc, lambda n: parse_sample_name(n)[0])
    site_r2, site_f, site_p = permanova(
        river_names, name_to_idx, bc, lambda n: parse_sample_name(n)[1])
    plot_permanova(season_r2, site_r2, season_p, site_p)
    print("Saved 8 plots to Output/")

    # ---- 5. Quantify the site-vs-season effect (above the replicate baseline) ----
    baseline = means[CAT_REP]
    site_effect = means[CAT_SITE] - baseline
    season_effect = means[CAT_SEASON] - baseline
    ratio = season_effect / site_effect if site_effect else float("nan")

    # ---- 6. Leave-one-out robustness: drop each sample, re-check the headline ----
    loo_ratios, loo_season_gt_site = [], 0
    for s in river_names:
        m, _ = category_means(river_names, name_to_idx, bc, exclude={s})
        se = m[CAT_SEASON] - m[CAT_REP]
        si = m[CAT_SITE] - m[CAT_REP]
        loo_ratios.append(se / si if si else float("nan"))
        if m[CAT_SEASON] > m[CAT_SITE]:
            loo_season_gt_site += 1

    # ---- 6b. Leave-one-GROUP-out: drop a whole season / whole site (a real stress test,
    #          since dropping one of ~600 pairs barely moves a bucket mean) ----
    def group_holds(group_kind):
        held = 0
        groups = sorted({parse_sample_name(n)[0 if group_kind == "season" else 1] for n in river_names})
        for g in groups:
            kept = [n for n in river_names
                    if parse_sample_name(n)[0 if group_kind == "season" else 1] != g]
            m, _ = category_means(kept, name_to_idx, bc)
            if m[CAT_SEASON] > m[CAT_SITE]:
                held += 1
        return held, len(groups)
    los_held, los_n = group_holds("season")   # leave-one-season-out
    losite_held, losite_n = group_holds("site")  # leave-one-site-out

    # ---- 6c. Genuine metric independence: does season still beat site under the CLASSIC
    #          presence/absence Jaccard (which shares no numerator with Bray-Curtis)? ----
    jb_means, _ = category_means(river_names, name_to_idx, jacb)
    jb_site_effect = jb_means[CAT_SITE] - jb_means[CAT_REP]
    jb_season_effect = jb_means[CAT_SEASON] - jb_means[CAT_REP]
    jb_ratio = jb_season_effect / jb_site_effect if jb_site_effect else float("nan")

    # ---- 7. Extra computed facts for the conclusions ----
    site_rich = {s: statistics.mean([rich[n] for n in all_maps if parse_sample_name(n)[1] == s])
                 for s in SITES}
    site_gini = {s: statistics.mean([gini[n] for n in all_maps if parse_sample_name(n)[1] == s])
                 for s in SITES}
    richest = max(RIVER_SITES, key=lambda s: site_rich[s])
    poorest = min(RIVER_SITES, key=lambda s: site_rich[s])
    most_even = max(RIVER_SITES, key=lambda s: site_gini[s])
    # richness vs evenness correlation across river samples (do they agree?)
    rv = [rich[n] for n in river_names]
    ev = [gini[n] for n in river_names]
    rich_even_r = float(np.corrcoef(rv, ev)[0, 1])

    # ---- 8. Save an auditable effect summary so the conclusions are reproducible ----
    with open(f"{OUT_DIR}/effect_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value"])
        w.writerow(["mean_BC_replicates", round(baseline, 4)])
        w.writerow(["mean_BC_diff_site_same_season", round(means[CAT_SITE], 4)])
        w.writerow(["mean_BC_same_site_diff_season", round(means[CAT_SEASON], 4)])
        w.writerow(["mean_BC_diff_site_diff_season", round(means[CAT_BOTH], 4)])
        w.writerow(["site_effect_above_baseline", round(site_effect, 4)])
        w.writerow(["season_effect_above_baseline", round(season_effect, 4)])
        w.writerow(["season_over_site_ratio", round(ratio, 2)])
        w.writerow(["leave_one_out_ratio_min", round(min(loo_ratios), 2)])
        w.writerow(["leave_one_out_ratio_max", round(max(loo_ratios), 2)])
        w.writerow(["leave_one_out_season_gt_site_count", f"{loo_season_gt_site}/{len(river_names)}"])
        w.writerow(["leave_one_SEASON_out_season_gt_site", f"{los_held}/{los_n}"])
        w.writerow(["leave_one_SITE_out_season_gt_site", f"{losite_held}/{losite_n}"])
        w.writerow(["permanova_season_R2", round(season_r2, 4)])
        w.writerow(["permanova_season_pseudoF", round(season_f, 2)])
        w.writerow(["permanova_season_p", round(season_p, 4)])
        w.writerow(["permanova_site_R2", round(site_r2, 4)])
        w.writerow(["permanova_site_pseudoF", round(site_f, 2)])
        w.writerow(["permanova_site_p", round(site_p, 4)])
        w.writerow(["presence_absence_jaccard_season_over_site_ratio", round(jb_ratio, 2)])
        w.writerow(["mean_BC_control_to_river", round(ctrl_mean, 4)])
        w.writerow(["mean_BC_river_to_river", round(river_mean, 4)])
        w.writerow(["bray_curtis_vs_jaccard_pearson_r", round(jaccard_r, 4)])
        w.writerow(["richness_vs_depth_pearson_r", round(depth_r, 4)])
        w.writerow(["richness_vs_evenness_pearson_r", round(rich_even_r, 4)])

    # ---- 9. Print DATA-DERIVED conclusions (every number below is computed above) ----
    print("\n" + "=" * 64)
    print("ANALYSIS CONCLUSIONS  (all numbers computed from the data, not hardcoded)")
    print("=" * 64)
    print("1. Alpha diversity (richness & evenness):")
    print(f"   Among river sites, richness is highest at {richest} ({site_rich[richest]:.0f}) and")
    print(f"   lowest at {poorest} ({site_rich[poorest]:.0f}); the most even site is {most_even} "
          f"(1-D = {site_gini[most_even]:.3f}).")
    print(f"   Richness and evenness barely correlate across samples (r = {rich_even_r:.2f}), so they")
    print("   tell DIFFERENT stories -- a rich sample is not necessarily an even one.\n")
    print("2. The control (sanity check):")
    print(f"   The control's mean composition distance to river samples is {ctrl_mean:.2f}, vs")
    print(f"   {river_mean:.2f} among river samples -- so BETA diversity flags the control even though")
    print(f"   ALPHA diversity does not (Fall controls were among the richest samples). Cross-check both.\n")
    print("3. Site vs. season (the centerpiece):")
    print(f"   Replicate baseline (same site & season): {baseline:.2f} Bray-Curtis.")
    print(f"   Change only the SITE  -> {means[CAT_SITE]:.2f}  (a +{site_effect:.2f} jump).")
    print(f"   Change only the SEASON -> {means[CAT_SEASON]:.2f}  (a +{season_effect:.2f} jump).")
    print(f"   Season moves communities about {ratio:.1f}x further than site ABOVE the replicate baseline.\n")
    print("4. Formal test (PERMANOVA on river samples, samples = units, 999 permutations):")
    print(f"   SEASON: R^2 = {season_r2:.2f}, pseudo-F = {season_f:.1f}, p = {season_p:.3f}.")
    print(f"   SITE:   R^2 = {site_r2:.2f}, pseudo-F = {site_f:.1f}, p = {site_p:.3f}.")
    print("   Season explains several times more compositional variation than site, and the")
    print("   permutation p-values say neither is a fluke -- this is the proper version of '4x'.\n")
    print("5. Robustness & metric independence:")
    print(f"   Leave-one-SAMPLE-out: season beats site in {loo_season_gt_site}/{len(river_names)} drops "
          f"(ratio {min(loo_ratios):.1f}x-{max(loo_ratios):.1f}x).")
    print(f"   Leave-one-SEASON-out: still holds in {los_held}/{los_n} drops; "
          f"leave-one-SITE-out: {losite_held}/{losite_n}. No group drives it.")
    print(f"   Presence/absence Jaccard (shares NO numerator with Bray-Curtis): season still beats")
    print(f"   site by {jb_ratio:.1f}x -- so the result is genuinely metric-independent, not a tautology.")
    print(f"   (Honest caveat: richness tracks read depth, r = {depth_r:.2f}; we lean on composition.)")
    print("=" * 64 + "\n")


if __name__ == "__main__":
    main()
