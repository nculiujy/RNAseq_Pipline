"""
5_2_circos.py  —  Circos-style circular plot for DEG chromosomal distribution
Usage:
    python 5_2_circos.py --deg <DEG_result.csv> --bed <mRNA.bed>
                         --outdir <output_dir> [--padj 0.05] [--lfc 1]
"""
import argparse, os, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--deg",          required=True, help="DEG result CSV")
    p.add_argument("--bed",          required=True, help="mRNA BED12 file for gene coordinates")
    p.add_argument("--outdir",       required=True, help="Output directory")
    p.add_argument("--pvalue_type",  default="padj", choices=["pvalue","padj"])
    p.add_argument("--pvalue_cut",   type=float, default=0.05)
    p.add_argument("--lfc",          type=float, default=1.0)
    return p.parse_args()

def load_gene_coords(bed_path):
    """Return dict: gene_id -> (chrom, midpoint)"""
    coords = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            coords[name.upper()] = (chrom, (start + end) / 2)
    return coords

def chrom_sort_key(c):
    c = c.replace("Chr", "").replace("chr", "")
    try: return (0, int(c))
    except ValueError: return (1, c)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Load DEG
    deg = pd.read_csv(args.deg)
    sig = deg[(deg[args.pvalue_type] < args.pvalue_cut) & (deg["log2FoldChange"].abs() >= args.lfc)].copy()
    sig["gene_upper"] = sig["SYMBOL"].str.upper()
    print(f"Significant DEGs: {len(sig)}")

    # Load BED coordinates
    coords = load_gene_coords(args.bed)

    # Map DEGs to chromosomes — try gene_id (AGI locus) first, then SYMBOL
    chrom_up   = defaultdict(list)
    chrom_down = defaultdict(list)
    chrom_sizes = defaultdict(int)

    for _, row in sig.iterrows():
        gene = row["gene_upper"]
        # Also try gene_id column (AGI locus like AT1G01010)
        gene_id_upper = str(row.get("gene_id", "")).upper() if "gene_id" in sig.columns else ""
        hit = coords.get(gene) or coords.get(gene_id_upper)
        if not hit: continue
        chrom, mid = hit
        if row["log2FoldChange"] > 0:
            chrom_up[chrom].append(mid)
        else:
            chrom_down[chrom].append(mid)

    # Chromosome sizes from BED
    with open(args.bed) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 3: continue
            chrom, end = parts[0], int(parts[2])
            if end > chrom_sizes[chrom]:
                chrom_sizes[chrom] = end

    chroms = sorted(chrom_sizes.keys(), key=chrom_sort_key)
    if not chroms:
        print("No chromosome data found. Exiting.")
        return

    # ── Draw Circos ───────────────────────────────────────────────────────────
    n_chr = len(chroms)
    GAP = 0.04          # gap between chromosomes in radians
    total_len = sum(chrom_sizes[c] for c in chroms)
    total_arc = 2 * np.pi - n_chr * GAP

    # Compute angle ranges per chromosome
    chr_angles = {}   # chrom -> (start_angle, end_angle)
    cur = 0.0
    for c in chroms:
        span = (chrom_sizes[c] / total_len) * total_arc
        chr_angles[c] = (cur, cur + span)
        cur += span + GAP

    def pos_to_angle(chrom, pos):
        s, e = chr_angles[chrom]
        return s + (pos / chrom_sizes[chrom]) * (e - s)

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_axis_off()

    R_CHR   = 1.00   # outer chromosome ring
    R_INNER = 0.85   # inner edge of chromosome band
    R_UP    = 0.80   # up-regulated scatter ring
    R_DOWN  = 0.72   # down-regulated scatter ring

    CHR_COLORS = plt.cm.Set3(np.linspace(0, 1, n_chr))

    for i, c in enumerate(chroms):
        s, e = chr_angles[c]
        theta = np.linspace(s, e, 200)
        # Chromosome band
        ax.fill_between(theta, R_INNER, R_CHR, color=CHR_COLORS[i], alpha=0.85, zorder=2)
        # Label
        mid_angle = (s + e) / 2
        lbl = c.replace("Chr", "Chr ").replace("chr", "chr ")
        ax.text(mid_angle, R_CHR + 0.07, lbl, ha="center", va="center",
                fontsize=8, fontweight="bold",
                rotation=np.degrees(mid_angle) - 90 if mid_angle < np.pi else np.degrees(mid_angle) + 90)

        # Up-regulated dots
        for pos in chrom_up.get(c, []):
            a = pos_to_angle(c, pos)
            r = R_UP + np.random.uniform(-0.03, 0.03)
            ax.scatter(a, r, color="#E74C3C", s=12, alpha=0.7, zorder=3)

        # Down-regulated dots
        for pos in chrom_down.get(c, []):
            a = pos_to_angle(c, pos)
            r = R_DOWN + np.random.uniform(-0.03, 0.03)
            ax.scatter(a, r, color="#3498DB", s=12, alpha=0.7, zorder=3)

    # Ring labels
    ax.text(0, R_UP,   "UP",   ha="center", va="center", fontsize=7, color="#E74C3C")
    ax.text(0, R_DOWN, "DOWN", ha="center", va="center", fontsize=7, color="#3498DB")

    # Legend
    patches = [
        mpatches.Patch(color="#E74C3C", label=f"Up-regulated (n={sum(len(v) for v in chrom_up.values())})"),
        mpatches.Patch(color="#3498DB", label=f"Down-regulated (n={sum(len(v) for v in chrom_down.values())})"),
    ]
    ax.legend(handles=patches, loc="lower center", bbox_to_anchor=(0.5, -0.05),
              frameon=False, fontsize=9)

    prefix = os.path.splitext(os.path.basename(args.deg))[0].replace("_DEG_result", "")
    plt.title(f"DEG Chromosomal Distribution\n{prefix}", fontsize=11, pad=20)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(args.outdir, f"{prefix}_circos.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved Circos plot: {prefix}_circos.png/pdf")

if __name__ == "__main__":
    main()
