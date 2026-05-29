"""
5_2_circos.py  —  Circos-style circular plot for DEG chromosomal distribution
Style: colored chromosome arcs + gene tick labels + inner dot track + Bezier links
Usage:
    python 5_2_circos.py --deg <DEG_result.csv> --bed <mRNA.bed>
                         --outdir <output_dir> [--pvalue_type padj] [--pvalue_cut 0.05] [--lfc 1]
                         [--max_genes 100]
"""
import argparse, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patheffects as pe
from collections import defaultdict

CHR_COLORS = {
    "Chr1": "#E74C3C", "Chr2": "#27AE60", "Chr3": "#3498DB",
    "Chr4": "#F39C12", "Chr5": "#9B59B6",
    "1": "#E74C3C", "2": "#27AE60", "3": "#3498DB",
    "4": "#F39C12", "5": "#9B59B6",
}
DEFAULT_COLOR = "#95A5A6"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--deg",          required=True)
    p.add_argument("--bed",          required=True)
    p.add_argument("--outdir",       required=True)
    p.add_argument("--pvalue_type",  default="padj", choices=["pvalue","padj"])
    p.add_argument("--pvalue_cut",   type=float, default=0.05)
    p.add_argument("--lfc",          type=float, default=1.0)
    p.add_argument("--max_genes",    type=int,   default=100,
                   help="Max DEGs to show links for (top by |lfc|)")
    return p.parse_args()

def load_gene_coords(bed_path):
    """Return dict: locus_id_upper -> (chrom, midpoint)"""
    coords = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            mid = (start + end) / 2
            coords[name.upper()] = (chrom, mid)
            locus = name.split(".")[0].upper()
            if locus not in coords:
                coords[locus] = (chrom, mid)
    return coords

def chrom_sort_key(c):
    s = c.replace("Chr","").replace("chr","")
    try: return (0, int(s))
    except ValueError: return (1, s)

def angle_for(chrom, pos, chr_angles, chrom_sizes):
    s, e = chr_angles[chrom]
    return s + (pos / chrom_sizes[chrom]) * (e - s)

def polar_to_xy(r, theta):
    return r * np.cos(theta), r * np.sin(theta)

def bezier_link(ax, a1, a2, r_inner, color, alpha=0.4, lw=0.6):
    """Draw a cubic Bezier ribbon between two angles on the circle."""
    x1, y1 = polar_to_xy(r_inner, a1)
    x2, y2 = polar_to_xy(r_inner, a2)
    # Control points pulled toward center
    ctrl_r = r_inner * 0.15
    cx1, cy1 = polar_to_xy(ctrl_r, a1)
    cx2, cy2 = polar_to_xy(ctrl_r, a2)
    verts = [(x1,y1),(cx1,cy1),(cx2,cy2),(x2,y2)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path = Path(verts, codes)
    patch = mpatches.FancyArrowPatch.__new__(mpatches.FancyArrowPatch)
    from matplotlib.patches import PathPatch
    pp = PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=alpha, zorder=2)
    ax.add_patch(pp)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    deg = pd.read_csv(args.deg)
    sig = deg[(deg[args.pvalue_type] < args.pvalue_cut) &
              (deg["log2FoldChange"].abs() >= args.lfc)].copy()
    print(f"Significant DEGs: {len(sig)}")
    if sig.empty:
        print("No DEGs. Exiting.")
        return

    # Sort by |lfc| and take top max_genes
    sig = sig.reindex(sig["log2FoldChange"].abs().sort_values(ascending=False).index)
    sig = sig.head(args.max_genes)

    coords = load_gene_coords(args.bed)

    # Map each DEG to (chrom, pos)
    gene_pos = {}
    for _, row in sig.iterrows():
        gid = str(row.get("gene_id","")).upper()
        sym = str(row.get("SYMBOL","")).upper()
        hit = coords.get(gid) or coords.get(sym)
        if hit:
            gene_pos[row.name] = (hit[0], hit[1], row["log2FoldChange"])

    if not gene_pos:
        print("No DEGs could be mapped to BED coordinates. Exiting.")
        return

    # Chromosome sizes from BED
    chrom_sizes = defaultdict(int)
    with open(args.bed) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 3: continue
            chrom, end = parts[0], int(parts[2])
            if end > chrom_sizes[chrom]:
                chrom_sizes[chrom] = end

    chroms = sorted(chrom_sizes.keys(), key=chrom_sort_key)

    # Angle layout (in standard math coords, 0=right, CCW)
    GAP_DEG = 3.0
    GAP = np.radians(GAP_DEG)
    total_len = sum(chrom_sizes[c] for c in chroms)
    total_arc = 2*np.pi - len(chroms)*GAP

    chr_angles = {}
    cur = np.pi/2  # start at top
    for c in chroms:
        span = (chrom_sizes[c]/total_len)*total_arc
        chr_angles[c] = (cur, cur - span)
        cur = cur - span - GAP

    # ── Draw ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)

    R_OUT   = 1.20   # outer edge of chromosome arc
    R_IN    = 1.10   # inner edge of chromosome arc
    R_TICK  = 1.08   # tick marks for genes
    R_DOT   = 1.04   # dot track (beige band inner edge)
    R_BAND_OUT = 1.09
    R_BAND_IN  = 0.98
    R_LINK  = 0.97   # Bezier link attachment radius

    # Draw beige inner band
    theta_all = np.linspace(0, 2*np.pi, 500)
    ax.fill_between(np.cos(theta_all)*R_BAND_IN, np.sin(theta_all)*R_BAND_IN,
                    np.cos(theta_all)*R_BAND_OUT, np.sin(theta_all)*R_BAND_OUT,
                    color="#D4C5A9", alpha=0.5, zorder=0)
    # Actually use a proper annulus
    from matplotlib.patches import Wedge, Circle
    band = mpatches.Annulus((0,0), R_BAND_OUT, R_BAND_OUT-R_BAND_IN,
                             color="#D4C5A9", alpha=0.6, zorder=0)
    ax.add_patch(band)

    # Draw chromosome arcs and labels
    for c in chroms:
        s, e = chr_angles[c]
        color = CHR_COLORS.get(c, DEFAULT_COLOR)
        theta = np.linspace(e, s, 300)
        # Outer arc
        ax.plot(np.cos(theta)*R_OUT, np.sin(theta)*R_OUT, color=color, lw=4, solid_capstyle="butt", zorder=3)
        # Inner arc boundary
        ax.plot(np.cos(theta)*R_IN, np.sin(theta)*R_IN, color=color, lw=1, alpha=0.4, zorder=3)
        # End caps
        for ang in (s, e):
            ax.plot([np.cos(ang)*R_IN, np.cos(ang)*R_OUT],
                    [np.sin(ang)*R_IN, np.sin(ang)*R_OUT],
                    color=color, lw=1.5, zorder=3)
        # Chromosome label
        mid_a = (s+e)/2
        lx, ly = np.cos(mid_a)*1.32, np.sin(mid_a)*1.32
        lbl = c.replace("Chr","").replace("chr","")
        ax.text(lx, ly, lbl, ha="center", va="center", fontsize=11,
                fontweight="bold", color=color)

    # Draw DEG ticks and gene labels on the chromosome ring
    mapped_genes = []
    for idx, (row_idx, (chrom, pos, lfc)) in enumerate(gene_pos.items()):
        a = angle_for(chrom, pos, chr_angles, chrom_sizes)
        # Tick
        ax.plot([np.cos(a)*R_IN, np.cos(a)*R_TICK],
                [np.sin(a)*R_IN, np.sin(a)*R_TICK],
                color="#333333", lw=0.8, zorder=4)
        # Gene label (outside the arc)
        lx, ly = np.cos(a)*1.14, np.sin(a)*1.14
        rot = np.degrees(a)
        if np.cos(a) < 0:
            rot += 180
        row_data = sig.loc[row_idx]
        label = str(row_data.get("SYMBOL", row_data.get("gene_id", "")))
        ax.text(lx, ly, label, ha="left" if np.cos(a)>=0 else "right",
                va="center", fontsize=4.5, rotation=rot, rotation_mode="anchor",
                color="#333333", zorder=5)
        # Dot on band
        dot_color = "#E74C3C" if lfc > 0 else "#3498DB"
        ax.scatter(np.cos(a)*R_DOT, np.sin(a)*R_DOT,
                   s=10, color=dot_color, zorder=5, linewidths=0)
        mapped_genes.append((a, lfc))

    # Draw Bezier links between all pairs of DEGs
    n = len(mapped_genes)
    for i in range(n):
        for j in range(i+1, n):
            a1, lfc1 = mapped_genes[i]
            a2, lfc2 = mapped_genes[j]
            # Color by whether both are same direction
            if lfc1 > 0 and lfc2 > 0:
                lc = "#E74C3C"
            elif lfc1 < 0 and lfc2 < 0:
                lc = "#3498DB"
            else:
                lc = "#9B59B6"
            bezier_link(ax, a1, a2, R_LINK, lc, alpha=0.25, lw=0.5)

    # Legend
    patches = [
        mpatches.Patch(color="#E74C3C", label="Up-regulated links"),
        mpatches.Patch(color="#3498DB", label="Down-regulated links"),
        mpatches.Patch(color="#9B59B6", label="Mixed links"),
        mpatches.Patch(color="#E74C3C", label="Up DEG dot"),
        mpatches.Patch(color="#3498DB", label="Down DEG dot"),
    ]
    ax.legend(handles=patches, loc="lower center", bbox_to_anchor=(0.5,-0.02),
              ncol=3, frameon=False, fontsize=8)

    prefix = os.path.splitext(os.path.basename(args.deg))[0].replace("_DEG_result","")
    plt.title(f"DEG Circos Plot — {prefix}\n(top {len(mapped_genes)} DEGs by |log2FC|)",
              fontsize=11, pad=10)
    plt.tight_layout()
    for ext in ("png","pdf"):
        plt.savefig(os.path.join(args.outdir, f"{prefix}_circos.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved Circos plot: {prefix}_circos.png/pdf  ({len(mapped_genes)} genes, {n*(n-1)//2} links)")

if __name__ == "__main__":
    main()
