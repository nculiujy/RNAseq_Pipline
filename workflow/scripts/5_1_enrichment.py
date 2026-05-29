"""
5_1_enrichment.py  —  GO & KEGG enrichment analysis for DEG results
Usage:
    python 5_1_enrichment.py --deg <DEG_result.csv> --species <tair|human|mouse>
                             --outdir <output_dir> [--padj 0.05] [--lfc 1]
"""
import argparse, os, sys, gzip, shutil, urllib.request
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--deg",     required=True, help="DEG result CSV (from 4_1_RNAseq_DEseq.R)")
    p.add_argument("--species", required=True, choices=["tair","human","mouse"],
                   help="Species: tair | human | mouse")
    p.add_argument("--outdir",  required=True, help="Output directory")
    p.add_argument("--padj",    type=float, default=0.05, help="DEG padj threshold")
    p.add_argument("--lfc",     type=float, default=1.0,  help="|log2FC| threshold")
    p.add_argument("--top_n",   type=int,   default=10,   help="Top N terms per namespace")
    return p.parse_args()

# ── Download helpers ──────────────────────────────────────────────────────────
def _download(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f)

def ensure_file(path, url, gz=False):
    if os.path.exists(path):
        return
    print(f"Downloading {os.path.basename(path)} …")
    tmp = path + ".gz" if gz else path
    _download(url, tmp)
    if gz:
        with gzip.open(tmp, "rb") as fi, open(path, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        os.remove(tmp)

# ── GO enrichment ─────────────────────────────────────────────────────────────
OBO_URL = "https://release.geneontology.org/2024-01-17/ontology/go-basic.obo"
GAF_URLS = {
    "tair":  "https://release.geneontology.org/2024-01-17/annotations/tair.gaf.gz",
    "human": "https://release.geneontology.org/2024-01-17/annotations/goa_human.gaf.gz",
    "mouse": "https://release.geneontology.org/2024-01-17/annotations/mgi.gaf.gz",
}

def parse_obo(path):
    terms, cur = {}, {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                if cur.get("id"): terms[cur["id"]] = cur
                cur = {}
            elif line.startswith("id: GO:"): cur["id"] = line[4:]
            elif line.startswith("name: "):  cur["name"] = line[6:]
            elif line.startswith("namespace: "): cur["ns"] = line[11:]
            elif line.startswith("is_obsolete: true"): cur["obsolete"] = True
    if cur.get("id"): terms[cur["id"]] = cur
    return {k: v for k, v in terms.items() if not v.get("obsolete")}

def parse_gaf(path, obo):
    gene2go, go2genes = {}, {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("!"): continue
            parts = line.strip().split("\t")
            if len(parts) < 10: continue
            for col in (1, 2):          # DB_Object_ID and Symbol
                gene = parts[col].upper()
                go_id = parts[4]
                if go_id not in obo: continue
                gene2go.setdefault(gene, set()).add(go_id)
                go2genes.setdefault(go_id, set()).add(gene)
    return gene2go, go2genes

def run_go_fisher(study_genes, gene2go, go2genes, obo):
    NS_MAP = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}
    population = set(gene2go.keys())
    study_in_pop = study_genes & population
    pop_n, study_n = len(population), len(study_in_pop)
    rows = []
    for go_id, ann_genes in go2genes.items():
        k = len(study_in_pop & ann_genes)
        if k < 2: continue
        K = len(ann_genes & population)
        table = [[k, study_n - k], [K - k, pop_n - study_n - K + k]]
        _, p = fisher_exact(table, alternative="greater")
        term = obo[go_id]
        rows.append({"GO_ID": go_id, "Term": term["name"],
                     "Namespace": NS_MAP.get(term["ns"], term["ns"]),
                     "Count": k, "pvalue": p})
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df["p.adjust"] = padj
    df["-log10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    return df.sort_values("pvalue")

def plot_go(result, top_n, outdir, prefix):
    sig = result[result["pvalue"] < 0.05].copy()
    if sig.empty:
        print("No significant GO terms (p<0.05). Skipping GO plot.")
        return
    frames = [sig[sig["Namespace"] == ns].nsmallest(top_n, "pvalue")
              for ns in ["BP", "MF", "CC"]]
    plot_df = pd.concat([f for f in frames if not f.empty])
    plot_df["ns_order"] = plot_df["Namespace"].map({"BP": 0, "MF": 1, "CC": 2})
    plot_df = plot_df.sort_values(["ns_order", "-log10p"], ascending=[True, False]).iloc[::-1]

    COLOR = {"BP": "#E8A838", "MF": "#6BAED6", "CC": "#31A354"}
    LABEL = {"BP": "Biological process", "MF": "Molecular function", "CC": "Cellular component"}
    fig, ax = plt.subplots(figsize=(14, max(6, len(plot_df) * 0.42)))
    colors = [COLOR[ns] for ns in plot_df["Namespace"]]
    ax.barh(range(len(plot_df)), plot_df["-log10p"].values,
            color=colors, edgecolor="none", height=0.7)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["Term"].values, fontsize=9)
    ax.set_xlabel("−log10(P-value)", fontsize=10)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)
    patches = [mpatches.Patch(color=COLOR[ns], label=LABEL[ns]) for ns in ["BP", "MF", "CC"]]
    ax.legend(handles=patches, title="Term type", frameon=False,
              loc="lower right", bbox_to_anchor=(1.35, 0.0))
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(outdir, f"{prefix}_GO_barplot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved GO plot: {prefix}_GO_barplot.png/pdf")

# ── KEGG enrichment ───────────────────────────────────────────────────────────
KEGG_ORG = {"tair": "ath", "human": "hsa", "mouse": "mmu"}

def fetch_kegg_pathway_genes(org, cache_dir):
    """Fetch KEGG pathway→gene mapping via KEGG REST API."""
    list_file = os.path.join(cache_dir, f"kegg_{org}_pathways.txt")
    if not os.path.exists(list_file):
        print(f"Fetching KEGG pathway list for {org} …")
        url = f"https://rest.kegg.jp/list/pathway/{org}"
        _download(url, list_file)

    pathways = {}
    with open(list_file) as f:
        for line in f:
            pid = line.split("\t")[0].strip()   # e.g. path:ath00010
            pathways[pid] = line.split("\t")[1].strip() if "\t" in line else pid

    pw2genes = {}
    for pid, name in pathways.items():
        gene_file = os.path.join(cache_dir, f"{pid.replace(':', '_')}_genes.txt")
        if not os.path.exists(gene_file):
            url = f"https://rest.kegg.jp/link/genes/{pid}"
            try:
                _download(url, gene_file)
            except Exception:
                open(gene_file, "w").close()
        genes = set()
        with open(gene_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    # gene id like "ath:AT1G01010" → "AT1G01010"
                    genes.add(parts[1].split(":")[-1].upper())
        if genes:
            pw2genes[pid] = {"name": name, "genes": genes}
    return pw2genes

def run_kegg_fisher(study_genes, pw2genes):
    all_kegg_genes = set()
    for v in pw2genes.values():
        all_kegg_genes |= v["genes"]
    population = all_kegg_genes
    study_in_pop = study_genes & population
    pop_n, study_n = len(population), len(study_in_pop)
    if study_n == 0:
        return pd.DataFrame()
    rows = []
    for pid, info in pw2genes.items():
        ann_genes = info["genes"]
        k = len(study_in_pop & ann_genes)
        if k < 2: continue
        K = len(ann_genes & population)
        table = [[k, study_n - k], [K - k, pop_n - study_n - K + k]]
        _, p = fisher_exact(table, alternative="greater")
        rows.append({"Pathway_ID": pid, "Term": info["name"],
                     "Count": k, "pvalue": p, "GeneRatio": k / study_n})
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    _, padj, _, _ = multipletests(df["pvalue"], method="fdr_bh")
    df["p.adjust"] = padj
    df["-log10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    return df.sort_values("pvalue")

def plot_kegg(result, top_n, outdir, prefix):
    sig = result[result["pvalue"] < 0.05].nsmallest(top_n * 3, "pvalue").copy()
    if sig.empty:
        print("No significant KEGG pathways (p<0.05). Skipping KEGG plot.")
        return
    sig = sig.sort_values("-log10p").tail(min(top_n * 3, len(sig)))

    fig, ax = plt.subplots(figsize=(12, max(6, len(sig) * 0.45)))
    sc = ax.scatter(sig["-log10p"], range(len(sig)),
                    c=sig["p.adjust"], cmap="RdYlGn_r",
                    s=sig["Count"] * 20, vmin=0, vmax=0.05,
                    edgecolors="grey", linewidths=0.5)
    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(sig["Term"].values, fontsize=9)
    ax.set_xlabel("−log10(P-value)", fontsize=10)
    ax.set_title("KEGG Pathway Enrichment", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("p.adjust", fontsize=9)
    # size legend
    for sz in [5, 20, 50]:
        ax.scatter([], [], s=sz * 20, c="grey", alpha=0.6, label=f"n={sz}")
    ax.legend(title="Gene count", frameon=False, loc="lower right",
              bbox_to_anchor=(1.35, 0.0))
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(outdir, f"{prefix}_KEGG_dotplot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved KEGG plot: {prefix}_KEGG_dotplot.png/pdf")

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cache_dir = os.path.join(args.outdir, ".cache")
    os.makedirs(cache_dir, exist_ok=True)

    # Load DEG
    deg = pd.read_csv(args.deg)
    # Column names from 4_1_RNAseq_DEseq.R output: SYMBOL, log2FoldChange, padj
    sig_mask = (deg["padj"] < args.padj) & (deg["log2FoldChange"].abs() >= args.lfc)
    study_genes = set(deg.loc[sig_mask, "SYMBOL"].dropna().str.upper())
    print(f"Significant DEGs: {len(study_genes)}")
    if not study_genes:
        print("No DEGs passed threshold. Exiting.")
        sys.exit(0)

    prefix = os.path.splitext(os.path.basename(args.deg))[0].replace("_DEG_result", "")

    # ── GO ────────────────────────────────────────────────────────────────────
    obo_path = os.path.join(cache_dir, "go-basic.obo")
    gaf_path = os.path.join(cache_dir, f"{args.species}.gaf")
    ensure_file(obo_path, OBO_URL)
    ensure_file(gaf_path, GAF_URLS[args.species], gz=True)

    obo = parse_obo(obo_path)
    gene2go, go2genes = parse_gaf(gaf_path, obo)
    go_result = run_go_fisher(study_genes, gene2go, go2genes, obo)
    if not go_result.empty:
        go_result.to_csv(os.path.join(args.outdir, f"{prefix}_GO_enrichment.csv"), index=False)
        plot_go(go_result, args.top_n, args.outdir, prefix)

    # ── KEGG ──────────────────────────────────────────────────────────────────
    org = KEGG_ORG.get(args.species)
    if org:
        pw2genes = fetch_kegg_pathway_genes(org, cache_dir)
        kegg_result = run_kegg_fisher(study_genes, pw2genes)
        if not kegg_result.empty:
            kegg_result.to_csv(os.path.join(args.outdir, f"{prefix}_KEGG_enrichment.csv"), index=False)
            plot_kegg(kegg_result, args.top_n, args.outdir, prefix)
    else:
        print(f"No KEGG organism code for species '{args.species}'. Skipping KEGG.")

if __name__ == "__main__":
    main()
