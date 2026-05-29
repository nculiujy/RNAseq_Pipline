import os, sys

SPECIES_MAP = {"TAIR": "tair", "homo": "human", "mm": "mouse"}

rule enrichment:
    input:
        deseq_flag = "result/{species}/4_DEseq/DEseq_finished.txt"
    output:
        flag = "result/{species}/5_enrichment/enrichment_finished.txt"
    params:
        outdir      = "result/{species}/5_enrichment",
        script_1    = "workflow/scripts/5_1_enrichment.py",
        script_2    = "workflow/scripts/5_2_circos.py",
        pvalue_type = config.get("deg_pvalue_type",   "padj"),
        pvalue_cut  = config.get("deg_pvalue_cutoff", 0.05),
        lfc         = config.get("deg_lfc_cutoff",    1.0),
    log:
        "logs/{species}_enrichment.log"
    run:
        import glob, subprocess

        species   = wildcards.species
        sp_code   = SPECIES_MAP.get(species, "tair")
        deseq_dir = os.path.join("result", species, "4_DEseq")
        os.makedirs(params.outdir, exist_ok=True)

        # Find BED file from config
        bed_file = None
        for proj in config["projects"]:
            if proj["species"] == species:
                bed_file = proj["gtf_file"].replace(".gtf", ".bed")
                break

        deg_files = glob.glob(os.path.join(deseq_dir, "*_DEG_result.csv"))
        for deg in deg_files:
            base = os.path.splitext(os.path.basename(deg))[0].replace("_DEG_result", "")
            sub  = os.path.join(params.outdir, base)
            os.makedirs(sub, exist_ok=True)

            # 5_1: GO + KEGG enrichment
            subprocess.run([
                sys.executable, params.script_1,
                "--deg", deg, "--species", sp_code,
                "--outdir", sub,
                "--pvalue_type", params.pvalue_type,
                "--pvalue_cut",  str(params.pvalue_cut),
                "--lfc",         str(params.lfc)
            ], check=True)

            # 5_2: Circos plot (requires BED)
            if bed_file and os.path.exists(bed_file):
                subprocess.run([
                    sys.executable, params.script_2,
                    "--deg", deg, "--bed", bed_file,
                    "--outdir", sub,
                    "--pvalue_type", params.pvalue_type,
                    "--pvalue_cut",  str(params.pvalue_cut),
                    "--lfc",         str(params.lfc)
                ], check=True)
            else:
                print(f"BED file not found ({bed_file}), skipping Circos plot.")

        with open(output.flag, "w") as f:
            f.write(f"Enrichment finished for {species}\n")
