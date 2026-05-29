import os
import json
import pandas as pd

def get_deseq_config(metadata_csv):
    if not os.path.exists(metadata_csv):
        return "{}"
    df = pd.read_csv(metadata_csv)
    config_dict = {}
    for _, row in df.iterrows():
        config_dict[str(row['sample_name'])] = str(row['group'])
    # format as string with single quotes for command line
    return json.dumps(config_dict).replace('"', "'")

def get_latest_merge_dir(species):
    flag_file = os.path.join("result", species, "3_Align_Filter", "Merge_finished.txt")
    if not os.path.exists(flag_file):
        return None
    with open(flag_file, "r") as f:
        line = f.read().strip()
        if "Output: " in line:
            return line.split("Output: ")[1]
    return None

rule DEseq_analysis:
    input:
        merge_flag = "result/{species}/3_Align_Filter/Merge_finished.txt",
        metadata = "config/RNAseq_metadata.csv"
    output:
        flag = "result/{species}/4_DEseq/DEseq_finished.txt"
    params:
        outdir       = "result/{species}/4_DEseq",
        script_4_1   = "workflow/scripts/4_1_RNAseq_DEseq.R",
        script_4_2   = "workflow/scripts/4_2_DEseq_result_plot_pipline.R",
        pvalue_type  = config.get("deg_pvalue_type",   "padj"),
        pvalue_cut   = config.get("deg_pvalue_cutoff", 0.05),
        lfc          = config.get("deg_lfc_cutoff",    1.0)
    log:
        "logs/{species}_DEseq.log"
    run:
        import glob
        import subprocess
        
        species = wildcards.species
        merge_dir = get_latest_merge_dir(species)
        if not merge_dir:
            raise Exception("Cannot find merge directory from Merge_finished.txt")
            
        config_str = get_deseq_config(input.metadata)
        
        # Ensure output directory exists
        os.makedirs(params.outdir, exist_ok=True)
        
        # Find Count matrices (there might be multiple if multiple annotations are used, e.g. mRNA, eRNA, etc.)
        all_count_matrices = glob.glob(os.path.join(merge_dir, "*_Count_matrix.csv"))
        # Filter out Transcript level count matrix, only run DESeq2 on Gene level
        count_matrices = [f for f in all_count_matrices if "Transcript" not in os.path.basename(f)]
        hisat2_log = os.path.join(merge_dir, f"{species}_hisat2_result.log")
        
        for count_matrix in count_matrices:
            # e.g. TAIR_mRNA_TAIR10_stringtie_Count_matrix.csv
            base_name = os.path.basename(count_matrix).replace("_Count_matrix.csv", "")
            
            deg_out = os.path.join(params.outdir, f"{base_name}_DEG_result.csv")
            
            # 1. Run DESeq2
            cmd1 = [
                "Rscript", params.script_4_1,
                "--input", count_matrix,
                "--config", config_str,
                "--output", deg_out
            ]
            print(f"Running DESeq2 for {base_name}...")
            subprocess.run(cmd1, check=True)
            
            # 2. Run Plotting
            outdir_table = os.path.join(params.outdir, base_name, "tables")
            outdir_plot = os.path.join(params.outdir, base_name, "plots")
            os.makedirs(outdir_table, exist_ok=True)
            os.makedirs(outdir_plot, exist_ok=True)
            
            # Need bed file. We can infer it from species or config.
            # But where is the bed file?
            # For TAIR10, maybe we don't have a bed file, or we should generate one from GTF?
            # 4_2_DEseq_result_plot_pipline.R requires --bed
            # Let's get the GTF path from config and generate a simple BED if not exists.
            gtf_file = config["projects"][0]["gtf_file"] # Wait, this is hardcoded to project 0.
            # Better find the correct project:
            proj_gtf = None
            for proj in config["projects"]:
                if proj["species"] == species:
                    proj_gtf = proj["gtf_file"]
                    break
            
            bed_file = proj_gtf.replace(".gtf", ".bed")
            if not os.path.exists(bed_file):
                # Generate BED from GTF
                # format: chr start end gene_id score strand
                # grep -v "^#" gtf | awk '$3=="gene"' | awk '{print $1"\t"$4"\t"$5"\t"$10"\t.\t"$7}' 
                # gene_id might be $10, let's be careful.
                print(f"Generating BED file from {proj_gtf}...")
                with open(proj_gtf, "r") as f_in, open(bed_file, "w") as f_out:
                    for line in f_in:
                        if line.startswith("#"): continue
                        parts = line.split('\t')
                        if parts[2] == 'transcript' or parts[2] == 'gene':
                            # extract gene_id
                            attrs = parts[8]
                            gene_id = "unknown"
                            if 'gene_id "' in attrs:
                                gene_id = attrs.split('gene_id "')[1].split('"')[0]
                            # bed format: chr start end gene_id score strand
                            # bed is 0-indexed, gtf is 1-indexed. start-1
                            f_out.write(f"{parts[0]}\t{int(parts[3])-1}\t{parts[4]}\t{gene_id}\t.\t{parts[6]}\n")
            
            cmd2 = [
                "Rscript", params.script_4_2,
                "--log", hisat2_log,
                "--deg", deg_out,
                "--outdir_table", outdir_table,
                "--outdir_plot", outdir_plot,
                "--bed", bed_file,
                "--pvalue_type", str(params.pvalue_type),
                "--pvalue_cut",  str(params.pvalue_cut),
                "--lfc",         str(params.lfc)
            ]
            print(f"Running Plots for {base_name}...")
            subprocess.run(cmd2, check=True)
            
        with open(output.flag, "w") as f:
            f.write(f"DEseq finished for {species}\n")
