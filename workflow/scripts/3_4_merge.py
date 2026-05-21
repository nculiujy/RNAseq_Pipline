
import os
import glob
import argparse
from collections import defaultdict
import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.4: Merge Quantification Results")
    parser.add_argument("--inputdir", required=True, help="Directory containing quantification results")
    parser.add_argument("--outputdir", required=True, help="Output directory for merged matrices")
    parser.add_argument("--filter_csv", required=True, help="CSV file containing QC results (alignment_quality.csv)")
    # 移除 species 参数，改为自动检测
    parser.add_argument("--gtf_base", help="GTF annotation base directory")
    return parser.parse_args()

def get_species_from_path(path):
    path_lower = path.lower()
    if "homo" in path_lower:
        return "human"
    elif "mouse" in path_lower:
        return "mouse"
    elif "tair" in path_lower:
        return "TAIR"
    return None

def get_annotation_dirs(species, gtf_base_dir):
    """Reuse annotation structure from Quant script to know which directories to look for"""
    if species == 'human':
        return [
            "mRNA/genecode/stringtie",
            "eRNA/EnhancerAtlas/stringtie",
            "eRNA/Ensembl/stringtie",
            "eRNA/FANTOM5/stringtie",
            "lncRNA/GENCODE/stringtie",
            "miRNA/miRBase/stringtie",
            "miRNA/MirGeneDB/stringtie",
            "lncRNA/NONCODE/stringtie"
        ]
    elif species == 'mouse':
        return [
            "mRNA/genecode/stringtie",
            "eRNA/EnhancerAtlas/stringtie",
            "eRNA/Ensembl/stringtie",
            "eRNA/FANTOM5/stringtie",
            "lncRNA/GENCODE/stringtie",
            "miRNA/miRBase/stringtie",
            "miRNA/MirGeneDB/stringtie",
            "lncRNA/NONCODE/stringtie"
        ]
    elif species == 'TAIR':
        return [
            "mRNA/TAIR10/stringtie"
        ]
    return []

def load_qc_pass_samples(qc_csv):
    """Load sample IDs that passed QC"""
    if not os.path.exists(qc_csv):
        print(f"Warning: QC file {qc_csv} not found. Proceeding with all samples.")
        return None
    
    try:
        passed_samples = set()
        with open(qc_csv, 'r') as f:
            lines = f.readlines()
            if not lines: return None
            
            headers = lines[0].strip().split(',')
            try:
                sid_idx = headers.index('Sample_ID')
                pass_idx = headers.index('Passed')
            except ValueError:
                # Fallback to no-QC-filtering if headers don't match
                return None
                
            for line in lines[1:]:
                parts = line.strip().split(',')
                # Handle possible missing 'Passed' column gracefully
                if len(parts) > max(sid_idx, pass_idx) and parts[pass_idx] == 'Yes':
                    passed_samples.add(parts[sid_idx])
                    
        # If no samples passed but we have lines, maybe QC file format changed, skip filtering
        if len(passed_samples) == 0:
            print(f"Warning: No passing samples found in QC file (or failed to parse). Proceeding with all samples.")
            return None
            
        print(f"Loaded {len(passed_samples)} passing samples from QC file.")
        return passed_samples
    except Exception as e:
        print(f"Error reading QC file: {e}")
        return None

def merge_expression_matrices(inputdir, outputdir, gtf_base, passed_samples):
    """
    Iterate through all subdirectories to find gene_abund.tab
    Detect species from path.
    """
    
    # We need to scan inputdir to find "homo" or "mouse" directories first
    # Or just search all gene_abund.tab and group by species/annotation type
    
    # Let's search for all gene_abund.tab files
    search_pattern = os.path.join(inputdir, "**", "gene_abund.tab")
    all_files = glob.glob(search_pattern, recursive=True)
    
    # Group files by (species, annotation_type)
    # File path structure: .../inputdir/homo/.../mRNA/genecode/stringtie/SampleID/gene_abund.tab
    
    grouped_files = defaultdict(list)
    
    for fpath in all_files:
        # Infer species from full path
        rel_path = os.path.relpath(fpath, inputdir)
        species = get_species_from_path(fpath) # use fpath to ensure species is matched
        
        if not species:
            continue
            
        # Infer annotation type from path
        # Check against known annotation dirs
        anno_dirs = get_annotation_dirs(species, gtf_base)
        matched_anno = None
        for anno in anno_dirs:
            if anno in rel_path:
                matched_anno = anno
                break
        
        if matched_anno:
            grouped_files[(species, matched_anno)].append(fpath)

    # Extract species from grouped_files
    species_set = {sp for sp, _ in grouped_files.keys()}
    if not species_set:
        print("No valid files found for merge.")
        return
    # Assuming one species per inputdir
    species = list(species_set)[0]

    # Combine hisat2 logs
    os.makedirs(outputdir, exist_ok=True)
    combined_log = os.path.join(outputdir, f"{species}_hisat2_result.log")
    with open(combined_log, "w") as out_f:
        # Search for QC_results.log
        log_pattern = os.path.join(inputdir, "**", "QC_results.log")
        log_files = glob.glob(log_pattern, recursive=True)
        for log_f in log_files:
            # Try to extract sample id from path
            # path is like result/TAIR/3_Align_Filter/SampleID/hisat2file/SampleID/QC_results.log
            # or result/TAIR/3_Align_Filter/SampleID/hisat2file/QC_results.log
            parts = log_f.split(os.sep)
            # find index of 3_Align_Filter
            try:
                idx = parts.index("3_Align_Filter")
                sample_id = parts[idx + 1]
            except ValueError:
                sample_id = "unknown"
                
            if passed_samples is not None and sample_id not in passed_samples:
                continue
                
            out_f.write(f"===== {sample_id} =====\n")
            with open(log_f, "r") as in_f:
                out_f.write(in_f.read())
            out_f.write("\n")
    print(f"  Combined HISAT2 logs into {combined_log}")

    # Process each group
    for (species, anno_subpath), files in grouped_files.items():
        print(f"Processing {species} - {anno_subpath} ({len(files)} files)")
        
        tpm_data = defaultdict(dict)
        sample_ids = set()
        
        safe_name = anno_subpath.replace('/', '_')
        
        # Prepare prepDE.py input
        os.makedirs(outputdir, exist_ok=True)
        prepde_input = os.path.join(outputdir, f"{species}_{safe_name}_prepDE_in.txt")
        with open(prepde_input, 'w') as pf:
            for fpath in files:
                sample_id = os.path.basename(os.path.dirname(fpath))
                
                if passed_samples is not None and sample_id not in passed_samples:
                    continue
                    
                sample_ids.add(sample_id)
                
                # write GTF path for prepDE.py
                gtf_path = os.path.join(os.path.dirname(fpath), "transcripts.gtf")
                if os.path.exists(gtf_path):
                    pf.write(f"{sample_id}\t{gtf_path}\n")
                
                try:
                    # Read file manually to avoid pandas/numpy dependency issues
                    with open(fpath, 'r') as f:
                        lines = f.readlines()
                        if not lines: continue
                        
                        headers = lines[0].strip().split('\t')
                        try:
                            gid_idx = headers.index('Gene ID')
                            tpm_idx = headers.index('TPM')
                        except ValueError:
                            continue
                            
                        for line in lines[1:]:
                            parts = line.strip().split('\t')
                            if len(parts) > max(gid_idx, tpm_idx):
                                gene_id = parts[gid_idx]
                                try:
                                    tpm = float(parts[tpm_idx])
                                    tpm_data[gene_id][sample_id] = tpm
                                except ValueError:
                                    pass
                except Exception as e:
                    print(f"  Error reading {fpath}: {e}")
        
        if not tpm_data:
            continue
            
        # Manually construct matrix and write to CSV
        all_genes = sorted(tpm_data.keys())
        all_samples = sorted(list(sample_ids))
        
        # Output TPM Matrix
        tpm_out_file = os.path.join(outputdir, f"{species}_{safe_name}_TPM_matrix.csv")
        os.makedirs(os.path.dirname(tpm_out_file), exist_ok=True)
        
        with open(tpm_out_file, 'w') as f:
            f.write("Gene_ID," + ",".join(all_samples) + "\n")
            for gene in all_genes:
                row = [gene]
                for sample in all_samples:
                    val = tpm_data[gene].get(sample, 0.0)
                    row.append(str(val))
                f.write(",".join(row) + "\n")
                
        print(f"  Saved TPM matrix to {tpm_out_file}")
        
        # Run prepDE.py to get Raw Counts
        count_out_file = os.path.join(outputdir, f"{species}_{safe_name}_Count_matrix.csv")
        transcript_out_file = os.path.join(outputdir, f"{species}_{safe_name}_Transcript_Count_matrix.csv")
        
        cmd = [
            "python", "workflow/scripts/3_5_prepDE.py",
            "-i", prepde_input,
            "-g", count_out_file,
            "-t", transcript_out_file
        ]
        
        try:
            print(f"  Running prepDE.py for {species} - {anno_subpath} to generate raw counts...")
            import subprocess
            subprocess.run(cmd, check=True)
            print(f"  Saved Count matrix to {count_out_file}")
        except Exception as e:
            print(f"  Failed to run prepDE.py: {e}")

def main():
    args = parse_args()
    
    if not args.gtf_base:
        args.gtf_base = "/home/public_software_annotation"

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_outputdir = args.outputdir.rstrip(os.sep)
    final_outputdir = f"{base_outputdir}_{timestamp}"
    print(f"Output directory updated with timestamp: {final_outputdir}")
    
    passed_samples = load_qc_pass_samples(args.filter_csv)
    
    merge_expression_matrices(args.inputdir, final_outputdir, args.gtf_base, passed_samples)
    
    parent_dir = os.path.dirname(args.outputdir)
    flag_file = os.path.join(parent_dir, "Merge_finished.txt")
    
    with open(flag_file, "w") as f:
        f.write(f"Merge finished at {timestamp}. Output: {final_outputdir}\n")
    print(f"Created finish flag at {flag_file}")

if __name__ == "__main__":
    main()
