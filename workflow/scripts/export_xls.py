import pandas as pd
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Convert CSV matrices to .xls format (tab-separated)")
    parser.add_argument("--count_csv", required=True, help="Input raw count matrix (CSV)")
    parser.add_argument("--tpm_csv", required=True, help="Input TPM matrix (CSV)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()

def main():
    args = parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Convert Count matrix -> readcount.xls
    if os.path.exists(args.count_csv):
        df_count = pd.read_csv(args.count_csv)
        # Assuming the first column is the ID column, rename it to 'id' if needed to match common .xls formats
        # or leave it as 'gene_id'. Let's rename the first column to 'gene_id' just in case.
        col_name = df_count.columns[0]
        if col_name != 'gene_id':
            df_count.rename(columns={col_name: 'gene_id'}, inplace=True)
        
        out_count = os.path.join(args.outdir, "readcount.xls")
        df_count.to_csv(out_count, sep='\t', index=False)
        print(f"Generated {out_count}")
    else:
        print(f"Error: {args.count_csv} not found.")

    # 2. Convert TPM matrix -> fpkm_sample.xls (or tpm_sample.xls, but user requested fpkm_sample.xls)
    # Note: Although the user asked for fpkm_sample.xls, StringTie outputs TPM. We will save the TPM data into the requested filename.
    if os.path.exists(args.tpm_csv):
        df_tpm = pd.read_csv(args.tpm_csv)
        col_name = df_tpm.columns[0]
        if col_name != 'gene_id':
            df_tpm.rename(columns={col_name: 'gene_id'}, inplace=True)
            
        out_tpm = os.path.join(args.outdir, "fpkm_sample.xls")
        df_tpm.to_csv(out_tpm, sep='\t', index=False)
        print(f"Generated {out_tpm} (Note: Contains TPM values as calculated by the pipeline)")
    else:
        print(f"Error: {args.tpm_csv} not found.")

if __name__ == "__main__":
    main()
