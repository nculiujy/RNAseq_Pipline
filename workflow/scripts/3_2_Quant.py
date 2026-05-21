
import os
import argparse
import subprocess
import glob
from concurrent.futures import ThreadPoolExecutor

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.2: StringTie Quantification")
    parser.add_argument("--inputdir", required=True, help="Directory containing hisat2 output")
    parser.add_argument("--outputdir", required=True, help="Output directory")
    # 移除 species 参数，改为自动检测
    parser.add_argument("--threads", type=int, default=4, help="Threads per task")
    parser.add_argument("--maxparallel", type=int, default=4, help="Max parallel tasks")
    parser.add_argument("--gtf_base", help="GTF annotation base directory")
    parser.add_argument("--gtf_file", type=str, default="", help="Specific GTF file to use (overrides gtf_base auto-detection)")
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

def get_annotations(species, gtf_base_dir):
    annotations = []
    if species == 'human':
        annotations = [
            {"gff": f"{gtf_base_dir}/genomeanno/gencode.v44.mRNA.annotation.gtf", "dir": "mRNA/genecode/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/EnhancerAtlasv2.0_eRNA.gtf", "dir": "eRNA/EnhancerAtlas/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/Ensemblv110_eRNA.gtf", "dir": "eRNA/Ensembl/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/FANTOM5_eRNA.gtf", "dir": "eRNA/FANTOM5/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/GENCODEv44_lncRNA.gtf", "dir": "lncRNA/GENCODE/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/miRBasev22.1_miRNA.gtf", "dir": "miRNA/miRBase/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/MirGeneDBv2.1_miRNA.gtf", "dir": "miRNA/MirGeneDB/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/human/NONCODEv6_lncRNA.gtf", "dir": "lncRNA/NONCODE/stringtie"}
        ]
    elif species == 'mouse':
        annotations = [
            {"gff": f"{gtf_base_dir}/genomeanno/gencode.vM25.mRNA.annotation.gtf", "dir": "mRNA/genecode/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/mouse/EnhancerAltasv2.0_eRNA.gtf", "dir": "eRNA/EnhancerAtlas/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/mouse/Ensemblv110_eRNA.gtf", "dir": "eRNA/Ensembl/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/mouse/FANTOM5_eRNA.gtf", "dir": "eRNA/FANTOM5/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/mouse/GENCODEvM25_lncRNA.gtf", "dir": "lncRNA/GENCODE/stringtie"},
            {"gff": f"{gtf_base_dir}/Reference/RNAanno/mouse/miRBasev22.1_miRNA.gtf", "dir": "miRNA/miRBase/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/mouse/MirGeneDBv2.1_miRNA.gtf", "dir": "miRNA/MirGeneDB/stringtie"},
            {"gff": f"{gtf_base_dir}/ncRNAanno/mouse/NONCODEv6_lncRNA.gtf", "dir": "lncRNA/NONCODE/stringtie"}
        ]
    elif species == 'TAIR':
        annotations = [
            {"gff": "/home/jyliu/GBM_machine_learn/RNAseq_Pipline/workflow/resources/TAIR_reference/TAIR10.gtf", "dir": "mRNA/TAIR10/stringtie"}
        ]
    return annotations

def run_stringtie(hisat2_path, args, rel_dir, sample_id, species_path):
    if args.gtf_file:
        annotations = [{"gff": args.gtf_file, "dir": "mRNA/TAIR10/stringtie"}] # 默认当作一个通用的目录名，或者你可以自己定义规则
        species = "Custom"
    else:
        # 自动推断物种
        species = get_species_from_path(species_path)
        if not species:
            print(f"[{sample_id}] Warning: Could not determine species from path {species_path}. Skipping.")
            return False
            
        annotations = get_annotations(species, args.gtf_base)
    
    print(f"[{sample_id}] Starting StringTie quantification (Species: {species})...")
    dedup_bam = os.path.join(hisat2_path, f"{sample_id}.dedup.bam")

    if not os.path.exists(dedup_bam):
        return False

    for anno in annotations:
        # outputdir/rel_dir/anno_dir/sample_id
        out_dir = os.path.join(args.outputdir, rel_dir, anno["dir"], sample_id)
        os.makedirs(out_dir, exist_ok=True)
        out_gtf = os.path.join(out_dir, "transcripts.gtf")
        gene_abund = os.path.join(out_dir, "gene_abund.tab")
        
        if not os.path.exists(out_gtf):
            stringtie_cmd = [
                "stringtie", "-p", str(args.threads), "-e", "-B",
                "-G", anno["gff"],
                "-A", gene_abund,
                "-o", out_gtf,
                dedup_bam
            ]
            try:
                subprocess.run(stringtie_cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"[{sample_id}] StringTie failed ({anno['dir']}): {e}")
    return True

def main():
    args = parse_args()
    
    if not args.gtf_base:
        args.gtf_base = "/home/public_software_annotation"
    
    tasks = []
    # 递归查找 hisat2file 目录
    for root, dirs, files in os.walk(args.inputdir):
        if "hisat2file" in dirs:
            rel_dir = os.path.relpath(root, args.inputdir)
            species_path = root # 使用绝对路径判断物种
            hisat2_base = os.path.join(root, "hisat2file")
            for sid in os.listdir(hisat2_base):
                if os.path.isdir(os.path.join(hisat2_base, sid)):
                    tasks.append((os.path.join(hisat2_base, sid), args, rel_dir, sid, species_path))

    with ThreadPoolExecutor(max_workers=args.maxparallel) as executor:
        futures = [executor.submit(run_stringtie, *t) for t in tasks]
        for f in futures:
            f.result()

    with open(os.path.join(args.outputdir, "Quant_finished.txt"), "w") as f:
        f.write("Quantification finished.\n")

if __name__ == "__main__":
    main()
