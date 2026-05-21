
import os
import argparse
import subprocess
import glob
import re
from concurrent.futures import ThreadPoolExecutor

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.1: HISAT2 Alignment")
    parser.add_argument("--inputdir", required=True, help="Input directory containing .clean.fastq files")
    parser.add_argument("--outputdir", required=True, help="Output directory")
    # 移除 species 参数，改为自动检测
    parser.add_argument("--threads", type=int, default=4, help="Threads per task")
    parser.add_argument("--maxparallel", type=int, default=4, help="Max parallel tasks")
    parser.add_argument("--hisat2_index", type=str, default="", help="Path to HISAT2 index")
    parser.add_argument("--picard_jar", type=str, default="workflow/env/picard.jar", help="Path to picard.jar")
    return parser.parse_args()

def get_species_from_path(path):
    """
    根据路径判断物种。
    假设路径包含 "homo", "mouse" 或 "TAIR"。
    """
    path_lower = path.lower()
    if "homo" in path_lower:
        return "human"
    elif "mouse" in path_lower:
        return "mouse"
    elif "tair" in path_lower:
        return "TAIR"
    return None

def get_index_dir(species):
    if species == 'human':
        return "/home/public_software_annotation/genomeanno/humanhisat2Index/GRCh38"
    elif species == 'mouse':
        return "/home/public_software_annotation/genomeanno/mousehisat2Index/GRCm38"
    elif species == 'TAIR':
        return "/home/jyliu/GBM_machine_learn/RNAseq_Pipline/workflow/resources/TAIR_reference/TAIR10"
    return None

def run_hisat2(sample_id, fq1, fq2, args, rel_dir, species_path):
    # 根据物种获取索引
    index_dir = args.hisat2_index
    species = get_species_from_path(species_path)
    
    if not index_dir:
        # 如果没有传入，尝试自动推断（向后兼容）
        if not species:
            print(f"[{sample_id}] Warning: Could not determine species from path {species_path}. Skipping.")
            return False
        index_dir = get_index_dir(species)

    if not index_dir:
        print(f"[{sample_id}] Error: HISAT2 index directory not provided and could not be inferred.")
        return False
        
    species_display = species if species else "Custom"

    print(f"[{sample_id}] Starting HISAT2 alignment (Species: {species_display})...")
    
    # 保持输出目录结构：outputdir/rel_dir/hisat2file/sample_id
    # 注意：这里的 rel_dir 已经是移除了 rawdata 的干净路径
    hisat2_dir = os.path.join(args.outputdir, rel_dir, "hisat2file", sample_id)
    os.makedirs(hisat2_dir, exist_ok=True)
    
    bam_file = os.path.join(hisat2_dir, f"{sample_id}.sorted.bam")
    dedup_bam = os.path.join(hisat2_dir, f"{sample_id}.dedup.bam")
    qc_log = os.path.join(hisat2_dir, "QC_results.log")

    if os.path.exists(dedup_bam):
        print(f"[{sample_id}] Dedup BAM already exists. Skipping.")
        return True

    sam_file = os.path.join(hisat2_dir, "accepted_hits.sam")
    
    # HISAT2 command
    hisat2_cmd = [
        "hisat2", "-x", index_dir, "-p", str(args.threads), "--dta",
        "--rg-id", sample_id, "--rg", f"SM:{sample_id}"
    ]
    if fq2:
        hisat2_cmd.extend(["-1", fq1, "-2", fq2])
    else:
        hisat2_cmd.extend(["-U", fq1])
    
    hisat2_cmd.extend(["-S", sam_file])
    
    try:
        with open(qc_log, "w") as log:
            subprocess.run(hisat2_cmd, check=True, stderr=log)
        
        subprocess.run(f"samtools view -bS {sam_file} | samtools sort -@ {args.threads} -o {bam_file}", shell=True, check=True)
        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
        
        if os.path.exists(sam_file):
            os.remove(sam_file)

        picard_jar = args.picard_jar
        dedup_cmd = [
            "java", "-Xmx15g", "-jar", picard_jar, "MarkDuplicates",
            f"I={bam_file}", f"O={dedup_bam}",
            f"METRICS_FILE={hisat2_dir}/{sample_id}.metrics",
            "REMOVE_DUPLICATES=true", "ASSUME_SORT_ORDER=coordinate"
        ]
        subprocess.run(dedup_cmd, check=True)
        subprocess.run(f"samtools index {dedup_bam}", shell=True, check=True)

        if os.path.exists(bam_file):
            os.remove(bam_file)
            os.remove(bam_file + ".bai")

    except subprocess.CalledProcessError as e:
        print(f"[{sample_id}] Alignment failed: {e}")
        return False
    
    return True

def main():
    args = parse_args()
            
    # 递归搜索所有子目录
    tasks = []
    for root, dirs, files in os.walk(args.inputdir):
        for f in files:
            if f.endswith((".clean.fastq", ".clean.fq", ".clean.fastq.gz", ".clean.fq.gz")):
                # 获取相对路径目录
                rel_dir = os.path.relpath(root, args.inputdir)
                
                # rel_dir 已经是相对路径，例如 homo/GSE... 或 TAIR/cleandata
                # 为了防止 rel_dir 中没有物种名称，直接用 root 来判断
                species_path = root
                
                match = re.search(r"(.+)_([12])\.clean\.(fastq|fq)(\.gz)?$", f)
                if match:
                    sid = match.group(1)
                    read_type = match.group(2)
                    ext = match.group(3)
                    gz = match.group(4) or ""
                    if read_type == "1":
                        # 检查配对
                        fq1 = os.path.join(root, f)
                        fq2 = os.path.join(root, f.replace(f"_1.clean.{ext}{gz}", f"_2.clean.{ext}{gz}"))
                        if not os.path.exists(fq2):
                            fq2 = None
                        tasks.append((sid, fq1, fq2, args, rel_dir, species_path))

    print(f"Found {len(tasks)} samples to process.")
    
    def process_task(task):
        sid, fq1, fq2, args, rel_dir, species_path = task
        return run_hisat2(sid, fq1, fq2, args, rel_dir, species_path)

    with ThreadPoolExecutor(max_workers=args.maxparallel) as executor:
        futures = [executor.submit(process_task, t) for t in tasks]
        for f in futures:
            f.result()

    with open(os.path.join(args.outputdir, "Align_finished.txt"), "w") as f:
        f.write("Alignment finished.\n")

if __name__ == "__main__":
    main()
