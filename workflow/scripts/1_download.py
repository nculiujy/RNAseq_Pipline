
import os
import sys
import subprocess
import concurrent.futures
import shutil
import glob
import time

def get_srr_jobs(rawdata_dir):
    print(f"[Info] 扫描目录: {rawdata_dir}")
    jobs = []
    folders = []
    # 递归扫描所有子目录寻找 SRR.txt
    for root, dirs, files in os.walk(rawdata_dir):
        if 'SRR.txt' in files:
            # 获取相对于 rawdata_dir 的路径，例如 "homo/rawdata/GSE263588"
            folder = os.path.relpath(root, rawdata_dir)
            srr_txt = os.path.join(root, 'SRR.txt')
            print(f"[Info] 检测到SRR子目录: {folder}，包含SRR.txt")
            folders.append(folder)
            with open(srr_txt) as f:
                srr_list = [line.strip() for line in f if line.strip()]
            print(f"[Info] {folder} 中需下载SRA编号: {srr_list}")
            jobs += [(folder, srr) for srr in srr_list]
    print(f"[Info] 需下载的总Jobs: {len(jobs)}，涉及文件夹: {folders}")
    return jobs, folders

def download_sra(result_dir, folder, srr):
    # 下载到 result_dir/folder 下
    target_dir = os.path.join(result_dir, folder)
    os.makedirs(target_dir, exist_ok=True)
    
    sra_path = os.path.join(target_dir, f"{srr}.sra")
    if os.path.exists(sra_path):
        print(f"[Skip] {folder}: {srr}.sra 已存在，跳过下载")
        return True
    print(f"[Start] {folder}: 正在下载 {srr} ...")
    cmd = [
        "prefetch",
        "--max-size", "30GB",
        "-O", target_dir,
        srr
    ]
    print(f"[Debug] 执行命令: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        print(result.stdout)
        print(f"[Success] {folder}: {srr} 下载完成")
        
        # prefetch 可能会创建 srr 文件夹，也可能直接下载文件，需要处理
        srr_dir = os.path.join(target_dir, srr)
        srr_file_in_dir = os.path.join(srr_dir, f"{srr}.sra")
        
        if os.path.exists(srr_file_in_dir):
             shutil.move(srr_file_in_dir, sra_path)
             try:
                os.rmdir(srr_dir)
             except:
                pass

        return result.returncode == 0
    except Exception as e:
        print(f"[Error] {folder}: {srr} 下载失败: {e}")
        return False

def get_failed_jobs(result_dir, jobs):
    """返回未下载成功的(SRR, folder)列表"""
    failed = []
    for folder, srr in jobs:
        sra_path = os.path.join(result_dir, folder, f"{srr}.sra")
        if not os.path.exists(sra_path):
            failed.append((folder, srr))
    return failed

def retry_failed_jobs(result_dir, failed_jobs, max_workers):
    if not failed_jobs:
        print("[Info] 没有需要重试的任务")
        return []
    print(f"[Info] 准备重试 {len(failed_jobs)} 个失败的SRR")
    def job_func(args):
        folder, srr = args
        print(f"[Retry] 再次尝试下载: {folder} - {srr}")
        return download_sra(result_dir, folder, srr)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(job_func, failed_jobs))
    
    still_failed = []
    for i, (folder, srr) in enumerate(failed_jobs):
        sra_path = os.path.join(result_dir, folder, f"{srr}.sra")
        if not os.path.exists(sra_path):
            still_failed.append((folder, srr))
    return still_failed

def decompress_all_sra(result_dir, folders, max_workers=8):
    print("[Step] 开始解压所有SRA文件")
    sra_files = []
    for folder in folders:
        folder_path = os.path.join(result_dir, folder)
        sra_files += glob.glob(os.path.join(folder_path, "*.sra"))
    print(f"[Info] 共需要解压 {len(sra_files)} 个SRA文件")

    def decompress_one(sra_file):
        folder_path = os.path.dirname(sra_file)
        cmd = [
            "fasterq-dump",
            "--split-files",
            "--threads", "4",
            "-O", folder_path,
            sra_file
        ]
        print(f"[Decompress] {os.path.basename(sra_file)} -> {folder_path}")
        try:
            subprocess.run(cmd, check=True)
            print(f"[Success] {os.path.basename(sra_file)} 解压完成")
            try:
                os.remove(sra_file)
                print(f"[Cleanup] 删除原SRA文件: {sra_file}")
            except Exception as del_e:
                print(f"[Warning] 删除 {sra_file} 失败: {del_e}")
            return True
        except Exception as e:
            print(f"[Error] {os.path.basename(sra_file)} 解压失败: {e}")
            return False

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(decompress_one, sra_files))

    failed = [f for f, r in zip(sra_files, results) if not r]
    return failed

def main():
    print("="*40)
    print("[Step] 开始执行下载流程")
    if len(sys.argv) != 3:
        print("[Error] 用法: python 1_download.py rawdata_folder_path result_folder_path")
        sys.exit(1)
    
    rawdata_dir = sys.argv[1] # rawdata 根目录
    result_dir = sys.argv[2]  # 下载结果目录
    
    print(f"[Info] 原始数据配置目录: {rawdata_dir}")
    print(f"[Info] 结果输出目录: {result_dir}")
    
    if not os.path.exists(rawdata_dir):
        print(f"[Error] 目录不存在: {rawdata_dir}")
        sys.exit(1)
        
    os.makedirs(result_dir, exist_ok=True)
    
    jobs, folders = get_srr_jobs(rawdata_dir)
    print("[Step] 进入下载阶段")
    max_workers = int(os.environ.get("DOWNLOAD_THREADS", "8"))

    def job_func(args):
        folder, srr = args
        return download_sra(result_dir, folder, srr)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        list(executor.map(job_func, jobs))

    # === 自动重试 ===
    retry_times = 2
    failed_jobs = get_failed_jobs(result_dir, jobs)
    for cycle in range(retry_times):
        if not failed_jobs:
            break
        print(f"[RetryPhase] 第{cycle+1}次重试 {len(failed_jobs)} 个任务")
        still_failed = retry_failed_jobs(result_dir, failed_jobs, max_workers)
        if not still_failed:
            print("[Info] 所有重试任务已成功")
            break
        failed_jobs = still_failed
        time.sleep(3)

    # 下载全部成功才能解压
    if not get_failed_jobs(result_dir, jobs):
        decompress_failed = decompress_all_sra(result_dir, folders, max_workers=8)
        if not decompress_failed:
            finished_path = os.path.join(result_dir, "finished.txt")
            with open(finished_path, "w") as f:
                f.write("done\n")
            print(f"[Done] 工作全部完成，生成标志文件: {finished_path}")
        else:
            print(f"[Error] 解压失败，详情请检查日志")
    else:
        print(f"[Error] 下载未完全成功")
    print("="*40)

if __name__ == "__main__":
    main()
