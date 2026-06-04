#!/bin/bash
# 新服务器环境部署脚本
# 用法: bash workflow/env/setup_new_server.sh
set -e

echo "=== Step 1: 创建 conda 环境 ==="
conda env create -f environment.yml
echo "conda 环境创建完成"

echo ""
echo "=== Step 2: 安装 R 包 ==="
conda run -n RNAseq_Pipline Rscript workflow/env/install_R_packages.R

echo ""
echo "=== Step 3: 验证关键工具 ==="
conda run -n RNAseq_Pipline bash -c "
  hisat2 --version | head -1
  samtools --version | head -1
  stringtie --version
  snakemake --version
  python --version
"

echo ""
echo "=== 部署完成，激活环境: conda activate RNAseq_Pipline ==="
