RNAseq_Pipline：基于 Snakemake 的 RNA-seq 自动化分析流程
=======================================================================

## 简介

**RNAseq_Pipline** 是一套基于 **Snakemake** 构建的模块化、可配置的 RNA-seq 数据分析流程。该流程支持多物种（拟南芥 TAIR10、人类 GRCh38、小鼠 mm 等）、多实验项目的并行处理，只需修改配置文件即可灵活切换分析模式。

### 主要功能

- 自动化比对：基于 HISAT2 完成双端测序数据比对
- BAM 处理：SAMtools 排序、Picard 去重、质量控制
- 转录本定量：使用 StringTie 进行转录本组装和定量分析
- 表达矩阵生成：自动生成基因和转录本的 Count 和 TPM 矩阵
- 差异表达分析：基于 DESeq2 进行差异基因表达分析
- 可视化：生成火山图、PCA 图、样本距离热图等
- 多项目支持：单个配置文件中可同时定义多个物种/实验项目，独立控制各模块开关

---

## 分析流程

```
原始 FASTQ 数据
      │
      ▼
  [可选] 数据下载 (1_download)
   SRA-tools 下载
      │
      ▼
  [可选] 质控 (2_QC)
   FastQC / fastp
      │
      ▼
  Step1: 比对 & 过滤 (3_Align_Filter)
   HISAT2 → SAMtools → Picard 去重
      │
      ▼
  Step2: 转录本定量 (3_Align_Filter)
   StringTie 组装与定量
      │
      ▼
  Step3: 表达矩阵合并 (3_Align_Filter)
   prepDE.py 生成 Count/TPM 矩阵
      │
      ▼
  Step4: 差异表达分析 (4_DEseq)
   DESeq2 差异分析 + 可视化
   火山图、PCA、热图等
```

---

## 目录结构

```
RNAseq_Pipline/
├── Snakefile                     # 流程入口，读取 config 并汇总所有目标文件
├── README.md                     # 项目说明文档
├── .gitignore                    # Git 忽略规则（大文件、结果目录等）
├── environment.yml               # Conda 环境配置文件
│
├── config/                       # 配置中心（用户修改区域）
│   ├── config.yaml               # 全局参数与各项目/模块开关
│   └── RNAseq_metadata.csv       # 样本分组信息（实验组 vs 对照组）
│
├── workflow/                     # 核心流程逻辑
│   ├── anno/                     # 基因组注释文件（本地准备，不上传）
│   │   ├── TAIR10/               # 拟南芥 TAIR10 HISAT2 索引和注释
│   │   │   ├── TAIR10.*.ht2      # HISAT2 索引文件
│   │   │   ├── TAIR10.fa         # 基因组序列
│   │   │   ├── TAIR10.gtf        # 基因注释
│   │   │   ├── TAIR10.mRNA.gtf   # mRNA 注释
│   │   │   └── README.md         # 注释文件说明
│   │   ├── homo/                 # 人类基因组注释
│   │   └── mm/                   # 小鼠基因组注释
│   │
│   ├── resources/                # 原始测序数据（本地准备，不上传）
│   │   ├── TAIR/                 # 拟南芥实验数据
│   │   │   └── fastqfile/        # FASTQ 文件目录
│   │   ├── homo/                 # 人类实验数据
│   │   └── mm/                   # 小鼠实验数据
│   │
│   ├── scripts/                  # 分析脚本
│   │   ├── 1_download.py         # SRA 数据下载
│   │   ├── 2_QC.pl               # 质控脚本
│   │   ├── 3_1_Align.pl          # HISAT2 比对
│   │   ├── 3_2_Filter.pl         # BAM 过滤和去重
│   │   ├── 3_3_Quant.pl          # StringTie 定量
│   │   ├── 3_4_Merge.pl          # 表达矩阵合并
│   │   └── 4_DEseq.R             # DESeq2 差异分析
│   │
│   └── profiles/                 # Snakemake 配置文件
│
├── result/                       # 分析结果输出（自动生成，不上传）
│   └── {species}/
│       └── 3_Align_Filter/       # 比对、过滤、定量结果
│           ├── {sample}/         # 各样本结果
│           │   ├── hisat2file/   # BAM 文件和质控结果
│           │   └── mRNA/         # StringTie 定量结果
│           └── Matrices_*/       # 表达矩阵
│       └── 4_DEseq/              # 差异表达分析结果
│           ├── plots/            # 可视化图表
│           └── tables/           # 差异基因列表
│
├── logs/                         # 各步骤日志（自动生成，不上传）
└── benchmarks/                   # 性能基准测试（自动生成，不上传）
```

---

## 安装与依赖

### 软件依赖

| 软件 | 用途 |
|------|------|
| Snakemake | 流程管理 |
| HISAT2 | 短序列比对 |
| SAMtools | BAM 文件处理 |
| Picard | 去重 |
| StringTie | 转录本组装与定量 |
| prepDE.py | 表达矩阵生成 |
| FastQC / fastp | 质控（可选） |
| **R 4.0+** | DESeq2、ggplot2 等 |
| **Python 3** | pandas、subprocess |
| **Perl** | 各步骤调度脚本 |

### 环境配置（一键安装）

项目根目录提供了完整的 [`environment.yml`](environment.yml)，可一键还原所有 Python/工具依赖：

```bash
conda env create -f environment.yml
conda activate RNAseq_Pipline
```

### R 包安装（版本锁定）

R 包需在 conda 环境激活后单独安装。[`workflow/env/install_R_packages.R`](workflow/env/install_R_packages.R) 锁定了精确版本（DESeq2 1.44.0、ggplot2 3.5.1 等），确保跨服务器结果一致：

```bash
conda activate RNAseq_Pipline
Rscript workflow/env/install_R_packages.R
```

> **一键部署新服务器**：也可直接运行 [`workflow/env/setup_new_server.sh`](workflow/env/setup_new_server.sh)，自动完成 conda 环境创建 + R 包安装 + 工具验证：
> ```bash
> bash workflow/env/setup_new_server.sh
> ```

### Picard 安装

Picard 未包含在 conda 环境中，需单独下载 JAR 文件：

```bash
cd workflow/env
wget https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar
# 验证安装
java -jar picard.jar --version
```

下载后路径默认为 `workflow/env/picard.jar`，与 [`config/config.yaml`](config/config.yaml) 中的默认配置一致。如需使用其他路径，修改配置文件中的 `picard_jar` 字段即可。详见 [`workflow/env/README.md`](workflow/env/README.md)。

---

## 快速开始

### 1. 克隆仓库

```bash
git clone https://github.com/nculiujy/RNAseq_Pipline.git
cd RNAseq_Pipline
```

### 2. 准备参考基因组索引

将 HISAT2 索引文件和基因组注释文件放置在 `workflow/anno/` 对应子目录下：

```
workflow/anno/TAIR10/TAIR10.*.ht2
workflow/anno/TAIR10/TAIR10.fa
workflow/anno/TAIR10/TAIR10.gtf
workflow/anno/TAIR10/TAIR10.mRNA.gtf
```

### 3. 准备原始数据

将 `.fastq.gz` 或 `.fq.gz` 双端测序文件放置在 `workflow/resources/` 对应物种目录下，命名格式为：

```
{样本名}_1.clean.fq.gz
{样本名}_2.clean.fq.gz
```

### 4. 修改配置文件

**[`config/config.yaml`](config/config.yaml)** — 全局参数与项目列表：

```yaml
picard_dir: "/path/to/picard.jar"   # Picard 安装路径
threads: 8                           # 并行线程数

projects:
  - species: "TAIR"
    experiment: "MyExperiment"
    rawdata_dir: "workflow/resources/TAIR/fastqfile"
    index_dir: "workflow/anno/TAIR10/TAIR10"
    gtf_file: "workflow/anno/TAIR10/TAIR10.mRNA.gtf"
    modules:
      1_download: false
      2_QC: false
      3_Align_Filter: true
      4_DEseq: true
```

**[`config/RNAseq_metadata.csv`](config/RNAseq_metadata.csv)** — 样本分组信息：

```csv
sample,group
sample1,control
sample2,control
sample3,treatment
sample4,treatment
```

> `sample` 列填写样本文件名前缀（不含 `_1.clean.fq.gz` 后缀），`group` 列填写分组信息。

### 5. 运行流程

```bash
# 试运行（查看将执行哪些规则，不实际运行）
snakemake -n

# 正式运行（使用 8 个核心）
snakemake -c 8 --rerun-incomplete

# 强制重新运行所有步骤
snakemake -c 8 --forceall
```

---

## 配置说明

### 模块开关

每个 project 可以独立控制各分析模块：

| 模块键 | 功能 | 说明 |
|--------|------|------|
| `1_download` | 数据下载 | 从 SRA 下载原始数据 |
| `2_QC` | 质量控制 | FastQC / fastp |
| `3_Align_Filter` | 比对与定量 | HISAT2 比对 + Picard 去重 + StringTie 定量 |
| `4_DEseq` | 差异分析 | DESeq2 差异表达分析与可视化 |

设为 `true` 开启，`false` 关闭。

### 支持物种

| `species` 配置值 | 对应生物 | 基因组版本 |
|-----------------|---------|-----------|
| `TAIR` | 拟南芥 | TAIR10 |
| `homo` | 人类 | GRCh38 |
| `mm` | 小鼠 | GRCm39 |

---

## 输出结果

### 3_Align_Filter 目录

| 文件类型 | 说明 |
|---------|------|
| `*.dedup.bam` | 去重后的 BAM 文件 |
| `*.dedup.bam.bai` | BAM 索引文件 |
| `*.metrics` | Picard 去重统计 |
| `QC_results.log` | 质控统计结果 |
| `transcripts.gtf` | StringTie 组装的转录本 |
| `gene_abund.tab` | 基因表达丰度 |
| `*_Count_matrix.csv` | 基因 Count 矩阵 |
| `*_TPM_matrix.csv` | 基因 TPM 矩阵 |

### 4_DEseq 目录

| 子目录/文件 | 说明 |
|-----------|------|
| `*_DEG_result.csv` | 差异表达基因完整结果表 |
| `plots/mRNA_volcano.pdf` | 火山图 |
| `plots/sample_PCA.pdf` | PCA 主成分分析图 |
| `plots/sample_dist_heatmap.pdf` | 样本距离热图 |
| `plots/sample_boxplot.pdf` | 样本表达量箱线图 |
| `tables/UP_genes_name.csv` | 上调基因列表 |
| `tables/DOWN_genes_name.csv` | 下调基因列表 |

---

## 注意事项

1. **大文件不上传**：原始 FASTQ 数据、HISAT2 索引（`.ht2`）、基因组 FASTA 和 GTF 文件体积过大，已在 `.gitignore` 中排除，需本地自行准备。
2. **Picard 路径**：请在 [`config.yaml`](config/config.yaml) 中修改 `picard_dir` 为本机实际安装路径。
3. **路径中文字符**：如果 `rawdata_dir` 包含中文字符，请确保系统 locale 支持 UTF-8（`export LANG=zh_CN.UTF-8`）。
4. **多项目并行**：在 `projects` 列表中添加多个条目，Snakemake 会自动并行处理，互不干扰。
5. **样本命名**：样本名称应避免使用特殊字符，建议使用字母、数字、下划线和连字符。

---

## 常见问题

### Q1: 如何添加新的物种？

1. 在 `workflow/anno/` 下创建新物种目录
2. 准备 HISAT2 索引、基因组序列和 GTF 注释文件
3. 在 `config/config.yaml` 中添加新的 project 配置
4. 在 `workflow/resources/` 下创建对应的数据目录

### Q2: 如何只运行部分样本？

修改 [`config/RNAseq_metadata.csv`](config/RNAseq_metadata.csv)，只保留需要分析的样本行。

### Q3: 如何调整差异分析的阈值？

在 [`workflow/scripts/4_DEseq.R`](workflow/scripts/4_DEseq.R) 中修改 `padj` 和 `log2FoldChange` 的阈值。

---

## 引用软件

- [Snakemake](https://snakemake.readthedocs.io/)
- [HISAT2](http://daehwankimlab.github.io/hisat2/)
- [SAMtools](http://www.htslib.org/)
- [Picard](https://broadinstitute.github.io/picard/)
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [DESeq2](https://bioconductor.org/packages/DESeq2/)

---

## 许可证

本项目采用 MIT 许可证。详见 LICENSE 文件。

---

## 联系方式

如有问题或建议，请通过 GitHub Issues 联系。
