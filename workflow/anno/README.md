# HISAT2 注释文件 (Index) 下载与配置指南

本目录 (`workflow/anno`) 用于统一存放 RNA-seq 流程中 HISAT2 比对步骤所需的参考基因组索引文件和基因注释文件。

根据流程的设计，目前主要支持以下三类核心物种。本教程将指导您如何获取并配置这三个物种的注释文件。

---

## 1. 拟南芥 (TAIR / TAIR10)

拟南芥的参考基因组和注释文件支持两种来源：**TAIR 官方 GFF3**（推荐，注释最权威）和 **Ensembl Plants GTF**。两种方式最终都需要构建 HISAT2 索引。

---

### 方式 1：使用 TAIR 官方 GFF3 文件（推荐）

TAIR 官方提供的 GFF3 文件注释最为权威，但需要先将其转换为 GTF 格式才能用于 HISAT2 和 StringTie。

```bash
# 1. 进入注释文件存放目录
cd workflow/anno/TAIR10

# 2. 下载 TAIR 官方基因组 FASTA 文件
wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas -O TAIR10.fa

# 3. 下载 TAIR 官方 GFF3 注释文件
wget "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff" -O TAIR10_genes.gff

# 4. 将 GFF3 转换为 GTF 格式（需要 gffread，已包含在 conda 环境中）
gffread TAIR10_genes.gff -T -o TAIR10.gtf

# 5. 提取 mRNA 类型的注释进行定量（StringTie 使用）
grep -E '^#|gene_id' TAIR10.gtf | awk '$3 == "exon" || $3 == "transcript" || $3 == "gene"' > TAIR10.mRNA.gtf

# 6. 使用 hisat2-build 构建索引（需要几分钟）
hisat2-build TAIR10.fa TAIR10
```

> **注意**：TAIR 官方 GFF3 中染色体名称为 `Chr1`、`Chr2` 等，而 Ensembl 版本为 `1`、`2` 等。请确保 FASTA 和 GFF3 的染色体命名一致。

---

### 方式 2：使用 Ensembl Plants GTF 文件

```bash
# 1. 进入注释文件存放目录
cd workflow/anno/TAIR10

# 2. 从 Ensembl Plants 下载基因组 FASTA
wget http://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# 3. 下载 GTF 注释文件
wget http://ftp.ensemblgenomes.org/pub/plants/release-58/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz

# 4. 解压并重命名
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.58.gtf.gz
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10.fa
mv Arabidopsis_thaliana.TAIR10.58.gtf TAIR10.gtf

# 5. 提取 mRNA 类型的注释进行定量
awk '$3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "UTR"' TAIR10.gtf | grep -E 'gene_biotype "protein_coding"|transcript_biotype "protein_coding"' > TAIR10.mRNA.gtf

# 6. 使用 hisat2-build 构建索引
hisat2-build TAIR10.fa TAIR10
```

---

构建完成后，该目录下会生成以 `TAIR10` 为前缀的 `.ht2` 索引文件（共 8 个）。

### 目录文件清单

- `TAIR10.fa` - 未压缩的完整参考基因组序列文件
- `TAIR10.gtf` - 完整的基因注释文件
- `TAIR10.mRNA.gtf` - 仅包含 mRNA (protein_coding) 类型的基因注释文件，专门用于 StringTie 的基因/转录本定量
- `TAIR10.*.ht2` - 用于 HISAT2 进行序列比对的索引文件（共 8 个文件）

### 配置文件设置

在 `config/config.yaml` 中配置：

```yaml
projects:
  - species: "TAIR"
    experiment: "MyExperiment"
    rawdata_dir: "workflow/resources/TAIR/fastqfile"
    index_dir: "workflow/anno/TAIR10/TAIR10"
    gtf_file: "workflow/anno/TAIR10/TAIR10.mRNA.gtf"
    modules:
      3_Align_Filter: true
      4_DEseq: true
```

---

## 2. 人类 (homo / GRCh38)

对于人类基因组，我们推荐使用 Illumina 或 Johns Hopkins University 提供的官方预编译 HISAT2 索引，这样可以节省大量的本地编译时间。

### 下载与解压步骤

在终端中依次运行以下命令：

```bash
# 1. 进入注释文件存放目录
cd workflow/anno

# 2. 创建并进入人类基因组目录
mkdir -p homo
cd homo

# 3. 下载官方提供的 GRCh38 (hg38) 预编译 HISAT2 索引包
# 注意：这个文件较大（约 4.2 GB），下载可能需要一些时间
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# 4. 解压文件
tar -xzf grch38_genome.tar.gz

# 5. 下载 GTF 注释文件
wget http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
mv Homo_sapiens.GRCh38.110.gtf GRCh38.gtf

# 6. 提取 mRNA 注释（可选，根据需要）
awk '$3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "UTR"' GRCh38.gtf | grep -E 'gene_biotype "protein_coding"|transcript_biotype "protein_coding"' > GRCh38.mRNA.gtf

# 7. 删除压缩包以节省空间
rm grch38_genome.tar.gz
```

### 配置文件设置

在 `config/config.yaml` 中配置：

```yaml
projects:
  - species: "homo"
    experiment: "MyExperiment"
    rawdata_dir: "workflow/resources/homo/fastqfile"
    index_dir: "workflow/anno/homo/grch38/genome"
    gtf_file: "workflow/anno/homo/GRCh38.mRNA.gtf"
    modules:
      3_Align_Filter: true
      4_DEseq: true
```

---

## 3. 小鼠 (mm / GRCm39)

与人类基因组类似，小鼠基因组也提供了官方预编译的 HISAT2 索引文件，可直接下载使用。

### 下载与解压步骤

在终端中依次运行以下命令：

```bash
# 1. 进入注释文件存放目录
cd workflow/anno

# 2. 创建并进入小鼠基因组目录
mkdir -p mm
cd mm

# 3. 下载官方提供的 GRCm39 (mm39) 预编译 HISAT2 索引包
# 注意：这个文件较大（约 3.8 GB），下载可能需要一些时间
wget https://genome-idx.s3.amazonaws.com/hisat/grcm39_genome.tar.gz

# 4. 解压文件
tar -xzf grcm39_genome.tar.gz

# 5. 下载 GTF 注释文件
wget http://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz
gunzip Mus_musculus.GRCm39.110.gtf.gz
mv Mus_musculus.GRCm39.110.gtf GRCm39.gtf

# 6. 提取 mRNA 注释（可选，根据需要）
awk '$3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "UTR"' GRCm39.gtf | grep -E 'gene_biotype "protein_coding"|transcript_biotype "protein_coding"' > GRCm39.mRNA.gtf

# 7. 删除压缩包以节省空间
rm grcm39_genome.tar.gz
```

### 配置文件设置

在 `config/config.yaml` 中配置：

```yaml
projects:
  - species: "mm"
    experiment: "MyExperiment"
    rawdata_dir: "workflow/resources/mm/fastqfile"
    index_dir: "workflow/anno/mm/grcm39/genome"
    gtf_file: "workflow/anno/mm/GRCm39.mRNA.gtf"
    modules:
      3_Align_Filter: true
      4_DEseq: true
```

---

## 常见问题

### Q1: 如何验证索引文件是否正确构建？

对于 HISAT2 索引，应该有 8 个 `.ht2` 文件（编号从 1 到 8）：

```bash
ls -lh workflow/anno/TAIR10/TAIR10.*.ht2
# 应该看到：
# TAIR10.1.ht2
# TAIR10.2.ht2
# ...
# TAIR10.8.ht2
```

### Q2: 下载速度太慢怎么办？

可以使用国内镜像源或使用 `aria2c` 等多线程下载工具：

```bash
# 使用 aria2c 加速下载
aria2c -x 16 -s 16 https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
```

### Q3: 如何选择合适的 GTF 注释版本？

- 建议使用与参考基因组版本匹配的 GTF 注释文件
- 对于 RNA-seq 定量分析，通常使用 `mRNA.gtf`（仅包含 protein_coding 基因）
- 如果需要分析 lncRNA 等非编码 RNA，使用完整的 GTF 文件

### Q4: 磁盘空间不足怎么办？

- 构建索引后可以删除原始的 FASTA 文件（但建议保留备份）
- 解压后立即删除压缩包
- 只下载和构建当前需要分析的物种索引

---

## 总结

无论您分析的是哪种物种，都需要在 `config/config.yaml` 中正确设置以下参数：

1. **`species`** - 决定了下游分析所调用的物种特异性参数
2. **`index_dir`** - HISAT2 索引文件的**路径前缀**（注意：**不要**包含 `.1.ht2` 等后缀）
3. **`gtf_file`** - 基因注释文件的完整路径，用于 StringTie 定量和 DESeq2 分析

配置示例：

```yaml
# 拟南芥
index_dir: "workflow/anno/TAIR10/TAIR10"
gtf_file: "workflow/anno/TAIR10/TAIR10.mRNA.gtf"

# 人类
index_dir: "workflow/anno/homo/grch38/genome"
gtf_file: "workflow/anno/homo/GRCh38.mRNA.gtf"

# 小鼠
index_dir: "workflow/anno/mm/grcm39/genome"
gtf_file: "workflow/anno/mm/GRCm39.mRNA.gtf"
```

---

## 参考资源

- [HISAT2 官方网站](http://daehwankimlab.github.io/hisat2/)
- [Ensembl Plants](http://plants.ensembl.org/)
- [Ensembl](http://www.ensembl.org/)
- [HISAT2 预编译索引下载](https://daehwankimlab.github.io/hisat2/download/)
