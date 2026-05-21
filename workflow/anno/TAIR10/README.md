# TAIR10 参考基因组与注释文件

本目录包含了用于拟南芥 (*Arabidopsis thaliana*，TAIR10版本) RNA-seq 分析的参考基因组、HISAT2 索引以及基因注释文件。

## 数据来源

本目录下的文件主要来源于 Ensembl Plants (Release 58)。

- **参考基因组 (FASTA):** `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz`
- **基因注释 (GTF):** `Arabidopsis_thaliana.TAIR10.58.gtf.gz`

## 如何重新生成这些文件

如果您需要从头开始重新生成或下载这些文件，可以参考以下步骤：

### 1. 下载参考基因组和 GTF 注释文件
```bash
wget http://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget http://ftp.ensemblgenomes.org/pub/plants/release-58/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz
```

### 2. 解压并重命名文件
为了方便后续流程调用，我们将解压后的文件重命名：
```bash
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.58.gtf.gz
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10.fa
mv Arabidopsis_thaliana.TAIR10.58.gtf TAIR10.gtf
```

### 3. 构建 HISAT2 索引
确保您的环境中已安装了 `hisat2` 工具，然后运行以下命令构建比对索引：
```bash
hisat2-build TAIR10.fa TAIR10
```
这将在当前目录下生成一系列名为 `TAIR10.1.ht2`, `TAIR10.2.ht2` ... `TAIR10.8.ht2` 的索引文件。

### 4. 提取 mRNA 类型的注释进行定量
为了只针对 mRNA 进行表达量定量（过滤掉 ncRNA 等），需要从原始的 GTF 文件中单独提取包含蛋白编码基因或转录本的信息：
```bash
awk '$3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "UTR"' TAIR10.gtf | grep -E 'gene_biotype "protein_coding"|transcript_biotype "protein_coding"' > TAIR10.mRNA.gtf
```
*(注意：不同的注释文件可能在 `biotype` 字段的命名上有所不同，提取时需根据实际 GTF 的属性列调整 grep 关键字。针对当前使用的 GTF，可直接根据实际情况提取)*

## 目录文件清单

- `TAIR10.fa` - 未压缩的完整参考基因组序列文件。
- `TAIR10.gtf` - 完整的基因注释文件。
- `TAIR10.mRNA.gtf` - 仅包含 mRNA (protein_coding) 类型的基因注释文件，专门用于 StringTie 的基因/转录本定量。
- `TAIR10.*.ht2` - 用于 HISAT2 进行序列比对的索引文件。