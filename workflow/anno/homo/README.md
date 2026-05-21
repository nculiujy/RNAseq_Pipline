# 人类 (Homo sapiens) 参考基因组与注释文件

本目录应包含用于人类 (Homo sapiens) RNA-seq 分析的参考基因组、HISAT2 索引以及基因注释文件。

## 如何生成这些文件 (示例流程)

您可以参考以下通用流程从 Ensembl 或 GENCODE 获取并构建必要的文件。

### 1. 下载参考基因组和 GTF 注释文件
*(以 GENCODE v44 为例)*
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
```

### 2. 解压文件
```bash
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v44.annotation.gtf.gz
mv GRCh38.primary_assembly.genome.fa homo.fa
mv gencode.v44.annotation.gtf homo.gtf
```

### 3. 构建 HISAT2 索引
*(构建人类基因组索引需要大量内存和时间)*
```bash
hisat2-build homo.fa homo_index
```
这将生成一系列名为 `homo_index.1.ht2` 到 `homo_index.8.ht2` 的索引文件。

### 4. 提取 mRNA 类型的注释进行定量
为了只针对 mRNA 进行表达量定量，提取 `protein_coding` 相关信息：
```bash
awk '$3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "UTR"' homo.gtf | grep -E 'gene_type "protein_coding"|transcript_type "protein_coding"' > homo.mRNA.gtf
```
*(注意：GENCODE 中的属性名通常为 `gene_type`，Ensembl 通常为 `gene_biotype`，请根据实际下载的 GTF 格式进行调整)*

## 预期目录文件清单

- `homo.fa` - 未压缩的完整参考基因组序列文件。
- `homo.gtf` - 完整的基因注释文件。
- `homo.mRNA.gtf` - 仅包含 mRNA 类型的基因注释文件。
- `homo_index.*.ht2` - HISAT2 索引文件。