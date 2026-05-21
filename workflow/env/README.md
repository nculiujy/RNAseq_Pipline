# 软件环境配置指南

本目录 (`workflow/env`) 用于存放流程所需的外部软件包，主要是 **Picard**（无法通过 conda 直接安装的工具）。

---

## Picard 安装

Picard 是一个用于处理高通量测序数据的 Java 工具集，本流程使用它进行 BAM 文件去重（MarkDuplicates）。

### 下载步骤

```bash
# 进入软件目录
cd workflow/env

# 下载 Picard JAR 文件（推荐使用稳定版本）
wget https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

# 验证下载是否完整
java -jar picard.jar --version
```

下载完成后，`picard.jar` 应位于 `workflow/env/picard.jar`，与 [`config/config.yaml`](../../config/config.yaml) 中的默认路径一致：

```yaml
picard_jar: "workflow/env/picard.jar"
```

### 使用其他版本或路径

如果已在其他位置安装了 Picard，只需修改 [`config/config.yaml`](../../config/config.yaml) 中的路径：

```yaml
picard_jar: "/path/to/your/picard.jar"
```

### 系统要求

- **Java 8 或更高版本**（推荐 Java 11+）
- 验证 Java 版本：`java -version`
- 建议分配至少 **8 GB 内存**（脚本默认使用 `-Xmx15g`）

---

## 注意事项

- `picard.jar` 文件体积约 60 MB，已在 `.gitignore` 中排除，不会上传到 GitHub
- 每次克隆仓库后需重新下载 `picard.jar`
