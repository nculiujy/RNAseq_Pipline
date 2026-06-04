#!/usr/bin/env Rscript
# R包安装脚本 - 锁定版本，用于跨服务器环境复现
# 在新服务器conda环境激活后运行: Rscript workflow/env/install_R_packages.R

options(repos = c(
  CRAN     = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
  CRANfull = "https://cran.r-project.org/"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

# 配置Bioconductor镜像
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# ---- 锁定版本的包列表 ----
cran_pkgs <- list(
  optparse     = "1.7.5",
  yaml         = "2.3.10",
  ggplot2      = "3.5.1",
  dplyr        = "1.1.4",
  RColorBrewer = "1.1-3",
  pheatmap     = "1.0.12",
  ggrepel      = "0.9.5"
)
bioc_pkgs <- list(
  DESeq2 = "1.44.0"
)

install_exact <- function(pkg, ver, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    installed_ver <- as.character(packageVersion(pkg))
    if (installed_ver == ver) {
      message(pkg, " ", ver, " already installed, skipping.")
      return(invisible(NULL))
    }
    message(pkg, " version mismatch: installed=", installed_ver, " required=", ver)
  }
  if (bioc) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    # 尝试从CRAN archive安装精确版本
    url <- paste0("https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/Archive/",
                  pkg, "/", pkg, "_", ver, ".tar.gz")
    tryCatch(
      install.packages(url, repos = NULL, type = "source"),
      error = function(e) install.packages(pkg)  # fallback到最新版
    )
  }
}

for (pkg in names(cran_pkgs)) install_exact(pkg, cran_pkgs[[pkg]])
for (pkg in names(bioc_pkgs)) install_exact(pkg, bioc_pkgs[[pkg]], bioc = TRUE)

# 验证
cat("\n=== 安装结果验证 ===\n")
all_pkgs <- c(names(cran_pkgs), names(bioc_pkgs))
for (pkg in all_pkgs) {
  ver <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "NOT INSTALLED")
  cat(sprintf("%-15s %s\n", pkg, ver))
}
