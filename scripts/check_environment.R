#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2L && args[[1]] == "--config") args[[2]] else "config.yml"
if (!file.exists(config_path)) stop("找不到配置文件：", config_path, call. = FALSE)
if (!requireNamespace("yaml", quietly = TRUE)) stop("缺少 yaml 包。", call. = FALSE)
config <- yaml::read_yaml(config_path)
packages <- c(
  "DESeq2", "ggplot2", "clusterProfiler", "GSVA", "limma", "readxl", "pheatmap", "patchwork",
  "ggrepel", "ggalluvial", "openxlsx"
)
# Avoid requiring the full pipeline just to run this small check.
packages <- unique(packages)
runtime_config <- config$runtime
if (is.null(runtime_config) || !length(runtime_config)) runtime_config <- list()
ggplot2_lib_loc <- runtime_config$ggplot2_lib_loc
if (!is.null(ggplot2_lib_loc) && length(ggplot2_lib_loc) && !dir.exists(ggplot2_lib_loc)) {
  stop("runtime.ggplot2_lib_loc 不存在或不是目录：", ggplot2_lib_loc, call. = FALSE)
}
available <- vapply(packages, function(package) {
  if (identical(package, "ggplot2") && !is.null(ggplot2_lib_loc) && length(ggplot2_lib_loc)) {
    return(requireNamespace(package, lib.loc = ggplot2_lib_loc, quietly = TRUE))
  }
  requireNamespace(package, quietly = TRUE)
}, logical(1))
print(data.frame(package = packages, available = available, row.names = NULL))
if (!is.null(ggplot2_lib_loc) && length(ggplot2_lib_loc)) {
  message("ggplot2 将从指定库检查并加载：", ggplot2_lib_loc)
}
if (!all(available)) {
  message("缺失的包仅会阻止对应步骤：DESeq2/ggplot2 用于差异分析与图；clusterProfiler 用于富集；GSVA/limma 用于 GSVA；readxl/pheatmap/patchwork 用于 GO-BP 分类后 GSVA 绘图；ggrepel 用于分类 GSEA 气泡图；ggalluvial/openxlsx 用于分类后趋势桑基图/Excel 汇总。")
  quit(status = 1L)
}
message("R 环境满足完整 pipeline 的依赖。")
