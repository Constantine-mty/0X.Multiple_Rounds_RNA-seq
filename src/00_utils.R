# Shared helpers for the configuration-driven multi-round RNA-seq workflow.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    abort("缺少 R 包 '", package, "'。请先安装它后再运行此步骤。")
  }
}

load_configured_ggplot2 <- function(config = list()) {
  runtime <- config$runtime %||% list()
  library_path <- runtime$ggplot2_lib_loc %||% NULL
  if (is.null(library_path) || !length(library_path)) {
    require_package("ggplot2")
    return(invisible(getNamespaceInfo(asNamespace("ggplot2"), "path")))
  }
  library_path <- as.character(library_path)
  if (length(library_path) != 1L || is.na(library_path) || !nzchar(library_path) || !dir.exists(library_path)) {
    abort("runtime.ggplot2_lib_loc 不存在或不是目录：", paste(library_path, collapse = ", "))
  }
  requested_path <- normalizePath(file.path(library_path, "ggplot2"), mustWork = FALSE)
  if (!dir.exists(requested_path)) {
    abort("runtime.ggplot2_lib_loc 中找不到 ggplot2：", requested_path)
  }
  if ("ggplot2" %in% loadedNamespaces()) {
    loaded_path <- normalizePath(getNamespaceInfo(asNamespace("ggplot2"), "path"), mustWork = FALSE)
    if (!identical(loaded_path, requested_path)) {
      abort(
        "ggplot2 已从另一位置加载：", loaded_path,
        "；要求使用：", requested_path,
        "。请在新的 R 会话中运行，避免两个版本同时存在。"
      )
    }
    return(invisible(loaded_path))
  }
  tryCatch(
    loadNamespace("ggplot2", lib.loc = library_path),
    error = function(error) abort("无法从 runtime.ggplot2_lib_loc 加载 ggplot2：", conditionMessage(error))
  )
  loaded_path <- normalizePath(getNamespaceInfo(asNamespace("ggplot2"), "path"), mustWork = FALSE)
  if (!identical(loaded_path, requested_path)) {
    abort("ggplot2 实际加载位置与 runtime.ggplot2_lib_loc 不一致：", loaded_path)
  }
  message("ggplot2 已从指定库加载：", loaded_path)
  invisible(loaded_path)
}

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    abort("无法创建目录：", path)
  }
  invisible(path)
}

resolve_path <- function(project_root, path) {
  if (is.null(path) || !nzchar(path)) abort("配置中的路径不能为空。")
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(project_root, path)
}

safe_path_component <- function(x) {
  x <- trimws(as.character(x))
  # Keep non-ASCII labels (for example Chinese group names), but never allow a path separator.
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[[:cntrl:]]+", "_", x)
  x <- gsub("(^_+|_+$)", "", x)
  if (!nzchar(x)) abort("组名或比较名不能只包含特殊字符。")
  x
}

read_tsv <- function(path, required_columns = character()) {
  if (!file.exists(path)) abort("找不到输入文件：", path)
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext == "csv") "," else "\t"
  out <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  missing <- setdiff(required_columns, colnames(out))
  if (length(missing)) {
    abort(basename(path), " 缺少必要列：", paste(missing, collapse = ", "))
  }
  out
}

write_tsv <- function(x, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = row.names,
    col.names = TRUE, na = ""
  )
  invisible(path)
}

write_rds <- function(x, path) {
  ensure_dir(dirname(path))
  saveRDS(x, path)
  invisible(path)
}

as_logical_flag <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

read_counts_matrix <- function(path, gene_column) {
  counts_df <- read_tsv(path, required_columns = gene_column)
  if (ncol(counts_df) < 2L) abort("count 矩阵至少需要一列基因 ID 和一列样本。")
  genes <- as.character(counts_df[[gene_column]])
  if (anyNA(genes) || any(!nzchar(genes))) abort("count 矩阵中存在空基因 ID。")
  if (anyDuplicated(genes)) abort("count 矩阵中存在重复基因 ID。")
  value_columns <- setdiff(colnames(counts_df), gene_column)
  if (anyDuplicated(value_columns)) abort("count 矩阵中存在重复样本列名。")
  values <- as.matrix(counts_df[, value_columns, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  if (any(!is.finite(values))) abort("count 矩阵含有非数值或缺失值。")
  if (any(values < 0)) abort("count 矩阵不能含负值。")
  if (any(abs(values - round(values)) > 1e-8)) {
    abort("DESeq2 只接受原始整数 count；检测到非整数，请不要输入 TPM/FPKM/CPM 或变换后的矩阵。")
  }
  rownames(values) <- genes
  values
}

read_metadata <- function(path, sample_column, group_column, cycle_column) {
  meta <- read_tsv(path, c(sample_column, group_column, cycle_column))
  for (column in c(sample_column, group_column, cycle_column)) {
    meta[[column]] <- trimws(as.character(meta[[column]]))
    if (anyNA(meta[[column]]) || any(!nzchar(meta[[column]]))) {
      abort("metadata 的 ", column, " 列含有空值。")
    }
  }
  if (anyDuplicated(meta[[sample_column]])) {
    duplicated_samples <- unique(meta[[sample_column]][duplicated(meta[[sample_column]])])
    abort("metadata 中存在重复样本：", paste(duplicated_samples, collapse = ", "))
  }
  meta
}

align_counts_and_metadata <- function(counts, metadata, sample_column) {
  missing_samples <- setdiff(metadata[[sample_column]], colnames(counts))
  if (length(missing_samples)) {
    abort("metadata 中这些样本不在 count 矩阵列名中：", paste(missing_samples, collapse = ", "))
  }
  unused_counts <- setdiff(colnames(counts), metadata[[sample_column]])
  if (length(unused_counts)) {
    message("提示：count 矩阵中有 ", length(unused_counts), " 个样本未在 metadata 中出现，将不参与分析。")
  }
  counts[, metadata[[sample_column]], drop = FALSE]
}

read_contrasts <- function(path) {
  contrasts <- read_tsv(path, c("group_id", "comparison_id", "numerator", "denominator"))
  if (!"enabled" %in% colnames(contrasts)) contrasts$enabled <- TRUE
  for (column in c("group_id", "comparison_id", "numerator", "denominator")) {
    contrasts[[column]] <- trimws(as.character(contrasts[[column]]))
    if (anyNA(contrasts[[column]]) || any(!nzchar(contrasts[[column]]))) {
      abort("contrasts 的 ", column, " 列含有空值。")
    }
  }
  contrasts$enabled <- vapply(contrasts$enabled, as_logical_flag, logical(1))
  duplicate_key <- paste(contrasts$group_id, contrasts$comparison_id, sep = "\r")
  if (anyDuplicated(duplicate_key)) abort("同一 group_id 内的 comparison_id 必须唯一。")
  contrasts
}

relative_result_path <- function(...) {
  file.path(...)
}

read_manifest <- function(project_root) {
  path <- file.path(project_root, "data", "pipeline", "de_manifest.tsv")
  read_tsv(path, c("group_id", "comparison_id", "status", "result_file"))
}

get_config_value <- function(config, name, default = NULL) {
  config[[name]] %||% default
}
