# Shared helpers for interactive, classified-GSEA notebooks. No work runs on source().

gsea_nb_read_classified <- function(path) {
  if (!file.exists(path)) stop("找不到分类后 GSEA 文件：", path, call. = FALSE)
  out <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  out <- out[, !grepl("^(Unnamed: ?[0-9]*|\\.\\.\\.[0-9]+|)$", names(out)), drop = FALSE]
  required <- c("Group", "Comparison", "ID", "Description", "NES", "p.adjust", "core_enrichment")
  missing <- setdiff(required, names(out))
  if (length(missing)) stop("分类后 GSEA 缺少列：", paste(missing, collapse = ", "), call. = FALSE)
  if (!"Direction" %in% names(out)) {
    if ("GSEA_Direction" %in% names(out)) out$Direction <- out$GSEA_Direction
    else out$Direction <- ifelse(out$NES >= 0, "Positive", "Negative")
  }
  out$Direction <- gsea_nb_direction(out$Direction)
  out
}

gsea_nb_direction <- function(x) {
  x <- as.character(x)
  lookup <- c(Positive = "Positive", Negative = "Negative", Up = "Positive",
              Down = "Negative", up = "Positive", down = "Negative")
  value <- unname(lookup[x])
  bad <- unique(x[is.na(value) & !is.na(x)])
  if (length(bad)) stop("未知 GSEA 方向：", paste(bad, collapse = ", "), call. = FALSE)
  value
}

gsea_nb_core_genes <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) return(character())
  genes <- trimws(strsplit(as.character(x), "/", fixed = TRUE)[[1]])
  unique(genes[nzchar(genes)])
}

gsea_nb_require <- function(data, columns, label = "数据") {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " 缺少列：", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

gsea_nb_safe_name <- function(x) {
  x <- gsub("[^[:alnum:]_.-]+", "_", as.character(x))
  if (!nzchar(x) || x %in% c(".", "..")) stop("不安全的输出名称。", call. = FALSE)
  x
}

gsea_nb_read_matrix <- function(path, gene_column = "Gene") {
  if (!file.exists(path)) stop("找不到表达矩阵：", path, call. = FALSE)
  data <- readr::read_tsv(path, show_col_types = FALSE, name_repair = "minimal")
  if (sum(names(data) == gene_column) != 1L) stop("表达矩阵缺少唯一的基因列：", gene_column, call. = FALSE)
  genes <- as.character(data[[gene_column]])
  samples <- names(data)[names(data) != gene_column]
  if (!length(samples) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("表达矩阵需要唯一且非空的 Gene / Sample 名称。", call. = FALSE)
  }
  mat <- as.matrix(data[, samples, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat)) || any(mat < 0)) stop("表达矩阵必须是有限的非负数。", call. = FALSE)
  rownames(mat) <- genes
  mat
}
