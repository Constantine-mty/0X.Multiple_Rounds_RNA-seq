# Interactive RNA-seq quality control. Nothing in this file runs on source().

qc_read_matrix <- function(path, gene_column = "Gene") {
  if (!file.exists(path)) stop("文件不存在：", path, call. = FALSE)
  tab <- readr::read_tsv(path, show_col_types = FALSE, name_repair = "minimal")
  if (sum(names(tab) == gene_column) != 1L || ncol(tab) < 2L) {
    stop("矩阵需要基因列 ", gene_column, " 和至少一个样本列：", path, call. = FALSE)
  }
  genes <- as.character(tab[[gene_column]])
  samples <- names(tab)[names(tab) != gene_column]
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("基因或样本 ID 不能为空、重复。", call. = FALSE)
  }
  raw <- as.matrix(tab[, samples, drop = FALSE])
  suppressWarnings(storage.mode(raw) <- "numeric")
  if (any(!is.finite(raw)) || any(raw < 0)) stop("矩阵必须是有限的非负数。", call. = FALSE)
  rownames(raw) <- genes
  raw
}

qc_prepare_data <- function(counts, metadata, tpm = NULL, fpkm = NULL,
                            groups = NULL, cycles = NULL,
                            sample_column = "Sample", group_column = "Group",
                            cycle_column = "Cycle") {
  if (!is.data.frame(metadata)) stop("metadata 必须是 data.frame。", call. = FALSE)
  required <- c(sample_column, group_column, cycle_column)
  if (!all(required %in% names(metadata))) {
    stop("metadata 缺少列：", paste(setdiff(required, names(metadata)), collapse = ", "), call. = FALSE)
  }
  meta <- as.data.frame(metadata, stringsAsFactors = FALSE)
  for (column in required) {
    meta[[column]] <- as.character(meta[[column]])
    if (anyNA(meta[[column]]) || any(!nzchar(meta[[column]]))) {
      stop("metadata 列 ", column, " 不允许空值。", call. = FALSE)
    }
  }
  if (anyDuplicated(meta[[sample_column]])) stop("metadata Sample 重复。", call. = FALSE)
  if (!is.null(groups)) meta <- meta[meta[[group_column]] %in% groups, , drop = FALSE]
  if (!is.null(cycles)) meta <- meta[meta[[cycle_column]] %in% cycles, , drop = FALSE]
  if (!nrow(meta)) stop("筛选后没有样本。", call. = FALSE)
  group_order <- if (is.null(groups)) unique(meta[[group_column]]) else groups
  cycle_order <- if (is.null(cycles)) unique(meta[[cycle_column]]) else cycles
  meta <- meta[order(match(meta[[group_column]], group_order),
                     match(meta[[cycle_column]], cycle_order), meta[[sample_column]]), , drop = FALSE]
  rownames(meta) <- meta[[sample_column]]
  validate <- function(x, label, integer_only = FALSE) {
    if (is.null(x)) return(NULL)
    x <- as.matrix(x)
    if (is.null(rownames(x)) || is.null(colnames(x)) || anyNA(rownames(x)) ||
        anyNA(colnames(x)) || any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
        anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
      stop(label, " 必须有非空且唯一的行/列名。", call. = FALSE)
    }
    suppressWarnings(storage.mode(x) <- "numeric")
    if (any(!is.finite(x)) || any(x < 0)) stop(label, " 必须是有限的非负数。", call. = FALSE)
    if (integer_only && any(x != round(x))) stop(label, " 必须是原始整数 counts。", call. = FALSE)
    missing <- setdiff(meta[[sample_column]], colnames(x))
    if (length(missing)) stop(label, " 缺少 metadata 样本：", paste(missing, collapse = ", "), call. = FALSE)
    x[, meta[[sample_column]], drop = FALSE]
  }
  counts <- validate(counts, "counts", TRUE)
  tpm <- validate(tpm, "TPM")
  fpkm <- validate(fpkm, "FPKM")
  if (!is.null(tpm) && !identical(rownames(counts), rownames(tpm))) {
    if (!setequal(rownames(counts), rownames(tpm))) stop("counts 和 TPM 基因集合不一致。", call. = FALSE)
    tpm <- tpm[rownames(counts), , drop = FALSE]
  }
  if (!is.null(fpkm) && !identical(rownames(counts), rownames(fpkm))) {
    if (!setequal(rownames(counts), rownames(fpkm))) stop("counts 和 FPKM 基因集合不一致。", call. = FALSE)
    fpkm <- fpkm[rownames(counts), , drop = FALSE]
  }
  list(counts = counts, tpm = tpm, fpkm = fpkm, metadata = meta,
       columns = list(sample = sample_column, group = group_column, cycle = cycle_column),
       group_order = group_order, cycle_order = cycle_order)
}

qc_calculate_metrics <- function(qc, tpm_thresholds = c(0.1, 1),
                                 fpkm_thresholds = c(0.1, 1),
                                 mitochondrial_pattern = "^MT-",
                                 ribosomal_pattern = "^(RPL|RPS)") {
  counts <- qc$counts
  totals <- colSums(counts)
  out <- qc$metadata
  out$total_counts <- as.numeric(totals)
  out$log10_total_counts <- log10(out$total_counts + 1)
  out$detected_genes_count_gt0 <- colSums(counts > 0)
  add_thresholds <- function(mat, prefix, thresholds) {
    if (is.null(mat)) return(invisible(NULL))
    for (threshold in thresholds) {
      if (!is.finite(threshold) || threshold < 0) stop("检测阈值必须是非负数。", call. = FALSE)
      label <- gsub("\\.", "p", format(threshold, trim = TRUE))
      out[[paste0("detected_genes_", prefix, "_gt", label)]] <<- colSums(mat > threshold)
    }
  }
  add_thresholds(qc$tpm, "tpm", tpm_thresholds)
  add_thresholds(qc$fpkm, "fpkm", fpkm_thresholds)
  fraction <- function(pattern) {
    idx <- grepl(pattern, rownames(counts))
    if (!any(idx)) warning("没有匹配 ", pattern, " 的基因；该比例设为 NA。", call. = FALSE)
    if (any(totals == 0)) warning("存在 library size 为 0 的样本；比例设为 NA。", call. = FALSE)
    if (!any(idx)) return(rep(NA_real_, ncol(counts)))
    value <- colSums(counts[idx, , drop = FALSE]) / totals * 100
    value[totals == 0] <- NA_real_
    as.numeric(value)
  }
  out$mitochondrial_pct <- fraction(mitochondrial_pattern)
  out$ribosomal_pct <- fraction(ribosomal_pattern)
  rownames(out) <- NULL
  out
}

qc_transform_counts <- function(qc, method = c("vst", "log_cpm"),
                                min_count = 10, min_samples = 3L) {
  method <- match.arg(method)
  counts <- qc$counts
  min_samples <- as.integer(min_samples)
  if (min_samples < 1L) stop("min_samples 至少为 1。", call. = FALSE)
  keep <- rowSums(counts >= min_count) >= min_samples
  if (sum(keep) < 2L) stop("通过过滤的基因少于 2 个；请调整 min_count/min_samples。", call. = FALSE)
  counts <- counts[keep, , drop = FALSE]
  if (any(colSums(counts) == 0)) stop("过滤后有零测序深度样本。", call. = FALSE)
  if (method == "log_cpm") {
    return(log2(t(t(counts) / colSums(counts) * 1e6) + 1))
  }
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("VST 需要安装 DESeq2；或手动选择 method='log_cpm'。", call. = FALSE)
  coldata <- data.frame(row.names = colnames(counts), qc$metadata, check.names = FALSE)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  as.matrix(SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)))
}

qc_plot_metric <- function(metrics, metric = "total_counts", qc,
                           group_colors = NULL, cycle_colors = NULL) {
  if (!metric %in% names(metrics) || !is.numeric(metrics[[metric]])) {
    stop("metrics 中不存在数值列：", metric, call. = FALSE)
  }
  group <- qc$columns$group; cycle <- qc$columns$cycle
  df <- metrics
  df$.group <- factor(as.character(df[[group]]), levels = qc$group_order)
  df$.cycle <- factor(as.character(df[[cycle]]), levels = qc$cycle_order)
  df$.group_cycle <- factor(paste(df$.group, df$.cycle, sep = " / "),
                            levels = as.vector(t(outer(qc$group_order, qc$cycle_order, paste, sep = " / "))))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .group_cycle, y = .data[[metric]], fill = .cycle)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.7, na.rm = TRUE) +
    ggplot2::geom_jitter(ggplot2::aes(color = .group), width = 0.12, height = 0, size = 2.5, na.rm = TRUE) +
    ggplot2::labs(x = NULL, y = metric, color = "Group", fill = "Cycle") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (metric == "total_counts") p <- p + ggplot2::scale_y_continuous(labels = scales::label_scientific())
  if (!is.null(group_colors)) p <- p + ggplot2::scale_color_manual(values = group_colors)
  if (!is.null(cycle_colors)) p <- p + ggplot2::scale_fill_manual(values = cycle_colors)
  p
}

qc_plot_pca <- function(transformed, qc, group_colors = NULL,
                        label_samples = TRUE, top_genes = NULL) {
  mat <- transformed[, qc$metadata[[qc$columns$sample]], drop = FALSE]
  if (!is.null(top_genes)) {
    if (top_genes < 2L) stop("top_genes 至少为 2。", call. = FALSE)
    variance <- apply(mat, 1, stats::var)
    mat <- mat[order(variance, decreasing = TRUE)[seq_len(min(top_genes, nrow(mat)))], , drop = FALSE]
  }
  if (ncol(mat) < 3L) stop("PCA 至少需要 3 个样本。", call. = FALSE)
  fit <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  variance <- fit$sdev^2 / sum(fit$sdev^2) * 100
  if (ncol(fit$x) < 2L || !all(is.finite(variance[1:2]))) stop("PCA 无法得到两个非零主成分。", call. = FALSE)
  df <- qc$metadata
  df$PC1 <- fit$x[, 1]; df$PC2 <- fit$x[, 2]
  df$.group <- factor(df[[qc$columns$group]], levels = qc$group_order)
  df$.cycle <- factor(df[[qc$columns$cycle]], levels = qc$cycle_order)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = .group, shape = .cycle)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]),
                  color = "Group", shape = "Cycle") + ggplot2::theme_bw(base_size = 11)
  if (label_samples) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[qc$columns$sample]]), max.overlaps = Inf, size = 3)
    } else {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = .data[[qc$columns$sample]]), vjust = -0.7, size = 3)
    }
  }
  if (!is.null(group_colors)) p <- p + ggplot2::scale_color_manual(values = group_colors)
  p
}

qc_plot_correlation <- function(transformed, qc, annotation_colors = NULL,
                                cluster_samples = TRUE, show_sample_names = TRUE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("相关性热图需要 pheatmap。", call. = FALSE)
  samples <- qc$metadata[[qc$columns$sample]]
  mat <- transformed[, samples, drop = FALSE]
  if (ncol(mat) < 2L) stop("相关性热图至少需要 2 个样本。", call. = FALSE)
  cors <- stats::cor(mat, use = "pairwise.complete.obs")
  if (any(!is.finite(cors))) stop("相关性矩阵存在 NA；请检查恒定表达样本。", call. = FALSE)
  annotation <- qc$metadata[, c(qc$columns$group, qc$columns$cycle), drop = FALSE]
  names(annotation) <- c("Group", "Cycle")
  rownames(annotation) <- samples
  annotation$Group <- factor(annotation$Group, levels = qc$group_order)
  annotation$Cycle <- factor(annotation$Cycle, levels = qc$cycle_order)
  pheatmap::pheatmap(cors, annotation_col = annotation, annotation_row = annotation,
                     annotation_colors = annotation_colors, cluster_rows = cluster_samples,
                     cluster_cols = cluster_samples, show_rownames = show_sample_names,
                     show_colnames = show_sample_names, border_color = NA, silent = TRUE)
}

qc_housekeeping_table <- function(qc, genes = c("ACTB", "GAPDH", "B2M", "RPLP0", "HPRT1",
                                               "TBP", "PPIA", "PGK1", "YWHAZ", "SDHA", "RPL13A", "TUBB")) {
  if (is.null(qc$tpm)) stop("Housekeeping QC 需要 TPM。", call. = FALSE)
  present <- intersect(unique(genes), rownames(qc$tpm))
  if (!length(present)) stop("指定 housekeeping genes 均不在 TPM 中。", call. = FALSE)
  missing <- setdiff(unique(genes), present)
  if (length(missing)) warning("Housekeeping genes 未匹配：", paste(missing, collapse = ", "), call. = FALSE)
  mat <- log2(qc$tpm[present, , drop = FALSE] + 1)
  df <- data.frame(Sample = rep(colnames(mat), each = nrow(mat)),
                   Gene = rep(rownames(mat), times = ncol(mat)),
                   log2_tpm = as.vector(mat), stringsAsFactors = FALSE)
  meta <- qc$metadata
  meta$Sample <- meta[[qc$columns$sample]]
  merge(df, meta, by = "Sample", sort = FALSE)
}

qc_plot_housekeeping <- function(qc, genes = c("ACTB", "GAPDH", "B2M", "RPLP0", "HPRT1",
                                              "TBP", "PPIA", "PGK1", "YWHAZ", "SDHA", "RPL13A", "TUBB"),
                                 summary = c("per_gene", "mean", "median"), group_colors = NULL) {
  summary <- match.arg(summary)
  df <- qc_housekeeping_table(qc, genes)
  group <- qc$columns$group; cycle <- qc$columns$cycle
  if (summary != "per_gene") {
    by_sample <- split(df$log2_tpm, df$Sample)
    value <- vapply(by_sample, if (summary == "mean") base::mean else stats::median, numeric(1))
    df <- qc$metadata
    df$Sample <- df[[qc$columns$sample]]
    df$log2_tpm <- unname(value[df$Sample])
  }
  df$.group <- factor(df[[group]], levels = qc$group_order)
  df$.cycle <- factor(df[[cycle]], levels = qc$cycle_order)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .cycle, y = log2_tpm, color = .group)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = interaction(.group, .cycle)),
                          position = ggplot2::position_dodge(width = 0.7), outlier.shape = NA, width = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 1.5) +
    ggplot2::labs(x = "Cycle", y = "log2(TPM + 1)", color = "Group") + ggplot2::theme_bw(base_size = 11)
  if (summary == "per_gene") p <- p + ggplot2::facet_wrap(~Gene, scales = "free_y")
  if (!is.null(group_colors)) p <- p + ggplot2::scale_color_manual(values = group_colors)
  p
}
