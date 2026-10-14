# DESeq2 stage: fit one model for each requested within-group cycle contrast.

build_design_formula <- function(metadata, cycle_column, design_covariates) {
  design_columns <- unique(c(design_covariates, cycle_column))
  missing <- setdiff(design_columns, colnames(metadata))
  if (length(missing)) {
    abort("设计公式所需的 metadata 列不存在：", paste(missing, collapse = ", "))
  }
  variable_covariates <- design_covariates[vapply(
    design_covariates,
    function(x) length(unique(metadata[[x]])) > 1L,
    logical(1)
  )]
  dropped <- setdiff(design_covariates, variable_covariates)
  if (length(dropped)) {
    message("提示：当前比较中协变量无变化，已从模型移除：", paste(dropped, collapse = ", "))
  }
  stats::as.formula(paste("~", paste(c(variable_covariates, cycle_column), collapse = " + ")))
}

# write_deseq_qc <- function(dds, result, metadata, cycle_column, figure_dir, comparison_label) {
#   require_package("ggplot2")
#   ensure_dir(figure_dir)
# 
#   grDevices::pdf(file.path(figure_dir, "ma_plot.pdf"), width = 7, height = 6)
#   try(DESeq2::plotMA(result, main = comparison_label), silent = TRUE)
#   grDevices::dev.off()
# 
#   if (ncol(dds) >= 3L && nrow(dds) >= 2L) {
#     grDevices::pdf(file.path(figure_dir, "pca.pdf"), width = 7, height = 6)
#     try({
#       vsd <- DESeq2::vst(dds, blind = FALSE)
#       print(DESeq2::plotPCA(vsd, intgroup = cycle_column) + ggplot2::ggtitle(comparison_label))
#     }, silent = TRUE)
#     grDevices::dev.off()
#   }
# }

run_one_deseq2_contrast <- function(
    counts,
    metadata,
    contrast_row,
    input_config,
    deseq_config,
    project_root) {
  require_package("DESeq2")

  sample_column <- input_config$sample_column
  group_column <- input_config$group_column
  cycle_column <- input_config$cycle_column
  group_id <- contrast_row$group_id
  comparison_id <- contrast_row$comparison_id
  numerator <- contrast_row$numerator
  denominator <- contrast_row$denominator

  result_rel <- relative_result_path(
    "data", "de_results", safe_path_component(group_id),
    safe_path_component(comparison_id), "deseq2_results.tsv"
  )
  normalized_rel <- relative_result_path(
    "data", "de_results", safe_path_component(group_id),
    safe_path_component(comparison_id), "normalized_counts.tsv"
  )
  metadata_rel <- relative_result_path(
    "data", "de_results", safe_path_component(group_id),
    safe_path_component(comparison_id), "metadata_used.tsv"
  )
  figure_dir <- file.path(
    project_root, "figures", "de_results", safe_path_component(group_id),
    safe_path_component(comparison_id)
  )

  group_meta <- metadata[metadata[[group_column]] == group_id, , drop = FALSE]
  pair_meta <- group_meta[group_meta[[cycle_column]] %in% c(numerator, denominator), , drop = FALSE]
  n_numerator <- sum(pair_meta[[cycle_column]] == numerator)
  n_denominator <- sum(pair_meta[[cycle_column]] == denominator)
  min_group_samples <- as.integer(deseq_config$min_group_samples %||% 2L)

  manifest_row <- data.frame(
    group_id = group_id,
    comparison_id = comparison_id,
    numerator = numerator,
    denominator = denominator,
    status = "skipped",
    reason = NA_character_,
    n_numerator = n_numerator,
    n_denominator = n_denominator,
    n_genes_tested = NA_integer_,
    result_file = result_rel,
    normalized_counts_file = normalized_rel,
    metadata_file = metadata_rel,
    stringsAsFactors = FALSE
  )

  if (!nrow(group_meta)) {
    manifest_row$reason <- "group_id 不存在于 metadata"
    return(manifest_row)
  }
  if (identical(numerator, denominator)) {
    manifest_row$reason <- "numerator 与 denominator 不能相同"
    return(manifest_row)
  }
  if (n_numerator < min_group_samples || n_denominator < min_group_samples) {
    manifest_row$reason <- paste0(
      "生物学重复不足：", numerator, "=", n_numerator, "；", denominator, "=", n_denominator,
      "；每组至少需要 ", min_group_samples
    )
    return(manifest_row)
  }

  pair_meta <- pair_meta[order(pair_meta[[cycle_column]], pair_meta[[sample_column]]), , drop = FALSE]
  pair_meta[[cycle_column]] <- factor(pair_meta[[cycle_column]], levels = c(denominator, numerator))
  rownames(pair_meta) <- pair_meta[[sample_column]]
  pair_counts <- counts[, pair_meta[[sample_column]], drop = FALSE]
  design <- build_design_formula(pair_meta, cycle_column, deseq_config$design_covariates %||% character())
  design_matrix <- stats::model.matrix(design, data = pair_meta)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    manifest_row$reason <- paste0("设计矩阵不满秩（轮次与协变量可能完全混杂）：", deparse(design))
    return(manifest_row)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pair_counts, colData = pair_meta, design = design)
  min_count <- as.numeric(deseq_config$min_count %||% 10)
  min_samples <- as.integer(deseq_config$min_samples %||% 3L)
  keep_genes <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
  dds <- dds[keep_genes, ]
  if (!nrow(dds)) {
    manifest_row$reason <- "没有基因通过低表达过滤"
    return(manifest_row)
  }

  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  result <- DESeq2::results(
    dds,
    contrast = c(cycle_column, numerator, denominator),
    alpha = as.numeric(deseq_config$alpha %||% 0.05)
  )
  result_df <- as.data.frame(result)
  result_df$Gene <- rownames(result_df)
  result_df <- result_df[, c("Gene", setdiff(colnames(result_df), "Gene")), drop = FALSE]
  result_df <- result_df[order(is.na(result_df$padj), result_df$padj), , drop = FALSE]

  write_tsv(result_df, file.path(project_root, result_rel))
  normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  normalized_counts$Gene <- rownames(normalized_counts)
  normalized_counts <- normalized_counts[, c("Gene", setdiff(colnames(normalized_counts), "Gene")), drop = FALSE]
  write_tsv(normalized_counts, file.path(project_root, normalized_rel))
  write_tsv(pair_meta, file.path(project_root, metadata_rel), row.names = TRUE)
  # write_deseq_qc(dds, result, pair_meta, cycle_column, figure_dir, paste(numerator, "vs", denominator))

  manifest_row$status <- "success"
  manifest_row$n_genes_tested <- nrow(result_df)
  manifest_row
}

run_deseq2_stage <- function(config, project_root) {
  require_package("DESeq2")
  require_package("ggplot2")
  input_config <- config$input
  deseq_config <- config$deseq2 %||% list()
  counts <- read_counts_matrix(resolve_path(project_root, input_config$counts_file), input_config$gene_column)
  metadata <- read_metadata(
    resolve_path(project_root, input_config$metadata_file), input_config$sample_column,
    input_config$group_column, input_config$cycle_column
  )
  counts <- align_counts_and_metadata(counts, metadata, input_config$sample_column)
  contrasts <- read_contrasts(resolve_path(project_root, input_config$contrasts_file))
  contrasts <- contrasts[contrasts$enabled, , drop = FALSE]
  if (!nrow(contrasts)) abort("contrasts 文件中没有 enabled=TRUE 的比较。")

  rows <- lapply(seq_len(nrow(contrasts)), function(i) {
    message("DESeq2: ", contrasts$group_id[i], " / ", contrasts$comparison_id[i])
    tryCatch(
      run_one_deseq2_contrast(counts, metadata, contrasts[i, , drop = FALSE], input_config, deseq_config, project_root),
      error = function(e) data.frame(
        group_id = contrasts$group_id[i], comparison_id = contrasts$comparison_id[i],
        numerator = contrasts$numerator[i], denominator = contrasts$denominator[i],
        status = "failed", reason = conditionMessage(e),
        n_numerator = NA_integer_, n_denominator = NA_integer_, n_genes_tested = NA_integer_,
        result_file = NA_character_, normalized_counts_file = NA_character_, metadata_file = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  })
  manifest <- do.call(rbind, rows)
  write_tsv(manifest, file.path(project_root, "data", "pipeline", "de_manifest.tsv"))
  message("DESeq2 完成：", sum(manifest$status == "success"), "/", nrow(manifest), " 个比较成功。")
  manifest
}
