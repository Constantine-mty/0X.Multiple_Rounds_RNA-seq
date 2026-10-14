# GSVA stage: score project-defined pathways, then test cycle differences within each group.

read_gsva_expression_matrix <- function(path, gene_column) {
  expression_df <- read_tsv(path, required_columns = gene_column)
  if (ncol(expression_df) < 2L) abort("GSVA 表达矩阵至少需要一列基因 ID 和一列样本。")
  genes <- as.character(expression_df[[gene_column]])
  if (anyNA(genes) || any(!nzchar(genes))) abort("GSVA 表达矩阵中存在空基因 ID。")
  if (anyDuplicated(genes)) abort("GSVA 表达矩阵中存在重复基因 ID。")
  sample_columns <- setdiff(colnames(expression_df), gene_column)
  if (anyDuplicated(sample_columns)) abort("GSVA 表达矩阵中存在重复样本列名。")
  values <- as.matrix(expression_df[, sample_columns, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  if (any(!is.finite(values))) abort("GSVA 表达矩阵含有非数值或缺失值。")
  if (any(values < 0)) abort("GSVA 表达矩阵不能含负值。")
  rownames(values) <- genes
  values
}

align_matrix_to_metadata <- function(matrix, metadata, sample_column, label) {
  missing_samples <- setdiff(metadata[[sample_column]], colnames(matrix))
  if (length(missing_samples)) {
    abort(label, " 中缺少 metadata 样本：", paste(missing_samples, collapse = ", "))
  }
  unused_samples <- setdiff(colnames(matrix), metadata[[sample_column]])
  if (length(unused_samples)) {
    message("提示：", label, " 中有 ", length(unused_samples), " 个样本未在 metadata 中出现，将不参与 GSVA。")
  }
  matrix[, metadata[[sample_column]], drop = FALSE]
}

build_gsva_gene_sets <- function(pathways, expressed_genes, gsva_config) {
  all_sets <- split(pathways$term2gene$gene, pathways$term2gene$term)
  all_sets <- lapply(all_sets, function(genes) intersect(unique(genes), expressed_genes))
  min_size <- as.integer(gsva_config$min_gs_size %||% 10L)
  max_size <- as.integer(gsva_config$max_gs_size %||% 300L)
  keep <- lengths(all_sets) >= min_size & lengths(all_sets) <= max_size
  gene_sets <- all_sets[keep]
  if (!length(gene_sets)) {
    abort("没有通路在表达矩阵中拥有 ", min_size, "–", max_size, " 个基因，无法进行 GSVA。")
  }
  gene_sets
}

run_gsva_scoring <- function(expression_log, gene_sets, gsva_config) {
  require_package("GSVA")
  kcdf <- gsva_config$kcdf %||% "Gaussian"
  min_size <- as.integer(gsva_config$min_gs_size %||% 10L)
  max_size <- as.integer(gsva_config$max_gs_size %||% 300L)
  if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    parameters <- GSVA::gsvaParam(
      exprData = expression_log,
      geneSets = gene_sets,
      kcdf = kcdf,
      minSize = min_size,
      maxSize = max_size
    )
    scores <- GSVA::gsva(parameters)
  } else {
    scores <- GSVA::gsva(
      expr = expression_log,
      gset.idx.list = gene_sets,
      method = "gsva",
      kcdf = kcdf,
      min.sz = min_size,
      max.sz = max_size,
      verbose = FALSE
    )
  }
  scores <- as.matrix(scores)
  if (!nrow(scores) || !ncol(scores)) abort("GSVA 未返回有效的通路打分矩阵。")
  scores
}

pathway_description <- function(pathway_ids, pathways) {
  names_lookup <- setNames(pathways$term2name$name, pathways$term2name$term)
  description <- unname(names_lookup[pathway_ids])
  description[is.na(description)] <- pathway_ids[is.na(description)]
  description
}

gsva_scores_table <- function(scores, pathways) {
  out <- as.data.frame(scores, check.names = FALSE)
  out$Pathway <- rownames(scores)
  out$Description <- pathway_description(out$Pathway, pathways)
  out[, c("Pathway", "Description", setdiff(colnames(out), c("Pathway", "Description"))), drop = FALSE]
}

write_gsva_scoring_outputs <- function(
    scores,
    pathways,
    gene_sets,
    expressed_genes,
    project_root,
    output_dir_rel,
    input_manifest) {
  score_rel <- file.path(output_dir_rel, "gsva_scores.tsv")
  score_rds_rel <- file.path(output_dir_rel, "gsva_scores.rds")
  gene_set_rel <- file.path(output_dir_rel, "gene_sets_used.tsv")
  input_manifest_rel <- file.path(output_dir_rel, "input_manifest.tsv")
  write_rds(scores, file.path(project_root, score_rds_rel))
  write_tsv(gsva_scores_table(scores, pathways), file.path(project_root, score_rel))
  used_mapping <- pathways$term2gene[
    pathways$term2gene$term %in% names(gene_sets) & pathways$term2gene$gene %in% expressed_genes,
    , drop = FALSE
  ]
  used_mapping$Description <- pathway_description(used_mapping$term, pathways)
  write_tsv(used_mapping, file.path(project_root, gene_set_rel))
  write_tsv(input_manifest, file.path(project_root, input_manifest_rel))
  list(
    score_file = score_rel,
    score_rds_file = score_rds_rel,
    gene_sets_file = gene_set_rel,
    input_manifest_file = input_manifest_rel
  )
}

build_gsva_design <- function(metadata, cycle_column, cycle_levels, design_covariates) {
  metadata$.gsva_cycle <- factor(as.character(metadata[[cycle_column]]), levels = cycle_levels)
  if (anyNA(metadata$.gsva_cycle)) abort("GSVA 模型中出现未定义的轮次。")
  missing_covariates <- setdiff(design_covariates, colnames(metadata))
  if (length(missing_covariates)) {
    abort("GSVA 设计公式所需的 metadata 列不存在：", paste(missing_covariates, collapse = ", "))
  }
  varying_covariates <- design_covariates[vapply(
    design_covariates,
    function(column) length(unique(metadata[[column]])) > 1L,
    logical(1)
  )]
  dropped_covariates <- setdiff(design_covariates, varying_covariates)
  if (length(dropped_covariates)) {
    message("提示：当前 GSVA 检验中协变量无变化，已从模型移除：", paste(dropped_covariates, collapse = ", "))
  }
  formula <- stats::as.formula(paste("~ 0 +", paste(c(".gsva_cycle", varying_covariates), collapse = " + ")))
  design <- stats::model.matrix(formula, data = metadata)
  cycle_columns <- setNames(paste0(".gsva_cycle", make.names(cycle_levels)), cycle_levels)
  if (!all(cycle_columns %in% colnames(design))) {
    abort("无法在 GSVA 设计矩阵中定位轮次列。")
  }
  list(metadata = metadata, formula = formula, design = design, cycle_columns = cycle_columns)
}

classify_gsva_result <- function(result, gsva_config) {
  padj_cutoff <- as.numeric(gsva_config$padj_cutoff %||% 0.05)
  logfc_cutoff <- as.numeric(gsva_config$logfc_cutoff %||% 0)
  result$Direction <- ifelse(
    !is.na(result$adj.P.Val) & result$adj.P.Val <= padj_cutoff & result$logFC > logfc_cutoff,
    "Up",
    ifelse(
      !is.na(result$adj.P.Val) & result$adj.P.Val <= padj_cutoff & result$logFC < -logfc_cutoff,
      "Down", "NS"
    )
  )
  result
}

extract_gsva_limma_result <- function(fit, coefficient, pathways, gsva_config) {
  result <- limma::topTable(
    fit,
    coef = coefficient,
    number = Inf,
    adjust.method = gsva_config$p_adjust_method %||% "BH",
    sort.by = "P"
  )
  result$Pathway <- rownames(result)
  result$Description <- pathway_description(result$Pathway, pathways)
  result <- result[, c("Pathway", "Description", setdiff(colnames(result), c("Pathway", "Description"))), drop = FALSE]
  classify_gsva_result(result, gsva_config)
}

run_one_gsva_pairwise <- function(
    scores,
    metadata,
    contrast_row,
    input_config,
    gsva_config,
    pathways,
    project_root,
    result_root_rel = file.path("data", "gsva_results"),
    scoring_scope = "global") {
  sample_column <- input_config$sample_column
  group_column <- input_config$group_column
  cycle_column <- input_config$cycle_column
  group_id <- contrast_row$group_id
  comparison_id <- contrast_row$comparison_id
  numerator <- contrast_row$numerator
  denominator <- contrast_row$denominator
  result_dir_rel <- file.path(result_root_rel, safe_path_component(group_id), safe_path_component(comparison_id))
  result_rel <- file.path(result_dir_rel, "limma_results.tsv")
  significant_rel <- file.path(result_dir_rel, "limma_significant.tsv")
  fit_rel <- file.path(result_dir_rel, "limma_fit.rds")
  group_metadata <- metadata[metadata[[group_column]] == group_id, , drop = FALSE]
  pair_metadata <- group_metadata[group_metadata[[cycle_column]] %in% c(numerator, denominator), , drop = FALSE]
  n_numerator <- sum(pair_metadata[[cycle_column]] == numerator)
  n_denominator <- sum(pair_metadata[[cycle_column]] == denominator)
  min_group_samples <- as.integer(gsva_config$min_group_samples %||% 2L)
  manifest <- data.frame(
    scoring_scope = scoring_scope,
    group_id = group_id,
    comparison_id = comparison_id,
    numerator = numerator,
    denominator = denominator,
    status = "skipped",
    reason = NA_character_,
    n_numerator = n_numerator,
    n_denominator = n_denominator,
    n_pathways = NA_integer_,
    n_up = NA_integer_,
    n_down = NA_integer_,
    result_file = result_rel,
    significant_file = significant_rel,
    fit_file = fit_rel,
    stringsAsFactors = FALSE
  )
  if (!nrow(group_metadata)) {
    manifest$reason <- "group_id 不存在于 metadata"
    return(manifest)
  }
  if (identical(numerator, denominator)) {
    manifest$reason <- "numerator 与 denominator 不能相同"
    return(manifest)
  }
  if (n_numerator < min_group_samples || n_denominator < min_group_samples) {
    manifest$reason <- paste0("生物学重复不足：", numerator, "=", n_numerator, "；", denominator, "=", n_denominator)
    return(manifest)
  }
  pair_metadata <- pair_metadata[order(pair_metadata[[cycle_column]], pair_metadata[[sample_column]]), , drop = FALSE]
  rownames(pair_metadata) <- pair_metadata[[sample_column]]
  score_pair <- scores[, pair_metadata[[sample_column]], drop = FALSE]
  design_info <- build_gsva_design(
    pair_metadata, cycle_column, c(denominator, numerator), gsva_config$design_covariates %||% character()
  )
  if (qr(design_info$design)$rank < ncol(design_info$design)) {
    manifest$reason <- paste0("设计矩阵不满秩（轮次与协变量可能完全混杂）：", deparse(design_info$formula))
    return(manifest)
  }
  contrast_vector <- setNames(rep(0, ncol(design_info$design)), colnames(design_info$design))
  contrast_vector[design_info$cycle_columns[[numerator]]] <- 1
  contrast_vector[design_info$cycle_columns[[denominator]]] <- -1
  contrast_matrix <- matrix(contrast_vector, ncol = 1L, dimnames = list(colnames(design_info$design), comparison_id))
  fit <- limma::lmFit(score_pair, design_info$design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit)
  result <- extract_gsva_limma_result(fit, comparison_id, pathways, gsva_config)
  significant <- result[result$Direction != "NS", , drop = FALSE]
  write_tsv(result, file.path(project_root, result_rel))
  write_tsv(significant, file.path(project_root, significant_rel))
  write_rds(
    list(fit = fit, design = design_info$design, contrast = contrast_matrix, metadata = design_info$metadata, result = result),
    file.path(project_root, fit_rel)
  )
  manifest$status <- "success"
  manifest$n_pathways <- nrow(result)
  manifest$n_up <- sum(result$Direction == "Up")
  manifest$n_down <- sum(result$Direction == "Down")
  manifest
}

calculate_gsva_cycle_means <- function(scores, group_metadata, input_config, pathways) {
  sample_column <- input_config$sample_column
  cycle_column <- input_config$cycle_column
  cycle_levels <- unique(as.character(group_metadata[[cycle_column]]))
  mean_matrix <- sapply(cycle_levels, function(cycle) {
    samples <- group_metadata[[sample_column]][group_metadata[[cycle_column]] == cycle]
    rowMeans(scores[, samples, drop = FALSE])
  })
  mean_matrix <- as.matrix(mean_matrix)
  if (is.null(colnames(mean_matrix))) colnames(mean_matrix) <- cycle_levels
  out <- as.data.frame(mean_matrix, check.names = FALSE)
  out$Pathway <- rownames(scores)
  out$Description <- pathway_description(out$Pathway, pathways)
  out[, c("Pathway", "Description", cycle_levels), drop = FALSE]
}

run_one_gsva_specific <- function(
    scores,
    metadata,
    group_id,
    input_config,
    gsva_config,
    pathways,
    project_root,
    result_root_rel = file.path("data", "gsva_results"),
    scoring_scope = "global") {
  sample_column <- input_config$sample_column
  group_column <- input_config$group_column
  cycle_column <- input_config$cycle_column
  result_dir_rel <- file.path(result_root_rel, safe_path_component(group_id))
  result_rel <- file.path(result_dir_rel, "specific_results.tsv")
  peak_rel <- file.path(result_dir_rel, "specific_peak_significant.tsv")
  means_rel <- file.path(result_dir_rel, "cycle_means.tsv")
  fit_rel <- file.path(result_dir_rel, "specific_limma_fit.rds")
  group_metadata <- metadata[metadata[[group_column]] == group_id, , drop = FALSE]
  manifest <- data.frame(
    scoring_scope = scoring_scope, group_id = group_id, status = "skipped", reason = NA_character_, n_cycles = NA_integer_,
    n_pathways = NA_integer_, n_peak_specific = NA_integer_, result_file = result_rel,
    peak_significant_file = peak_rel, cycle_means_file = means_rel, fit_file = fit_rel,
    stringsAsFactors = FALSE
  )
  if (!nrow(group_metadata)) {
    manifest$reason <- "group_id 不存在于 metadata"
    return(manifest)
  }
  cycle_counts <- table(group_metadata[[cycle_column]])
  min_group_samples <- as.integer(gsva_config$min_group_samples %||% 2L)
  usable_cycles <- names(cycle_counts)[cycle_counts >= min_group_samples]
  if (length(usable_cycles) < 2L) {
    manifest$reason <- "至少需要两个轮次且每轮次满足最小重复数，才能进行特异性检验"
    return(manifest)
  }
  group_metadata <- group_metadata[group_metadata[[cycle_column]] %in% usable_cycles, , drop = FALSE]
  group_metadata <- group_metadata[order(group_metadata[[cycle_column]], group_metadata[[sample_column]]), , drop = FALSE]
  rownames(group_metadata) <- group_metadata[[sample_column]]
  score_group <- scores[, group_metadata[[sample_column]], drop = FALSE]
  design_info <- build_gsva_design(
    group_metadata, cycle_column, usable_cycles, gsva_config$design_covariates %||% character()
  )
  if (qr(design_info$design)$rank < ncol(design_info$design)) {
    manifest$reason <- paste0("设计矩阵不满秩（轮次与协变量可能完全混杂）：", deparse(design_info$formula))
    return(manifest)
  }
  contrast_matrix <- sapply(usable_cycles, function(target) {
    vector <- setNames(rep(0, ncol(design_info$design)), colnames(design_info$design))
    vector[design_info$cycle_columns[[target]]] <- 1
    others <- setdiff(usable_cycles, target)
    vector[design_info$cycle_columns[others]] <- -1 / length(others)
    vector
  })
  contrast_matrix <- as.matrix(contrast_matrix)
  rownames(contrast_matrix) <- colnames(design_info$design)
  colnames(contrast_matrix) <- paste0(usable_cycles, "_specific")
  fit <- limma::lmFit(score_group, design_info$design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit)
  specific_results <- do.call(rbind, lapply(seq_along(usable_cycles), function(index) {
    result <- extract_gsva_limma_result(fit, colnames(contrast_matrix)[index], pathways, gsva_config)
    result$Specific_Cycle <- usable_cycles[index]
    result
  }))
  cycle_means <- calculate_gsva_cycle_means(scores, group_metadata, input_config, pathways)
  mean_matrix <- as.matrix(cycle_means[match(specific_results$Pathway, cycle_means$Pathway), usable_cycles, drop = FALSE])
  storage.mode(mean_matrix) <- "numeric"
  specific_index <- match(specific_results$Specific_Cycle, usable_cycles)
  peak_index <- apply(mean_matrix, 1, function(values) {
    if (anyNA(values) || sum(values == max(values)) != 1L) return(NA_integer_)
    which.max(values)
  })
  specific_results$Peak_Cycle <- ifelse(is.na(peak_index), NA_character_, usable_cycles[peak_index])
  specific_results$is_peak_specific <- !is.na(peak_index) & peak_index == specific_index
  peak_significant <- specific_results[
    specific_results$Direction != "NS" & specific_results$is_peak_specific,
    , drop = FALSE
  ]
  write_tsv(cycle_means, file.path(project_root, means_rel))
  write_tsv(specific_results, file.path(project_root, result_rel))
  write_tsv(peak_significant, file.path(project_root, peak_rel))
  write_rds(
    list(fit = fit, design = design_info$design, contrast = contrast_matrix, metadata = design_info$metadata,
         cycle_means = cycle_means, result = specific_results),
    file.path(project_root, fit_rel)
  )
  manifest$status <- "success"
  manifest$n_cycles <- length(usable_cycles)
  manifest$n_pathways <- nrow(scores)
  manifest$n_peak_specific <- nrow(peak_significant)
  manifest
}

write_gsva_pairwise_summary_figure <- function(
    manifest,
    project_root,
    figure_root_rel = file.path("figures", "gsva_results")) {
  successful <- manifest[manifest$status == "success", , drop = FALSE]
  if (!nrow(successful)) return(invisible(NULL))
  for (group_id in unique(successful$group_id)) {
    group_results <- successful[successful$group_id == group_id, , drop = FALSE]
    counts <- rbind(Up = group_results$n_up, Down = -group_results$n_down)
    figure_dir <- file.path(project_root, figure_root_rel, safe_path_component(group_id))
    ensure_dir(figure_dir)
    grDevices::pdf(file.path(figure_dir, "pairwise_significant_pathway_counts.pdf"), width = 8, height = 5)
    graphics::barplot(
      counts, beside = TRUE, col = c("#B2182B", "#2166AC"),
      names.arg = group_results$comparison_id, ylab = "Pathway count (Down shown as negative)",
      main = paste0("GSVA significant pathways: ", group_id)
    )
    graphics::abline(h = 0, lty = 2)
    graphics::legend("topright", legend = c("Up", "Down"), fill = c("#B2182B", "#2166AC"), bty = "n")
    grDevices::dev.off()
  }
  invisible(NULL)
}

run_gsva_stage <- function(config, project_root) {
  if (is.null(config$gsva)) abort("配置文件缺少 gsva 节。")
  gsva_config <- config$gsva
  if (!isTRUE(gsva_config$enabled %||% TRUE)) {
    message("GSVA 已在配置中禁用。")
    return(invisible(NULL))
  }
  require_package("GSVA")
  require_package("limma")
  input_config <- config$input
  if (is.null(gsva_config$tpm_file) || !nzchar(gsva_config$tpm_file)) {
    abort("gsva 节缺少 tpm_file。")
  }
  counts <- read_counts_matrix(resolve_path(project_root, input_config$counts_file), input_config$gene_column)
  metadata <- read_metadata(
    resolve_path(project_root, input_config$metadata_file), input_config$sample_column,
    input_config$group_column, input_config$cycle_column
  )
  counts <- align_matrix_to_metadata(counts, metadata, input_config$sample_column, "count 矩阵")
  tpm <- read_gsva_expression_matrix(resolve_path(project_root, gsva_config$tpm_file), input_config$gene_column)
  tpm <- align_matrix_to_metadata(tpm, metadata, input_config$sample_column, "TPM 矩阵")
  common_genes <- intersect(rownames(counts), rownames(tpm))
  if (!length(common_genes)) abort("count 与 TPM 矩阵没有共同的基因 ID。")
  counts <- counts[common_genes, , drop = FALSE]
  tpm <- tpm[common_genes, , drop = FALSE]
  count_cutoff <- as.numeric(gsva_config$count_filter_min %||% 10)
  min_samples <- as.integer(gsva_config$count_filter_min_samples %||% 3L)
  keep_genes <- rowSums(counts > count_cutoff) >= min_samples
  if (!any(keep_genes)) abort("没有基因通过 GSVA 的 counts 过滤。")
  pseudocount <- as.numeric(gsva_config$log_pseudocount %||% 1)
  if (pseudocount <= 0) abort("gsva.log_pseudocount 必须大于 0。")
  expression_log <- log2(tpm[keep_genes, , drop = FALSE] + pseudocount)
  pathways <- load_pathway_definitions(config$enrichment %||% list(), project_root)
  gene_sets <- build_gsva_gene_sets(pathways, rownames(expression_log), gsva_config)
  scores <- run_gsva_scoring(expression_log, gene_sets, gsva_config)
  score_rel <- file.path("data", "gsva", "gsva_scores.tsv")
  score_rds_rel <- file.path("data", "gsva", "gsva_scores.rds")
  gene_set_rel <- file.path("data", "gsva", "gene_sets_used.tsv")
  input_manifest_rel <- file.path("data", "gsva", "input_manifest.tsv")
  write_rds(scores, file.path(project_root, score_rds_rel))
  write_tsv(gsva_scores_table(scores, pathways), file.path(project_root, score_rel))
  used_mapping <- pathways$term2gene[pathways$term2gene$term %in% names(gene_sets) & pathways$term2gene$gene %in% rownames(expression_log), , drop = FALSE]
  used_mapping$Description <- pathway_description(used_mapping$term, pathways)
  write_tsv(used_mapping, file.path(project_root, gene_set_rel))
  input_manifest <- data.frame(
    scoring_scope = "global",
    scoring_group = NA_character_,
    counts_file = resolve_path(project_root, input_config$counts_file),
    tpm_file = resolve_path(project_root, gsva_config$tpm_file),
    n_common_genes = length(common_genes), n_genes_after_count_filter = sum(keep_genes),
    n_gene_sets = length(gene_sets), n_samples = ncol(scores),
    stringsAsFactors = FALSE
  )
  write_tsv(input_manifest, file.path(project_root, input_manifest_rel))
  contrasts <- read_contrasts(resolve_path(project_root, input_config$contrasts_file))
  contrasts <- contrasts[contrasts$enabled, , drop = FALSE]
  if (!nrow(contrasts)) abort("contrasts 文件中没有 enabled=TRUE 的比较。")
  pair_rows <- lapply(seq_len(nrow(contrasts)), function(index) {
    message("GSVA pairwise: ", contrasts$group_id[index], " / ", contrasts$comparison_id[index])
    tryCatch(
      run_one_gsva_pairwise(scores, metadata, contrasts[index, , drop = FALSE], input_config, gsva_config, pathways, project_root),
      error = function(error) data.frame(
        scoring_scope = "global",
        group_id = contrasts$group_id[index], comparison_id = contrasts$comparison_id[index],
        numerator = contrasts$numerator[index], denominator = contrasts$denominator[index],
        status = "failed", reason = conditionMessage(error), n_numerator = NA_integer_, n_denominator = NA_integer_,
        n_pathways = NA_integer_, n_up = NA_integer_, n_down = NA_integer_, result_file = NA_character_,
        significant_file = NA_character_, fit_file = NA_character_, stringsAsFactors = FALSE
      )
    )
  })
  pair_manifest <- do.call(rbind, pair_rows)
  write_tsv(pair_manifest, file.path(project_root, "data", "pipeline", "gsva_manifest.tsv"))
  write_gsva_pairwise_summary_figure(pair_manifest, project_root)
  if (isTRUE(gsva_config$run_specific_tests %||% TRUE)) {
    group_ids <- unique(contrasts$group_id)
    specific_rows <- lapply(group_ids, function(group_id) {
      message("GSVA specific: ", group_id)
      tryCatch(
        run_one_gsva_specific(scores, metadata, group_id, input_config, gsva_config, pathways, project_root),
      error = function(error) data.frame(
          scoring_scope = "global", group_id = group_id, status = "failed", reason = conditionMessage(error), n_cycles = NA_integer_,
          n_pathways = NA_integer_, n_peak_specific = NA_integer_, result_file = NA_character_,
          peak_significant_file = NA_character_, cycle_means_file = NA_character_, fit_file = NA_character_,
          stringsAsFactors = FALSE
        )
      )
    })
    write_tsv(do.call(rbind, specific_rows), file.path(project_root, "data", "pipeline", "gsva_specific_manifest.tsv"))
  }
  groupwise <- run_gsva_groupwise_stage(
    expression_log = expression_log,
    metadata = metadata,
    contrasts = contrasts,
    input_config = input_config,
    gsva_config = gsva_config,
    pathways = pathways,
    gene_sets = gene_sets,
    project_root = project_root,
    input_summary = input_manifest
  )
  message(
    "GSVA 完成：全样本打分分支 ", sum(pair_manifest$status == "success"), "/", nrow(pair_manifest),
    " 个组内比较成功；分组独立打分分支 ", groupwise$n_success, "/", groupwise$n_total, " 个组内比较成功。"
  )
  pair_manifest
}

run_gsva_groupwise_stage <- function(
    expression_log,
    metadata,
    contrasts,
    input_config,
    gsva_config,
    pathways,
    gene_sets,
    project_root,
    input_summary) {
  if (!isTRUE(gsva_config$run_groupwise_scoring %||% TRUE)) {
    message("GSVA 分组独立打分分支已在配置中禁用。")
    return(list(n_success = 0L, n_total = 0L))
  }

  sample_column <- input_config$sample_column
  group_column <- input_config$group_column
  group_ids <- unique(contrasts$group_id)
  score_rows <- list()
  pair_rows <- list()
  specific_rows <- list()
  pair_index <- 1L
  score_index <- 1L
  specific_index <- 1L

  for (group_id in group_ids) {
    group_contrasts <- contrasts[contrasts$group_id == group_id, , drop = FALSE]
    group_metadata <- metadata[metadata[[group_column]] == group_id, , drop = FALSE]
    score_dir_rel <- file.path("data", "gsva", "by_group", safe_path_component(group_id))
    result_root_rel <- file.path("data", "gsva_results", "by_group")

    tryCatch({
      if (!nrow(group_metadata)) abort("group_id 不存在于 metadata")
      group_expression <- expression_log[, group_metadata[[sample_column]], drop = FALSE]
      group_scores <- run_gsva_scoring(group_expression, gene_sets, gsva_config)
      group_input_manifest <- input_summary
      group_input_manifest$scoring_scope <- "within_group"
      group_input_manifest$scoring_group <- group_id
      group_input_manifest$n_samples <- ncol(group_scores)
      group_input_manifest$n_gene_sets <- nrow(group_scores)
      group_artifacts <- write_gsva_scoring_outputs(
        scores = group_scores,
        pathways = pathways,
        gene_sets = gene_sets,
        expressed_genes = rownames(group_expression),
        project_root = project_root,
        output_dir_rel = score_dir_rel,
        input_manifest = group_input_manifest
      )
      score_rows[[score_index]] <- data.frame(
        scoring_scope = "within_group", group_id = group_id, status = "success", reason = NA_character_,
        n_samples = ncol(group_scores), n_pathways = nrow(group_scores),
        score_file = group_artifacts$score_file, score_rds_file = group_artifacts$score_rds_file,
        gene_sets_file = group_artifacts$gene_sets_file, input_manifest_file = group_artifacts$input_manifest_file,
        stringsAsFactors = FALSE
      )
      score_index <- score_index + 1L

      for (contrast_index in seq_len(nrow(group_contrasts))) {
        message("GSVA within-group scoring: ", group_id, " / ", group_contrasts$comparison_id[contrast_index])
        pair_rows[[pair_index]] <- tryCatch(
          run_one_gsva_pairwise(
            scores = group_scores,
            metadata = metadata,
            contrast_row = group_contrasts[contrast_index, , drop = FALSE],
            input_config = input_config,
            gsva_config = gsva_config,
            pathways = pathways,
            project_root = project_root,
            result_root_rel = result_root_rel,
            scoring_scope = "within_group"
          ),
          error = function(error) data.frame(
            scoring_scope = "within_group", group_id = group_id,
            comparison_id = group_contrasts$comparison_id[contrast_index],
            numerator = group_contrasts$numerator[contrast_index], denominator = group_contrasts$denominator[contrast_index],
            status = "failed", reason = conditionMessage(error), n_numerator = NA_integer_, n_denominator = NA_integer_,
            n_pathways = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
            result_file = NA_character_, significant_file = NA_character_, fit_file = NA_character_, stringsAsFactors = FALSE
          )
        )
        pair_index <- pair_index + 1L
      }

      if (isTRUE(gsva_config$run_specific_tests %||% TRUE)) {
        message("GSVA within-group specific: ", group_id)
        specific_rows[[specific_index]] <- tryCatch(
          run_one_gsva_specific(
            scores = group_scores,
            metadata = metadata,
            group_id = group_id,
            input_config = input_config,
            gsva_config = gsva_config,
            pathways = pathways,
            project_root = project_root,
            result_root_rel = result_root_rel,
            scoring_scope = "within_group"
          ),
          error = function(error) data.frame(
            scoring_scope = "within_group", group_id = group_id, status = "failed", reason = conditionMessage(error),
            n_cycles = NA_integer_, n_pathways = NA_integer_, n_peak_specific = NA_integer_,
            result_file = NA_character_, peak_significant_file = NA_character_, cycle_means_file = NA_character_,
            fit_file = NA_character_, stringsAsFactors = FALSE
          )
        )
        specific_index <- specific_index + 1L
      }
    }, error = function(error) {
      score_rows[[score_index]] <<- data.frame(
        scoring_scope = "within_group", group_id = group_id, status = "failed", reason = conditionMessage(error),
        n_samples = NA_integer_, n_pathways = NA_integer_, score_file = NA_character_, score_rds_file = NA_character_,
        gene_sets_file = NA_character_, input_manifest_file = NA_character_, stringsAsFactors = FALSE
      )
      score_index <<- score_index + 1L
      for (contrast_index in seq_len(nrow(group_contrasts))) {
        pair_rows[[pair_index]] <<- data.frame(
          scoring_scope = "within_group", group_id = group_id,
          comparison_id = group_contrasts$comparison_id[contrast_index],
          numerator = group_contrasts$numerator[contrast_index], denominator = group_contrasts$denominator[contrast_index],
          status = "failed", reason = paste0("分组独立打分失败：", conditionMessage(error)),
          n_numerator = NA_integer_, n_denominator = NA_integer_, n_pathways = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
          result_file = NA_character_, significant_file = NA_character_, fit_file = NA_character_, stringsAsFactors = FALSE
        )
        pair_index <<- pair_index + 1L
      }
      if (isTRUE(gsva_config$run_specific_tests %||% TRUE)) {
        specific_rows[[specific_index]] <<- data.frame(
          scoring_scope = "within_group", group_id = group_id, status = "failed",
          reason = paste0("分组独立打分失败：", conditionMessage(error)), n_cycles = NA_integer_,
          n_pathways = NA_integer_, n_peak_specific = NA_integer_, result_file = NA_character_,
          peak_significant_file = NA_character_, cycle_means_file = NA_character_, fit_file = NA_character_, stringsAsFactors = FALSE
        )
        specific_index <<- specific_index + 1L
      }
    })
  }

  score_manifest <- do.call(rbind, score_rows)
  pair_manifest <- do.call(rbind, pair_rows)
  write_tsv(score_manifest, file.path(project_root, "data", "pipeline", "gsva_groupwise_score_manifest.tsv"))
  write_tsv(pair_manifest, file.path(project_root, "data", "pipeline", "gsva_groupwise_manifest.tsv"))
  write_gsva_pairwise_summary_figure(
    pair_manifest,
    project_root,
    figure_root_rel = file.path("figures", "gsva_results", "by_group")
  )
  if (isTRUE(gsva_config$run_specific_tests %||% TRUE) && length(specific_rows)) {
    write_tsv(do.call(rbind, specific_rows), file.path(project_root, "data", "pipeline", "gsva_groupwise_specific_manifest.tsv"))
  }
  list(
    n_success = sum(pair_manifest$status == "success"),
    n_total = nrow(pair_manifest)
  )
}
