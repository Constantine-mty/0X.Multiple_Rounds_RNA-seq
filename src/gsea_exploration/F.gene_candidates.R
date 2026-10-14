# Candidate genes from significant, classified GSEA leading-edge genes.
# Source F.common.R before this file.

gm_extract_core <- function(gsea, group, comparison, category, direction,
                            category_column = "Category_Merged") {
  gsea_nb_require(gsea, c("Group", "Comparison", "Direction", "ID", "Description",
                          "core_enrichment", category_column), "分类 GSEA")
  direction <- gsea_nb_direction(direction)
  if (length(direction) != 1L || is.na(direction)) stop("请指定单个方向。", call. = FALSE)
  rows <- gsea[as.character(gsea$Group) == group &
                 as.character(gsea$Comparison) == comparison &
                 as.character(gsea$Direction) == direction &
                 !is.na(gsea[[category_column]]) &
                 as.character(gsea[[category_column]]) %in% category, , drop = FALSE]
  if (!nrow(rows)) stop("指定的 Group × Comparison × Category × Direction 没有显著通路。", call. = FALSE)
  rows <- rows[!duplicated(as.character(rows$ID)), , drop = FALSE]
  parsed <- lapply(rows$core_enrichment, gsea_nb_core_genes)
  links <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    if (!length(parsed[[i]])) return(NULL)
    data.frame(Gene = parsed[[i]], ID = as.character(rows$ID[[i]]),
               Description = as.character(rows$Description[[i]]), stringsAsFactors = FALSE)
  }))
  if (is.null(links) || !nrow(links)) stop("所选通路均没有可解析的 core_enrichment 基因。", call. = FALSE)
  links <- unique(links)
  counts <- sort(table(links$Gene), decreasing = TRUE)
  candidates <- data.frame(Gene = names(counts), n_core_pathways = as.integer(counts),
                           core_pathway_fraction = as.numeric(counts) / nrow(rows),
                           n_selected_pathways = nrow(rows), stringsAsFactors = FALSE)
  candidates$pathway_ids <- vapply(candidates$Gene, function(gene) {
    paste(links$ID[links$Gene == gene], collapse = ";")
  }, character(1))
  list(pathways = rows, core_links = links, candidates = candidates,
       selection = list(group = group, comparison = comparison,
                        category = category, category_column = category_column,
                        direction = direction))
}

gm_add_metrics <- function(extraction, de, tpm, metadata, contrasts,
                           gene_column = "Gene", sample_column = "Sample",
                           group_column = "Group", cycle_column = "Cycle") {
  gsea_nb_require(de, c(gene_column, "log2FoldChange", "padj"), "DESeq2")
  gsea_nb_require(metadata, c(sample_column, group_column, cycle_column), "metadata")
  gsea_nb_require(contrasts, c("group_id", "comparison_id", "numerator", "denominator"), "contrasts")
  sel <- extraction$selection
  match_contrast <- as.character(contrasts$group_id) == sel$group &
    as.character(contrasts$comparison_id) == sel$comparison
  if (sum(match_contrast) != 1L) stop("contrasts 中必须唯一匹配当前 Group + Comparison。", call. = FALSE)
  contrast <- contrasts[match_contrast, , drop = FALSE]
  numerator <- as.character(contrast$numerator[[1]])
  denominator <- as.character(contrast$denominator[[1]])
  group_meta <- metadata[as.character(metadata[[group_column]]) == sel$group, , drop = FALSE]
  num_samples <- as.character(group_meta[[sample_column]][as.character(group_meta[[cycle_column]]) == numerator])
  den_samples <- as.character(group_meta[[sample_column]][as.character(group_meta[[cycle_column]]) == denominator])
  if (!length(num_samples) || !length(den_samples)) stop("metadata 中缺少比较对应的轮次样本。", call. = FALSE)
  if (anyDuplicated(metadata[[sample_column]])) stop("metadata Sample 重复。", call. = FALSE)
  tpm <- as.matrix(tpm)
  if (is.null(rownames(tpm)) || is.null(colnames(tpm)) || anyDuplicated(rownames(tpm))) {
    stop("TPM 必须具有唯一的基因行名和样本列名。", call. = FALSE)
  }
  missing_samples <- setdiff(c(num_samples, den_samples), colnames(tpm))
  if (length(missing_samples)) stop("TPM 缺少样本：", paste(missing_samples, collapse = ", "), call. = FALSE)
  if (anyDuplicated(as.character(de[[gene_column]]))) stop("DESeq2 结果的 Gene 重复。", call. = FALSE)
  out <- extraction$candidates
  de_index <- match(out$Gene, as.character(de[[gene_column]]))
  out$log2FoldChange <- as.numeric(de$log2FoldChange[de_index])
  out$padj <- as.numeric(de$padj[de_index])
  if ("baseMean" %in% names(de)) out$baseMean <- as.numeric(de$baseMean[de_index])
  tpm_index <- match(out$Gene, rownames(tpm))
  column_mean <- function(samples) {
    value <- rep(NA_real_, nrow(out))
    found <- which(!is.na(tpm_index))
    if (length(found)) value[found] <- rowMeans(tpm[tpm_index[found], samples, drop = FALSE])
    value
  }
  out$numerator_mean_tpm <- column_mean(num_samples)
  out$denominator_mean_tpm <- column_mean(den_samples)
  out$pair_mean_tpm <- column_mean(c(num_samples, den_samples))
  high_samples <- if (sel$direction == "Positive") num_samples else den_samples
  out$high_condition_mean_tpm <- column_mean(high_samples)
  out$high_condition_detected_fraction <- NA_real_
  found <- which(!is.na(tpm_index))
  if (length(found)) {
    out$high_condition_detected_fraction[found] <- rowMeans(
      tpm[tpm_index[found], high_samples, drop = FALSE] > 1
    )
  }
  out$matches_gsea_direction <- !is.na(out$log2FoldChange) &
    if (sel$direction == "Positive") out$log2FoldChange > 0 else out$log2FoldChange < 0
  out$Group <- sel$group
  out$Comparison <- sel$comparison
  out$Direction <- sel$direction
  out$Category <- paste(sel$category, collapse = ";")
  out
}

gm_rank_candidates <- function(metrics, min_mean_tpm = 1,
                               rank_by = c("frequency", "log2fc"),
                               expression_metric = "high_condition_mean_tpm") {
  rank_by <- match.arg(rank_by)
  if (length(expression_metric) != 1L || is.na(expression_metric) || !nzchar(expression_metric)) {
    stop("expression_metric 必须是单个表达量列名。", call. = FALSE)
  }
  gsea_nb_require(metrics, c("Gene", "Direction", "n_core_pathways",
                             "log2FoldChange", expression_metric), "候选基因表")
  if (length(min_mean_tpm) != 1L || !is.finite(min_mean_tpm) || min_mean_tpm < 0) {
    stop("min_mean_tpm 必须是非负的有限数值。", call. = FALSE)
  }
  out <- as.data.frame(metrics, stringsAsFactors = FALSE)
  direction <- gsea_nb_direction(out$Direction)
  fc <- as.numeric(out$log2FoldChange)
  expression <- as.numeric(out[[expression_metric]])
  frequency <- as.numeric(out$n_core_pathways)
  direction_ok <- is.finite(fc) &
    ((direction == "Positive" & fc > 0) | (direction == "Negative" & fc < 0))
  keep <- !is.na(direction_ok) & direction_ok & is.finite(expression) &
    expression >= min_mean_tpm & is.finite(frequency) & frequency > 0
  out <- out[keep, , drop = FALSE]
  if (rank_by == "frequency") {
    order_index <- order(-frequency[keep], -abs(fc[keep]),
                         -expression[keep], as.character(out$Gene))
  } else {
    order_index <- order(-abs(fc[keep]), -frequency[keep],
                         -expression[keep], as.character(out$Gene))
  }
  out <- out[order_index, , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}
