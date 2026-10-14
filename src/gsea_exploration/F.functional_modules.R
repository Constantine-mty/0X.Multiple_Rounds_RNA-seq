# Multi-resolution semantic clustering of classified GO BP GSEA terms.
# Source F.common.R before this file. No analysis runs on source().

gmc_selection_keys <- function(gsea,
                               mode = c("group_comparison_direction", "group_comparison"),
                               groups = NULL, comparisons = NULL, directions = NULL) {
  mode <- match.arg(mode)
  gsea_nb_require(gsea, c("Group", "Comparison", "Direction"), "分类 GSEA")
  if (mode == "group_comparison" && !is.null(directions)) {
    stop("group_comparison 模式合并两个方向；请将 directions 设为 NULL。", call. = FALSE)
  }
  rows <- as.data.frame(gsea[, c("Group", "Comparison", "Direction"), drop = FALSE])
  rows$Group <- as.character(rows$Group)
  rows$Comparison <- as.character(rows$Comparison)
  rows$Direction <- gsea_nb_direction(rows$Direction)
  rows <- rows[!is.na(rows$Group) & nzchar(rows$Group) &
                 !is.na(rows$Comparison) & nzchar(rows$Comparison) &
                 !is.na(rows$Direction), , drop = FALSE]
  if (!is.null(groups)) rows <- rows[rows$Group %in% groups, , drop = FALSE]
  if (!is.null(comparisons)) rows <- rows[rows$Comparison %in% comparisons, , drop = FALSE]
  if (!is.null(directions)) rows <- rows[rows$Direction %in% gsea_nb_direction(directions), , drop = FALSE]
  if (!nrow(rows)) stop("没有符合选择条件的 GSEA 分支。", call. = FALSE)
  if (mode == "group_comparison_direction") {
    keys <- unique(rows)
    names(keys)[names(keys) == "Direction"] <- "Analysis_Direction"
  } else {
    keys <- unique(rows[, c("Group", "Comparison"), drop = FALSE])
    keys$Analysis_Direction <- "Combined"
  }
  keys$key <- make.unique(paste(keys$Group, keys$Comparison,
                                keys$Analysis_Direction, sep = "__"))
  rownames(keys) <- NULL
  keys
}

gmc_prepare_terms <- function(gsea, group, comparison, direction = NULL, max_terms = Inf,
                              categories = NULL, category_column = "Category_Merged",
                              include_background = TRUE) {
  gsea_nb_require(gsea, c("Group", "Comparison", "Direction", "ID", "Description", "p.adjust", "NES", category_column), "分类 GSEA")
  if (length(group) != 1L || is.na(group) || !nzchar(group)) stop("group 必须是单个非空名称。", call. = FALSE)
  if (length(comparison) != 1L || is.na(comparison) || !nzchar(comparison)) stop("comparison 必须是单个非空名称。", call. = FALSE)
  if (!is.null(direction)) {
    direction <- gsea_nb_direction(direction)
    if (length(direction) != 1L || is.na(direction)) stop("direction 必须是单个有效方向。", call. = FALSE)
  }
  if (!is.infinite(max_terms) && (!is.finite(max_terms) || max_terms < 1)) stop("max_terms 至少为 1。", call. = FALSE)
  keep <- !is.na(gsea$Group) & as.character(gsea$Group) == group &
    !is.na(gsea$Comparison) & as.character(gsea$Comparison) == comparison
  if (!is.null(direction)) keep <- keep & !is.na(gsea$Direction) &
    gsea_nb_direction(gsea$Direction) == direction
  rows <- as.data.frame(gsea[keep, , drop = FALSE])
  rows$Direction <- gsea_nb_direction(rows$Direction)
  if (!is.null(categories)) rows <- rows[as.character(rows[[category_column]]) %in% categories, , drop = FALSE]
  if (!include_background) {
    gsea_nb_require(rows, "Background", "分类 GSEA")
    rows <- rows[!tolower(as.character(rows$Background)) %in% c("true", "t", "1"), , drop = FALSE]
  }
  rows <- rows[!is.na(rows$ID) & grepl("^GO:[0-9]{7}$", as.character(rows$ID)), , drop = FALSE]
  if (!nrow(rows)) stop("所选 Group × Comparison × Direction 范围没有可用的 GO ID。", call. = FALSE)
  rows$.padj_order <- suppressWarnings(as.numeric(rows$p.adjust))
  rows$.padj_order[!is.finite(rows$.padj_order)] <- Inf
  rows$.nes_order <- abs(suppressWarnings(as.numeric(rows$NES)))
  rows$.nes_order[!is.finite(rows$.nes_order)] <- -Inf
  rows <- rows[order(rows$.padj_order, -rows$.nes_order, as.character(rows$ID)), , drop = FALSE]
  counts <- table(as.character(rows$ID))
  sources <- vapply(split(as.character(rows$Direction), as.character(rows$ID)),
                    function(x) paste(intersect(c("Positive", "Negative"), unique(x)),
                                      collapse = ";"), character(1))
  direction_counts <- vapply(split(as.character(rows$Direction), as.character(rows$ID)),
                             function(x) length(unique(x)), integer(1))
  rows <- rows[!duplicated(as.character(rows$ID)), , drop = FALSE]
  rows$Analysis_Direction <- if (is.null(direction)) "Combined" else direction
  rows$Source_Directions <- unname(sources[as.character(rows$ID)])
  rows$n_enrichment_rows <- as.integer(counts[as.character(rows$ID)])
  rows$n_directions <- unname(direction_counts[as.character(rows$ID)])
  if (is.finite(max_terms)) rows <- head(rows, as.integer(max_terms))
  rows$.padj_order <- NULL
  rows$.nes_order <- NULL
  rownames(rows) <- NULL
  rows
}

gmc_semantic_similarity <- function(terms, orgdb = "org.Hs.eg.db", ont = "BP",
                                    method = "Rel", semdata = NULL) {
  gsea_nb_require(terms, "ID", "GO terms")
  ids <- unique(as.character(terms$ID))
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) stop("GO ID 无效或重复。", call. = FALSE)
  # Loading Rcpp first registers the C-callables used by GOSemSim's IC methods.
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("语义相似度计算需要 Rcpp 包。", call. = FALSE)
  if (is.null(semdata)) semdata <- GOSemSim::godata(annoDb = orgdb, ont = ont, keytype = "ENTREZID")
  raw <- GOSemSim::goSim(ids, ids, semData = semdata, measure = method)
  sim <- matrix(as.numeric(raw), nrow = length(ids), ncol = length(ids),
                dimnames = list(ids, ids))
  # Unsupported GO terms are excluded, never silently treated as similar.
  valid <- is.finite(diag(sim))
  omitted <- data.frame(ID = ids[!valid],
                        reason = rep("no_semantic_information", sum(!valid)))
  sim <- sim[valid, valid, drop = FALSE]
  if (!nrow(sim)) stop("所有 GO term 都缺少语义信息；请检查 OrgDb 和 ontology。", call. = FALSE)
  sim <- (sim + t(sim)) / 2
  missing_pairs <- sum(!is.finite(sim[upper.tri(sim)]))
  if (missing_pairs) {
    # Undefined cross-branch similarities are maximally distant; report count.
    sim[!is.finite(sim)] <- 0
  }
  sim[sim < 0] <- 0
  sim[sim > 1] <- 1
  diag(sim) <- 1
  list(similarity = sim, terms = terms[match(rownames(sim), as.character(terms$ID)), , drop = FALSE],
       omitted = omitted, n_missing_pairs_as_zero = missing_pairs,
       settings = list(orgdb = orgdb, ont = ont, method = method))
}

gmc_cluster_multi <- function(similarity, terms, distance_cuts = c(0.9, 0.7, 0.5),
                              linkage = "complete") {
  gsea_nb_require(terms, c("Group", "Comparison", "Analysis_Direction",
                           "ID", "Description", "p.adjust"), "GO terms")
  sim <- as.matrix(similarity)
  if (!is.numeric(sim) || !nrow(sim) || nrow(sim) != ncol(sim) ||
      is.null(rownames(sim)) || !identical(rownames(sim), colnames(sim)) ||
      any(!is.finite(sim)) || any(sim < 0 | sim > 1) ||
      any(abs(sim - t(sim)) > 1e-8)) stop("similarity 需要 0–1 的对称具名矩阵。", call. = FALSE)
  if (!length(distance_cuts) || any(!is.finite(distance_cuts)) ||
      any(distance_cuts < 0 | distance_cuts > 1) || anyDuplicated(distance_cuts)) {
    stop("distance_cuts 必须是互不重复的 0–1 距离阈值。", call. = FALSE)
  }
  if (!all(rownames(sim) %in% as.character(terms$ID))) stop("terms 缺少相似度矩阵中的 GO ID。", call. = FALSE)
  terms <- as.data.frame(terms[match(rownames(sim), as.character(terms$ID)), , drop = FALSE])
  if (anyDuplicated(as.character(terms$ID)) ||
      nrow(unique(terms[, c("Group", "Comparison", "Analysis_Direction"), drop = FALSE])) != 1L) {
    stop("一次聚类只能包含一个 Group × Comparison × Analysis_Direction 分支，且 GO ID 不能重复。", call. = FALSE)
  }
  scores <- -log10(pmax(as.numeric(terms$p.adjust), .Machine$double.xmin))
  scores[!is.finite(scores)] <- 0
  distance <- 1 - sim
  distance[distance < 0] <- 0
  tree <- if (nrow(sim) > 1L) stats::hclust(stats::as.dist(distance), method = linkage) else NULL
  cuts <- list()
  summaries <- list()
  for (i in seq_along(distance_cuts)) {
    cut <- distance_cuts[[i]]
    label <- paste0("d", format(cut, scientific = FALSE, trim = TRUE))
    cluster <- if (is.null(tree)) stats::setNames(1L, rownames(sim)) else stats::cutree(tree, h = cut)
    # Re-number clusters by representative significance for stable heatmap order.
    members <- split(seq_along(cluster), cluster)
    reps <- vapply(members, function(ix) ix[order(-scores[ix], rownames(sim)[ix])[[1]]], integer(1))
    score_by_cluster <- stats::setNames(scores[reps], names(reps))
    cluster_order <- names(sort(score_by_cluster, decreasing = TRUE))
    ordered_id <- match(as.character(cluster), cluster_order)
    rep_idx <- reps[as.character(cluster)]
    within_uniqueness <- vapply(seq_along(cluster), function(j) {
      peers <- which(cluster == cluster[[j]] & seq_along(cluster) != j)
      if (!length(peers)) return(1)
      1 - mean(sim[j, peers])
    }, numeric(1))
    assignment <- data.frame(Group = as.character(terms$Group),
                             Comparison = as.character(terms$Comparison),
                             Analysis_Direction = as.character(terms$Analysis_Direction),
                             distance_cut = cut, cluster = ordered_id,
                             ID = rownames(sim), Description = as.character(terms$Description),
                             parent = rownames(sim)[rep_idx],
                             parentTerm = as.character(terms$Description[rep_idx]),
                             score = scores,
                             termUniqueness = 1 - (rowSums(sim) - 1) / pmax(1, nrow(sim) - 1),
                             termUniquenessWithinCluster = within_uniqueness,
                             termDispensability = ifelse(seq_along(cluster) == rep_idx, 0,
                                                         sim[cbind(seq_along(cluster), rep_idx)]),
                             stringsAsFactors = FALSE)
    for (field in intersect(c("Direction", "Source_Directions", "n_enrichment_rows",
                              "n_directions", "Category", "Category_Merged",
                              "Background", "NES", "p.adjust"), names(terms))) {
      assignment[[field]] <- terms[[field]]
    }
    cuts[[label]] <- assignment
    summary <- aggregate(ID ~ Group + Comparison + Analysis_Direction +
                           distance_cut + cluster + parent + parentTerm,
                         data = assignment, FUN = length)
    names(summary)[names(summary) == "ID"] <- "n_terms"
    summary$representative_padj <- terms$p.adjust[match(summary$parent, as.character(terms$ID))]
    summaries[[label]] <- summary[order(summary$cluster), , drop = FALSE]
  }
  assignments <- do.call(rbind, cuts)
  cluster_summary <- do.call(rbind, summaries)
  rownames(assignments) <- NULL
  rownames(cluster_summary) <- NULL
  list(terms = terms, similarity = sim, tree = tree,
       assignments = assignments, cluster_summary = cluster_summary,
       by_cut = cuts, distance_cuts = distance_cuts)
}

gmc_get_cut <- function(result, distance_cut) {
  name <- paste0("d", format(distance_cut, scientific = FALSE, trim = TRUE))
  if (!name %in% names(result$by_cut)) stop("结果中不存在 distance_cut：", distance_cut, call. = FALSE)
  result$by_cut[[name]]
}

gmc_plot_heatmap <- function(result, distance_cut, max_plot_terms = 120,
                             show_labels = TRUE, palette = "blue_white") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("热图需要 pheatmap 包。", call. = FALSE)
  assignment <- gmc_get_cut(result, distance_cut)
  if (!is.finite(max_plot_terms) || max_plot_terms < 1) stop("max_plot_terms 至少为 1。", call. = FALSE)
  assignment <- assignment[order(assignment$cluster, -assignment$score, assignment$ID), , drop = FALSE]
  # Include one representative from each module before filling with extra terms.
  module_rows <- split(seq_len(nrow(assignment)), factor(assignment$cluster,
                           levels = unique(assignment$cluster)))
  round_robin <- unlist(lapply(seq_len(max(lengths(module_rows))), function(position) {
    vapply(module_rows, function(indices) if (length(indices) >= position) indices[[position]] else NA_integer_,
           integer(1))
  }), use.names = FALSE)
  selected <- head(round_robin[!is.na(round_robin)], as.integer(max_plot_terms))
  assignment <- assignment[sort(selected), , drop = FALSE]
  ids <- assignment$ID
  sim <- result$similarity[ids, ids, drop = FALSE]
  label <- if (show_labels) paste0(ids, "  ", substr(assignment$Description, 1, 55)) else ids
  gaps <- which(diff(assignment$cluster) != 0)
  clusters <- factor(assignment$cluster, levels = unique(assignment$cluster))
  annotation <- data.frame(Module = clusters, row.names = ids)
  if (is.null(palette)) palette <- "blue_white" # 兼容旧调用
  if (length(palette) == 1L) {
    palette <- switch(palette,
      blue_white = grDevices::colorRampPalette(c("#FFFFFF", "#DCE9F2", "#5B8DB8", "#173A63"))(100),
      viridis = {
        if (!requireNamespace("viridisLite", quietly = TRUE))
          stop("palette='viridis' 需要 viridisLite 包。", call. = FALSE)
        viridisLite::viridis(100)
      },
      stop("palette 请选择 'blue_white'、'viridis'，或提供至少两种颜色的向量。", call. = FALSE)
    )
  }
  if (length(palette) < 2L || anyNA(palette)) stop("palette 需要至少两种有效颜色。", call. = FALSE)
  tryCatch(grDevices::col2rgb(palette), error = function(e)
    stop("palette 包含无效颜色：", conditionMessage(e), call. = FALSE))
  scope <- paste(unique(assignment$Group), unique(assignment$Comparison),
                 unique(assignment$Analysis_Direction), sep = " / ")
  plot <- pheatmap::pheatmap(sim, color = palette, breaks = seq(0, 1, length.out = length(palette) + 1),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            labels_row = label, labels_col = ids,
                            annotation_row = annotation, annotation_col = annotation,
                            gaps_row = gaps, gaps_col = gaps,
                            show_colnames = nrow(sim) <= 40L,
                            show_rownames = show_labels, fontsize_row = 6,
                            fontsize_col = 6, border_color = NA, silent = TRUE,
                            main = paste0(scope, " | Semantic GO similarity | distance cut ", distance_cut,
                                          " | shown ", nrow(sim), "/", nrow(result$similarity), " GO terms"))
  attr(plot, "plotted_ids") <- ids
  attr(plot, "palette") <- palette
  plot
}

gmc_plot_scatter <- function(result, distance_cut, label_representatives = TRUE) {
  assignment <- gmc_get_cut(result, distance_cut)
  if (nrow(assignment) < 2L) stop("散点图至少需要两个 GO term。", call. = FALSE)
  distance <- 1 - result$similarity
  distance[distance < 0] <- 0
  coords <- stats::cmdscale(stats::as.dist(distance), k = min(2L, nrow(assignment) - 1L))
  coords <- as.matrix(coords)
  plot_df <- assignment
  plot_df$x <- coords[, 1]
  plot_df$y <- if (ncol(coords) >= 2L) coords[, 2] else 0
  plot_df$Module <- factor(plot_df$cluster)
  plot_df$label <- ifelse(label_representatives & plot_df$ID == plot_df$parent,
                          substr(plot_df$Description, 1, 35), "")
  scope <- paste(unique(plot_df$Group), unique(plot_df$Comparison),
                 unique(plot_df$Analysis_Direction), sep = " / ")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = Module, size = score)) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(title = paste0(scope, " | GO functional modules | distance cut ", distance_cut),
                  subtitle = paste(nrow(plot_df), "GO terms; representative labels"),
                  x = "PCoA 1", y = "PCoA 2", size = "−log10 adjusted P") +
    ggplot2::theme_bw(base_size = 11)
  if (label_representatives) {
    if (requireNamespace("ggrepel", quietly = TRUE))
      p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = label), max.overlaps = Inf, size = 2.8, show.legend = FALSE)
    else p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = 2.8, check_overlap = TRUE, show.legend = FALSE)
  }
  p
}
