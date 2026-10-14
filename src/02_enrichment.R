# Enrichment stage: all ORA and GSEA analyses use the same project-supplied pathway definitions.

read_term2gene <- function(path) {
  term2gene <- read_tsv(path, c("term", "gene"))
  term2gene <- term2gene[, c("term", "gene"), drop = FALSE]
  term2gene$term <- trimws(as.character(term2gene$term))
  term2gene$gene <- trimws(as.character(term2gene$gene))
  if (anyNA(term2gene$term) || anyNA(term2gene$gene) ||
      any(!nzchar(term2gene$term)) || any(!nzchar(term2gene$gene))) {
    abort("TERM2GENE 含有空的 term 或 gene。")
  }
  term2gene <- unique(term2gene)
  if (!nrow(term2gene)) abort("TERM2GENE 为空。")
  term2gene
}

read_term2name <- function(path) {
  term2name <- read_tsv(path, c("term", "name"))
  term2name <- term2name[, c("term", "name"), drop = FALSE]
  term2name$term <- trimws(as.character(term2name$term))
  term2name$name <- trimws(as.character(term2name$name))
  if (anyNA(term2name$term) || anyNA(term2name$name) ||
      any(!nzchar(term2name$term)) || any(!nzchar(term2name$name))) {
    abort("TERM2NAME 含有空的 term 或 name。")
  }
  duplicated_terms <- duplicated(term2name$term) | duplicated(term2name$term, fromLast = TRUE)
  if (any(duplicated_terms)) {
    conflicting <- vapply(split(term2name$name[duplicated_terms], term2name$term[duplicated_terms]),
      function(x) length(unique(x)) > 1L, logical(1)
    )
    if (any(conflicting)) abort("TERM2NAME 中同一 term 对应了多个不同名称。")
  }
  term2name[!duplicated(term2name$term), , drop = FALSE]
}

load_pathway_definitions <- function(enrichment_config, project_root) {
  required <- c("term2gene_file", "term2name_file")
  missing <- required[vapply(required, function(x) is.null(enrichment_config[[x]]) || !nzchar(enrichment_config[[x]]), logical(1))]
  if (length(missing)) {
    abort("enrichment 节缺少：", paste(missing, collapse = ", "), "。")
  }
  term2gene_path <- resolve_path(project_root, enrichment_config$term2gene_file)
  term2name_path <- resolve_path(project_root, enrichment_config$term2name_file)
  term2gene <- read_term2gene(term2gene_path)
  term2name <- read_term2name(term2name_path)
  unnamed_terms <- setdiff(unique(term2gene$term), unique(term2name$term))
  if (length(unnamed_terms)) {
    warning("TERM2NAME 缺少 ", length(unnamed_terms), " 个 TERM2GENE term 的名称；结果中将保留 term ID。", call. = FALSE)
  }
  list(
    term2gene = term2gene,
    term2name = term2name,
    genes = unique(term2gene$gene),
    term2gene_path = normalizePath(term2gene_path, mustWork = TRUE),
    term2name_path = normalizePath(term2name_path, mustWork = TRUE)
  )
}

write_pathway_reference <- function(pathways, project_root) {
  reference <- data.frame(
    reference = c("TERM2GENE", "TERM2NAME"),
    file = c(pathways$term2gene_path, pathways$term2name_path),
    md5 = unname(tools::md5sum(c(pathways$term2gene_path, pathways$term2name_path))),
    n_rows = c(nrow(pathways$term2gene), nrow(pathways$term2name)),
    n_terms = c(length(unique(pathways$term2gene$term)), length(unique(pathways$term2name$term))),
    stringsAsFactors = FALSE
  )
  write_tsv(reference, file.path(project_root, "data", "pipeline", "pathway_reference.tsv"))
}

enrichment_dataframe <- function(result) {
  if (is.null(result)) return(data.frame())
  as.data.frame(result)
}

filter_gsea_dataframe <- function(result, enrichment_config) {
  table <- enrichment_dataframe(result)
  if (!nrow(table)) return(table)
  pvalue_cutoff <- as.numeric(enrichment_config$pvalue_cutoff %||% 0.05)
  # This reproduces clusterProfiler::GSEA's reporting rule while leaving the RDS unfiltered.
  table <- table[!is.na(table$pvalue) & table$pvalue <= pvalue_cutoff, , drop = FALSE]
  table <- table[!is.na(table$p.adjust) & table$p.adjust <= pvalue_cutoff, , drop = FALSE]
  table
}

run_pathway_ora <- function(genes, universe, pathways, enrichment_config) {
  genes <- intersect(unique(stats::na.omit(as.character(genes))), pathways$genes)
  universe <- intersect(unique(stats::na.omit(as.character(universe))), pathways$genes)
  min_gs_size <- as.integer(enrichment_config$min_gs_size %||% 10L)
  if (length(genes) < min_gs_size || !length(universe)) return(NULL)
  clusterProfiler::enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = pathways$term2gene,
    TERM2NAME = pathways$term2name,
    pAdjustMethod = enrichment_config$p_adjust_method %||% "BH",
    pvalueCutoff = as.numeric(enrichment_config$pvalue_cutoff %||% 0.05),
    qvalueCutoff = as.numeric(enrichment_config$qvalue_cutoff %||% 0.2),
    minGSSize = min_gs_size,
    maxGSSize = as.integer(enrichment_config$max_gs_size %||% 500L)
  )
}

run_pathway_gsea <- function(de_result, pathways, enrichment_config) {
  rank_column <- enrichment_config$rank_column %||% "stat"
  required <- c("Gene", rank_column)
  missing <- setdiff(required, colnames(de_result))
  if (length(missing)) abort("GSEA 所需的 DE 列不存在：", paste(missing, collapse = ", "))
  ranking <- data.frame(
    Gene = as.character(de_result$Gene),
    Rank = suppressWarnings(as.numeric(de_result[[rank_column]])),
    stringsAsFactors = FALSE
  )
  ranking <- ranking[!is.na(ranking$Gene) & nzchar(ranking$Gene) & is.finite(ranking$Rank), , drop = FALSE]
  ranking <- ranking[ranking$Gene %in% pathways$genes, , drop = FALSE]
  ranking <- ranking[order(-abs(ranking$Rank)), , drop = FALSE]
  ranking <- ranking[!duplicated(ranking$Gene), , drop = FALSE]
  gene_list <- ranking$Rank
  names(gene_list) <- ranking$Gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  min_gs_size <- as.integer(enrichment_config$min_gs_size %||% 10L)
  if (length(gene_list) < min_gs_size) {
    return(list(result = NULL, ranking = ranking, reason = "TERM2GENE 可匹配的排序基因数少于 min_gs_size"))
  }
  result <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = pathways$term2gene,
    TERM2NAME = pathways$term2name,
    minGSSize = min_gs_size,
    maxGSSize = as.integer(enrichment_config$max_gs_size %||% 500L),
    # Keep every tested pathway in the gseaResult RDS.  Filtering is applied only to the TSV export.
    pvalueCutoff = 1,
    pAdjustMethod = enrichment_config$p_adjust_method %||% "BH",
    verbose = FALSE,
    seed = TRUE
  )
  list(result = result, ranking = ranking, reason = NA_character_)
}

run_one_enrichment <- function(manifest_row, enrichment_config, trend_config, project_root, pathways) {
  de_result <- read_tsv(file.path(project_root, manifest_row$result_file), "Gene")
  required <- c("Gene", "log2FoldChange", "padj")
  missing <- setdiff(required, colnames(de_result))
  if (length(missing)) abort("DE 结果缺少富集所需列：", paste(missing, collapse = ", "))
  fc_cutoff <- log2(as.numeric(trend_config$fc_cutoff %||% 2))
  padj_cutoff <- as.numeric(trend_config$padj_cutoff %||% 0.05)
  direction <- ifelse(
    !is.na(de_result$padj) & de_result$padj <= padj_cutoff & de_result$log2FoldChange >= fc_cutoff,
    "Up",
    ifelse(!is.na(de_result$padj) & de_result$padj <= padj_cutoff & de_result$log2FoldChange <= -fc_cutoff, "Down", "NotSig")
  )
  universe <- unique(stats::na.omit(as.character(de_result$Gene[!is.na(de_result$padj)])))
  if (!length(universe)) universe <- unique(stats::na.omit(as.character(de_result$Gene)))
  mapped_universe <- intersect(universe, pathways$genes)
  n_up_mapped <- length(intersect(de_result$Gene[direction == "Up"], pathways$genes))
  n_down_mapped <- length(intersect(de_result$Gene[direction == "Down"], pathways$genes))
  n_de_mapped <- length(intersect(de_result$Gene[direction %in% c("Up", "Down")], pathways$genes))
  common <- data.frame(
    group_id = manifest_row$group_id,
    comparison_id = manifest_row$comparison_id,
    n_universe = length(universe),
    n_universe_mapped = length(mapped_universe),
    n_up_mapped = n_up_mapped,
    n_down_mapped = n_down_mapped,
    n_de_mapped = n_de_mapped,
    enrichment_directory = file.path("data", "enrichment_results", safe_path_component(manifest_row$group_id), safe_path_component(manifest_row$comparison_id)),
    ora_up_rds_file = NA_character_,
    ora_down_rds_file = NA_character_,
    ora_up_down_rds_file = NA_character_,
    gsea_rds_file = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!length(mapped_universe)) {
    return(cbind(common, status = "skipped_no_mapped_genes", reason = "DE 基因与 TERM2GENE 没有重叠；请检查 Gene ID 类型。"))
  }
  output_dir <- file.path(project_root, common$enrichment_directory)
  ensure_dir(output_dir)
  common$ora_up_rds_file <- file.path(common$enrichment_directory, "ora_up.rds")
  common$ora_down_rds_file <- file.path(common$enrichment_directory, "ora_down.rds")
  common$ora_up_down_rds_file <- file.path(common$enrichment_directory, "ora_up_down.rds")
  common$gsea_rds_file <- file.path(common$enrichment_directory, "gsea_go_bp.rds")
  outputs <- list(
    ora_up = run_pathway_ora(de_result$Gene[direction == "Up"], universe, pathways, enrichment_config),
    ora_down = run_pathway_ora(de_result$Gene[direction == "Down"], universe, pathways, enrichment_config),
    ora_up_down = run_pathway_ora(de_result$Gene[direction %in% c("Up", "Down")], universe, pathways, enrichment_config)
  )
  # Save the untouched objects first: their @result slot retains every tested pathway.
  write_rds(outputs$ora_up, file.path(project_root, common$ora_up_rds_file))
  write_rds(outputs$ora_down, file.path(project_root, common$ora_down_rds_file))
  write_rds(outputs$ora_up_down, file.path(project_root, common$ora_up_down_rds_file))
  for (name in names(outputs)) write_tsv(enrichment_dataframe(outputs[[name]]), file.path(output_dir, paste0(name, ".tsv")))
  gsea <- run_pathway_gsea(de_result, pathways, enrichment_config)
  write_rds(gsea$result, file.path(project_root, common$gsea_rds_file))
  write_tsv(gsea$ranking, file.path(output_dir, "gsea_ranking.tsv"))
  write_tsv(filter_gsea_dataframe(gsea$result, enrichment_config), file.path(output_dir, "gsea_go_bp.tsv"))
  cbind(common, status = "success", reason = gsea$reason)
}

run_enrichment_stage <- function(config, project_root, manifest = NULL) {
  enrichment_config <- config$enrichment %||% list()
  if (!isTRUE(enrichment_config$enabled %||% TRUE)) {
    message("富集分析已在配置中禁用。")
    return(invisible(NULL))
  }
  require_package("clusterProfiler")
  manifest <- manifest %||% read_manifest(project_root)
  successful <- manifest[manifest$status == "success", , drop = FALSE]
  if (!nrow(successful)) abort("没有成功的 DE 比较可用于富集分析。")
  pathways <- load_pathway_definitions(enrichment_config, project_root)
  write_pathway_reference(pathways, project_root)
  rows <- lapply(seq_len(nrow(successful)), function(i) {
    message("Enrichment: ", successful$group_id[i], " / ", successful$comparison_id[i])
    tryCatch(
      run_one_enrichment(
        successful[i, , drop = FALSE], enrichment_config, config$trend %||% list(), project_root, pathways
      ),
      error = function(e) data.frame(
        group_id = successful$group_id[i], comparison_id = successful$comparison_id[i],
        n_universe = NA_integer_, n_universe_mapped = NA_integer_, n_up_mapped = NA_integer_,
        n_down_mapped = NA_integer_, n_de_mapped = NA_integer_,
        enrichment_directory = NA_character_, ora_up_rds_file = NA_character_, ora_down_rds_file = NA_character_,
        ora_up_down_rds_file = NA_character_, gsea_rds_file = NA_character_,
        status = "failed", reason = conditionMessage(e), stringsAsFactors = FALSE
      )
    )
  })
  output <- do.call(rbind, rows)
  write_tsv(output, file.path(project_root, "data", "pipeline", "enrichment_manifest.tsv"))
  message("富集分析完成：", sum(output$status == "success"), "/", nrow(output), " 个比较成功。")
  output
}
