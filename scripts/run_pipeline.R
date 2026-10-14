#!/usr/bin/env Rscript

# Configuration-driven entry point for multi-round, within-group RNA-seq analysis.

parse_cli <- function(args) {
  output <- list(config = "config.yml", stage = "all", dry_run = FALSE)
  index <- 1L
  while (index <= length(args)) {
    arg <- args[[index]]
    if (arg %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/run_pipeline.R --config config.yml [--stage all|de|enrichment|gsva] [--dry-run]\n")
      quit(status = 0L)
    }
    if (arg == "--dry-run") {
      output$dry_run <- TRUE
      index <- index + 1L
      next
    }
    if (arg %in% c("--config", "--stage")) {
      if (index == length(args)) stop(arg, " 需要一个值。", call. = FALSE)
      output[[sub("^--", "", arg)]] <- args[[index + 1L]]
      index <- index + 2L
      next
    }
    stop("未知参数：", arg, call. = FALSE)
  }
  if (identical(output$stage, "trend")) {
    stop("趋势统计已移至 GSEA 分类后：请完成分类后运行 scripts/run_post_category_trends.R。", call. = FALSE)
  }
  if (!output$stage %in% c("all", "de", "enrichment", "gsva")) {
    stop("--stage 只能是 all、de、enrichment 或 gsva。", call. = FALSE)
  }
  output
}

script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
template_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source(file.path(template_root, "src", "00_utils.R"))
source(file.path(template_root, "src", "01_deseq2.R"))
source(file.path(template_root, "src", "02_enrichment.R"))
source(file.path(template_root, "src", "03_gsva.R"))

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  require_package("yaml")
  config_file <- normalizePath(cli$config, mustWork = FALSE)
  if (!file.exists(config_file)) abort("找不到配置文件：", config_file)
  config <- yaml::read_yaml(config_file)
  # Load the configured ggplot2 before any plotting dependency (for example
  # patchwork) can resolve the default library-tree version.
  load_configured_ggplot2(config)
  if (is.null(config$input)) abort("配置文件缺少 input 节。")
  project_root <- normalizePath(
    resolve_path(dirname(config_file), config$project_root %||% "."),
    mustWork = FALSE
  )
  if (!dir.exists(project_root)) abort("project_root 不存在：", project_root)
  required_input <- c("counts_file", "metadata_file", "contrasts_file", "gene_column", "sample_column", "group_column", "cycle_column")
  missing_input <- setdiff(required_input, names(config$input))
  if (length(missing_input)) abort("input 节缺少：", paste(missing_input, collapse = ", "))

  if (cli$dry_run) {
    for (name in c("counts_file", "metadata_file", "contrasts_file")) {
      path <- resolve_path(project_root, config$input[[name]])
      if (!file.exists(path)) abort("dry-run 找不到 ", name, "：", path)
    }
    if (cli$stage %in% c("all", "enrichment") && isTRUE((config$enrichment %||% list())$enabled %||% TRUE)) {
      for (name in c("term2gene_file", "term2name_file")) {
        path <- resolve_path(project_root, config$enrichment[[name]])
        if (!file.exists(path)) abort("dry-run 找不到 ", name, "：", path)
      }
    }
    if (cli$stage == "gsva" && is.null(config$gsva)) abort("配置文件缺少 gsva 节。")
    if (cli$stage %in% c("all", "gsva") && !is.null(config$gsva) && isTRUE(config$gsva$enabled %||% TRUE)) {
      path <- resolve_path(project_root, config$gsva$tpm_file)
      if (!file.exists(path)) abort("dry-run 找不到 gsva.tpm_file：", path)
    }
    message("配置与所需输入路径检查通过；dry-run 未执行统计分析。")
    return(invisible(NULL))
  }

  de_manifest <- NULL
  if (cli$stage %in% c("all", "de")) de_manifest <- run_deseq2_stage(config, project_root)
  if (cli$stage %in% c("all", "enrichment")) run_enrichment_stage(config, project_root, de_manifest)
  if (cli$stage == "gsva" && is.null(config$gsva)) abort("配置文件缺少 gsva 节。")
  if (cli$stage %in% c("all", "gsva") && !is.null(config$gsva) && isTRUE(config$gsva$enabled %||% TRUE)) {
    run_gsva_stage(config, project_root)
  }
  message("Pipeline 运行结束。")
}

main()
