#!/usr/bin/env Rscript

# Re-draw DE, trend and GSEA summary figures without re-running statistics.

parse_result_plot_cli <- function(args) {
  output <- list(config = "config.yml")
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("-h", "--help")) {
      cat("Usage: Rscript scripts/run_result_plots.R --config config.yml\n")
      quit(status = 0L)
    }
    if (argument == "--config") {
      if (index == length(args)) stop("--config 需要一个值。", call. = FALSE)
      output$config <- args[[index + 1L]]
      index <- index + 2L
      next
    }
    stop("未知参数：", argument, call. = FALSE)
  }
  output
}

script_argument <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_argument[grep("^--file=", script_argument)][1])
template_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source(file.path(template_root, "src", "00_utils.R"))
source(file.path(template_root, "src", "08_result_plots.R"))

main <- function() {
  cli <- parse_result_plot_cli(commandArgs(trailingOnly = TRUE))
  require_package("yaml")
  config_file <- normalizePath(cli$config, mustWork = FALSE)
  if (!file.exists(config_file)) abort("找不到配置文件：", config_file)
  config <- yaml::read_yaml(config_file)
  load_configured_ggplot2(config)
  project_root <- normalizePath(
    resolve_path(dirname(config_file), config$project_root %||% "."),
    mustWork = FALSE
  )
  if (!dir.exists(project_root)) abort("project_root 不存在：", project_root)
  run_result_plot_stage(config, project_root)
}

main()
