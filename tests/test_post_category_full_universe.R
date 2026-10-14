source("src/00_utils.R")
source("src/06_post_category_trends.R")

project <- tempfile("post-category-full-universe-")
dir.create(project)
on.exit(unlink(project, recursive = TRUE), add = TRUE)

methods::setClass("PostCategoryMockGseaResult", slots = c(result = "data.frame"))
write_raw <- function(comparison_id, ids) {
  directory <- file.path(project, "data", "enrichment_results", "CQ", comparison_id)
  dir.create(directory, recursive = TRUE)
  object <- methods::new("PostCategoryMockGseaResult", result = data.frame(ID = ids))
  saveRDS(object, file.path(directory, "gsea_go_bp.rds"))
}
write_raw("C1", c("A", "B", "C", "E"))
write_raw("C2", c("A", "B", "C", "D", "E"))
write_raw("C3", c("A", "B", "C", "D", "E"))
invisible(methods::removeClass("PostCategoryMockGseaResult"))

classified <- data.frame(
  Group = "CQ", Comparison = c("C1", "C1", "C2", "C3"),
  Feature_ID = c("A", "C", "B", "D"),
  Description = c("A", "C", "B", "D"),
  State = c("Positive", "Negative", "Negative", "Positive"),
  stringsAsFactors = FALSE
)
comparison_ids <- c("C1", "C2", "C3")
sig_pool_before <- post_category_build_gsea_table(classified, "CQ", comparison_ids, character())
full <- post_category_build_full_gsea_plot_table(classified, "CQ", comparison_ids, project)
sig_pool_after <- post_category_build_gsea_table(classified, "CQ", comparison_ids, character())
stopifnot(identical(sig_pool_before, sig_pool_after))
stopifnot(nrow(sig_pool_before) == 4L, nrow(full$table) == 5L)
stopifnot(identical(full$table$Feature_ID, c("A", "B", "C", "D", "E")))
stopifnot(identical(full$table$Direction_C1, c("Positive", "NS", "Negative", "NS", "NS")))
stopifnot(identical(full$table$Direction_C2, c("NS", "Negative", "NS", "NS", "NS")))
stopifnot(identical(full$table$Direction_C3, c("NS", "NS", "NS", "Positive", "NS")))
stopifnot(identical(full$counts$n_tested, c(4L, 5L, 5L)))
stopifnot(identical(full$counts$n_significant, c(2L, 1L, 1L)))
stopifnot(identical(full$counts$n_nonsig, c(3L, 4L, 4L)))
stopifnot(!"n_not_tested" %in% names(full$counts))
stopifnot(!any(grepl("NotTested", full$table$Full_Trend, fixed = TRUE)))
stopifnot(all(full$counts$n_significant + full$counts$n_nonsig == full$counts$n_union))

long <- post_category_make_sankey_long(full$table, comparison_ids)
labels <- post_category_sankey_labels(
  long, comparison_ids, c("Positive", "NS", "Negative"), "NS",
  c(Positive = "Positive", NS = "NotSig", Negative = "Negative"),
  c(Positive = "+", Negative = "-")
)
stopifnot(labels$n[labels$Comparison == "C1" & labels$State == "NS"] == 3L)
stopifnot(!"NotTested" %in% labels$State)
stopifnot(!grepl("NA", labels$label[labels$Comparison == "C2" & labels$State == "Negative"]))

wrong <- rbind(classified, data.frame(
  Group = "CQ", Comparison = "C1", Feature_ID = "Z", Description = "Z",
  State = "Positive", stringsAsFactors = FALSE
))
stopifnot(inherits(try(post_category_build_full_gsea_plot_table(
  wrong, "CQ", comparison_ids, project
), silent = TRUE), "try-error"))

cat("post-category full-universe checks passed\n")
