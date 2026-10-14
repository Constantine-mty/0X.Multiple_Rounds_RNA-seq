
# 表达矩阵和metadata因各个项目上游来源不一致，无法标准化，只能手动整理；
# 首先规范轮次命名，原先Non，R1，R2，R3转换为R0，R1，R2，R3...

# ================== Prepare exp counts ================

load("~/Project/11.piRNA-Discovery/2026.07-新论文路线/object/00.ZHJ_Exp_Meta.RData")

exp.count <- exp.gene$Counts

# cq.exp <- exp.comuts |>
#   dplyr::select(dplyr::starts_with("CQ_"))
# 
# cq.exp <- cq.exp[, meta_cq$Run]

exp.count <- exp.count[, meta$Run]

counts_out <- exp.count |>
  tibble::rownames_to_column("Gene")

write.table(
  counts_out,
  file = "./data/counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


exp.tpm <- exp.gene$TPM
exp.tpm <- exp.tpm[, meta$Run]


tpm_out <- exp.tpm |>
  tibble::rownames_to_column("Gene")


write.table(
  tpm_out,
  file = "./data/tpm.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)




# ================== Prepare metadata ================


meta <- meta |>
  dplyr::select(
    Sample = Run,
    Group = CAR,
    Cycle,
    Rep
  ) |>
  dplyr::mutate(
    Cycle = dplyr::recode(Cycle, "Non" = "R0")
  )



write.table(
  meta,
  file = "./data/metadata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)




# ================== Prepare Contrasts ================


cycles <- c("R0", "R1", "R2", "R3")

pairs <- combn(cycles, 2, simplify = FALSE)

contrasts <- tidyr::crossing(
  group_id = c("CQ", "N396"),
  pair = pairs
) |>
  dplyr::mutate(
    denominator = purrr::map_chr(pair, 1),
    numerator   = purrr::map_chr(pair, 2),
    comparison_id = paste0(numerator, "_vs_", denominator),
    enabled = TRUE
  ) |>
  dplyr::select(
    group_id,
    comparison_id,
    numerator,
    denominator,
    enabled
  )

write.table(
  contrasts,
  "./data/contrasts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



# ================== Prepare Trend group ================


trend_sets <- tidyr::crossing(
  group_id = c("CQ", "N396"),
  trend_id = c("adjacent")
) |>
  dplyr::mutate(
    comparison_ids = dplyr::case_when(
      trend_id == "adjacent" ~ "R1_vs_R0;R2_vs_R1;R3_vs_R2"
    )
  )


write.table(
  trend_sets,
  "./data/trend_sets.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)




# ================== Prepare Pathways ================


library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(dplyr)


packageVersion("org.Hs.eg.db")
packageVersion("GO.db")

metadata(org.Hs.eg.db)

AnnotationDbi::metadata(org.Hs.eg.db) %>%
  dplyr::filter(grepl("GO|EG|DATE", name, ignore.case = TRUE))



# GO all(+ ancestor GO terms)
gobp_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  keytype = "SYMBOL",
  columns = c("GOALL", "ONTOLOGYALL")
) %>%
  dplyr::filter(ONTOLOGYALL == "BP") %>%
  dplyr::filter(!is.na(GOALL)) %>%
  dplyr::distinct(SYMBOL, GOALL) %>%
  dplyr::rename(
    Gene = SYMBOL,
    GO_ID = GOALL
  )


# 15517
go_name_df <- AnnotationDbi::select(
  GO.db,
  keys = unique(gobp_df$GO_ID),
  keytype = "GOID",
  columns = c("TERM", "ONTOLOGY")
) %>%
  filter(ONTOLOGY == "BP") %>%
  distinct(GOID, TERM)


gobp_df <- gobp_df %>%
  left_join(
    go_name_df,
    by = c("GO_ID" = "GOID")
  )


# gobp_gene_sets <- gobp_df %>%
#   dplyr::select(TERM, Gene) %>%
#   distinct() %>%
#   split(x = .$Gene, f = .$TERM)


TERM2GENE <- gobp_df %>%
  dplyr::select(
    term = GO_ID,
    gene = Gene
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(term, gene)

TERM2NAME <- gobp_df %>%
  dplyr::select(
    term = GO_ID,
    name = TERM
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(term)

head(TERM2GENE)
head(TERM2NAME)


readr::write_tsv(
  TERM2GENE,
  "./data/GO_BP_TERM2GENE.tsv"
)

readr::write_tsv(
  TERM2NAME,
  "./data/GO_BP_TERM2NAME.tsv"
)



# not GO ALL
# gobp_df2 <- AnnotationDbi::select(
#   org.Hs.eg.db,
#   keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
#   keytype = "SYMBOL",
#   columns = c("GO", "ONTOLOGY")
# ) %>%
#   dplyr::filter(ONTOLOGY == "BP") %>%
#   dplyr::filter(!is.na(GO)) %>%
#   dplyr::distinct(SYMBOL, GO) %>%
#   dplyr::rename(
#     Gene = SYMBOL,
#     GO_ID = GO
#   )
# 
# 
# go_name_df2 <- AnnotationDbi::select(
#   GO.db,
#   keys = unique(gobp_df2$GO_ID),
#   keytype = "GOID",
#   columns = c("TERM", "ONTOLOGY")
# ) %>%
#   filter(ONTOLOGY == "BP") %>%
#   distinct(GOID, TERM)
# 
# 
# gobp_df2 <- gobp_df2 %>%
#   left_join(
#     go_name_df2,
#     by = c("GO_ID" = "GOID")
#   )




de.check <- read_tsv("./data/pipeline/de_manifest.tsv")
head(de.check)





