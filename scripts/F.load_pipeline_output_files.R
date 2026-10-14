


# =============================== 1. DESeq2 results =============================== 

read_deseq2_results <- function(
    comparison.df,
    result_dir,
    filename = "deseq2_results.tsv",
    enabled_only = TRUE
) {
  
  # 检查 comparison.df 必需字段
  required_cols <- c("group_id", "comparison_id")
  missing_cols <- setdiff(required_cols, colnames(comparison.df))
  
  if (length(missing_cols) > 0) {
    stop(
      "comparison.df 缺少必要列: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # 根据 enabled 过滤
  comp <- comparison.df
  
  if (enabled_only && "enabled" %in% colnames(comp)) {
    comp <- comp[!is.na(comp$enabled) & comp$enabled, , drop = FALSE]
  }
  
  # 初始化结果
  res <- list()
  
  for (i in seq_len(nrow(comp))) {
    
    group_id      <- comp$group_id[[i]]
    comparison_id <- comp$comparison_id[[i]]
    
    file_path <- file.path(
      result_dir,
      group_id,
      comparison_id,
      filename
    )
    
    if (!file.exists(file_path)) {
      stop(
        "找不到 DESeq2 结果文件:\n",
        file_path
      )
    }
    
    # 初始化第一层 group
    if (is.null(res[[group_id]])) {
      res[[group_id]] <- list()
    }
    
    # 第二层 comparison
    res[[group_id]][[comparison_id]] <- readr::read_tsv(
      file_path,
      show_col_types = FALSE
    ) |>
      tibble::as_tibble()
  }
  
  res
}









# =============================== 2. Trend results =============================== 

read_trend_results <- function(
    trend.df,
    result_dir,
    filename = "trend_summary.tsv",
    enabled_only = TRUE
) {
  
  # 检查 comparison.df 必需字段
  required_cols <- c("group_id", "trend_id")
  missing_cols <- setdiff(required_cols, colnames(trend.df))
  
  if (length(missing_cols) > 0) {
    stop(
      "trend.df 缺少必要列: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # 根据 enabled 过滤
  comp <- trend.df
  
  if (enabled_only && "enabled" %in% colnames(comp)) {
    comp <- comp[!is.na(comp$enabled) & comp$enabled, , drop = FALSE]
  }
  
  # 初始化结果
  res <- list()
  
  for (i in seq_len(nrow(comp))) {
    
    group_id      <- comp$group_id[[i]]
    trend_id <- comp$trend_id[[i]]
    
    file_path <- file.path(
      result_dir,
      group_id,
      trend_id,
      filename
    )
    
    if (!file.exists(file_path)) {
      stop(
        "找不到 DESeq2 结果文件:\n",
        file_path
      )
    }
    
    # 初始化第一层 group
    if (is.null(res[[group_id]])) {
      res[[group_id]] <- list()
    }
    
    # 第二层 trend
    res[[group_id]][[trend_id]] <- readr::read_tsv(
      file_path,
      show_col_types = FALSE
    ) |>
      tibble::as_tibble()
  }
  
  res
}




# =============================== 3. Enrichment significant results =============================== 


read_enrichment_results <- function(
    comparison_df,
    enrichment_dir = "./data/enrichment_results",
    gsea_file = "gsea_go_bp.tsv",
    ora_files = c(
      Up      = "ora_up.tsv",
      Down    = "ora_down.tsv",
      Up_Down = "ora_up_down.tsv"
    ),
    enabled_only = TRUE
) {
  
  # ------------------------------------------------------------
  # 0. 检查 contrasts
  # ------------------------------------------------------------
  
  required_cols <- c(
    "group_id",
    "comparison_id",
    "enabled"
  )
  
  missing_cols <- setdiff(required_cols, colnames(comparison_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "comparison_df 缺少以下列: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (enabled_only) {
    comparison_df <- comparison_df |>
      dplyr::filter(enabled)
  }
  
  # 防止重复定义
  duplicated_contrasts <- comparison_df |>
    dplyr::count(group_id, comparison_id) |>
    dplyr::filter(n > 1)
  
  if (nrow(duplicated_contrasts) > 0) {
    stop("comparison_df 中存在重复的 group_id + comparison_id")
  }
  
  
  # ------------------------------------------------------------
  # 1. 内部读取函数
  # ------------------------------------------------------------
  
  read_one_tsv <- function(group, comparison, filename) {
    
    file_path <- file.path(
      enrichment_dir,
      group,
      comparison,
      filename
    )
    
    if (!file.exists(file_path)) {
      warning("文件不存在，跳过: ", file_path)
      return(NULL)
    }
    
    readr::read_tsv(
      file_path,
      show_col_types = FALSE,
      progress = FALSE
    ) |>
      tibble::as_tibble()
  }
  
  
  # ============================================================
  # GSEA
  # ============================================================
  
  # ------------------------------------------------------------
  # 2. GSEA list
  #
  # gsea_list
  # ├── CQ
  # │   ├── R1_vs_R0
  # │   ├── R2_vs_R0
  # │   └── ...
  # └── N396
  # ------------------------------------------------------------
  
  gsea_list <- split(
    comparison_df,
    comparison_df$group_id
  ) |>
    purrr::imap(
      function(group_df, group_name) {
        
        comparisons <- group_df$comparison_id
        
        result <- purrr::map(
          comparisons,
          ~ {
            df <- read_one_tsv(
              group = group_name,
              comparison = .x,
              filename = gsea_file
            )
            
            if (is.null(df)) {
              return(NULL)
            }
            
            if (!"NES" %in% colnames(df)) {
              stop(
                "GSEA结果缺少 NES 列: ",
                group_name, "/", .x
              )
            }
            
            df |>
              dplyr::mutate(
                Direction = dplyr::if_else(
                  NES > 0,
                  "Positive",
                  "Negative"
                ),
                .after = NES
              )
          }
        )
        
        names(result) <- comparisons
        
        result
      }
    )
  
  
  # ------------------------------------------------------------
  # 3. GSEA 合并为一个 dataframe
  # ------------------------------------------------------------
  
  gsea_df <- purrr::imap_dfr(
    gsea_list,
    function(comparison_list, group_name) {
      
      purrr::imap_dfr(
        comparison_list,
        function(df, comparison_name) {
          
          if (is.null(df)) {
            return(NULL)
          }
          
          df |>
            dplyr::mutate(
              Group = group_name,
              Comparison = comparison_name,
              .before = 1
            )
        }
      )
    }
  )
  
  
  # ============================================================
  # ORA
  # ============================================================
  
  # ------------------------------------------------------------
  # 4. ORA list
  #
  # ora_list
  # ├── CQ
  # │   ├── Up
  # │   │   ├── R1_vs_R0
  # │   │   └── ...
  # │   ├── Down
  # │   └── Up_Down
  # └── N396
  # ------------------------------------------------------------
  
  ora_list <- split(
    comparison_df,
    comparison_df$group_id
  ) |>
    purrr::imap(
      function(group_df, group_name) {
        
        comparisons <- group_df$comparison_id
        
        purrr::imap(
          ora_files,
          function(filename, set_name) {
            
            result <- purrr::map(
              comparisons,
              ~ read_one_tsv(
                group = group_name,
                comparison = .x,
                filename = filename
              )
            )
            
            names(result) <- comparisons
            
            result
          }
        )
      }
    )
  
  
  # ------------------------------------------------------------
  # 5. ORA 合并为一个 dataframe
  # ------------------------------------------------------------
  
  ora_df <- purrr::imap_dfr(
    ora_list,
    function(set_list, group_name) {
      
      purrr::imap_dfr(
        set_list,
        function(comparison_list, set_name) {
          
          purrr::imap_dfr(
            comparison_list,
            function(df, comparison_name) {
              
              if (is.null(df)) {
                return(NULL)
              }
              
              df |>
                dplyr::mutate(
                  Group = group_name,
                  Comparison = comparison_name,
                  Sets = set_name,
                  .before = 1
                )
            }
          )
        }
      )
    }
  )
  
  
  # ============================================================
  # 6. 返回
  # ============================================================
  
  list(
    gsea = list(
      list = gsea_list,
      df = gsea_df
    ),
    ora = list(
      list = ora_list,
      df = ora_df
    )
  )
}






# =============================== 4. GSVA limma for contrast results =============================== 


read_gsva.limma_results <- function(
    comparison.df,
    result_dir,
    filename = "limma_results.tsv",
    enabled_only = TRUE
) {
  
  # 检查 comparison.df 必需字段
  required_cols <- c("group_id", "comparison_id")
  missing_cols <- setdiff(required_cols, colnames(comparison.df))
  
  if (length(missing_cols) > 0) {
    stop(
      "comparison.df 缺少必要列: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # 根据 enabled 过滤
  comp <- comparison.df
  
  if (enabled_only && "enabled" %in% colnames(comp)) {
    comp <- comp[!is.na(comp$enabled) & comp$enabled, , drop = FALSE]
  }
  
  # 初始化结果
  res <- list()
  
  for (i in seq_len(nrow(comp))) {
    
    group_id      <- comp$group_id[[i]]
    comparison_id <- comp$comparison_id[[i]]
    
    file_path <- file.path(
      result_dir,
      group_id,
      comparison_id,
      filename
    )
    
    if (!file.exists(file_path)) {
      stop(
        "找不到 DESeq2 结果文件:\n",
        file_path
      )
    }
    
    # 初始化第一层 group
    if (is.null(res[[group_id]])) {
      res[[group_id]] <- list()
    }
    
    # 第二层 comparison
    res[[group_id]][[comparison_id]] <- readr::read_tsv(
      file_path,
      show_col_types = FALSE
    ) |>
      tibble::as_tibble()
  }
  
  res
}





# 使用示例 --------------------------------

# 见Notebook/Check_pipeline_output.rmd










