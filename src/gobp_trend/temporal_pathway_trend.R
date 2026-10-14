# Temporal pathway/category trend analysis
# Dependencies: dplyr, tidyr, tibble, stringr, rlang
#
# Do not attach these packages.  This module is often sourced after
# Bioconductor packages (for example S4Vectors) that mask dplyr::first(),
# dplyr::select(), and related functions.  All non-base calls below therefore
# use explicit namespaces.
for (package in c("dplyr", "tidyr", "tibble", "stringr", "rlang")) {
  if (!requireNamespace(package, quietly=TRUE)) {
    stop("Temporal trend analysis 缺少 R 包 '", package, "'。", call.=FALSE)
  }
}
.data <- rlang::.data

`%||%` <- function(x, y) if (is.null(x)) y else x

prepare_annotation <- function(annotation, group = NULL,
                               group_col = "Group", id_col = "ID",
                               description_col = "Description",
                               category_col = "Category",
                               extra_cols = c("Category_CN", "Category_Merged")) {
  req <- c(id_col, description_col, category_col)
  if (!is.null(group)) req <- c(group_col, req)
  miss <- setdiff(req, colnames(annotation)); if(length(miss)) stop("Missing annotation columns: ", paste(miss, collapse=", "))
  x <- annotation
  if (!is.null(group)) x <- x |> dplyr::filter(.data[[group_col]] == group)
  keep <- unique(c(if(!is.null(group)) group_col, id_col, description_col, category_col, intersect(extra_cols, names(x))))
  x <- x |> dplyr::select(dplyr::all_of(keep)) |> dplyr::distinct()
  conflicts <- x |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(id_col, category_col)))) |>
    dplyr::count(.data[[id_col]], name="n_category") |>
    dplyr::filter(.data$n_category > 1)
  if(nrow(conflicts)) warning(nrow(conflicts), " pathway IDs map to >1 category. Resolve before interpretation.")
  x
}

extract_score_subset <- function(score_matrix, annotation, id_col="ID") {
  if(is.null(rownames(score_matrix))) stop("score_matrix must have pathway IDs as row names.")
  ids <- unique(annotation[[id_col]])
  matched <- intersect(ids, rownames(score_matrix)); missing <- setdiff(ids, rownames(score_matrix))
  anno <- annotation |>
    dplyr::filter(.data[[id_col]] %in% matched) |>
    dplyr::distinct(.data[[id_col]], .keep_all=TRUE) |>
    dplyr::arrange(match(.data[[id_col]], matched))
  mat <- score_matrix[anno[[id_col]], , drop=FALSE]
  list(score=mat, annotation=anno, missing_ids=missing,
       matching=data.frame(n_selected=length(ids), n_matched=length(matched), n_missing=length(missing)))
}

make_readable_score <- function(x, description_col="Description") {
  m <- x$score; rownames(m) <- make.unique(x$annotation[[description_col]]); m
}

prepare_time_from_metadata <- function(score_matrix, metadata, group_name,
                                       sample_column="Sample", group_column="Group",
                                       cycle_column="Cycle", cycle_order=NULL) {
  # Select and order score columns with the project metadata, rather than parsing
  # time points from sample names.  This also makes the function safe to use with
  # a whole-cohort score matrix: only samples from group_name are retained.
  if (is.null(colnames(score_matrix)) || any(!nzchar(colnames(score_matrix)))) {
    stop("score_matrix 必须具有非空的样本列名。", call. = FALSE)
  }
  needed <- c(sample_column, group_column, cycle_column)
  missing <- setdiff(needed, colnames(metadata))
  if (length(missing)) {
    stop("metadata 缺少必要列：", paste(missing, collapse=", "), call. = FALSE)
  }

  meta <- as.data.frame(metadata, stringsAsFactors=FALSE)
  for (column in needed) {
    meta[[column]] <- trimws(as.character(meta[[column]]))
    if (anyNA(meta[[column]]) || any(!nzchar(meta[[column]]))) {
      stop("metadata 的 ", column, " 列含有空值。", call. = FALSE)
    }
  }
  if (anyDuplicated(meta[[sample_column]])) {
    duplicates <- unique(meta[[sample_column]][duplicated(meta[[sample_column]])])
    stop("metadata 中存在重复样本：", paste(duplicates, collapse=", "), call. = FALSE)
  }

  group_meta <- meta[meta[[group_column]] == as.character(group_name), , drop=FALSE]
  if (!nrow(group_meta)) {
    stop("metadata 中不存在 Group='", group_name, "'。", call. = FALSE)
  }
  available <- group_meta[[sample_column]] %in% colnames(score_matrix)
  if (!any(available)) {
    stop("Group '", group_name, "' 的 metadata 样本均不在 score_matrix 列名中。", call. = FALSE)
  }
  missing_scores <- group_meta[[sample_column]][!available]
  if (length(missing_scores)) {
    warning(
      "Group '", group_name, "' 有 ", length(missing_scores),
      " 个 metadata 样本不在 score_matrix 中，将不参与趋势分析：",
      paste(missing_scores, collapse=", "), call.=FALSE
    )
  }
  unmatched_scores <- setdiff(colnames(score_matrix), meta[[sample_column]])
  if (length(unmatched_scores)) {
    warning(
      "score_matrix 中有 ", length(unmatched_scores),
      " 个样本未见于 metadata，已忽略：", paste(unmatched_scores, collapse=", "),
      call.=FALSE
    )
  }

  group_meta <- group_meta[available, , drop=FALSE]
  observed_cycles <- unique(group_meta[[cycle_column]])
  if (is.null(cycle_order) || !length(cycle_order)) {
    # The metadata row order is the default chronological order.  Projects with
    # another desired order should set cycle_order explicitly in config.yml.
    time_levels <- observed_cycles
  } else {
    cycle_order <- trimws(as.character(unlist(cycle_order, use.names=FALSE)))
    if (anyNA(cycle_order) || any(!nzchar(cycle_order)) || anyDuplicated(cycle_order)) {
      stop("cycle_order 必须是无重复、非空的轮次名称向量。", call. = FALSE)
    }
    absent <- setdiff(observed_cycles, cycle_order)
    if (length(absent)) {
      stop("cycle_order 缺少 metadata 中实际出现的轮次：", paste(absent, collapse=", "), call. = FALSE)
    }
    time_levels <- cycle_order[cycle_order %in% observed_cycles]
  }

  group_meta$.temporal_time <- factor(group_meta[[cycle_column]], levels=time_levels, ordered=TRUE)
  group_meta <- group_meta[order(group_meta$.temporal_time, group_meta[[sample_column]]), , drop=FALSE]
  sample_time_map <- stats::setNames(as.character(group_meta[[cycle_column]]), group_meta[[sample_column]])

  list(
    score_matrix=score_matrix[, names(sample_time_map), drop=FALSE],
    sample_metadata=group_meta,
    sample_time_map=sample_time_map,
    time_levels=time_levels
  )
}

score_to_long <- function(x, group_name, time_map=NULL,
                          id_col="ID", sample_parser=NULL,
                          sample_time_map=NULL, time_levels=NULL) {
  # Legacy time_map: named regular-expression patterns -> ordered time labels.
  # New code should pass sample_time_map built from metadata instead.
  z <- as.data.frame(x$score) |>
    tibble::rownames_to_column(id_col) |>
    tidyr::pivot_longer(-dplyr::all_of(id_col), names_to="Sample", values_to="Score") |>
    dplyr::left_join(x$annotation, by=id_col) |>
    dplyr::mutate(Group=group_name)
  if (!is.null(sample_time_map)) {
    if (is.null(names(sample_time_map)) || any(!nzchar(names(sample_time_map)))) {
      stop("sample_time_map 必须是以样本名为名称的命名向量。", call. = FALSE)
    }
    unassigned <- setdiff(unique(z$Sample), names(sample_time_map))
    if (length(unassigned)) {
      stop("这些 GSVA 样本无法由 metadata 匹配轮次：", paste(unassigned, collapse=", "), call. = FALSE)
    }
    z$Time <- unname(sample_time_map[z$Sample])
    time_levels <- time_levels %||% unique(unname(sample_time_map))
  } else if(!is.null(sample_parser)) {
    z$Time <- sample_parser(z$Sample)
  } else {
    if (is.null(time_map) || !length(time_map)) {
      stop("请提供 sample_metadata（推荐）或旧式 time_map。", call. = FALSE)
    }
    z$Time <- NA_character_
    for(p in names(time_map)) z$Time[stringr::str_detect(z$Sample, p)] <- unname(time_map[p])
    time_levels <- time_levels %||% unique(unname(time_map))
  }
  z$Time <- factor(z$Time, levels=time_levels, ordered=TRUE)
  if(any(is.na(z$Time))) warning(sum(is.na(z$Time)), " samples could not be assigned to a time point.")
  z
}

transform_scores <- function(df, zscore=TRUE, id_col="ID") {
  if(!zscore) return(df |> dplyr::mutate(Analysis_Score=Score))
  df |>
    dplyr::group_by(.data$Group, .data[[id_col]]) |>
    dplyr::mutate(Analysis_Score=ifelse(is.finite(stats::sd(.data$Score, na.rm=TRUE)) && stats::sd(.data$Score, na.rm=TRUE)>0,
                                        (.data$Score-mean(.data$Score,na.rm=TRUE))/stats::sd(.data$Score,na.rm=TRUE), 0)) |>
    dplyr::ungroup()
}

pathway_trajectory <- function(df, id_col="ID", description_col="Description", category_col="Category") {
  extras <- intersect(c("Category_CN","Category_Merged"), names(df))
  gr <- c("Group", category_col, extras, id_col, description_col, "Time")
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(gr))) |>
    dplyr::summarise(
      mean_raw_score=mean(.data$Score,na.rm=TRUE), mean_score=mean(.data$Analysis_Score,na.rm=TRUE),
      sd_score=stats::sd(.data$Analysis_Score,na.rm=TRUE), n_sample=dplyr::n(), .groups="drop"
    )
}

category_cas <- function(pathway_df, category_col="Category") {
  extras <- intersect(c("Category_CN","Category_Merged"), names(pathway_df))
  gr <- c("Group", category_col, extras, "Time")
  pathway_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(gr))) |>
    dplyr::summarise(
      CAS=median(.data$mean_score,na.rm=TRUE), CAS_mean=mean(.data$mean_score,na.rm=TRUE),
      n_pathway=dplyr::n_distinct(.data$ID), Q25=stats::quantile(.data$mean_score,.25,na.rm=TRUE),
      Q75=stats::quantile(.data$mean_score,.75,na.rm=TRUE), .groups="drop"
    )
}

add_delta <- function(df, value_col, by_cols, id_name="Interval") {
  df |>
    dplyr::arrange(dplyr::across(dplyr::all_of(by_cols)), .data$Time) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by_cols))) |>
    dplyr::mutate(
      Previous_Time=dplyr::lag(.data$Time), Delta=.data[[value_col]]-dplyr::lag(.data[[value_col]]),
      !!id_name := ifelse(!is.na(.data$Previous_Time), paste0(as.character(.data$Previous_Time),"_to_",as.character(.data$Time)), NA_character_)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$Delta))
}

suggest_thresholds <- function(pathway_delta, category_delta,
                               pathway_quantile=.25, category_quantile=.25,
                               amplitude_values=NULL, amplitude_quantile=.25,
                               dcs_default=.60) {
  # These are scale-adaptive operational suggestions, NOT significance cutoffs.
  ep <- unname(quantile(abs(pathway_delta$Delta), pathway_quantile, na.rm=TRUE))
  ec <- unname(quantile(abs(category_delta$Delta), category_quantile, na.rm=TRUE))
  ea <- if(!is.null(amplitude_values)) unname(quantile(amplitude_values, amplitude_quantile, na.rm=TRUE)) else NA_real_
  list(epsilon_pathway=ep, epsilon_category=ec, amplitude_cutoff=ea, dcs_cutoff=dcs_default,
       method=paste0("Q", pathway_quantile*100,"(|pathway delta|), Q",category_quantile*100,
                     "(|category delta|), Q",amplitude_quantile*100,"(amplitude); DCS empirical default"))
}

direction_metrics <- function(pathway_delta, category_delta, epsilon_pathway, epsilon_category, category_col="Category") {
  pd <- pathway_delta |>
    dplyr::rename(Delta_pathway=Delta) |>
    dplyr::mutate(Pathway_Direction=dplyr::case_when(
      .data$Delta_pathway > epsilon_pathway ~ "Up", .data$Delta_pathway < -epsilon_pathway ~ "Down", TRUE ~ "Stable"
    ))
  cd <- category_delta |>
    dplyr::rename(Delta_CAS=Delta) |>
    dplyr::select(Group, dplyr::all_of(category_col), Interval, Delta_CAS) |>
    dplyr::mutate(CAS_Direction=dplyr::case_when(
      .data$Delta_CAS > epsilon_category ~ "Up", .data$Delta_CAS < -epsilon_category ~ "Down", TRUE ~ "Stable"
    ))
  pd |>
    dplyr::left_join(cd, by=c("Group",category_col,"Interval")) |>
    dplyr::group_by(.data$Group, .data[[category_col]], .data$Interval, .data$Delta_CAS, .data$CAS_Direction) |>
    dplyr::summarise(
      n_total=dplyr::n(), n_up=sum(.data$Pathway_Direction=="Up"),
      n_down=sum(.data$Pathway_Direction=="Down"), n_stable=sum(.data$Pathway_Direction=="Stable"),
      P_up=.data$n_up/.data$n_total, P_down=.data$n_down/.data$n_total, P_stable=.data$n_stable/.data$n_total,
      NDI=(.data$n_up-.data$n_down)/.data$n_total,
      DCS_all=mean(.data$Pathway_Direction==.data$CAS_Direction),
      DCS_directional=dplyr::case_when(
        dplyr::first(.data$CAS_Direction)=="Up" ~ mean(.data$Pathway_Direction=="Up"),
        dplyr::first(.data$CAS_Direction)=="Down" ~ mean(.data$Pathway_Direction=="Down"),
        TRUE ~ mean(.data$Pathway_Direction=="Stable")
      ),
      .groups="drop"
    )
}

global_metrics <- function(cas, category_col="Category") {
  extras <- intersect(c("Category_CN","Category_Merged"), names(cas)); gr <- c("Group",category_col,extras)
  cas |>
    dplyr::arrange(dplyr::across(dplyr::all_of(gr)), .data$Time) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(gr))) |>
    dplyr::summarise(
      n_pathway=dplyr::first(.data$n_pathway),
      Trend_rho=if(dplyr::n_distinct(.data$Time)>=3) stats::cor(as.numeric(.data$Time),.data$CAS,method="spearman",use="complete.obs") else NA_real_,
      Amplitude=max(.data$CAS,na.rm=TRUE)-min(.data$CAS,na.rm=TRUE),
      Peak_Time=as.character(.data$Time[which.max(.data$CAS)]),
      Trough_Time=as.character(.data$Time[which.min(.data$CAS)]),
      Start_CAS=dplyr::first(.data$CAS), End_CAS=dplyr::last(.data$CAS),
      Net_Change=dplyr::last(.data$CAS)-dplyr::first(.data$CAS), .groups="drop"
    )
}

classify_patterns <- function(cas, direction, global, epsilon_category, dcs_cutoff=.60, amplitude_cutoff=0, category_col="Category") {
  # General classification uses sign sequence; named biological labels are most interpretable for 4 ordered time points.
  wide_delta <- add_delta(cas, "CAS", c("Group",category_col)) |>
    dplyr::select(Group, dplyr::all_of(category_col), Interval, Delta) |>
    dplyr::mutate(sign=dplyr::case_when(
      .data$Delta>epsilon_category~1L, .data$Delta< -epsilon_category~-1L, TRUE~0L
    )) |>
    dplyr::group_by(.data$Group,.data[[category_col]]) |>
    dplyr::summarise(
      Sign_Pattern=paste(.data$sign,collapse=","), n_up=sum(.data$sign==1),
      n_down=sum(.data$sign==-1), n_flat=sum(.data$sign==0), .groups="drop"
    )
  cons <- direction |>
    dplyr::group_by(.data$Group,.data[[category_col]]) |>
    dplyr::summarise(
      Mean_DCS=mean(.data$DCS_directional,na.rm=TRUE), Min_DCS=min(.data$DCS_directional,na.rm=TRUE),
      Mean_abs_NDI=mean(abs(.data$NDI),na.rm=TRUE), .groups="drop"
    )
  out <- global |>
    dplyr::left_join(wide_delta,by=c("Group",category_col)) |>
    dplyr::left_join(cons,by=c("Group",category_col)) |>
    dplyr::mutate(Pattern=dplyr::case_when(
      .data$Amplitude <= amplitude_cutoff ~ "Stable",
      .data$Mean_DCS < dcs_cutoff ~ "Heterogeneous",
      .data$n_down==0 & .data$n_up>=2 ~ "Sustained activation",
      .data$n_up==0 & .data$n_down>=2 ~ "Sustained suppression",
      .data$Sign_Pattern %in% c("1,-1,-1","1,-1,0") ~ "Early activation",
      .data$Sign_Pattern %in% c("0,1,-1","1,1,-1") ~ "Transient activation",
      .data$Sign_Pattern %in% c("-1,0,1","-1,1,1","0,0,1") ~ "Late activation",
      .data$Sign_Pattern %in% c("-1,1,1","-1,1,0") ~ "Early suppression",
      .data$Sign_Pattern %in% c("0,-1,1","-1,-1,1") ~ "Transient suppression",
      .data$Sign_Pattern %in% c("1,0,-1","0,0,-1") ~ "Late suppression",
      TRUE ~ "Fluctuating / biphasic"
    ))
  out
}

run_temporal_pathway_analysis <- function(score_matrix, annotation, group_name, time_map=NULL,
                                          sample_metadata=NULL,
                                          sample_column="Sample", group_column="Group", cycle_column="Cycle",
                                          cycle_order=NULL,
                                          zscore=TRUE,
                                          thresholds=list(mode="auto", pathway_quantile=.25, category_quantile=.25,
                                                          amplitude_quantile=.25, dcs_cutoff=.60,
                                                          epsilon_pathway=NULL, epsilon_category=NULL, amplitude_cutoff=NULL),
                                          id_col="ID", description_col="Description", category_col="Category") {
  metadata_time <- NULL
  if (!is.null(sample_metadata)) {
    metadata_time <- prepare_time_from_metadata(
      score_matrix=score_matrix, metadata=sample_metadata, group_name=group_name,
      sample_column=sample_column, group_column=group_column, cycle_column=cycle_column,
      cycle_order=cycle_order
    )
    score_matrix <- metadata_time$score_matrix
  }
  anno <- prepare_annotation(annotation, group_name, id_col=id_col,description_col=description_col,category_col=category_col)
  sub <- extract_score_subset(score_matrix, anno, id_col)
  long <- score_to_long(
    sub, group_name, time_map=time_map, id_col=id_col,
    sample_time_map=metadata_time$sample_time_map %||% NULL,
    time_levels=metadata_time$time_levels %||% NULL
  ) |> transform_scores(zscore,id_col)
  pathway <- pathway_trajectory(long,id_col,description_col,category_col)
  cas <- category_cas(pathway,category_col)
  pdelta <- add_delta(pathway,"mean_score",c("Group",category_col,id_col))
  cdelta <- add_delta(cas,"CAS",c("Group",category_col))
  global <- global_metrics(cas,category_col)

  mode <- thresholds$mode %||% "auto"
  auto <- suggest_thresholds(pdelta,cdelta,
    thresholds$pathway_quantile %||% .25, thresholds$category_quantile %||% .25,
    global$Amplitude, thresholds$amplitude_quantile %||% .25, thresholds$dcs_cutoff %||% .60)
  used <- auto
  if(mode=="manual") {
    used$epsilon_pathway <- thresholds$epsilon_pathway %||% stop("manual mode requires epsilon_pathway")
    used$epsilon_category <- thresholds$epsilon_category %||% stop("manual mode requires epsilon_category")
    used$amplitude_cutoff <- thresholds$amplitude_cutoff %||% stop("manual mode requires amplitude_cutoff")
    used$dcs_cutoff <- thresholds$dcs_cutoff %||% .60
    used$method <- "manual"
  }
  direction <- direction_metrics(pdelta,cdelta,used$epsilon_pathway,used$epsilon_category,category_col)
  patterns <- classify_patterns(cas,direction,global,used$epsilon_category,used$dcs_cutoff,used$amplitude_cutoff,category_col)
  summary <- patterns |> dplyr::arrange(.data$Pattern, dplyr::desc(.data$Amplitude))

  diagnostics <- list(
    abs_pathway_delta=quantile(abs(pdelta$Delta),c(.1,.25,.5,.75,.9),na.rm=TRUE),
    abs_category_delta=quantile(abs(cdelta$Delta),c(.1,.25,.5,.75,.9),na.rm=TRUE),
    amplitude=quantile(global$Amplitude,c(.1,.25,.5,.75,.9),na.rm=TRUE),
    dcs=quantile(direction$DCS_directional,c(.1,.25,.5,.75,.9),na.rm=TRUE),
    matching=sub$matching, missing_ids=sub$missing_ids)

  structure(list(settings=list(
                 group=group_name, zscore=zscore,
                 time_source=if (is.null(metadata_time)) "sample_name_pattern" else "metadata",
                 time_map=if (is.null(metadata_time)) time_map else metadata_time$sample_time_map,
                 time_levels=if (is.null(metadata_time)) unique(unname(time_map)) else metadata_time$time_levels,
                 metadata_columns=if (is.null(metadata_time)) NULL else list(
                   sample=sample_column, group=group_column, cycle=cycle_column
                 ),
                 thresholds=used),
                 annotation=sub$annotation,score_subset=sub$score,score_readable=make_readable_score(sub,description_col),
                 long=long,pathway_trajectory=pathway,CAS=cas,pathway_delta=pdelta,category_delta=cdelta,
                 direction_metrics=direction,global_metrics=global,summary=summary,diagnostics=diagnostics),
            class="temporal_pathway_analysis")
}

print.temporal_pathway_analysis <- function(x, ...) {
  cat("Temporal pathway analysis\n")
  cat("Group:",x$settings$group," | Z-score:",x$settings$zscore,"\n")
  cat("Pathways:",nrow(x$annotation)," | Categories:",dplyr::n_distinct(x$annotation$Category),"\n")
  cat("Thresholds used:\n"); print(x$settings$thresholds[c("epsilon_pathway","epsilon_category","amplitude_cutoff","dcs_cutoff","method")])
  cat("\nPattern counts:\n"); print(table(x$summary$Pattern))
  invisible(x)
}
