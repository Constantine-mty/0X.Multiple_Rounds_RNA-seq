# Temporal pathway/category trend analysis
# Dependencies: dplyr, tidyr, tibble, stringr, purrr

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr); library(purrr)
})

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
  if (!is.null(group)) x <- x %>% filter(.data[[group_col]] == group)
  keep <- unique(c(if(!is.null(group)) group_col, id_col, description_col, category_col, intersect(extra_cols, names(x))))
  x <- x %>% select(all_of(keep)) %>% distinct()
  conflicts <- x %>% distinct(across(all_of(c(id_col, category_col)))) %>% count(.data[[id_col]], name="n_category") %>% filter(n_category > 1)
  if(nrow(conflicts)) warning(nrow(conflicts), " pathway IDs map to >1 category. Resolve before interpretation.")
  x
}

extract_score_subset <- function(score_matrix, annotation, id_col="ID") {
  if(is.null(rownames(score_matrix))) stop("score_matrix must have pathway IDs as row names.")
  ids <- unique(annotation[[id_col]])
  matched <- intersect(ids, rownames(score_matrix)); missing <- setdiff(ids, rownames(score_matrix))
  anno <- annotation %>% filter(.data[[id_col]] %in% matched) %>% distinct(.data[[id_col]], .keep_all=TRUE) %>% arrange(match(.data[[id_col]], matched))
  mat <- score_matrix[anno[[id_col]], , drop=FALSE]
  list(score=mat, annotation=anno, missing_ids=missing,
       matching=data.frame(n_selected=length(ids), n_matched=length(matched), n_missing=length(missing)))
}

make_readable_score <- function(x, description_col="Description") {
  m <- x$score; rownames(m) <- make.unique(x$annotation[[description_col]]); m
}

score_to_long <- function(x, group_name, time_map,
                          id_col="ID", sample_parser=NULL) {
  # time_map: named character vector; names are regex patterns, values are ordered time labels.
  # Example c("_Non_"="R0", "_R1_"="R1", "_R2_"="R2", "_R3_"="R3")
  z <- as.data.frame(x$score) %>% rownames_to_column(id_col) %>%
    pivot_longer(-all_of(id_col), names_to="Sample", values_to="Score") %>%
    left_join(x$annotation, by=id_col) %>% mutate(Group=group_name)
  if(!is.null(sample_parser)) {
    z$Time <- sample_parser(z$Sample)
  } else {
    z$Time <- NA_character_
    for(p in names(time_map)) z$Time[str_detect(z$Sample, p)] <- unname(time_map[p])
  }
  time_levels <- unique(unname(time_map))
  z$Time <- factor(z$Time, levels=time_levels, ordered=TRUE)
  if(any(is.na(z$Time))) warning(sum(is.na(z$Time)), " samples could not be assigned to a time point.")
  z
}

transform_scores <- function(df, zscore=TRUE, id_col="ID") {
  if(!zscore) return(df %>% mutate(Analysis_Score=Score))
  df %>% group_by(Group, .data[[id_col]]) %>%
    mutate(Analysis_Score=ifelse(is.finite(sd(Score, na.rm=TRUE)) && sd(Score, na.rm=TRUE)>0,
                                 (Score-mean(Score,na.rm=TRUE))/sd(Score,na.rm=TRUE), 0)) %>% ungroup()
}

pathway_trajectory <- function(df, id_col="ID", description_col="Description", category_col="Category") {
  extras <- intersect(c("Category_CN","Category_Merged"), names(df))
  gr <- c("Group", category_col, extras, id_col, description_col, "Time")
  df %>% group_by(across(all_of(gr))) %>% summarise(
    mean_raw_score=mean(Score,na.rm=TRUE), mean_score=mean(Analysis_Score,na.rm=TRUE),
    sd_score=sd(Analysis_Score,na.rm=TRUE), n_sample=n(), .groups="drop")
}

category_cas <- function(pathway_df, category_col="Category") {
  extras <- intersect(c("Category_CN","Category_Merged"), names(pathway_df))
  gr <- c("Group", category_col, extras, "Time")
  pathway_df %>% group_by(across(all_of(gr))) %>% summarise(
    CAS=median(mean_score,na.rm=TRUE), CAS_mean=mean(mean_score,na.rm=TRUE),
    n_pathway=n_distinct(ID), Q25=quantile(mean_score,.25,na.rm=TRUE), Q75=quantile(mean_score,.75,na.rm=TRUE), .groups="drop")
}

add_delta <- function(df, value_col, by_cols, id_name="Interval") {
  df %>% arrange(across(all_of(by_cols)), Time) %>% group_by(across(all_of(by_cols))) %>%
    mutate(Previous_Time=lag(Time), Delta=.data[[value_col]]-lag(.data[[value_col]]),
           !!id_name := ifelse(!is.na(Previous_Time), paste0(as.character(Previous_Time),"_to_",as.character(Time)), NA_character_)) %>%
    ungroup() %>% filter(!is.na(Delta))
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
  pd <- pathway_delta %>% rename(Delta_pathway=Delta) %>% mutate(Pathway_Direction=case_when(
    Delta_pathway > epsilon_pathway ~ "Up", Delta_pathway < -epsilon_pathway ~ "Down", TRUE ~ "Stable"))
  cd <- category_delta %>% rename(Delta_CAS=Delta) %>% select(Group, all_of(category_col), Interval, Delta_CAS) %>%
    mutate(CAS_Direction=case_when(Delta_CAS > epsilon_category ~ "Up", Delta_CAS < -epsilon_category ~ "Down", TRUE ~ "Stable"))
  pd %>% left_join(cd, by=c("Group",category_col,"Interval")) %>%
    group_by(Group, .data[[category_col]], Interval, Delta_CAS, CAS_Direction) %>% summarise(
      n_total=n(), n_up=sum(Pathway_Direction=="Up"), n_down=sum(Pathway_Direction=="Down"), n_stable=sum(Pathway_Direction=="Stable"),
      P_up=n_up/n_total, P_down=n_down/n_total, P_stable=n_stable/n_total,
      NDI=(n_up-n_down)/n_total,
      DCS_all=mean(Pathway_Direction==CAS_Direction),
      DCS_directional=case_when(
        first(CAS_Direction)=="Up" ~ mean(Pathway_Direction=="Up"),
        first(CAS_Direction)=="Down" ~ mean(Pathway_Direction=="Down"),
        TRUE ~ mean(Pathway_Direction=="Stable")), .groups="drop")
}

global_metrics <- function(cas, category_col="Category") {
  extras <- intersect(c("Category_CN","Category_Merged"), names(cas)); gr <- c("Group",category_col,extras)
  cas %>% arrange(across(all_of(gr)), Time) %>% group_by(across(all_of(gr))) %>% summarise(
    n_pathway=first(n_pathway), Trend_rho=if(n_distinct(Time)>=3) cor(as.numeric(Time),CAS,method="spearman",use="complete.obs") else NA_real_,
    Amplitude=max(CAS,na.rm=TRUE)-min(CAS,na.rm=TRUE), Peak_Time=as.character(Time[which.max(CAS)]),
    Trough_Time=as.character(Time[which.min(CAS)]), Start_CAS=first(CAS), End_CAS=last(CAS), Net_Change=last(CAS)-first(CAS), .groups="drop")
}

classify_patterns <- function(cas, direction, global, epsilon_category, dcs_cutoff=.60, amplitude_cutoff=0, category_col="Category") {
  # General classification uses sign sequence; named biological labels are most interpretable for 4 ordered time points.
  wide_delta <- add_delta(cas, "CAS", c("Group",category_col)) %>% select(Group,all_of(category_col),Interval,Delta) %>%
    mutate(sign=case_when(Delta>epsilon_category~1L, Delta< -epsilon_category~-1L, TRUE~0L)) %>%
    group_by(Group,.data[[category_col]]) %>% summarise(Sign_Pattern=paste(sign,collapse=","), n_up=sum(sign==1),n_down=sum(sign==-1),n_flat=sum(sign==0),.groups="drop")
  cons <- direction %>% group_by(Group,.data[[category_col]]) %>% summarise(
    Mean_DCS=mean(DCS_directional,na.rm=TRUE), Min_DCS=min(DCS_directional,na.rm=TRUE), Mean_abs_NDI=mean(abs(NDI),na.rm=TRUE),.groups="drop")
  out <- global %>% left_join(wide_delta,by=c("Group",category_col)) %>% left_join(cons,by=c("Group",category_col)) %>%
    mutate(Pattern=case_when(
      Amplitude <= amplitude_cutoff ~ "Stable",
      Mean_DCS < dcs_cutoff ~ "Heterogeneous",
      n_down==0 & n_up>=2 ~ "Sustained activation",
      n_up==0 & n_down>=2 ~ "Sustained suppression",
      Sign_Pattern %in% c("1,-1,-1","1,-1,0") ~ "Early activation",
      Sign_Pattern %in% c("0,1,-1","1,1,-1") ~ "Transient activation",
      Sign_Pattern %in% c("-1,0,1","-1,1,1","0,0,1") ~ "Late activation",
      Sign_Pattern %in% c("-1,1,1","-1,1,0") ~ "Early suppression",
      Sign_Pattern %in% c("0,-1,1","-1,-1,1") ~ "Transient suppression",
      Sign_Pattern %in% c("1,0,-1","0,0,-1") ~ "Late suppression",
      TRUE ~ "Fluctuating / biphasic"))
  out
}

run_temporal_pathway_analysis <- function(score_matrix, annotation, group_name, time_map,
                                          zscore=TRUE,
                                          thresholds=list(mode="auto", pathway_quantile=.25, category_quantile=.25,
                                                          amplitude_quantile=.25, dcs_cutoff=.60,
                                                          epsilon_pathway=NULL, epsilon_category=NULL, amplitude_cutoff=NULL),
                                          id_col="ID", description_col="Description", category_col="Category") {
  anno <- prepare_annotation(annotation, group_name, id_col=id_col,description_col=description_col,category_col=category_col)
  sub <- extract_score_subset(score_matrix, anno, id_col)
  long <- score_to_long(sub, group_name, time_map, id_col=id_col) %>% transform_scores(zscore,id_col)
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
  summary <- patterns %>% arrange(Pattern,desc(Amplitude))

  diagnostics <- list(
    abs_pathway_delta=quantile(abs(pdelta$Delta),c(.1,.25,.5,.75,.9),na.rm=TRUE),
    abs_category_delta=quantile(abs(cdelta$Delta),c(.1,.25,.5,.75,.9),na.rm=TRUE),
    amplitude=quantile(global$Amplitude,c(.1,.25,.5,.75,.9),na.rm=TRUE),
    dcs=quantile(direction$DCS_directional,c(.1,.25,.5,.75,.9),na.rm=TRUE),
    matching=sub$matching, missing_ids=sub$missing_ids)

  structure(list(settings=list(group=group_name,zscore=zscore,time_map=time_map,thresholds=used),
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
