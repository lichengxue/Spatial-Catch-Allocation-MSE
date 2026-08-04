#### ============================================================
#### Catch allocation MSE figures: six main OM scenarios plus all OM scenarios
####
#### Main-text OM scenarios:
####   NBL, MBH, MUL, MUH, NBLD, NULD
####
#### Figure 3: Catch, relative difference from Equal allocation
#### Figure 4: SSB, relative difference from Equal allocation
####
#### Other outputs for six main OMs and all OMs:
####   AAV catch, AAV SSB, P(overfishing), P(overfished)
####   Summed scaled performance heatmaps for system, region 1, region 2
####
#### This version rebuilds performance2_*.RDS files once so that risk
#### probabilities, AAV catch, and AAV SSB are calculated across all
#### feedback years.
#### After those files are rebuilt, set RUN_PERFORMANCE_EXTRACTION <- FALSE
#### for faster replotting.
#### ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(scales)
  library(patchwork)
})

## ============================================================
## 0) User settings
## ============================================================

WORK_DIR <- "C:/Users/liche/Desktop/Catch_Paper"
setwd(WORK_DIR)

RUN_PERFORMANCE_EXTRACTION <- FALSE

MODEL_LEVELS <- c("Equal", "Gear-Equal", "Fleet-Catch", "Index-Equal", "Index-Gear")
BASELINE_MODEL <- "Equal"
REGION_LEVELS <- c("R1", "R2", "G")

## Raw file codes used in Mods_*.RDS and performance2_*.RDS.
## These are recoded below to the final scenario names used in plots.
SCENARIOS_MAIN_RAW <- c("NB", "MBH", "MU", "MUH", "NBC", "NUC")
SCENARIOS_MAIN <- c("NBL", "MBH", "MUL", "MUH", "NBLD", "NULD")

SCENARIOS_ALL_RAW <- c(
  "NB", "NBH", "NU", "NUH", "NBI", "NBO", "NBC", "NUC", "NBHS",
  "MB", "MBH", "MU", "MUH", "MBI", "MBO", "MBHS"
)
SCENARIOS_ALL <- c(
  "NBL", "NBH", "NUL", "NUH", "NBLI", "NBLO", "NBLD", "NULD", "NBH-BH",
  "MBL", "MBH", "MUL", "MUH", "MBLI", "MBLO", "MBH-BH"
)

## Factor levels for horizontal plots. ggplot draws the last y-level at the top,
## so reversing here makes the top-to-bottom order read as SCENARIOS_MAIN.
SCENARIOS_MAIN_Y <- rev(SCENARIOS_MAIN)
SCENARIOS_ALL_Y <- rev(SCENARIOS_ALL)

STRATEGY_COLORS <- c(
  "Equal"       = "#2874A6",
  "Gear-Equal"  = "#F39C12",
  "Fleet-Catch" = "#27AE60",
  "Index-Equal" = "#E74C3C",
  "Index-Gear"  = "#8E44AD"
)

REGION_LABELS <- c(
  "R1" = "Region 1",
  "R2" = "Region 2",
  "G"  = "System"
)

MODEL_LABELS_SHORT <- c(
  "Equal"       = "Equal",
  "Gear-Equal"  = "Gear-Equal",
  "Fleet-Catch" = "Fleet-Catch",
  "Index-Equal" = "Index-Equal",
  "Index-Gear"  = "Index-Gear"
)

## Relative catch and SSB are defined relative to Equal allocation in
## extract_performance_df_from_mods(). Equal is therefore exactly zero.
## Keeping Equal out of relative figures reduces clutter; the dashed zero line
## is the Equal reference. Change this to FALSE if you want Equal points shown.
DROP_EQUAL_IN_RELATIVE_FIGURES <- TRUE

## Optional fixed x-limits for Figure 3 and Figure 4.
## Leave NULL for automatic symmetric limits.
CATCH_RELATIVE_XLIM <- NULL
SSB_RELATIVE_XLIM <- NULL

SAVE_PDF <- TRUE

MAIN_FIG_WIDTH <- 14
MAIN_FIG_HEIGHT <- 18
ALL_OM_FIG_WIDTH <- 15
ALL_OM_FIG_HEIGHT <- 24

## Larger vertical separation among the five allocation strategies within
## each OM row. This is most helpful for the six-OM main-text figures.
STRATEGY_DODGE_WIDTH <- 0.82

## ============================================================
## 1) Scenario recoding
## ============================================================

recode_scenario <- function(x) {
  dplyr::recode(
    as.character(x),
    "NB"  = "NBL",
    "MB"  = "MBL",
    "NU"  = "NUL",
    "MU"  = "MUL",
    "NBI" = "NBLI",
    "MBI" = "MBLI",
    "NBO" = "NBLO",
    "MBO" = "MBLO",
    "NBC" = "NBLD",
    "NUC" = "NULD",
    "NBHS" = "NBH-BH",
    "MBHS" = "MBH-BH",
    .default = as.character(x)
  )
}

stopifnot(identical(recode_scenario(SCENARIOS_MAIN_RAW), SCENARIOS_MAIN))
stopifnot(identical(recode_scenario(SCENARIOS_ALL_RAW), SCENARIOS_ALL))

## ============================================================
## 2) Performance extractor
## ============================================================

extract_performance_df_from_mods <- function(mods,
                                             is.nsim,
                                             use.n.years.first = 10,
                                             use.n.years.last = 10,
                                             start.years = 36,
                                             new_model_names = NULL) {

  make_ref_mat <- function(ref, nyear, ncol_needed) {
    if (is.null(dim(ref))) {
      ref <- as.numeric(ref)
      if (length(ref) == 1) {
        return(matrix(ref, nrow = nyear, ncol = ncol_needed))
      }
      if (length(ref) == ncol_needed) {
        return(matrix(ref, nrow = nyear, ncol = ncol_needed, byrow = TRUE))
      }
      if (length(ref) == nyear) {
        return(matrix(ref, nrow = nyear, ncol = ncol_needed))
      }
      return(matrix(rep(ref, length.out = nyear * ncol_needed),
                    nrow = nyear, ncol = ncol_needed))
    }

    mat <- as.matrix(ref)
    if (nrow(mat) == nyear && ncol(mat) == ncol_needed) return(mat)
    if (nrow(mat) == 1 && ncol(mat) == ncol_needed) {
      return(matrix(as.numeric(mat[1, ]),
                    nrow = nyear, ncol = ncol_needed, byrow = TRUE))
    }
    if (ncol(mat) == 1 && nrow(mat) == ncol_needed) {
      return(matrix(as.numeric(mat[, 1]),
                    nrow = nyear, ncol = ncol_needed, byrow = TRUE))
    }
    matrix(rep(as.numeric(mat), length.out = nyear * ncol_needed),
           nrow = nyear, ncol = ncol_needed)
  }

  calc_aav_mat <- function(X_mat, idx) {
    idx <- sort(unique(idx))
    if (length(idx) < 2) return(rep(NA_real_, ncol(X_mat)))

    Xw <- X_mat[idx, , drop = FALSE]
    mu <- colMeans(Xw, na.rm = TRUE)
    dX <- abs(Xw[-1, , drop = FALSE] - Xw[-nrow(Xw), , drop = FALSE])
    out <- colMeans(dX, na.rm = TRUE) / mu
    out[!is.finite(out)] <- NA_real_
    out
  }

  if (!is.nsim) stop("Only is.nsim = TRUE is implemented in this script.")

  n_models <- length(mods[[1]])
  n_reps <- length(mods)
  results_list <- vector("list", 0)

  for (r in seq_len(n_reps)) {
    mod_list <- mods[[r]]
    model_labels <- paste0("Model", seq_len(n_models))
    ts_list <- vector("list", n_models)
    names(ts_list) <- model_labels

    for (m in seq_len(n_models)) {
      rep <- mod_list[[m]]$om$rep
      input <- mod_list[[m]]$om$input$data

      n_fleets <- input$n_fleets[1]
      n_regions <- input$n_regions[1]
      nyear <- nrow(rep$SSB)
      yrs <- seq_len(nyear)

      ssb_reg <- as.matrix(rep$SSB)
      ssb_system <- rowSums(ssb_reg)

      catch_fleet <- as.matrix(rep$pred_catch)
      if (!is.null(input$fleet_regions)) {
        fleet_regions <- input$fleet_regions
        catch_reg <- sapply(seq_len(n_regions), function(k) {
          idx <- which(fleet_regions == k)
          if (length(idx) == 0) return(rep(0, nrow(catch_fleet)))
          rowSums(catch_fleet[, idx, drop = FALSE])
        })
      } else {
        catch_reg <- catch_fleet
      }
      catch_system <- rowSums(catch_reg)

      Fbar_mat <- as.matrix(rep$Fbar)
      Fbar_reg <- Fbar_mat[, (n_fleets + 1):(n_fleets + n_regions), drop = FALSE]
      Fbar_system <- Fbar_mat[, n_fleets + n_regions + 1]

      SSB_mat <- cbind(ssb_reg, ssb_system)
      SSBref <- exp(rep$log_SSB_FXSPR)
      SSBref_mat <- make_ref_mat(SSBref, nyear, n_regions + 1)
      SSB_ratio <- SSB_mat / SSBref_mat

      F_mat <- cbind(Fbar_reg, Fbar_system)
      Fref <- exp(rep$log_Fbar_XSPR)
      Fref_mat <- make_ref_mat(Fref, nyear, n_regions + 1)
      F_ratio <- F_mat / Fref_mat

      ts_list[[m]] <- list(
        years = yrs,
        n_regions = n_regions,
        ssb_reg = ssb_reg,
        ssb_system = ssb_system,
        catch_reg = catch_reg,
        catch_system = catch_system,
        SSB_ratio = SSB_ratio,
        F_ratio = F_ratio
      )
    }

    yrs <- ts_list[[1]]$years
    nY <- length(yrs)

    idx_short <- start.years:(start.years + use.n.years.first - 1)
    idx_short <- idx_short[idx_short >= 1 & idx_short <= nY]

    idx_long <- (nY - use.n.years.last + 1):nY
    idx_long <- idx_long[idx_long >= 1 & idx_long <= nY]

    idx_feedback <- start.years:nY
    idx_feedback <- idx_feedback[idx_feedback >= 1 & idx_feedback <= nY]

    base_ts <- ts_list[[1]]

    base_short_reg_ssb <- colMeans(base_ts$ssb_reg[idx_short, , drop = FALSE], na.rm = TRUE)
    base_long_reg_ssb <- colMeans(base_ts$ssb_reg[idx_long, , drop = FALSE], na.rm = TRUE)
    base_short_system_ssb <- mean(base_ts$ssb_system[idx_short], na.rm = TRUE)
    base_long_system_ssb <- mean(base_ts$ssb_system[idx_long], na.rm = TRUE)

    base_short_reg_catch <- colMeans(base_ts$catch_reg[idx_short, , drop = FALSE], na.rm = TRUE)
    base_long_reg_catch <- colMeans(base_ts$catch_reg[idx_long, , drop = FALSE], na.rm = TRUE)
    base_short_system_catch <- mean(base_ts$catch_system[idx_short], na.rm = TRUE)
    base_long_system_catch <- mean(base_ts$catch_system[idx_long], na.rm = TRUE)

    for (m in seq_len(n_models)) {
      lab <- model_labels[m]
      ts <- ts_list[[m]]
      n_regions <- ts$n_regions
      region_names <- paste0("R", seq_len(n_regions))

      short_reg_ssb <- colMeans(ts$ssb_reg[idx_short, , drop = FALSE], na.rm = TRUE)
      short_system_ssb <- mean(ts$ssb_system[idx_short], na.rm = TRUE)
      short_reg_catch <- colMeans(ts$catch_reg[idx_short, , drop = FALSE], na.rm = TRUE)
      short_system_catch <- mean(ts$catch_system[idx_short], na.rm = TRUE)

      long_reg_ssb <- colMeans(ts$ssb_reg[idx_long, , drop = FALSE], na.rm = TRUE)
      long_system_ssb <- mean(ts$ssb_system[idx_long], na.rm = TRUE)
      long_reg_catch <- colMeans(ts$catch_reg[idx_long, , drop = FALSE], na.rm = TRUE)
      long_system_catch <- mean(ts$catch_system[idx_long], na.rm = TRUE)

      rel_ssb_short_reg <- (short_reg_ssb - base_short_reg_ssb) / base_short_reg_ssb
      rel_ssb_long_reg <- (long_reg_ssb - base_long_reg_ssb) / base_long_reg_ssb
      rel_ssb_short_system <- (short_system_ssb - base_short_system_ssb) / base_short_system_ssb
      rel_ssb_long_system <- (long_system_ssb - base_long_system_ssb) / base_long_system_ssb

      rel_catch_short_reg <- (short_reg_catch - base_short_reg_catch) / base_short_reg_catch
      rel_catch_long_reg <- (long_reg_catch - base_long_reg_catch) / base_long_reg_catch
      rel_catch_short_system <- (short_system_catch - base_short_system_catch) / base_short_system_catch
      rel_catch_long_system <- (long_system_catch - base_long_system_catch) / base_long_system_catch

      p_overfish_all <- colMeans(ts$F_ratio[idx_feedback, , drop = FALSE] > 1, na.rm = TRUE)
      p_overfished_all <- colMeans(ts$SSB_ratio[idx_feedback, , drop = FALSE] < 0.5, na.rm = TRUE)

      ssb_all_mat <- cbind(ts$ssb_reg, ts$ssb_system)
      catch_all_mat <- cbind(ts$catch_reg, ts$catch_system)

      aav_ssb_all <- calc_aav_mat(ssb_all_mat, idx_feedback)
      aav_catch_all <- calc_aav_mat(catch_all_mat, idx_feedback)

      results_list[[length(results_list) + 1L]] <- data.frame(
        Realization = r,
        Model = lab,
        Period = "short",
        Scale = "regional",
        Region = region_names,
        SSB = as.numeric(short_reg_ssb),
        Catch = as.numeric(short_reg_catch),
        rel_SSB = as.numeric(rel_ssb_short_reg),
        rel_Catch = as.numeric(rel_catch_short_reg),
        AAV_SSB = NA_real_,
        AAV_Catch = NA_real_,
        p_overfishing = NA_real_,
        p_overfished = NA_real_,
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1L]] <- data.frame(
        Realization = r,
        Model = lab,
        Period = "short",
        Scale = "system",
        Region = "G",
        SSB = as.numeric(short_system_ssb),
        Catch = as.numeric(short_system_catch),
        rel_SSB = as.numeric(rel_ssb_short_system),
        rel_Catch = as.numeric(rel_catch_short_system),
        AAV_SSB = NA_real_,
        AAV_Catch = NA_real_,
        p_overfishing = NA_real_,
        p_overfished = NA_real_,
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1L]] <- data.frame(
        Realization = r,
        Model = lab,
        Period = "long",
        Scale = "regional",
        Region = region_names,
        SSB = as.numeric(long_reg_ssb),
        Catch = as.numeric(long_reg_catch),
        rel_SSB = as.numeric(rel_ssb_long_reg),
        rel_Catch = as.numeric(rel_catch_long_reg),
        AAV_SSB = NA_real_,
        AAV_Catch = NA_real_,
        p_overfishing = NA_real_,
        p_overfished = NA_real_,
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1L]] <- data.frame(
        Realization = r,
        Model = lab,
        Period = "long",
        Scale = "system",
        Region = "G",
        SSB = as.numeric(long_system_ssb),
        Catch = as.numeric(long_system_catch),
        rel_SSB = as.numeric(rel_ssb_long_system),
        rel_Catch = as.numeric(rel_catch_long_system),
        AAV_SSB = NA_real_,
        AAV_Catch = NA_real_,
        p_overfishing = NA_real_,
        p_overfished = NA_real_,
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1L]] <- data.frame(
        Realization = r,
        Model = lab,
        Period = "all",
        Scale = "regional",
        Region = region_names,
        SSB = NA_real_,
        Catch = NA_real_,
        rel_SSB = NA_real_,
        rel_Catch = NA_real_,
        AAV_SSB = as.numeric(aav_ssb_all[1:n_regions]),
        AAV_Catch = as.numeric(aav_catch_all[1:n_regions]),
        p_overfishing = as.numeric(p_overfish_all[1:n_regions]),
        p_overfished = as.numeric(p_overfished_all[1:n_regions]),
        stringsAsFactors = FALSE
      )

      results_list[[length(results_list) + 1L]] <- data.frame(
        Realization = r,
        Model = lab,
        Period = "all",
        Scale = "system",
        Region = "G",
        SSB = NA_real_,
        Catch = NA_real_,
        rel_SSB = NA_real_,
        rel_Catch = NA_real_,
        AAV_SSB = as.numeric(aav_ssb_all[n_regions + 1]),
        AAV_Catch = as.numeric(aav_catch_all[n_regions + 1]),
        p_overfishing = as.numeric(p_overfish_all[n_regions + 1]),
        p_overfished = as.numeric(p_overfished_all[n_regions + 1]),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- bind_rows(results_list)
  rownames(out) <- NULL

  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(out$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    out$Model <- factor(
      out$Model,
      levels = paste0("Model", seq_along(new_model_names)),
      labels = new_model_names
    )
  }

  out
}

## ============================================================
## 3) Optional rebuild of performance2_*.RDS for all scenarios
## ============================================================

if (RUN_PERFORMANCE_EXTRACTION) {
  mods_files <- paste0("Mods_", SCENARIOS_ALL_RAW, ".RDS")
  missing_mods <- mods_files[!file.exists(mods_files)]
  if (length(missing_mods) > 0) {
    stop("Missing Mods files:\n", paste(missing_mods, collapse = "\n"))
  }

  for (i in seq_along(mods_files)) {
    cat("Processing", mods_files[i], "\n")
    mods_i <- readRDS(mods_files[i])

    ## Keep the five catch allocation strategies used in the manuscript.
    mods_i <- lapply(mods_i, function(x) {
      if (length(x) >= length(MODEL_LEVELS)) x[seq_along(MODEL_LEVELS)] else x
    })

    perf_i <- extract_performance_df_from_mods(
      mods = mods_i,
      is.nsim = TRUE,
      start.years = 36,
      use.n.years.first = 10,
      use.n.years.last = 10,
      new_model_names = MODEL_LEVELS
    )

    perf_i$Scenario <- SCENARIOS_ALL_RAW[i]
    saveRDS(perf_i, file = paste0("performance2_", SCENARIOS_ALL_RAW[i], ".RDS"))
    rm(mods_i, perf_i)
    gc()
  }
}

## ============================================================
## 4) Load performance data
## ============================================================

load_performance <- function(scenario_raw, scenario_levels) {
  perf_files <- paste0("performance2_", scenario_raw, ".RDS")
  missing_perf <- perf_files[!file.exists(perf_files)]
  if (length(missing_perf) > 0) {
    stop("Missing performance files:\n", paste(missing_perf, collapse = "\n"))
  }

  bind_rows(lapply(seq_along(perf_files), function(i) {
    df <- readRDS(perf_files[i])
    if (!"Scenario" %in% names(df)) df$Scenario <- scenario_raw[i]
    df$Scenario_file <- scenario_raw[i]
    df
  })) %>%
    mutate(
      Scenario_raw = as.character(Scenario_file),
      Scenario = recode_scenario(Scenario_raw),
      Scenario = factor(Scenario, levels = scenario_levels),
      Model = factor(as.character(Model), levels = MODEL_LEVELS),
      Region = factor(as.character(Region), levels = REGION_LEVELS),
      Period = factor(as.character(Period), levels = c("short", "long", "all"))
    ) %>%
    filter(
      Scenario %in% scenario_levels,
      Model %in% MODEL_LEVELS,
      Region %in% REGION_LEVELS,
      Period %in% c("short", "long", "all")
    )
}

add_relative_ratio_columns <- function(perf) {
  perf %>%
  mutate(
    Catch_ratio_to_equal = 1 + rel_Catch,
    SSB_ratio_to_equal = 1 + rel_SSB,
    log_Catch_ratio_to_equal = if_else(Catch_ratio_to_equal > 0,
                                       log(Catch_ratio_to_equal), NA_real_),
    log_SSB_ratio_to_equal = if_else(SSB_ratio_to_equal > 0,
                                     log(SSB_ratio_to_equal), NA_real_)
  )
}

perf_main <- load_performance(SCENARIOS_MAIN_RAW, SCENARIOS_MAIN) %>%
  add_relative_ratio_columns()

perf_all_om <- load_performance(SCENARIOS_ALL_RAW, SCENARIOS_ALL) %>%
  add_relative_ratio_columns()

perf_all <- perf_main

cat("Loaded", nrow(perf_main), "rows for main scenarios:",
    paste(SCENARIOS_MAIN, collapse = ", "), "\n")
cat("Loaded", nrow(perf_all_om), "rows for all OM scenarios:",
    paste(SCENARIOS_ALL, collapse = ", "), "\n")

## ============================================================
## 5) Shared plot helpers
## ============================================================

theme_manuscript <- function() {
  theme_bw() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 13, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
      plot.title = element_text(size = 17, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

make_scenario_bands <- function(scenario_levels, fill = "grey95") {
  idx <- seq_along(scenario_levels)
  band_levels <- scenario_levels[idx %% 2 == 1]
  data.frame(
    Scenario = factor(band_levels, levels = scenario_levels),
    xmin = -Inf,
    xmax = Inf,
    ymin = as.numeric(factor(band_levels, levels = scenario_levels)) - 0.5,
    ymax = as.numeric(factor(band_levels, levels = scenario_levels)) + 0.5,
    fill = fill
  )
}

save_figure <- function(plot, filename, width = 14, height = 9, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
  if (SAVE_PDF) {
    pdf_name <- sub("\\.png$", ".pdf", filename)
    ggsave(pdf_name, plot, width = width, height = height)
  }
  invisible(plot)
}

auto_symmetric_xlim <- function(sum_dat, pad = 1.08, minimum = 0.05, include_iqr = TRUE) {
  vals <- if (include_iqr) c(sum_dat$Q1, sum_dat$Q3, sum_dat$Median) else sum_dat$Median
  max_abs <- suppressWarnings(max(abs(vals), na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- minimum
  c(-max_abs * pad, max_abs * pad)
}

## ============================================================
## 6) Median and IQR plots
## ============================================================

make_metric_summary <- function(perf,
                                value_col,
                                metric_prefix,
                                scenario_levels_y = SCENARIOS_MAIN_Y,
                                period_keep = c("short", "long"),
                                period_suffix = c(
                                  short = "ST",
                                  long = "LT",
                                  all = "All"
                                ),
                                drop_equal = FALSE) {
  model_keep <- if (drop_equal) setdiff(MODEL_LEVELS, BASELINE_MODEL) else MODEL_LEVELS
  metric_levels <- paste0(metric_prefix, "_", unname(period_suffix[period_keep]))

  perf %>%
    filter(
      Model %in% model_keep,
      as.character(Period) %in% period_keep
    ) %>%
    mutate(
      Scenario = factor(as.character(Scenario), levels = scenario_levels_y),
      Model = factor(as.character(Model), levels = model_keep),
      Region = factor(as.character(Region), levels = REGION_LEVELS),
      Metric = paste0(metric_prefix, "_", period_suffix[as.character(Period)]),
      Value = .data[[value_col]]
    ) %>%
    filter(!is.na(Value), !is.na(Metric)) %>%
    mutate(Metric = factor(Metric, levels = metric_levels)) %>%
    group_by(Scenario, Region, Model, Metric) %>%
    summarise(
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_metric_iqr <- function(sum_dat,
                            filename,
                            xlab,
                            metric_labels,
                            xlim = NULL,
                            auto_xlim = TRUE,
                            facet_scales = "fixed",
                            percent_axis = FALSE,
                            ratio_axis = FALSE,
                            add_zero_line = FALSE,
                            show_iqr = TRUE,
                            clip_to_xlim = FALSE,
                            median_segment_from = NULL,
                            median_segment_linewidth = 1.15,
                            median_segment_alpha = 0.28,
                            width = 14,
                            height = 9,
                            title = NULL) {
  scen_levels <- levels(sum_dat$Scenario)
  model_levels <- levels(droplevels(sum_dat$Model))
  sum_dat <- sum_dat %>%
    mutate(
      ## ggplot dodges horizontal y-axis groups from bottom to top. Reversing
      ## only the dodge factor makes the visible top-to-bottom strategy order
      ## match the legend order.
      Model_dodge = factor(as.character(Model), levels = rev(model_levels))
    )
  band_df <- make_scenario_bands(scen_levels)
  dodge <- position_dodge(width = STRATEGY_DODGE_WIDTH)

  if (is.null(xlim) && auto_xlim && facet_scales == "fixed") {
    xlim <- auto_symmetric_xlim(sum_dat, include_iqr = show_iqr)
  }
  x_range <- if (!is.null(xlim)) diff(xlim) else NA_real_

  if (clip_to_xlim && !is.null(xlim)) {
    sum_dat <- sum_dat %>%
      mutate(
        Median_raw = Median,
        Q1_raw = Q1,
        Q3_raw = Q3,
        Median = pmin(pmax(Median, xlim[1]), xlim[2]),
        Q1 = pmin(pmax(Q1, xlim[1]), xlim[2]),
        Q3 = pmin(pmax(Q3, xlim[1]), xlim[2]),
        clipped_low = Median_raw < xlim[1] | Q1_raw < xlim[1] | Q3_raw < xlim[1],
        clipped_high = Median_raw > xlim[2] | Q1_raw > xlim[2] | Q3_raw > xlim[2]
      )
  } else {
    sum_dat <- sum_dat %>%
      mutate(
        clipped_low = FALSE,
        clipped_high = FALSE
      )
  }

  p <- ggplot(sum_dat, aes(y = Scenario, color = Model, group = Model_dodge)) +
    geom_rect(
      data = band_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey95",
      color = NA
    )

  if (show_iqr) {
    p <- p +
      geom_errorbarh(
        aes(xmin = Q1, xmax = Q3),
        position = dodge,
        height = 0.22,
        linewidth = 0.35,
        alpha = 0.35
      )
  }

  if (!is.null(median_segment_from)) {
    p <- p +
      geom_errorbarh(
        aes(
          xmin = median_segment_from,
          xmax = Median
        ),
        position = dodge,
        height = 0,
        linewidth = median_segment_linewidth,
        alpha = median_segment_alpha
      )
  }

  p <- p +
    geom_point(
      aes(x = Median),
      position = dodge,
      size = 2.7
    )

  if (clip_to_xlim && !is.null(xlim)) {
    p <- p +
      geom_text(
        data = sum_dat %>% filter(clipped_low),
        aes(
          x = xlim[1] + 0.015 * x_range,
          label = "<"
        ),
        position = dodge,
        fontface = "bold",
        size = 3.5,
        show.legend = FALSE
      ) +
      geom_text(
        data = sum_dat %>% filter(clipped_high),
        aes(
          x = xlim[2] - 0.015 * x_range,
          label = ">"
        ),
        position = dodge,
        fontface = "bold",
        size = 3.5,
        show.legend = FALSE
      )
  }

  p <- p +
    facet_grid(
      Metric ~ Region,
      scales = facet_scales,
      labeller = labeller(
        Metric = metric_labels,
        Region = REGION_LABELS
      )
    ) +
    scale_color_manual(
      values = STRATEGY_COLORS[model_levels],
      breaks = model_levels,
      labels = MODEL_LABELS_SHORT[model_levels],
      name = "Strategy"
    ) +
    labs(title = title, x = xlab, y = "Operating model scenario") +
    theme_manuscript()

  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }

  if (add_zero_line) {
    p <- p +
      geom_vline(xintercept = 0, linetype = "dashed",
                 linewidth = 0.45, color = "grey35")
  }

  if (ratio_axis) {
    ratio_breaks <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10, 25, 50, 100, 250, 500)
    p <- p +
      scale_x_continuous(
        breaks = log(ratio_breaks),
        labels = paste0(ratio_breaks, "x")
      )
  } else if (percent_axis) {
    p <- p + scale_x_continuous(labels = percent_format(accuracy = 1))
  }

  save_figure(p, filename, width = width, height = height)
}

## ============================================================
## 7) Main text figures: relative catch and relative SSB
## ============================================================

catch_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "rel_Catch",
  metric_prefix = "Catch",
  drop_equal = DROP_EQUAL_IN_RELATIVE_FIGURES
)

plot_metric_iqr(
  sum_dat = catch_sum,
  filename = "Figure_3_Catch_Relative_6OM.png",
  xlab = "Relative difference from Equal allocation",
  metric_labels = c("Catch_ST" = "Catch ST", "Catch_LT" = "Catch LT"),
  xlim = CATCH_RELATIVE_XLIM,
  percent_axis = TRUE,
  add_zero_line = TRUE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

ssb_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "rel_SSB",
  metric_prefix = "SSB",
  drop_equal = DROP_EQUAL_IN_RELATIVE_FIGURES
)

plot_metric_iqr(
  sum_dat = ssb_sum,
  filename = "Figure_4_SSB_Relative_6OM.png",
  xlab = "Relative difference from Equal allocation",
  metric_labels = c("SSB_ST" = "SSB ST", "SSB_LT" = "SSB LT"),
  xlim = SSB_RELATIVE_XLIM,
  auto_xlim = FALSE,
  facet_scales = "free_x",
  percent_axis = TRUE,
  add_zero_line = TRUE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

## Absolute-scale versions for cases where relative SSB becomes unstable
## because the Equal baseline is very small.
catch_abs_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "Catch",
  metric_prefix = "Catch",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = catch_abs_sum,
  filename = "Figure4_Catch_Absolute_6OM.png",
  xlab = "Catch",
  metric_labels = c("Catch_ST" = "Catch ST", "Catch_LT" = "Catch LT"),
  auto_xlim = FALSE,
  facet_scales = "free_x",
  percent_axis = FALSE,
  add_zero_line = FALSE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

ssb_abs_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "SSB",
  metric_prefix = "SSB",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = ssb_abs_sum,
  filename = "Figure5_SSB_Absolute_6OM.png",
  xlab = "SSB",
  metric_labels = c("SSB_ST" = "SSB ST", "SSB_LT" = "SSB LT"),
  auto_xlim = FALSE,
  facet_scales = "free_x",
  percent_axis = FALSE,
  add_zero_line = FALSE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

## ============================================================
## 8) Supporting metric figures
## ============================================================

aav_catch_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "AAV_Catch",
  metric_prefix = "AAV_Catch",
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = aav_catch_sum,
  filename = "Figure_S2_AAV_Catch_6OM.png",
  xlab = "Average annual variation in catch",
  metric_labels = c("AAV_Catch_All" = "AAV catch"),
  xlim = c(0, max(aav_catch_sum$Median, na.rm = TRUE) * 1.08),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

aav_ssb_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "AAV_SSB",
  metric_prefix = "AAV_SSB",
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = aav_ssb_sum,
  filename = "Figure_S3_AAV_SSB_6OM.png",
  xlab = "Average annual variation in SSB",
  metric_labels = c("AAV_SSB_All" = "AAV SSB"),
  xlim = c(0, max(aav_ssb_sum$Median, na.rm = TRUE) * 1.08),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

pof_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "p_overfishing",
  metric_prefix = "POF",
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = pof_sum,
  filename = "Figure_S4_P_overfishing_6OM.png",
  xlab = "Probability of overfishing",
  metric_labels = c("POF_All" = "P(overfishing)"),
  xlim = c(0, 1),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

pofd_sum <- make_metric_summary(
  perf = perf_all,
  value_col = "p_overfished",
  metric_prefix = "POFD",
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = pofd_sum,
  filename = "Figure_S5_P_overfished_6OM.png",
  xlab = "Probability of being overfished",
  metric_labels = c("POFD_All" = "P(overfished)"),
  xlim = c(0, 1),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  show_iqr = FALSE,
  median_segment_from = 0,
  width = MAIN_FIG_WIDTH,
  height = MAIN_FIG_HEIGHT
)

## ============================================================
## 8b) All-OM supporting figures
## ============================================================

catch_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "rel_Catch",
  metric_prefix = "Catch",
  scenario_levels_y = SCENARIOS_ALL_Y,
  drop_equal = DROP_EQUAL_IN_RELATIVE_FIGURES
)

plot_metric_iqr(
  sum_dat = catch_sum_all_om,
  filename = "AllOM_Catch_Relative.png",
  xlab = "Relative difference from Equal allocation",
  metric_labels = c("Catch_ST" = "Catch ST", "Catch_LT" = "Catch LT"),
  xlim = CATCH_RELATIVE_XLIM,
  percent_axis = TRUE,
  add_zero_line = TRUE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

ssb_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "rel_SSB",
  metric_prefix = "SSB",
  scenario_levels_y = SCENARIOS_ALL_Y,
  drop_equal = DROP_EQUAL_IN_RELATIVE_FIGURES
)

plot_metric_iqr(
  sum_dat = ssb_sum_all_om,
  filename = "AllOM_SSB_Relative.png",
  xlab = "Relative difference from Equal allocation",
  metric_labels = c("SSB_ST" = "SSB ST", "SSB_LT" = "SSB LT"),
  xlim = SSB_RELATIVE_XLIM,
  auto_xlim = FALSE,
  facet_scales = "free_x",
  percent_axis = TRUE,
  add_zero_line = TRUE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

catch_abs_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "Catch",
  metric_prefix = "Catch",
  scenario_levels_y = SCENARIOS_ALL_Y,
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = catch_abs_sum_all_om,
  filename = "Figure_S6_Catch_AllOM.png",
  xlab = "Catch",
  metric_labels = c("Catch_ST" = "Catch ST", "Catch_LT" = "Catch LT"),
  auto_xlim = FALSE,
  facet_scales = "free_x",
  percent_axis = FALSE,
  add_zero_line = FALSE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

ssb_abs_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "SSB",
  metric_prefix = "SSB",
  scenario_levels_y = SCENARIOS_ALL_Y,
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = ssb_abs_sum_all_om,
  filename = "Figure_S7_SSB_AllOM.png",
  xlab = "SSB",
  metric_labels = c("SSB_ST" = "SSB ST", "SSB_LT" = "SSB LT"),
  auto_xlim = FALSE,
  facet_scales = "free_x",
  percent_axis = FALSE,
  add_zero_line = FALSE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

aav_catch_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "AAV_Catch",
  metric_prefix = "AAV_Catch",
  scenario_levels_y = SCENARIOS_ALL_Y,
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = aav_catch_sum_all_om,
  filename = "Figure_S8_AAV_Catch_AllOM.png",
  xlab = "Average annual variation in catch",
  metric_labels = c("AAV_Catch_All" = "AAV catch"),
  xlim = c(0, max(aav_catch_sum_all_om$Q3, na.rm = TRUE) * 1.08),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

aav_ssb_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "AAV_SSB",
  metric_prefix = "AAV_SSB",
  scenario_levels_y = SCENARIOS_ALL_Y,
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = aav_ssb_sum_all_om,
  filename = "Figure_S9_AAV_SSB_AllOM.png",
  xlab = "Average annual variation in SSB",
  metric_labels = c("AAV_SSB_All" = "AAV SSB"),
  xlim = c(0, max(aav_ssb_sum_all_om$Q3, na.rm = TRUE) * 1.08),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

pof_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "p_overfishing",
  metric_prefix = "POF",
  scenario_levels_y = SCENARIOS_ALL_Y,
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = pof_sum_all_om,
  filename = "Figure_S10_P_overfishing_AllOM.png",
  xlab = "Probability of overfishing",
  metric_labels = c("POF_All" = "P(overfishing)"),
  xlim = c(0, 1),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

pofd_sum_all_om <- make_metric_summary(
  perf = perf_all_om,
  value_col = "p_overfished",
  metric_prefix = "POFD",
  scenario_levels_y = SCENARIOS_ALL_Y,
  period_keep = "all",
  drop_equal = FALSE
)

plot_metric_iqr(
  sum_dat = pofd_sum_all_om,
  filename = "Figure_S11_P_overfished_AllOM.png",
  xlab = "Probability of being overfished",
  metric_labels = c("POFD_All" = "P(overfished)"),
  xlim = c(0, 1),
  percent_axis = FALSE,
  add_zero_line = FALSE,
  width = ALL_OM_FIG_WIDTH,
  height = ALL_OM_FIG_HEIGHT
)

## ============================================================
## 9) Summed scaled performance scores
## ============================================================

METRIC_LEVELS <- c(
  "Catch_ST",
  "SSB_ST",
  "Catch_LT",
  "SSB_LT",
  "AAV_Catch_All",
  "AAV_SSB_All",
  "p_overfishing_All",
  "p_overfished_All"
)

METRIC_LABELS <- c(
  "Catch_ST" = "Catch ST",
  "SSB_ST" = "SSB ST",
  "Catch_LT" = "Catch LT",
  "SSB_LT" = "SSB LT",
  "AAV_Catch_All" = "AAV catch",
  "AAV_SSB_All" = "AAV SSB",
  "p_overfishing_All" = "P(overfishing)",
  "p_overfished_All" = "P(overfished)"
)

GOOD_METRICS <- c("Catch_ST", "SSB_ST", "Catch_LT", "SSB_LT")
SUMMED_SCORE_MAX <- length(METRIC_LEVELS)

score_metric_0_1 <- function(x, metric_name) {
  x <- as.numeric(x)
  metric_name <- as.character(metric_name)[1]

  if (all(is.na(x))) return(rep(NA_real_, length(x)))

  rng <- suppressWarnings(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) {
    if (metric_name %in% GOOD_METRICS) {
      return(rep(1, length(x)))
    }
    if (grepl("^p_", metric_name)) {
      p <- x[which(!is.na(x))[1]]
      p <- max(0, min(1, p))
      return(rep(1 - p, length(x)))
    }
    return(rep(1, length(x)))
  }

  if (metric_name %in% GOOD_METRICS) {
    return((x - min(x, na.rm = TRUE)) / rng)
  }

  (max(x, na.rm = TRUE) - x) / rng
}

make_metric_long <- function(perf, scenario_levels) {
  perf %>%
    mutate(
      Period_chr = as.character(Period)
    ) %>%
    transmute(
      Scenario,
      Realization,
      Region,
      Model,
      Catch_ST = if_else(Period_chr == "short", rel_Catch, NA_real_),
      Catch_LT = if_else(Period_chr == "long", rel_Catch, NA_real_),
      SSB_ST = if_else(Period_chr == "short", rel_SSB, NA_real_),
      SSB_LT = if_else(Period_chr == "long", rel_SSB, NA_real_),
      AAV_Catch_All = if_else(Period_chr == "all", AAV_Catch, NA_real_),
      AAV_SSB_All = if_else(Period_chr == "all", AAV_SSB, NA_real_),
      p_overfishing_All = if_else(Period_chr == "all", p_overfishing, NA_real_),
      p_overfished_All = if_else(Period_chr == "all", p_overfished, NA_real_)
    ) %>%
    pivot_longer(
      cols = all_of(METRIC_LEVELS),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value)) %>%
    mutate(
      Scenario = factor(as.character(Scenario), levels = scenario_levels),
      Region = factor(as.character(Region), levels = REGION_LEVELS),
      Model = factor(as.character(Model), levels = MODEL_LEVELS),
      Metric = factor(Metric, levels = METRIC_LEVELS)
    )
}

make_scores <- function(perf, scenario_levels) {
  scores_long <- make_metric_long(perf, scenario_levels) %>%
    group_by(Scenario, Realization, Region, Model, Metric) %>%
    summarise(Value_real = median(Value, na.rm = TRUE), .groups = "drop") %>%
    group_by(Scenario, Region, Model, Metric) %>%
    summarise(Value = mean(Value_real, na.rm = TRUE), .groups = "drop") %>%
    group_by(Scenario, Region, Metric) %>%
    mutate(Score = score_metric_0_1(Value, Metric)) %>%
    ungroup()

  score_totals <- scores_long %>%
    group_by(Scenario, Region, Model) %>%
    summarise(Summed_score = sum(Score, na.rm = TRUE), .groups = "drop")

  list(scores_long = scores_long, score_totals = score_totals)
}

plot_score_heatmap <- function(region_keep,
                               scores_long,
                               score_totals,
                               scenario_levels,
                               filename,
                               title,
                               ncol = 3,
                               width = 12.5,
                               height = 10.5,
                               summed_score_ylim = c(0, SUMMED_SCORE_MAX)) {
  dat <- scores_long %>%
    filter(Region == region_keep) %>%
    mutate(
      Scenario = factor(as.character(Scenario), levels = scenario_levels),
      Metric = factor(as.character(Metric), levels = rev(METRIC_LEVELS)),
      Model = factor(as.character(Model), levels = MODEL_LEVELS)
    )

  totals <- score_totals %>%
    filter(Region == region_keep) %>%
    mutate(
      Scenario = factor(as.character(Scenario), levels = scenario_levels),
      Model = factor(as.character(Model), levels = MODEL_LEVELS)
    )

  scen_levels <- scenario_levels

  one_panel <- function(sc) {
    d_sc <- dat %>% filter(Scenario == sc)
    t_sc <- totals %>%
      filter(Scenario == sc) %>%
      mutate(Winner = Summed_score == max(Summed_score, na.rm = TRUE))

    p_bar <- ggplot(t_sc, aes(x = Model, y = Summed_score, fill = Winner)) +
      geom_col(color = "black", linewidth = 0.2) +
      scale_fill_manual(
        values = c("FALSE" = "grey82", "TRUE" = "#2A9D8F"),
        guide = "none"
      ) +
      scale_x_discrete(labels = MODEL_LABELS_SHORT) +
      scale_y_continuous(
        limits = summed_score_ylim,
        breaks = seq(summed_score_ylim[1], summed_score_ylim[2], by = 3),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(title = sc, x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(2, 4, 0, 4)
      )

    p_tile <- ggplot(d_sc, aes(x = Model, y = Metric, fill = Score)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_viridis_c(name = "Score", limits = c(0, 1)) +
      scale_x_discrete(labels = MODEL_LABELS_SHORT) +
      scale_y_discrete(labels = METRIC_LABELS) +
      labs(x = "Strategy", y = "Metric") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 13, face = "bold"),
        panel.grid = element_blank(),
        plot.margin = margin(0, 4, 2, 4)
      )

    p_bar / p_tile + plot_layout(heights = c(1.25, 5.25))
  }

  p <- wrap_plots(lapply(scen_levels, one_panel), ncol = ncol, guides = "collect") +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17))
    ) &
    theme(
      legend.position = "right",
      legend.text = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 14, face = "bold")
    )

  save_figure(p, filename, width = width, height = height)
}

scores_main <- make_scores(perf_main, SCENARIOS_MAIN)
scores_all_om <- make_scores(perf_all_om, SCENARIOS_ALL)

plot_score_heatmap(
  region_keep = "G",
  scores_long = scores_main$scores_long,
  score_totals = scores_main$score_totals,
  scenario_levels = SCENARIOS_MAIN,
  filename = "Figure_5_Summed_Scaled_Performance_System_6OM.png",
  title = "Summed scaled performance, system level"
)

plot_score_heatmap(
  region_keep = "R1",
  scores_long = scores_main$scores_long,
  score_totals = scores_main$score_totals,
  scenario_levels = SCENARIOS_MAIN,
  filename = "Figure_6_Summed_Scaled_Performance_Region1_6OM.png",
  title = "Summed scaled performance, region 1"
)

plot_score_heatmap(
  region_keep = "R2",
  scores_long = scores_main$scores_long,
  score_totals = scores_main$score_totals,
  scenario_levels = SCENARIOS_MAIN,
  filename = "Figure_7_Summed_Scaled_Performance_Region2_6OM.png",
  title = "Summed scaled performance, region 2"
)

plot_score_heatmap(
  region_keep = "G",
  scores_long = scores_all_om$scores_long,
  score_totals = scores_all_om$score_totals,
  scenario_levels = SCENARIOS_ALL,
  filename = "Figure_S12_Summed_Scaled_Performance_System_AllOM.png",
  title = "Summed scaled performance, system level, all OM scenarios",
  ncol = 4,
  width = 15,
  height = 17
)

plot_score_heatmap(
  region_keep = "R1",
  scores_long = scores_all_om$scores_long,
  score_totals = scores_all_om$score_totals,
  scenario_levels = SCENARIOS_ALL,
  filename = "Figure_S13_Summed_Scaled_Performance_Region1_AllOM.png",
  title = "Summed scaled performance, region 1, all OM scenarios",
  ncol = 4,
  width = 15,
  height = 17
)

plot_score_heatmap(
  region_keep = "R2",
  scores_long = scores_all_om$scores_long,
  score_totals = scores_all_om$score_totals,
  scenario_levels = SCENARIOS_ALL,
  filename = "Figure_S14_Summed_Scaled_Performance_Region2_AllOM.png",
  title = "Summed scaled performance, region 2, all OM scenarios",
  ncol = 4,
  width = 15,
  height = 17
)

cat(
  "DONE.\n",
  "Main scenarios: ", paste(SCENARIOS_MAIN, collapse = ", "), "\n",
  "All OM scenarios: ", paste(SCENARIOS_ALL, collapse = ", "), "\n",
  "Figure 3 and Figure 4 use relative differences from Equal allocation.\n",
  "Output folder: ", getwd(), "\n",
  sep = ""
)
