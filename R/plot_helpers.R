# Plotting helpers for manuscript figures.
# Scoring aggregation is defined in reproduce_figures.R, using raw means.
MODEL_LEVELS <- c("Equal", "Gear-Equal", "Fleet-Catch", "Index-Equal", "Index-Gear")

BASELINE_MODEL <- "Equal"

REGION_LEVELS <- c("R1", "R2", "G")

SCENARIOS_MAIN_RAW <- c("NB", "MBH", "MU", "MUH", "NBC", "NUC")

SCENARIOS_MAIN <- c("NBL", "MBH", "MUL", "MUH", "NBLD", "NULD")

SCENARIOS_ALL_RAW <- c("NB", "NBH", "NU", "NUH", "NBI", "NBO", "NBC", "NUC", "NBHS", "MB", "MBH", "MU",
    "MUH", "MBI", "MBO", "MBHS")

SCENARIOS_ALL <- c("NBL", "NBH", "NUL", "NUH", "NBLI", "NBLO", "NBLD", "NULD", "NBH-BH", "MBL", "MBH",
    "MUL", "MUH", "MBLI", "MBLO", "MBH-BH")

SCENARIOS_MAIN_Y <- rev(SCENARIOS_MAIN)

SCENARIOS_ALL_Y <- rev(SCENARIOS_ALL)

STRATEGY_COLORS <- c(Equal = "#2874A6", `Gear-Equal` = "#F39C12", `Fleet-Catch` = "#27AE60", `Index-Equal` = "#E74C3C",
    `Index-Gear` = "#8E44AD")

REGION_LABELS <- c(R1 = "Region 1", R2 = "Region 2", G = "System")

MODEL_LABELS_SHORT <- c(Equal = "Equal", `Gear-Equal` = "Gear-Equal", `Fleet-Catch` = "Fleet-Catch",
    `Index-Equal` = "Index-Equal", `Index-Gear` = "Index-Gear")

STRATEGY_DODGE_WIDTH <- 0.82

recode_scenario <- function(x) {
    dplyr::recode(as.character(x), NB = "NBL", MB = "MBL", NU = "NUL", MU = "MUL", NBI = "NBLI", MBI = "MBLI",
        NBO = "NBLO", MBO = "MBLO", NBC = "NBLD", NUC = "NULD", NBHS = "NBH-BH", MBHS = "MBH-BH", .default = as.character(x))
}

theme_manuscript <- function() {
    theme_bw() + theme(legend.position = "right", axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"), axis.title.x = element_text(size = 16,
            face = "bold", margin = margin(t = 10)), axis.title.y = element_text(size = 16, face = "bold",
            margin = margin(r = 10)), plot.title = element_text(size = 17, face = "bold"), strip.text = element_text(size = 14,
            face = "bold"), legend.text = element_text(size = 13, face = "bold"), legend.title = element_text(size = 14,
            face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

make_scenario_bands <- function(scenario_levels, fill = "grey95") {
    idx <- seq_along(scenario_levels)
    band_levels <- scenario_levels[idx%%2 == 1]
    data.frame(Scenario = factor(band_levels, levels = scenario_levels), xmin = -Inf, xmax = Inf, ymin = as.numeric(factor(band_levels,
        levels = scenario_levels)) - 0.5, ymax = as.numeric(factor(band_levels, levels = scenario_levels)) +
        0.5, fill = fill)
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
    vals <- if (include_iqr)
        c(sum_dat$Q1, sum_dat$Q3, sum_dat$Median)
    else sum_dat$Median
    max_abs <- suppressWarnings(max(abs(vals), na.rm = TRUE))
    if (!is.finite(max_abs) || max_abs == 0)
        max_abs <- minimum
    c(-max_abs * pad, max_abs * pad)
}

make_metric_summary <- function(perf, value_col, metric_prefix, scenario_levels_y = SCENARIOS_MAIN_Y,
    period_keep = c("short", "long"), period_suffix = c(short = "ST", long = "LT", all = "All"), drop_equal = FALSE) {
    model_keep <- if (drop_equal)
        setdiff(MODEL_LEVELS, BASELINE_MODEL)
    else MODEL_LEVELS
    metric_levels <- paste0(metric_prefix, "_", unname(period_suffix[period_keep]))
    perf %>% filter(Model %in% model_keep, as.character(Period) %in% period_keep) %>% mutate(Scenario = factor(as.character(Scenario),
        levels = scenario_levels_y), Model = factor(as.character(Model), levels = model_keep), Region = factor(as.character(Region),
        levels = REGION_LEVELS), Metric = paste0(metric_prefix, "_", period_suffix[as.character(Period)]),
        Value = .data[[value_col]]) %>% filter(!is.na(Value), !is.na(Metric)) %>% mutate(Metric = factor(Metric,
        levels = metric_levels)) %>% group_by(Scenario, Region, Model, Metric) %>% summarise(Median = median(Value,
        na.rm = TRUE), Q1 = quantile(Value, 0.25, na.rm = TRUE), Q3 = quantile(Value, 0.75, na.rm = TRUE),
        .groups = "drop")
}

plot_metric_iqr <- function(sum_dat, filename, xlab, metric_labels, xlim = NULL, auto_xlim = TRUE, facet_scales = "fixed",
    percent_axis = FALSE, ratio_axis = FALSE, add_zero_line = FALSE, show_iqr = TRUE, clip_to_xlim = FALSE,
    median_segment_from = NULL, median_segment_linewidth = 1.15, median_segment_alpha = 0.28, width = 14,
    height = 9, title = NULL) {
    scen_levels <- levels(sum_dat$Scenario)
    model_levels <- levels(droplevels(sum_dat$Model))
    sum_dat <- sum_dat %>% mutate(Model_dodge = factor(as.character(Model), levels = rev(model_levels)))
    band_df <- make_scenario_bands(scen_levels)
    dodge <- position_dodge(width = STRATEGY_DODGE_WIDTH)
    if (is.null(xlim) && auto_xlim && facet_scales == "fixed") {
        xlim <- auto_symmetric_xlim(sum_dat, include_iqr = show_iqr)
    }
    x_range <- if (!is.null(xlim))
        diff(xlim)
    else NA_real_
    if (clip_to_xlim && !is.null(xlim)) {
        sum_dat <- sum_dat %>% mutate(Median_raw = Median, Q1_raw = Q1, Q3_raw = Q3, Median = pmin(pmax(Median,
            xlim[1]), xlim[2]), Q1 = pmin(pmax(Q1, xlim[1]), xlim[2]), Q3 = pmin(pmax(Q3, xlim[1]), xlim[2]),
            clipped_low = Median_raw < xlim[1] | Q1_raw < xlim[1] | Q3_raw < xlim[1], clipped_high = Median_raw >
                xlim[2] | Q1_raw > xlim[2] | Q3_raw > xlim[2])
    }
    else {
        sum_dat <- sum_dat %>% mutate(clipped_low = FALSE, clipped_high = FALSE)
    }
    p <- ggplot(sum_dat, aes(y = Scenario, color = Model, group = Model_dodge)) + geom_rect(data = band_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey95",
        color = NA)
    if (show_iqr) {
        p <- p + geom_errorbarh(aes(xmin = Q1, xmax = Q3), position = dodge, height = 0.22, linewidth = 0.35,
            alpha = 0.35)
    }
    if (!is.null(median_segment_from)) {
        p <- p + geom_errorbarh(aes(xmin = median_segment_from, xmax = Median), position = dodge, height = 0,
            linewidth = median_segment_linewidth, alpha = median_segment_alpha)
    }
    p <- p + geom_point(aes(x = Median), position = dodge, size = 2.7)
    if (clip_to_xlim && !is.null(xlim)) {
        p <- p + geom_text(data = sum_dat %>% filter(clipped_low), aes(x = xlim[1] + 0.015 * x_range,
            label = "<"), position = dodge, fontface = "bold", size = 3.5, show.legend = FALSE) + geom_text(data = sum_dat %>%
            filter(clipped_high), aes(x = xlim[2] - 0.015 * x_range, label = ">"), position = dodge,
            fontface = "bold", size = 3.5, show.legend = FALSE)
    }
    p <- p + facet_grid(Metric ~ Region, scales = facet_scales, labeller = labeller(Metric = metric_labels,
        Region = REGION_LABELS)) + scale_color_manual(values = STRATEGY_COLORS[model_levels], breaks = model_levels,
        labels = MODEL_LABELS_SHORT[model_levels], name = "Strategy") + labs(title = title, x = xlab,
        y = "Operating model scenario") + theme_manuscript()
    if (!is.null(xlim)) {
        p <- p + coord_cartesian(xlim = xlim)
    }
    if (add_zero_line) {
        p <- p + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, color = "grey35")
    }
    if (ratio_axis) {
        ratio_breaks <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10, 25, 50, 100, 250, 500)
        p <- p + scale_x_continuous(breaks = log(ratio_breaks), labels = paste0(ratio_breaks, "x"))
    }
    else if (percent_axis) {
        p <- p + scale_x_continuous(labels = percent_format(accuracy = 1))
    }
    save_figure(p, filename, width = width, height = height)
}

METRIC_LEVELS <- c("Catch_ST", "SSB_ST", "Catch_LT", "SSB_LT", "AAV_Catch_All", "AAV_SSB_All", "p_overfishing_All",
    "p_overfished_All")

METRIC_LABELS <- c(Catch_ST = "Catch ST", SSB_ST = "SSB ST", Catch_LT = "Catch LT", SSB_LT = "SSB LT",
    AAV_Catch_All = "AAV catch", AAV_SSB_All = "AAV SSB", p_overfishing_All = "P(overfishing)", p_overfished_All = "P(overfished)")

GOOD_METRICS <- c("Catch_ST", "SSB_ST", "Catch_LT", "SSB_LT")

SUMMED_SCORE_MAX <- length(METRIC_LEVELS)

AVERAGE_SCORE_MAX <- 1

score_metric_0_1 <- function(x, metric_name) {
    x <- as.numeric(x)
    metric_name <- as.character(metric_name)[1]
    if (all(is.na(x)))
        return(rep(NA_real_, length(x)))
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
        return((x - min(x, na.rm = TRUE))/rng)
    }
    (max(x, na.rm = TRUE) - x)/rng
}

plot_score_heatmap <- function(region_keep, scores_long, score_totals, scenario_levels, filename, title,
    ncol = 3, width = 12.5, height = 10.5, score_col = "Summed_score", score_ylim = c(0, SUMMED_SCORE_MAX),
    score_breaks = seq(score_ylim[1], score_ylim[2], by = 3)) {
    dat <- scores_long %>% filter(Region == region_keep) %>% mutate(Scenario = factor(as.character(Scenario),
        levels = scenario_levels), Metric = factor(as.character(Metric), levels = rev(METRIC_LEVELS)),
        Model = factor(as.character(Model), levels = MODEL_LEVELS))
    totals <- score_totals %>% filter(Region == region_keep) %>% mutate(Scenario = factor(as.character(Scenario),
        levels = scenario_levels), Model = factor(as.character(Model), levels = MODEL_LEVELS))
    scen_levels <- scenario_levels
    n_panels <- length(scen_levels)
    n_rows <- ceiling(n_panels/ncol)
    one_panel <- function(sc, panel_index) {
        panel_row <- ceiling(panel_index/ncol)
        panel_col <- ((panel_index - 1)%%ncol) + 1
        is_bottom_row <- panel_row == n_rows
        is_left_col <- panel_col == 1
        d_sc <- dat %>% filter(Scenario == sc)
        t_sc <- totals %>% filter(Scenario == sc) %>% mutate(Winner = .data[[score_col]] == max(.data[[score_col]],
            na.rm = TRUE))
        p_bar <- ggplot(t_sc, aes(x = Model, y = .data[[score_col]], fill = Winner)) + geom_col(color = "black",
            linewidth = 0.2) + scale_fill_manual(values = c(`FALSE` = "grey82", `TRUE` = "#2A9D8F"),
            guide = "none") + scale_x_discrete(labels = MODEL_LABELS_SHORT) + scale_y_continuous(limits = score_ylim,
            breaks = score_breaks, expand = expansion(mult = c(0, 0.02))) + labs(title = sc, x = NULL,
            y = NULL) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold",
            size = 13), axis.text.x = element_blank(), axis.text.y = if (is_left_col)
            element_text(size = 10)
        else element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            plot.margin = margin(2, 4, 0, 4))
        p_tile <- ggplot(d_sc, aes(x = Model, y = Metric, fill = Score)) + geom_tile(color = "white",
            linewidth = 0.25) + scale_fill_viridis_c(name = "Score", limits = c(0, 1)) + scale_x_discrete(labels = MODEL_LABELS_SHORT) +
            scale_y_discrete(labels = METRIC_LABELS) + labs(x = "Strategy", y = "Metric") + theme_minimal() +
            theme(axis.text.x = if (is_bottom_row) {
                element_text(size = 11, angle = 45, hjust = 1, face = "bold")
            }
            else {
                element_blank()
            }, axis.text.y = if (is_left_col)
                element_text(size = 11)
            else element_blank(), axis.title.x = if (is_bottom_row) {
                element_text(size = 13, face = "bold")
            }
            else {
                element_blank()
            }, axis.title.y = if (is_left_col) {
                element_text(size = 13, face = "bold")
            }
            else {
                element_blank()
            }, panel.grid = element_blank(), plot.margin = margin(0, 4, 2, 4))
        p_bar/p_tile + plot_layout(heights = c(1.25, 5.25))
    }
    p <- wrap_plots(lapply(seq_along(scen_levels), function(i) one_panel(scen_levels[i], i)), ncol = ncol,
        guides = "collect") + plot_annotation(title = title, theme = theme(plot.title = element_text(hjust = 0.5,
        face = "bold", size = 17))) & theme(legend.position = "right", legend.text = element_text(size = 13,
        face = "bold"), legend.title = element_text(size = 14, face = "bold"))
    save_figure(p, filename, width = width, height = height)
}
