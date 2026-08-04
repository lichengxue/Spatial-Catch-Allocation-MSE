library(dplyr)
library(tidyr)
library(ggplot2)

setwd("C:/Users/liche/Desktop/Catch_Paper")

mods_file <- "Mods_NB.RDS"   # NBL raw scenario
scenario_label <- "NBL"
realization_id <- 1

strategy_names <- c("Equal", "Gear-Equal", "Fleet-Catch", "Index-Equal", "Index-Gear")

strategy_cols <- c(
  "Equal"       = "#2874A6",
  "Gear-Equal"  = "#F39C12",
  "Fleet-Catch" = "#27AE60",
  "Index-Equal" = "#E74C3C",
  "Index-Gear"  = "#8E44AD"
)

feedback_start <- 36
fallback_start_year <- 1989

mods <- readRDS(mods_file)
mod_list <- mods[[realization_id]][seq_along(strategy_names)]

extract_om_trajectory <- function(mod_obj, strategy_name) {
  rep <- mod_obj$om$rep
  input <- mod_obj$om$input$data
  
  n_regions <- input$n_regions[1]
  nyear <- nrow(rep$SSB)
  
  years <- if (!is.null(input$years)) {
    input$years
  } else {
    seq(fallback_start_year, by = 1, length.out = nyear)
  }
  
  ssb_reg <- as.matrix(rep$SSB)
  colnames(ssb_reg) <- paste0("Region ", seq_len(n_regions))
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
  
  catch_reg <- as.matrix(catch_reg)
  colnames(catch_reg) <- paste0("Region ", seq_len(n_regions))
  catch_system <- rowSums(catch_reg)
  
  ssb_df <- as.data.frame(cbind(ssb_reg, System = ssb_system)) %>%
    mutate(Year_index = seq_len(nyear), Year = years, Metric = "SSB")
  
  catch_df <- as.data.frame(cbind(catch_reg, System = catch_system)) %>%
    mutate(Year_index = seq_len(nyear), Year = years, Metric = "Catch")
  
  bind_rows(ssb_df, catch_df) %>%
    pivot_longer(
      cols = c(starts_with("Region"), System),
      names_to = "Spatial_scale",
      values_to = "Value"
    ) %>%
    mutate(
      Scenario = scenario_label,
      Realization = realization_id,
      Strategy = strategy_name,
      Period = if_else(Year_index < feedback_start, "Historical", "Feedback")
    )
}

traj_df <- bind_rows(lapply(seq_along(mod_list), function(i) {
  extract_om_trajectory(mod_list[[i]], strategy_names[i])
}))

plot_df <- traj_df %>%
  mutate(
    Strategy = factor(Strategy, levels = strategy_names),
    Spatial_scale = factor(Spatial_scale, levels = c("Region 1", "Region 2", "System")),
    Metric = factor(Metric, levels = c("Catch", "SSB"))
  )

historical_df <- plot_df %>%
  filter(Year_index < feedback_start, Strategy == "Equal")

feedback_df <- plot_df %>%
  filter(Year_index >= feedback_start)

feedback_x <- unique(plot_df$Year[plot_df$Year_index == feedback_start])[1]

p <- ggplot() +
  geom_line(
    data = historical_df,
    aes(x = Year, y = Value, group = interaction(Metric, Spatial_scale)),
    color = "black",
    linewidth = 0.85,
    alpha = 0.95
  ) +
  geom_vline(
    xintercept = feedback_x,
    linetype = "dashed",
    color = "gray45",
    linewidth = 0.6
  ) +
  geom_line(
    data = feedback_df,
    aes(x = Year, y = Value, color = Strategy),
    linewidth = 0.85,
    alpha = 0.72
  ) +
  facet_grid(Metric ~ Spatial_scale, scales = "free_y") +
  scale_color_manual(values = strategy_cols, drop = FALSE) +
  labs(
    x = "Year",
    y = NULL,
    color = "Strategy"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray88", linewidth = 0.35)
  )

print(p)

ggsave(
  filename = paste0("NBL_realization_", realization_id, "_OM_Catch_SSB_Trajectory_actual_years.png"),
  plot = p,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  filename = paste0("NBL_realization_", realization_id, "_OM_Catch_SSB_Trajectory_actual_years.pdf"),
  plot = p,
  width = 11,
  height = 7
)