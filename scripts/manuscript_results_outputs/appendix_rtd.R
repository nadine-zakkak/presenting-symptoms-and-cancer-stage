## --------------------------------
## Script name: appendix_rtd.R
## Author: Nadine Zakkak
## Purpose: Prep figures for sensitivity stage rtd-symptom analysis
## ---------------------------------
# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "appendix_rtd.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "flextable", "scales", "patchwork")
sapply(pckgs, function(pckg)  { load_package(pckg) })

run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_for_tables.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")
sx_fig_rename <- read.delim("./lookup_files/sx_figures_rename.txt")

# Locations where to read data from 
results_loc_main <- configs$main_configs$results_loc
results_generated_dt_main   <- configs$main_configs$results_generated_dt
results_loc_main_suppress <- str_glue("{results_loc_main}suppressed/")

results_loc <- paste0(configs$sa_configs$results_loc, str_split(sa_rtd_loc, pattern = "/")[[1]][2], "/")
results_loc_suppress <- str_glue("{results_loc}/suppressed/")
results_loc_dt <- configs$sa_configs$results_generated_dt

# Directory to save utputs
save_dir    <- str_glue("./tables_figs/manuscript/")
if(!dir.exists(save_dir)) dir.create(save_dir)
save_dir_suppress         <- str_glue("{save_dir}suppressed/")

# Load data ------
N_tumour_by_sex_symptom_rtd_suppress <- read.csv(str_glue("{results_loc_suppress}n_tumour_by_sex_symptom_rtd_{results_loc_dt}.csv"))
stage_bin_by_sex_symptom_rtd_suppress <- read.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_rtd_{results_loc_dt}.csv"))
sx_order <- readRDS(str_glue("{results_loc_main}order_stage_by_sex_symptom_{results_generated_dt_main}.rds"))
order_short <- setNames(sx_fig_rename$New, sx_fig_rename$Original)
order_m       <- sx_order[["Men"]]
final_m       <- order_short[order_m]
final_m       <- c("N/A", "N/K", setdiff(final_m, c("N/A", "N/K")))
order_w       <- sx_order[["Women"]]
final_w       <- order_short[order_w]
final_w       <- c("N/A", "N/K", setdiff(final_w, c("N/A", "N/K")))


prep_data <- function(data) {
  prepped_data <-
    data |> 
    label_sex() |>
    left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom") |>
    arrange(desc(symptom_short))
  
  # only keep relevant symptoms
  ## Men
  prepped_data_m <- prepped_data |> filter(sex_lb == "Men", symptom_short %in% final_m)
  ## Women
  prepped_data_w <- prepped_data |> filter(sex_lb == "Women", symptom_short %in% final_w)
  
  # join back together
  prepped_data_all <- prepped_data_m |> bind_rows(prepped_data_w)
  
  prepped_data_all |>
    mutate(symptom_short = factor(symptom_short, levels = c("N/A", "N/K", setdiff(prepped_data$symptom_short, c("N/A", "N/K")))),
           typereferral_recode = factor(typereferral_recode, levels = referral_order)
    ) |>
    filter(N_tumour >= 100) |>
    mutate(y_label = str_glue("N={N_tumour}"))
}

plot_sx_by_all_rtd <- function(data) {
  data |>
    ggplot() +
    geom_bar(aes(y    = prop100, 
                 x    = symptom_short, 
                 fill = reorder(typereferral_recode , desc(typereferral_recode))),
             stat  = "identity", 
             width = .5) +
    ylim(-12, 101) +
    geom_text(aes(x = symptom_short, y = -12, label = y_label),
              size   = 3.5, 
              hjust  = 0, 
              colour = "black") +
    scale_fill_manual(values = rtd_colours) +
    guides(fill = guide_legend(reverse = T, title = "Route to diagnosis")) +
    coord_flip() +
    facet_grid(~sex_lb) +
    theme_minimal() +
    theme(text               = element_text(size = 14, colour = "black"),
          axis.text.y        = element_text(size = 14, colour = "black"),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    xlab("Symptom") +
    ylab("Proportion (%)")
}

plot_data_rtd_by_sx  <- N_tumour_by_sex_symptom_rtd_suppress  |> prep_data()
p_rtd_sx <- plot_data_rtd_by_sx  |> plot_sx_by_all_rtd()
ggsave(plot     = p_rtd_sx, 
       path     = save_dir_suppress,
       filename = str_glue("App_rtd_sx_{run_date}.png"),
       width    = 38, 
       height   = 35, 
       unit     = "cm")

# Stratify by stage -----
stage_bin_by_sex_symptom_suppress <- read.csv(str_glue("{results_loc_main_suppress}stage_bin_by_sex_symptom_{results_generated_dt_main}.csv"))

stage_bin_by_sex_symptom_suppress <- 
  stage_bin_by_sex_symptom_suppress |> 
  label_sex() |>
  label_stage(binary = T) |>
  mutate(typereferral_recode = "All") |>
  select(sex_lb, symptom, stage_bin_lb, prop100_suppressed)
stage_bin_by_sex_symptom_suppress <- stage_bin_by_sex_symptom_suppress |>  mutate(typereferral_recode = "All")

# only visualise if N>=100
plot_data <- 
  stage_bin_by_sex_symptom_suppress |>
  bind_rows(
    stage_bin_by_sex_symptom_rtd_suppress |>
      filter(suppressed_N >= 100) |>
      label_stage(binary = T) |>
      select(sex_lb, symptom, stage_bin_lb, prop100_suppressed, typereferral_recode)
  ) 

plot_data <-
  plot_data |> 
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")

plot_rtd_sex_sx_stg_fct <- function(data, filter_by_sex, sx_order) {
  data_filtered <-   data |>
    filter(sex_lb == filter_by_sex) |>
    filter(stage_bin_lb  == "1-3") |>
    group_by(symptom_short) |>
    filter(n_distinct(typereferral_recode) > 2) |>
    ungroup()
  
  data_filtered |>
    mutate(symptom_short = factor(symptom_short, levels = rev(sx_order))) |>
    mutate(typereferral_recode = factor(typereferral_recode, levels = c("All", referral_order))) |>
    arrange(symptom_short, typereferral_recode) |>
    ggplot() +
    geom_bar(aes(y    = prop100_suppressed, 
                 x    = typereferral_recode),
             stat     = "identity", 
             colour   = "#4292C6",
             fill     = NA,
             width    = .5) +
    facet_wrap( ~ symptom_short, nrow = 3, scales = "free") +
    ylim(0, 100) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 100))
    theme_classic() +
    theme(text               = element_text(size = 14, colour = "black"),
          axis.text.y        = element_text(size = 15, colour = "black"),
          axis.text.x        = element_text(size = 15, angle = 45, hjust = 1),
          strip.text         = element_text(size = 14, colour = "black")) +
    xlab("Route to diagnosis") +
    ylab("% Stage 1-3") +
    labs(title = filter_by_sex)
}

p_rtd_stage_w <- plot_data |> plot_rtd_sex_sx_stg_fct("Women", final_w)
ggsave(plot     = p_rtd_stage_w, 
       path     = save_dir_suppress,
       filename = str_glue("App_rtd_sx_w_{run_date}.png"),
       width    = 45, 
       height   = 25, 
       unit     = "cm")
p_rtd_stage_m <- plot_data |> plot_rtd_sex_sx_stg_fct("Men", final_m)
ggsave(plot     = p_rtd_stage_m, 
       path     = save_dir_suppress,
       filename = str_glue("App_rtd_sx_m_{run_date}.png"),
       width    = 45, 
       height   = 25, 
       unit     = "cm")

