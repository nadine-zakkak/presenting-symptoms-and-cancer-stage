## Script name: main_figures.R
##
## Author: Nadine Zakkak
## 
## Purpose: Produce figures of symptom combinations for main manuscript

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "All main figures"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "flextable", "patchwork", "tidytext", "ggh4x")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)

# helper functions and global variables
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_functions_for_figures.R")
source("./scripts/common_scripts/00_global_vars.R")

run_date                  <- readRDS("./run_date.rds")
sx_fig_rename             <- read.delim("./lookup_files/sx_figures_rename.txt")
results_generated_dt_main <- configs$main_configs$results_generated_dt
results_loc_main          <- configs$main_configs$results_loc
results_loc_main_suppress <- str_glue("{results_loc_main}suppressed/")
results_loc_sa            <- paste0(configs$sa_configs$results_loc, "prostate_cancer/")
results_loc_sa_suppress   <- str_glue("{results_loc_sa}suppressed/")
results_generated_dt_sa   <- configs$sa_configs$results_generated_dt
save_dir                  <- str_glue("./tables_figs/manuscript/")
if(!dir.exists(save_dir)) dir.create(save_dir)
save_dir_suppress         <- str_glue("{save_dir}suppressed/")
if(!dir.exists(save_dir_suppress)) dir.create(save_dir_suppress)

# Stage distribution by symptom -----
## Without suppression -----
print(paste0("*****", current_script, ": Stage by symptom visuals*****"))

stage_by_sex_symptom <- read.csv(str_glue("{results_loc_main}stage_by_sex_symptom_{results_generated_dt_main}.csv"))
sx_order <- readRDS(str_glue("{results_loc_main}order_stage_by_sex_symptom_{results_generated_dt_main}.rds"))

stage_by_sex_symptom <- 
  stage_by_sex_symptom |> 
  label_sex() |>
  label_stage(binary = F) |>
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")

w <- plot_stage_sx(stage_by_sex_symptom, sx_order, "Women", stage_colours)
ggsave(plot     = w, 
       path     = save_dir,
       filename = str_glue("Fig2_stage_symptom_women_{run_date}.png"),
       width    = 30, 
       height   = 35, 
       unit     = "cm")

m <- plot_stage_sx(stage_by_sex_symptom, sx_order, "Men", stage_colours)
ggsave(plot     = m, 
       path     = save_dir,
       filename = str_glue("Fig3_stage_symptom_men_{run_date}.png"),
       width    = 30, 
       height   = 35, 
       unit     = "cm")

print(paste0("*****", current_script, ": Stage by symptom and cancer site visuals*****"))
## With suppression -----
print(paste0("*****", current_script, ": Stage by symptom visuals*****"))

stage_by_sex_symptom_suppress <- read.csv(str_glue("{results_loc_main_suppress}stage_by_sex_symptom_{results_generated_dt_main}.csv"))
sx_order <- readRDS(str_glue("{results_loc_main}order_stage_by_sex_symptom_{results_generated_dt_main}.rds"))

stage_by_sex_symptom_suppress <-
  stage_by_sex_symptom_suppress |>
  select(sex, symptom, stage_recode,
         n_tumour   = suppressed_n_tidy,
         N_tumour   = suppressed_N,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)

stage_by_sex_symptom_suppress <- 
  stage_by_sex_symptom_suppress |> 
  label_sex() |>
  label_stage(binary = F) |>
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")

w_suppress <- plot_stage_sx(stage_by_sex_symptom_suppress, sx_order, "Women", stage_colours)
ggsave(plot     = w_suppress, 
       path     = save_dir_suppress,
       filename = str_glue("Fig2_stage_symptom_women_{run_date}.png"),
       width    = 30, 
       height   = 35, 
       unit     = "cm")

m_suppress <- plot_stage_sx(stage_by_sex_symptom_suppress, sx_order, "Men", stage_colours)
ggsave(plot     = m_suppress, 
       path     = save_dir_suppress,
       filename = str_glue("Fig3_stage_symptom_men_{run_date}.png"),
       width    = 30, 
       height   = 35, 
       unit     = "cm")

print(paste0("*****", current_script, ": Stage by symptom and cancer site visuals*****"))
# Stage distribution by symptom, stratified by cancer site -----
## Without suppression -----
stage_bin_by_symptom_site_exemplar <- read.csv(str_glue("{results_loc_main}stage_bin_by_symptom_site_exemplar_{results_generated_dt_main}.csv"))

p_all <- plot_stage_sx_site(stage_bin_by_symptom_site_exemplar, y_indent = 1.5, stagebin_colours)
p_redflags    <- p_all$redflags
p_nonspecific <- p_all$nonspecific
p_general     <- p_all$general
p_combined <- 
  p_redflags$Lung + p_redflags$Colon + p_redflags$Rectum + p_redflags$Oesophagus + p_redflags$Kidney + p_redflags$Pancreas +
  p_nonspecific$Lung + p_nonspecific$Colon + p_nonspecific$Rectum + p_nonspecific$Oesophagus + p_nonspecific$Kidney + p_nonspecific$Pancreas +
  p_general$rf + p_general$nonsp +
  plot_layout(design = c(
    area(t = 1, l = 1),
    area(t = 1, l = 2),
    area(t = 1, l = 3),
    area(t = 1, l = 4),
    area(t = 1, l = 5),
    area(t = 1, l = 6),
    area(t = 2, l = 1, b = 3),
    area(t = 2, l = 2, b = 3),
    area(t = 2, l = 3, b = 3),
    area(t = 2, l = 4, b = 3),
    area(t = 2, l = 5, b = 3),
    area(t = 2, l = 6, b = 3),
    area(t = 1, l = 1, b = 1, r = 6),
    area(t = 2, l = 1, b = 3, r = 6)
  )
  )

ggsave(plot     = p_combined, 
       path     = save_dir,
       filename = str_glue("Fig4_exemplar_rf_nonsp_{run_date}.png"), 
       width    = 35, 
       height   = 15, 
       unit     = "cm")

## With suppression -----
stage_bin_by_symptom_site_exemplar_suppress <- read.csv(str_glue("{results_loc_main_suppress}stage_bin_by_symptom_site_exemplar_{results_generated_dt_main}.csv"))

stage_bin_by_symptom_site_exemplar_suppress <-
  stage_bin_by_symptom_site_exemplar_suppress |>
  select(symptom, cancer_site_desc, stage_bin,
         n_tumour = suppressed_n_tidy,
         N_tumour = suppressed_N,
         prop100 = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)

p_all_suppress <- plot_stage_sx_site(stage_bin_by_symptom_site_exemplar_suppress, y_indent = 1.5, stagebin_colours)
p_redflags_suppress    <- p_all_suppress$redflags
p_nonspecific_suppress <- p_all_suppress$nonspecific
p_general_suppress     <- p_all_suppress$general
p_combined_suppress <- 
  p_redflags_suppress$Lung + p_redflags_suppress$Colon + p_redflags_suppress$Rectum + p_redflags_suppress$Oesophagus + p_redflags_suppress$Kidney + p_redflags_suppress$Pancreas +
  p_nonspecific_suppress$Lung + p_nonspecific_suppress$Colon + p_nonspecific_suppress$Rectum + p_nonspecific_suppress$Oesophagus + p_nonspecific_suppress$Kidney + p_nonspecific_suppress$Pancreas +
  p_general_suppress$rf + p_general_suppress$nonsp +
  plot_layout(design = c(
    area(t = 1, l = 1),
    area(t = 1, l = 2),
    area(t = 1, l = 3),
    area(t = 1, l = 4),
    area(t = 1, l = 5),
    area(t = 1, l = 6),
    area(t = 2, l = 1, b = 3),
    area(t = 2, l = 2, b = 3),
    area(t = 2, l = 3, b = 3),
    area(t = 2, l = 4, b = 3),
    area(t = 2, l = 5, b = 3),
    area(t = 2, l = 6, b = 3),
    area(t = 1, l = 1, b = 1, r = 6),
    area(t = 2, l = 1, b = 3, r = 6)
  )
  )

ggsave(plot     = p_combined_suppress, 
       path     = save_dir_suppress,
       filename = str_glue("Fig4_exemplar_rf_nonsp_{run_date}.png"), 
       width    = 35, 
       height   = 15, 
       unit     = "cm")

# Blood test visuals ----
## without suppression ----
stage_bin_by_sex_symptom_bt  <- read.csv(str_glue("{results_loc_main}stage_bin_by_sex_symptom_bt_{results_generated_dt_main}.csv")) 
t_test                       <- read.csv(str_glue("{results_loc_main}t_test_bt_{results_generated_dt_main}.csv")) 

# Read in results from supplementary analysis (no prostate)
stage_bin_by_sex_symptom_bt_sa  <- read.csv(paste0(results_loc_sa, "stage_bin_by_sex_symptom_bt_", results_generated_dt_sa, ".csv")) 

# Only keep rows with "significant" p-values
vis_df <- 
  stage_bin_by_sex_symptom_bt |>
  inner_join(t_test, by = c('sex_lb', 'symptom')) |>
  filter(p <= 0.05) |>
  mutate(Analysis = 'All cancers') |>
  bind_rows(
    stage_bin_by_sex_symptom_bt_sa |>
      inner_join(t_test, by = c('sex_lb', 'symptom')) |>
      filter(p <= 0.05) |>
      mutate(Analysis = 'All cancers excl. prostate cancer')
  )

vis_df <-
  vis_df |>
  select(sex_lb, symptom, blood_test, stage_bin, n_tumour, N_tumour, prop100, prop100_lb, prop100_ub, p, Analysis) |>
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")


# women - plot distribution of binary stage by blood test and symptom
w_opt2 <- plot_stage_sx_bt(vis_df, "Women", stagebin_colours)
ggsave(plot     = w_opt2,
       path     = save_dir,
       filename = str_glue("Fig5_stage_bin_bt_women_{run_date}.png"),
       width    = 26, 
       height   = 15, 
       unit     = "cm")

# men - plot distribution of binary stage by blood test and symptom
m_opt2 <- plot_stage_sx_bt(vis_df, "Men", stagebin_colours)
ggsave(plot     = m_opt2, 
       path     = save_dir,
       filename = str_glue("Fig6_stage_bin_bt_men_{run_date}.png"),
       width    = 26, 
       height   = 10, 
       unit = "cm")

## with suppression ----
stage_bin_by_sex_symptom_bt_suppress  <- read.csv(str_glue("{results_loc_main_suppress}stage_bin_by_sex_symptom_bt_{results_generated_dt_main}.csv")) 
t_test                       <- read.csv(str_glue("{results_loc_main}t_test_bt_{results_generated_dt_main}.csv")) 

# Read in results from supplementary analysis (no prostate)
stage_bin_by_sex_symptom_bt_sa_suppress  <- read.csv(str_glue("{results_loc_sa_suppress}stage_bin_by_sex_symptom_bt_{results_generated_dt_sa}.csv")) 

stage_bin_by_sex_symptom_bt_suppress <-
  stage_bin_by_sex_symptom_bt_suppress |>
  select(sex, sex_lb, symptom, inv, blood_test, stage_bin,
         n_tumour   = suppressed_n_tidy,
         N_tumour   = suppressed_N,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)

stage_bin_by_sex_symptom_bt_sa_suppress <-
  stage_bin_by_sex_symptom_bt_sa_suppress |>
  select(sex, sex_lb, symptom, inv, blood_test, stage_bin,
         n_tumour   = suppressed_n_tidy,
         N_tumour   = suppressed_N,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)

# Only keep rows with "significant" p-values
vis_df_suppress <- 
  stage_bin_by_sex_symptom_bt_suppress |>
  inner_join(t_test, by = c('sex_lb', 'symptom')) |>
  filter(p <= 0.05) |>
  mutate(Analysis = 'All cancers') |>
  mutate(Analysis = if_else(sex_lb == "Men", paste("(a)", Analysis), Analysis)) |>
  bind_rows(
    stage_bin_by_sex_symptom_bt_sa_suppress |>
      inner_join(t_test, by = c('sex_lb', 'symptom')) |>
      filter(p <= 0.05) |>
      mutate(Analysis = 'All cancers excl. prostate cancer') |>
      mutate(Analysis = if_else(sex_lb == "Men", paste("(b)", Analysis), Analysis))
  )

vis_df_suppress <-
  vis_df_suppress |>
  select(sex_lb, symptom, blood_test, stage_bin, n_tumour, N_tumour, prop100, prop100_lb, prop100_ub, p, Analysis) |>
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")

# women - plot distribution of binary stage by blood test and symptom
w_opt2_suppress <- plot_stage_sx_bt(vis_df_suppress, "Women", stagebin_colours)
ggsave(plot     = w_opt2_suppress,
       path     = save_dir_suppress,
       filename = str_glue("Fig5_stage_bin_bt_women_{run_date}.png"),
       width    = 26, 
       height   = 15, 
       unit     = "cm")

# men - plot distribution of binary stage by blood test and symptom
m_opt2_suppress <- plot_stage_sx_bt(vis_df_suppress, "Men", stagebin_colours)
ggsave(plot     = m_opt2_suppress, 
       path     = save_dir_suppress,
       filename = str_glue("Fig6_stage_bin_bt_men_{run_date}.png"),
       width    = 26, 
       height   = 10, 
       unit = "cm")

# End script ----
print(str_glue("***** {current_script}: files saved at {save_dir} ******"))
print(str_glue("***** {current_script}: detach packages ******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
