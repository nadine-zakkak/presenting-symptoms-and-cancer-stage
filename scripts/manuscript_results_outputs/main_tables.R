## Script name: main_figures.R
##
## Author: Nadine Zakkak
## 
## Purpose: Produce tables
##  of symptom combinations for main manuscript

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "All main tables"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "flextable", "tidytext")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)

# helper functions and global variables
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_functions_for_tables.R")
source("./scripts/common_scripts/00_global_vars.R")

run_date                  <- readRDS("./run_date.rds")
sx_fig_rename             <- read.delim("./lookup_files/sx_figures_rename.txt")
results_generated_dt_main <- configs$main_configs$results_generated_dt
results_loc_main          <- configs$main_configs$results_loc
results_loc_main_suppress <- str_glue("{results_loc_main}suppressed/")
save_dir                  <- str_glue("./tables_figs/manuscript/")
if(!dir.exists(save_dir)) dir.create(save_dir)
save_dir_suppress         <- str_glue("{save_dir}suppressed/")
if(!dir.exists(save_dir_suppress)) dir.create(save_dir_suppress)

print(paste0("*****", current_script, ": Charactestics by all patients*****"))
## Descriptive - Characteristics of all patients -----
# load data
summary_total <- read.csv(str_glue("{results_loc_main}summary_all_{results_generated_dt_main}.csv"))
# Total nuber of patients
total_numb_pat <- 
  data.frame(
    var   = "Total",
    n_patients = unique(summary_total$N_pat)
  )

# format and tidy table to export into word
summary_total |>
  mutate(across(starts_with("prop100"), ~format_numb(.x))) |>
  format_ci() |>
  bind_rows(total_numb_pat) |>
  mutate(var = factor(var),
         var = factor(var, levels = c("Total", "Sex", "Age", "IMD"))) |>
  arrange(var) |>
  # select(var, val, `n patients` = n_patients, `% (95% CI)` = `%_(95% CI)`) |>
  select(var, val, `n patients` = n_patients, `%` = `prop100`) |>
  flextable() |>
  set_header_labels(var = "", 
                    val = "") |>
  merge_v(j = 1) |>
  valign(i = NULL, j = 1,    valign = "top",    part = "body") |>
  add_table_theme() |>
  save_as_docx(path = str_glue("{save_dir}Table1_all_{run_date}.docx"))

rm(summary_total, total_numb_pat)

# Repeat for suppressed values
# load data
summary_total_suppress <- read.csv(str_glue("{results_loc_main_suppress}summary_all_{results_generated_dt_main}.csv"))
summary_total_suppress <-
  summary_total_suppress |>
  select(var, val,
         n_patients = suppressed_n_tidy,
         N_patients = suppressed_N,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)
# Total number of patients
total_numb_pat <- 
  data.frame(
    var       = "Total",
    n_patients = unique(summary_total_suppress$N_patients)
  )

# format and tidy table to export into word
summary_total_suppress |>
  mutate(across(starts_with("prop100"), ~format_numb(.x))) |>
  format_ci() |>
  bind_rows(total_numb_pat) |>
  mutate(var = factor(var),
         var = factor(var, levels = c("Total", "Sex", "Age", "IMD"))) |>
  arrange(var) |>
  # select(var, val, `n patients` = n_patients, `% (95% CI)` = `%_(95% CI)`) |>
  select(var, val, `n patients` = n_patients, `%` = `prop100`) |>
  flextable() |>
  set_header_labels(var = "", 
                    val = "") |>
  merge_v(j = 1) |>
  valign(i = NULL, j = 1,    valign = "top",    part = "body") |>
  add_table_theme() |>
  save_as_docx(path = str_glue("{save_dir_suppress}Table1_all_{run_date}.docx"))

# End script ----
print(str_glue("***** {current_script}: files saved at {save_dir} ******"))
print(str_glue("***** {current_script}: detach packages ******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
