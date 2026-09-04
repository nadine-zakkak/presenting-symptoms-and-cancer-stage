## --------------------------------
## Script name: analysis.R
## Author: Nadine Zakkak
## Purpose: Analysing NCDA data for sensitivity stage rtd-symptom analysis
## ---------------------------------

source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (RtD)-2_descriptive.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse")
sapply(pckgs, function(pckg) { load_package(pckg) })

run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")

# Locations where to read data from and write results to
data_sa_loc <- configs$sa_configs$data_loc
data_sa_dt  <- configs$sa_configs$last_generated_dt

results_loc <- paste0(configs$sa_configs$results_loc, str_split(sa_rtd_loc, pattern = "/")[[1]][2], "/")
results_loc_suppress <- str_glue("{results_loc}/suppressed/")

# Load data 
sx_long  <- readRDS(paste0(data_sa_loc, "data_sx_long_rtd_",  data_sa_dt, ".rds"))

## find the proportion of patients with each stage for sex and symptom and referral type
#' symptom_inv_stage_count
#' Find the proportion of stage stratified by
#' sex, symptom and investigation type
#'
#' @param data 
#' @param stage_type: specify which stage variable is required 
#'
#' @returns: summarised data frame with number and proportion of tumours 
#' for each of the stage level by sex, symptom and investigation type
symptom_rtd_stage_count <- function(data, stage_type) {
  data |>
    group_by(sex, symptom, typereferral_recode, {{stage_type}}) |>
    summarise(n_tumour = sum(present), .groups = "drop_last") |>
    mutate(N_tumour = sum(n_tumour)) |>
    mutate(prop100    = 100*n_tumour/N_tumour,
           prop100_lb = 100*calculate_ci(n_tumour, N_tumour)$lower,
           prop100_ub = 100*calculate_ci(n_tumour, N_tumour)$upper) |>
    label_sex() |>
    ungroup()
}

stage_bin_by_sex_symptom_rtd <- 
  sx_long |> 
  select(sex, symptom, present, typereferral_recode, stage_bin) |>
  symptom_rtd_stage_count(stage_type = stage_bin)
stage_bin_by_sex_symptom_rtd |> write.csv(str_glue("{results_loc}stage_bin_by_sex_symptom_rtd_{run_date}.csv"), row.names = F)

## Suppress values
stage_bin_by_sex_symptom_rtd_suppress <- stage_bin_by_sex_symptom_rtd |> group_by(sex_lb, symptom, typereferral_recode) |> suppress_values_with_prop(n_var = n_tumour) |> ungroup()
stage_bin_by_sex_symptom_rtd_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_rtd_{run_date}.csv"), row.names = F)

# build it back-up to get rtd by symptom without stage using the un(suppressed) information with stage stratification
summarise_total_tumours <- function(data, N) {
  data |> 
    distinct(sex, symptom, typereferral_recode, {{N}}) |>
    rename(n_tumour = {{N}}) |>
    group_by(sex, symptom) |> 
    mutate(N_tumour = sum(n_tumour),
           prop100    = 100*n_tumour/N_tumour
    ) |>
    ungroup()
}
N_tumour_by_sex_symptom_rtd <- stage_bin_by_sex_symptom_rtd |> summarise_total_tumours(N_tumour)
write.csv(N_tumour_by_sex_symptom_rtd, 
          str_glue("{results_loc}n_tumour_by_sex_symptom_rtd_{run_date}.csv"), 
          row.names = F)

N_tumour_by_sex_symptom_rtd_suppress <- stage_bin_by_sex_symptom_rtd_suppress |> summarise_total_tumours(suppressed_N)
write.csv(N_tumour_by_sex_symptom_rtd_suppress, 
          str_glue("{results_loc_suppress}n_tumour_by_sex_symptom_rtd_{run_date}.csv"), 
          row.names = F)

print(str_glue("***** {current_script}: deleting data and detaching packages ******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
          
          
