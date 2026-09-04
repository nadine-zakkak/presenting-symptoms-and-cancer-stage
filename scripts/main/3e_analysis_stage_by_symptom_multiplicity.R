## --------------------------------
## Script name: analysis_stage_by_symptom_multiplicity.R
## Author: Nadine Zakkak
## Purpose: Main analysis of stages of cancer by symptom multiplicity (symptom count)
## QA: Matt Barclay
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Main-3e_analysis_stage_by_symptom_multiplicity.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "broom")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)
run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")

# Locations where to read data from and write results to
data_loc    <- configs$main_configs$data_loc
data_dt     <- configs$main_configs$last_generated_dt
results_loc <- configs$main_configs$results_loc
results_loc_suppress <- str_glue("{results_loc}suppressed/")

# Load data into R ------
print(paste0("*****", current_script, ": Loading data*****"))
data  <- readRDS(paste0(data_loc, "data_wide_",  data_dt, ".rds"))

print(str_glue("***** {current_script}: Describing stage distribution by symptom multiplicity******"))
# Set the symptoms with N/A or N/K to be "-1" counts
data <- 
  data |>
  mutate(numb_sx_update = ifelse(
    if_any(c(`symptom_N/A`, `symptom_N/K`), function(x) x == 1), 
    -1, 
    numb_sx)
  )

# Describe stage distribution by symptom ----
total_tumour_count <- function(data, numb_sx_multiplicity) {
  data |> 
    group_by({{numb_sx_multiplicity}}, sex) |> 
    summarise(N_tumour = n()) |>
    ungroup()
}

## MBQA - factorisation
symptom_numb_stage_count <- function(data, numb_sx_multiplicity, stage_type) {
  data |>
    group_by(sex, {{numb_sx_multiplicity}}, {{stage_type}}) |>
    summarise(n_tumour = n()) |> 
    inner_join(N_tumour_by_sex_symptom, by = c("sex", numb_sx_multiplicity |> substitute() |> deparse())) |>
    mutate(prop100    = 100*n_tumour/N_tumour,
           prop100_lb = 100*calculate_ci(n_tumour, N_tumour)$lower,
           prop100_ub = 100*calculate_ci(n_tumour, N_tumour)$upper) |>
    label_sex() |>
    ungroup()
}

# Get the total number of tumours by number of symptoms and sex
N_tumour_by_sex_symptom <- data |> total_tumour_count(numb_sx_multiplicity = numb_sx_update)

## Find the proportion of patients with each stage for sex and symptom multiplicity
stage_by_sex_symptom_numb     <- data |> symptom_numb_stage_count(numb_sx_multiplicity = numb_sx_update, stage_type = stage_recode)
stage_bin_by_sex_symptom_numb <- data |> symptom_numb_stage_count(numb_sx_multiplicity = numb_sx_update, stage_type = stage_bin)

stage_by_sex_symptom_numb     |> write.csv(str_glue("{results_loc}stage_by_sex_numb_symptom_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_numb |> write.csv(str_glue("{results_loc}stage_bin_by_sex_numb_symptom_{run_date}.csv"), row.names = F)

## suppress values
stage_by_sex_symptom_numb_suppress <- 
  stage_by_sex_symptom_numb |> 
  group_by(sex_lb, numb_sx_update) |> 
  suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_bin_by_sex_symptom_numb_suppress <- 
  stage_bin_by_sex_symptom_numb |> 
  group_by(sex_lb, numb_sx_update) |> 
  suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)

stage_by_sex_symptom_numb_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_numb_symptom_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_numb_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_numb_symptom_{run_date}.csv"), row.names = F)

# categorised into 1 symptom vs 2+ vs 3+ symptoms -----
# create new flag
data <- 
  data |>
  mutate(
    numb_sx_group2 = ifelse(numb_sx_update >= 2, "2+", ifelse(numb_sx_update == -1, "N/A or N/K", as.character(numb_sx_update))),
    numb_sx_group3 = ifelse(numb_sx_update >= 3, "3+", ifelse(numb_sx_update == -1, "N/A or N/K", as.character(numb_sx_update)))
  )

# Get the total number of tumours by number of symptoms and sex
N_tumour_by_sex_symptom <- 
  N_tumour_by_sex_symptom |>
  bind_rows(data |> total_tumour_count(numb_sx_multiplicity = numb_sx_group2)) |>
  bind_rows(data |> total_tumour_count(numb_sx_multiplicity = numb_sx_group3))

## Find the proportion of patients with each stage for sex and symptom multiplicity (2+)
stage_by_sex_symptom_categ_two     <- data |> symptom_numb_stage_count(numb_sx_multiplicity = numb_sx_group2, stage_type = stage_recode)
stage_bin_by_sex_symptom_categ_two <- data |> symptom_numb_stage_count(numb_sx_multiplicity = numb_sx_group2, stage_type = stage_bin)

stage_by_sex_symptom_categ_two     <- stage_by_sex_symptom_categ_two     |> select(-c(numb_sx_update, numb_sx_group3))
stage_bin_by_sex_symptom_categ_two <- stage_bin_by_sex_symptom_categ_two |> select(-c(numb_sx_update, numb_sx_group3))

stage_by_sex_symptom_categ_two     |> write.csv(str_glue("{results_loc}stage_by_sex_numb_categ2_symptom_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_categ_two |> write.csv(str_glue("{results_loc}stage_bin_by_sex_numb_categ2_symptom_{run_date}.csv"), row.names = F)

## suppress values
stage_by_sex_symptom_categ_two_suppress     <- stage_by_sex_symptom_categ_two     |> group_by(sex_lb, numb_sx_group2) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_bin_by_sex_symptom_categ_two_suppress <- stage_bin_by_sex_symptom_categ_two |> group_by(sex_lb, numb_sx_group2) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)

stage_by_sex_symptom_categ_two_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_numb_categ2_symptom_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_categ_two_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_numb_categ2_symptom_{run_date}.csv"), row.names = F)

## Find the proportion of patients with each stage for sex and symptom multiplicity (3+)
stage_by_sex_symptom_categ_three     <- data |> symptom_numb_stage_count(numb_sx_multiplicity = numb_sx_group3, stage_type = stage_recode)
stage_bin_by_sex_symptom_categ_three <- data |> symptom_numb_stage_count(numb_sx_multiplicity = numb_sx_group3, stage_type = stage_bin)

stage_by_sex_symptom_categ_three     <- stage_by_sex_symptom_categ_three     |> select(-c(numb_sx_update, numb_sx_group2))
stage_bin_by_sex_symptom_categ_three <- stage_bin_by_sex_symptom_categ_three |> select(-c(numb_sx_update, numb_sx_group2))

stage_by_sex_symptom_categ_three     |> write.csv(str_glue("{results_loc}stage_by_sex_numb_categ3_symptom_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_categ_three |> write.csv(str_glue("{results_loc}stage_bin_by_sex_numb_categ3_symptom_{run_date}.csv"), row.names = F)

## suppress values
stage_by_sex_symptom_categ_three_suppress     <- stage_by_sex_symptom_categ_three     |> group_by(sex_lb, numb_sx_group3) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_bin_by_sex_symptom_categ_three_suppress <- stage_bin_by_sex_symptom_categ_three |> group_by(sex_lb, numb_sx_group3) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)

stage_by_sex_symptom_categ_three_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_numb_categ3_symptom_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_categ_three_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_numb_categ3_symptom_{run_date}.csv"), row.names = F)

# update json file -----
# update date of generation of results
configs$main_configs$results_generated_dt <- run_date
# Set main analysis to false
configs$main_configs$analysis_regenerate <- FALSE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")
# detach packages -----
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
