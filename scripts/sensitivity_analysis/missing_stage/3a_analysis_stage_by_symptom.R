## --------------------------------
## Script name: analysis_stage_by_symptom.R
## Author: Nadine Zakkak
## Purpose: Sensitivity (incl. missing stage) analysis of proportion of stages of cancer by symptom presentation
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (missing stage)-3a_analysis_stage_by_symptom.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)
run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")

# Locations where to read data from and write results to
data_loc    <- configs$sa_configs$data_loc
data_dt     <- configs$sa_configs$last_generated_dt
results_loc <- paste0(configs$sa_configs$results_loc, str_split(sa_miss_loc, pattern = "/")[[1]][2], "/")
results_loc_suppress <- str_glue("{results_loc}/suppressed/")

# Load data into R ------
print(paste0("*****", current_script, ": Loading data*****"))
sx_long  <- readRDS(paste0(data_loc, "data_sx_long_miss_",  data_dt, ".rds"))

print(str_glue("***** {current_script}: Describing stage distribution by symptom******"))
# Describe stage distribution by symptom ----
# Get the total number of tumours by symptom and sex
N_tumour_by_sex_symptom <- 
  sx_long |> 
  group_by(sex, symptom) |> 
  summarise(N_tumour = sum(present))

## Find the proportion of patients with each stage for sex and symptom groups
stage_by_sex_symptom     <- sx_long |> symptom_stage_count(stage_type = stage_recode)
stage_bin_by_sex_symptom <- sx_long |> symptom_stage_count(stage_type = stage_bin)

stage_by_sex_symptom     |> write.csv(str_glue("{results_loc}stage_by_sex_symptom_{run_date}.csv"), row.names = F)
stage_bin_by_sex_symptom |> write.csv(str_glue("{results_loc}stage_bin_by_sex_symptom_{run_date}.csv"), row.names = F)

# Suppress values
stage_by_sex_symptom_suppress     <- stage_by_sex_symptom     |> group_by(sex_lb, symptom) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_bin_by_sex_symptom_suppress <- stage_bin_by_sex_symptom |> group_by(sex_lb, symptom) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_by_sex_symptom_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_symptom_{run_date}.csv"), row.names = F)
stage_bin_by_sex_symptom_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_{run_date}.csv"), row.names = F)

print(paste0("*****", current_script, ": detaching packages and deleting data******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
