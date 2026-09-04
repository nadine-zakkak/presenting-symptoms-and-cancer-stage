## --------------------------------
## Script name: preprocessing.R
## Author: Nadine Zakkak
## Purpose: Preprocessing NCDA data for sensitivity stage (including missing)-symptom analysis
## exclude patients with tumours that didn't have any eligible symptoms
## (this is based on the preprocessing in the main analysis)
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (missing stage)-1d_preprocessing.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "xlsx")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)
run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")

# Locations where to read data from and write results to
data_loc    <- configs$preprocess$data_loc
data_dt     <- configs$preprocess$last_generated_dt

data_loc_main <- configs$main_configs$data_loc
data_dt_main  <- configs$main_configs$last_generated_dt

# Load data into R ------
print(str_glue("***** {current_script}: Loading data *****"))
cohort   <- readRDS(paste0(configs$preprocess$data_loc, "data_filtered_", configs$preprocess$last_generated_dt, ".rds"))
data     <- readRDS(paste0(configs$preprocess$data_loc, "data_wide_",     configs$preprocess$last_generated_dt, ".rds"))
sx_long  <- readRDS(paste0(configs$preprocess$data_loc, "data_sx_long_",  configs$preprocess$last_generated_dt, ".rds")) 
inv_long <- readRDS(paste0(configs$preprocess$data_loc, "data_inv_long_", configs$preprocess$last_generated_dt, ".rds")) 
bt_long  <- readRDS(paste0(configs$preprocess$data_loc, "data_bt_long_",  configs$preprocess$last_generated_dt, ".rds")) 

# load in excluded tumours based on absence of eligible symptoms
tumours_to_exclude <- 
  readRDS(paste0(data_loc_main, "excl_tumours_no_elig_sx_", data_dt_main, ".rds")) |>
  distinct(patid, tumourid)

cohort   <- cohort   |> anti_join(tumours_to_exclude, by = c("tumourid", "patid"))
data     <- data     |> anti_join(tumours_to_exclude, by = c("tumourid", "patid"))
sx_long  <- sx_long  |> anti_join(tumours_to_exclude, by = c("tumourid", "patid"))
inv_long <- inv_long |> anti_join(tumours_to_exclude, by = c("tumourid", "patid"))
bt_long  <- bt_long  |> anti_join(tumours_to_exclude, by = c("tumourid", "patid"))

# save data
cohort   |> saveRDS(paste0(configs$sa_configs$data_loc, "data_filtered_miss_",  run_date, ".rds"))
data     |> saveRDS(paste0(configs$sa_configs$data_loc, "data_wide_miss_",      run_date, ".rds"))
sx_long  |> saveRDS(paste0(configs$sa_configs$data_loc, "data_sx_long_miss_",   run_date, ".rds"))
inv_long |> saveRDS(paste0(configs$sa_configs$data_loc, "data_inv_long_miss_",  run_date, ".rds"))
bt_long  |> saveRDS(paste0(configs$sa_configs$data_loc, "data_bt_long_miss_",   run_date, ".rds"))

# Set main analysis to TRUE 
configs$sa_configs$analysis_regenerate <- TRUE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")

print(paste0("*****", current_script, ": deleting data and detaching packages******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
