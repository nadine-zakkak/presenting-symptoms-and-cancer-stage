## --------------------------------
## Script name: preprocessing.R
## Author: Nadine Zakkak
## Purpose: Preprocessing NCDA data for sensitivity stage rtd-symptom analysis
## ---------------------------------

source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (RtD)-1d_preprocessing.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse")
sapply(pckgs, function(pckg) { load_package(pckg) })

run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")

# Locations where to read data from and write results to
data_orig_loc    <- configs$preprocess$data_loc
data_orig_dt     <- configs$preprocess$last_generated_dt
data_loc         <- configs$main_configs$data_loc
data_dt          <- configs$main_configs$last_generated_dt

data_sa_loc <- configs$sa_configs$data_loc

# Load data into R ------
print(paste0("*****", current_script, ": Loading data*****"))
sx_long  <- readRDS(paste0(data_loc, "data_sx_long_",  data_dt, ".rds"))

# load in data from initial pre-processing to get RtD
data_orig  <- readRDS(paste0(data_orig_loc, "data_filtered_", data_orig_dt, ".rds"))
data_orig <- 
  data_orig |>
  select(patid, tumourid,
         typereferral, typereferral_desc)

# add in information on type of referral and any other additional sign or test
print(paste0("*****", current_script, ": Adding in relevant information*****"))
sx_long <- sx_long |> left_join(data_orig, by = c("patid", "tumourid"))

# referral recode to EP, Non-EP and Unknown
sx_long <- 
  sx_long |>
  mutate(
    typereferral_recode = case_when(
      typereferral %in% c(0, 1, 6, 2, 3, 7, 8)   ~ "Non-EP",
      typereferral %in% c(4)   ~ "EP",
      typereferral %in% c(999) ~ "Unknown"
    )
  )

sx_long  |> saveRDS(paste0(data_sa_loc, "data_sx_long_rtd_",   run_date, ".rds"))

# Set main analysis to TRUE 
configs$sa_configs$analysis_regenerate <- TRUE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")

print(paste0("*****", current_script, ": deleting data and detaching packages******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
