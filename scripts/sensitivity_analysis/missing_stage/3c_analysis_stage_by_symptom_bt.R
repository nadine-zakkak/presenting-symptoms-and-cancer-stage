## --------------------------------
## Script name: analysis_stage_by_symptom_bt.R
## Author: Nadine Zakkak
## Purpose: Sensitivity (incl. missing stage) analysis of stages of cancer stratified 
## by symptom and blood test
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (missing stage)-3c_analysis_stage_by_symptom_bt.R"

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
inv_long <- readRDS(paste0(data_loc, "data_inv_long_miss_", data_dt, ".rds"))

print(str_glue("***** {current_script}: Describing stage distribution by symptom and inv type s******"))
# Further stratify by investigation type ----
# Transform data into long symptom-investigation format
sx_inv_long <-
  sx_long |>
  rename(present_sx = present) |>
  left_join(inv_long |>  select(patid, tumourid, inv, present_inv = present), 
            by           = c("patid", "tumourid"), 
            multiple     = "all",
            relationship = "many-to-many")

# Get the total number of tumours by symptom, sex and investigation types
N_tumour_by_sex_symptom_invt <- 
  sx_inv_long |> 
  group_by(sex, symptom, inv, present_inv) |> 
  summarise(N_tumour = sum(present_sx)) |> 
  ungroup()

## find the proportion of patients with each stage for sex and symptom and investigation groups
stage_by_sex_symptom_inv     <- sx_inv_long |> symptom_inv_stage_count(stage_type = stage_recode)
stage_bin_by_sex_symptom_inv <- sx_inv_long |> symptom_inv_stage_count(stage_type = stage_bin)

stage_by_sex_symptom_inv     |> write.csv(str_glue("{results_loc}stage_by_sex_symptom_inv_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_inv |> write.csv(str_glue("{results_loc}stage_bin_by_sex_symptom_inv_{run_date}.csv"), row.names = F)

## Suppress values
stage_by_sex_symptom_inv_suppress     <- stage_by_sex_symptom_inv     |> group_by(sex_lb, symptom, inv, present_inv) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour) 
stage_bin_by_sex_symptom_inv_suppress <- stage_bin_by_sex_symptom_inv |> group_by(sex_lb, symptom, inv, present_inv) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)

stage_by_sex_symptom_inv_suppress   |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_symptom_inv_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_inv_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_inv_{run_date}.csv"), row.names = F)

gc()

print(str_glue("***** {current_script}: Describing stage distribution by symptom and blood test ******"))
# Further stratify by blood test status -----
# Get the total number of tumours by symptom, sex and blood test presence
N_tumour_by_sex_symptom_bt <-
  sx_long |> 
  group_by(sex, symptom, blood_test = `inv_Blood test`) |> 
  summarise(N_tumour = sum(present)) |> 
  ungroup()

## find the proportion of patients with each stage for sex and symptom and blood test presence
stage_sx_bt_report <- function(data, stage_type = "stage") {
  data |>
    filter(inv == "Blood test") |>
    rename(blood_test = present_inv)
}

stage_by_sex_symptom_bt     <- stage_by_sex_symptom_inv     |> stage_sx_bt_report()
stage_bin_by_sex_symptom_bt <- stage_bin_by_sex_symptom_inv |> stage_sx_bt_report()

stage_by_sex_symptom_bt     |> write.csv(str_glue("{results_loc}stage_by_sex_symptom_bt_{run_date}.csv"), row.names = F)
stage_bin_by_sex_symptom_bt |> write.csv(str_glue("{results_loc}stage_bin_by_sex_symptom_bt_{run_date}.csv"), row.names = F)

## Suppress values
stage_by_sex_symptom_bt_suppress     <- stage_by_sex_symptom_bt     |> group_by(sex_lb, symptom, inv, blood_test) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour) 
stage_bin_by_sex_symptom_bt_suppress <- stage_bin_by_sex_symptom_bt |> group_by(sex_lb, symptom, inv, blood_test) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)

stage_by_sex_symptom_bt_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_symptom_bt_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_bt_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_bt_{run_date}.csv"), row.names = F)

# End script ----
print(str_glue("***** {current_script}: end {as.character(Sys.time())} ******"))
print(str_glue("***** {current_script}: detaching packages and deleting data ******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
