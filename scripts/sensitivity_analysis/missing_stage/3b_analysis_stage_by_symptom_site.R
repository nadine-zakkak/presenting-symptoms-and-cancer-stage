## --------------------------------
## Script name: analysis_stage_by_symptom_site.R
## Author: Nadine Zakkak
## Purpose: Sensitivity (incl. missing stage) analysis of proportion of stages of cancer stratified 
## by 6 cancer sites and relevant alarm and non-alarm symptoms
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (missing stage)-3b_analysis_stage_by_symptom_site.R"

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
sx_long <- readRDS(str_glue("{data_loc}data_sx_long_miss_{data_dt}.rds"))

# Use exemplar cancer sites and symptoms---- 
print(paste0("*****", current_script, ": Analysis using exemplar red flag and non-specific symptoms and cancer sites******"))
# red-flag + 1 non-specific symptom -----
cancer_sites <- names(cancer_redflag)
symptoms     <- c(nonspecific_symptoms, unname(unlist(cancer_redflag)))
stage_bin_by_symptom_site_exemplar <-
  sx_long |> 
  find_stage_prop_by_site_and_sx(symptoms = symptoms, cancer_sites = cancer_sites, stage_type = stage_bin)
stage_by_symptom_site_exemplar <-
  sx_long |> 
  find_stage_prop_by_site_and_sx(symptoms = symptoms, cancer_sites = cancer_sites, stage_type = stage_recode)
stage_bin_by_symptom_site_exemplar |> write.csv(str_glue("{results_loc}stage_bin_by_symptom_site_exemplar_{run_date}.csv"), row.names = F)
stage_by_symptom_site_exemplar     |> write.csv(str_glue("{results_loc}stage_by_symptom_site_exemplar_{run_date}.csv"), row.names = F)

# Suppress values
stage_bin_by_symptom_site_exemplar_suppress <- stage_bin_by_symptom_site_exemplar |> group_by(symptom, cancer_site_desc) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_by_symptom_site_exemplar_suppress     <- stage_by_symptom_site_exemplar     |> group_by(symptom, cancer_site_desc) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_bin_by_symptom_site_exemplar_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_symptom_site_exemplar_{run_date}.csv"), row.names = F)
stage_by_symptom_site_exemplar_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_symptom_site_exemplar_{run_date}.csv"),     row.names = F)

print(paste0("*****", current_script, ": detach packages******"))
# detach packages
sapply(pckgs, function(pckg)  { detach_package(pckg) })
rm(list = ls()); gc()

