## --------------------------------
## Script name: analysis.R
## Author: Nadine Zakkak
## Purpose: 
## Analysis of NCDA data for sensitivity (excl. prostate cancer) 
## stage-symptom  in men
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "SA (prostate cancer)-2_analysis.R"

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
data_loc    <- configs$sa_configs$data_loc
data_dt     <- configs$sa_configs$last_generated_dt
results_loc <- paste0(configs$sa_configs$results_loc, str_split(sa_prst_loc, pattern = "/")[[1]][2], "/")
results_loc_suppress <- str_glue("{results_loc}/suppressed/")

# Load data into R ------
print(paste0("*****", current_script, ": Loading data*****"))
sx_long  <- readRDS(paste0(data_loc, "data_sx_long_",  data_dt, ".rds"))
inv_long <- readRDS(paste0(data_loc, "data_inv_long_", data_dt, ".rds"))
sx_fig_rename <- read.delim("./lookup_files/sx_figures_rename.txt")

# Transform data into long symptom-investigation format
sx_inv_long <-
  sx_long |>
  rename(present_sx = present) |>
  left_join(inv_long |>  select(patid, tumourid, inv, present_inv = present), 
            by           = c("patid", "tumourid"), 
            multiple     = "all",
            relationship = "many-to-many")

rm(sx_long, inv_long)

# Get the total number of tumours by symptom, sex and investigation types
N_tumour_by_sex_symptom_invt <- 
  sx_inv_long |> 
  group_by(sex, symptom, inv, present_inv) |> 
  summarise(N_tumour = sum(present_sx)) |> 
  ungroup()

## Using detailed stage, find the proportion of patients with each stage for sex and symptom and investigation groups
stage_by_sex_symptom_inv     <- sx_inv_long |> symptom_inv_stage_count(stage_type = stage_recode)
stage_bin_by_sex_symptom_inv <- sx_inv_long |> symptom_inv_stage_count(stage_type = stage_bin)

stage_by_sex_symptom_inv     |> write.csv(str_glue("{results_loc}stage_by_sex_symptom_inv_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_inv |> write.csv(str_glue("{results_loc}stage_bin_by_sex_symptom_inv_{run_date}.csv"), row.names = F)

## Suppress values
stage_by_sex_symptom_inv_suppress     <- stage_by_sex_symptom_inv     |> group_by(sex_lb, symptom, inv, present_inv) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour) 
stage_bin_by_sex_symptom_inv_suppress <- stage_bin_by_sex_symptom_inv |> group_by(sex_lb, symptom, inv, present_inv) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)

stage_by_sex_symptom_inv_suppress   |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_symptom_inv_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_inv_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_inv_{run_date}.csv"), row.names = F)

# Further stratify by blood test status -----
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

# What was the absolute percentage points difference in stage distriubtion -----
t_test <- read.csv(paste0(configs$main_configs$results_loc, "t_test_bt_", configs$main_configs$results_generated_dt, ".csv")) 

stage_bin_by_sex_symptom_bt <- 
  stage_bin_by_sex_symptom_bt |>
  select(-inv) |>
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")

# Only keep rows with "significant" p-values from main analysis
vis_df <- 
  stage_bin_by_sex_symptom_bt |>
  inner_join(t_test|> filter(p <= 0.05), by = c('sex_lb', 'symptom'))

single_stg <- 
  vis_df |>
  filter(stage_bin == '4') |>
  pivot_wider(id_cols     = c(sex, symptom, stage_bin),
              names_from  = blood_test,
              values_from = prop100) |>
  mutate(abs_diff = abs(`1` - `0`)) 

diff_range <-
  single_stg |>
  group_by(sex) |>
  summarise(min_diff = min(abs_diff),
            max_diff = max(abs_diff))
write.csv(diff_range, paste0(results_loc, "diff_perc_bt_", run_date, ".csv"), row.names = F)

# End script ----
print(paste0("*****", current_script, ": updating json file******"))
# update json file
# update date of generation of results
configs$sa_configs$results_generated_dt <- run_date
# Set sensitivity analysis to false
configs$sa_configs$analysis_regenerate <- FALSE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")

print(paste0("*****", current_script, ": detaching packages and deleting data******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
