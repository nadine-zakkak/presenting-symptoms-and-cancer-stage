## --------------------------------
## Script name: preprocessing.R
## Author: Nadine Zakkak
## Purpose: 
## Preprocessing NCDA data for sensitivity (excl. prostate cancer) 
## stage-symptom analysis in men
## (this is based on the preprocessing in the main analysis)
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "SA (prostate cancer)-1d_preprocessing.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("dplyr", "tidyr")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)
run_date <- readRDS("./run_date.rds")


# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")

print(paste0("*****", current_script, ": Loading data*****"))
# Load data into R ------
data     <- readRDS(paste0(configs$preprocess$data_loc, "data_wide_",     configs$preprocess$last_generated_dt, ".rds"))
sx_long  <- readRDS(paste0(configs$preprocess$data_loc, "data_sx_long_",  configs$preprocess$last_generated_dt, ".rds"))
inv_long <- readRDS(paste0(configs$preprocess$data_loc, "data_inv_long_", configs$preprocess$last_generated_dt, ".rds"))
cohort_counts <- read.csv(paste0(configs$preprocess$data_loc, "cohort_flowchart_", configs$preprocess$last_generated_dt, ".csv"))

# Limit to only men ---
print(paste0("*****", current_script, ": Excluding women from analysis*****"))
data <- data |> filter(sex != 2)
cohort_counts <-  cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. women"))

# Exclude missing stage ----
print(paste0("*****", current_script, ": Excluding missing stage*****"))
data <- data |> filter(stage_recode != -99)
cohort_counts <-  cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. missing stage"))

# Exclude tumours with no eligible symptoms -----
print(paste0("*****", current_script, ": Excluding tumours with no eligible symptoms*****"))
tumours_to_exclude <- readRDS(paste0(configs$main_configs$data_loc, "excl_tumours_no_elig_sx_", configs$main_configs$last_generated_dt, ".rds"))
data <- data |> anti_join(tumours_to_exclude, by = c("tumourid", "patid"))
cohort_counts <-  cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. tumours with no eligible symptoms"))

# Exclude prostate cancer ----
print(paste0("*****", current_script, ": Excluding prostate cancer in men*****"))
data <- data |> filter(cancer_site_desc != 'Prostate')
cohort_counts <-  cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. prostate cancer"))

# helper function to update the dataframes
update_df <- function(df, data){
  df |> inner_join(data |> select(patid, tumourid), by = c("patid", "tumourid"))
}

# Update all dataframes with filtered data
sx_long  <- update_df(sx_long, data)
inv_long <- update_df(inv_long, data)

# count the difference in patient and tumour counts
cohort_counts <- 
  cohort_counts |>
  mutate(patients_diff = as.numeric(as.character(n_patients)) - lag(as.numeric(as.character(n_patients))),
         tumours_diff  = as.numeric(as.character(n_tumours))  - lag(as.numeric(as.character(n_tumours)))) |>
  select(description, n_patients, patients_diff, n_tumours, tumours_diff)

print(paste0("*****", current_script, ": Saving data and updating json file*****"))
# Save data for sensitifity analysis
results_loc <- paste0(configs$sa_configs$results_loc, str_split(sa_prst_loc, pattern = "/")[[1]][2], "/")
sx_long  |> saveRDS(paste0(configs$sa_configs$data_loc, "data_sx_long_",  run_date, ".rds"))
inv_long |> saveRDS(paste0(configs$sa_configs$data_loc, "data_inv_long_", run_date, ".rds"))
cohort_counts |> write.csv(paste0(results_loc, "cohort_flowchart_prostate_", run_date, ".csv"), row.names = F)

# update json file
# Update main analysis data generated data
configs$sa_configs$last_generated_dt <- run_date
# Set main analysis preprocessing to False
configs$sa_configs$regenerate <- FALSE
# Set main analysis to TRUE 
configs$sa_configs$analysis_regenerate <- TRUE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")

print(paste0("*****", current_script, ": deleting data and detaching packages******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()

