## --------------------------------
## Script name: preprocessing.R
## Author: Nadine Zakkak
## Purpose: Preprocessing NCDA data for main stage-symptom analysis
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
source("./scripts/common_scripts/00_global_vars.R")
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_formatting.R")
current_script <- "Main-1d_preprocessing.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse")
sapply(pckgs, function(pckg) { load_package(pckg) })

options(dplyr.summarise.inform = F)

print(paste0("*****", current_script, ": Loading data*****"))
# Load data into R ------
run_date <- readRDS("./run_date.rds")

cohort   <- readRDS(paste0(configs$preprocess$data_loc, "data_filtered_",  configs$preprocess$last_generated_dt, ".rds"))
data     <- readRDS(paste0(configs$preprocess$data_loc, "data_wide_",     configs$preprocess$last_generated_dt, ".rds"))
sx_long  <- readRDS(paste0(configs$preprocess$data_loc, "data_sx_long_",  configs$preprocess$last_generated_dt, ".rds"))
bt_long  <- readRDS(paste0(configs$preprocess$data_loc, "data_bt_long_",  configs$preprocess$last_generated_dt, ".rds"))
inv_long <- readRDS(paste0(configs$preprocess$data_loc, "data_inv_long_", configs$preprocess$last_generated_dt, ".rds"))
cohort_counts <- read.csv(paste0(configs$preprocess$data_loc, "cohort_flowchart_", configs$preprocess$last_generated_dt, ".csv"))

# Exclude missing stage ----
print(paste0("*****", current_script, ": Excluding missing stage*****"))
data <- data |> filter(stage_recode != -99)
cohort_counts <-  cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. missing stage"))

# Update 'symptom long' dataframes with filtered cohort
# helper function to update the dataframes
update_df <- function(df, data){
  df |> inner_join(data |> select(patid, tumourid), by = c("patid", "tumourid"))
}

sx_long  <- update_df(sx_long, data)

# Exclude tumours that only presented with symptoms that had < 100 observations in analysis sample, by gender ------
## Step 1: find out which symptoms are going to be excluded from our analysis
N_tumour_by_sex_symptom <- 
  sx_long |> 
  group_by(sex, symptom) |> 
  summarise(N_tumour = sum(present))

## Using detailed stage, find the proportion of patients with each stage for sex and symptom groups
stage_by_sex_symptom <- sx_long |> symptom_stage_count(stage_type = stage_recode)

print(paste0("*****", current_script, ": Finding which symptoms were excluded******"))
# Find which symptoms were excluded ----
# Remove discordant sex-symptom records from excluded symptom list
stage_by_sex_symptom <-
  stage_by_sex_symptom |>
  filter(!((sex_lb   == "Men"   & symptom %in% discord_sx_by_sex$Men) |
             (sex_lb == "Women" & symptom %in% discord_sx_by_sex$Women))
  )
# Which symptoms were excluded based on the 100 tumours threshold
symptoms_excluded_by_sex <- 
  stage_by_sex_symptom |> 
  distinct(sex_lb, symptom, N_tumour) |> 
  filter(N_tumour < 100) |>
  arrange(sex_lb, desc(N_tumour)) |>
  select(sex = sex_lb, symptom, N_tumour)

## Step 2: exclude tumours in patients that only presented with those symptoms
# now extend this back to patient level
only_sx <- sx_long |> select(patid, tumourid, sex, symptom, present, multiple_sx, stage_recode)
only_sx <- only_sx |> filter(present == 1) |> label_sex()

# add flag if sex-symptom discordant
tumour_exclude <-
  only_sx |>
  inner_join(symptoms_excluded_by_sex, 
             by = c("sex_lb" = "sex", "symptom"))

# now join back with flag to all data
only_sx <-
  only_sx |>
  left_join(tumour_exclude |> mutate(sx_exclude = T),
            by = c("patid", "tumourid", "sex", "symptom", "present", "multiple_sx", "stage_recode", "sex_lb")) |>
  mutate(sx_exclude = ifelse(is.na(sx_exclude), F, sx_exclude))

# now we will only exclude tumour if all sx_exclude = TRUE
only_sx <-
  only_sx |>
  group_by(patid, tumourid) |>
  mutate(exclude = all(sx_exclude == T)) |>
  ungroup()

# tumours to remain
tumours_to_keep <- only_sx |> filter(!exclude)
tumours_to_keep <- tumours_to_keep |> distinct(patid, tumourid)
# join back to sx_long so we get all symptoms whether present or not
data_updated <-
  data |>
  inner_join(tumours_to_keep, 
             by           = c("patid", "tumourid"), 
             relationship = "many-to-one")

cohort_counts <-  cohort_counts |> bind_rows(data_updated |> count_tumours_and_patients("excl. tumours with no eligible symptoms"))

# to be excluded
tumours_to_exclude <- only_sx |> filter(exclude) 
tumours_to_exclude |> summarise(n_tumours = n_distinct(tumourid),
                                n_patid   = n_distinct(patid)) 

tumours_to_exclude  |> saveRDS(paste0(configs$main_configs$data_loc, "excl_tumours_no_elig_sx_", run_date, ".rds"))

# count the difference in patient and tumour counts
cohort_counts <- 
  cohort_counts |>
  mutate(patients_diff = as.numeric(as.character(n_patients)) - lag(as.numeric(as.character(n_patients))),
         tumours_diff  = as.numeric(as.character(n_tumours))  - lag(as.numeric(as.character(n_tumours)))) |>
  select(description, n_patients, patients_diff, n_tumours, tumours_diff)

# update all dataframes with new cohort
sx_long  <- update_df(sx_long, data_updated)
bt_long  <- update_df(bt_long, data_updated)
inv_long <- update_df(inv_long, data_updated)
cohort   <- update_df(cohort, data_updated)

print(paste0("*****", current_script, ": Saving data and updating json file*****"))
# Save data for main analysis
cohort   |> saveRDS(paste0(configs$main_configs$data_loc, "data_",          run_date, ".rds"))
data     |> saveRDS(paste0(configs$main_configs$data_loc, "data_wide_",     run_date, ".rds"))
sx_long  |> saveRDS(paste0(configs$main_configs$data_loc, "data_sx_long_",  run_date, ".rds"))
bt_long  |> saveRDS(paste0(configs$main_configs$data_loc, "data_bt_long_",  run_date, ".rds"))
inv_long |> saveRDS(paste0(configs$main_configs$data_loc, "data_inv_long_", run_date, ".rds"))
cohort_counts |> write.csv(paste0(configs$main_configs$results_loc, "cohort_flowchart_", run_date, ".csv"), row.names = F)

# update json file
# Update main analysis data generated data
configs$main_configs$last_generated_dt   <- run_date
# Set main analysis preprocessing to False
configs$main_configs$regenerate          <- FALSE
# Set main analysis to TRUE 
configs$main_configs$analysis_regenerate <- TRUE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")

print(paste0("*****", current_script, ": deleting data and detaching packages******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()

