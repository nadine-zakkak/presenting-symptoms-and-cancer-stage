## --------------------------------
## Script name: descriptive.R
## Author: Nadine Zakkak
## Purpose: Sensitivity analysis (missing stage) descriptive analysis
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Sensitivity analysis (missing stage)-2_descriptvie.R"

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
data_loc    <- configs$sa_configs$data_loc
data_dt     <- configs$sa_configs$last_generated_dt
results_loc <- paste0(configs$sa_configs$results_loc, str_split(sa_miss_loc, pattern = "/")[[1]][2], "/")
results_loc_suppress <- str_glue("{results_loc}/suppressed/")

# Load data into R ------
print(str_glue("***** {current_script}: Loading data *****"))
cohort    <- readRDS(paste0(data_loc, "data_filtered_miss_", data_dt, ".rds"))
data      <- readRDS(paste0(data_loc, "data_wide_miss_",     data_dt, ".rds"))
sx_long   <- readRDS(paste0(data_loc, "data_sx_long_miss_",  data_dt, ".rds")) 

print(str_glue("***** {current_script}: finding number of tumours by cancer sites ******"))
# Number and % of tumours ----
numb_tumours_site <- data |> find_prop_tumours_by_site()
write.csv(numb_tumours_site, str_glue("{results_loc}numb_tumours_site_{run_date}.csv"), row.names = F)
# No need to suppress values as we only use cancer sites with N>=100 tumours

print(str_glue("***** {current_script}: running multiple tumour summaries ******"))
# How often do we have patients with more than 1 tumour? -----
# Find the number of patients that had 1,2,3,... tumours
n_pat_by_tumour_numb <- data |> find_n_patients_by_tumour_count()
write.csv(n_pat_by_tumour_numb, str_glue("{results_loc}n_pat_by_tumour_numb_{run_date}.csv"), row.names = F)
# No need to suppress values

# Characteristics of cohort by the number of tumours ----
# flag patients if they had more than 1 tumour
data_tumour_count <- data |> flag_multiple_tumours(cohort_data = cohort)

# Get total number of patients by multiple tumour flag
N_pat_by_mult_tumour <- 
  data_tumour_count |> 
  group_by(mult_tumour) |> 
  summarise(N_patients = n_distinct(patid))

summary_by_mult_tumour <- 
  data_tumour_count |>
  count_by_mult_tumour(variable = sex) |>
  bind_rows(count_by_mult_tumour(variable = imd)) |>
  bind_rows(count_by_mult_tumour(variable = age)) |>
  inner_join(N_pat_by_mult_tumour, by = "mult_tumour") |>
  mutate(prop100    = 100*n_patients/N_patients,
         prop100_lb = 100*calculate_ci(n_patients, N_patients)$lower,
         prop100_ub = 100*calculate_ci(n_patients, N_patients)$upper) |>
  select(mult_tumour, var, val, starts_with("n", ignore.case = T), starts_with("prop100"))
summary_by_mult_tumour |> write.csv(str_glue("{results_loc}summary_by_mult_tumour_{run_date}.csv"), row.names = F)

## Suppress values
summary_by_mult_tumour_suppress <- summary_by_mult_tumour |> group_by(mult_tumour, var) |> suppress_values_with_prop(n_var = n_patients, N_var = N_patients)
summary_by_mult_tumour_suppress |> write.csv(str_glue("{results_loc_suppress}summary_by_mult_tumour_{run_date}.csv"), row.names = F)

print(paste0("*****", current_script, ": characterising cohort******"))
# Characterise cohort -----
# Get total number of patients
N_pat_total <- 
  N_pat_by_mult_tumour |> 
  summarise(N_patients = sum(N_patients)) |> 
  pull(N_patients)

# Summarise all cohort by sex, imd and age
summary_all <- 
  data_tumour_count |>
  count_any(variable = sex) |>
  bind_rows(count_any(variable = imd)) |>
  bind_rows(count_any(variable = age)) |>
  mutate(N_patients = N_pat_total,
         prop100    = 100*n_patients/N_patients,
         prop100_lb = 100*calculate_ci(n_patients, N_patients)$lower,
         prop100_ub = 100*calculate_ci(n_patients, N_patients)$upper) |>
  select(var, val, starts_with("n", ignore.case = T), starts_with("prop100"))
summary_all |> write.csv(str_glue("{results_loc}summary_all_{run_date}.csv"), row.names = F)

## Suppress values
summary_all_suppress <- summary_all |> group_by(var) |> suppress_values_with_prop(n_var = n_patients, N_var = N_patients)
summary_all_suppress |> write.csv(str_glue("{results_loc_suppress}summary_all_{run_date}.csv"), row.names = F)

print(str_glue("***** {current_script}: finding frequency of number of presenting symptoms per patient ******"))
# Find frequency of number of presenting symptoms per patient ----
excel_file_name      <- str_glue("{results_loc}numb_presenting_symptoms_{run_date}.xlsx")
# delete excel file if it exists or will run into errors because of same sheet naming
if(file.exists(excel_file_name)) file.remove(excel_file_name)

presenting_symptoms <- data |> filter(!(`symptom_N/A` == 1 | `symptom_N/K` == 1))

numb_presenting_sx <- presenting_symptoms |> count_tumour_by_sx_numb(numb_sx)
numb_presenting_sx |> write_to_excel(sheet_name = "individual", excel_file_name)

numb_presenting_sx_2plus <-
  presenting_symptoms |>
  mutate(new_numb_sx = ifelse(numb_sx == 1, 1, "2+")) |>
  count_tumour_by_sx_numb(new_numb_sx)
numb_presenting_sx_2plus |> write_to_excel(sheet_name = "1, 2+", excel_file_name)

numb_presenting_sx_3plus <-
  presenting_symptoms |>
  mutate(new_numb_sx = ifelse(numb_sx %in% c(1, 2), numb_sx, "3+")) |>
  count_tumour_by_sx_numb(new_numb_sx)
numb_presenting_sx_3plus |> write_to_excel(sheet_name = "1, 2, 3+", excel_file_name)

## No need to suppress values

# Number of symptoms included and excluded from general analyses -----
N_tumour_by_sex_symptom <-
  sx_long |>
  group_by(sex, symptom) |>
  summarise(N_tumour = sum(present)) |>
  ungroup()

# Categorise the number of symptoms into 2 groups >= 100 occurrences and <100 occurrences
n_symptoms_by_sex <- 
  N_tumour_by_sex_symptom |>
  filter(N_tumour > 0) |>
  mutate(desc = ifelse(N_tumour >= 100, "Freq >= 100", "Freq > 0 & Freq < 100")) |>
  group_by(sex, desc) |>
  summarise(n = n()) |>
  bind_rows(N_tumour_by_sex_symptom |> 
              filter(N_tumour > 0) |>
              group_by(sex) |> 
              summarise(n = n()) |> 
              mutate(desc = "All symptoms")) |>
  label_sex() |>
  pivot_wider(id_cols = sex_lb, names_from = desc, values_from = n) |>
  select(sex = sex_lb, `All symptoms`, everything())
write.csv(n_symptoms_by_sex, str_glue("{results_loc}n_symptoms_by_sex_{run_date}.csv"), row.names = F)
## No need to suppress values

# Find number of eligible symptoms in common in men and women and unique to each ----
numb_elig_symptoms <- N_tumour_by_sex_symptom |> find_n_elig_sx(cut_off = 100)
write.csv(numb_elig_symptoms, str_glue("{results_loc}numb_elig_symptoms_{run_date}.csv"), row.names = F)
## No need to suppress values

print(str_glue("***** {current_script}: deleting data and detaching packages ******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
