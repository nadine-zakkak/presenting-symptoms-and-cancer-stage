## Script name: exclusions.R
##
## Author: Nadine Zakkak
## 
## Purpose: Preprocessing NCDA data for stage-symptom analysis: apply inclusion/exclusion criteria
##
## QA: Matt Barclay

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
source("./scripts/common_scripts/00_functions_analysis.R")

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse")
sapply(pckgs, function(pckg) { load_package(pckg) })
options(dplyr.summarise.inform = F)

# Set up working space and load data -----------------
current_script <- "1b_exclusions.R"
run_date <- readRDS("./run_date.rds")

print(paste0("*****", current_script, ": loading data******"))
data    <- readRDS(paste0(configs$preprocess$data_loc, "data_clean_", run_date, ".rds"))
data_sx <- readRDS(paste0(paste0(configs$preprocess$data_loc, "data_sx_clean_", run_date, ".rds")))

print(paste0("*****", current_script, ": Incl./Excl. steps******"))
# Apply exclusions and count -----------------------------------
cohort_counts <- data |> count_tumours_and_patients("original")

## Excl screen-detected cancers ----
data <- data |> filter(screen_detected != 1 & typereferral != 5)
cohort_counts <- cohort_counts |> bind_rows(data |> count_tumours_and_patients("Not via screening"))

## Excl. patients aged <25 ----
data <- data |> filter(tenyearageband != "<25") # exc 645 ROWS, 641 patients
cohort_counts <- cohort_counts |> bind_rows(data |> count_tumours_and_patients("Age >=25")) 

## Excl. diagnostic Interval < 0 or > 730 days ----
data <- data |> filter(is.na(diagnostic_interval) | between(diagnostic_interval, 0, 730)) 
cohort_counts <- cohort_counts |> bind_rows(data |> count_tumours_and_patients("diagnostic interval (0, 730) days or missing"))

## Excl. discordant sex - site -----
data <- data |> filter(!((sex == 1 & cancer_group == "Gynaecological") | (sex == 2 & cancer_group == "Prostate and other male organs"))) 
cohort_counts <- cohort_counts |> bind_rows(data |> count_tumours_and_patients("concordant sex-site"))

## Excl. discordant sex - symptom ----
# get patid of patients with discordant sex-symptom record
discord_patid <- 
  data_sx |> 
  filter(
    (sex == 1 & symptom_group == "Female specific") |
      (sex == 2 & symptom_group == "Male specific")) |>
  select(patid) |> 
  pull()

data <- data |> filter(!(patid %in% discord_patid))
cohort_counts <- cohort_counts |> bind_rows(data |> count_tumours_and_patients("concordant sex-symptom"))

## Excl. cancer sites of CNS, Multiple myeloma and of ICD-10 codes C91, C92, C95 relating to Leukaemia ----
data <-
  data |>
  filter(!(grepl("leukaemia", cancer_site_desc)               |
             cancer_site_desc %in% c("CNS", "Multiple myeloma") |
             icd10_3dig       %in% c("C91", "C92", "C95")))
cohort_counts <-  cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. Leukaemia, multiple myeloma or CNS related cancer"))

## Excl. tumours with recorded stage 0 of cancer (all breast cancer) ----
data <- data |> filter(stage_recode != 0) 
cohort_counts <- cohort_counts |> bind_rows(data |> count_tumours_and_patients("excl. stage 0 cancer"))

# Find the differences in patient and tumour counts across the different steps
cohort_counts <- 
  cohort_counts |>
  mutate(patients_diff = as.numeric(as.character(n_patients)) - lag(as.numeric(as.character(n_patients))),
         tumours_diff  = as.numeric(as.character(n_tumours))  - lag(as.numeric(as.character(n_tumours)))) |>
  select(description, n_patients, patients_diff, n_tumours, tumours_diff)

print(paste0("*****", current_script, ": saving data and updating json file******"))
# Save data -----
data          |> saveRDS(paste0(configs$preprocess$data_loc, "data_filtered_", run_date, ".rds"))
cohort_counts |> write.csv(paste0(configs$preprocess$data_loc, "cohort_flowchart_", run_date, ".csv"), row.names = F)

print(paste0("*****", current_script, ": deleting data and detaching packages******"))
# detach packages ----
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
