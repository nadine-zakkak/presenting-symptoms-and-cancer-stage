## Script name: create_analysis_data.R
##
## Author: Nadine Zakkak
## 
## Purpose: Preprocessing NCDA data for stage-symptom analysis: create the analysis data
##
## QA: Matt Barclay

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
source("./scripts/common_scripts/00_functions_analysis.R")

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "fastDummies")
sapply(pckgs, function(pckg) { load_package(pckg) })
options(dplyr.summarise.inform = F)

# Set up working space and load data -----------------
current_script <- "1c_create_analysis_data.R"
run_date <- readRDS("./run_date.rds")

print(paste0("*****", current_script, ": loading data******"))
data    <- readRDS(paste0(configs$preprocess$data_loc, "data_filtered_", run_date, ".rds"))
data_sx <- readRDS(paste0(paste0(configs$preprocess$data_loc, "data_sx_clean_", run_date, ".rds")))

# Import lookups ----------------------------------------------------------
sx_groups  <- read.delim("./lookup_files/symptom_groups.txt", sep = " ")
inv_lookup <- read.delim("./lookup_files/investigations.txt", colClasses = c("numeric", "character", "character"))
bt_lookup  <- read.delim(file = "./lookup_files/bloodtst_desc.txt")

# Update data_sx to only include filtered cohort ----
data_sx <- data_sx |> inner_join(data |> distinct(patid, tumourid), by = c("patid", "tumourid"))

print(paste0("*****", current_script, ": Investigations and blood tests ******"))
# Investigations and blood tests -----
# each investigation on separate row
# Add in blood test descriptions
# Add in investigation description

# blood tests are in investigations (Numbers 1 to 10)
data_invest <-
  data |> 
  select(patid, tumourid, investigations) |>
  separate_rows(investigations, sep = " ") |>
  mutate(investigations = as.numeric(investigations)) |>
  left_join(bt_lookup, by = c("investigations" = "blood_test")) |>
  left_join(inv_lookup |> distinct(investigation_numb, investigation_type), by = c("investigations" = "investigation_numb"))

# Create long format data with investigations, symptoms and blood tests
data_long <-
  data_sx |>
  left_join(data_invest,
            by           = c('patid', 'tumourid'), 
            relationship = 'many-to-many')

# convert to wide format including all symptoms, investigation types and blood tests
data_wide <-
  data_long |>
  fastDummies::dummy_cols(
    select_columns          = c("symptom_recod_desc", "investigation_type", "blood_test_desc"),
    remove_selected_columns = T
  )

# keep only description of blood tests
data_wide <- data_wide |> select(-c(matches("^blood_test_desc_NA", perl = F)))

# For each tumour, flag if patient had had any of the symptoms, investigations and blood tests (detailed)
# Find the total number of symptoms per tumour
# Rename columns
# NOTE: NUMBER OF SYMPTOMS CAN BE 1 FOR PATIENTS WITH N/A AND N/K SYMPTOMS
data_wide <-
  data_wide |>
  group_by(patid, tumourid, sex) |>
  summarise(across(starts_with("symptom_recod_desc"), ~max(.x)),
            across(starts_with("investigation_type"), ~max(.x)),
            across(starts_with("blood_test_desc"),    ~max(.x))) |>
  ungroup()

data_wide <-
  data_wide |>
  mutate(numb_sx   = rowSums(across(starts_with("symptom_recod_desc")))) |>
  rename_with(.fn  = ~gsub("symptom_recod_desc_", "symptom_", .x),
              cols = starts_with("symptom_recod_desc_")) |>
  rename_with(.fn  = ~gsub("investigation_type_", "inv_", .x),
              cols = starts_with("investigation_type")) |>
  rename_with(.fn  = ~gsub("blood_test_desc_", "bt_", .x),
              cols = starts_with("blood_test_desc_"))

# Add in patient and cancer info ----
data_wide <-
  data_wide |>
  left_join(data |> select(patid,
                           tumourid,
                           pseudo_gp_code,
                           tenyearageband,
                           quintile_2019,
                           ca_reg_ethnicity_group,
                           site_icd10_o2,
                           cancer_site_desc,
                           cancer_group,
                           stage_recode,
                           stage_bin,
                           placepresentation_desc,
                           diagnostic_interval),
            by = c("patid", "tumourid")) |>
  as_tibble()

# Flag if patient/tumour had multiple symptoms
data_wide <- data_wide |> mutate(multiple_sx = numb_sx > 1)

print(paste0("*****", current_script, ": converting to long format dfs ******"))
## long format dfs ---
convert_to_long_format <- function(data, prefix) {
  update_prefix <- substr(prefix, 1, (nchar(prefix) - 1))
  data |>
    pivot_longer(starts_with(prefix),
                 names_prefix = prefix,
                 names_to     = update_prefix,
                 values_to    = "present")
}
# based on symptom
sx_long  <- data_wide |> convert_to_long_format("symptom_")
# based on blood tests
bt_long  <- data_wide |> convert_to_long_format("bt_")
# based on investigations
inv_long <- data_wide |> convert_to_long_format("inv_")

# Add symptom super groups ----
sx_long <- sx_long |> left_join(sx_groups, by = c("symptom"))

print(paste0("*****", current_script, ": saving data and updating json file******"))
# Save data -----
data_wide     |> saveRDS(paste0(configs$preprocess$data_loc, "data_wide_",        run_date, ".rds"))
sx_long       |> saveRDS(paste0(configs$preprocess$data_loc, "data_sx_long_",     run_date, ".rds"))
bt_long       |> saveRDS(paste0(configs$preprocess$data_loc, "data_bt_long_",     run_date, ".rds"))
inv_long      |> saveRDS(paste0(configs$preprocess$data_loc, "data_inv_long_",    run_date, ".rds"))

# update json file ----
# update last preprocess generated date
configs$preprocess$last_generated_dt     <- run_date
# Set preprocess regenerate to false
configs$preprocess$regenerate            <- FALSE
# Set all downstream processes to true
configs$main_configs$regenerate          <- TRUE
configs$main_configs$analysis_regenerate <- TRUE
configs$sa_configs$regenerate            <- TRUE
configs$sa_configs$analysis_regenerate   <- TRUE
configs_json <- toJSON(configs)
write(configs_json, "./scripts/common_scripts/00_configs.json")

print(paste0("*****", current_script, ": deleting data and detaching packages******"))
# detach packages ----
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()


