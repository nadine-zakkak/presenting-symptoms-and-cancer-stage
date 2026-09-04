## Script name: preprocessing.R
##
## Author: Nadine Zakkak
## 
## Purpose: Preprocessing NCDA data for stage-symptom analysis
##
## QA: Matt Barclay

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse")
sapply(pckgs, function(pckg) { load_package(pckg) })
options(dplyr.summarise.inform = F)

# Set up working space and data -----------------
current_script <- "1a_preprocessing.R"
run_date <- readRDS("./run_date.rds")

print(paste0("*****", current_script, ": loading data******"))
data <- read.csv("../../Data_DONOTEDIT/Original_files/ODR1920_196_A2_Yoryos_Lyratzopoulos_NCDA_2018_010620/ODR1920_196_A2_Yoryos_Lyratzopoulos_NCDA_2018_010620.csv", stringsAsFactors = FALSE)

# change all column names to lower case and remove prefix
colnames(data) <- tolower(colnames(data))
data <- data |> rename_with(.fn = ~gsub("ncda_", "", .x), cols = starts_with("ncda_"))

# keep a selected number of columns
data <- 
  data |> 
  select(
    # IDs
    tumour_pseudoid, 
    patient_pseudoid, 
    pseudo_gp_code,
    # Demographics
    sex, 
    fiveyearageband, 
    quintile_2019,
    ca_reg_ethnicity_group,
    # Tumour info
    site_icd10_o2, 
    morph_coded, 
    morph_icd10_o2, 
    behaviour_icd10_o2, 
    behaviour_coded,
    grade, 
    stage_best, 
    stage_best_system,
    # Presentation info
    presentsymptom, 
    presentsymptom_desc, 
    presentsignortest, 
    presentsignortest_desc,
    placepresentation, 
    placepresentation_desc,
    # Investigations and referrals
    investigations, 
    investigations_desc,
    typereferral, 
    typereferral_desc,
    # Intervals and delays
    diagnostic_interval, 
    patient_interval, 
    pci, 
    pci_alt, 
    ref_to_spec_cons,
    delaypresentation, 
    delayafterreferral,
    # Hx of comorbidity
    chrl_tot_27_03, 
    chrl_tot_78_06,
    # Screening info
    screening_participation, 
    screening_participationd,
    screen_detected, 
    screen_detected_desc,
    comorbidity, 
    comorbidity_desc,
    # Hx of consultations info
    consultations, 
    multiconsults, 
    multiconsultsd, 
    consultations_unknown
  ) |> 
  rename(
    tumourid     = tumour_pseudoid,
    patid        = patient_pseudoid,
    stage        = stage_best,
    symptom      = presentsymptom,
    symptom_desc = presentsymptom_desc
  )

# Import lookups ----------------------------------------------------------
sx_recode        <- read.delim("./lookup_files/sx_recode.txt",      header = TRUE, sep = " ")
sx_recode_desc   <- read.delim("./lookup_files/sx_recode_desc.txt", header = TRUE)
sx_groups        <- read.delim("./lookup_files/symptom_groups.txt", sep = " ")

# THIS CODELIST WAS BASED ON THE AVAILABLE DATA - DOESN'T CONTAIN ALL ICD-10 CODES **
cancer_groups <- read.delim("./lookup_files/cancer_site_groups.txt") 

# Data cleaning -----------------------------------------------------------
# manipulate data i.e. age groups, stage cleaning
print(str_glue("***** {current_script}: Converting to 10 year age bands and recoding stage******"))
## Convert to 10 year age bands ---- 
data <-
  data |> 
  mutate(
    tenyearageband = factor(case_when(
      fiveyearageband %in% c("0-4", "5-9", "10-14", "15-19", "20-24") ~ "<25",
      fiveyearageband %in% c("25-29")                                 ~ "25-29",
      fiveyearageband %in% c("30-34", "35-39")                        ~ "30-39",
      fiveyearageband %in% c("40-44", "45-49")                        ~ "40-49",
      fiveyearageband %in% c("50-54", "55-59")                        ~ "50-59",
      fiveyearageband %in% c("60-64", "65-69")                        ~ "60-69",
      fiveyearageband %in% c("70-74", "75-79")                        ~ "70-79",
      fiveyearageband %in% c("80-84", "85-89")                        ~ "80-89",
      fiveyearageband == "90+"                                        ~ "90+",
      TRUE                                                            ~ NA
    ))
  )

# Add cancer sites description and super groups  ---------- 
# change C649 to C64
data <- data |> mutate(site_icd10_o2 = ifelse(site_icd10_o2 == "C649", "C64", site_icd10_o2))
data <- data |> left_join(cancer_groups, by = c("site_icd10_o2" = "icd10_4dig"))

## Recode stage ----
data <- data |> mutate(stage_recode = ifelse(stage %in% c("?", "U" , "", "A", "B", "C"), -99, stage))
#only extract number
data <- data |> mutate(stage_recode = parse_number(stage_recode))

## Set unknown primary cancer sites to stage IV ----
data <- data |> mutate(stage_recode = ifelse(cancer_site_desc == "Unknown primary", 4, stage_recode))

## Dichotomise stage -----
data <- 
  data |>
  mutate(
    stage_bin = case_when(
      stage_recode %in% c(1, 2, 3) ~ "1-3",
      stage_recode == 4            ~ "4",
      stage_recode == -99          ~ "Missing",
      TRUE                         ~ "Unknown")
  )

# Recode symptom (following MK - double checked) ----
# each symptom on separate row
data_sx <- data |> separate_rows(symptom, sep = " ") 
# recode symptoms
data_sx <- data_sx |> left_join(sx_recode |> mutate(symptom = as.character(sx_old)) |> select(symptom, symptom_recod = sx_new), by = 'symptom')
# add description of symptom (per row)
data_sx <- data_sx |> left_join(sx_recode_desc |> rename(symptom_recod_desc = sx_desc), by = c("symptom_recod" = "sx")) 
# only keep required columns
data_sx <- 
  data_sx |> 
  select(
    patid, tumourid, sex,
    symptom_recod, symptom_recod_desc
  ) |>
  distinct()
# Add symptom group descriptions
data_sx <- data_sx |> left_join(sx_groups, by = c("symptom_recod_desc" = "symptom"))

# Save data ------
data    |> saveRDS(paste0(configs$preprocess$data_loc, "data_clean_",    run_date, ".rds"))
data_sx |> saveRDS(paste0(configs$preprocess$data_loc, "data_sx_clean_", run_date, ".rds"))

# Unload packages and clean up environment ----
print(str_glue("*****{current_script}: deleting data and detaching packages******"))
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
