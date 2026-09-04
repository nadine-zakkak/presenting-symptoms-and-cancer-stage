## --------------------------------
## Script name: analysis_stage_by_symptom_bt.R
## Author: Nadine Zakkak
## Purpose: Main analysis of stages of cancer stratified 
## by symptom and blood test
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Main-3c_analysis_stage_by_symptom_bt.R"

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
data_loc    <- configs$main_configs$data_loc
data_dt     <- configs$main_configs$last_generated_dt
results_loc <- configs$main_configs$results_loc
results_loc_suppress <- str_glue("{results_loc}suppressed/")

# Load data into R ------
print(paste0("*****", current_script, ": Loading data*****"))
sx_long  <- readRDS(paste0(data_loc, "data_sx_long_",  data_dt, ".rds"))
inv_long <- readRDS(paste0(data_loc, "data_inv_long_", data_dt, ".rds"))
sx_fig_rename <- read.delim("./lookup_files/sx_figures_rename.txt")

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
rm(inv_long)

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
stage_by_sex_symptom_inv_suppress     <- stage_by_sex_symptom_inv     |> group_by(sex_lb, symptom, inv, present_inv) |> suppress_values_with_prop(n_var = n_tumour) 
stage_bin_by_sex_symptom_inv_suppress <- stage_bin_by_sex_symptom_inv |> group_by(sex_lb, symptom, inv, present_inv) |> suppress_values_with_prop(n_var = n_tumour)

stage_by_sex_symptom_inv_suppress   |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_symptom_inv_{run_date}.csv"),     row.names = F)
stage_bin_by_sex_symptom_inv_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_inv_{run_date}.csv"), row.names = F)

gc()

print(str_glue("***** {current_script}: Describing stage distribution by symptom and blood test ******"))
# Further stratify by blood test status and perform chi-squared test if applicable -----
# Get the total number of tumours by symptom, sex and blood test presence
## find the proportion of patients with each stage for sex and symptom and blood test presence
# MBQA - factorisation
stage_sx_bt_report <- function(data) {
  data |>
    filter(inv == "Blood test") |> rename(blood_test = present_inv)
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


# Permutation t-test if applicable -----
print(str_glue("***** {current_script}: permutation test ({n_permut} permutations) for stage by blood test******"))
# Get the total number of tumours by symptom, sex and blood test presence
N_tumour_by_sex_symptom_bt <-
  sx_long |> 
  group_by(sex, symptom, blood_test = `inv_Blood test`) |> 
  summarise(N_tumour = sum(present)) |> 
  ungroup()

# Flag if valid for chi-squared/t-test test: if at least 100 tumours exist in the sex-symptom group
# For each sex, symptom and blood test presence group, flag if there were 100 tumours and more
# then, regardless of blood test presence, flag the sex and symptom group as eligible if both blood test presence groups had 100 or more tumours
# otherwise flag as ineligible
sex_sx_bt_elig <- 
  N_tumour_by_sex_symptom_bt |>
  mutate(t_test = N_tumour >= 100) |>
  group_by(sex, symptom) |>
  summarise(t_test = as.logical(min(t_test))) |>
  label_sex()

# Filter data to only include sex-symptoms groups eligible for t-test and only keep required columns
t_data <- 
  sx_long |>
  left_join(sex_sx_bt_elig, 
            by = c('sex', 'symptom')) |>
  filter(t_test, present == 1) |>
  label_sex() |>
  mutate(stage_advanced = case_when(stage_bin == "4" ~ 1, stage_bin == "1-3" ~ 0)) |>
  select(patid, tumourid, sex_lb, symptom, blood_test = `inv_Blood test`, stage_advanced)

set.seed(1433401)
# Run permutation test using t-test
t_all <- run_permutation_test(t_data, n_permut, method = "t-test")
# Extract "base" run
t0 <- t_all[[1]]
# Tidy up results in a dataframe and 
# check if t-statistic of permutation was more extreme than base case (t <= t0 or t>= t0)
t_all_df <-  
  do.call(rbind, t_all[-1]) |> 
  full_join(t0 |> select(-permutation, t0 = t), by = c('sex_lb', 'symptom')) |> 
  mutate(extreme = abs(t) >= abs(t0))
# Calculate permutation test p-value by finding the mean number of tests that were more extreme than base case
t_p_vals <- 
  t_all_df |> 
  group_by(sex_lb, symptom) |> 
  summarise(p = mean(extreme)) |>
  mutate(signif    = p <= 0.05) |>
  mutate(method    = paste0("Permutation test with ", n_permut, " permutations"),
         statistic = "t-statistic") |>
  arrange(sex_lb, p)
t_p_vals |> write.csv(str_glue("{results_loc}t_test_bt_{run_date}.csv"), row.names = F)

## How many sex-symptom stratum were eligible to use to study blood test associations with stage ----
n_elig_symptoms_bt <- 
  t_p_vals |>
  distinct(sex_lb, symptom) |>
  count(sex_lb)
n_elig_symptoms_bt |> write.csv(str_glue("{results_loc}n_elig_symptoms_bt_{run_date}.csv"), row.names = F)

## What was the absolute percentage points difference in stage distribution -----
stage_bin_by_sex_symptom_bt <- 
  stage_bin_by_sex_symptom_bt |>
  select(-inv) |>
  left_join(sx_fig_rename |> select(symptom = Original, symptom_short = New), by = "symptom")
# Only keep rows with "significant" p-values
vis_df <- 
  stage_bin_by_sex_symptom_bt |>
  inner_join(t_p_vals, by = c('sex_lb', 'symptom')) |>
  filter(p <= 0.05) 

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
diff_range |> write.csv(str_glue("{results_loc}diff_perc_bt_{run_date}.csv"), row.names = F)

# End script ----
print(str_glue("***** {current_script}: end {as.character(Sys.time())} ******"))
print(str_glue("***** {current_script}: detaching packages and deleting data ******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
