## --------------------------------
## Script name: analysis_stage_by_symptom.R
## Author: Nadine Zakkak
## Purpose: Main analysis of proportion of stages of cancer by symptom presentation
## QA: Matt Barclay
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Main-3a_analysis_stage_by_symptom.R"

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

print(str_glue("***** {current_script}: Describing stage distribution by symptom******"))
# Describe stage distribution by symptom ----
# Get the total number of tumours by symptom and sex
N_tumour_by_sex_symptom <- 
  sx_long |> 
  group_by(sex, symptom) |> 
  summarise(N_tumour = sum(present))

## Using detailed stage, find the proportion of patients with each stage for sex and symptom groups
stage_by_sex_symptom     <- sx_long |> symptom_stage_count(stage_type = stage_recode)
stage_bin_by_sex_symptom <- sx_long |> symptom_stage_count(stage_type = stage_bin)
stage_by_sex_symptom     |> write.csv(str_glue("{results_loc}stage_by_sex_symptom_{run_date}.csv"), row.names = F)
stage_bin_by_sex_symptom |> write.csv(str_glue("{results_loc}stage_bin_by_sex_symptom_{run_date}.csv"), row.names = F)

# Suppress values
stage_by_sex_symptom_suppress     <- stage_by_sex_symptom     |> group_by(sex_lb, symptom) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_bin_by_sex_symptom_suppress <- stage_bin_by_sex_symptom |> group_by(sex_lb, symptom) |> suppress_values_with_prop(n_var = n_tumour, N_var = N_tumour)
stage_by_sex_symptom_suppress     |> write.csv(str_glue("{results_loc_suppress}stage_by_sex_symptom_{run_date}.csv"), row.names = F)
stage_bin_by_sex_symptom_suppress |> write.csv(str_glue("{results_loc_suppress}stage_bin_by_sex_symptom_{run_date}.csv"), row.names = F)

# Get order of symptoms to be used in all corresponding figure
sx_order <- 
  stage_bin_by_sex_symptom |> 
  filter(stage_bin == "1-3" & N_tumour >= 100) |> 
  arrange(sex_lb, prop100) |>
  select(sex_lb, symptom)

# convert to a list
sx_order_list <- as.list(split(sx_order$symptom, sx_order$sex_lb))
# move N/A and N/K to the beginning of the ordered symptoms in both men and women
sx_order_list <- lapply(sx_order_list, function(x) c("N/A", "N/K", setdiff(x, c("N/A", "N/K"))))
saveRDS(sx_order_list, str_glue("{results_loc}order_stage_by_sex_symptom_{run_date}.rds"))

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
symptoms_excluded_by_sex |> write.csv(str_glue("{results_loc}symptoms_excluded_by_sex_{run_date}.csv"), row.names = F)

print(paste0("*****", current_script, ": detaching packages and deleting data******"))
# detach packages
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
