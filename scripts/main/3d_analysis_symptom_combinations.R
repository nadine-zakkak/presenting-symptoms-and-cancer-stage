## --------------------------------
## Script name: analysis_symptom_combinations.R
## Author: Nadine Zakkak
## Purpose: Main analysis of stages of cancer by symptom combinations
## ---------------------------------

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Main-3d_analysis_symptom_combinations.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "broom", "stringr", "janitor", "marginaleffects")
sapply(pckgs, function(pckg)  { load_package(pckg) })

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

# Load data into R ------
print(str_glue("*****{current_script}: Loading data *****"))
cohort    <- readRDS(paste0(data_loc, "data_",          data_dt, ".rds"))
data      <- readRDS(paste0(data_loc, "data_wide_",     data_dt, ".rds"))
data_long <- readRDS(paste0(data_loc, "data_sx_long_",  data_dt, ".rds"))

print(str_glue("*****{current_script}: Preprocessing data******"))
# Preprocess ----
# === === === === ===  === === === === ===  === === === === ===  === === === === ===  === === === === === 
# Exclude unknown primary from analysis as it is all advanced stage
model_data <-
  data |>
  select(sex,
         tenyearageband,
         cancer_site_desc,
         stage_bin,
         multiple_sx,
         starts_with("symptom_")
  ) |>
  filter(cancer_site_desc != "Unknown primary")

all_symptoms <- colnames(model_data |> select(starts_with("symptom")))

model_data <-
  model_data |>
  mutate(advanced = stage_bin == "4",
         sex      = factor(sex, levels = c("1", "2")),
         age      = factor(case_when(
           tenyearageband %in% c("25-29", "30-39") ~ "<40",
           TRUE ~ tenyearageband),
           levels = c("60-69",
                      "<40",
                      "40-49",
                      "50-59",
                      "70-79",
                      "80-89",
                      "90+")
         ),
         across(starts_with("symptom"), ~factor(.x, levels = c("0", "1")))
  ) |>
  select(-stage_bin)

print(str_glue("***** {current_script}: analysis ******"))
# === === === === ===  === === === === ===  === === === === ===  === === === === ===  === === === === === 
# Include interaction terms between symptoms -----
print(paste0("*****", current_script, ": Count pairwise symptom combinations******"))
# === === === === ===  === === === === ===  === === === === ===  === === === === ===  === === === === === 
count_pairwise_symptoms <- function(data, symptom1, symptom2){
  # only count rows where both symptoms are present
  data |>
    filter(!!sym(symptom1) == 1 &
             !!sym(symptom2) == 1) |>
    group_by(sex, cancer_site_desc) |>
    summarise(n = n()) |>
    mutate(symptom1 = symptom1,
           symptom2 = symptom2) |>
    select(sex, cancer_site_desc, symptom1, symptom2, n) |>
    ungroup()
}
all_symptoms <- colnames(data |> select(starts_with("symptom")))

## assumes combination symptom1-symptom2 is different than symptom2-symptom1 
all_pairs <- 
  expand.grid(symptom1 = all_symptoms, symptom2 = all_symptoms) |> 
  filter(symptom1 != symptom2) |>
  arrange(symptom1)

## unique pair of combinations: symptom1-symptom2 is equal to symptom2-symptom1 
unique_pairs <- combn(all_symptoms, 2) |> t() |> as.data.frame() |> set_names(c("symptom1", "symptom2"))

# actual potential combinations ----
# 1) exclude all N/A and N/K combinations
# 2) exclude all discord sex-specific symptom combinations
discord_sx_by_sex_comb <- expand.grid(symptom1 = discord_sx_by_sex$Women, symptom2 = discord_sx_by_sex$Men) |> 
  bind_rows(expand.grid(symptom1 = discord_sx_by_sex$Men, symptom2 = discord_sx_by_sex$Women)) |>
  mutate(across(c(symptom1, symptom2), function (x) paste0("symptom_", x)))

unique_pairs_actual <-
  unique_pairs |>
  filter(!(symptom1 %in% c("symptom_N/A", "symptom_N/K")) & !(symptom2 %in% c("symptom_N/A", "symptom_N/K"))) |>
  anti_join(discord_sx_by_sex_comb)

all_pairs_actual <-
  all_pairs |>
  filter(!(symptom1 %in% c("symptom_N/A", "symptom_N/K")) & !(symptom2 %in% c("symptom_N/A", "symptom_N/K"))) |>
  anti_join(discord_sx_by_sex_comb)

viable_symptom_comb_sex_site <- data.frame()
unique_viable_symptom_comb_sex_site <- data.frame()
for(row in 1:nrow(all_pairs_actual)) {
  if(row %% (10^(nchar(nrow(all_pairs_actual))-2)) == 0) cat(str_glue("{row}/ "))
  
  # Case 1: symptom1-symptom2 is different than symptom2-symptom1
  symptom1 <- as.character(all_pairs_actual$symptom1[row])
  symptom2 <- as.character(all_pairs_actual$symptom2[row])
  viable_symptom_comb_sex_site <- 
    viable_symptom_comb_sex_site |> 
    bind_rows(data |> 
                filter(multiple_sx) |>
                count_pairwise_symptoms(symptom1, symptom2))
  
  # Case 2: symptom1-symptom2 is the same as symptom2-symptom1
  if(row <= nrow(unique_pairs_actual)) {
    symptom1 <- as.character(unique_pairs_actual$symptom1[row])
    symptom2 <- as.character(unique_pairs_actual$symptom2[row])
    unique_viable_symptom_comb_sex_site <- 
      unique_viable_symptom_comb_sex_site |> 
      bind_rows(data |> 
                  filter(multiple_sx) |>
                  count_pairwise_symptoms(symptom1, symptom2))
  }
}

cat("\n")

unique_viable_symptom_comb_all <-
  unique_viable_symptom_comb_sex_site |>
  filter(n > 0) |>
  group_by(symptom1, symptom2) |>
  summarise(N = sum(n)) |> 
  ungroup()

unique_viable_symptom_comb_sex_site <-
  unique_viable_symptom_comb_sex_site |>
  mutate(cutoff_hundr = n >= 100)

cutoff_hundr_df <-
  unique_viable_symptom_comb_sex_site |>
  group_by(cutoff_hundr) |>
  summarise(freq = n())

unique_viable_symptom_all <-
  data_long |>
  group_by(sex, symptom, cancer_site_desc) |>
  summarise(N_tumour = sum(present)) |> 
  mutate(cutoff_hundr = N_tumour >= 100) |>
  label_sex() |>
  ungroup()

viable_symptom_comb_sex_site        |> write.csv(str_glue("{results_loc}n_symptom_combinations_{run_date}.csv"),                       row.names = F)
unique_viable_symptom_comb_sex_site |> write.csv(str_glue("{results_loc}n_unique_symptom_combinations_{run_date}.csv"),                row.names = F)
cutoff_hundr_df                     |> write.csv(str_glue("{results_loc}n_unique_symptom_combinations_cutoff_hundred_{run_date}.csv"), row.names = F)

model_data_new <- model_data |> janitor::clean_names()
lkup <- data.frame(old_sx = model_data |> colnames(),
                   new_sx = model_data_new |> colnames()
)
model_data <- model_data_new ; rm(model_data_new)

unique_viable_symptom_comb_sex_site  <- 
  unique_viable_symptom_comb_sex_site |> 
  left_join(lkup, by = c("symptom1" = "old_sx")) |>
  rename(symptom1_old = symptom1,
         symptom1     = new_sx) |>
  left_join(lkup, by = c("symptom2" = "old_sx")) |>
  rename(symptom2_old = symptom2,
         symptom2     = new_sx)

# Fit model only stratified by sex and cancer site -----
symptom_comb_alt <-
  unique_viable_symptom_comb_sex_site |>
  mutate(symptom = str_glue("`{symptom1}`*`{symptom2}`")) |>
  bind_rows(
    unique_viable_symptom_all |>
      mutate(symptom = str_glue("symptom_{symptom}")) |>
      left_join(lkup, by = c("symptom" = "old_sx")) |>
      rename(symptom_old = symptom,
             symptom     = new_sx)
  ) |>
  filter(cutoff_hundr) |>
  group_by(sex, cancer_site_desc) |>
  summarise(symptoms = paste(symptom, collapse = " + ")) |>
  filter(grepl("\\*", symptoms)) |>
  ungroup() |> 
  mutate(sex = factor(sex, levels = c("1", "2")))

# Data to be used for modelling
data_comb_alt <-
  model_data |>
  inner_join(symptom_comb_alt, 
             by       = c("sex", "cancer_site_desc"), 
             multiple = "all")

model_df <-
  data_comb_alt |>
  group_by(sex, cancer_site_desc, symptoms) |>
  nest() |>
  mutate(glm_mod = map(data,  ~glm(as.formula(paste0("advanced ~ age + ", symptoms))
                                   , data = .x,
                                   family = binomial(link = "logit")))  ) |>
  ungroup() |>
  select(sex, cancer_site_desc, glm_mod)

# tidy the model results to to get the OR estimates for each of the terms along with their 95% CI and p-values
model_all_tidy <-
  model_df |>
  group_by(sex, cancer_site_desc) |>
  mutate(estimates = map(glm_mod, ~tidy(.x, conf.int = T, exponentiate = TRUE))) |>
  select(sex, cancer_site_desc, estimates) |>
  unnest(estimates) |>
  select(sex, cancer_site_desc, term = term,
         OR = estimate, lower = conf.low, upper = conf.high, p = p.value) |>
  ungroup()
model_all_tidy |> write.csv(str_glue("{results_loc}model_new_symptom_combinations_{run_date}.csv"), row.names = F)

# find the significant interaction terms
signif_interaction_terms <-
  model_df |>
  group_by(sex, cancer_site_desc) |>
  mutate(estimates = map(glm_mod, ~tidy(.x, conf.int = T, exponentiate = TRUE))) |>
  select(sex, cancer_site_desc, estimates) |>
  unnest(estimates) |>
  select(sex, cancer_site_desc, term = term,
         OR = estimate, lower = conf.low, upper = conf.high, p = p.value) |>
  ungroup() |>
  filter(p < 0.05) |>
  filter(grepl(":", term)) |>
  mutate(signif_interaction = T)

# find non significant interaction terms
# update on 6/Mar/26 -- to find OR of both symptoms = 1 when they don't have a significant interaction which is excluded in hypothesis
non_signif_interaction_terms <-
  model_df |>
  group_by(sex, cancer_site_desc) |>
  mutate(estimates = map(glm_mod, ~tidy(.x, conf.int = T, exponentiate = TRUE))) |>
  select(sex, cancer_site_desc, estimates) |>
  unnest(estimates) |>
  select(sex, cancer_site_desc, term = term,
         OR = estimate, lower = conf.low, upper = conf.high, p = p.value) |>
  ungroup() |>
  filter(p >= 0.05) |>
  filter(grepl(":", term)) |>
  mutate(signif_interaction = F)

# set-up the hypothesis definition to get the OR when symptom combinations = 1
sx_to_test <- 
  signif_interaction_terms |>
  group_by(sex, cancer_site_desc) |>
  mutate(symptom1 = str_split_i(term, ":", i = 1),
         symptom2 = str_split_i(term, ":", i = 2)) |>
  mutate(hypothesis = paste0(symptom1, " + ", symptom2, " + `", term, "` = 0")) |>
  ungroup()

sx_to_test_no_interaction <-
  non_signif_interaction_terms |>
  bind_rows(signif_interaction_terms) |>
  group_by(sex, cancer_site_desc) |>
  mutate(symptom1 = str_split_i(term, ":", i = 1),
         symptom2 = str_split_i(term, ":", i = 2)) |>
  mutate(hypothesis = paste0(symptom1, " + ", symptom2, " = 0")) |>
  ungroup()

all_sx_to_test <-
  sx_to_test |> 
  bind_rows(sx_to_test_no_interaction)

# symptoms to be included in final plot because had significant interaction terms
sx_to_include <- 
  sx_to_test |> 
  distinct(sex, cancer_site_desc, symptom1) |> 
  rename(symptom = symptom1) |>
  bind_rows(  sx_to_test |> 
                distinct(sex, cancer_site_desc, symptom2) |> 
                rename(symptom = symptom2) ) |>
  distinct(sex, cancer_site_desc, symptom) |>
  mutate(symptom = str_replace_all(symptom, c(symptom_ = "", `1` = "")))

# find odds ratio of required symptom combinations using hypotheses function from marginal effects
pairwise_estimates <-
  model_df |>
  # update on 6/Mar/26 -- to find OR of both symptoms = 1 for 
  # (1) have significant interaction which is 
  # a. included in hypothesis
  # b. excluded in hypothesis
  # (2) don't have significant interaction which is excluded in hypothesis
  inner_join(all_sx_to_test |> select(sex, cancer_site_desc, signif_interaction, hypothesis),
             by = c("sex", "cancer_site_desc"), relationship = "one-to-many") |>
  group_by(sex, cancer_site_desc, signif_interaction, hypothesis) |>
  mutate(hypothesis_test = map(glm_mod, ~marginaleffects::hypotheses(.x, hypothesis = hypothesis, vcov = "HC3")
  )) |>
  rename(hypothesis_set = hypothesis) |>
  unnest(hypothesis_test, names_repair = "universal")  |>
  ungroup()

pairwise_estimates <- 
  pairwise_estimates |>
  mutate(across(c("estimate", "conf.low", "conf.high"), exp))  |>
  arrange(sex, cancer_site_desc, desc(signif_interaction), hypothesis_set)

pairwise_estimates <-
  pairwise_estimates |>
  select(-c(std.error, statistic, s.value, hypothesis_set, glm_mod)) |>
  rename(OR    = estimate,
         lower = conf.low,
         upper = conf.high,
         p     = p.value) |>
  mutate(sx_label = 
           paste0(
             str_split_i(hypothesis, pattern = "\\+", i = 1), " and ",
             str_split_i(hypothesis, pattern = "\\+", i = 2)
           ) |> 
           str_replace_all(c(symptom_ = "", `1` = "", `= 0` = "", "_" = " ", "=0" = "")) |>
           str_to_sentence() |>
           trimws() 
  )
pairwise_estimates |> write.csv(str_glue("{results_loc}model_symptom_pairwise_estimates_{run_date}.csv"), row.names = F)

# get the odds ratios of the symptoms from the original model, excluding the interaction terms
all_sx_model_tidy <-
  model_df |>
  group_by(sex, cancer_site_desc) |>
  mutate(estimates = map(glm_mod, ~tidy(.x, conf.int = T, exponentiate = TRUE))) |>
  select(sex, cancer_site_desc, estimates) |>
  unnest(estimates) |>
  select(sex, cancer_site_desc, term = term,
         OR = estimate, lower = conf.low, upper = conf.high, p = p.value) |>
  ungroup() |>
  filter(!grepl(":", term) & grepl("symptom_", term))

# only leave symptoms that were included in at least one significant combination
model_tidy_to_plot <- 
  all_sx_model_tidy |> 
  mutate(sx_label = str_replace_all(term, c(symptom_ = "", `1` = ""))) |>
  inner_join(sx_to_include, by = c("sx_label" = "symptom", "sex" = "sex", "cancer_site_desc" = "cancer_site_desc"))

# prepare the data to plot by extracting the symptom labels and preparing the OR with the 95% CI
plot_data <-
  pairwise_estimates |>
  filter(signif_interaction) |>
  mutate(with_interaction = grepl(":", hypothesis)) |>
  mutate(sx_label = 
           paste0(
             str_split_i(hypothesis, pattern = "\\+", i = 1), " and ",
             str_split_i(hypothesis, pattern = "\\+", i = 2)
           ) |> 
           str_replace_all(c(symptom_ = "", `1` = "", `= 0` = "", `=0` = "")) |>
           trimws()
  ) |>
  ungroup() |>
  select(sex, cancer_site_desc, sx_label, OR, lower, upper, with_interaction) |>
  bind_rows(model_tidy_to_plot |> select(sex, cancer_site_desc, sx_label, OR, lower, upper))

# tidy up symptom labels
plot_data <-
  plot_data |>
  mutate(sx_label = str_replace_all(sx_label, "_", " ")
  ) |>
  mutate(sx_label = str_to_sentence(sx_label)) |>
  mutate(sx_label = str_replace_all(sx_label, "Breast lump mass", "Breast lump/mass"))
plot_data |> write.csv(str_glue("{results_loc}model_symptom_combintion_to_plot_{run_date}.csv"), row.names = F)

print(str_glue("***** {current_script}: detaching packages and deleting data ******"))
# detach packages ----
sapply(pckgs, function(pckg) { detach_package(pckg) })
rm(list = ls()); gc()
