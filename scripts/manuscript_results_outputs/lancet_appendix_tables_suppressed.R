## Script name: lancet_appendix_tables_suppressed.R
##
## Author: Nadine Zakkak
## 
## Purpose: Produce tables for the appendix 
## with the exception of symptom combinations analysiss

# Initial setup ----
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "Appendix tables"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "officer", "flextable", "xlsx")
sapply(pckgs, function(pckg) { load_package(pckg) } )
options(dplyr.summarise.inform = F)

# helper functions and global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_for_tables.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")
run_date <- readRDS("./run_date.rds")

# Locations where to read data from
data_loc_main <- configs$main_configs$data_loc
data_dt_main  <- configs$main_configs$last_generated_dt

data_loc_sa_miss <- configs$preprocess$data_loc
data_dt_sa_miss  <- configs$preprocess$last_generated_dt

results_generated_dt_main <- configs$main_configs$results_generated_dt
results_loc_main          <- configs$main_configs$results_loc
results_loc_main_suppress <- str_glue("{results_loc_main}suppressed/")

results_generated_dt_sa_miss <- configs$sa_configs$results_generated_dt
results_loc_sa_miss          <- paste0(configs$sa_configs$results_loc, "missing_stage/")
results_loc_sa_miss_suppress <- str_glue("{results_loc_sa_miss}suppressed/")

results_generated_dt_sa_pc <- configs$sa_configs$results_generated_dt
results_loc_sa_pc          <- paste0(configs$sa_configs$results_loc, "prostate_cancer/")
results_loc_sa_pc_suppress <- str_glue("{results_loc_sa_pc}suppressed/")

save_dir    <- "./tables_figs/manuscript/suppressed"
if(!dir.exists(save_dir)) dir.create(save_dir)

# Table 1: All symptoms and flag for included symptoms ------------------
table_numb <- 1

sx_long <- readRDS(paste0(data_loc_main, "data_sx_long_", data_dt_main, ".rds"))
symptoms_excluded_by_gender <- read.csv(paste0(results_loc_main, "symptoms_excluded_by_sex_", results_generated_dt_main, ".csv"))

all_sx_by_gender <-
  sx_long |>
  distinct(symptom_group, Symptom = symptom) |>
  mutate(Men    = symptom_group != "Female specific",
         Women  = symptom_group != "Male specific") |>
  select(-symptom_group) |>
  arrange(desc(Men), desc(Women), Symptom) |>
  mutate(`Available in data for` = case_when(
    Men == T & Women == T ~ "M, W",
    Men == T              ~ "M",
    Women == T            ~ "W"
  )) |>
  select(-Men, -Women)

symptoms_included_by_gender <-
  symptoms_excluded_by_gender |>
  select(Gender = sex, Symptom = symptom) |>
  mutate(include = F) |>
  spread(Gender, include, fill = T) |>
  full_join(all_sx_by_gender, by = c("Symptom")) |>
  mutate(
    Men = case_when(
      is.na(Men) & grepl("M.*", `Available in data for`) ~ T,
      is.na(Men) ~ F,
      TRUE ~ Men
    ),
    Women = case_when(
      is.na(Women) & grepl("W.*", `Available in data for`, ignore.case = T) ~ T,
      is.na(Women) ~ F,
      TRUE ~ Women
    )
  ) |>
  mutate(`Included in analysis for` = case_when(
    Men == T & Women == T ~ "M, W",
    Men == T              ~ "M",
    Women == T            ~ "W",
    TRUE                  ~ ".."
  )) |>
  select(-c(Men, Women, `Available in data for`))

app_tb1 <- 
  all_sx_by_gender |>
  left_join(
    symptoms_included_by_gender,
    by = "Symptom"
  )

app_tb1 <-
  app_tb1 |>
  flextable() |>
  add_table_theme() |>
  set_caption(str_glue("Appendix Table {table_numb}. All available presenting symptoms in the original dataset."))

rm(sx_long, symptoms_excluded_by_gender)
table_numb <- table_numb + 1

# Table 2: List all cancer sites included in analysis ----
data <- readRDS(paste0(data_loc_main, "data_", data_dt_main, ".rds"))
app_tb2 <-
  data |>
  distinct(cancer_group, `Cancer Site` = cancer_site_desc, site_icd10_o2) |>
  mutate(Men   = cancer_group != "Gynaecological",
         Women = cancer_group != "Prostate and other male organs") |>
  select(-cancer_group) |>
  arrange(`Cancer Site`, site_icd10_o2) |>
  group_by(`Cancer Site`, Men, Women) |>
  summarise(`ICD-10` = paste(site_icd10_o2, collapse = ", "),
            .groups = "drop") |>
  mutate(Gender = case_when(
    Men == T & Women == T ~ "M, W",
    Men == T              ~ "M",
    Women == T            ~ "W"
  )) |>
  arrange(desc(Men), desc(Women), `Cancer Site`) |>
  select(-Men, -Women) |>
  left_join(
    data |>
      rename(`Cancer Site` = cancer_site_desc) |>
      group_by(`Cancer Site`) |>
      summarise(`Tumour count` = n(), .groups = "drop") |>
      mutate(`Tumour count` = prettyNum(`Tumour count`, big.mark = ",")),
    by = c("Cancer Site")
  )
app_tb2 <-
  app_tb2 |>
  flextable() |>
  add_table_theme() |>
  set_caption(stringr::str_glue("Appendix Table {table_numb}. Cancer sites analysed in the study."))
rm(data)

table_numb <- table_numb + 1

# Table 3: Characteristics by tumour multiplicity ----
print(paste0("*****", current_script, ": Charactestics by tumour multiplicity*****"))
# Load data
summary_by_mult_tumour <- read.csv(paste0(results_loc_main_suppress, "summary_by_mult_tumour_", results_generated_dt_main, ".csv"))

# Get total number of patients by tumour multiplicity
total_numb_pat_by_mult <- 
  summary_by_mult_tumour |> 
  select(mult_tumour, n_patients = N_patients) |> 
  distinct() |> 
  mutate(var   = "Total",
         n_patients = as.character(n_patients))

summary_by_mult_tumour <-
  summary_by_mult_tumour |>
  select(mult_tumour, var, val,
         n_patients = suppressed_n_tidy,
         N_patients = suppressed_N,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)

app_tb3 <-
  summary_by_mult_tumour |>
  mutate(across(starts_with("prop100"), ~format_numb(.x))) |>
  format_ci() |>
  bind_rows(total_numb_pat_by_mult) |>
  mutate(n_patients = prettyNum(n_patients, big.mark = ",")) |>
  pivot_wider(id_cols     = c(var, val),
              names_from  = mult_tumour,
              values_from = c(n_patients, `%_(95% CI)`),
              names_vary  = "slowest") |>
  mutate(var = if_else(var == "Sex", "Gender", var)) |> 
  mutate(var = factor(var),
         var = factor(var, levels = c("Total", "Gender", "Age", "IMD"))) |>
  arrange(var) |>
  # mutate(across(starts_with("n_patients"), ~format(.x |> as.numeric(), big.mark = ","))) |>
  flextable() |>
  add_header_row(values = c("", "Single Tumour", "Multiple Tumours"), colwidths = c(2,2,2)) |>
  set_header_labels(
    var                = "",
    val                = "",
    n_patients_FALSE        = "n patients",
    n_patients_TRUE         = "n patients",
    `%_(95% CI)_FALSE` = "column % (95% CI)",
    `%_(95% CI)_TRUE`  = "column % (95% CI)") |>
  merge_v(j = 1) |>
  # align( i = 1,    j = NULL, align = "center", part = "header") |>
  valign(i = NULL, j = 1,    valign = "top",   part = "body") |>
  add_table_theme() |>
  set_caption(stringr::str_glue("Appendix Table {table_numb}. Patient sample characteristics by tumour multiplicity"))
rm(summary_by_mult_tumour, total_numb_pat_by_mult)

table_numb <- table_numb + 1

# Table 4: Stage distribution by symptom, main analysis -----
sx_order                 <- readRDS(paste0(results_loc_main, "order_stage_by_sex_symptom_", results_generated_dt_main, ".rds"))
stage_by_sex_symptom     <- read.csv(paste0(results_loc_main_suppress, "stage_by_sex_symptom_", results_generated_dt_main, ".csv"))
stage_by_sex_symptom <- stage_by_sex_symptom |> filter(N_tumour >= 100)

stage_by_sex_symptom <-
  stage_by_sex_symptom |>
  select(sex, symptom, stage_recode,
         N_tumour   = suppressed_N,
         n_tumour   = suppressed_n_tidy,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)

app_tb4 <- stage_by_sex_symptom |> as.data.frame() |> prepare_table_stage_by_symptom(sx_order)
app_tb4 <-
  app_tb4 |>
  mutate(
    prop_ci_Stage1 = str_glue("{`Proportion (%) Stage 1`} {`(95% CI) Stage 1`}"),
    prop_ci_Stage2 = str_glue("{`Proportion (%) Stage 2`} {`(95% CI) Stage 2`}"),
    prop_ci_Stage3 = str_glue("{`Proportion (%) Stage 3`} {`(95% CI) Stage 3`}"),
    prop_ci_Stage4 = str_glue("{`Proportion (%) Stage 4`} {`(95% CI) Stage 4`}")
  ) |>
  select(-c(contains("Proportion"), contains("(95% CI)"))) |>
  select(Gender, Symptom, `Total tumour count`,
         contains("Stage 1"), prop_ci_Stage1,
         contains("Stage 2"), prop_ci_Stage2,
         contains("Stage 3"), prop_ci_Stage3,
         contains("Stage 4"), prop_ci_Stage4)
app_tb4 <-
  app_tb4 |>
  mutate(across(starts_with("Tumour count"), ~format(.x |> as.numeric(), big.mark = ","))) |>
  flextable() |>
  add_header_row(values = c("",
                            "Stage 1",
                            "Stage 2",
                            "Stage 3",
                            "Stage 4"),
                 colwidths = c(3,
                               2,
                               2,
                               2,
                               2)) |>
  set_header_labels(
    `Tumour count Stage 1` = "Tumour count",
    prop_ci_Stage1         = "% (95% CI)",
    `Tumour count Stage 2` = "Tumour count",
    prop_ci_Stage2         = "% (95% CI)",
    `Tumour count Stage 3` = "Tumour count",
    prop_ci_Stage3         = "% (95% CI)",
    `Tumour count Stage 4` = "Tumour count",
    prop_ci_Stage4         = "% (95% CI)") |>
  add_table_theme() |>
  set_caption(stringr::str_glue("Appendix Table {table_numb}. Cancer stage by symptom, stratified by gender."))
rm(sx_order, stage_by_sex_symptom)

table_numb <- table_numb + 1

# Table 5: Stage distribution by symptom and cancer site, main analysis-----
stage_by_symptom_site_exemplar <- read.csv(str_glue("{results_loc_main_suppress}stage_bin_by_symptom_site_exemplar_{results_generated_dt_main}.csv"))

stage_by_symptom_site_exemplar <-
  stage_by_symptom_site_exemplar |>
  select(symptom, cancer_site_desc, stage_bin,
         N_tumour   = suppressed_N,
         n_tumour   = suppressed_n_tidy,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)
app_tb5 <- stage_by_symptom_site_exemplar |> prepare_table_stage_by_symptom_exemplar(binary = T)
app_tb5 <-
  app_tb5 |>
  mutate(
    `prop_ci_Stage1-3` = str_glue("{`Proportion (%) Stage 1-3`} {`(95% CI) Stage 1-3`}"),
    prop_ci_Stage4 = str_glue("{`Proportion (%) Stage 4`} {`(95% CI) Stage 4`}")
  ) |>
  select(-c(contains("Proportion"), contains("(95% CI)"))) |>
  select(Symptom, `Cancer Site`, `Total tumour count`,
         contains("Stage 1"), `prop_ci_Stage1-3`,
         contains("Stage 4"), prop_ci_Stage4)
app_tb5 <-
  app_tb5 |>
  flextable() |>
  add_header_row(values = c("",
                            "Stage 1-3",
                            "Stage 4"),
                 colwidths = c(3,
                               2,
                               2)) |>
  set_header_labels(
    `Tumour count Stage 1-3` = "Tumour count",
    `prop_ci_Stage1-3`       = "% (95% CI)",
    `Tumour count Stage 4` = "Tumour count",
    prop_ci_Stage4         = "% (95% CI)") |>
  add_table_theme() |>
  set_caption(stringr::str_glue("Appendix Table {table_numb}. Cancer stage by symptom, stratified by cancer site"))
rm(stage_by_symptom_site_exemplar)

table_numb <- table_numb + 1

# Table 6: Stage distribution by blood test, main analysis -----
stage_bin_by_sex_symptom_bt  <- read.csv(paste0(results_loc_main_suppress, "stage_bin_by_sex_symptom_bt_", results_generated_dt_main, ".csv"))
t_test                       <- read.csv(paste0(results_loc_main, "t_test_bt_",  results_generated_dt_main, ".csv"))
N_tumour_sex_symptom <-
  stage_bin_by_sex_symptom_bt |>
  distinct(sex_lb, symptom, N_tumour) |>
  group_by(sex_lb, symptom) |>
  summarise(Total = sum(N_tumour))
stage_bin_by_sex_symptom_bt <-
  stage_bin_by_sex_symptom_bt |>
  inner_join(t_test |> filter(signif), by = c('sex_lb', 'symptom'))

stage_bin_by_sex_symptom_bt <-
  stage_bin_by_sex_symptom_bt |>
  select(sex_lb, symptom, blood_test, stage_bin,
         N_tumour   = suppressed_N,
         n_tumour   = suppressed_n_tidy,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)
app_tb6 <-
  stage_bin_by_sex_symptom_bt |> prepare_table_bt(N_tumour_sex_symptom, "Women", binary = T)  |>
  bind_rows(
    stage_bin_by_sex_symptom_bt |> prepare_table_bt(N_tumour_sex_symptom, "Men", binary = T)
  )
app_tb6 <-
  app_tb6 |>
  mutate(
    prop_ci_Stage1_3 = str_glue("{`Proportion (%) Stage 1-3`} {`(95% CI) Stage 1-3`}"),
    prop_ci_Stage4 = str_glue("{`Proportion (%) Stage 4`} {`(95% CI) Stage 4`}")
  ) |>
  select(-c(contains("Proportion"), contains("(95% CI)"))) |>
  select(Gender, Symptom, `Blood test flag`, `Total tumour count`,
         contains("Stage 1"), prop_ci_Stage1_3,
         contains("Stage 4"), prop_ci_Stage4)
app_tb6 <-
  app_tb6 |>
  flextable() |>
  add_header_row(values = c("",
                            "Stage 1-3",
                            "Stage 4"),
                 colwidths = c(4,
                               2,
                               2)) |>
  set_header_labels(
    `Tumour count Stage 1-3` = "Tumour count",
    prop_ci_Stage1_3         = "% (95% CI)",
    `Tumour count Stage 4`   = "Tumour count",
    prop_ci_Stage4           = "% (95% CI)") |>
  add_table_theme() |>
  set_caption(stringr::str_glue("Appendix Table {table_numb}. Cancer stage by symptom, stratified by blood test."))
rm(stage_bin_by_sex_symptom_bt, t_test, N_tumour_sex_symptom)

table_numb <- table_numb + 1

# Table 7: Stage distribution by blood test, sensitivity analysis excluding prostate cancer-----
stage_bin_by_sex_symptom_bt  <- read.csv(paste0(results_loc_sa_pc_suppress, "stage_bin_by_sex_symptom_bt_", results_generated_dt_sa_pc, ".csv"))
t_test                       <- read.csv(paste0(results_loc_main, "t_test_bt_",  results_generated_dt_main, ".csv"))
N_tumour_sex_symptom <-
  stage_bin_by_sex_symptom_bt |>
  distinct(sex_lb, symptom, N_tumour) |>
  group_by(sex_lb, symptom) |>
  summarise(Total = sum(N_tumour))
stage_bin_by_sex_symptom_bt <-
  stage_bin_by_sex_symptom_bt |>
  inner_join(t_test |> filter(signif), by = c('sex_lb', 'symptom'))

stage_bin_by_sex_symptom_bt <-
  stage_bin_by_sex_symptom_bt |>
  select(sex_lb, symptom, blood_test, stage_bin,
         N_tumour   = suppressed_N,
         n_tumour   = suppressed_n_tidy,
         prop100    = prop100_suppressed,
         prop100_lb = prop100_lb_min,
         prop100_ub = prop100_ub_max)
app_tb7 <-
  stage_bin_by_sex_symptom_bt |> prepare_table_bt(N_tumour_sex_symptom, "Women", binary = T)  |>
  bind_rows(
    stage_bin_by_sex_symptom_bt |> prepare_table_bt(N_tumour_sex_symptom, "Men", binary = T)
  )
app_tb7 <-
  app_tb7 |>
  mutate(
    prop_ci_Stage1_3 = str_glue("{`Proportion (%) Stage 1-3`} {`(95% CI) Stage 1-3`}"),
    prop_ci_Stage4 = str_glue("{`Proportion (%) Stage 4`} {`(95% CI) Stage 4`}")
  ) |>
  select(-c(contains("Proportion"), contains("(95% CI)"))) |>
  select(Symptom, `Blood test flag`, `Total tumour count`,
         contains("Stage 1"), prop_ci_Stage1_3,
         contains("Stage 4"), prop_ci_Stage4)
app_tb7 <-
  app_tb7 |>
  flextable() |>
  add_header_row(values = c("",
                            "Stage 1-3",
                            "Stage 4"),
                 colwidths = c(3,
                               2,
                               2)) |>
  set_header_labels(
    `Tumour count Stage 1-3` = "Tumour count",
    prop_ci_Stage1_3         = "% (95% CI)",
    `Tumour count Stage 4`   = "Tumour count",
    prop_ci_Stage4           = "% (95% CI)") |>
  add_table_theme() |>
  set_caption(stringr::str_glue("Appendix Table {table_numb}. Cancer stage, after excluding prostate cancer, by symptom in men, stratified by blood test."))
rm(stage_bin_by_sex_symptom_bt, t_test, N_tumour_sex_symptom)

table_numb <- table_numb + 1

# Add to word doc ------
doc <- read_docx()

for(numb in 1:(table_numb-1)) {
  table_var <- paste0("app_tb", numb)
  cat(str_glue("Writing {table_var}...."))
  doc <- body_add_flextable(doc, eval(as.name(table_var)))
}

doc <- body_set_default_section(doc, prop_section(page_size = page_size(orient = "landscape")))

print(doc, target = str_glue("{save_dir}/lancet_appendix_tables_{run_date}.docx"))
