## Script name: appendix_symptom_combinations.R
##
## Author: Nadine Zakkak
## 
## Purpose: Produce figures and tables of symptom combinations for the appendix

# Initial set-up ------
source("./scripts/common_scripts/00_base_setup.R")
current_script <- "appendix_symptom_combinations.R"

# Load packages, suppressing warnings so easier when batch run
pckgs <- c("tidyverse", "flextable", "scales", "patchwork")
sapply(pckgs, function(pckg)  { load_package(pckg) })

options(dplyr.summarise.inform = F)
run_date <- readRDS("./run_date.rds")

# Load helper functions & global variables
source("./scripts/common_scripts/00_functions_analysis.R")
source("./scripts/common_scripts/00_functions_for_tables.R")
source("./scripts/common_scripts/00_functions_formatting.R")
source("./scripts/common_scripts/00_global_vars.R")

# Locations where to read data from 
results_generated_dt_main <- configs$main_configs$results_generated_dt
results_loc_main          <- configs$main_configs$results_loc
results_loc_main_suppress <- str_glue("{results_loc_main}suppressed/")

save_dir    <- str_glue("./tables_figs/manuscript/")
if(!dir.exists(save_dir)) dir.create(save_dir)

table_numb <- 1

## Pairwise symptom combinations -----
unique_viable_symptom_comb_sex_site <- read.csv(str_glue("{results_loc_main}n_unique_symptom_combinations_{results_generated_dt_main}.csv"))
unique_viable_symptom_comb_sex_site_tab <- 
  unique_viable_symptom_comb_sex_site |>
  filter(cutoff_hundr) |>
  label_sex() |>
  mutate(symptom1 = str_replace(symptom1, "symptom_", ""),
         symptom2 = str_replace(symptom2, "symptom_", "")) |>
  mutate(across(c(symptom1, symptom2), .fns = function(x) str_replace(x, "^LUTS", "Lower Urinary Tract Symptoms"))) |>
  mutate(`Symptom combination` = str_glue("{symptom1} - {symptom2}")) |>
  select(Gender = sex_lb, `Cancer Site` = cancer_site_desc, `Symptom combination`) |>
  arrange(desc(Gender), `Cancer Site`, `Symptom combination`)
# add in gender-only rows
women_end <- unique_viable_symptom_comb_sex_site_tab |> mutate(row_numb = row_number()) |> filter(Gender == "Women") |> tail(n = 1) |> pull(row_numb)
unique_viable_symptom_comb_sex_site_tab <-
  unique_viable_symptom_comb_sex_site_tab |>
  add_row(Gender = "Women", .before = 1) |>
  add_row(Gender = "Men", .after = women_end + 1)

# tidy and save as word document
unique_viable_symptom_comb_sex_site_tab |>
  flextable() |>
  merge_h(i = c(1, women_end)) |>
  merge_v(j = c(1, 2)) |>
  valign(j = c(1, 2), valign = "top") |>
  add_table_theme() |>
  save_as_docx(path = paste0(save_dir, "App1_Table", table_numb, "_symptom_combinations_", run_date, ".docx"))

table_numb <- table_numb + 1

rm(unique_viable_symptom_comb_sex_site, unique_viable_symptom_comb_sex_site_tab)

# pairwise symptom combination model results ------
or_pairwise <-  read.csv(str_glue("{results_loc_main}model_symptom_pairwise_estimates_{results_generated_dt_main}.csv"))
model_tidy <- read.csv(str_glue("{results_loc_main}model_new_symptom_combinations_{results_generated_dt_main}.csv"))

sx_terms <-
  model_tidy |>
  label_sex() |>
  select(-sex) |>
  select(Gender = sex_lb, everything()) |>
  arrange(Gender) |>
  filter(grepl("symptom", term))
interaction_terms <- sx_terms |> filter(grepl("\\:", term))

# extract OR for people with pairwise interactions and reverse back tidy symptoms
or_interaction <-
  or_pairwise |>
  filter(grepl("and", sx_label) & grepl(":", hypothesis)) |>
  mutate(sx1 = str_split_i(sx_label, " and ", i = 1),
         sx2 = str_split_i(sx_label, " and ", i = 2)
  ) |>
  mutate(across(c(sx1, sx2), function(x) {
    x <- x |> trimws() |> tolower() |> str_replace_all(" ", "_") |> str_replace_all("/", "_")
    paste0("symptom_", x, "1")
  })) |>
  mutate(sx_interaction = str_glue("{sx1}:{sx2}")) |>
  label_sex() |>
  select(-sex) |>
  select(Gender = sex_lb, everything())

or_combination <-
  or_pairwise |>
  filter(grepl("and", sx_label) & !grepl(":", hypothesis)) |>
  mutate(sx1 = str_split_i(sx_label, " and ", i = 1),
         sx2 = str_split_i(sx_label, " and ", i = 2)
  ) |>
  mutate(across(c(sx1, sx2), function(x) {
    x <- x |> trimws() |> tolower() |> str_replace_all(" ", "_") |> str_replace_all("/", "_")
    paste0("symptom_", x, "1")
  })) |>
  mutate(sx_interaction = str_glue("{sx1}:{sx2}")) |>
  label_sex() |>
  select(-sex) |>
  select(Gender = sex_lb, everything())

tab_interaction <-
  interaction_terms |>
  mutate(
    sx1 = str_split_i(term, "\\:", i = 1),
    sx2 = str_split_i(term, "\\:", i = 2)
  ) |>
  left_join(
    sx_terms |> select(Gender, cancer_site_desc, term, OR1 = OR, lower1 = lower, upper1 = upper),
    by = c("Gender", "cancer_site_desc", "sx1" = "term")
  ) |>
  left_join(
    sx_terms |> select(Gender, cancer_site_desc, term, OR2 = OR, lower2 = lower, upper2 = upper),
    by = c("Gender", "cancer_site_desc", "sx2" = "term")
  ) |>
  left_join(
    or_interaction |> select(Gender, cancer_site_desc, sx_interaction, OR_comb_inter = OR, lower_comb_inter = lower, upper_comb_inter = upper),
    by = c("Gender", "cancer_site_desc", "term" = "sx_interaction")
  ) |>
  left_join(
    or_combination |> select(Gender, cancer_site_desc, sx_interaction, OR_comb = OR, lower_comb = lower, upper_comb = upper),
    by = c("Gender", "cancer_site_desc", "term" = "sx_interaction")
  ) |>
  mutate(sx_label = paste(
    sx1, "and", sx2
  ) |>
    str_remove_all("symptom_") |>
    str_remove_all("1") |>
    str_replace_all("_", " ") |>
    str_to_sentence()
  )
tab_interaction <-
  tab_interaction |>
  mutate(OR_interact = if_else(p < 0.05 & !is.na(OR_comb_inter), OR, NA),
         lower_interact = if_else(p < 0.05 & !is.na(OR_comb_inter), lower, NA),
         upper_interact = if_else(p < 0.05 & !is.na(OR_comb_inter), upper, NA)
  ) |>
  select(Gender, cancer_site_desc, sx_label,
         OR1, lower1, upper1,
         OR2, lower2, upper2,
         p_interact = p,
         OR_interact, lower_interact, upper_interact,
         OR_comb_inter, lower_comb_inter, upper_comb_inter,
         OR_comb, lower_comb, upper_comb)
tab_interaction <-
  tab_interaction |>
  mutate(across(c(starts_with("OR"), starts_with("lower"), starts_with("upper")), function(x) x |> format_numb()),
         p_interact = p_interact |> format_p_val(lancet = T)) |>
  mutate(p_interact = p_interact |> str_remove("p=?")) |>
  mutate(`OR1 (95% CI)` = str_glue("{OR1} ({lower1}, {upper1})"),
         `OR2 (95% CI)` = str_glue("{OR2} ({lower2}, {upper2})"),
         `OR_inter (95% CI)` = if_else(is.na(OR_interact), NA, str_glue("{OR_interact} ({lower_interact}, {upper_interact})")),
         `OR_comb_inter (95% CI)` = if_else(is.na(OR_comb_inter), NA, str_glue("{OR_comb_inter} ({lower_comb_inter}, {upper_comb_inter})")),
         `OR_comb (95% CI)` = if_else(is.na(OR_comb), NA, str_glue("{OR_comb} ({lower_comb}, {upper_comb})"))
  ) |>
  mutate(sx_label = sx_label |> str_replace_all("luts nocturia frequency hesitancy urgency retention", "LUTS")) |>
  select(-c(OR1, lower1, upper1,
            OR2, lower2, upper2,
            OR_interact, lower_interact, upper_interact,
            OR_comb_inter, lower_comb_inter, upper_comb_inter,
            OR_comb, lower_comb, upper_comb)) |>
  select(Gender,
         `Cancer site` = cancer_site_desc,
         `Pairwise combination` = sx_label,
         starts_with("OR1"),
         starts_with("OR2"),
         p_interact,
         starts_with("OR_inter"),
         starts_with("OR_comb_inter"),
         starts_with("OR_comb")
  ) |>
  arrange(desc(Gender)) |>
  flextable() |>
  merge_v(j = c(1, 2)) |>
  valign(j = c(1, 2), valign = "top") |>
  add_table_theme() |>
  save_as_docx(path = stringr::str_glue("{save_dir}/App1_Table{table_numb}_OR_sx_comb_{run_date}.docx"))

table_numb <- table_numb + 1

rm(or_pairwise, sx_terms, interaction_terms, or_interaction)

# all symptoms model results ------
only_sx_terms <-
  model_tidy |>
  filter(grepl("symptom", term) & !grepl(":", term))

only_sx_terms <-
  only_sx_terms |>
  mutate(sx_label = term |>
           str_remove_all("symptom_") |>
           str_remove_all("1") |>
           str_replace_all("_", " ") |>
           str_to_sentence()
  ) |>
  label_sex() |>
  select(-sex) |>
  select(Gender = sex_lb, everything())

only_sx_terms <-
  only_sx_terms |>
  mutate(across(c(starts_with("OR"), starts_with("lower"), starts_with("upper")), function(x) x |> format_numb()),
         p_tidy = p |> format_p_val(lancet = T)) |>
  mutate(p_tidy = p_tidy |> str_remove("p=?"),
         `OR (95% CI)` = str_glue("{OR} ({lower}, {upper})")
  ) |>
  mutate(sx_label = case_when(
    sx_label == "luts nocturia frequency hesitancy urgency retention" ~ "LUTS",
    sx_label == "N a" ~ "N/A*",
    sx_label == "N k" ~ "N/K*",
    TRUE ~ sx_label)
  )

only_sx_terms |>
  select(Gender,
         `Cancer site` = cancer_site_desc,
         Symptom = sx_label,
         `OR (95% CI)`,
         p = p_tidy
  ) |>
  arrange(desc(Gender)) |>
  flextable() |>
  merge_v(j = c(1, 2)) |>
  valign(j = c(1, 2), valign = "top") |>
  add_table_theme() |>
  save_as_docx(path = stringr::str_glue("{save_dir}/App1_Table{table_numb}_OR_all_sx_{run_date}.docx"))

rm(model_tidy, only_sx_terms)

# figure -----
or_plot_data <-  read.csv(str_glue("{results_loc_main}model_symptom_combintion_to_plot_{results_generated_dt_main}.csv")) 

w_plot_data <- or_plot_data |> filter(sex == 2)
w_plot_data
w_plot_data <-
  w_plot_data |>
  mutate(sx_label = factor(sx_label, levels = rev(c(
    "Breast lump/mass",
    "Breast pain",
    "Nipple changes",
    "Cough",
    "Chest pain",
    "Dyspnoea",
    "Breast lump/mass and breast pain",
    "Breast lump/mass and nipple changes",
    "Cough and chest pain",
    "Cough and dyspnoea"
  ))))

w_or <-
  w_plot_data |>
  select(sex, cancer_site_desc, sx_label, with_interaction, OR, lower, upper) |>
  filter(sex == 2) |>
  mutate(sex_lb = case_when(sex == 1 ~ "Men", sex == 2 ~ "Women", TRUE ~ "Other")) |>
  mutate(with_interaction_label = if_else(is.na(with_interaction), "", 
                                          if_else(with_interaction, "with interaction", 
                                                  "without interaction"))) |>
  filter(!is.na(sx_label)) |>
  arrange(cancer_site_desc) |>
  ggplot(aes(x = OR, xmin = lower, xmax = upper, y = sx_label, colour = with_interaction_label, label = with_interaction_label))  +
  geom_point(show.legend = F, position = position_dodge(.7)) +
  geom_errorbar(width = .15, show.legend = F, position = position_dodge(.7)) + 
  geom_vline(xintercept = 1, colour = "darkgrey", linetype = "dashed") +
  geom_text(vjust = -0.5, show.legend = F, position = position_dodge(.7)) +
  facet_grid(cancer_site_desc ~ sex_lb, scales = "free_y") +
  scale_colour_manual(values = c("black", "darkblue", "brown")) +
  scale_x_log10(breaks = c(2^seq(-7, -1, by = 2), 0, 2^seq(1, 3, by = 2)), limits = c(2^-7, 2^1),
                labels = c(paste0("1/", 2^seq(7, 1, by = -2)), "1", 2^(seq(1, 3, by = 2)))
  ) +
  labs(colour = "Cancer Site", x = "Odds ratio", y = "") +
  theme_minimal() +
  theme(
    text             = element_text(size = 14, colour = "black"),
    panel.background = element_rect(fill = NA, colour = "black")
  ) 

ggsave(plot     = w_or,
       path     = save_dir,
       filename = str_glue("App1_Fig1a_sx_combinations_women_{run_date}.png"),
       width    = 27, 
       height   = 18, 
       unit     = "cm")

m_plot_data <- or_plot_data |> filter(sex == 1)

m_plot_data
m_plot_data <-
  m_plot_data |>
  mutate(sx_label = factor(sx_label, levels = rev(c(
    "Cough",
    "Haemoptysis",
    "Cough and haemoptysis"
  ))))

m_or <-
  m_plot_data |>
  select(sex, cancer_site_desc, sx_label, with_interaction, OR, lower, upper) |>
  filter(sex == 1) |>
  mutate(sex_lb = case_when(sex == 1 ~ "Men", sex == 2 ~ "Women", TRUE ~ "Other")) |>
  mutate(with_interaction_label = if_else(is.na(with_interaction), "", 
                                          if_else(with_interaction, "with interaction", 
                                                  "without interaction"))) |>
  filter(!is.na(sx_label)) |>
  arrange(cancer_site_desc) |>
  ggplot(aes(x = OR, xmin = lower, xmax = upper, y = sx_label, colour = with_interaction_label, label = with_interaction_label))  +
  geom_point(show.legend = F, position = position_dodge(.5)) +
  geom_errorbar(width = .15, show.legend = F, position = position_dodge(.5)) + 
  geom_vline(xintercept = 1, colour = "darkgrey", linetype = "dashed") +
  geom_text(vjust = -0.5, show.legend = F, , position = position_dodge(.5)) +
  facet_grid(cancer_site_desc ~ sex_lb, scales = "free_y") +
  scale_colour_manual(values = c("black", "darkblue", "brown")) +
  scale_x_log10(breaks = c(2^seq(-3, -1, by = 2), 0, 2^seq(1, 3, by = 2)), limits = c(2^-3, 2^1),
                labels = c(paste0("1/", 2^seq(3, 1, by = -2)), "1", 2^(seq(1, 3, by = 2)))
  ) +
  labs(colour = "Cancer Site", x = "Odds ratio", y = "") +
  theme_minimal() +
  theme(
    text             = element_text(size = 14, colour = "black"),
    panel.background = element_rect(fill = NA, colour = "black")
    
  ) 
ggsave(plot     = m_or,
       path     = save_dir,
       filename = str_glue("App1_Fig1b_sx_combinations_men_{run_date}.png"),
       width    = 26, 
       height   = 13, 
       unit     = "cm")

combined_or <- ((w_or + ggtitle("a)")) / (m_or + ggtitle("b)"))) + 
  plot_layout(heights = c(2, 1))

ggsave(plot     = combined_or,
       path     = save_dir,
       filename = str_glue("App1_Fig1_sx_combinations_{run_date}.png"),
       width    = 27, 
       height   = 40, 
       unit     = "cm")
