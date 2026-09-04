# Helper functions to produce final tables 

# Add table theme (using flextable - need to make sure it is loaded)
add_table_theme <- function(table) {
  table |> 
    fontsize(part = "all", size = 9) |>
    bold(part     = "header") |>
    font(part     = "all", fontname =  "Calibri") |>
    align(align = 'center', part = 'header') |>
    set_table_properties(layout = "autofit") 
}


# Stage distribution by symptom ----
prepare_table_stage_by_symptom <- function(data, df_sx_order, binary = F) {
  # function to produce tables of stage distribution by symptom
  data <- 
    data |> 
    label_sex() |>
    label_stage(binary = binary) |>
    filter(N_tumour >= 100)
  
  if(binary) {
    data <- data |> rename(stage_lb = stage_bin_lb)
  }
  
  # get symptom order for each sex
  order_men <- df_sx_order[["Men"]]
  order_men <- c("N/A", "N/K", setdiff(order_men, c("N/A", "N/K")))
  
  order_women <- df_sx_order[["Women"]]
  order_women <- c("N/A", "N/K", setdiff(order_women, c("N/A", "N/K")))
  
  prep_table <- function(data) {
    data |>
      select(Sex = sex_lb, Symptom, N_tumour, stage_lb, n_tumour, starts_with("prop100")) |>
      mutate(across(starts_with("prop100"), ~format_numb(.x))) |>
      format_ci(with_prop = F) |>
      rename(Gender  = Sex,
             `Proportion (%)`     = prop100,
             `Total tumour count` = N_tumour,
             `Tumour count`       = n_tumour) |>
      pivot_wider(id_cols     = c(Gender, Symptom, `Total tumour count`), 
                  names_from  = stage_lb, 
                  values_from = c(`Tumour count`, `Proportion (%)`, `(95% CI)`), 
                  names_vary  = "slowest", 
                  names_glue  = "{.value} Stage {stage_lb}")
  }
  
  w <-
    data |>
    filter(sex_lb == "Women") |> 
    rename(Symptom = symptom) |>
    filter(Symptom %in% order_women) |>
    prep_table() |>
    mutate(Symptom = factor(Symptom, levels = order_women)) |>
    arrange(desc(Symptom)) 
  
  m <-
    data |>
    filter(sex_lb == "Men") |>
    rename(Symptom = symptom) |>
    filter(Symptom %in% order_men) |>
    prep_table() |>
    mutate(Symptom = factor(Symptom, levels = order_men)) |>
    arrange(desc(Symptom)) |>
    select(-starts_with("prop100"))
  
  w |> bind_rows(m)
}

# Exemplar cancer site and symptoms ----
prepare_table_stage_by_symptom_exemplar <- function(data, binary = F) {
  # function to produce tables of stage distribution by symptom for specific sites and symptoms
  # Tidy up the data and present in wide format to be written to excel file 
  
  cancer_sites <- names(cancer_redflag)
  cancer_sites_alarm_sx <- 
    cancer_redflag |>
    bind_rows() |> 
    pivot_longer(cols      = everything(), 
                 names_to  = "cancer_site", 
                 values_to = "symptom") |> 
    distinct() |> 
    arrange(cancer_site, symptom)
  
  # for organ-specific alarm symptoms, only keep relevant cancer sites
  # Filter data to specific cancer site and order symptoms
  data_alarm <- 
    data |>
    inner_join(cancer_sites_alarm_sx, by = c("cancer_site_desc" = "cancer_site", "symptom" = "symptom")) |>
    mutate(symptom = factor(symptom)) |>
    arrange(symptom) 
  
  # for vague symptoms, keep all required cancer sites and all vague symptoms
  data_vague <-
    data |>
    filter(cancer_site_desc %in% cancer_sites &
             symptom %in% nonspecific_symptoms) |>
    mutate(symptom = factor(symptom)) |>
    arrange(symptom) 
  
  data_updated <- data_alarm |> bind_rows(data_vague)
  
  data_summary <-
    data_updated |> 
    label_stage(binary = binary)
  
  if(binary) {
    data_summary <- data_summary |> rename(stage_lb = stage_bin_lb)
  }
  
  data_summary <- data_summary |> arrange(stage_lb)
  
  data_summary |>
    select(Symptom = symptom, 'Cancer Site' = cancer_site_desc, N_tumour, stage_lb, n_tumour, starts_with("prop100")) |>
    mutate(across(starts_with("prop100"), ~format_numb(.x))) |>
    format_ci(with_prop = F) |>
    rename(`Proportion (%)`     = prop100,
           `Total tumour count` = N_tumour,
           `Tumour count`       = n_tumour) |>
    pivot_wider(id_cols     = c(Symptom, `Cancer Site`, `Total tumour count`), 
                names_from  = c(stage_lb), 
                values_from = c(`Tumour count`, `Proportion (%)`, `(95% CI)`), 
                names_vary  = "slowest", 
                names_glue  = "{.value} Stage {stage_lb}") |>
    mutate(Symptom = factor(Symptom, levels = c(unique(unlist(cancer_redflag)), 
                                                "Abdominal pain (NOS)",
                                                "Nausea and/or vomiting",
                                                "Weight loss")),
           `Cancer Site` = factor(`Cancer Site`, levels = names(cancer_redflag))) |>
    arrange(Symptom, `Cancer Site`, `Proportion (%) Stage 4`)
}

# Table of proportion of stage by blood test status ----
prepare_table_bt <- function(data, N_tumour_df, filter_by_sex, binary = T) {
  # function to produce tables of advanced stage proportion in symptoms by blood test status
  
  data <- 
    data |>
    label_stage(binary = binary) |>
    mutate(across(starts_with("prop100"), ~format_numb(.x))) |>
    format_ci(with_prop = F) |>
    mutate(`Blood test flag` = ifelse(blood_test, "Blood test", "No blood test")) |>
    rename(`Proportion (%)`     = prop100,
           `Total tumour count` = N_tumour,
           `Tumour count`       = n_tumour) |>
    pivot_wider(id_cols     = c(sex_lb, symptom, `Blood test flag`, `Total tumour count`),
                names_from  = c(stage_bin_lb), 
                values_from = c(`Tumour count`, `Proportion (%)`, `(95% CI)`), 
                names_glue  = ("{.value} Stage {stage_bin_lb}"),
                names_vary  = "slowest",
                names_sort  = F
    ) |>
    left_join(N_tumour_sex_symptom, by = c("sex_lb", "symptom")) |>
    select(sex_lb, symptom, Total, everything()) |>
    filter(Total >= 100)
  
  data |>
    filter(sex_lb == filter_by_sex) |>
    arrange(desc(Total)) |>
    mutate(symptom = factor(symptom, levels = c(setdiff(unique(symptom), c("N/A", "N/K")), c("N/A", "N/K")))) |>
    arrange(symptom) |>
    select(-Total) |>
    rename(Gender     = sex_lb, 
           Symptom = symptom)
}


