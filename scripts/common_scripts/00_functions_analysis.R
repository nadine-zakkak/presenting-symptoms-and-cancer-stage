#' Counts number of patients and records
#'
#' @param df: data (dataframe)
#' @param step: step of eligibility criteria (string)
#'
#' @return: dataframe with the symptom description, number of records, number of unique patients and description of step of eligibility criteria
count_tumours_and_patients <- function(df, step = NA) {
  df <-
    df |>
    summarise(n_tumours  = n_distinct(tumourid),
              n_patients = n_distinct(patid)) 
  
  if(!is.na(step)) {
    df <- df |>
      mutate(description = step) |>
      ungroup()
  }
  return(df)
}

#' calculate_ci (95% CI) - from Matt
#'
#' @param n - number of samples (numerator)
#' @param N - total number of samples (denominator)
#'
#' @return: a list with 2 values: "lower" = lower CI, "upper" = upper CI
calculate_ci <- function(n, N) {
  lower = (1/(1+(qnorm(0.975)^2)/N))*((n/N)+(qnorm(0.975)^2)/(2*N)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/N)))*sqrt((n/N)*(1-n/N)/N + (qnorm(0.975)^2)/(4*(N^2)))
  upper = (1/(1+(qnorm(0.975)^2)/N))*((n/N)+(qnorm(0.975)^2)/(2*N)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/N)))*sqrt((n/N)*(1-n/N)/N + (qnorm(0.975)^2)/(4*(N^2)))
  ls <- list(lower, upper)
  names(ls) <- c("lower", "upper")
  return(ls)
}

#' find_prop_tumours_by_site
#' Finds the proportion of tumours stratified by cancer site
#'
#'
#' @param df with cancer_site_desc and tumourid columns 
#'
#' @returns: summarised dataframe with number of tumours per cancer site 
find_prop_tumours_by_site <- function(df) {
  df |>
    group_by(cancer_site_desc) |>
    summarise(n = n_distinct(tumourid)) |>
    ungroup() |>
    mutate(prop = 100*n/sum(n))
}

#' find_n_patients_by_tumour_count
#' Calculates the number of patients by tumour multiplicity
#'
#'
#' @param df with patid and tumourid columns
#'
#' @returns: summarised dataframe with number of tumours per patient 
find_n_patients_by_tumour_count <- function(df) {
  data |>
    group_by(patid) |>
    summarise(n_tumours = n_distinct(tumourid)) |>
    group_by(n_tumours) |>
    summarise(n_pat = n_distinct(patid)) |>
    mutate(total_pat = sum(n_pat))
}

#' flag_multiple_tumours
#' Add flag if patients had multiple tumours
#'
#' @param df: dataframe with patid, tumourid and sex columns
#' @param cohort_data: dataframe with demographic info (patid, tumourid, imd and age)
#'
#' @returns: dataframe with column flagging if patient had multiple tumours
flag_multiple_tumours <- function(df, cohort_data) {
  # flag patients if they had more than 1 tumour
  df |> 
    add_count(patid) |>
    mutate(mult_tumour = n > 1) |>
    distinct(patid, tumourid, sex, mult_tumour) |>
    inner_join(cohort_data |> 
                 select(patid, tumourid, 
                        imd = quintile_2019, 
                        age = tenyearageband), 
               by = c("patid", "tumourid"))
}

#' Helper function for patients' characteristics summaries
#'
#' @param data : dataframe to tidy 
#' @param variable: which variable to tidy (options are sex, imd and age)
#'
#' @returns
cleaning_var <- function(data, variable) {
  if (variable == "sex") {
    data |> label_sex() |> mutate(val = sex_lb, var = "Sex")
  }
  else if (variable == "imd") {
    data |> mutate(val = imd, var = "IMD")
  }
  else if (variable == "age") {
    data |> mutate(val = age, var = "Age")
  }
}

#' count_by_mult_tumour
#  Count the number of patients and tumours by tumour multiplicity and variable of interest
#'
#' @param data 
#' @param variable: variable to group by, to count number of patients and tumours
#'
#' @returns: summarised data frame with number of patients and tumours by
#' tumour multiplicity and variable of interest
count_by_mult_tumour <- function(data = data_tumour_count, variable) {
  variable_str <- variable |> substitute() |> deparse()
  data |>
    group_by(mult_tumour, {{variable}}) |>
    count_tumours_and_patients()|>
    cleaning_var(variable_str)
}

#' #' count_by_mult_tumour
#  Count the number of patients and tumours by the variable of interest
#'
#' @param data 
#' @param variable: variable to group by, to count number of patients and tumours 
#'
#' @returns: summarised data frame with number of patients and tumours by
#' variable of interest
count_any <- function(data = data_tumour_count, variable) {
  variable_str <- variable |> substitute() |> deparse()
  data |>
    group_by({{variable}}) |>
    count_tumours_and_patients() |>
    cleaning_var(variable_str)
}

#' symptom_stage_count
#' Find the proportion of stage stratified by sex and symptom
#'
#' @param data 
#' @param stage_type: stage var to use (binary or detailed)? 
#'
#' @returns: summarised data frame with number and proportion of tumours
#' for each of the stage level by sex and symptom
symptom_stage_count <- function(data, stage_type) {
  data |>
    group_by(sex, symptom, {{stage_type}}) |>
    summarise(n_tumour = sum(present)) |>
    inner_join(N_tumour_by_sex_symptom, by = c("sex", "symptom")) |>
    mutate(
      prop100    = 100*n_tumour/N_tumour,
      prop100_lb = 100*calculate_ci(n_tumour, N_tumour)$lower,
      prop100_ub = 100*calculate_ci(n_tumour, N_tumour)$upper
    ) |>
    label_sex() |>
    ungroup()
}

#' symptom_inv_stage_count
#' Find the proportion of stage stratified by
#' sex, symptom and investigation type
#'
#' @param data 
#' @param stage_type: specify which stage variable is required 
#'
#' @returns: summarised data frame with number and proportion of tumours 
#' for each of the stage level by sex, symptom and investigation type
symptom_inv_stage_count <- function(data, stage_type) {
  data |>
    group_by(sex, symptom, inv, present_inv, {{stage_type}}) |>
    summarise(n_tumour = sum(present_sx)) |> 
    inner_join(N_tumour_by_sex_symptom_invt, 
               by = c("sex", "symptom", "inv", "present_inv")) |>
    mutate(prop100    = 100*n_tumour/N_tumour,
           prop100_lb = 100*calculate_ci(n_tumour, N_tumour)$lower,
           prop100_ub = 100*calculate_ci(n_tumour, N_tumour)$upper) |>
    label_sex() |>
    ungroup()
}

#' find_stage_prop_by_site_and_sx
#' Find the proportion of dichotomised stage stratified by
#' symptom and cancer site
#' for filtered symptoms and cancer sites
#'
#' @param data 
#' @param symptoms 
#' @param cancer_sites 
#'
#' @returns: summarised data frame with number and proportion of tumours 
#' for each of the stage level by symptom and cancer site
find_stage_prop_by_site_and_sx <- function(data, symptoms, cancer_sites, stage_type) { 
  
  data <- data |> filter(symptom %in% symptoms,
                         cancer_site_desc %in% cancer_sites)
  
  N_tumour_by_symptom_site_exemplar <-
    data |> 
    group_by(symptom, cancer_site_desc) |>
    summarise(N_tumour = sum(present)) |>
    arrange(cancer_site_desc, desc(N_tumour)) |>
    ungroup()
  
  # Using stage, find the proportion of tumours with each stage for symptom and cancer site groups
  data |> 
    group_by(symptom, cancer_site_desc, {{stage_type}}) |> 
    summarise(n_tumour = sum(present)) |>
    ungroup() |>
    left_join(N_tumour_by_symptom_site_exemplar, by = c("symptom", "cancer_site_desc")) |>
    mutate(prop100     = 100*n_tumour/N_tumour,
           prop100_lb  = 100*calculate_ci(n_tumour, N_tumour)$lower,
           prop100_ub  = 100*calculate_ci(n_tumour, N_tumour)$upper) |>
    arrange(symptom, desc(prop100)) |>
    ungroup()
}

#' find_n_elig_sx
#' Find the number of eligible symptoms stratified by sex:
#' eligible symptoms common to both, only men, only women
#'
#' @param df 
#' @param cut_off: minimum number of tumours to be eligible 
#'
#' @returns: dataframe with the number of eligible symptoms common 
#' to both men and women and each separately
find_n_elig_sx <- function(df, cut_off = 100) { 
  symptoms_both <-   
    df |>
    filter(N_tumour >= cut_off) |>
    count(symptom) |>
    filter(n > 1) |>
    mutate(sex_lb = "Both") |>
    select(-n)
  
  symptoms_both |>
    bind_rows(
      df |>
        filter(N_tumour >= cut_off) |>
        anti_join(symptoms_both |> distinct(symptom), by = "symptom") |>
        label_sex() |>
        distinct(symptom, sex_lb) 
    ) |>
    group_by(sex_lb) |>
    mutate(numb_symptoms = row_number()) |>
    ungroup() |>
    select(sex_lb, symptom, numb_symptoms)
}

#' count_tumour_by_sx_numb
#' Count the number of tumours by symptom multiplicity
#'
#' @param df 
#' @param numb_sx_var: symptom multplicity column to group by 
#'
#' @returns: summarised dataframe with the number of tumours by symptom count
count_tumour_by_sx_numb <- function(df, numb_sx_var) {
  df  |>
    group_by({{numb_sx_var}}) |> 
    summarise(n_tumour = n()) |> 
    mutate(prop = 100*n_tumour/sum(n_tumour))
}


#' run_permutation_test
#' Function to run permutation tests of a t-test of advanced stage and blood test presenced
#' stratified by sex and symptom
#' 
#' @param data: dataframe that includes data to be used for the t-test
#' should include sex_lb, tumourid, stage_advanced (0/1) and blood_test (0/1)
#' @param n_permut : number of permutations
#' @param method : method to be used, default it t-test, no other options are currently available
#'
#' @returns: a list with the outputs of the t-tests
run_permutation_test <- function(data, n_permut, method = c("t-test")) {
  # very primitive checks to make sure "input" is reasonable
  if(!method %in% c("t-test") | n_permut < 0) {
    print("Invalid method or number of permutations")
  }
  else {
    run_test <- function(data, permut_numb, method) {
      
      # function to randomly assign blood test occurrence
      shuffle_labels <- function(data) {
        data |>
          group_by(tumourid) |>
          mutate(blood_test = sample(0:1, size = 1)) |>
          ungroup()
      }
      
      # shuffle blood test labels if not "base" round of permutation
      if(permut_numb > 0) data <- shuffle_labels(data)
      
      if(method == "t-test") {
        stats_df <- 
          data |>
          group_by(sex_lb, symptom) |>
          nest() |>
          mutate(t_stat = map(data, ~t.test(stage_advanced ~ blood_test, data = .x)),
                 stat   = map(t_stat, glance)) |>
          unnest(stat) |>
          select(-data, -t_stat) |>
          ungroup() |>
          mutate(permutation = permut_numb) |>
          select(sex_lb, symptom, t = statistic, permutation)
      }
      if(permut_numb %% (n_permut/10) == 0) print(paste0("******** permutations complete:", permut_numb, " ********"))
      
      return(stats_df)
    }
    return(lapply(0:n_permut, FUN = function(x) run_test(data, x, method)))
  }
}
