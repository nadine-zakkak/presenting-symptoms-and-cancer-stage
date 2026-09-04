# Formatting functions ----
# Helper function for suppression of numbers -----
# Code for suppression rules from Matt
# NZ added max_n and min_n
#' suppress_counts
#' Suppress counts based on restrictions on the data
#'
#' @param table: dataframe that contains values to be suppressed 
#' @param suppress_n: maximum value allowed before suppression
#' @param round_to: what value to round to 
#' @param n_var: column name of values to be suppressed 
#'
#' @returns: same table with 4 additional columns:
#' rules (suppressed or rounded, depending which was applied)
#' suppressed_n (the new value to be used)
#' min_n (the minimum possible value it could be based on the rules applied)
#' max_n (the maximum possible value it could be based on the rules applied)
suppress_counts <- function(table, suppress_n = 10, round_to = 5, n_var) {
  if(round_to*round(suppress_n/round_to) < suppress_n) {
    issue <- round_to*round(suppress_n/round_to)
    max_n_round_down <- issue + floor(round_to/2)
    print(str_glue("caution: rounded values may becomes less than suppressed. it is suppressed at {suppress_n}, but {max_n_round_down} will round to {issue}"))
  }
  table <-
    table |>
    mutate(suppressed = {{n_var}} < suppress_n) |>
    mutate(rounded = max(suppressed) |> as.logical()) |>
    mutate(rules = case_when(
      suppressed ~ "Suppressed",
      rounded    ~ "Rounded",
      TRUE       ~ ""
    ))
  
  table <-
    table |>
    mutate(
      suppressed_n = case_when(
        rules == "Suppressed" ~ suppress_n - 1,
        rules == "Rounded"    ~ round_to*round({{n_var}}/round_to, digits = 0),
        TRUE ~ {{n_var}}
      )
    )  |>
    # set min n for
    ## suppressed to be 0
    ## rounded to be the min number it could have been given what it has been rounded to, except if it has been rounded to the suppressed value then set as 0
    ## otherwise keep as it is
    mutate(min_n = case_when(
      rules == "Suppressed" ~ 0,
      rules == "Rounded"    ~ if_else(suppressed_n < suppress_n, 0, suppressed_n - floor(round_to/2)),
      TRUE ~ {{n_var}}
    )) |>
    # set max n for
    ## suppressed to be the number to be the max suppressed value
    ## rounded to be the max number it could have been given what it has been rounded to, except if it has been rounded to the suppressed value then set as max suppressed value
    ## otherwise keep as it is
    mutate(max_n = case_when(
      rules == "Suppressed" ~ suppress_n - 1,
      rules == "Rounded"    ~ if_else(suppressed_n < suppress_n, suppressed_n, suppressed_n + floor(round_to/2)),
      TRUE ~ {{n_var}}
    )) |>
    # fix an negative values that could arise and set to 0
    mutate(max_n = if_else(max_n < 0, 0, max_n),
           min_n = if_else(min_n < 0, 0, min_n)
    )
  
  table |> select(-suppressed, -rounded)
}

#' calculate_suppressed_total
#' Find the total values based on the suppression rules applied
#'
#' @param table: data that includes the suppressed value columns 
#' @param suppressed_n_var: column name that contains the suppressed value 
#' @param min_n_var: column name contains the minimum possible value it could have been 
#' @param max_n_var: column name contains the maximum possible value it could have been 
#'
#' @returns: same table with 3 additional columns:
#' 1. suppressed_N: total count based on the suppresed values
#' 2. min_N: minimum total count based on the minimum possible values
#' 3. max_N: maximum total count based on the maximum possible values
calculate_suppressed_total <- function(table, suppressed_n_var, min_n_var, max_n_var) {
  table |>
    mutate(suppressed_N = sum({{suppressed_n_var}}),
           min_N        = sum({{min_n_var}}),
           max_N        = sum({{max_n_var}}))
}

#' calculate_suppressed_prop_with_ci
#' Calculate the proportions and 95% CI using the new suppressed values
#'
#' @param table: dataframe that includes all suppressed value columns 
#' @param n_var: column name of the original values 
#' @param suppressed_n_var: column name of the suppressed values 
#' @param min_n_var: column name of the minimum suppressed values
#' @param max_n_var : column name of the maximum suppressed values
#' @param suppressed_N_var : column name of the total count of suppressed values
#' @param min_N_var : column name of the total count of the least suppressed values
#' @param max_N_var : column name of the total count of the highest suppresed values
#'
#' @returns : same dataframe with 3 additional columns
#' 1. prop100_suppressed (proportion based onthe suppressed values)
#' 2. prop100_lb_min (proportion based on the least possible 95% LB)
#' 3. prop100_ub_max (proportion based on the max possible 95% UB)
calculate_suppressed_prop_with_ci <- function(table, 
                                              n_var,
                                              suppressed_n_var, 
                                              min_n_var, 
                                              max_n_var, 
                                              suppressed_N_var,
                                              min_N_var,
                                              max_N_var) {
  # calculate the proportion using the new suppressed value as numerator
  # calculate the lower 95% CI using the minimum possible value the number could have been (min numerator and max denominator)
  # calculate the upper 95% CI using the maximum possible value the number could have been (max numerator and min denominator)
  table |> 
    mutate(
      prop100_suppressed = 100*{{suppressed_n_var}}/{{suppressed_N_var}},
      prop100_lb_min = 100*calculate_ci(n = {{min_n_var}}, N = {{max_N_var}})$lower,
      prop100_ub_max = 100*calculate_ci(n = {{max_n_var}}, N = {{min_N_var}})$upper
    ) |>
    mutate(prop100_lb_min = if_else(is.na(prop100_lb_min) | prop100_lb_min < 0, 0, prop100_lb_min),
           prop100_ub_max = if_else(is.na(prop100_ub_max) | prop100_ub_max > 100, 100, prop100_ub_max))
}

#' tidy_suppressed_vals
#' Tidy the suppressed values
#'
#' @param table : dataframe that includes the suppressed count column  
#' @param suppressed_n_var : column name that includes the suppressed counts
#'
#' @returns : table with new column suppressed_n_tidy that has tidied counts
tidy_suppressed_vals <- function(table, suppressed_n_var) {
  
  table |>
    mutate(suppressed_n_tidy = if_else(rules == "Suppressed", 
                                       ({{suppressed_n_var}} + 1), 
                                       {{suppressed_n_var}}) |>
             as.character()
    ) |>
    mutate(suppressed_n_tidy = if_else(rules == "Suppressed", paste0("<", suppressed_n_tidy), suppressed_n_tidy))
}

#' suppress_values_with_prop
#' Main function that calls all other helper function to 
#' find the suppressed counts, along with their proportions and 95% CI
#'
#' @param table : dataframe that includes the column that has values to be suppressed
#' @param suppress_n : maximum value allowed before suppression
#' @param round_to : what value to round to 
#' @param n_var : column name of values to be suppressed 
#'
#' @returns : same dataframe with all additional columns relating to suppression and proportions
#' suppressed_n, min_n, max_n, suppressed_n_tidy
#' suppressed_N, min_N, max_N
#' prop100_suppressed, prop100_lb_min, prop100_ub_max
suppress_values_with_prop <- function(table, suppress_n = 10, round_to = 5, n_var) {
  table |>
    suppress_counts(suppress_n = suppress_n, 
                    round_to   = round_to, 
                    n_var      = {{n_var}}) |>
    calculate_suppressed_total(suppressed_n_var = suppressed_n, 
                               min_n_var        = min_n, 
                               max_n_var        = max_n) |>
    calculate_suppressed_prop_with_ci(n_var = {{n_var}},
                                      suppressed_n_var = suppressed_n, 
                                      min_n_var = min_n, 
                                      max_n_var = max_n, 
                                      suppressed_N_var = suppressed_N, 
                                      min_N_var = min_N, 
                                      max_N_var = max_N) |>
    tidy_suppressed_vals(suppressed_n_var = suppressed_n)
}

#' format_numb
#'
#' @param numb - number to format 
#'
#' @return formatted number to include 3 digits. if <10: 2 decimal places, (10,100): 1 decimal place and >100 0 decimal places
#'
#' @examples: format_numb(0.128) --> 0.13, format_numb(90.4736) --> 90.5, format_numb(123.2737) --> 123
format_numb <- function(numb) {
  numb = ifelse(numb < 10,
                sprintf(numb, fmt='%.2f'),
                ifelse(numb < 100,
                       sprintf(numb, fmt = '%.1f'),
                       sprintf(numb, fmt = '%1.0f')))
}

#' format_ci
#'format prop and CI into a single variable in the format prop (prop_lb, prop_ub)
#' or (prop_lb, prop_ub)
#' 
#' @param data : dataframe with the 3 columns prop100, prop100_lb, prop100_ub
#' @param with_prop : to include the proportion?
#'
#' @returns : same dataframe with 1 additional column %_(95% CI) or (95% CI)
format_ci <- function(data, with_prop = T) {
  if(with_prop) {
    data |>
      mutate(`%_(95% CI)` = ifelse(is.na(prop100), 
                                   "-",
                                   str_glue("{prop100}% ({prop100_lb}, {prop100_ub})"))
      )
  } else {
    data |>
      mutate(`(95% CI)` = ifelse(is.na(prop100), 
                                 "-",
                                 str_glue("({prop100_lb}, {prop100_ub})"))
      )
  }
}

#' format_p_val
#' Format the p-values in a tidy format for publishing
#'
#' @param p_val : the p-value to be formatted
#' @param lancet : whether to use the lancet formatting:
#' p < 0.0001 and if larger, than format to 2 significant figures
#'
#' @returns the formatted p-value
format_p_val <- function(p_val, lancet = F) {
  if(!lancet) {
    ifelse(p_val < 0.001, "p<0.001", paste0("p=", sprintf(p_val, fmt = '%.3f')))
  } else {
    ifelse(p_val < 0.0001, "p<0.0001", paste0("p=", p_val |> signif(digits = 2) |> formatC(digits = 2, format = "fg", flag = "#")))
  }
}

#' label_sex
#' Add labels to Men and Women
#' assumes 1 = Men and 2 = Women and there are no other options in this instance
#' 
#' @param data : dataframe that includes the column 'sex' as numeric 
#'
#' @returns : same dataframe with additional column 'sex_lb'
label_sex <- function(data) {
  data |> 
    mutate(sex_lb = case_when(sex == 1 ~ "Men", sex == 2 ~ "Women"))
}

#' label_stage
#' Add labels to stage to be consistent throughout
#' 
#' @param data : dataframe that includes the relevant stage column (stage_bin or stage_recode)
#' assumes stage_bin includes levels 1-3, 4 and -99
#' and stage_recode includes levels 1, 2, 3, 4 and -99
#' @param binary : T/F set to T if stage to be labelled is dichotomised
#'
#' @returns : same dataframe with an additional column stage_bin_lb or stage_lb
label_stage <- function(data, binary = T) {
  if(binary){
    data |>
      mutate(stage_bin_lb = factor(stage_bin, levels = c("1-3", "4", "-99"), labels = c("1-3", "4", "Missing")))
  }
  else{
    data |>
      mutate(stage_lb = factor(stage_recode, levels = c(1:4, -99), labels = c(1:4, "Missing")))
  }
}

#' write_to_excel
#'
#' @param df : table to be written to excel sheet
#' @param sheet_name : sheet name in excel
#' @param excel_file_name : excel file name to save to
#'
#' @returns : nothing. data written to specified excel file
write_to_excel <- function(df, sheet_name, excel_file_name){
  if(is_tibble(df)) df <- df |> as.data.frame()
  write.xlsx(df, 
             file      = excel_file_name, 
             sheet     = sheet_name,
             row.names = F, 
             append    = T
  )
}
