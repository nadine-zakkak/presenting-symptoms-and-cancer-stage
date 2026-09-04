# Script to run through entire analysis based on Configs file

# Open log file in rw format ("a+b")
log_file <- file("run-complete-analysis_v1.log", open = "a+b")
# sink to log file all outputs (print statements) and all error/warning messages
sink(log_file, append = T, split = T, type = "output")
sink(log_file, append = T,  type ="message")

print("-----------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------")
print(paste0("Start date: ", as.character(Sys.time())))
print(R.version$version.string)

run_date <- as.character(Sys.Date())
saveRDS(run_date, "./run_date.rds")

source("./scripts/common_scripts/00_base_setup.R")

# Print configs file with its current conditions
print(paste0("Preprocess regenerate: "            , configs$preprocess$regenerate))
print(paste0("Main preprocess regenerate: "       , configs$main_configs$regenerate))
print(paste0("Main analysis regenerate: "         , configs$main_configs$analysis_regenerate))
print(paste0("Sensitivity preprocess regenerate: ", configs$sa_configs$regenerate))
print(paste0("Sensitivity analysis regenerate: "  , configs$sa_configs$analysis_regenerate))

# Common ----
{
  ## update common preprocessing, if needed ----
  if(configs$preprocess$regenerate) {
    print(paste0("*****", current_script, ": sourcing 1a_preprocessing.R, 1b_exclusions.R and 1c_create_analysis_data.R*****"))
    source("./scripts/common_scripts/1a_preprocessing.R")
    source("./scripts/common_scripts/1b_exclusions.R")
    source("./scripts/common_scripts/1c_create_analysis_data.R")
    source("./scripts/common_scripts/00_base_setup.R")
  }
}

# Main analysis ----
{
  ## update main analysis preprocessing, if needed  -----
  if(configs$main_configs$regenerate) {
    print(paste0("*****", current_script, ": sourcing ", main_loc, "1d_preprocessing.R*****"))
    source(paste0("./scripts/", main_loc, "/1d_preprocessing.R"))
    source("./scripts/common_scripts/00_base_setup.R")
  }
  
  ## update main analysis, if needed  -----
  if(configs$main_configs$analysis_regenerate) {
    print(paste0("*****", current_script, ": sourcing ",  main_loc, "/2_descriptive.R*****"))
    source(paste0("./scripts/", main_loc, "/2_descriptive.R"))
    source("./scripts/common_scripts/00_base_setup.R")

    print(paste0("*****", current_script, ": sourcing ",  main_loc, "/3a_analysis_stage_by_symptom.R*****"))
    source(paste0("./scripts/", main_loc, "/3a_analysis_stage_by_symptom.R"))
    source("./scripts/common_scripts/00_base_setup.R")

    print(paste0("*****", current_script, ": sourcing ",  main_loc, "/3b_analysis_stage_by_symptom_site.R*****"))
    source(paste0("./scripts/", main_loc, "/3b_analysis_stage_by_symptom_site.R"))
    source("./scripts/common_scripts/00_base_setup.R")

    print(paste0("*****", current_script, ": sourcing ",  main_loc, "/3c_analysis_stage_by_symptom_bt.R*****"))
    source(paste0("./scripts/", main_loc, "/3c_analysis_stage_by_symptom_bt.R"))
    source("./scripts/common_scripts/00_base_setup.R")

    # print("**** SKIPPING SYMPTOM ANALYSIS COMBINATIONS AS NO CHANGES IN THIS SCRIPT")
    print("**** RUNNING SYMPTOM ANALYSIS COMBINATIONS AS WAS SKIPPED IN PREVIOUS RUN")
    print(paste0("*****", current_script, ": sourcing ",  main_loc, "/3d_analysis_symptom_combinations.R*****"))
    source(paste0("./scripts/", main_loc, "/3d_analysis_symptom_combinations.R"))
    source("./scripts/common_scripts/00_base_setup.R")

    print(paste0("*****", current_script, ": sourcing ",  main_loc, "/3e_analysis_stage_by_symptom_multiplicity.R*****"))
    source(paste0("./scripts/", main_loc, "/3e_analysis_stage_by_symptom_multiplicity.R"))
    source("./scripts/common_scripts/00_base_setup.R")
  }
}

# Sensitivity analyses ---- 
{
  ## update sensitivity analysis preprocessing, if needed  -----
  if(configs$sa_configs$regenerate) {
    print(paste0("*****", current_script, ": sourcing ", sa_miss_loc, "/1d_preprocessing.R*****"))
    source(paste0("./scripts/", sa_miss_loc, "/1d_preprocessing.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    print(paste0("*****", current_script, ": sourcing ", sa_prst_loc, "/1d_preprocessing.R*****"))
    source(paste0("./scripts/", sa_prst_loc, "/1d_preprocessing.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    print(paste0("*****", current_script, ": sourcing ", sa_rtd_loc, "/1d_preprocessing.R*****"))
    source(paste0("./scripts/", sa_rtd_loc, "/1d_preprocessing.R"))
    source("./scripts/common_scripts/00_base_setup.R")
  }
  
  ## update sensitivity analysis, if needed  -----
  if(configs$sa_configs$analysis_regenerate) {
    # missing stage analysis
    print(paste0("*****", current_script, ": sourcing ", sa_miss_loc, "/2_descriptive.R*****"))
    source(paste0("./scripts/", sa_miss_loc, "/2_descriptive.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    print(paste0("*****", current_script, ": sourcing ", sa_miss_loc, "/3a_analysis_stage_by_symptom.R*****"))
    source(paste0("./scripts/", sa_miss_loc, "/3a_analysis_stage_by_symptom.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    print(paste0("*****", current_script, ": sourcing ", sa_miss_loc, "/3b_analysis_stage_by_symptom_site.R*****"))
    source(paste0("./scripts/", sa_miss_loc, "/3b_analysis_stage_by_symptom_site.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    print(paste0("*****", current_script, ": sourcing ", sa_miss_loc, "/3c_analysis_stage_by_symptom_bt.R*****"))
    source(paste0("./scripts/", sa_miss_loc, "/3c_analysis_stage_by_symptom_bt.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    # excl prostate cancer
    print(paste0("*****", current_script, ": sourcing ", sa_prst_loc, "/2_analysis.R*****"))
    source(paste0("./scripts/", sa_prst_loc, "/2_analysis.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
    # by RtD
    print(paste0("*****", current_script, ": sourcing ", sa_rtd_loc, "/2_analysis.R*****"))
    source(paste0("./scripts/", sa_prst_loc, "/2_analysis.R"))
    source("./scripts/common_scripts/00_base_setup.R")
    
  }
}

# Table figures and manuscripts ---- 
{
  print(paste0("*****", current_script, ": sourcing ", manuscript_figs_loc, "/main_figures*****"))
  source(paste0("./scripts/", manuscript_figs_loc, "/main_figures.R"))
  source("./scripts/common_scripts/00_base_setup.R")  
  
  print(paste0("*****", current_script, ": sourcing ", manuscript_figs_loc, "/main_tables*****"))
  source(paste0("./scripts/", manuscript_figs_loc, "/main_tables.R"))
  source("./scripts/common_scripts/00_base_setup.R") 
  
  print(paste0("*****", current_script, ": sourcing ", manuscript_figs_loc, "/lancet_appendix_tables_suppressed*****"))
  source(paste0("./scripts/", manuscript_figs_loc, "/lancet_appendix_tables_suppressed.R"))
  source("./scripts/common_scripts/00_base_setup.R")  
  
  print(paste0("*****", current_script, ": sourcing ", manuscript_figs_loc, "/appendix_symptom_combinations*****"))
  source(paste0("./scripts/", manuscript_figs_loc, "/appendix_symptom_combinations.R"))
  source("./scripts/common_scripts/00_base_setup.R")  
  
  print(paste0("*****", current_script, ": sourcing ", manuscript_figs_loc, "/appendix_rtd*****"))
  source(paste0("./scripts/", manuscript_figs_loc, "/appendix_rtd.R"))
  source("./scripts/common_scripts/00_base_setup.R")
  
  print(paste0("*****", current_script, ": sourcing ", manuscript_figs_loc, "/text outputs*****"))
  knitr::knit("./scripts/manuscript_results_outputs/text_ouputs.Rmd", output = "./scripts/manuscript_results_outputs/text_outputs.docx")
  source("./scripts/common_scripts/00_base_setup.R")
}

print(paste0("End date: ", as.character(Sys.time())))
# save all outputs to log file
sink(append = T)
sink(append = T, type = "message")

