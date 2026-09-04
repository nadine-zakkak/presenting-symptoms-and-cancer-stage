# function used to load packages current session while suppressing warnings and start up messages
load_package <- function(pckg) {
  suppressWarnings(suppressPackageStartupMessages(library(pckg, character.only = T, warn.conflicts = F, verbose = F, quietly = T)))
}

# function to detach packages from current session
detach_package <- function(pckg) {
  suppressWarnings(detach(paste0("package:", pckg), character.only = T))
}

load_package("jsonlite")
load_package("RColorBrewer")
load_package("stringr")

current_script <- "run_analysis.R"
main_loc       <- "main"
sa_miss_loc    <- "sensitivity_analysis/missing_stage"
sa_prst_loc    <- "sensitivity_analysis/prostate_cancer"
sa_rtd_loc    <- "sensitivity_analysis/rtd"
manuscript_figs_loc <-  "./manuscript_results_outputs"
configs        <- fromJSON("./scripts/common_scripts/00_configs.json")

