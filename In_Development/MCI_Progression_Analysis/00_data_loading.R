################################################################################
# 00_data_loading.R
# Load and define data structure
################################################################################

library(tidyverse)

#' Load Alamar biomarker data
#' 
#' @param filepath Path to filtered_combined_post_QC.csv
#' @return List with data and metadata
load_alamar_data <- function(filepath) {
  
  # Load data
  dat <- read_csv(filepath, show_col_types = FALSE)
  
  # Define biomarker columns (columns 26 onwards based on header inspection)
  biomarker_cols <- names(dat)[26:ncol(dat)]
  
  # Define covariate columns
  covariate_cols <- c("Record_ID", "Group", "Ethnicity", "Race", "sex", 
                      "age_at_subject", "CDX", "Case_Control", "AOO", 
                      "Years_Onset", "APOE.geno", "BMI", "Site", 
                      "Group_Race_Ethnicity")
  
  # Return structured list
  list(
    data = dat,
    biomarker_cols = biomarker_cols,
    covariate_cols = covariate_cols,
    n_biomarkers = length(biomarker_cols),
    n_samples = nrow(dat)
  )
}

#' Get biomarker matrix
#' 
#' @param dat_list Output from load_alamar_data
#' @return Matrix of biomarkers only
get_biomarker_matrix <- function(dat_list) {
  dat_list$data %>%
    select(all_of(dat_list$biomarker_cols)) %>%
    as.matrix()
}

#' Get covariates dataframe
#' 
#' @param dat_list Output from load_alamar_data
#' @return Dataframe of covariates
get_covariates <- function(dat_list) {
  dat_list$data %>%
    select(any_of(dat_list$covariate_cols))
}
