
# Author: Zoe Dyson (zoe.dyson@lshtm.ac.uk)
# Title: NARMS_functions.R
# Date: 2024-10-17


# Package management
# List required packages
packages <- c("tidyverse")

# Load packages & install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# combine_sal_narms: Combine Salmonella NARMS AST files
#   Parameters: 
#     1. user_path = path to NARMS files in excel format (one per serovar)
#     2. file_pattern = sting to identify files (unique prefix or suffix)
#        (default = "Salmonella_NARMS")
#
#   Returns: data frame of combined AST data
#
combine_narms <- function(user_path="AST_data", file_pattern="Salmonella_NARMS"){
  
  # list NARMS files
  narms_files <- list.files(user_path)[grep(file_pattern,list.files(user_path))]
  
  # create empty data frame
  combined_narms <- NULL
  
  # combine NARMS files
  for (file in narms_files){
    temp_file <- readxl::read_excel(paste0(user_path,"/",file))
    combined_narms <- bind_rows(combined_narms, temp_file)
  }
  
  # return combined data frame
  return(as.data.frame(combined_narms))
}
