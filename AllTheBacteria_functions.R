
# Author: Zoe Dyson (zoe.dyson@lshtm.ac.uk)
# Title: AllTheBacteria_functions.R
# Date: 17/10/2024


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

# atb_amrfp_filter_by_taxa: Extract specific taxa from the AllTheBacteria
# AMRFinderPlus (AMRFP) data (from https://osf.io/zgexh, accessed 17/10/24)
# (via browser, click three dots, download)
#
#   Parameters: 
#     1. user_taxa = character vector of taxa to be extracted (default = Salmonella)
#     2. atb_amrfp_results = tab delimited file containing complete AMRFP calls
#     3. atb_species_results = tab delimited file containing taxonomic assignments
#
#   Returns: data frame of AMRFP hits for specific taxa
#
atb_amrfp_filter_by_taxa <- function(user_taxa=c("Salmonella"), atb_armfp_results_path="AllTheBacteria/AMRFP_results.tsv.gz", atb_armfp_qc_status_path="AllTheBacteria/ATB_AMRFP_status.tsv.gz", atb_species_results_path="AllTheBacteria/ATB_species_calls.tsv.gz"){
  
  # load files
  atb_species_results <- read_tsv(atb_species_results_path)
  atb_armfp_qc_status <- read_tsv(atb_armfp_qc_status_path)
    
  # create an empty vector for sample accessions
  combined_sample_accessions <- NULL
  
  # get complete list of all relevant sample accessions
  for (taxa in user_taxa){
    
    # Get sample accessions for all taxa
    temp_sample_accessions <- atb_species_results %>%
      filter(grepl(taxa, Species)) %>%
      select(Sample) %>%
      unique() %>%
      pull()
    
    # Add accessions to character vector
    combined_sample_accessions <- c(temp_sample_accessions, combined_sample_accessions)
  }
  
   # nested callback function to read in only those lines containing selected taxa (for chunking)
  f <- function(x, pos) subset(x, Name %in% combined_sample_accessions)
  
  # read lines for selected taxa  
  selected_atb_amrfp <- read_tsv_chunked(file=atb_armfp_results_path, callback=DataFrameCallback$new(f), chunk_size = 10000)
  
  selected_atb_amrfp <- selected_atb_amrfp %>%
    left_join(atb_species_results, by=c("Name"="Sample")) %>% 
    left_join(atb_armfp_qc_status, by=c("Name"="sample"))
  
  # return data frame for selected taxa
  return(selected_atb_amrfp)
}
