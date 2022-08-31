# Helper script
# R function to make dataset ready for genind
# 
# Femke Batsleer
#
# April 2021
#
#

#allele_data: data with structure of for each locus 2 columns, and a column with Sample (name of sample)
#pop_data: gives vector or dataframe with the population-names 

library(dplyr)
library(tibble)

data4genind <- function(locus_data){
  #convert sample column into rownames
  locus_df <- locus_data %>% remove_rownames %>% column_to_rownames(var="Sample")
  #concatenate the columns pairwise
  n2_microsats <- ncol(locus_df) #number of microsats
  cols <- split(names(locus_df), sub("_.", "", names(locus_df)))#get correct names for microsats
  
  #alleles data
  data_conc <- locus_df %>% bind_cols(as.data.frame(sapply(names(cols), function(col) {
    do.call(paste, c(locus_df[cols[[col]]], sep = "/"))}))) %>%
    dplyr::select(-1:-n2_microsats/2) #only select the concatenated data
  #names of samples are in the rownames of the data
  
  return(data_conc)
}