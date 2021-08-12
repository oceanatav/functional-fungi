###Pipeline V1###

library(tidyverse)
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("phyloseq")
library(phyloseq)
#library(ggplot)


####
#Function: Return matching species names#
####

ref_ID <- readRDS("master_ID_df.rds")

#match_names will compare species binomials assigned in your data with an expanding reference list (ref_ID). You input a corresponding OTU and binomial list ("OTU" and "ID") and it will spit out any overlapping OTU with our existing reference list. Only reports exact matches. For instance, Amanita muscaria clone. J459 will not match with Amanita muscaria. match_names also binds the rows of your input data to the masterlist, so is intended to be used concurrently. For the moment renaming the overlapping OTU's has to happen manually.

match_names <- function(OTU, ID, ref_ID){
  input_df = data.frame(OTU,ID)
  input_df = unique(input_df)
  input_df$ref_name = rep(0)
  for(i in 1:nrow(input_df)) {
    for(j in 1:nrow(ref_ID)) {
      if(input_df$ID[i] == ref_ID$ID[j]) {
        input_df$ref_name[i] = ref_ID$OTU[j]
        print(paste("Replace ID:", input_df$ID[i], "with", ref_ID$OTU[j]))}}}
  ref_ID <- bind_rows(ref_ID, input_df) %>% distinct(.)
 # print(unique(input_df$ref_name))
  ref_ID <<- ref_ID
}

#match_names(test$Ecto_spp, test$Ecto_ID, ref_ID)
#erlandson_OTU_combined <- erlandson_OTU %>% select(order(colnames(.))) %>% 
 # rename("17" = Cortinarius1,
  #       "3" = Cenococcum_geophilum)

####








