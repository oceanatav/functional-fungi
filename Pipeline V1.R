###Pipeline V1###

library(tidyverse)
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("phyloseq")
library(phyloseq)
#library(ggplot)

####
#Read in your data
####

    ##Tejon
tejon_sample <- readRDS("NMDS_4_input/tejon_enviro_list.rds")
tejon_OTUs <- readRDS("NMDS_4_input/tejon_OTU_matrix.rds")
tejon_traits <- readRDS("NMDS_4_input/tejon_trait_matrix.rds")

    ##Erlandson
erl_sample <- readRDS("NMDS_4_input/erl_enviro_list.rds")
erl_OTUs <- readRDS("NMDS_4_input/erl_OTU_matrix.rds")
erl_traits <- readRDS("NMDS_4_input/erl_trait_matrix.rds")

    ##Mendocino
mendo_sample <- readRDS("NMDS_4_input/mendo_enviro_list.rds")
mendo_OTUs <- readRDS("NMDS_4_input/mendo_OTU_matrix.rds")
mendo_traits <- readRDS("NMDS_4_input/mendo_trait_matrix.rds")

    ##Richard
rich_sample <- readRDS("NMDS_4_input/rich_enviro_list.rds")
rich_OTUs <- readRDS("NMDS_4_input/rich_OTU_matrix.rds")
rich_traits <- readRDS("NMDS_4_input/rich_trait_matrix.rds")

####
#Function: Clean up and conslidate similar spp names in ref_ID
####
ref_ID_raw <- readRDS("master_ID_df.rds")
ref_ID <- ref_ID_raw

clean_up_ref <- function(ref_ID) {
  for(i in 1:nrow(ref_ID)){
    #ten_percent <- 0.1*nchar(ref_ID$ID[i])
    if(grepl("sp.", ref_ID$ID[i],ignore.case = T) == FALSE) {
    fuzzy_matches <- agrep(ref_ID$ID[i], x = ref_ID$ID, max =(0.001*nchar(ref_ID$ID[i])), ignore.case = TRUE) #returns the indices of fuzzy matches
    
    shortest_name <- min(fuzzy_matches) #I want to pull the shortest name from this and then make that a row in a new, "cleaned up" ref_ID
    
    #I also need to homogenize anything with matching OTU's
    NEW_LIST =rbind(ref_ID[fuzzy_matches,])
    
  #we should exclude anything with "sp." at the end from this analysis huh
    #that way we retain the diffierent OTU's
    }
    }
  new_df <- NEW_LIST
  fuzz_test <<- fuzzy_matches
}

clean_up_ref(ref_ID)
cenge <- agrep(ref_ID$ID[1], x = ref_ID$ID, max = 2, ignore.case = TRUE)
cenge_short <- min(nchar(ref_ID[cenge,]))
tomentella <- agrep(ref_ID$ID[187], x = ref_ID$ID, max = 1, ignore.case = TRUE)
 
r1 <- ref_ID_raw %>% 
ref_ID_for_tejon <- ref_ID_raw %>% filter(., !grepl('OTU_', OTU)) #this is only for tejon data-since I re-annotated species names, removing the previous annotations will help with duplicate ID's (hopefully)

####
#Function: Return matching species names#
####

#match_names will compare species binomials assigned in your data with an expanding reference list (ref_ID). You input a corresponding OTU and binomial list ("OTU" and "ID") and it will spit out any overlapping OTU with our existing reference list. Only reports exact matches. For instance, Amanita muscaria clone. J459 will not match with Amanita muscaria. match_names also binds the rows of your input data to the masterlist, so is intended to be used concurrently. For the moment renaming the overlapping OTU's has to happen manually.

match_names <- function(OTU, ID, ref_ID, site){
  input_df = data.frame(OTU,ID)
  input_df = unique(input_df)
  input_df$ref_name = rep(0)
  for(i in 1:nrow(input_df)) {
    for(j in 1:nrow(ref_ID)) {
      if(input_df$ID[i] == ref_ID$ID[j]) {
        input_df$ref_name[i] = ref_ID$OTU[j]
       # print(paste("Replace ID:", input_df$ID[i], "with", ref_ID$OTU[j]))
      }
      if(input_df$ref_name[i] == 0){
       input_df$ref_name[i] = input_df$OTU[i]
      }
    }}
output_df <- input_df
 # output_df$change <- output_df$OTU != output_df$ref_name
 # output_df$new_name <- rep(0)
  #for(i in 1:nrow(output_df)){
   # if(output_df$ref_name[i] != "0"){
    #  output_df$new_name[i] = output_df$ref_name[i]}
    #else{
    #  output_df$new_name[i] = output_df$OTU
    #}
  #}
  output_df <<- output_df
  #ref_ID <- bind_rows(ref_ID, input_df) %>% distinct(.) %>% mutate(study = rep(site))
 # print(unique(input_df$ref_name))
  #ref_ID <<- ref_ID
}

#match_names(test$Ecto_spp, test$Ecto_ID, ref_ID)
#erlandson_OTU_combined <- erlandson_OTU %>% select(order(colnames(.))) %>% 
 # rename("17" = Cortinarius1,
  #       "3" = Cenococcum_geophilum)

##Run Tejon and Richard datasets (all others are already included)

    #Tejon
tejon_spp<-readRDS("NMDS_4_input/tejon_spp_list.rds") %>% unite(ID, Genus, Species, sep=" ") %>% rename(OTU = OUT_MatchwMendo)
match_names(tejon_spp$OTU, tejon_spp$ID, ref_ID, "tejon")

    #Richard








