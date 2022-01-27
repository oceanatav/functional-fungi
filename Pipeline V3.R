###Pipeline V3###

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
#Clean up and consolidate similar spp names in ref_ID####
####
# here is Laura's building a table # 01/18/2022
clean_ref <- function(ref_ID) {
  
  ref_ID_noNAs = ref_ID %>% filter(ID != "NA")

output = tibble()

for(i in 1:nrow(ref_ID_noNAs)) {
  approx_match <- agrep(ref_ID_noNAs$ID[i], x = ref_ID_noNAs$ID, 
                        max =(0.001*nchar(ref_ID_noNAs$ID[i])), ignore.case = TRUE)
  print(as.vector(ref_ID_noNAs[approx_match[1],]))
  output = rbind(as.vector(ref_ID_noNAs[approx_match[1],]), output)
  
}

new_ref_ID <<- output %>% distinct()

}
clean_ref(ref_ID)
####
#Return matching species names####
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
  output_df <<- output_df %>% mutate(site = rep(site))
  ref_ID <- bind_rows(ref_ID, input_df) %>% distinct(.) #%>% mutate(study = rep(site))
  # print(unique(input_df$ref_name))
  ref_ID_newer <<- ref_ID
}

#match_names(test$Ecto_spp, test$Ecto_ID, ref_ID)
#erlandson_OTU_combined <- erlandson_OTU %>% select(order(colnames(.))) %>% 
# rename("17" = Cortinarius1,
#       "3" = Cenococcum_geophilum)

####Input Tejon and Richard datasets to ref_ID (all others are already included) ####

#Tejon
tejon_spp<-readRDS("NMDS_4_input/tejon_spp_list.rds") %>% unite(ID, Genus, Species, sep=" ") %>% rename(OTU = OUT_MatchwMendo)
match_names(tejon_spp$OTU, tejon_spp$ID, new_ref_ID, "tejon")
tejon_match_names_output <- output_df


#Richard
richard_raw<- read_delim("Richard et al OTU table.csv", ";", escape_double = FALSE, trim_ws = TRUE)
richard_spp<-richard_raw %>% select(Taxon) %>%  
transmute(OTU = str_replace(Taxon, " ", "_"),
          ID = Taxon)
clean_ref(ref_ID_newer)
match_names(richard_spp$OTU, richard_spp$ID, new_ref_ID, "richard")
richard_match_names_output <- output_df
