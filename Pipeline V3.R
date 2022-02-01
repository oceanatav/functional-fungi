###Pipeline V3###

library(tidyverse)

#library(ggplot)
#### 1- Functions ####

####
#a.Clean up and consolidate similar spp names in ref_ID####


ref_ID<-readRDS(file = "master_ID_df.rds")

# here is Laura's building a table # 01/18/2022
clean_ref <- function(reference_list_OTU_ID) {
  
  ref_ID_noNAs = reference_list_OTU_ID %>% filter(ID != "NA")

output = tibble()

for(i in 1:nrow(ref_ID_noNAs)) {
  approx_match <- agrep(ref_ID_noNAs$ID[i], x = ref_ID_noNAs$ID, 
                        max =(0.001*nchar(ref_ID_noNAs$ID[i])), ignore.case = TRUE)
  print(as.vector(ref_ID_noNAs[approx_match[1],]))
  output = rbind(as.vector(ref_ID_noNAs[approx_match[1],]), output)
  
}

ref_ID_new <<- output %>% distinct()

}


clean_ref(ref_ID)

####
#b.Return matching species names####
####

#match_names will compare species binomials assigned in your data with an expanding reference list (ref_ID). You input a corresponding OTU and binomial list ("OTU" and "ID") and it will spit out any overlapping OTU with our existing reference list. Only reports exact matches. For instance, Amanita muscaria clone. J459 will not match with Amanita muscaria. match_names also binds the rows of your input data to the masterlist, so is intended to be used concurrently. For the moment renaming the overlapping OTU's has to happen manually.


#https://stackoverflow.com/questions/48834536/r-saving-the-values-from-a-for-loop-in-a-vector-or-list/48834776 

match_names <- function(OTU, ID, ref_ID, site){
  input_df = data.frame(OTU,ID)
  input_df = unique(input_df)
  input_df$ref_name = rep(0)
  for(i in 1:nrow(input_df)) {
    for(j in 1:nrow(ref_ID)) {
      if(input_df$ID[i] == ref_ID$ID[j]) {
        input_df$ref_name[i] = ref_ID$OTU[j]
        print(paste("Replace ID:", input_df$ID[i], "with", ref_ID$OTU[j]))
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

####c.Input Tejon and Richard datasets to ref_ID (all others are already included) ####

#Tejon
tejon_spp<-readRDS("NMDS_4_input/tejon_spp_list.rds") %>% unite(ID, Genus, Species, sep=" ") %>% rename(OTU = OUT_MatchwMendo)
match_names(tejon_spp$OTU, tejon_spp$ID, ref_ID_new, "tejon")
tejon_match_names_output <- output_df
clean_ref(ref_ID_newer)


#Richard
richard_raw<- read_delim("Richard et al OTU table.csv", ";", escape_double = FALSE, trim_ws = TRUE)
richard_spp<-richard_raw %>% select(Taxon) %>%  
transmute(OTU = str_replace(Taxon, " ", "_"),
          ID = Taxon)

match_names(richard_spp$OTU, richard_spp$ID, ref_ID_new, "richard")
richard_match_names_output <- output_df
clean_ref(ref_ID_newer)

####d. Prepare ref_ID to replace all OTU names in future combined matrix####

newest_ref_ID <- clean_ref(ref_ID_newer)

newest_ref_ID$ref_name = ifelse(is.na(newest_ref_ID$ref_name), paste(newest_ref_ID$OTU),paste(newest_ref_ID$ref_name));newest_ref_ID

#tail(newest_ref_ID)

#### 2- Load in data ####

##Tejon
tejon_sample <- readRDS("NMDS_4_input/tejon_enviro_list.rds")
tejon_OTUs <- readRDS("NMDS_4_input/tejon_OTU_matrix.rds")
tejon_traits <- readRDS("NMDS_4_input/tejon_trait_matrix.rds")

##Erlandson
erl_sample <- readRDS("NMDS_4_input/erlandson_enviro_list.rds")
erl_OTUs <- readRDS("NMDS_4_input/erlandson_OTU_matrix.rds")
erl_traits <- readRDS("NMDS_4_input/erlandson_trait_matrix.rds")

##Mendocino
mendo_sample <- readRDS("NMDS_4_input/mendo_enviro_list.rds")
mendo_OTUs <- readRDS("NMDS_4_input/mendo_OTU_matrix.rds")
mendo_traits <- readRDS("NMDS_4_input/mendo_trait_matrix.rds")

##Richard
rich_sample <- readRDS("NMDS_4_input/richard_enviro_list.rds")
rich_OTUs <- readRDS("NMDS_4_input/richard_OTU_matrix.rds")
rich_traits <- readRDS("NMDS_4_input/richard_trait_matrix.rds")


####3- Combine Data #### 

####a. Update OTU Names with New Consolidated OTUs ####

#This has to happen in the datasets individually to reduce issues with too many OTUS in the combined matrix

update_OTUS <- function(input_dataset, newest_reference_list) {
  for(i in 1:ncol(input_dataset)){
    for(j in 1:nrow(newest_reference_list)){
      if()
    }
  }
}

####b. Combine Matrices ####
#Combo OTU

##Make columns and rows
###Columns

OTU_list <- unique(c(colnames(mendo_OTUs), colnames(tejon_OTUs), colnames(erl_OTUs),colnames(rich_OTUs)))
###Rows
tree_names <- c(rownames(mendo_OTUs), rownames(tejon_OTUs), rownames(erl_OTUs), rownames(rich_OTUs))
#tree_list <- 1:(dim(tejon_OTU_matrix)[1]+dim(mendo_OTU_matrix)[1] + dim(erlandson_OTU_matrix)[1])

##Make empty matrix the with our columns and rows

combo_OTU_matrix <- matrix(rep(0,length(tree_names)*length(OTU_list)),nrow=length(tree_names),ncol=length(OTU_list))

colnames(combo_OTU_matrix) <- OTU_list
rownames(combo_OTU_matrix) <- tree_names
combo_OTU_matrix <- combo_OTU_matrix %>% replace(is.na(.), 0)

## Now to fill in the values

for(i in 1:dim(combo_OTU_matrix)[2]){
  OTU <- colnames(combo_OTU_matrix)[i]
  
  # Find OTU data for tejon
  colnum <- 0
  for(j in 1:ncol(tejon_OTUs)){
    if(colnames(tejon_OTUs)[j]==OTU){
      colnum <- j 	}	}	
  if(colnum==0){ 	tejon.abun <- rep(0,nrow(tejon_OTUs)) 	} 
  if(colnum!=0){ 	tejon.abun <- tejon_OTUs[,colnum] }
  
  # Find OTU data for mendocino	
  colnum <- 0
  for(j in 1:dim(mendo_OTUs)[2]){
    if(colnames(mendo_OTUs)[j]==OTU){
      colnum <- j  } 	}	
  if(colnum==0){	mendo.abun <- rep(0,nrow(mendo_OTUs)) 	} 
  if(colnum!=0){	mendo.abun <- mendo_OTUs[,colnum]	}
  
  #Find OTU data for erlandson (Oceana)
  
  colnum <- 0
  for(j in 1:ncol(erl_OTUs)) {
    if(colnames(erl_OTUs)[j]==OTU){
      colnum <- j }	}
  if(colnum == 0){ erlandson.abun <- rep(0, nrow(erl_OTUs))}
  if(colnum != 0){erlandson.abun <- erl_OTUs[,colnum]}
  
  #Find OTU data for richard (Oceana)
  
  colnum <- 0
  for(j in 1:ncol(rich_OTUs)) {
    if(colnames(rich_OTUs)[j]==OTU){
      colnum <- j }	}
  if(colnum == 0){ richard.abun <- rep(0, nrow(rich_OTUs))}
  if(colnum != 0){richard.abun <- rich_OTUs[,colnum]}
  
  combo_OTU_matrix[,i] <- c(tejon.abun, mendo.abun, erlandson.abun, richard.abun)
  
}

# In the above for loop Holly creates an object (colnum). then loops through the columns of a site-specific matrix (eg tejon, mendocino, erlandson's study) and sets the object equal to that value (column name/OTU). if a column name matching that in the site-specific matrix is found, then colnum is not equal to zero and the data from that column (rows) are saved to an object (X.abun). ...
#All OTU matrices are relativized within-study. All observations of frequency are summed across a tree and divided by the total abundance per tree. 

####c. Combine Sample Data####




####4-Phyloseq Object ####

#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("phyloseq")
library(phyloseq)



####5-Ordinate####




####6-Figures####


