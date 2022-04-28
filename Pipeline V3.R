###Pipeline V3###

library(tidyverse)


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
  old_name <- c()
  new_name <- c()
  for(i in 1:nrow(input_df)) {
    for(j in 1:nrow(ref_ID)) {
      if(input_df$ID[i] == ref_ID$ID[j]) {
        input_df$ref_name[i] = ref_ID$OTU[j]
        print(paste("Replace ID:", input_df$ID[i], "with", ref_ID$OTU[j]))
      
        old_name <- c(old_name, input_df$OTU[i])
        new_name <- c(new_name, ref_ID$OTU[j])
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
 # df_replace <<- data.frame(old_name,new_name) %>% mutate(site=rep(site))
  #print(df)
  ref_ID <- bind_rows(ref_ID, input_df) %>% distinct(.) #%>% mutate(study = rep(site))
  # print(unique(input_df$ref_name))
  ref_ID_newer <<- ref_ID
  
  original_OTU<<-old_name
  new_OTU <<-new_name
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
erl_sample <- readRDS("NMDS_4_input/erlandson_enviro_list.rds") %>% rename(tree=Plspp)
erl_OTUs <- readRDS("NMDS_4_input/erlandson_OTU_matrix.rds")
erl_traits <- readRDS("NMDS_4_input/erlandson_trait_matrix.rds")

##Mendocino
mendo_sample <- readRDS("NMDS_4_input/mendo_enviro_list.rds") %>% rename(tree=Tree) %>% mutate(tree=as.character(tree))
mendo_OTUs <- readRDS("NMDS_4_input/mendo_OTU_matrix.rds")
mendo_traits <- readRDS("NMDS_4_input/mendo_trait_matrix.rds") %>% rownames_to_column(var="tree") %>% mutate(tree=as.character(tree)) %>% column_to_rownames(var="tree")

##Richard
rich_sample <- readRDS("NMDS_4_input/richard_enviro_list.rds") %>% rename(tree=Tree)
rich_OTUs <- readRDS("NMDS_4_input/richard_OTU_matrix.rds")
rich_traits <- readRDS("NMDS_4_input/richard_trait_matrix.rds") %>% t()



####3- Combine Data #### 

####a. Update OTU Names with New Consolidated OTUs ####

#This is going to happen in the datasets individually for my peace of mind.

## Update Tejon ##


replace_names <- tejon_match_names_output %>% filter(ref_name != OTU)

old_col <- c(colnames(tejon_OTUs))
legacy_old_col <- old_col
for(j in 1:nrow(replace_names)) {
  for(i in 1:length(old_col)) {
    if(old_col[i] %in% replace_names$OTU[j]){
     print("yep")
      print(replace_names$ref_name[j])
      
  old_col<<- replace(old_col, i, replace_names$ref_name[j])
   
    }
  
  }
}
#out_colnames == legacy_old_col
tej_comp<- data.frame(legacy_old_col, old_col) %>% filter(legacy_old_col != old_col)
tej_comp2 <- tej_comp


# Now we have to update the OTU list

tej_rep_names <- tejon_OTUs %>%  t() %>% as_tibble(., rownames = "OTU")
for(i in 1:nrow(tej_rep_names)){
  for(j in 1:nrow(tej_comp2)){
  if(tej_rep_names$OTU[i] == tej_comp2$legacy_old_col[j]){
    tej_rep_names$OTU[i] = tej_comp2$old_col[j]
    tej_comp2$legacy_old_col[j] = "DONE"
  }
}

}

tej_resum <- tej_rep_names %>% group_by(OTU) %>% summarise(across(everything(), sum)) %>% column_to_rownames(var="OTU")

tej_OTU_up <- tej_resum %>% t()

## Repeat for Richard ##
## 

replace_names_rich <- richard_match_names_output %>% filter(ref_name != OTU)

old_col <- c(colnames(rich_OTUs))
legacy_old_col <- old_col
for(j in 1:nrow(replace_names_rich)) {
  for(i in 1:length(old_col)) {
    if(old_col[i] %in% replace_names_rich$OTU[j]){
      print("yep")
      print(replace_names_rich$ref_name[j])
      
      old_col<<- replace(old_col, i, replace_names_rich$ref_name[j])
      
    }
    
  }
}
#out_colnames == legacy_old_col
rich_comp<- data.frame(legacy_old_col, old_col) %>% filter(legacy_old_col != old_col)
rich_comp2 <- rich_comp


# Now we have to update the OTU list

rich_rep_names <- rich_OTUs %>%  t() %>% as_tibble(., rownames = "OTU")
for(i in 1:nrow(rich_rep_names)){
  for(j in 1:nrow(rich_comp2)){
    if(rich_rep_names$OTU[i] == rich_comp2$legacy_old_col[j]){
      rich_rep_names$OTU[i] = rich_comp2$old_col[j]
      rich_comp2$legacy_old_col[j] = "DONE"
    }
  }
  
}

rich_resum <- rich_rep_names %>% group_by(OTU) %>% summarise(across(everything(), sum)) %>% column_to_rownames(var="OTU")

rich_OTU_up <- rich_resum %>% t()
#check dim, no change
#damn, so nothing had to be updated? whaaat a waste of time 


####

####b. Combine Matrices ####
#Combo OTU

##Make columns and rows
###Columns

OTU_list <- unique(c(colnames(mendo_OTUs), colnames(tej_OTU_up), colnames(erl_OTUs),colnames(rich_OTU_up)))
###Rows
tree_names <- c(rownames(mendo_OTUs), rownames(tej_OTU_up), rownames(erl_OTUs), rownames(rich_OTU_up))
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
  for(j in 1:ncol(tej_OTU_up)){
    if(colnames(tej_OTU_up)[j]==OTU){
      colnum <- j 	}	}	
  if(colnum==0){ 	tejon.abun <- rep(0,nrow(tej_OTU_up)) 	} 
  if(colnum!=0){ 	tejon.abun <- tej_OTU_up[,colnum] }
  
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
  for(j in 1:ncol(rich_OTU_up)) {
    if(colnames(rich_OTU_up)[j]==OTU){
      colnum <- j }	}
  if(colnum == 0){ richard.abun <- rep(0, nrow(rich_OTU_up))}
  if(colnum != 0){richard.abun <- rich_OTU_up[,colnum]}
  
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

####a. Create PHYLOSEQ objects for OTU NMDS ####
# NOT using taxa table ftm
#https://joey711.github.io/phyloseq/import-data.html

##Sample table


com_sam_df <- bind_rows(mendo_sample, tejon_sample,  erl_sample, rich_sample) %>% replace(is.na(.), 0) %>% column_to_rownames(var="tree")

com_sam_df %>% group_by(enviro) %>% summarise(sum = n()) 
com_sam <- sample_data(com_sam_df)

##OTU table

com_otu <- otu_table(combo_OTU_matrix, taxa_are_rows = F) 

## Phyloseq object

com_ps <- phyloseq(com_otu,com_sam)


####b.Phyloseq objects for trait NMDS ####

com_trait_df <-  bind_rows(mendo_traits, tejon_traits,  erl_traits, data.frame(rich_traits)) %>% replace(is.na(.), 0)# %>% column_to_rownames(var="tree")


com_trait <- otu_table(com_trait_df, taxa_are_rows = F)

dat <- c("long_rhizo", "short_norhizo" ,  "unknown_unknown", "medium_rhizo","medium_unknown" ,"unknown_rhizo" ,  "short_rhizo"  ,   "short_unknown" ,  "unknown_norhizo")
rnames <- c("long_rhizo", "short_norhizo" ,  "unknown_unknown", "medium_rhizo","medium_unknown" ,"unknown_rhizo" ,  "short_rhizo"  ,   "short_unknown" ,  "unknown_norhizo")
cnames <- "Species"

com_trait_tax_df <- matrix(dat, nrow =9, ncol =1, dimnames=list(rnames,cnames))
com_trait_tax <- tax_table(com_trait_tax_df)

com_trait_sam <- com_trait_df %>% rownames_to_column(., var = "tree") %>% right_join(com_trait_tax, .) %>% column_to_rownames(., var = "tree")


com_trait_ps <- phyloseq(com_sam, com_trait_tax, com_trait)

####5-Ordinate####

####a.Taxonomy####


#Taxonomy distribution
com_ps_ord <- ordinate(com_ps, method = "PCoA", distance = "bray", trymax=50)
p_com = plot_ordination(com_ps, com_ps_ord, color="site")
print(p_com)

####b.Trait####

com_ps_ord_trait <- ordinate(com_trait_ps, method = "PCoA", distance = "bray", trymax=50)


####6-Figures####

###a. Taxa figure####

p_com = plot_ordination(com_ps, com_ps_ord, color="site")
print(p_com)

####b. Trait figure####

p_com_t = plot_ordination(com_trait_ps, com_ps_ord_trait, color="site", shape = "site", label = "enviro")
print(p_com_t)


plot_bar(com_trait_ps, fill = "Species")

#   interesting- weird spike CS and EA samples (richard)
#   otherwise all sum to Abund = 1
#   Maybe we should NOT relative abundance?
#   OR we have to fix weird relabund spike


## OTU heatmap

# plot_heatmap(com_ps, method = "PCoA", distance = "bray")
#   weird
#  " Transformation introduced infinite values in discrete y-axis "

## Trait heatmap

plot_heatmap(com_trait_ps, method = "PCoA", distance = "bray")

#   Warning message:
#   Transformation introduced infinite values in discrete y-axis 
#   
#   Unknowns clearly dominat
#   also short_norhizo


## Potential way to view by trait/extract eigenvalues
## Maybe through subsetting our data
## e.g.https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#read-the-data-and-create-phyloseq-objects 7.2
sample_variables(com_trait_ps)
#   "Pygmy"  "enviro" "site"   "mean" 
rank_names(com_trait_ps)
#   "Species"


short <- subset_taxa(com_trait_ps, Species == "short_norhizo" |Species ==  "short_rhizo" | Species == "short_unknown")
short

medium <- subset_taxa(com_trait_ps, Species == "medium_rhizo" |Species ==  "medium_unknown")

long <- subset_taxa(com_trait_ps, Species == "long_rhizo")

rhizo <- subset_taxa(com_trait_ps, Species == "long_rhizo" |Species ==  "medium_rhizo" | Species == "unknown_rhizo" | Species == "short_rhizo")

norhizo <- subset_taxa(com_trait_ps, Species == "short_norhizo" |Species ==  "unknown_norhizo")

unknown <- subset_taxa(com_trait_ps, Species == "unknown_unknown" |Species ==  "medium_unknown" | Species == "short_unknown")

short_ord <- ordinate(short, "PCoA", "bray")

####7- Export data ####

saveRDS(com_sam_df, "V3 Files/sample_dataframe.rds")
saveRDS(com_trait_df, "V3 Files/trait_dataframe.rds")
saveRDS(combo_OTU_matrix, "V3 Files/OTU_matrix.rds")

saveRDS(com_trait_tax_df, "V3 Files/trait_taxa_table.rds")

