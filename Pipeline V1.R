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

match_names <- function()
erlandson_names <- erlandson_metadata %>% select (Tree_site, Ecto_ID, Ecto_spp, Tip) %>% mutate(ref_name = numeric(nrow(.)))


for(i in 1:nrow(erlandson_names)) {
  for(j in 1:nrow(ref_ID)) {
    if(erlandson_names$Ecto_ID[i] == ref_ID$ID[j]) {
      erlandson_names$ref_name[i] = ref_ID$Species[j] 
      if(erlandson_names$ref_name[i] == 0) {
        erlandson_names$ref_name[i] = erlandson_names$Ecto_spp[i]
      }
      
    }
    
  }
}

unique(erlandson_names$ref_name)