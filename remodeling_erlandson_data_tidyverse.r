# Trying out a new pipeline for the Erlandson data.

library(tidyverse)
library(vegan)

# Goal: Get 1) a table of relative OTU abundances by tree,
# 2) a table of relative TRAIT abundances by tree.
# Push both of those through an NMDS
# Color code by site aridity

#setwd("/Users/laurabogar/Documents/2020-2021/Oceana research/TraitAnalysisforOceana")

erlandson_taxa = read_csv("ccectos_wtaxonomy.csv")

# Relative abundance by tree:

# Several root tips don't have tree-level identifiers, remove.
erlandson_taxa = erlandson_taxa %>% filter(!is.na(Ind))

# Group by tree individual?
totalabund = erlandson_taxa %>% group_by(Ind, Ecto_spp) %>% summarize(count = n())

relabund = totalabund %>% group_by(Ind) %>% summarize(totaltips = sum(count))
# No, there are a bunch of trees that only have a couple of root 

# Let's group by tree species and site (Plspp.)
totalabund = erlandson_taxa %>% group_by(Plspp, Ecto_spp) %>% summarize(count = n())
bysiteandtree = totalabund %>% group_by(Plspp) %>% summarize(totaltips = sum(count))

together = left_join(totalabund, bysiteandtree)

relabund = together %>% mutate(relabund = count/totaltips)

justrelabund = relabund %>% select(Tree_site = Plspp, Ecto_spp, relabund)

erlandson_otus = justrelabund %>% spread(Ecto_spp, relabund)
erlandson_otus[is.na(erlandson_otus)] = 0  # I'm sure there's a tidyr way to do this but I don't know it.

erlandson_otus_for_NMDS = as.matrix(erlandson_otus)
# Note: You will probably have to remove the Tree_site column before
# merging with other data (Tejon, Mendocino) for metaMDS().

# Now to repeat for the trait data!

# Note: I am not sure how Holly relativized the Tejon and Mendocino data
# so that they were comparable to each other. They weren't relative abundance, whatever
# was happening. I will ask. We should use whatever metric she used while pulling other
# data into the analysis.