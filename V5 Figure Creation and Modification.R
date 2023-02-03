#### V5 Data Vis with Vegan####

library(vegan)
library(tidyverse)

####1- Import Data####

#### All 4 datasets together

samples_raw <- readRDS("V4 Files/sample_dataframe.rds")
OTU_raw <- readRDS("V4 Files/OTU_matrix.rds")
traits_raw <- readRDS("V4 Files/trait_dataframe.rds")
traits_as_taxa_raw <- readRDS("V4 Files/trait_taxa_table.rds")

#### Using this tutorial: https://www.rpubs.com/RGrieger/545184

OTU <- as.data.frame(OTU_raw) %>% rownames_to_column("tree")
samples <- samples_raw %>% rownames_to_column("tree") %>% select(tree, enviro, site)

OTU_enviro <- OTU %>% left_join(., samples)

traits <- as.data.frame(traits_raw) %>% rownames_to_column("tree")

traits_enviro <- traits %>% left_join(., samples)

####2- Ordinate####

OTU_ord <- capscale(OTU_raw~1, distance="bray")
OTU_ord_fit <- envfit(OTU_ord ~ enviro, data=OTU_enviro, perm=999)
OTU_ord_fit
plot(OTU_ord_fit)

scores(OTU_ord,display="sites", tidy=TRUE) #or display=�species�
biplot(OTU_ord)
plot(OTU_ord)


traits_ord <- capscale(traits_raw~1, distance="bray")

ordiplot (traits_ord, display = 'sp', type = 'n')
orditorp (traits_ord, display = 'sp')
# clearly short and rhizo separqate, also some between long, short, and rhizo. noerhizo not on this graph, "unknown" very far (probably skewing data)




scores(traits_ord,display="species", tidy=TRUE) #or display=�species�
biplot(traits_ord)

traits_enviro_4_ordiplot <- traits_enviro %>% 
  mutate(enviro_code = if_else(enviro == "arid", 1, 2),
         site_code=rep(0),
         site_code=if_else(site=="mendocino", 1, site_code),
         site_code=if_else(site=="tejon", 2, site_code),
         site_code=if_else(site=="cedar_creek", 3, site_code),
         site_code=if_else(site=="montpellier", 4,site_code))

# plot(traits_ord)
ordiplot(traits_ord, type = "n", arrows=TRUE)
# ordilabel (traits_ord, display = 'sp', cex = 0.25)
points (traits_ord, col = traits_enviro_4_ordiplot$site_code, pch=traits_enviro_4_ordiplot$enviro_code)


legend.vector <- c("mendocino", "tejon", "cedar_creek", "montpellier")

plot(traits_ord, type="n", main="Traits", ylim=c(-5,5)) |>
  points("sites", pch=traits_enviro_4_ordiplot$enviro_code, col=traits_enviro_4_ordiplot$site_code) |>
  text("species", arrows = TRUE, length=0.05, col="blue")

# Making a Legend
#https://stackoverflow.com/questions/38753832/need-help-making-a-legend-for-a-plot-using-the-vegan-package-in-r

pchv <- 1:2
colv <- 1:4

ordiplot(traits_ord, display = 'sites', type = 'n', cex=.75, main="Traits as Species")
with(traits_enviro_4_ordiplot, points(traits_ord, col=colv[traits_enviro_4_ordiplot$site_code], pch=pchv[traits_enviro_4_ordiplot$enviro_code]))
with(traits_enviro_4_ordiplot, legend(x="bottomleft", legend=levels(traits_enviro_4_ordiplot$site_code), col=colv, pch=pchv))



traits_ord_fit <- envfit(traits_ord ~ enviro, data=traits_enviro, perm=999)
traits_ord_fit
plot(traits_ord_fit)

# https://www.davidzeleny.net/anadat-r/doku.php/en:ordiagrams_examples

#to get the scores of the vector arrows:
scores(traits_ord_fit, "vectors")
plot(traits_ord_fit)

####
data(varespec)
vare.pca <- prcomp(varespec)
scores(vare.pca, choices=c(1,2))

