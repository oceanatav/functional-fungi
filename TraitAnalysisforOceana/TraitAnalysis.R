#path <- "~/Documents/Fungi Textbooks/TraitAnalysisforOceana"
#setwd(path)
library(vegan)
################################
## species matrix function #####
################################

#this code takes an input dataframe and makes it into
#a matrix that can be used for community similarity 

speciesMatrix<-function(data,species,abundance,sample.ID){
	
	speciesNames<-unique(data[[species]])
	sampleNames<-unique(data[[sample.ID]])
	
	speciesSummary<-aggregate(data[[abundance]],by=list(data[[species]]),FUN=sum,na.rm=TRUE)
	
	siteSummary<-aggregate(data[[abundance]], by=list(data[[sample.ID]]), FUN=sum, na.rm=TRUE)

	speciesByPlot<-aggregate(data[[abundance]], by=list(data[[sample.ID]],data[[species]]), FUN=sum, na.rm=TRUE)

	names(speciesByPlot)=c("sample","species","abundance")

	sppMat<-matrix(nrow=length(sampleNames),ncol=length(speciesNames))
	sppMat[is.na(sppMat)]=0	

	rownames(sppMat)<-sampleNames
	colnames(sppMat)<-speciesNames

	colNumbers<-match(speciesByPlot$species,speciesNames)
	rowNumbers<-match(speciesByPlot$sample,sampleNames)

for(i in 1:length(colNumbers)) {
	
	sppMat[rowNumbers[i],colNumbers[i]]=speciesByPlot$abundance[i]
	
	}

return(sppMat)


}

##END FUNCTION


######################################
## Summary Table w/ Two Treatments ###
######################################

summary.table2<-function(dataset,response,treatments){	
mean.table<-aggregate(dataset[[response]],by=list(dataset[[treatments[1]]],dataset[[treatments[2]]]),mean,na.rm=TRUE)

std.table<-aggregate(dataset[[response]],by=list(dataset[[treatments[1]]],dataset[[treatments[2]]]),sd,na.rm=TRUE)

n.table<-aggregate(dataset[[response]],by=list(dataset[[treatments[1]]],dataset[[treatments[2]]]),length)

se.vector<-std.table$x/sqrt(n.table$x)

final.table<-cbind(mean.table,n.table$x,std.table$x,se.vector)
names(final.table)<-c("Group1","Group2","Mean","n","SD","SE")

return(final.table)

}


##END FUNCTION

################################
## Summar Table w/ 1 Treatment #
################################

summary.table1<-function(dataset,response,treatments){	
mean.table<-aggregate(dataset[[response]],by=list(dataset[[treatments]]),mean,na.rm=TRUE)

std.table<-aggregate(dataset[[response]],by=list(dataset[[treatments]]),sd,na.rm=TRUE)

n.table<-aggregate(dataset[[response]],by=list(dataset[[treatments]]),length)

se.vector<-std.table$x/sqrt(n.table$x)

final.table<-cbind(mean.table,n.table$x,std.table$x,se.vector)

names(final.table)<-c("Group1","Mean","n","SD","SE")

return(final.table)

}

##END FUNCTION

#################################
##Geodetic Distance from Matrix##
#################################
#Here is my program to compute geodetic inter-site distance matrix. Input requires *vectors* in degrees. 
#From http://www.biostat.umn.edu/~sudiptob/Software/distonearth.R

geodetic.distance.matrix <- function (long,lat) {		

NSITES <- length(lat)
R <- 6371

latitude <- lat	
longitude <- long	
latlong <- cbind(latitude, longitude)*pi/180	

d <- matrix(nrow=NSITES, ncol=NSITES)	

for(i in 1:(NSITES-1)) {	
d[i,i] <- 1.0		
	for (j in (i+1):NSITES) {					
	d[i,j] <- sin(latlong[i,1]) * sin(latlong[j,1]) + 			cos(latlong[i,1]) * cos(latlong[j,1])	*				cos(abs(latlong[i,2] - latlong[j,2]))
	d[j,i] <- d[i,j]
	}	
}
d[NSITES, NSITES] <- 1.0
d <- R*acos(d)
d

}
#####END FUNCTION



#####################################
##Extract p-value from Linear Model##
#####################################

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}
#####END FUNCTION



TejonData <- function(){}
setwd("~/GoogleSync/Proposals/DOE/BER_2019/TraitData")
tejon.trees <- read.csv('fungal-abundance.csv')
tejon.ecto <- read.csv('ecto_OTUs.csv')
tejon.ecto.traits <- read.csv('ecto_OTU_traits.csv')
dim(tejon.trees)
dim(tejon.ecto)
head(tejon.ecto$OTU_ID)
head(colnames(tejon.trees))

tejon.trees <- tejon.trees[tejon.trees$site!='intermediate',]
tejon.trees.ecto <- as.data.frame(cbind(tejon.trees$tree,tejon.trees$site))

# Subset out just the Ecto data
for(i in 1:dim(tejon.ecto)[1]){
	OTU <- tejon.ecto$OTU_ID[i]
	
	for(j in 1:dim(tejon.trees)[2]){
		if(colnames(tejon.trees)[j]==OTU){
			colnum <- j
		}
	}
	
	tejon.trees.ecto <- cbind(tejon.trees.ecto,tejon.trees[,colnum])
}

tejon.trees.ecto <- cbind(rep('tejon',dim(tejon.trees.ecto)[1]),tejon.trees.ecto)
tail(tejon.trees.ecto)
colnames(tejon.trees.ecto) <- c('site','tree','enviro',as.character(tejon.ecto$OUT_MatchwMendo))
write.csv(tejon.trees.ecto,'tejon_trees_ecto.csv')

# Make a table that groups by traits, not OTUs
head(tejon.ecto.traits)
levels(tejon.ecto.traits$trait)

tejon.tree.traits <- tejon.trees.ecto[,1:3]
tejon.tree.traits$short <- 0
tejon.tree.traits$long <- 0
tejon.tree.traits$medium <- 0
tejon.tree.traits$norhizo <- 0
tejon.tree.traits$rhizo <- 0

short <- tejon.ecto.traits[tejon.ecto.traits$type=='short',]$OTU_ID
medium <- tejon.ecto.traits[tejon.ecto.traits$type=='medium',]$OTU_ID
long <- tejon.ecto.traits[tejon.ecto.traits$type=='long distance',]$OTU_ID
rhizo <- tejon.ecto.traits[tejon.ecto.traits$type=='r_true',]$OTU_ID
norhizo <- tejon.ecto.traits[tejon.ecto.traits$type=='r_false',]$OTU_ID

for(i in 4:dim(tejon.trees.ecto)[2]){
	# for each fungus
	# determine whether it is a short-dist forager
	hold <- rep(0, dim(tejon.trees.ecto)[1])
	if(colnames(tejon.trees.ecto)[i] %in% short){ hold <- tejon.trees.ecto[,i] }
	tejon.tree.traits$short <- tejon.tree.traits$short + hold
		
	# determine whether it is a medium-dist forager
	hold <- rep(0, dim(tejon.trees.ecto)[1])
	if(colnames(tejon.trees.ecto)[i] %in% medium){ hold <- tejon.trees.ecto[,i] }
	tejon.tree.traits$medium <- tejon.tree.traits$medium + hold
		
	# determine whether it is a long-dist forager
	hold <- rep(0, dim(tejon.trees.ecto)[1])
	if(colnames(tejon.trees.ecto)[i] %in% long){ hold <- tejon.trees.ecto[,i] }
	tejon.tree.traits$long <- tejon.tree.traits$long + hold
	
	# determine whether it makes rhizomorphs
	hold <- rep(0, dim(tejon.trees.ecto)[1])
	if(colnames(tejon.trees.ecto)[i] %in% rhizo){ hold <- tejon.trees.ecto[,i] }
	tejon.tree.traits$rhizo <- tejon.tree.traits$rhizo + hold
		
	# determine whether it doesn't make rhizomorphs
	hold <- rep(0, dim(tejon.trees.ecto)[1])
	if(colnames(tejon.trees.ecto)[i] %in% norhizo){ hold <- tejon.trees.ecto[,i] }
	tejon.tree.traits$norhizo <- tejon.tree.traits$norhizo + hold
	
}



MendocinoData <- function(){}

#read in files
pygmy.raw<-read.csv("Curated_Assembly_3July2012.csv")
head(pygmy.raw) #looks at first 6 lines
names(pygmy.raw) #looks at variable names
pygmy.clean<-read.csv("pygmy-cleandata-wsoils.csv")
is.na(tail(pygmy.clean))

pygmy.crop<-pygmy.clean[order(pygmy.clean$SpeciesCode),]
sort(unique(pygmy.clean$SpeciesCode))
pygmy.crop[(470:490),(1:5)]
pygmy.crop<-pygmy.crop[(1:475),]
dim(pygmy.crop)



#############Condensing to core by species
names(pygmy.clean)
pygmy.core<-aggregate(pygmy.clean$No.Tips,list(pygmy.clean$Sample.ID,pygmy.clean$SpeciesCode,pygmy.clean$Tree,pygmy.clean$Core,pygmy.clean$Terrace,pygmy.clean$Pygmy,pygmy.clean$Trait.ForagingType,pygmy.clean$Trait.Hydrophobicity,pygmy.clean$Trait.Rhizomorphs),sum)
colnames(pygmy.core)<-c("Sample.ID","SpeciesCode","Tree","Core","Terrace","Pygmy","Trait.ForagingType","Trait.Hydrophobicity","Trait.Rhizomorphs","Tips")
head(pygmy.core)
pygmy.core<-cbind(pygmy.core,c(rep(1,dim(pygmy.core)[1])))
colnames(pygmy.core)<-c("Sample.ID","SpeciesCode","Tree","Core","Terrace","Pygmy","Trait.ForagingType","Trait.Hydrophobicity","Trait.Rhizomorphs","Tips","Spp")


#############Condensing to tree by species
names(pygmy.clean)
pygmy.tree<-aggregate(pygmy.clean$No.Tips,list(pygmy.clean$SpeciesCode,pygmy.clean$Tree,pygmy.clean$Terrace,pygmy.clean$Pygmy,pygmy.clean$Trait.ForagingType,pygmy.clean$Trait.Hydrophobicity,pygmy.clean$Trait.Rhizomorphs),sum)
colnames(pygmy.core)<-c("SpeciesCode","Tree","Terrace","Pygmy","Trait.ForagingType","Trait.Hydrophobicity","Trait.Rhizomorphs","Tips")
head(pygmy.tree)
pygmy.tree<-cbind(pygmy.tree,c(rep(1,dim(pygmy.tree)[1])))
colnames(pygmy.tree)<-c("SpeciesCode","Tree","Terrace","Pygmy","Trait.ForagingType","Trait.Hydrophobicity","Trait.Rhizomorphs","Tips","Spp")
tail(pygmy.tree)

tip_count_trait=tapply(tree.trait.matrix.cond[,1],terrace,FUN=sum)+tapply(tree.trait.matrix.cond[,2],terrace,FUN=sum)+tapply(tree.trait.matrix.cond[,3],terrace,FUN=sum)
tip_count = tapply(rowSums(tree.matrix),terrace,FUN=sum)
tip_count_trait/tip_count


tree.matrix<-speciesMatrix(pygmy.clean,"SpeciesCode","No.Tips","Tree")
dim(tree.matrix)  #Check to make sure there're 29 rows; if not...
#tree.matrix=tree.matrix[1:29,]

tree.trait.matrix<-cbind(speciesMatrix(pygmy.clean,"Trait.ForagingType","No.Tips","Tree"),speciesMatrix(pygmy.clean,"Trait.Hydrophobicity","No.Tips","Tree"),speciesMatrix(pygmy.clean,"Trait.Rhizomorphs","No.Tips","Tree"))
head(tree.trait.matrix)
tree.trait.matrix<-tree.trait.matrix[,c(1:2,4:9,12:13)]

#Condensed trait matrix to make all medium f-types into one column
tree.trait.matrix.cond<-cbind(tree.trait.matrix[,1:2],rowSums(tree.trait.matrix[,3:5]),tree.trait.matrix[,6:10])

#Proportion of tips w/ traits
trait_tip_count = tapply(tree.trait.matrix.cond[,1],terrace,FUN=sum)+tapply(tree.trait.matrix.cond[,2],terrace,FUN=sum)+tapply(tree.trait.matrix.cond[,3],terrace,FUN=sum)
tip_count = tapply(rowSums(tree.matrix),terrace,FUN=sum)
trait_tip_count/tip_count



ConsolidateMatrices <- function(){}
OTU_list <- unique(c(colnames(tree.matrix),colnames(tejon.trees.ecto)[4:120]))
length(OTU_list)

tree_list <- 1:(dim(tejon.trees.ecto)[1]+dim(tree.matrix)[1])

combo.matrix <- matrix(rep(0,length(tree_list)*length(OTU_list)),nrow=length(tree_list),ncol=length(OTU_list))

dim(combo.matrix)

colnames(combo.matrix) <- OTU_list
rownames(combo.matrix) <- tree_list

for(i in 1:dim(combo.matrix)[2]){
	OTU <- colnames(combo.matrix)[i]
	
	# Find OTU data for tejon
	colnum <- 0
	for(j in 1:dim(tejon.trees.ecto)[2]){
		if(colnames(tejon.trees.ecto)[j]==OTU){
			colnum <- j 	}	}	
	if(colnum==0){ 	tejon.abun <- rep(0,dim(tejon.trees.ecto)[1]) 	} 
	if(colnum!=0){ 	tejon.abun <- tejon.trees.ecto[,colnum] }

	# Find OTU data for mendocino	
	colnum <- 0
	for(j in 1:dim(tree.matrix)[2]){
		if(colnames(tree.matrix)[j]==OTU){
			colnum <- j  } 	}	
	if(colnum==0){	mendo.abun <- rep(0,dim(tree.matrix)[1]) 	} 
	if(colnum!=0){	mendo.abun <- tree.matrix[,colnum]	}
	
	combo.matrix[,i] <- c(tejon.abun,mendo.abun)
	
}


combo.trait.matrix <- matrix(rep(0,length(tree_list)*5),nrow=length(tree_list),ncol=5)
colnames(combo.trait.matrix) <- c('short','med','long','rhizo','norhizo')
rownames(combo.trait.matrix) <- tree_list

combo.trait.matrix <- as.data.frame(combo.trait.matrix)

head(tejon.tree.traits)
head(tree.trait.matrix.cond)
colnames(tree.trait.matrix.cond) <- c('short','long','med','contact','nohydro','hydro','norhizo','rhizo')
mendo.tree.traits <- as.data.frame(tree.trait.matrix.cond)
combo.trait.matrix$short <- c(tejon.tree.traits$short, mendo.tree.traits$short)
combo.trait.matrix$med <- c(tejon.tree.traits$medium,mendo.tree.traits$med)
combo.trait.matrix$long <- c(tejon.tree.traits$long,mendo.tree.traits$long)
combo.trait.matrix$rhizo <- c(tejon.tree.traits$rhizo,mendo.tree.traits$rhizo)
combo.trait.matrix$norhizo <- c(tejon.tree.traits$norhizo,mendo.tree.traits$norhizo)



NMDSPlots <- function(){}
require(vegan)
siteIDs <- c(rep(1,dim(tejon.trees.ecto)[1]),rep(2,dim(tree.matrix)[1]))

pygmy<-pygmy.clean$Pygmy[match(row.names(tree.matrix),pygmy.clean$Tree)]

aridIDs <- c(tejon.trees.ecto$enviro,pygmy)
for(i in 1:length(aridIDs)){if(aridIDs[i]==3){aridIDs[i]<-0}}

speciesNMDS <- metaMDS(as.matrix(combo.matrix),k=2,trymax=50)
traitNMDS <- metaMDS(as.matrix(combo.trait.matrix))

quartz(height=4.2,width=8.7)
par(mar=c(4,4,1,1),mfrow=c(1,2))
plot(speciesNMDS$points[,1],speciesNMDS$points[,2],pch=c(21,22)[unclass(siteIDs)],main='OTUs',bg=c('black','white')[as.factor(aridIDs)],xlim=c(-2,2.2),cex=1.5,xlab='NMDS 1',ylab='NMDS 2');#legend(-1.6,2.1,c('Arid, Tejon','Wet, Tejon','Arid, Mendo.','Wet, Mendo.'),pch=c(21,21,22,22),pt.cex=1.5,pt.bg=c('white','black','white','black'))


plot(traitNMDS$points[,1], traitNMDS$points[,2],pch=c(21,22)[unclass(siteIDs)],main='Traits',bg=c('black','white')[as.factor(aridIDs)],ylim=c(-1,1.2),cex=1.5,xlab='NMDS 1',ylab='NMDS 2');legend(-1.3,1.25,c('Arid, Tejon','Wet, Tejon','Arid, Mendo.','Wet, Mendo.'),pch=c(21,21,22,22),pt.cex=1.5,pt.bg=c('white','black','white','black'),bg='white'); #ordiellipse(traitNMDS,aridIDs)

