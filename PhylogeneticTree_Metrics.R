####Phylo metrics calculation#######
install.packages("picante")
library(picante)
library(vegan)
library(ape)
library(car)
library(ggplot2)
#read in tree

setwd("C:/Users/shann/OneDrive/Documents/Swei lab/R")
phylo<-read.tree("tree_squirrel.newick")

class(phylo)
plot(phylo)
#this is consensus tree

#abundance data
abund<-read.csv("community.csv", header=TRUE, row.names=2)

#make sure community and tree match
matched <- picante::match.phylo.comm(phy = phylo, comm = abund)
#keeping matched community data
abund<-matched$comm

#calculate Faiths PD###
phylo.pd<-pd(abund,phylo)  #SR= species richness
head(phylo.pd)
write.csv(phylo.pd, "~/Swei lab/R/Faith_Pd.csv") #saving

boxplot(phylo.pd$PD~all$Site)
plot(phylo.pd$PD~phylo.pd$SR, xlab="Species richness", ylab = "Faith's PD")

#MPD and MNTD####
#convert phylogeny to a distance matrix for mpd
phy.dist<-cophenetic(phylo)

#calculate ses.mpd (standaridzed mpd)
#use taxa.labels which randomly shuffles tip labels
sesmpd <- ses.mpd(abund, phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, 
                  runs = 999)
head(sesmpd)

write.csv(sesmpd, "~/Swei lab/R/mpd_taxalabels.csv")

#####

###calculate mntd####
#with taxa labels
sesmntd <- ses.mntd(abund, phy.dist, null.model = "taxa.labels", abundance.weighted = FALSE, 
                    runs = 999)

write.csv(sesmntd, "~/Swei lab/R/mntd_taxa.csv")

######
