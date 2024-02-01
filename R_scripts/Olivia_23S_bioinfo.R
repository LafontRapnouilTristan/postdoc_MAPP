############################################################################################################################
#### Project: Olivia's project
#### Title:   Bioinfo on 23S data
#### Author:  Vincent Jassey (vincent.jassey@univ-tlse3.fr)
#### Date:    22 January 2024
#### -----------------------------------------------------------------------------------------------------------------------


#packages ----
source('https://bioconductor.org/biocLite.R')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("Biobase")
BiocManager::install("rhdf5")
BiocManager::install("multtest")
BiocManager::install("phyloseq")
BiocManager::install("ggtree")

library(ggplot2); library(plyr); library(reshape2); library(grid); library(gridExtra); library(ape); library(picante)
library(seqinr); library(phytools); library(scales); library(nlme); library(vegan); library(RColorBrewer)
library("BiocManager"); library("Biobase"); library("rhdf5"); library("devtools"); library(phyloseq)
library(ggtree); library(tidyjson); library(RColorBrewer);library(tidyverse); library(viridis); library(randomForest); library(caret); library(cluster)
library(tidyverse); library(cluster);  library(parallelDist); library(boral);  library(corrplot); library(igraph); library(car); library(fossil); library(ROCR); library(doMC)
  #library(divesitree); 
library(caret)
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R" )
library(dplyr); library(tidyr)

#source("nested_donut_plot.R")
#source("test.a.R")

###### DO IT only the first time- pre-treatment of the fastafile ###################

#23S 
fastafile<-read.FASTA("23S_Affiliation_postprocess__sequences_Olivia.fasta")
V<-names(fastafile)
for (i in 1:length(V)) { V[i]<- substr(V[i], 1, nchar(V[i])-13)}
names(fastafile)<-V
write.dna(fastafile, "23S_Olivia_correctnames.fasta", format = "fasta", nbcol = -1, colsep = "")

################### 23S bioinfo #################################################################################################### ----

####
#### Import du biom 23S
####

#Import and rename IDs
Biom23 <-import_biom("23S_Olivia_Biom.biom1", refseqfilename='23S_Olivia_correctnames.fasta')
colnames(tax_table(Biom23)) <-c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Do it only the first time
noms<-as.vector(sample_names(Biom23))
id<- as.vector(substr(noms, 1, nchar(noms)-18)) #remove the last 18 characters to have shorter names
sample_names(Biom23)<-id
write.table(as.data.frame(id), "ids_23s.csv", sep = ";")

#Create a .txt table with metadata
metadata23S <- read.delim("metadata_23S_Oliviaids_23s.csv", sep = ";", dec = ",")
metadata23S <- metadata23S[match(sample_names(Biom23), metadata23S$id), ]
rownames(metadata23S)<-sample_names(Biom23)
sample_data(Biom23)<-metadata23S
sample_sums(Biom23)

# Correct the Biom for Domain and Kingdom Mise en place d'une colonne domaine (Prokaryote) dans la tax_table de 23 prokaryote
a23prokaryote=data.frame(tax_table(Biom23))
a23prokaryote[] <- lapply(a23prokaryote, as.character)
a23prokaryote[a23prokaryote=='k__unclassified']=as.character('k__Eukaryota')
a23prokaryote[a23prokaryote=='k__Bacteria']=as.character('k__Prokaryota')
tax_table(Biom23)=as.matrix(a23prokaryote)
colnames(tax_table(Biom23prokaryote)) <-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table(Biom23) <- tax_table(Biom23)[,2:8]
sum(sample_sums(Biom23)) # 1228909 sequences in total

#remove Extraction tests
Biom23_woTE <- subset_samples(Biom23, sample_names(Biom23) !="TE-MEC") #
Biom23_woTE  <- subset_samples(Biom23_woTE, sample_names(Biom23_woTE) !="TE-MED") # 
Biom23_woTE  <- subset_samples(Biom23_woTE, sample_names(Biom23_woTE) !="TE-OLC") #
Biom23_woTE  <- subset_samples(Biom23_woTE, sample_names(Biom23_woTE) !="TE-OLD") # 
Biom23_woTE  <- subset_samples(Biom23_woTE, sample_names(Biom23_woTE) !="TE-OMD") # 
Biom23_woTE  <- subset_samples(Biom23_woTE, sample_names(Biom23_woTE) !="inconnu") # 

sample_names(Biom23_woTE)

###
### Filtration des plantes (Class=Embryophyceae) - Before removing plantsfiltration : 1 228 909 sequences pour 33 249 OTUs
###

# Enlevement des sequences de plantes (Class=Embryophyceae  Domain = Eukaryota > Kingdom = Archaeplastida > Phylum = Streptophyta > Class = Embryophyceae)
Biom23_F = subset_taxa(Biom23_woTE, Kingdom!= "k__Chloroplastida" & Class!="c__Embryophyceae") # OTUs 14 092
sum(sample_sums(Biom23_F)) #After filtration 340 626 | 71% of sequences removed

###
### Rarefaction curve | cut at 13470
###
rarecurve(t(as.data.frame(otu_table(Biom23_F))), step=50, cex=0.5) #step est l'ecart des points sur l'axe x (samples)

###
### Correlation btwn total number of reads versus number of OTUs
###

nb.otu <- as.vector(colSums(!(otu_table(Biom23_F)==0))) #nbr otu per sample
nb.reads <- as.vector(colSums((otu_table(Biom23_F)))) #nbr reads per sample

plot(nb.otu~nb.reads)
summary(lm(nb.otu~nb.reads)) # strong correlation so we have to rarefy

# Cut OTUs following rarefaction
Biom23_F_cut <- rarefy_even_depth(Biom23_F, rngseed = 12345) # keep all samples
#4981 OTUs removed, 10% of OTUs removed

#taxonomy table 18S tropical RF
a23=as.data.frame(tax_table(Biom23_F_cut)) #Transformation taxonomie en data frame

####
#### Separation des bioms avant rerefaction (a refaire apres)
####

sum(sample_sums(Biom23_F_cut)) # 121 794 sequences 4981 OTUs

#Extract abundance matrix from the phyloseq object
otu_mat <- as(otu_table(Biom23_F_cut), "matrix")
# transpose
if(taxa_are_rows(Biom23_F_cut)){otu_mat <- t(otu_mat)}
# Coerce to data.frame
otu = as.data.frame(otu_mat)

####
#### Quick graphs to check the community composition ----
####

## Relative abundance per phyla (raw data)
tax_table(Biom23)
p <- plot_composition(Biom23_F_cut, taxaRank1 = "Kingdom", taxaSet1 = "k__Prokaryota",
                      taxaRank2 = "Phylum", fill = "Phylum", numberOfTaxa = 15)
p

## Relative abundance per phyla (grouped data)
groupBiom23S <- merge_samples(Biom23_F_cut, group = sample_data(Biom23_F_cut)$Code, fun = sum)
tax_table(groupBiom23S)

p <- plot_composition(groupBiom23S, taxaRank1 = "Kingdom", taxaSet1 = "k__Eukaryota",
                      taxaRank2 = "Phylum", fill = "Phylum", numberOfTaxa = 15)
p


df <- groupData23 %>% group_by(Kingdom, Code) %>% 
  summarise_at(c("Abundance"), sum)

theme_set(new = theme_classic())
ggplot(df, aes(fill=Kingdom, y=Abundance, x=Code)) + 
  geom_bar(position="fill", stat="identity")

################### Figure pie charts #################################################################################################### ----

###
### Nested pie charts
###

### Preprocessing of the data
# remove 'multi-affiliation' taxa and transform into relative abundance data
#groupData <- subset_taxa(groupBiom18photo_peat, Species != "Multi-affiliation")
#groupData_relab <- transform_sample_counts(groupData, function(OTU) OTU/sum(OTU))

# merge samples
#melted <- psmelt(groupData_relab)
#melted <- as.data.frame(melted)
#melted <- melted[,-c(4:8)]

# Summarise data by sum abundance by family
#df <- melted %>% group_by(Family, Order, Class, ecotype) %>% 
#  summarise(Abundance=sum(Abundance))
#df <- as.data.frame(df)

###
### 23S dataset
###

### Preprocessing of the data
#*****************************

# remove 'multi-affiliation' taxa and transform into relative abundance data
groupData <- subset_taxa(Biom23_F_cut, Phylum != "Multi-affiliation") # remove unknwon phylum
groupData_relab <- transform_sample_counts(groupData, function(x) x/sum(x)) # calculate relative abundances

# merge samples
melted <- psmelt(groupData_relab)

#check that all sample sum is equal to 1
melted %>%
  group_by(Sample) %>%
  summarise_at(c("Abundance"), sum)

#save results
groupData23 <- as.data.frame(melted)

### group RF and P data together to have siilar colors for Phylum

df <- groupData23 %>% group_by(Order, Class, Phylum) %>% 
  summarise_at(c("Abundance"), sum)

df <- as.data.frame(df)
df=df[order(df$Phylum, df$Abundance),]
length(unique(df$Phylum)) # to have the number of colors
magma(11) #
df$cols <- c("#000004FF", "#150E37FF", "#3B0F70FF", "#641A80FF", "#8C2981FF", "#B63679FF", "#DE4968FF", "#F76F5CFF", "#FE9F6DFF", "#FECE91FF",
             "#FCFDBFFF")[match(df$Phylum,sort(unique(df$Phylum)))]

#subset rainforest
df3 <- df %>% group_by(Phylum) %>% 
  summarise_at(c("Abundance"),sum)

metadata3 <- df %>%
  mutate(Tot = sum(Abundance)) %>% 
  group_by(Phylum, Class) %>% 
  mutate(CUM = cumsum(Abundance), DomSize = max(CUM)) %>%
  arrange(Class, Phylum)

df4 <- metadata3 %>%
  group_by(Phylum, Class) %>% 
  summarise(n=sum(Abundance)) %>%
  arrange(Class, Phylum)

df5 <- df4
df5$per <- df5$n/sum(df5$n)

theme_set(new = theme_classic())
plt <- ggplot() + geom_col(aes(x = 2, y = Abundance, fill = Phylum), 
                           data = df3, color = "black") + 
  geom_col(aes(x = 3, y = n, fill = Phylum), 
           data = df4, color = "black") + 
  geom_col(aes(x = 4, y = Abundance, fill = Phylum), 
           data = metadata3, color = "black") +
  xlim(0, 4.5) + labs(x = NULL, y = NULL) + 
  scale_fill_manual(values = c("#000004FF", "#150E37FF", "#3B0F70FF", "#641A80FF", "#8C2981FF", "#B63679FF", "#DE4968FF", "#F76F5CFF", "#FE9F6DFF", "#FECE91FF",
                               "#FCFDBFFF"))

plt23 <- plt + coord_polar(theta = "y") 
plt23

################### Oridantion plot ####################################################################################################

###
### NMDS analyses ----
###
library(parallelDist)

# NMDS analysis 23S
dist.23S <- parDist(as.matrix(otu), method = "bray")
sp23s = metaMDS(dist.23S, k = 5)
stressplot(sp23s)
plot(sp23s)

#extract site scores
sites.sp23s = sp23s$points
sites.sp23s = as.data.frame(sites.sp23s)

#prepare for plotting
srcs_traitcompsp<-cbind(sites.sp23s, familily = metadata23S$Code[1:53])
head(srcs_traitcompsp)
cent.tr.sp<-aggregate(cbind(MDS1, MDS2)~familily, data=srcs_traitcompsp, FUN=mean)
segs.tr.sp<-merge(srcs_traitcompsp, setNames(cent.tr.sp, c('familily','Axis1', 'Axis2')),
                  by='familily', sort=FALSE)


#plot
find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls <- ddply(srcs_traitcompsp,"familily", find_hull)

group_col<-c("#E69F00", "#FFCC66", "#3399CC", "#3366FF", "#FF6699", "#CC0033")
theme_set(new = theme_classic())
nmds23S <- ggplot(srcs_traitcompsp, aes(MDS1, MDS2, color =  familily, fill=familily)) +
  geom_point(aes(MDS1, y = MDS2), size = 5, alpha = I(1/2)) +
  geom_polygon(data=hulls, alpha = .25) +
  scale_fill_manual(values= group_col) +
  scale_color_manual(values= group_col) 
nmds23S

## MAKE AN ANOSIM behind to test for treatments and site effects

################### Co-occurrence network analyses ################# ----

###
### spread OTU data into a wide form
###

#23S group similar taxo by sample
groupData_23S_group <- groupData23 %>%
  group_by(Sample, OTU) %>%
  summarise_at(c("Abundance"), sum)

groupData_23S.wide <- as.data.frame(spread(groupData_23S_group %>% select(Sample, OTU, Abundance), OTU, Abundance, fill = 0))
row.names(groupData_23S.wide) <- groupData_23S.wide$Sample
groupData_23S.wide <- groupData_23S.wide %>% select(-Sample)

groupData_23S.widefil <- groupData_23S.wide[,colSums(groupData_23S.wide)>0] ### selecte min number of co variables

###
### PA data
###

groupData_23Spa <- decostand(groupData_23S.widefil, 'pa')

###
### Cor mat
###
source("test.a.R")
library(Hmisc); library(igraph)

rmat <- rcorr(as.matrix(groupData_23S.widefil))
View(rmat$P)

### distance matrix
res <- test.a(as.matrix(groupData_23Spa), nperm = 9) # do nperm = 999 ... can be bloody long :)
summary(res)
is.list(res)
#res.final=res
#save(res, file="Cluster_res.RData")
#load(file="Cluster_res.RData")

# Adjacency matrix from the "a" distance matrix 
adjm2 <- 1 - as.matrix(res$p.a.dist)

### Adjacency matrix from the Spearman rank correlation matrix
adjm3 <- cor(groupData_23S.widefil, method = "spearman")

#second method with p values
adjm3b=rcorr(as.matrix(groupData_23S.widefil), type = "spearman")
plot(adjm3b$r~adjm3b$P, pch = 20, col= ifelse(adjm3b$P <= 0.05,"red", "black"), xlab = "P-value", ylab = "rho") # this plot takes a lot of time too! 
abline(h = 0.45); abline(h = -0.45) # rho > 0.35 => all significant

# Only positive associations (rho > 0.2)
adjm2[adjm3b$P <= 0.05] <- 1
adjm2[adjm3b$P > 0.05] <- 0  # binary co-occurrence

# threshold 2, strength of the 'interaction'
adjm2[adjm2 < 0.4] <- 0 # check the plot
diag(adjm2) <- 0

# Only significang associations
adjm2b=adjm2
adjm2b[is.nan(adjm3b$P)] <- 0
adjm2b[adjm2==1 & adjm3b$P <= 0.05] <- 1
adjm2b[adjm2==1 & adjm3b$P > 0.05] <- 0  # binary co-occurrence

### Select an adjacency matrix
adjm <- adjm2

# Build graph
go <- graph_from_adjacency_matrix(adjm, 
                                  weighted = TRUE, 
                                  mode = "undirected")

#names
names23S <- as.data.frame(colnames(groupData_23Spa))
names23S$Code <- str_extract(names23S$`colnames(groupData_23Spa)`, "[a-z]{1,4}")
names23S$ab <- colSums(groupData_23S.widefil)


# Network structure detection: find densely connected subgraphs (modules or "communities") in a graph
#different options:
#wc <- cluster_walktrap(go)
wc <- cluster_edge_betweenness(go)
wc <- cluster_fast_greedy(go)
#wc <- cluster_leading_eigen(go)
wc <- cluster_spinglass(go)
#wc <- cluster_optimal(go)

set.seed(1)
plot(wc, go, vertex.label=NA, vertex.size=names23S$ab*2)


################### Richness and diversity maps #########


###
### 23S
###
spe <- otu
N0 = rowSums(spe > 0)				# Species richness
Chao1 <- apply(spe, 1, fossil::chao1) #richness using Chao1
H = diversity(spe)					# Shannon's entropy
N1 = exp(H)									# Shannon's diversity index
N2 = diversity(spe, "inv")	# Simpson's diversity index
J = H/log(N0)								# Pielou's evenness
E1 = N1 / N0								# Hill's evenness
E2 = N2 / N0								# Hill's evenness
simpson(spe)
div.23s = data.frame(N0=N0, Chao1 = Chao1, H = H, N1=N1, N2=N2, E1=E1, E2=E2, J=J, Chao1 = Chao1)


#quick stats
div.23s %>%
  summarise_at(c("N1"), c(mean, sd))

################### Responses at species (or whatever) level ######### ----

library(nlme); library(effectsize)

db_species <- groupData23 %>% group_by(Phylum, Sample) %>% 
  summarise_at(c("Abundance"), sum)
levels(as.factor(groupData23$Phylum))

env <- metadata23S[1:53,] %>% dplyr::select(id, Site, Treatment)
colnames(env) <- c("Sample", "Site", "Treatment")
species23 <- left_join(env, db_species, by="Sample")

#model to test responses
mod<-lm(Abundance ~ Site*Treatment, data = species23, na.action=na.omit) 
anova(mod); #summary(mod)
hist(mod$residuals)

#per site - loop
list_phylum23 <- split(species23, f = species23$Phylum)
species <- 1:11 # 4773 is your total number of species in the list to know it use this line: list_species23 <- split(species23, f = species23$OTU)

output <- lapply(X = species, FUN = function(x){
  mod <- lm(Abundance ~ Site*Treatment, data = list_phylum23[[x]], na.action=na.omit) 
  fvalueS <- anova(mod)$`F value`[1]
  fvalueT <- anova(mod)$`F value`[2]
  fvalueST <- anova(mod)$`F value`[3]
  pvalueS <- anova(mod)$`Pr(>F)`[1]
  pvalueT <- anova(mod)$`Pr(>F)`[2]
  pvalueST <- anova(mod)$`Pr(>F)`[3]
  esS <- effectsize::eta_squared(mod, partial = T)$Eta2_partial[1]
  esT <- effectsize::eta_squared(mod, partial = T)$Eta2_partial[2]
  esST <-effectsize::eta_squared(mod, partial = T)$Eta2_partial[3]
  lesS <- effectsize::eta_squared(mod, partial = T)$CI_low[1]
  lesT <- effectsize::eta_squared(mod, partial = T)$CI_low[2]
  lesST <- effectsize::eta_squared(mod, partial = T)$CI_low[3]
  uesS <- effectsize::eta_squared(mod, partial = T)$CI_high[1]
  uesT <- effectsize::eta_squared(mod, partial = T)$CI_high[2]
  uesST <-effectsize::eta_squared(mod, partial = T)$CI_high[3]
  ind.names <- c("esS", "lesS", "uesS", 'pvalueS', 'esT', 'lesT', 'uesT', 'pvalueT', 'lesST', 'uesST', 'pvalueST')
  ind.values <- c(unname(esS), lesS, uesS, pvalueS, unname(esT), lesT, uesT, pvalueT, lesST, uesST, pvalueST)
  es_mat <- as.data.frame(cbind(ind.names, ind.values))
  es_mat$name <- list_phylum23[[x]]$Phylum[1]
  return(es_mat)
})

output
es_species23 <- do.call(rbind, output)
es_species23 <- spread(es_species23, ind.names, ind.values) #only one species respond to P => Cluster 11 (Anabaenae, ES = -0.26)

TukeyHSD(aov(Abundance ~ Site*Treatment, data = list_phylum23[[2]], na.action=na.omit))

##plot
# Change the color by groups
es_classes$Domain = factor(es_classes$Domain, levels = c("prokaryotes", "eukaryotes"))
es_classes$name = factor(es_classes$name, levels = c("Acidobacteriae", "Bacteroidia", "Cyanobacteriia", "Planctomycetes", "Polyangia", "Vampirivibrionia", "Chlorophyceae",    
                                                     "Chrysophyceae", "Colpodea", "Dinophyceae", "Multi-affiliation", "Trebouxiophyceae", "Zygnemophyceae"))

esTplot <- ggplot(es_classes, aes(x = esT, y = forcats::fct_rev(factor(name)), xmin = esT-lesT, xmax = esT+uesT, color = Domain)) +
  geom_point() +
  geom_errorbarh(height=0) + xlim(-0.6, 0.6)
esPplot <- ggplot(es_classes, aes(x = esP, y = forcats::fct_rev(factor(name)), xmin = esP-lesP, xmax = esP+uesP, color = Domain)) +
  geom_point() + xlim(-0.6, 0.6) + geom_errorbarh(height=0)

gridExtra::grid.arrange(esTplot, esPplot, nrow=1)

###
#### co-occurrence networking using jSDM approach ----
###

## install JAGS first, then install package rjags and boral
## for reference, check Robroek et al. 2017 Nat Com. Sytiuk et al. 2022 Oikos; Pollock et al. 2014 Methods Ecolv. Evol. 

library(boral); library(corrplot)

#filter species by removing rare OTUs
hist(colSums(groupData_23S.widefil)) # remove 10% of data to avoid rare mtbs
species_filtered <- groupData_23S.widefil[,colSums(groupData_23S.widefil)>0.07] ### remove rare compounds (10% threshold)
hist(colSums(species_filtered)) # remove 10% of data to avoid rare mtbs    

# preprare env
env <- metadata23S[1:53,] %>% dplyr::select(pH, N, WT, Shade)
env$pH <- as.numeric(env$pH)
env$N <- as.numeric(env$N)
env$WT <- as.numeric(env$WT)
env$Shade <- as.numeric(env$Shade)


# w/o stem water content / w/o location effect THIS MODEL!!!
fit.lvm2 <- boral(y = species_filtered, X =env, lv.control = list(num.lv = 2, type = "independent", distmat = NULL), family = "normal", row.eff = "fixed", 
                       save.model = TRUE, calc.ics = F, mcmc.control = list(n.brunin =100, n.iterations = 1000, n.thin = 10, seed = 1), hypparams = c(10,10,10,30)) # To look at estimated parameter values
     
     
#$hpdintervals$lv.coefs # 95% credible intervals for model parameters.
     
## Dunn-Smyth residual plots  to check model assumption, outliers etc...
plot(fit.lvm2)
     
     
## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lvm2)
     
## Use corrplot package to plot residual correlations between species, e.g. possibly due to species interaction etc... ----
res.cors <- get.residual.cor(fit.lvm2, est = "median", prob = 0.95)
     
corrplot(res.cors$cor.upper, diag = F, type = "upper", title = "Residual correlations from LVM", method = "color", tl.srt = 45)
corrplot(res.cors$sig.cor, diag = F, type = "upper", title = "Residual correlations from LVM", method = "color", tl.srt = 45)

## Note correlations are similar to those from the GLMM, the main difference being stronger correlations for the
## two points at the bottom-right of the plot. Adding an extra latent variable to the model would fix this (num.lv=3). 
     
## Use corrplot package to plot correlations between species  DUE TO environmental variables ----
env.cors <- get.enviro.cor(fit.lvm2, est = "median", prob = 0.95)
     
corrplot(env.cors$cor.upper, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)
corrplot(env.cors$sig.cor, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)
     
     
### linkages with env parameters
model<-as.data.frame(fit.lvm2$X.coefs.mean)
write.table(model, file="jsdm_model_coefficients.csv", sep=";", quote=FALSE, col.names=NA)
     
vec<-c(1,2,1,4,4,4,5,5,4,2,1,2,2,3,3,3,2,2,2,4) # cluster assignations
model$cluster<-vec
     
clusters <- gather(model, key = "env", value = "value", -cluster)
     
theme_set(new = theme_classic())
ggplot(clusters, aes(x = env, y = value)) + 
geom_boxplot(aes(fill = as.factor(cluster))) + scale_y_continuous(name = "Coefficient", breaks = seq(-2, 2, 0.5),
                                                                         limits=c(-2,2)) +
  scale_x_discrete(name = "cluster") + coord_flip() + geom_hline(yintercept = 0, linetype = "dotted")
     
     
##########-----------------------------------------------------------------------------------------------
### (4) Co-occurrence network of residuals and environmental correlations among species ----
##########-----------------------------------------------------------------------------------------------

library(igraph)

# Build graph
cor.A = as.data.frame(env.cors$sig.cor) # only positive correlations

cor.net = cor.A
diag(cor.net) = 0
     
go <- graph_from_adjacency_matrix(as.matrix(cor.net), 
                                       weighted = TRUE, 
                                       mode = "undirected")
     
#par(mfrow=c(2,2), mar=c(1,1,1,1))
     
# Network structure detection: find densely connected subgraphs (modules or "communities") in a graph
#different options:
set.seed(23)
     
wc <- cluster_walktrap(go)
plot(wc, go, vertex.label=wc$names,vertex.size=4,vertex.label = 0.5,layout = layout_with_graphopt)
     
wc1 <- cluster_edge_betweenness(go)
plot(wc1, go, vertex.label=wc$names)
     
wc2 <- cluster_fast_greedy(go)
plot(wc2, go, vertex.label=wc$names)
   
wc3 <- cluster_leading_eigen(go)
plot(wc3, go, vertex.label=wc$names)
     
wc4 <- cluster_spinglass(go)
plot(wc4, go, vertex.label=wc$names)
     
wc5 <- cluster_optimal(go) # good one
plot(wc5, go, vertex.label=wc$names)
     
     
cor.B = as.data.frame(res.cors$sig.cor) # only positive correlations
#cor.A = as.data.frame(res.cors$cor.upper) # only positive correlations
cor.oc2 = cor.B
   
cor.net2 = cor.oc2
diag(cor.net2) = 0
   
go2 <- graph_from_adjacency_matrix(as.matrix(cor.net2), weighted = TRUE, 
                                        mode = "undirected")

#par(mfrow=c(2,2), mar=c(1,1,1,1)
# Network structure detection: find densely connected subgraphs (modules or "communities") in a graph
#different options

set.seed(23)
     
wc <- cluster_walktrap(go2)
plot(wc, go2, vertex.label=wc$names,vertex.size=4,vertex.label = 0.5,layout = layout_with_graphopt)
     
wc2 <- cluster_edge_betweenness(go2)
plot(wc2, go2, vertex.label=wc$names)
     
wc2 <- cluster_fast_greedy(go2)
plot(wc2, go2, vertex.label=wc$names)
     
wc3 <- cluster_leading_eigen(go2)
plot(wc3, go2, vertex.label=wc$names)

