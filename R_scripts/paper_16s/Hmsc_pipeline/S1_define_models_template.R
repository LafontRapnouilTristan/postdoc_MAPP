# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Original datafiles of the case study, placed in the data folder.

#	OUTPUT. Unfitted models, i.e., the list of Hmsc model(s) that have been defined
# but not fitted yet, stored in the file "models/unfitted_models.RData".
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(973)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# LOAD PACKAGES (BEGINNING)
##################################################################################################
EcophyCofog::Library(c("Hmsc","readr","dplyr","magrittr","tidyr"))

##################################################################################################
# LOAD PACKAGES (END)
##################################################################################################


##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "./R_scripts/paper_16s/Hmsc_pipeline"
modelDir = file.path(localDir, "models")

if(!dir.exists(modelDir)) dir.create(modelDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################


##################################################################################################
# READ AND EXPLORE THE DATA (BEGINNING)
##################################################################################################
# Write here the code needed to read in the data, and explore it by (with View, plot, hist, ...)
# to get and idea of it and ensure that the data are consistent
load("./R_scripts/paper_16s/outputs/meco_16s.RData")
load("./R_scripts/paper_16s/outputs/site_data_imputed.RDS")
site_data <- combine
rm(combine)

Y <- t(meco_16s_ns$otu_table)
rownames(Y)<-  gsub("[A-Z]|-|(?=_).*","",rownames(Y),perl = T)  # get site name (i.e., sample)
samp_to_keep <- rownames(Y)
Y <- Y[order(rownames(Y)),]


site_data$Location <- cumsum(!duplicated(site_data[c("X","Y")])) # ID for unique lon lat combination
length(unique(site_data$Location)) # 181 unique locations for 199 samples
length(unique(site_data$X)) # 181 unique Lon
length(unique(site_data$Y)) # 180 unique Lat

site_data %<>% mutate(Sample=as.factor(Sample),
                      Location=as.factor(Location))
rownames(site_data) <- site_data$Sample
rownames(as.data.frame(Y))==rownames(site_data) #check if same order
##################################################################################################
# READ AND EXPLORE THE DATA (END)
##################################################################################################

##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################
# Note that many of the components are optional 
# Here instructions are given at generic level, the details will depend on the model

# Organize the community data in the matrix Y
Y <- Y # to be explicit if not enough yet
# Organize the environmental data into a dataframe XData
XData <- site_data
XData$Depth <- rowSums(Y)
# Define the environmental model through XFormula
# XFormula = ~ ...
var_lt <- colnames(XData)[which(grepl("lt",colnames(XData)))]
var_st <- colnames(XData)[which(grepl("st",colnames(XData)))]
var_interest <- c("aet","pr","swe","tmmn","tmmx","vap","vpd","moist","ssm","elevation","evi","fpar","gpp","ghm","ph","water","soc")
var_lt <- var_lt[which(gsub("_lt","",var_lt)%in%var_interest)]
var_st <- var_st[which(gsub("_st","",var_st)%in%var_interest)]

XFormula_lt <- paste("XData$",var_lt,sep="", collapse = "+")
XFormula_st <- paste("XData$",var_st,sep="", collapse = "+")

XFormula_lt <- paste0("~ XData$ph +", XFormula_lt ,"+ XData$Depth")
XFormula_st <- paste0("~ XData$ph +", XFormula_st ,"+ XData$Depth")
XFormula_null <- "~ XData$Depth" # NULL model?

# Organize the trait data into a dataframe TrData

# No trait Data

# Set up a phylogenetic (or taxonomic tree) as myTree

# No phylogenetic tree cuz 16s is too weak?

# Define the studyDesign as a dataframe 
studyDesign <- data.frame(XData[,c("Location","Sample")]) # we have both locations sampled multiple times and unique samples. 

# Set up the random effects

# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
sData <- as.matrix(distinct(site_data[,c("X","Y")]))  # get spatial data removing duplicate rows
rownames(sData) <- unique(studyDesign$Location) # Rownames as locations unique ID
rL.location = HmscRandomLevel(sData = sData) # Spatial random effect to account for space induced patterns

# For another example, you may define the sample = sampling unit = row of matrix Y
# as a random effect, in case you are interested in co-occurrences at that level
# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
rL.Sample = HmscRandomLevel(units = levels(studyDesign$Sample)) # Sample random effect to account for sample identity
                                                                # Set priors???

# Use the Hmsc model constructor to define a model
# note that in the random effects the left-hand sides in the list (year, location, sample)
# refer to the columns of the studyDesign

# It is always a good idea to look at the model object, so type m to the console and press enter
# Look at the components of the model by exploring m$...

# You may define multiple models, e.g.
# alternative covariate selections
# alternative random effect choices
# alternative ways to model the data (e.g. lognormal Poisson versus hurdle model)
# models for different subsets of the data

XFormula_list <- c(XFormula_lt,XFormula_st,XFormula_null)
names(XFormula_list) <- c("XFormula_lt","XFormula_st","XFormula_null")
m_list <- NULL
for (i in 1:length(XFormula_list)){
    m <- Hmsc(Y=Y,
             XData = XData, 
             XFormula = as.formula(XFormula_list[i]),
             distr="lognormal poisson" ,
             studyDesign=studyDesign, 
             ranLevels=list(Sample=rL.Sample,Location = rL.location))

    m_list[[names(XFormula_list)[i]]] <- m
}


##################################################################################################
# SET UP THE MODEL (END)
##################################################################################################


##################################################################################################
# In the general case, we could have multiple models. We combine them into a list and given them names.
# COMBINING AND SAVING MODELS (START)
names(m_list) <- c("longterm","shortterm","null")
save(m_list, file = file.path(modelDir, "unfitted_models.RData"))
##################################################################################################
# COMBINING AND SAVING MODELS (END)
##################################################################################################


##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
##################################################################################################
for(i in 1:length(m_list)){
 print(i)
 sampleMcmc(m_list[[i]],samples=2)
}
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)
##################################################################################################
