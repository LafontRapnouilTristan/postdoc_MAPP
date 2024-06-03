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
EcophyCofog::Library(c("Hmsc","readr","dplyr","magrittr","tidyr","jsonify"))

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
site_data <- combined
rm(combined)

Y <- t(meco_16s_ns$otu_table)
rownames(Y)<-  gsub("[A-Z]|-|(?=_).*","",rownames(Y),perl = T)  # get site name (i.e., sample)
samp_to_keep <- rownames(Y)
Y <- Y[order(rownames(Y)),]

site_data$Y_jitd <- site_data$Y# Random noise (<1m lat) in coordinates to fit HMSC
site_data$Y_jitd[which(duplicated(site_data$Y))]    <- jitter(site_data$Y[which(duplicated(site_data$Y))],amount = 8.983112*10^-6)

which(duplicated(site_data$Y_jitd)) # check duplication

site_data$Location <- cumsum(!duplicated(site_data[c("X","Y_jitd")])) # ID for unique lon lat combination
length(unique(site_data$Location)) # 181 unique locations for 199 samples
length(unique(site_data$X)) # 181 unique Lon
length(unique(site_data$Y)) # 180 unique Lat

site_data %<>% mutate(Sample=as.factor(Sample),
                      Location=as.factor(Location))
rownames(site_data) <- site_data$Sample
rownames(as.data.frame(Y))==rownames(site_data) #check if same order


# Filter coms on prev and relab
prevalence <- c(0,1,2,3,4,5,10,15,20) # set prevalence threshold 0/5/10/15 and 20% of our 199 samples 
n_sites <- 199 # 199 samples
prevalence_nsites <- floor(n_sites*0.01*prevalence) # convert to number of sites
relab_thrshld <- c(0,1e-3,1e-5,1e-7)

com_matrices_list <- NULL
for(i in prevalence_nsites){
    matrix_tmp <- t(meco_16s_ns$otu_table)   
    
    matrix_tmp <- matrix_tmp[,which(colSums(matrix_tmp!=0)>i)] # remove cluster with prevalence < i
    for(j in relab_thrshld){
        matrix_tmp_relab <- matrix_tmp/rowSums(matrix_tmp) # compute relab
        matrix_tmp_relab[matrix_tmp_relab==0]  <- NA
        matrix_tmp <- matrix_tmp[,which(colMeans(matrix_tmp_relab,na.rm = T)>j)] # filter when mean relab of cluster is < threshold
        matrix_tmp <- matrix_tmp[rowSums(matrix_tmp)!=0,] # remove samples with no reads anymore?
        rownames(matrix_tmp) <- gsub("[A-Z]|-|(?=_).*","",rownames(matrix_tmp),perl = T) # rename
        com_matrices_list[[paste0("com_prev_",i,"_relab_",j)]] <- matrix_tmp # save com
    }
}

##################################################################################################
# READ AND EXPLORE THE DATA (END)
##################################################################################################

##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################
# Note that many of the components are optional 
# Here instructions are given at generic level, the details will depend on the model

# Organize the community data in the matrix Y
Y <- com_matrices_list$com_prev_39_relab_0.001 # to be explicit if not enough yet
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
studyDesign <- data.frame(Sample=XData$Sample) 
# Set up the random effects

# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
sData <- site_data[,c("X","Y_jitd")]  # get spatial data removing duplicate rows
names(sData) <- c("x-coodinate","y-coordinate")
rownames(sData) <- studyDesign$Sample # Rownames as locations unique ID
rL.spatial = HmscRandomLevel(sData = sData) # Spatial random effect to account for space induced patterns

# For another example, you may define the sample = sampling unit = row of matrix Y
# as a random effect, in case you are interested in co-occurrences at that level
# rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
# rL.Sample = HmscRandomLevel(units = levels(studyDesign$Sample)) # Sample random effect to account for sample identity
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
XData %<>%select(-c("sample_id","Country"))
XFormula_list <- c(XFormula_lt,XFormula_st,XFormula_null)
names(XFormula_list) <- c("XFormula_lt","XFormula_st","XFormula_null")
m_list <- NULL
for (i in 1:length(XFormula_list)){
    m <- Hmsc(Y=Y,
             XData = XData, 
             XFormula = as.formula(XFormula_list[i]),
             distr="lognormal poisson" ,
             studyDesign=studyDesign, 
             ranLevels=list(Sample = rL.spatial))

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

for(i in 1:length(m_list)){
    print(i)
    sampleMcmc(m_list[[i]],samples=2, updater = list(GammaEta = FALSE))
}


