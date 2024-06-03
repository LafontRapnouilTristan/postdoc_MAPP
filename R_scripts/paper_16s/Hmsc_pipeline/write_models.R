set.seed(973)

library("easypackages")
libraries("Hmsc","readr","dplyr","magrittr","tidyr","jsonify","ggplot2","patchwork")

getwd()


# SET DIRECTORIES 
localDir = "R_scripts/paper_16s/Hmsc_pipeline/"
modelDir = file.path(localDir, "models")

if(!dir.exists(modelDir)) dir.create(modelDir)


# READ AND EXPLORE THE DATA

load("R_scripts/paper_16s/outputs/meco_16s.RDS")
load("R_scripts/paper_16s/outputs/site_data_imputed.RDS")
site_data <- combined
rm(combined)

Y <- t(meco_16s_ns$otu_table)
rownames(Y)<-  gsub("[A-Z]|-|(?=_).*","",rownames(Y),perl = T)  # get site name (i.e., sample)
samp_to_keep <- rownames(Y)
Y <- Y[order(rownames(Y)),]

length(unique(site_data$Location)) # 181 unique locations for 199 samples
length(unique(site_data$X)) # 181 unique Lon
length(unique(site_data$Y)) # 180 unique Lat

site_data$Y_jitd <- site_data$Y# Random noise (<1m lat) in coordinates to fit HMSC
site_data$Y_jitd[which(duplicated(site_data$Y))] <- jitter(site_data$Y[which(duplicated(site_data$Y))],amount = 8.983112*10^-6)

which(duplicated(site_data$Y_jitd)) # check that we have no duplication anymore
length(unique(paste0(site_data$Y_jitd,site_data$X))) # 199

site_data$Location <- cumsum(!duplicated(site_data[c("X","Y_jitd")])) # ID for unique lon lat combination

site_data %<>% mutate(Sample=as.factor(Sample),
                      Location=as.factor(Location))
rownames(site_data) <- site_data$Sample
rownames(as.data.frame(Y))==rownames(site_data) #check if same order


# distribution of Transformed Y data 
test_trans <- 
Y%>%
    data.frame()%>%
    reshape2::melt()%>%
    ggplot(aes(x=value))+
    geom_histogram(bins = 100)+ 
    scale_y_log10()+
    ggtitle("Raw reads*sample distribution")+
    ylab("Log transformed counts")+
    xlab("reads count")+
Y%>%
    data.frame()%>%
    labdsv::hellinger()%>%
    reshape2::melt()%>%
    ggplot(aes(x=value))+
    geom_histogram(bins = 100)+ 
    scale_y_log10()+
    ggtitle("Hellinger reads*sample distribution")+
    ylab("Log transformed counts")+
    xlab("hellinger transformed reads")+
Y%>%
    data.frame()%>%
    mutate(across(everything(),.fns=~ifelse(.x==0,NA,.x)))%>%
    mutate(across(everything(),.fns=~log(.x)))%>%
    mutate(across(everything(),.fns=~scale(.x)))%>%
    reshape2::melt()%>%
    ggplot(aes(x=value))+
    geom_histogram(bins = 100)+
    ggtitle("Non null logged and scaled reads*sample distribution")
  

# Transform Y

Ypa <- 1*(Y>0) # convert Y matrix to presence-absence
Yabu <- Y
Yabu[Yabu == 0] <- NA # Convert zero values in zero-inflated Yabu columns to NA
Yabu <- log(Yabu)
Yabu <- scale(Yabu)
Yraw <- Y # rename to explicitly have Raw com
Yhell <- labdsv::hellinger(Y) # Hellinger transformation

# SET UP THE MODEL

# environmental data into a dataframe XData
XData <- site_data
XData$Depth <- log(rowSums(Y)) # logged sequencing depth as covariates

# Define the environmental model through XFormula

var_lt <- colnames(XData)[which(grepl("lt",colnames(XData)))]
var_st <- colnames(XData)[which(grepl("st",colnames(XData)))]

var_interest <- c("aet","pr","swe","tmmn","tmmx","vpd","elevation","evi","fpar","gpp","ph","soc")

var_lt <- var_lt[which(gsub("_lt","",var_lt)%in%var_interest)]
var_st <- var_st[which(gsub("_st","",var_st)%in%var_interest)]

XFormula_lt <- paste("XData$",var_lt,sep="", collapse = "+")
XFormula_st <- paste("XData$",var_st,sep="", collapse = "+")

XFormula_lt <- paste0("~ XData$ph +", XFormula_lt ,"+ XData$Depth")
XFormula_st <- paste0("~ XData$ph +", XFormula_st ,"+ XData$Depth")
XFormula_null <- "~ XData$Depth" # NULL model



# Define the studyDesign as a dataframe 
studyDesign <- data.frame(Coords=as.factor(1:nrow(Y))) 
# Set up the random effects

# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
xycoords <- site_data[,c("X","Y_jitd")]  # get spatial data removing duplicate rows
xycoords <- as.matrix(xycoords)
rownames(xycoords) <- studyDesign$Coords # Rownames as locations unique ID

rL.spatial = HmscRandomLevel(sData = xycoords) # Spatial random effect to account for space induced patterns


# Use the Hmsc model constructor to define a model
# note that in the random effects the left-hand sides in the list (year, location, sample)
# refer to the columns of the studyDesign

# models for different subsets of the data
XData %<>%select(-c("Sample","Country")) #remove country and sample_id from X
XFormula_list <- c(XFormula_lt,XFormula_st,XFormula_null) # list with my formulas
names(XFormula_list) <- c("XFormula_lt","XFormula_st","XFormula_null")
m_list <- NULL

m_list <- NULL
for(i in 1:length(XFormula_list)){
    
    # hurdle PA
    mpa <- Hmsc(Y=Ypa, XData = XData,  XFormula = as.formula(XFormula_list[i]), distr = "probit", studyDesign=studyDesign, ranLevels=list(Coords = rL.spatial))
    # Hurdle abundance
    mabu <- Hmsc(Y=Yabu, XData = XData,  XFormula = as.formula(XFormula_list[i]), distr = "normal", studyDesign=studyDesign, ranLevels=list(Coords = rL.spatial))
    # Classic model
    mclassic <- Hmsc(Y=Yraw, XData = XData,  XFormula = as.formula(XFormula_list[i]), distr = "lognormal poisson", studyDesign=studyDesign, ranLevels=list(Coords = rL.spatial))
    # Classic model with scaled Y
    mclassicscale <- Hmsc(Y=Yraw, YScale = T, XData = XData,  XFormula = as.formula(XFormula_list[i]), distr = "lognormal poisson", studyDesign=studyDesign, ranLevels=list(Coords = rL.spatial))
    # Classic model
    mhellinger <- Hmsc(Y=Yhell, XData = XData,  XFormula = as.formula(XFormula_list[i]), distr = "lognormal poisson", studyDesign=studyDesign, ranLevels=list(Coords = rL.spatial))
    # Classic model with scaled Y
    mhellingerscale <- Hmsc(Y=Yhell, YScale = T, XData = XData,  XFormula = as.formula(XFormula_list[i]), distr = "lognormal poisson", studyDesign=studyDesign, ranLevels=list(Coords = rL.spatial))
    
    m_list[[names(XFormula_list)[i]]][["mpa"]] <- mpa
    m_list[[names(XFormula_list)[i]]][["mabu"]] <- mabu
    m_list[[names(XFormula_list)[i]]][["mclassic"]] <- mclassic
    m_list[[names(XFormula_list)[i]]][["mclassicscale"]] <- mclassicscale
    m_list[[names(XFormula_list)[i]]][["mhellinger"]] <- mhellinger
    m_list[[names(XFormula_list)[i]]][["mhellingerscale"]] <- mhellingerscale
    
}

#  SAVING MODELS 
save(m_list, file = file.path(modelDir, "unfitted_models_V1.RData"))

# TESTING THAT MODELS FIT WITHOUT ERRORS
mod_type <- c("mpa","mabu","mclassic","mclassicscale","mhellinger","mhellingerscale")
test_mod <- NULL
for(i in 1:length(XFormula_list)){
    for(j in mod_type){
        print(paste0(names(XFormula_list)[i]," : ",j))
        m_list[[names(XFormula_list)[i]]][[j]] <- sampleMcmc(m_list[[i]][[j]],samples=2,verbose = 1)
    }
}
