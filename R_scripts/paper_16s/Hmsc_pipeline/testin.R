# testin' HMSC


library("easypackages")
libraries("Hmsc","readr","dplyr","magrittr","tidyr","jsonify","ggplot2","patchwork")

getwd()


# SET DIRECTORIES 
localDir = "R_scripts/paper_16s/Hmsc_pipeline/"
modelDir = file.path(localDir, "models")

if(!dir.exists(modelDir)) dir.create(modelDir)


# READ AND EXPLORE THE DATA

load("R_scripts/paper_16s/outputs/meco_16s_filtered.RData")
load("R_scripts/paper_16s/outputs/site_data_imputed.RDS")
site_data <- combined
rm(combined)

Y <- t(meco_16s_relabfiltered$otu_table)
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

# subset data
set.seed(973)
sites <- sample(1:199,50,replace = F)
clusters <- sample(1:822,100,replace = F)
Ysub <- Y[sites,clusters]
# distribution of Transformed Y data 
test_trans <- 
  Ysub%>%
  data.frame()%>%
  reshape2::melt()%>%
  ggplot(aes(x=value))+
  geom_histogram(bins = 100)+ 
  scale_y_log10()+
  ggtitle("Raw reads*sample distribution")+
  ylab("Log transformed counts")+
  xlab("reads count")+
  Ysub%>%
  data.frame()%>%
  labdsv::hellinger()%>%
  reshape2::melt()%>%
  ggplot(aes(x=value))+
  geom_histogram(bins = 100)+ 
  scale_y_log10()+
  ggtitle("Hellinger reads*sample distribution")+
  ylab("Log transformed counts")+
  xlab("hellinger transformed reads")+
  Ysub%>%
  data.frame()%>%
  mutate(across(everything(),.fns=~ifelse(.x==0,NA,.x)))%>%
  mutate(across(everything(),.fns=~log(.x)))%>%
  mutate(across(everything(),.fns=~scale(.x)))%>%
  reshape2::melt()%>%
  ggplot(aes(x=value))+
  geom_histogram(bins = 100)+
  ggtitle("Non null logged and scaled reads*sample distribution")


# Transform Y

Ypa <- 1*(Ysub>0) # convert Y matrix to presence-absence
Yabu <- Ysub
Yabu[Yabu == 0] <- NA # Convert zero values in zero-inflated Yabu columns to NA
Yabu <- log(Yabu)
Yabu <- scale(Yabu)
Yraw <- Ysub # rename to explicitly have Raw com
Yhell <- labdsv::hellinger(Ysub) # Hellinger transformation

# SET UP THE MODEL ############

# environmental data into a dataframe XData
XData <- site_data
XData <- XData[sites,]
XData$Depth <- log(rowSums(Ysub)) # logged sequencing depth as covariates

# Define the environmental model through XFormula

var_lt <- colnames(XData)[which(grepl("lt",colnames(XData)))]

var_interest <- c("pr","swe","tmmn","vpd","elevation","gpp","ph","soc")

var_lt <- var_lt[which(gsub("_lt","",var_lt)%in%var_interest)]

XFormula_lt <- paste("XData$",var_lt,sep="", collapse = "+")

XFormula_lt <- paste0("~ XData$ph +", XFormula_lt ,"+ XData$Depth")
XFormula_null <- "~ XData$Depth" # NULL model


# Define the studyDesign as a dataframe 
studyDesign <- data.frame(Coords=as.factor(1:nrow(Ysub))) 
# Set up the random effects

# Here locations.xy would be a matrix (one row per unique location)
# where row names are the levels of studyDesign$location,
# and the columns are the xy-coordinates
xycoords <- XData[,c("X","Y_jitd")]  # get spatial data removing duplicate rows
xycoords <- as.matrix(xycoords)
rownames(xycoords) <- studyDesign$Coords # Rownames as locations unique ID

rL.spatial = HmscRandomLevel(sData = xycoords) # Spatial random effect to account for space induced patterns
XData <- XData[,grep(paste(c(var_interest,"Depth"),collapse="|"),colnames(XData))]


# Use the Hmsc model constructor to define a model
# note that in the random effects the left-hand sides in the list (year, location, sample)
# refer to the columns of the studyDesign

# models for different subsets of the data
XFormula_list <- c(XFormula_lt,XFormula_null) # list with my formulas
names(XFormula_list) <- c("XFormula_lt","XFormula_null")
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

save(m_list, file = file.path(modelDir, "unfitted_models_flt_subset.RData"))

# TESTING THAT MODELS FIT WITHOUT ERRORS
mod_type <- c("mpa","mabu","mclassic","mclassicscale","mhellinger","mhellingerscale")

for(i in 1:length(XFormula_list)){
  for(j in mod_type){
    print(paste0(names(XFormula_list)[i]," : ",j))
    m_list[[names(XFormula_list)[i]]][[j]] <- sampleMcmc(m_list[[i]][[j]],samples=2,verbose = 1)
  }
}


##### SAMPLING ##########


nParallel = 3 #Default: nParallel = nChains
# set.seed(973) #reproductibility
# SET DIRECTORIES 
localDir = "R_scripts/paper_16s/Hmsc_pipeline/"
modelDir = file.path(localDir, "models")
load(file=file.path(modelDir,"unfitted_models_flt_subset.RData"))

samples_list = c(250)
thin_list = c(10)
nChains <- 4


# Classic Fit

mod_type <- c("mpa","mabu","mclassic","mclassicscale","mhellinger","mhellingerscale")
model_formu <- length(m_list)

for(i in 1:length(samples_list)){
  
  thin <- thin_list[i]
  samples <- samples_list[i]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  
  for (m_formu in 1:model_formu) {
    for(m_type in mod_type){
      
      mod_name <- paste0(names(m_list)[m_formu],"_",m_type)
      print(paste0("model -> ",mod_name))
      
      filename <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                           "_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           "_",mod_name,
                                           "_subset.Rdata",
                                           sep = ""))
      if(file.exists(filename)){
        print("model had been fitted already")
      } else {
        print(date())
        m <- m_list[[m_formu]][[m_type]]
        m <- sampleMcmc(m, 
                            samples = samples, 
                            thin=thin,
                            adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                            transient = ceiling(0.5*samples*thin),
                            nChains = nChains,
                            updater = list(GammaEta = FALSE),
                            verbose = 100) 
        m_list[[m_formu]][[m_type]] <- m
        gc()
        save(m,file=filename)
      }
    }
  }
}

# FIT WITH HMSC HPC


if(is.null(nParallel)) nParallel <- nChains



# Initiate sampling for HPC 
# Produced INIT files can be read into R and look like classic R hmsc output.
# They are just 'primed' fit for HPC to use.
for(i in 1:length(samples_list)){
  
  thin <- thin_list[i]
  samples <- samples_list[i]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  
  for (m_formu in 1:model_formu) {
    for(m_type in mod_type){
      
      mod_name <- paste0(names(m_list)[m_formu],"_",m_type)
      print(paste0("model -> ",mod_name))
      
      filename <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                           "_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           "_",mod_name,
                                           "_init_files.rds",
                                           sep = ""))
      if(file.exists(filename)){
        print("model had been fitted already")
      } else {
        print(date())
        m <- m_list[[m_formu]][[m_type]]
        minit <- sampleMcmc(m, 
                            samples = samples, 
                            thin=thin,
                            adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                            transient = ceiling(0.5*samples*thin),
                            nChains = nChains,
                            updater = list(GammaEta = FALSE),
                            engine = "HPC") 
        saveRDS(to_json(minit), file=filename)
        gc()
      }
    }
  }
}


# RUn with Pytorch HPS version.
python <- file.path("R_scripts/paper_16s/Hmsc_pipeline/hmsc-venv/", "Scripts", "python") 

verbose = 100
for(i in 1:length(samples_list)){
  thin <- thin_list[i]
  samples <- samples_list[i]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  
  for (m_formu in 1:model_formu) {
    for(m_type in mod_type){
      
      mod_name <- paste0(names(m_list)[m_formu],"_",m_type)
      print(paste0("model -> ",mod_name))
      
      filename <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                           "_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           "_",mod_name,
                                           "_init_files.rds",
                                           sep = ""))
      
      post_file_path <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                                 "_samples_", as.character(samples),
                                                 "_chains_",as.character(nChains),
                                                 "_",mod_name,
                                                 "_post_file.rds",
                                                 sep = "")) 
      
      if(file.exists(post_file_path)){
        print("model had been fitted already")
      } else {
        
        
        python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                                "--input", shQuote(filename),
                                "--output", shQuote(post_file_path),
                                "--samples", samples,
                                "--transient",  ceiling(0.5*samples*thin),
                                "--thin", thin,
                                "--verbose", verbose)
        cat(paste(shQuote(python), python_cmd_args), "\n")
        system2(python, python_cmd_args)
        gc()
      }
    }
  }
}


#### EVALUATION #####

localDir <- "R_scripts/paper_16s/Hmsc_pipeline/"
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
post_file_path <- list.files(modelDir,full.names = T)
post_file_path <- post_file_path[-c(13:15)]

nChains <- 4
importFromHPC <- NULL
postList <- NULL
fitTF <- NULL
varparts <- NULL
predcomputed <- NULL
predcomputedexp <- NULL
predcomputedreal <- NULL
nSamples = c(250)
thin = c(1000)
transient = ceiling(0.5*nSamples*thin)
nChains <- 4

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


  for(i in 1:length(post_file_path)){
    
    nameOf <- stringr::str_extract(post_file_path[i],"(?<=_)m[a-z]+(?=_)")
    if(grepl("null",post_file_path[i])){
      formu <- "XFormula_null"
    }else{
      formu <- "XFormula_lt"
    }
    fitTF[[formu]][[nameOf]] <- loadRData(post_file_path[[i]])
    
    varparts[[formu]][[nameOf]] <- computeVariancePartitioning(fitTF[[formu]][[nameOf]])
    predcomputedexp[[formu]][[nameOf]] <- computePredictedValues(fitTF[[formu]][[nameOf]],expected = T)
    predcomputedreal[[formu]][[nameOf]] <- computePredictedValues(fitTF[[formu]][[nameOf]],expected = F)
  }


varpart_list <- NULL
for(formu in names(varparts)){
  for(i in mod_type){
    
    varpart_list[[formu]][[i]] <- varparts[[formu]][[i]]$vals%>%
      as.data.frame()%>%
      mutate(var=rownames(.))%>%
      reshape2::melt()%>%
      ggplot(aes(x=variable,y=value,fill=var))+
      geom_bar(stat = "identity",position = "stack",width = 1)+
      scale_fill_manual(values=RColorBrewer::brewer.pal(12,"Paired"),name="Variables")+
      xlab("Clusters")+
      ylab("Variance partition (%)")+
      ggpubr::theme_classic2()+
      scale_y_continuous(expand = c(0,0))+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      ggtitle(i)
  }
}

varpart_lt <- varpart_list$XFormula_lt$mabu+varpart_list$XFormula_lt$mpa+varpart_list$XFormula_lt$mclassic+varpart_list$XFormula_lt$mclassicscale+varpart_list$XFormula_lt$mhellinger+varpart_list$XFormula_lt$mhellingerscale+plot_layout(guides="collect")
varpart_null <- varpart_list$XFormula_null$mabu+varpart_list$XFormula_null$mpa+varpart_list$XFormula_null$mclassic+varpart_list$XFormula_null$mclassicscale+varpart_list$XFormula_null$mhellinger+varpart_list$XFormula_null$mhellingerscale+plot_layout(guides="collect")
ggsave("R_scripts/paper_16s/Hmsc_pipeline/results/subset/varpart_lt.png",varpart_lt,width = 13,height = 7)
ggsave("R_scripts/paper_16s/Hmsc_pipeline/results/subset/varpart_null.png",varpart_null,width = 13,height = 7)




nm <- length(fitTF)

posterior_list <- NULL
post_tmp <- NULL

for(formu in names(fitTF)){
  for(i in mod_type){
    nr <- fitTF[[formu]][[i]]$nr
    mpost <- convertToCodaObject(fitTF[[formu]][[i]])
    
    for(p in 1:length(mpost)){
      
      if(length(mpost[[p]])!=nr){
        valuez <- mpost[[p]] 
      }else{
        valuez <- mpost[[p]][[1]]
      }
      
      post_tmp <- NULL
      for (c in 1:nChains){
        print(paste0(formu," - ",i,": ",names(mpost)[p]," - chain: ",c))
        tmp <- t(valuez[[c]])%>%
          as.data.frame()%>%
          mutate(postsamp = rownames(.))%>%
          tidyr::separate(col=postsamp,
                          into=c("var","cluster"),
                          sep=",")%>%
          mutate(var=gsub("[A-Z]\\[|\\(|\\)","",var),
                 cluster=sub("]","",cluster))%>%
          reshape2::melt(id.vars=c("var","cluster"))%>%
          mutate(param=names(mpost)[p])
        
        if(names(mpost)[p]=="Omega"){
          rm_sample <- sample(unique(tmp$cluster),100)
          tmp%<>%
            filter(cluster%in%rm_sample)
        }
        post_tmp <- rbind(post_tmp,tmp)
        
      }
      posterior_list[[formu]][[i]][[names(mpost)[p]]] <- post_tmp
      gc()
    }
    
  }
}



EvModFit <- NULL
for(formu in names(fitTF)){
  for(i in mod_type){
    print(paste0(formu,": ",i))
    EvModFit[[formu]][[i]] <- evaluateModelFit(fitTF[[formu]][[i]],predcomputed[[formu]][[i]])
  }
}


melted_modfit <- NULL
for(formu in names(fitTF)){

  for(i in mod_type){
    print(paste0(formu,": ",i))
    binded_melted_i <- EvModFit[[formu]][[i]]%>%
      bind_cols(.id = i)%>%
      reshape2::melt()
    binded_melted_i$sp <- fitTF[[formu]][[i]]$spNames
    melted_modfit[[formu]] <- rbind(melted_modfit[[formu]],binded_melted_i)
  }
}

vio_modfit <- NULL
for(formu in names(fitTF)){
  for(i in mod_type){
    print(paste0(formu,": ",i)) 
    vio_modfit[[formu]][[i]] <- 
      melted_modfit[[formu]]%>%
      filter(.id==i)%>%
      ggplot(aes(x=1,y=value))+
      geom_violin(fill="darkcyan")+
      facet_wrap(~variable, scales = "free_y")+
      ggtitle(i)
    ggsave(paste0(resultDir,"/subset/modfit/modfit_",formu,"_",i,".png"),vio_modfit[[formu]][[i]],height = 6,width = 5)
    
  }
}


gc()

caterpillar_list <- NULL
density_list <- NULL
for(formu in names(fitTF)){
  for(i in mod_type){
    print(paste0(formu,": ",i)) 
    nr <- fitTF[[formu]][[i]]$nr
  mpost <- convertToCodaObject(fitTF[[formu]][[i]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
  for(p in 1:length(names(mpost))){
    param <- names(mpost)[p]
    print(param)
    if(length(mpost[[param]])!=nr){
      dim_param <- dim(mpost[[param]][[1]]) 
    }else{
      dim_param <- dim(mpost[[param]][[1]][[1]])
    }
    
    tab_pooledchain <- NULL
    for (c in 1:nChains){
      
      if(length(mpost[[param]])!=nr){
        poolchain <- as.data.frame(mpost[[param]][[c]])
        
      }else{
        poolchain <- as.data.frame(mpost[[param]][[1]][[c]])
      }
      
      poolchain$chain <- paste0("Chain ",c)
      poolchain$X <- 1:samples_list
      tab_pooledchain <- rbind(tab_pooledchain,poolchain)
    }
    
    meltedframe <- reshape2::melt(tab_pooledchain,id.vars = c("chain","X"))
    rm_sample <- sample(unique(meltedframe$variable),ifelse(dim_param[2]>50,42,dim_param[2]/2))
    
    caterpillar_list[[formu]][[i]][[names(mpost)[p]]] <- 
      meltedframe %>% 
      filter(variable%in%rm_sample) %>%
      droplevels()%>% 
      mutate(chain=as.factor(chain)) %>% 
      ggplot(aes(x=X,y=value,color=chain))+
      geom_line(alpha=.7)+
      xlab("")+
      facet_wrap(~variable,scales = "free")+
      scale_color_manual(values=c("darkgreen","darkorange","cadetblue","coral1"))+
      ggtitle(names(mpost)[p])
    
    density_list[[formu]][[i]][[names(mpost)[p]]] <- 
      meltedframe %>% 
      filter(variable%in%rm_sample) %>%
      droplevels()%>% 
      mutate(chain=as.factor(chain)) %>% 
      ggplot(aes(x=value,color=chain))+
      geom_density()+
      xlab("")+
      facet_wrap(~variable,scales = "free")+
      scale_color_manual(values=c("darkgreen","darkorange","cadetblue","coral1"))+
      ggtitle(names(mpost)[p])
  }
}
}
  
  
for(formu in names(fitTF)){
  for(i in mod_type){
    for(j in 1:length(caterpillar_list[[formu]][[i]])){
      ggsave(paste0(resultDir,"/subset/caterpillar/","caterpillar_",formu,"_",i,"_",names(caterpillar_list[[formu]][[i]][j]),".png"),caterpillar_list[[formu]][[i]][[j]],width = 13,height = 7)
      ggsave(paste0(resultDir,"/subset/density/","density_",formu,"_",i,"_",names(density_list[[formu]][[i]][j]),".png"),density_list[[formu]][[i]][[j]],width = 13,height = 7)
      
    }
  }
}

gel_list <- NULL
psrf_list <- NULL
geltmp <- NULL
for(formu in names(fitTF)){
  for(i in mod_type){
    print(paste0(formu,": ",i)) 
    mpost <- convertToCodaObject(fitTF[[formu]][[i]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    nr <- fitTF[[formu]][[i]]$nr
    gc()
    for(p in 1:length(names(mpost))){
      param <- names(mpost)[p]
      print(param)
      if(length(mpost[[param]])!=nr){
        geltmp <- gelman.diag(mpost[[param]],multivariate = F)
        save(geltmp,file = paste0(resultDir,"/subset/gelman/gelman_",formu,"_",i,"_",param,".RData"))
        gc()
      }else{
        if(param=="Omega"){
          for(r in 1:nr){
            geltmp <- gelman.diag(mpost[[param]][[r]],multivariate = F)
            save(geltmp,file = paste0(resultDir,"/subset/gelman/gelman_",formu,"_",i,"_",param,"_",r,".RData"))
            gc()
          }
        }
        else{
          for(r in 1:nr){
            geltmp <- gelman.diag(mpost[[param]][[r]],multivariate = F)
            save(geltmp,file = paste0(resultDir,"/subset/gelman/gelman_",formu,"_",i,"_",param,"_",r,".RData"))
            gc()
          }
        }
      }
    }
  }
}



allgelman <- list.files(paste0(resultDir,"/subset/gelman/"),full.names = T)
psrf_plot <- NULL
for(formu in names(fitTF)){
  for(i in mod_type){
    print(paste0(formu,": ",i)) 
    psrf_data_plot <- NULL
    gelman_mod <- allgelman[grepl(paste0(formu,"_",i,"_"),allgelman)]
    for(p in 1:length(gelman_mod)){
      gelmantmp <- loadRData(gelman_mod[p])
      psrftmp <- gelmantmp$psrf
      param <- stringr::str_extract(gelman_mod[p],"(?<=gelman_).+(?=\\.)")
      dtmp <- data.frame(psrf=psrftmp[,1],
                         param=param)
      psrf_data_plot[[p]] <- dtmp
    }
    data_plot <-data.table::rbindlist(psrf_data_plot,use.names = T)
    psrf_plot[[formu]][[i]] <- data_plot%>%
      mutate(param=as.factor(param))%>%
      ggplot(aes(x=param,y=psrf,fill=param))+
      geom_violin()+
      facet_wrap(~param,scales = 'free',drop = T)
    ggsave(paste0(resultDir,"/subset/gelman/plots/gelmanplot_",formu,"_",i,".png"),plot=psrf_plot[[formu]][[i]],height = 7,width=10)
  }
}

supportLevel <- .95
omega_cor <- computeAssociations(fitTF$XFormula_lt$mabu)
toplot <- ((omega_cor[[1]]$support > supportLevel)
           + (omega_cor[[1]]$support < (1-supportLevel))
           > 0)* omega_cor[[1]]$mean

toplot%>%
  reshape2::melt()%>%
  ggplot(aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  scale_fill_gradient2(low="cornflowerblue",mid="white",high='darkorange3')

biplot_list <- NULL
factors <- c(1,2)

for(formu in names(fitTF)){
  for(i in mod_type){
  etaPost <- getPostEstimate(fitTF[[formu]][[i]],"Eta")
  lambdaPost <-getPostEstimate(fitTF[[formu]][[i]],"Lambda")
  
  scale1 <- abs(c(min(etaPost$mean[, factors[1]]), max(etaPost$mean[, 
                                                                    factors[1]])))
  scale2 <- abs(c(min(etaPost$mean[, factors[2]]), max(etaPost$mean[, 
                                                                    factors[2]])))
  scale1 <- min(scale1/abs(c(min(lambdaPost$mean[factors[1], 
  ]), max(lambdaPost$mean[factors[1], ]))))
  
  scale2 <- min(scale2/abs(c(min(lambdaPost$mean[factors[2], 
  ]), max(lambdaPost$mean[factors[2], ]))))
  
  scale <- min(scale1, scale2)
  
  lv_frame <- data.frame(LV1=c(etaPost$mean[, factors[1]],scale * lambdaPost$mean[factors[1], ]),
                         LV2=c(etaPost$mean[, factors[2]],scale * lambdaPost$mean[factors[2], ]),
                         scores=c(rep("eta_sites",fitTF[[formu]][[i]]$ny),
                                  rep("lambda_species",fitTF[[formu]][[i]]$ns)),
                         names=c(rownames(fitTF[[formu]][[i]]$X),fitTF[[formu]][[i]]$spNames))
 
  lv_frame%<>%
    left_join(mutate(as.data.frame(fitTF[[formu]][[i]]$X),
                     names=rownames(fitTF[[formu]][[i]]$X),
                     Depth=XData$Depth))
  
  name_var <- ifelse(grepl("null",formu),"Sequencing depth", "Min temperature (C°)")
  
  if(grepl("null",formu)){
    p <- lv_frame%>%
      ggplot(aes(x=LV1,y=LV2,shape=scores,color=`XData$Depth`,alpha=scores))
  }else{
    p <- lv_frame%>%
      ggplot(aes(x=LV1,y=LV2,shape=scores,color=`XData$tmmn_lt`,alpha=scores))
  }
  
  biplot_list[[formu]][[i]] <- p+
    geom_point()+
    scale_alpha_manual(values = c(1,.4),name="scores")+
    scale_shape_manual(values = c(16,17),name="scores")+
    scale_color_gradient2(low="yellow",
                          mid = "coral1",
                          high="darkmagenta",
                          midpoint = mean(ifelse(grepl("null",formu),lv_frame$`XData$Depth`,lv_frame$`XData$tmmn_lt`),na.rm = T),
                          na.value="aquamarine4",
                          name=name_var)+
    ggpubr::theme_classic2()+
    labs(x="Latent variable 1",y="Latent variable 2")+
    ggtitle(paste0(formu,": ",i))
  ggsave(paste0("R_scripts/paper_16s/Hmsc_pipeline/results/subset/biplots/biplot_",formu,"_",i,".png"),biplot_list[[formu]][[i]],width=6,height = 5)
}
}


### Diversity prediction

# yijprime <- sqrt(417/5000)
# (yijprime)^2*5000
# 
# yijprime <- sqrt(1/5000)
# (yijprime)^2*5000
# 
# 
# sqrt(1371/25139) #0.2333
alpha_pred_plot <- NULL
for(formu in names(fitTF)){
  for(i in mod_type){
    
    if(grepl("hellinger",i)){
      
      melted.df <- as.data.frame(predcomputedreal[[formu]][[i]])
      melted.df <-  mutate(melted.df,site=rownames(melted.df))
      # dplyr::select(c(site,grep("\\.1{1}$|\\.2{1}$|\\.3{1}$",colnames(.))))%>%
      melted.df <- reshape2::melt(melted.df)
      melted.df$cluster <-  lapply(strsplit(as.character(melted.df$variable),split = "\\."),"[",1)
      melted.df$mcmc_sample <-  lapply(strsplit(as.character(melted.df$variable),split = "\\."),"[",2)
      melted.df <- melted.df[,c(1,3,4,5)]
      melted.df <- pivot_wider(melted.df,names_from=cluster,values_from=value)
      melted.df <- cbind(data.frame(obs_alpha=rep(rowSums(fitTF[[formu]]$mclassicscale$Y),1000)),melted.df) 
      colsid <- colnames(melted.df)[grep("Cluster",colnames(melted.df))]
      df2use <- melted.df%>%
        group_by(mcmc_sample)%>%
        mutate(across(.cols=starts_with("Cluster"),.fns=~ifelse((.x^2)*obs_alpha<0.5,
                                                                0,
                                                                round((.x^2)*obs_alpha,digits=0))))%>%
        ungroup()%>%
        reshape2::melt(id.vars=c("obs_alpha","site","mcmc_sample"))
      names(df2use)[which(names(df2use)=="variable")] <- "cluster"
      rm(melted.df)
      gc()
    }else{
      df2use <- predcomputedreal[[formu]][[i]]%>%
        reshape2::melt()%>%
        rename(mcmc_sample=Var3,
               site=Var1,
               cluster=Var2)%>%
        mutate(value=ifelse(value<0.5,0,round(value,digits=0)))
    }
    
    mod_alpha_pred <- df2use%>%
      group_by(site,mcmc_sample)%>%
      summarise(nzeroz=sum(value==0),
                pred_alphadiv=sum(value>0))%>%
      ungroup()%>%
      mutate(input="predicted")
    
    
    
    # 
    # tset <- data.frame(alphadiv=rowSums(fitTF$mclassicscale$Y>0))%>%
    #   mutate(site=rownames(.),
    #          input="oberved")%>%
    #   full_join(mod_alpha_pred)
    # # tset should be 199199 row long
    # 
    
    minmax <- mod_alpha_pred%>%
      group_by(site)%>%
      summarise(minpred=min(pred_alphadiv),
                maxpred=max(pred_alphadiv))
    
    df2plot <- data.frame(alphadiv=rowSums(fitTF[[formu]]$mclassic$Y>0))%>%
      mutate(site=rownames(.),
             input="oberved")%>%
      full_join(minmax)
    
    alpha_pred_plot[[formu]][[i]] <- df2plot%>%
      arrange(desc(alphadiv))%>%
      mutate(site=forcats::fct_inorder(site))%>%
      ggplot()+
      geom_segment(aes(x=site,xend=site,y=minpred,yend=maxpred),color="black")+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90))+
      ylab("ASVs ricness")+
      geom_point(aes(x=site,y=alphadiv),fill="darkorange",size=2,shape=23,color="coral1")+
      ggtitle(i)
    ggsave(paste0("R_scripts/paper_16s/Hmsc_pipeline/results/subset/predalpha/predalpha_",i,".png"),alpha_pred_plot[[formu]][[i]],width=8,height = 5)
    
    
    
    rm(df2use,mod_alpha_pred,minmax)
    gc()
    
    
  }
}

df2use%>%
  group_by(mcmc_sample,site)%>%
  mutate(relab_pred=(value/sum(value))*100,
         relab_obs=(value_obs/sum(value_obs))*100)%>%
  ungroup()%>%
  mutate(delta_relab=abs(relab_pred-relab_obs))%>%
  group_by(site)%>%
  summarise(mean_delta=mean(delta_relab),
            sd_delta=sd(delta_relab))











test_trans2 <- 
  predcomputedreal$XFormula_lt$mabu%>%
  data.frame()%>%
  reshape2::melt()%>%
  ggplot(aes(x=value))+
  geom_histogram(bins = 100)+ 
  ggtitle("pred_mabu")+
  predcomputedexp$XFormula_lt$mhellinger%>%
  data.frame()%>%
  reshape2::melt()%>%
  ggplot(aes(x=value))+
  xlim(c(0,4))+
  geom_histogram(bins = 100)+ 
  ggtitle("pred_hell")+
  Ysub%>%
  data.frame()%>%
  mutate(across(everything(),.fns=~ifelse(.x==0,NA,.x)))%>%
  mutate(across(everything(),.fns=~log(.x)))%>%
  mutate(across(everything(),.fns=~scale(.x)))%>%
  reshape2::melt()%>%
  ggplot(aes(x=value))+
  geom_histogram(bins = 100)+
  ggtitle("Non null logged and scaled reads*sample distribution")







m <- fitTF$XFormula_lt$mabu
m$XFormula <-as.formula(paste(gsub("XData\\$","",m$XFormula),collapse = ""))
gdient <- constructGradient(m,focalVariable = "tmmn_lt")
predY <- predict(m,Gradient = gdient,expected = T)
plotGradient(m,gdient,pred=predY,measure = 'Y',index=20,showData=T)



m <- fitTF$XFormula_lt$mhellinger
m$XFormula <-as.formula(paste(gsub("XData\\$","",m$XFormula),collapse = ""))
gdient <- constructGradient(m,focalVariable = "tmmn_lt")
predY <- predict(m,Gradient = gdient,expected = T)
plotGradient(m,gdient,pred=predY,measure = 'Y',index=20,showData=T)
