# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Unfitted models.

#	OUTPUT. Fitted models, with fitting done for multiple RUNs:
# first short MCMC chains (to provide some results fast), and then with increasingly long MCMC chains
# (until MCMC convergence or computational limit is reached)
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
nParallel = 1 #Default: nParallel = nChains
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "./R_scripts/paper_16s/Hmsc_pipeline/"
modelDir = file.path(localDir, "models")
load(file=file.path(modelDir,"unfitted_models.RData"))
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

library(Hmsc)

nm <- length(m_list)
samples_list <- c(5,250,250,500)
thin_list <- c(1,1,10,50)
nChains <- 4
if(is.null(nParallel)) nParallel <- nChains

for(i in 1:length(samples_list)){
    set.seed(973)
    thin <- thin_list[i]
    samples <- samples_list[i]
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    filename <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
    if(file.exists(filename)){
        print("model had been fitted already")
    } else {
        print(date())
        for (mi in 1:nm) {
            print(paste0("model = ",names(m_list)[mi]))
            m <- m_list[[mi]]
            m <- sampleMcmc(m, 
                            samples = samples, 
                            thin=thin,
                            adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                            transient = ceiling(0.5*samples*thin),
                            nChains = nChains,
                            nParallel = nParallel) 
            m_list[[mi]] <- m
        }
    }
    save(m_list,file=filename)
}
