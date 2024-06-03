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


#################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
nParallel = 4 #Default: nParallel = nChains
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





##################################################################################################
# FIT WITH BASIC R HMSC
##################################################################################################
nm <- length(m_list)
samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nChains = 4
if(is.null(nParallel)) nParallel = nChains
Lst <- 1
while(Lst <= length(samples_list)){
    thin <- thin_list[Lst]
    samples <- samples_list[Lst]
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    filename <- file.path(modelDir,paste("model_thin_",as.character(thin),
                                         "_samples_",as.character(samples),
                                         '_chaines_',as.character(nChains),
                                         ".Rdata",sep = ""))
    if(file.exists(filename)){
        print("Model had been fitted already")
    } else {
      print(date())
      for(mod in 1:nm){
          print(paste0("Model : ",names(m_list)[mod]))
          m <- m_list[[mod]]
          m <- sampleMcmc(m,
                          samples = samples,
                          thin = thin,
                          adaptNf = rep(ceiling(0.4*samples*thin),m$nr),
                          transient = ceiling(0.5*samples*thin),
                          nChains = nChains,
                          nParallel = nParallel,
                          updater = list(GammaEta = F),
                          verbose = 100)
          m_list[[mod]] <- m
      }
        save(m_list,file = filename)
    }
    Lst <- Lst+1
}




##################################################################################################
# FIT WITH HMSC HPC
##################################################################################################

library(Hmsc)

nm <- length(m_list)
samples_list <- c(5)
thin_list <- c(1)
nChains <- 4
if(is.null(nParallel)) nParallel <- nChains

for(i in 1:length(samples_list)){
    set.seed(973)
    thin <- thin_list[i]
    samples <- samples_list[i]
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))


        for (mi in 1:nm) {
            print(paste0("model = ",names(m_list)[mi]))
            
            filename <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                                 "_samples_", as.character(samples),
                                                 "_chains_",as.character(nChains),
                                                 "_",names(m_list)[mi],
                                                 "_init_files.rds",
                                                 sep = ""))
            if(file.exists(filename)){
                print("model had been fitted already")
            } else {
                print(date())
            m <- m_list[[mi]]
            minit <- sampleMcmc(m, 
                            samples = samples, 
                            thin=thin,
                            adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                            transient = ceiling(0.5*samples*thin),
                            nChains = nChains,
                            updater = list(GammaEta = FALSE),
                            engine = "HPC") 
            saveRDS(to_json(minit), file=filename)
        }
    }
}

python <- file.path("./R_scripts/paper_16s/Hmsc_pipeline/hmsc-venv", "Scripts", "python") 

verbose = 100
for(i in 1:length(samples_list)){
    set.seed(973)
    thin <- thin_list[i]
    samples <- samples_list[i]
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    
    
    for (mi in 1:nm) {
        print(paste0("model = ",names(m_list)[mi]))
        
        filename <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             "_",names(m_list)[mi],
                                             "_init_files.rds",
                                             sep = ""))
        
        post_file_path <- file.path(modelDir,paste("models_thin_", as.character(thin),
                                                   "_samples_", as.character(samples),
                                                   "_chains_",as.character(nChains),
                                                   "_",names(m_list)[mi],
                                                   "_post_file.rds",
                                                   sep = "")) 
        
        python_cmd_args = paste("-m hmsc.run_gibbs_sampler",
                                "--input", shQuote(filename),
                                "--output", shQuote(post_file_path),
                                "--samples", samples,
                                "--transient",  ceiling(0.5*samples*thin),
                                "--thin", thin,
                                "--verbose", verbose)
        cat(paste(shQuote(python), python_cmd_args), "\n")
        system2(python, python_cmd_args)
    }
}
