nParallel = 10 #Default: nParallel = nChains
set.seed(973) #reproductibility
# SET DIRECTORIES 
localDir = "R_scripts/paper_16s/Hmsc_pipeline/"
modelDir = file.path(localDir, "models")
load(file=file.path(modelDir,"unfitted_models_V1.RData"))


# FIT WITH HMSC HPC

library(Hmsc)
library(jsonify)


samples_list = c(250)
thin_list = c(1000,10000)
nChains <- 4
if(is.null(nParallel)) nParallel <- nChains

mod_type <- c("mpa","mabu","mclassic","mclassicscale","mhellinger","mhellingerscale")
model_formu <- length(m_list)

# Initiate sampling for HPC 
# Produced INIT files can be read into R and look like classic R hmsc output.
# They are just 'primed' fit for HPC to use.
for(i in 1:length(samples_list)){
    
    thin <- thin_list[i]
    samples <- samples_list[i]
    print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
    
    
    for (m_formu in 1:model_formu) {
        for(m_type in mod_type){
            
            mod_name <- paste0(names(m_list)[m_formu]," : ",m_type)
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
            
            mod_name <- paste0(names(m_list)[m_formu]," : ",m_type)
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
