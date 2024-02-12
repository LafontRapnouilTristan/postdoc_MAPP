########## BORAL LVM ################

#jSDM
#HSMC
#MIXOMICS

load("output/microeco_obj.Rdata")
library(boral)
# https://github.com/Zanne-Lab/fungal_wood_endophytes/blob/d2ceab5cf14ff52894cc69f0aa9732f52940add8/code/boral_modfit.R#L4  INSPIRATION

# Setting up

set_mcmc_controls <- function(seed, runType){
    
    if(runType == "test"){
        n.burnin = 50
        n.iteration = 100
        n.thin = 1
    }
    
    if(runType == "real"){
        n.burnin = 100000 
        n.iteration = 200000
        n.thin = 10
    }
    
    mcmc.controls <- list(n.burnin = n.burnin,
                          n.iteration = n.iteration,
                          n.thin = n.thin,
                          seed = seed)
    
    return(mcmc.controls)
}




prior_controls <- function(){
    
    prior.controls = list(type = c("normal","normal","normal","uniform"),
                          hypparams = c(100, #OTU-specific intercepts and row effects
                                        20,  #lv coefs
                                        100, #X coefs
                                        20)) #dispersion
}

lv_controls <- function(){lv.control = list(num.lv = 2, 
                                            type = "independent", 
                                            distmat = NULL)
}
Y <- meco_23s$otu_table


# Model fit

fit <- boral(y = Y,
             family = "negative.binomial",
             lv.control = lv_controls(),
             row.eff = "random",
             save.model=T,
             mcmc.control = set_mcmc_controls(973,"test"),
             prior.control = prior_controls())


lvsplot(fit)

devtools::install_github("mbedward/ggboral")
library(ggboral)
gg_lvsplot(fit,include="objects")
gg_lvsplot(fit,include="attributes")


#### jSDM ########

# install.packages("jSDM")
library(jSDM)


mod_23s_jSDM_poisson <- jSDM_poisson_log(
    # Chains
    burnin=200, mcmc=1000, thin=1,
    # Response variable 
    count_data = t(Y[,!grepl("166|172|174|042|047",colnames(Y))]), 
    # Explanatory variable,
    site_data = na.omit(meco_23s$sample_table[c("X","Y")]),
    site_formula = ~.,
    # Model specification 
    n_latent=2, 
    site_effect="random",
    # Starting values
    alpha_start=0, beta_start=0,
    lambda_start=0, W_start=0,
    V_alpha=1, 
    # Priors
    shape_Valpha=0.1,
    rate_Valpha=0.1,
    mu_beta=0, V_beta=1,
    mu_lambda=0, V_lambda=1,
    # Various 
    seed=974, verbose=1)


### GDM #######

install.packages("gdm")
devtools::install_github('skiptoniam/bbgdm')
