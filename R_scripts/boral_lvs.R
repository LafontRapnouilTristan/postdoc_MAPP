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
library(microeco)

bioclim_predictors <- c("siteCol","X","Y","aet_lt","bulk","def_lt","evi_lt","gHM","gpp_lt","height",'lai_lt',"moist_lt","pdsi_lt","ph","pr_lt","ssm_lt","susm_lt","swe_lt","tmmn_lt","tmmx_lt","vpd_lt")

gdm_formated <- NULL
gdm_fitted <- NULL
for(i in c("16s","18s","23s")){
    
    meco_tmp <- get(paste0("meco_",i)) # get com datasets
    data_sat_tmp <- data_sat
    
    comunity_df_hellingered <- as.data.frame(labdsv::hellinger(t(meco_tmp$otu_table))) # extract community matrices and transform to hellinger
    
    comunity_df_hellingered$siteCol <-  gsub("[A-Z]|-|(?=_).*","",rownames(comunity_df_hellingered),perl = T)  # get site name (i.e., sample)
    comunity_df_hellingered %<>%  
        mutate(siteCol=ifelse(nchar(siteCol)<3,paste0("0",siteCol),siteCol))%>% # reformat 1 into 001 and 64 to 064 etc...
        mutate(siteCol=ifelse(nchar(siteCol)<3,paste0("0",siteCol),siteCol))
    data_sat_tmp %<>% rename(siteCol=Sample) # rename 
    
    gdm_biodata <- left_join(comunity_df_hellingered,data_sat_tmp[,c(3:5)]) %>% # create gdm_biodata with site X and Y (long and lat)
        relocate(siteCol,X,Y) %>% # relocate columns
        filter(!siteCol%in%c("116","119","166","042047","047042")) # remove problematic samples
    
    gdm_preddata <- data_sat_tmp %>%
        filter(siteCol%in%gdm_biodata$siteCol) %>% #remove sites not in biodata
        select(bioclim_predictors) # get only variables we are interested in
    
    gdm_formated_tmp <- gdm::formatsitepair(bioData = gdm_biodata, # create data in the gdm format
                                            bioFormat = 1,
                                            dist="horn",
                                            abundance = T,
                                            siteColumn = "siteCol",
                                            XColumn = "X",
                                            YColumn = "Y",
                                            predData = gdm_preddata)
    gdm_formated[[i]] <- gdm_formated_tmp # store formated data
    gdm_fitted[[i]] <- gdm::gdm(gdm_formated_tmp, # fit gdm
                                geo = T)
}
View(gdm_fitted$"16s")
plot(gdm_fitted$"18s")
plot(gdm_fitted$"23s")
plot(gdm_fitted$"16s")
View(gdm::plot.gdm)
# gdm with atlasr https://github.com/jiho/atlasr
save.image(file="ressource//env_gdm_testing.Rdata")
