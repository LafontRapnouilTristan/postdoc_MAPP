library(Hmsc)
library(jsonify)
library(magrittr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
localDir <- "R_scripts/paper_16s/Hmsc_pipeline/"
modelDir <- file.path(localDir, "out_models")
resultDir <- file.path(localDir, "results")
post_file_path <- list.files(modelDir,full.names = T)
load(file=file.path(file.path(localDir,"models"),"unfitted_models_flt_V1.RData"))
text.file <- file.path(resultDir,"/MCMC_convergence.txt")

nChains <- 4
importFromHPC <- NULL
postList <- NULL
fitTF <- NULL
varparts <- NULL
predcomputedexp <- NULL
predcomputedpred <- NULL

nSamples = c(250)
thin = c(1000)
transient = ceiling(0.5*nSamples*thin)
nChains <- 4


for(i in 1:length(post_file_path)){
  
  nameOf <- stringr::str_extract(post_file_path[i],"(?<=_)m[a-z]+(?=_)")
  
  importFromHPC[[nameOf]] <- from_json(readRDS(file = post_file_path[[i]])[[1]])
  postList[[nameOf]] <- importFromHPC[[nameOf]][1:nChains]
  cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nameOf]][[nChains+1]]))
  fitTF[[nameOf]] <- importPosteriorFromHPC(m_list$XFormula_lt[[nameOf]], postList[[nameOf]], nSamples, thin, transient)
  varparts[[nameOf]] <- computeVariancePartitioning(fitTF[[nameOf]])
  predcomputedexp[[nameOf]] <- computePredictedValues(fitTF[[nameOf]],expected = TRUE)
  predcomputedreal[[nameOf]] <- computePredictedValues(fitTF[[nameOf]],expected = FALSE)
}


varpart_lt <- NULL
for(i in names(varparts)){
  
  varpart_lt[[i]] <- varparts[[i]]$vals%>%
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

library(patchwork)
varpartz <- varpart_lt$mabu+varpart_lt$mpa+varpart_lt$mclassic+varpart_lt$mclassicscale+varpart_lt$mhellinger+varpart_lt$mhellingerscale+plot_layout(guides="collect")
ggsave("R_scripts/paper_16s/Hmsc_pipeline/results/varpartz.png",varpartz,width = 13,height = 7)

library(coda)



samples_list <- c(250)
thin_list <- c(1000)
nChains <- 4
nst <- length(thin_list)
Lst <- 1


thin = thin_list[Lst]
samples = samples_list[Lst]

models <- fitTF
nm <- length(models)

posterior_list <- NULL
post_tmp <- NULL
for(i in 1:nm){
  nr <- fitTF[[i]]$nr
  mpost <- convertToCodaObject(fitTF[[i]])
  
  for(p in 1:length(mpost)){
    
    if(length(mpost[[p]])!=nr){
      valuez <- mpost[[p]] 
    }else{
      valuez <- mpost[[p]][[1]]
    }
    
    post_tmp <- NULL
    for (c in 1:nChains){
      print(paste0(names(fitTF)[i],": ",names(mpost)[p]," - chain: ",c))
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
      rm_sample <- sample(unique(meltedframe$cluster),411)
      tmp%<>%
        filter(cluster%in%rm_sample)
      }
      post_tmp <- rbind(post_tmp,tmp)
      
    }
    posterior_list[[names(fitTF)[i]]][[names(mpost)[p]]] <- post_tmp
    gc()
  }
  
}




EvModFit <- NULL
for(i in 1:nm){
  print(names(models)[i])
  EvModFit[[names(models)[i]]] <- evaluateModelFit(models[[i]],predcomputed[[i]])
}
 
melted_modfit <- NULL
for(i in 1:nm){
  print(names(models)[i])
  binded_melted_i <-  EvModFit[[names(models)[i]]]%>%
    bind_cols(.id = names(models)[i])%>%
    reshape2::melt()
  binded_melted_i$sp <- models[[i]]$spNames
  melted_modfit <- rbind(melted_modfit,binded_melted_i)
}

vio_modfit <- NULL
for(i in 1:nm){
  print(names(models)[i])
  vio_modfit[[i]] <- 
    melted_modfit%>%
    filter(.id==names(models)[i])%>%
    ggplot(aes(x=1,y=value))+
    geom_violin(fill="darkcyan")+
    facet_wrap(~variable, scales = "free_y")+
    ggtitle(names(models)[i])
  ggsave(paste0(resultDir,"/modfit/modfit_",names(models)[i],".png"),vio_modfit[[i]],height = 6,width = 5)
  
}

gc()
 
caterpillar_list <- NULL
density_list <- NULL
for(j in 1:nm){
  print(names(models)[j])
  nr <- models[[j]]$nr
  mpost <- convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
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
    
    caterpillar_list[[mod_type[j]]][[names(mpost)[p]]] <- 
      meltedframe %>% 
      filter(variable%in%rm_sample) %>%
      droplevels()%>% 
      mutate(chain=as.factor(chain)) %>% 
      ggplot(aes(x=X,y=value,color=chain))+
      geom_line(alpha=.7)+
      xlab("")+
      facet_wrap(~variable,scales = "free_y")+
      scale_color_manual(values=c("darkgreen","darkorange","cadetblue","coral1"))+
      ggtitle(names(mpost)[p])
    
    density_list[[mod_type[j]]][[names(mpost)[p]]] <- 
      meltedframe %>% 
      filter(variable%in%rm_sample) %>%
      droplevels()%>% 
      mutate(chain=as.factor(chain)) %>% 
      ggplot(aes(x=value,color=chain))+
      geom_density()+
      xlab("")+
      facet_wrap(~variable,scales = "free_y")+
      scale_color_manual(values=c("darkgreen","darkorange","cadetblue","coral1"))+
      ggtitle(names(mpost)[p])
  }
}

for(i in 1:length(caterpillar_list)){
  for(j in 1:length(caterpillar_list[[i]])){
    ggsave(paste0(resultDir,"/caterpillar/","caterpillar_",names(caterpillar_list)[i],"_",names(caterpillar_list[[i]][j]),".png"),caterpillar_list[[i]][[j]],width = 13,height = 7)
    ggsave(paste0(resultDir,"/density/","density",names(density_list)[i],"_",names(density_list[[i]][j]),".png"),density_list[[i]][[j]],width = 13,height = 7)
    
  }
}

save(caterpillar_list,density_list,EvModFit,fitTF,importFromHPC,postList,predcomputed,varparts,varpart_lt,vio_modfit,file="somesavesoutHPC.Rdata")
rm(caterpillar_list,density_list,EvModFit,fitTF,importFromHPC,postList,predcomputed,varparts,varpart_lt,vio_modfit)
gc()

gel_list <- NULL
psrf_list <- NULL
geltmp <- NULL
for(j in 1:nm){
  print(names(models)[j])
  mpost <- convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
  nr <- models[[j]]$nr
  gc()
  for(p in 1:length(names(mpost))){
    param <- names(mpost)[p]
    print(param)
    if(length(mpost[[param]])!=nr){
      geltmp <- gelman.diag(mpost[[param]],multivariate = F)
      save(geltmp,file = paste0(resultDir,"/gelman/gelman_",names(models)[j],"_",param,".RData"))
      gc()
    }else{
      if(param=="Omega"){
        for(r in 1:nr){
          geltmp <- gelman.diag(mpost[[param]][[r]][,sample(1:822^2,10000)],multivariate = F)
          save(geltmp,file = paste0(resultDir,"/gelman/gelman_",names(models)[j],"_",param,"_",r,".RData"))
          gc()
        }
      }
      else{
        for(r in 1:nr){
          geltmp <- gelman.diag(mpost[[param]][[r]],multivariate = F)
          save(geltmp,file = paste0(resultDir,"/gelman/gelman_",names(models)[j],"_",param,"_",r,".RData"))
          gc()
        }
      }
    }
  }
}

loadRData <- function(filename){
  load(filename)
  get(ls()[ls() != "filename"])
}


allgelman <- list.files(paste0(resultDir,"/gelman/"),full.names = T)
psrf_plot <- NULL
for(i in 1:nm){
  print(names(models)[i])
  psrf_data_plot <- NULL
  gelman_mod <- allgelman[grepl(paste0(names(models)[i],"_"),allgelman)]
  for(p in 1:length(gelman_mod)){
  gelmantmp <- loadRData(gelman_mod[p])
  psrftmp <- gelmantmp$psrf
  param <- stringr::str_extract(gelman_mod[p],"(?<=gelman_).+(?=\\.)")
  dtmp <- data.frame(psrf=psrftmp[,1],
                     param=param)
  psrf_data_plot[[p]] <- dtmp
  }
  data_plot <-data.table::rbindlist(psrf_data_plot,use.names = T)
  psrf_plot[[names(models)[i]]] <- data_plot%>%
    mutate(param=as.factor(param))%>%
    ggplot(aes(x=param,y=psrf,fill=param))+
    geom_violin()+
    facet_wrap(~param,scales = 'free',drop = T)
  ggsave(paste0(resultDir,"/gelman/plots/gelmanplot_",names(models)[i],".png"),plot=psrf_plot[[names(models)[i]]],height = 7)
}


supportLevel <- .95
omega_cor <- computeAssociations(fitTF$mabu)
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

for(i in 1:nm){

etaPost <- getPostEstimate(fitTF[[i]],"Eta")
lambdaPost <-getPostEstimate(fitTF[[i]],"Lambda")

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
                       scores=c(rep("eta_sites",fitTF[[i]]$ny),
                                rep("lambda_species",fitTF[[i]]$ns)),
                       names=c(rownames(fitTF[[i]]$X),fitTF[[i]]$spNames))
lv_frame%<>%
  left_join(mutate(as.data.frame(fitTF[[i]]$X),names=rownames(fitTF[[i]]$X)))

biplot_list[[names(fitTF)[i]]] <- lv_frame%>%
  ggplot(aes(x=LV1,y=LV2,shape=scores,color=`XData$tmmx_lt`,alpha=scores))+
  geom_point()+
  scale_alpha_manual(values = c(1,.4),name="scores")+
  scale_shape_manual(values = c(16,17),name="scores")+
  scale_color_gradient2(low="yellow",
                       mid = "coral1",
                       high="darkmagenta",
                       midpoint = mean(lv_frame$`XData$tmmx_lt`,na.rm = T),
                       na.value="aquamarine4",
                       name="Max temperature (C°)")+
  ggpubr::theme_classic2()+
  labs(x="Latent variable 1",y="Latent variable 2")+
  ggtitle(names(fitTF)[i])
ggsave(paste0("R_scripts/paper_16s/Hmsc_pipeline/results/biplots/biplot_",names(fitTF)[i],".png"),biplot_list[[names(fitTF)[i]]],width=6,height = 5)
}



### Diversity prediction
load("somesavesoutHPC.Rdata")
rm(varpart_lt,varparts,importFromHPC,EvModFit,density_list,caterpillar_list,postList)
gc()
# yijprime <- sqrt(417/5000)
# (yijprime)^2*5000
# 
# yijprime <- sqrt(1/5000)
# (yijprime)^2*5000
# 
# 
# sqrt(1371/25139) #0.2333

library(dplyr)
library(ggplot2)
library(magritrr)

pred <- computePredictedValues(fitTF$mpa,expected = F)
pred %<>% reshape2::melt()%>%
  rename(mcmc_sample=Var3,
         site=Var1,
         cluster=Var2)

mod_alpha_pred <- pred%>%
  group_by(site,mcmc_sample)%>%
  summarise(nzeroz=sum(value==0),
            pred_alphadiv=sum(value>0))%>%
  ungroup()%>%
  mutate(input="predicted")

minmax <- mod_alpha_pred%>%
  group_by(site)%>%
  summarise(minpred=min(pred_alphadiv),
            maxpred=max(pred_alphadiv))

df2plot <- data.frame(alphadiv=rowSums(fitTF$mpa$Y>0))%>%
  mutate(site=rownames(.),
         input="oberved")%>%
  full_join(minmax)

mpa <- df2plot%>%
  arrange(desc(alphadiv))%>%
  mutate(site=forcats::fct_inorder(site))%>%
  ggplot()+
  geom_segment(aes(x=site,xend=site,y=minpred,yend=maxpred),color="black")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  ylab("ASVs ricness")+
  geom_point(aes(x=site,y=alphadiv),fill="darkorange",size=2,shape=23,color="coral1")+
  ggtitle("mpa")



pred <- computePredictedValues(fitTF$mclassic,expected = F)
pred %<>% reshape2::melt()%>%
  rename(mcmc_sample=Var3,
         site=Var1,
         cluster=Var2)

mod_alpha_pred <- pred%>%
  group_by(site,mcmc_sample)%>%
  summarise(nzeroz=sum(value==0),
            pred_alphadiv=sum(value>0))%>%
  ungroup()%>%
  mutate(input="predicted")

minmax <- mod_alpha_pred%>%
  group_by(site)%>%
  summarise(minpred=min(pred_alphadiv),
            maxpred=max(pred_alphadiv))

df2plot <- data.frame(alphadiv=rowSums(fitTF$mclassic$Y>0))%>%
  mutate(site=rownames(.),
         input="oberved")%>%
  full_join(minmax)

mclassic <- df2plot%>%
  arrange(desc(alphadiv))%>%
  mutate(site=forcats::fct_inorder(site))%>%
  ggplot()+
  geom_segment(aes(x=site,xend=site,y=minpred,yend=maxpred),color="black")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90))+
  ylab("ASVs ricness")+
  geom_point(aes(x=site,y=alphadiv),fill="darkorange",size=2,shape=23,color="coral1")+
  ggtitle("mclassic")


raw_melted <- fitTF$mclassic$Y%>%
  reshape2::melt()%>%
  rename(site=Var1,
         cluster=Var2,
         value_obs=value)

relab <- pred%>%
  left_join(raw_melted,by=c("site","cluster"))%>%
  group_by(mcmc_sample,site)%>%
  mutate(relab_pred=(value/sum(value))*100,
         relab_obs=(value_obs/sum(value_obs))*100)%>%
  ungroup()
  
relab_summarised <- relab%>%
  mutate(delta_relab=abs(relab_pred-relab_obs))%>%
  group_by(site,cluster)%>%
  summarise(mean_delta=mean(delta_relab),
            sd_delta=sd(delta_relab))

my_breaks = c(0, 10, 20, 30, 40, 50)
my_breaks <- plyr::round_any(exp(seq(log(1), log(50), length = 5)), 3)
relab_summarised%>%
  ggplot(aes(x=site,y=cluster,fill=mean_delta))+
  geom_hex()+
  scale_fill_gradient(low="cornflowerblue",high = "darkorange",trans="log10")







alpha_pred_plot <- NULL
for(i in names(predcomputed)){
  
  if(grepl("hellinger",i)){
    
    melted.df <- as.data.frame(predcomputed[[i]])
    melted.df <-  mutate(melted.df,site=rownames(melted.df))
    # dplyr::select(c(site,grep("\\.1{1}$|\\.2{1}$|\\.3{1}$",colnames(.))))%>%
    melted.df <- reshape2::melt(melted.df)
    melted.df$cluster <-  lapply(strsplit(as.character(melted.df$variable),split = "\\."),"[",1)
    melted.df$mcmc_sample <-  lapply(strsplit(as.character(melted.df$variable),split = "\\."),"[",2)
    melted.df <- melted.df[,c(1,3,4,5)]
    melted.df <- pivot_wider(melted.df,names_from=cluster,values_from=value)
    melted.df <- cbind(data.frame(obs_alpha=rep(rowSums(fitTF$mclassicscale$Y),1000)),melted.df) 
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
    df2use <- predcomputed[[i]]%>%
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
  
  df2plot <- data.frame(alphadiv=rowSums(fitTF$mclassic$Y>0))%>%
    mutate(site=rownames(.),
           input="oberved")%>%
    full_join(minmax)
  
  alpha_pred_plot[[i]] <- df2plot%>%
    arrange(desc(alphadiv))%>%
    mutate(site=forcats::fct_inorder(site))%>%
    ggplot()+
    geom_segment(aes(x=site,xend=site,y=minpred,yend=maxpred),color="black")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    ylab("ASVs ricness")+
    geom_point(aes(x=site,y=alphadiv),fill="darkorange",size=2,shape=23,color="coral1")+
    ggtitle(i)
  ggsave(paste0("R_scripts/paper_16s/Hmsc_pipeline/results/predalpha/predalpha_",i,".png"),alpha_pred_plot[[i]],width=8,height = 5)
  
  
  
  rm(df2use,mod_alpha_pred,minmax)
  gc()
  
  
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