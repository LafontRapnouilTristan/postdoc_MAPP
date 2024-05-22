color_pool <- c("#7f61d3","#84b83a","#be4eb2","#4db955","#d54272","#6abf87","#d14b33","#3fc1bf","#d99037",
                "#616bba","#b8af3d","#c387d1","#5a8233","#dd81a0","#37835d","#9a4769","#c6a96d","#619dd7")

dxy <- site_data_lt[,1:2] %>% set_rownames(site_data$Sample) %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()
color_xy <- as.numeric(as.factor(site_data$Country))[order.dendrogram(dxy)]

d_lt <- site_data_lt[,grepl("_lt",names(site_data_lt))] %>% set_rownames(site_data$Sample) %>% scale() %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()
color_lt <- as.numeric(as.factor(site_data$Country))[order.dendrogram(d_lt)]

d_st_before <- site_data_st_before_imputation_imputed[,grepl("_st",names(site_data_st_before_imputation_imputed))] %>% set_rownames(site_data$Sample) %>% scale() %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()
color_st_before <- as.numeric(as.factor(site_data$Country))[order.dendrogram(d_st_before)]

d_st_after <- site_data_st_after_imputation[,grepl("_st",names(site_data_st_after_imputation))] %>% set_rownames(site_data$Sample) %>% scale() %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()
color_st_after <- as.numeric(as.factor(site_data$Country))[order.dendrogram(d_st_after)]

dendro_list <- dendlist(dxy %>% 
                            set('labels_col',value = color_pool[color_xy])%>%
                            set("branches_lty", 1)%>%
                            set("branches_lwd",1)%>%
                            set("labels_cex",rep(1,nrow(site_data_lt))),
                        d_lt%>%
                            set('labels_col',value = color_pool[color_lt])%>%
                            set("branches_lty", 1)%>%
                            set("branches_lwd",1)%>%
                            set("labels_cex",rep(1,nrow(site_data_lt))) )

dendro_list  %>%
    untangle(method = "step1side") %>%
    tanglegram(common_subtrees_color_lines = T,
               margin_inner=3,
               main_left = "Lon/Lat distances",
               main_right = "Long term variable distances",
               lwd=.5,
               edge.lwd=.5,
               columns_width=c(5,1,5))
dev.copy(png,"R_scripts/paper_16s/Figures/tanglegram.png",width=14,height=16,units='in',res=800)
dev.off()



dendro_list <- dendlist(dxy %>% 
                            set('labels_col',value = color_pool[color_xy])%>%
                            set("branches_lty", 1)%>%
                            set("branches_lwd",1)%>%
                            set("labels_cex",rep(1,nrow(site_data_lt))),
                        d_st_before%>%
                            set('labels_col',value = color_pool[color_st_before])%>%
                            set("branches_lty", 1)%>%
                            set("branches_lwd",1)%>%
                            set("labels_cex",rep(1,nrow(site_data_st_before_imputation_imputed))) )

dendro_list  %>%
    untangle(method = "step1side") %>%
    tanglegram(common_subtrees_color_lines = T,
               margin_inner=3,
               main_left = "Lon/Lat distances",
               main_right = "Short term variable distances (imputed_before)",
               lwd=.5,
               edge.lwd=.5,
               columns_width=c(5,1,5))
dev.copy(png,"R_scripts/paper_16s/Figures/tanglegram_st_before.png",width=14,height=16,units='in',res=800)
dev.off()




dendro_list <- dendlist(dxy %>% 
                            set('labels_col',value = color_pool[color_xy])%>%
                            set("branches_lty", 1)%>%
                            set("branches_lwd",1)%>%
                            set("labels_cex",rep(1,nrow(site_data_lt))),
                        d_st_after%>%
                            set('labels_col',value = color_pool[color_st_after])%>%
                            set("branches_lty", 1)%>%
                            set("branches_lwd",1)%>%
                            set("labels_cex",rep(1,nrow(site_data_st_after_imputation))) )





dendro_list  %>%
    untangle(method = "step1side") %>%
    tanglegram(common_subtrees_color_lines = T,
               margin_inner=3,
               main_left = "Lon/Lat distances",
               main_right = "Long term variable distances (imputed_after)",
               lwd=.5,
               edge.lwd=.5,
               columns_width=c(5,1,5))
dev.copy(png,"R_scripts/paper_16s/Figures/tanglegram_st_after.png",width=14,height=16,units='in',res=800)
dev.off()




tree_list <- dendro_list  %>%
    untangle(method = "step1side")
tree_list[[1]] %>%
    as.ggdend()%>%
    ggplot(horiz = TRUE, theme = NULL)

tree_list[[2]] %>%
    as.ggdend()%>%
    ggplot(horiz = T) + scale_y_reverse(expand = c(0.2, 0)) 









dendro_list  %>%
    untangle(method = "step1side") %>%
    tanglegram(color_lines=ifelse( 1:(199*2)%in%distinct_edges(dxy,d_lt),"lightgrey","grey17"),
               common_subtrees_color_lines = F,
               margin_inner=1.8,
               main_left = "Lon/Lat distances",
               main_right = "Long term variable distances",
               lwd=.5,
               edge.lwd=.5,
               columns_width=c(5,1,5))




#removed from main script 

# Cluster sites

## On Var

pet and aet are highly correlated -> use aet
vap and vpd same -> use vpd
evi and ndvi -> use evi
fpar and lai -> use fpar
gpp and npp -> gpp
moist and water are suposed to be the same but water is at much finer scale (250m vs 9000)

Lon Lat clustering -> Not making any sens ain't it?
```{r}
site_XY_dend <- imputed_data_lt %>% 
    select(c(X,Y))
rownames(site_XY_dend) <- site_data$Sample
dist_site_XY <- vegdist(site_XY_dend,method="euclidean")
clust_site_XY <- hclust(dist_site_XY,method="ward.D2")
hcd_site_XY <- as.dendrogram(clust_site_XY)
dend_data_site_XY <- ggdendro::dendro_data(hcd_site_XY,type='rectangle')

dend_data_site_XY$labels %<>% rename(Sample=label)
dend_data_site_XY$labels %<>% left_join(site_data[,2:3])

dend_XY <- ggplot(dend_data_site_XY$segments)+
    geom_segment(aes(x=x,y=-sqrt(y),xend=xend,yend=-sqrt(yend)))+
    geom_text(data = dend_data_site_XY$labels,aes(x,-(y)+2,label=Sample,color=Country),size=3,angle=90)+
    theme_void()+
    scale_color_manual(values=c("#7f61d3",
                                "#84b83a",
                                "#be4eb2",
                                "#4db955",
                                "#d54272",
                                "#6abf87",
                                "#d14b33",
                                "#3fc1bf",
                                "#d99037",
                                "#616bba",
                                "#b8af3d",
                                "#c387d1",
                                "#5a8233",
                                "#dd81a0",
                                "#37835d",
                                "#9a4769",
                                "#c6a96d",
                                "#619dd7"))+
    ggtitle("Lon/Lat clustering")
```

```{r}
color_pool <- c("#7f61d3","#84b83a","#be4eb2","#4db955","#d54272","#6abf87","#d14b33","#3fc1bf","#d99037",
                "#616bba","#b8af3d","#c387d1","#5a8233","#dd81a0","#37835d","#9a4769","#c6a96d","#619dd7")

d_xy <- site_data_lt[,1:2] %>% set_rownames(site_data$Sample) %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()

d_lt <- imputed_data_lt[,grepl("_lt",names(imputed_data_lt))] %>% select(!grep("pet|vap|ndvi|lai|npp",colnames(.))) %>% set_rownames(site_data$Sample) %>% scale() %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()

d_st_before <- site_data_st_before_imputation_imputed[,grepl("_st",names(site_data_st_before_imputation_imputed))] %>% select(!grep("pet|vap|ndvi|lai|npp",colnames(.))) %>% set_rownames(site_data$Sample) %>% scale() %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()

d_st_after <- site_data_st_after_imputation[,grepl("_st",names(site_data_st_after_imputation))]  %>% select(!grep("pet|vap|ndvi|lai|npp",colnames(.))) %>% set_rownames(site_data$Sample) %>% scale() %>% dist() %>% hclust(method = "ward.D2") %>% as.dendrogram()

dendro_site <- NULL
entang_measure <- NULL
for(i in c("xy","lt","st_before","st_after")){
    tmp <- get(paste0("d_",i))
    tmp %<>% 
        # set('labels_col',value = color_pool[color_xy])%>%
        set("branches_lty", 1)%>%
        set("branches_lwd",1)%>%
        set("labels_cex",rep(1,nrow(site_data_lt)))%>%
        as.ggdend()
    tmp$labels %<>%
        mutate(Sample=label)%>%
        left_join(site_data[,1:2])%>%
        mutate(col=color_pool[as.numeric(as.factor(Country))])
    tmp$segments%<>%
        mutate(across(c(2,4),~log(1+.x)))
    p <- tmp%>%ggplot(labels = F)+
        ggtitle(paste0("Site clustering ",i))+
        theme(plot.title = element_text(face="bold"))+
        coord_fixed()+
        coord_flip()+
        scale_y_reverse()+
        ylim(c(max(tmp$segments$y),-0.9))+
        theme_classic2()
    
    dendro_site[[i]] <- p+ ggnewscale::new_scale_color()+
        geom_text(tmp$labels,mapping=aes(x=x,y=y-c(max(tmp$segments$y)*c(.05,.2)),
                                         label=label,colour=Country),
                  key_glyph = draw_key_rect)+
        scale_colour_manual(values=color_pool)+
        theme(axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())+
        ylab("log( 1 + euclidean distance)")
    
    entang_measure[[i]] <- dendlist(d_xy,get(paste0("d_",i))) %>% 
        untangle(method = "step1side") %>% entanglement()
}

dendro_site_plots <- (dendro_site$xy+dendro_site$lt)/
    (dendro_site$st_before+dendro_site$st_after)+
    plot_layout(guides="collect")

ggsave("Figures/dendro_sites.png",dendro_site_plots,width=16,height =21)


ggsave(as.data.frame(entang_measure)%>%
           select(!xy)%>%
           melt()%>%
           ggplot(aes(x=variable,y=value))+
           geom_bar(stat='identity')+
           theme_classic2()+
           scale_y_continuous(expand = c(0,0))+
           coord_flip()+
           theme(axis.line.y = element_blank())+
           xlab("Dataset used for clustering")+
           ylab("Entanglement measure between XY tree and the 3 imputated datasets"),
       filename="Figures/entanglement_site_trees.png")
```

```{r}
p_list_dend <- NULL
p_list_km <- NULL
best_part <- NULL
for(i in c("xy","lt","st_before","st_after")){
    
    data_tmp <-
        if(i == "xy"){
            site_data_lt[,1:2] 
        }else{
            if(i == "lt"){
                imputed_data_lt[,grepl("_lt",names(imputed_data_lt))]%>%
                    select(!grep("pet|vap|ndvi|lai|npp",colnames(.)))
            }else{
                if(i == "st_before"){
                    site_data_st_before_imputation_imputed %>%
                        select(!grep("pet|vap|ndvi|lai|npp",colnames(.)))
                }else{
                    site_data_st_after_imputation %>% 
                        select(!grep("pet|vap|ndvi|lai|npp",colnames(.)))
                }
            }
        }
    
    cah_res <- data_tmp %>% scale() %>% factoextra::hcut(hc_method = "ward.D2", hc_metric = "euclidean", stand = TRUE)
    km2_res <- data_tmp %>% scale() %>%kmeans(centers = 2)
    km3_res <- data_tmp %>% scale() %>% kmeans(centers=3)
    km4_res <- data_tmp %>% scale() %>%kmeans(centers=4)
    km5_res <- data_tmp %>% scale() %>% kmeans(centers=5)
    km6_res <- data_tmp %>% scale() %>% kmeans(centers=6)
    km7_res <- data_tmp %>% scale() %>% kmeans(centers=7)
    
    estim_clust_1 <- scale(data_tmp) %>%
        factoextra::fviz_nbclust( FUNcluster =factoextra::hcut,
                                  method = "silhouette",
                                  hc_method = "average",
                                  hc_metric = "euclidean",
                                  stand = TRUE)+
        ggtitle("Silhouette")
    estim_clust_2 <- scale(data_tmp) %>%
        factoextra::fviz_nbclust( FUNcluster =factoextra::hcut,
                                  method = "wss",
                                  hc_method = "average",
                                  hc_metric = "euclidean",
                                  stand = TRUE)+
        ggtitle("WSS")
    estim_clust_3 <- scale(data_tmp) %>%
        factoextra::fviz_nbclust( FUNcluster =factoextra::hcut, 
                                  method = "gap_stat",
                                  hc_method = "average",
                                  hc_metric = "euclidean",
                                  stand = TRUE)+
        ggtitle("Gap Stat")
    
    nb <- NbClust::NbClust(scale(data_tmp), distance = "euclidean", min.nc = 2,
                           max.nc = 10, method = "kmeans")
    
    nb_clust <- fviz_nbclust_fixed(nb)+
        ggtitle(i)
    
    km_plot2 <- factoextra::fviz_cluster(km2_res, scale(data_tmp),geom="point")+
        ggtitle("k2")+
        theme_classic2()
    km_plot3 <- factoextra::fviz_cluster(km3_res, scale(data_tmp),geom="point")+
        ggtitle("k3")+
        theme_classic2()
    km_plot4 <- factoextra::fviz_cluster(km4_res, scale(data_tmp),geom="point")+
        ggtitle("k4")+
        theme_classic2()
    km_plot5 <- factoextra::fviz_cluster(km5_res, scale(data_tmp),geom="point")+
        ggtitle("k5")+
        theme_classic2()
    km_plot6 <- factoextra::fviz_cluster(km6_res, scale(data_tmp),geom="point")+
        ggtitle("k6")+
        theme_classic2()
    km_plot7 <- factoextra::fviz_cluster(km7_res, scale(data_tmp),geom="point")+
        ggtitle("k7")+
        theme_classic2()
    
    
    p_list_dend[[i]] <- (nb_clust+estim_clust_1)/
        (estim_clust_2+estim_clust_3)
    
    p_list_km[[i]] <- (km_plot2+km_plot3+km_plot4)/
        (km_plot5+km_plot6+km_plot7) 
    
    best_part[[i]] <- nb$Best.partition
}
```

```{r}
for( i in 1:length(p_list_dend)){
    ggsave(p_list_dend[[i]],
           filename=paste0("Figures/cluster_diag_",names(p_list_dend)[i],".png"))
    
    ggsave(p_list_km[[i]],
           filename=paste0("Figures/cluster_plot_",names(p_list_dend)[i],".png"),
           width = 10 ,
           height= 10)
}
```

Map best partitions

```{r}
site_best_part <- meco_16s_ns$sample_table
site_best_part %<>% mutate(bp_lt =as.factor(ifelse(best_part$lt==2,1,2)),
                           bp_st_bef =as.factor(best_part$st_before),
                           bp_st_aft = as.factor(best_part$st_after),
                           bp_xy = as.factor(best_part$xy))

(sites <- st_as_sf(site_best_part, coords = c("X", "Y"), 
                   crs = 4326, agr = "constant"))

map_best_part <- (ggplot() +
    geom_sf(data = world_lambert, fill = "grey", color = "darkgrey", size = 0.1)+
    geom_sf(data = sites, size = 3, shape = 23,aes(fill = bp_lt),alpha=.7) +
    ggtitle("Site clustering long term") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+
    theme(panel.background = element_rect(fill = "azure"))+
    scale_fill_manual(values = c("darkgreen","darkorange3"),name="Cluster - best partition"))/
    (ggplot() +
    geom_sf(data = world_lambert, fill = "grey", color = "darkgrey", size = 0.1)+
    geom_sf(data = sites, size = 3, shape = 23,aes(fill = bp_st_bef),alpha=.7) +
    ggtitle("Site clustering short term before") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+
    theme(panel.background = element_rect(fill = "azure"))+
    scale_fill_manual(values = c("darkgreen","darkorange3"),name="Cluster - best partition"))/
    (ggplot() +
    geom_sf(data = world_lambert, fill = "grey", color = "darkgrey", size = 0.1)+
    geom_sf(data = sites, size = 3, shape = 23,aes(fill = bp_st_aft),alpha=.7) +
    ggtitle("Site clustering short term after") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+
    theme(panel.background = element_rect(fill = "azure"))+
    scale_fill_manual(values = c("darkgreen","darkorange3"),name="Cluster - best partition"))+
    plot_layout(guides="collect")
    
ggsave(map_best_part,filename="Figures/map_best_part.png",height=21)



ggplot() +
    geom_sf(data = world_lambert, fill = "grey", color = "darkgrey", size = 0.1)+
    geom_sf(data = sites, size = 3, shape = 23,aes(fill = bp_xy),alpha=.7) +
    ggtitle("Site clustering xy") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+
    theme(panel.background = element_rect(fill = "azure"))+
    scale_fill_manual(values = c("darkgreen","darkorange3"),name="Cluster - best partition")
    ```


aet
pr
swe
tmmn
tmmx # delta tminmax?
vpd
elevation
evi
fpar
gpp
ghm
ph
soc
moist

```{r}
data_sets <- c("site_data_lt",
               'imputed_data_lt',
               "site_data_st_after_imputation",
               "site_data_st_before_imputation_imputed")

target_var <- c("aet",
                "pr",
                "swe",
                "tmmn",
                "tmmx", # delta tminmax?
                "vpd",
                "height",
                "evi",
                "fpar",
                "gpp",
                "gHM",
                "ph",
                "soc",
                "moist")

for(j in data_sets){
    message(j)
    tmp <- get(j)
    names(tmp) <- gsub("_.t$","",names(tmp))

    pairs <- ggpairs(tmp,
                     columns=which(colnames(tmp)%in%target_var),
                     progress = F,
                     axisLabels = "none",
                     # LOWER TRIANGLE ELEMENTS: add line with smoothing; make points transparent and smaller
                     lower = list(continuous = function(...) 
                         ggally_smooth(..., colour="darkolivegreen3", alpha = 0.3, size=0.8)), 
                     # DIAGONAL ELEMENTS: histograms
                     diag = list(continuous = function(...) 
                         ggally_barDiag(..., fill="grey")),
                     
                     # UPPER TRIANGLE ELEMENTS: use fct. creating corr heatmap with sign stars
                     upper = list(continuous = cor_fun)) + 
        theme(strip.background = element_blank(), # remove color
              strip.text = element_text(size=6,angle=45), # change font and font size
              axis.line = element_line(colour = "grey"),
              # remove grid
              panel.grid.minor = element_blank(), )+  # remove smaller gridlines
        # panel.grid.major = element_blank()    # remove larger gridlines)
        ggtitle(j)+
        theme(plot.title = element_text(face="bold"))
    
    
    ggsave(plot=pairs,
           filename=paste0("Figures/pairs_slctd_",j,".png"))
}
```


## On Communities

```{r}
d_com <- meco_16s_ns$otu_table %>% t() %>% labdsv::hellinger() %>%set_rownames(gsub("[A-Z]|-|(?=_).*","",rownames(.),perl = T))  %>% vegan::vegdist(method = "horn") %>% hclust(method = "ward.D2") %>% as.dendrogram()

diss_com <- meco_16s_ns$otu_table %>% t() %>% labdsv::hellinger() %>%set_rownames(gsub("[A-Z]|-|(?=_).*","",rownames(.),perl = T))

assign("comu_dissimilarity",as.matrix(vegdist(diss_com,method = "horn")))
```


```{r}
tmp <- d_com
    tmp %<>% 
        # set('labels_col',value = color_pool[color_xy])%>%
        set("branches_lty", 1)%>%
        set("branches_lwd",1)%>%
        set("labels_cex",rep(1,nrow(site_data_lt)))%>%
        as.ggdend()
    tmp$labels %<>%
        mutate(Sample=label)%>%
        left_join(site_data[,1:2])%>%
        mutate(col=color_pool[as.numeric(as.factor(Country))])
    tmp$segments%<>%
        mutate(across(c(2,4),~log(1+.x)))
    p <- tmp%>%ggplot(labels = F)+
        ggtitle(paste0("Site clustering ","com"))+
        theme(plot.title = element_text(face="bold"))+
        coord_fixed()+
        coord_flip()+
        scale_y_reverse()+
        ylim(c(max(tmp$segments$y),-0.9))+
        theme_classic2()
    
    p <- p+ ggnewscale::new_scale_color()+
        geom_text(tmp$labels,mapping=aes(x=x,y=y-c(max(tmp$segments$y)*c(.05,.2)),
                                         label=label,colour=Country),
                  key_glyph = draw_key_rect)+
        scale_colour_manual(values=color_pool)+
        theme(axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())+
        ylab("log( 1 + euclidean distance)")
    
    entang_measure_dcom <- dendlist(d_xy,d_com) %>% 
        untangle(method = "step1side") %>% entanglement()

    entang_measure_dcom
    
    
    dendlist(d_lt,d_com) %>% 
        untangle(method = "step1side") %>% entanglement()
```

```{r}
rm(a,b,alpha_plot,alpha_plot_pct,b,c,cleaned_metablist,cor_reads_asvs,curve16s,df_ordi_16s,distrib_reads_and_asvs,donut_comp,eig,hcd_site_XY,km_plot2,km_plot3,km_plot4,km_plot5,km_plot6,km_plot7,km2_res,km3_res,km4_res,km5_res,km6_res,km7_res,meco_16s_singl,metab_list,ordi16s,p,p_list_dend,p_list_km,pairs,pca_site_biplot,pca_sites,pca_sites_after_st,pca_sites_before_st,pca_sites_lt,pca_tmp,physeq,physeqcl,plot_alpha_divz_sgnsg,plot_data,plot_list,plots_comp,protox,rare16s,scree_plot_list,singleton_plots,singletons,site_XY_dend,tab16s,tmp16s,tmp,var_missing_val)
gc()
```