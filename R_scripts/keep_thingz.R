
### PCOA


ordi23s <- ordinate(file2meco::meco2phyloseq(data_23s_hell), "PCoA", "horn")

expl1 <- ordi23s$values$Relative_eig[1]
expl2 <- ordi23s$values$Relative_eig[2]
expl3 <- ordi23s$values$Relative_eig[3]
df_ordi_23s <- data.frame(dim1=ordi23s$vectors[,1],
                          dim2=ordi23s$vectors[,2],
                          dim3=ordi23s$vectors[,3],
                          country=data_23s_hell$sample_table$Country
)


pcoa23s <- ggplot(df_ordi_23s,aes(x=dim1,y=dim2,color=country))+
    geom_point()+
    theme_classic()+
    ylab(paste0("PCoA 2 (",round(expl2,digits = 3)*100,"%)"))+
    xlab(paste0("PCoA 1 (",round(expl1,digits = 3)*100,"%)"))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle("MAPP_23s")+
    theme(plot.title = element_text(face="bold"))

ordi16s <- ordinate(file2meco::meco2phyloseq(data_16s_hell), "PCoA", "horn")

expl1 <- ordi16s$values$Relative_eig[1]
expl2 <- ordi16s$values$Relative_eig[2]
expl3 <- ordi16s$values$Relative_eig[3]
df_ordi_16s <- data.frame(dim1=ordi16s$vectors[,1],
                          dim2=ordi16s$vectors[,2],
                          dim3=ordi16s$vectors[,3],
                          country=data_16s_hell$sample_table$Country
)


pcoa16s <- ggplot(df_ordi_16s,aes(x=dim1,y=dim2,color=country))+
    geom_point()+
    theme_classic()+
    ylab(paste0("PCoA 2 (",round(expl2,digits = 3)*100,"%)"))+
    xlab(paste0("PCoA 1 (",round(expl1,digits = 3)*100,"%)"))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle("MAPP_16s")+
    theme(plot.title = element_text(face="bold"))


ordi18s <- ordinate(file2meco::meco2phyloseq(data_18s_hell), "PCoA", "horn")

expl1 <- ordi18s$values$Relative_eig[1]
expl2 <- ordi18s$values$Relative_eig[2]
expl3 <- ordi18s$values$Relative_eig[3]
df_ordi_18s <- data.frame(dim1=ordi18s$vectors[,1],
                          dim2=ordi18s$vectors[,2],
                          dim3=ordi18s$vectors[,3],
                          country=data_18s_hell$sample_table$Country
)


pcoa18s <- ggplot(df_ordi_18s,aes(x=dim1,y=dim2,color=country))+
    geom_point()+
    theme_classic()+
    ylab(paste0("PCoA 2 (",round(expl2,digits = 3)*100,"%)"))+
    xlab(paste0("PCoA 1 (",round(expl1,digits = 3)*100,"%)"))+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    ggtitle("MAPP_18s")+
    theme(plot.title = element_text(face="bold"))















hierarchy_path_kegg <- path_kegg_pi3[,c(1,3,4)] %>%
    separate(col=classification,into = paste0("c",1:4),sep = ";")

# 
# maincol <- RColorBrewer::brewer.pal(7,'Dark2')
# color_frame <- data.frame(c1=unique(hierarchy_path_kegg$c1),
#                           color_c1=maincol)
# 
# hierarchy_path_kegg %<>%
#     left_join(color_frame)
# 
# hierarchy_path_kegg%>%
#     arrange(c1,c2,c3,c4)%>%
#     ggplot()+
#     geom_rect(aes(ymin = cumsum(observation_sum)-observation_sum,
#                   ymax = cumsum(observation_sum),
#                   xmin=0,
#                   xmax=.5),
#                   fill=arrange(hierarchy_path_kegg,c1)$color_c1)+
#     geom_rect(aes(ymin = cumsum(observation_sum)-observation_sum,
#                   ymax = cumsum(observation_sum),
#                   xmin=.75,
#                   xmax=1),
#                   fill=arrange(hierarchy_path_kegg,c1)$color_c1)+
#     geom_rect(aes(ymin = cumsum(observation_sum)-observation_sum,
#                   ymax = cumsum(observation_sum),
#                   xmin=1.25,
#                   xmax=1.5),
#                   fill=arrange(hierarchy_path_kegg,c1)$color_c1)+
#     geom_rect(aes(ymin = cumsum(observation_sum)-observation_sum,
#                   ymax = cumsum(observation_sum),
#                   xmin=1.75,
#                   xmax=2),
#                   fill=arrange(hierarchy_path_kegg,c1)$color_c1)+
#     coord_polar(theta='y')+
#     theme_void()






selected <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")

results_nbclust <- NULL
for (i in 1:length(selected)) {
    
    results_nbclust[[i]] <- try(NbClust(data = t(meco_16s_ns$otu_table),
                                        diss = comu_dissimilarity,
                                        min.nc=2,
                                        max.nc=15, 
                                        method="kmeans",
                                        index=selected[i]))
    
}

nb_clust_com <- fviz_nbclust_fixed(nb_com)

best_part_com <- nb_com$Best.partition


site_best_part %<>% mutate(bp_com=as.factor(best_part_com))

(sites <- st_as_sf(site_best_part, coords = c("X", "Y"), 
                   crs = 4326, agr = "constant"))

map_com_best_part <- (ggplot() +
                          geom_sf(data = world_lambert, fill = "grey", color = "darkgrey", size = 0.1)+
                          geom_sf(data = sites, size = 3, shape = 23,aes(fill = bp_com),alpha=.7) +
                          ggtitle("Site clustering 16s communities") +
                          theme_minimal()+
                          theme(plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+
                          theme(panel.background = element_rect(fill = "azure"))+
                          scale_fill_manual(values = c("darkgreen","darkorange3"),name="Cluster - best partition"))
map_com_best_part
       
clust_heatmap <- site_best_part%>%
    select(grep("Sample|bp_",colnames(.)))%>%
    melt("Sample")%>%
    mutate(variable=fct_relevel(variable,"bp_xy","bp_lt","bp_com","bp_st_bef","bp_st_aft"))%>%
    ggplot(aes(x=variable,y=Sample,fill=value))+
    geom_tile()+
    scale_fill_manual(values = c("darkgreen","darkorange3"),name="Cluster - best partition")+
    scale_y_discrete(limits=rev)+
    theme_minimal()+
    site_data%>%
    select(c(Country,Sample))%>%
    ggplot(aes(x="Country",y=Sample,fill=Country))+
    geom_tile()+
    theme_void()+
    scale_fill_manual(values=c("#7f61d3",
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
                               "#619dd7"),
                      name="Country")+
    plot_layout(guides="collect",widths = c(5,1))


ggsave("Figures/cluster_heatmaps.png",clust_heatmap,height = 28.7/2.54)
 
site_best_part$bp_st_aft==site_best_part$bp_com
 