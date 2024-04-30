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