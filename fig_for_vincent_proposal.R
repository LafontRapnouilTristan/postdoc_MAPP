library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(microeco)
library(patchwork)
library(flextable)

agg.table.taxo <- function(tab, tax.lvl="genus", tax.table) {
  tax.table <- tax.table[match(rownames(tab),
                               rownames(tax.table)),]
  message(paste('Table aggregation to the', tax.lvl, "level."))
  #message('Please be sure that the ASV/OTU table and the taxonomy table are ordered the same way')
  if(nrow(tab) != nrow(tax.table)) stop("The ASV/OTU table and the taxonomy table do not have the same number of rows")
  tax <- tax.table[,grep(tax.lvl, colnames(tax.table), ignore.case = T)]
  tax[is.na(tax)] <- "Unknown"
  tab <- stats::aggregate(tab, by=list("taxo"=tax), FUN=sum)
  rownames(tab) <- tab[,1]  
  tab <- tab[,-1]
  return(tab)
}


load("output/metabarlists.Rdata")
rm(metab_lists)
cleaned_metablist <- list(MAPP_16s=cleaned_metablist$MAPP_16s,MAPP_18s=cleaned_metablist$MAPP_18s)
cleaned_metablist$MAPP_16s$samples%>%nrow #218 samples
cleaned_metablist$MAPP_18s$samples%>%nrow

sample_names_16s <- gsub("[A-Z]|-|(?=_).*","",cleaned_metablist$MAPP_16s$samples$sample_id,perl = T)
sample_names_18s <- gsub("[A-Z]|-|(?=_).*","",cleaned_metablist$MAPP_18s$samples$sample_id,perl = T)

setdiff(sample_names_16s,sample_names_18s)

biome_data <- read.csv("ressource/Site_edaphic_data/MAPP_ecoR_with_Ecoregions.csv")
site_coordinates <- read.csv("ressource/Site_edaphic_data/site_coord_MAPP.csv")  # get site coordinates
site_coordinates <- left_join(site_coordinates,biome_data,by = "Sample")

site_coordinates%<>%
  filter(!grepl("23[1-7]|239|240|241|042-047|047-042|246b|256a|87b|116|119|166",Sample))%>%
  mutate(Sample=ifelse(nchar(Sample)<2,paste0("00",Sample),
                       ifelse(nchar(Sample)<3,paste0("0",Sample),Sample)))

site_coordinates[which(site_coordinates$Sample=="87a"),]$Sample <- "087a"
site_coordinates[which(site_coordinates$Sample=="087a"),]$Sample 
site_coordinates[which(site_coordinates$Sample=="087a"),2:ncol(site_coordinates)] <- site_coordinates[which(site_coordinates$Sample=="087"),2:ncol(site_coordinates)]


site_coordinates[which(site_coordinates$Sample%in%c("020","021","022")),"BIOME_NAME"] <- "Boreal Forests/Taiga"
site_coordinates[which(site_coordinates$Sample%in%c("020","021","022")),"ECO_NAME"] <- "Scandinavian and Russian taiga"


site_coordinates[which(site_coordinates$Sample%in%c("092","093","094")),"BIOME_NAME"] <- "Tundra"
site_coordinates[which(site_coordinates$Sample%in%c("092","093","094")),"ECO_NAME"] <- "Scandinavian Montane Birch forest and grasslands"



world <-rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Subset the world data to include only the specified latitude range
cropping <- as(raster::extent(-180, 180, 20, 90), "SpatialPolygons")


world_subset <- raster::crop(sf::as_Spatial(world),cropping)
# Create a Lambert azimuthal equal-area projection centered on the North Pole
lambert_proj <- sf::st_crs("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
world_lambert <- sf::st_transform(sf::st_as_sf(world_subset), crs = lambert_proj)

biomes_color <- c("#F6D55C","#ED553B","#173F5F","#3CAEA3","#20639B")

for (i in c("16s","18s")){
  
  metab_tmp <- cleaned_metablist[[paste0("MAPP_",i)]]
  reads_tmp <- metab_tmp$reads
  samples_tmp <- metab_tmp$samples
  motus_tmp <- metab_tmp$motus

  rownames(reads_tmp) <- gsub("[A-Z]|-|(?=_).*","",rownames(reads_tmp),perl = T)
  rownames(samples_tmp) <- gsub("[A-Z]|-|(?=_).*","",rownames(samples_tmp),perl = T)
  samples_tmp%<>%
    mutate(Sample=rownames(.))
  
  samples_tmp%<>%
  filter(!grepl("23[1-7]|239|240|241|042-047|047-042|246b|256a|87b|116|119|166",sample_id))
  
  reads_tmp%<>%
    as.data.frame()%>%
    filter(!grepl("23[1-7]|239|240|241|042047|047042|246b|256a|87b|116|119|166",rownames(.)))
  
  reads_tmp <- reads_tmp[,which(colSums(reads_tmp)!=0)]
  
  motus_tmp%<>%
    mutate(asvs=rownames(.))
  motus_tmp%<>%filter(asvs%in%colnames(reads_tmp))
  
  samples_tmp%<>%
    mutate(Sample=ifelse(nchar(Sample)<2,paste0("00",Sample),
                                      ifelse(nchar(Sample)<3,paste0("0",Sample),Sample)))%>%
    left_join(site_coordinates)
  
 
  
  (sites <- sf::st_as_sf(samples_tmp, coords = c("X", "Y"), 
                         crs = 4326, agr = "constant"))
  
  # Plot the map
  map <- ggplot() +
    geom_sf(data = world_lambert, fill = "grey", color = "darkgrey", size = 0.1)+
    geom_sf(data = sites, size = 3, shape = 23,aes(fill = BIOME_NAME),alpha=.7) +
    ggtitle("Site Biomes") +
    theme(plot.title = element_text(hjust = 0.5, size = 16))+
    theme_minimal()+
    theme(panel.background = element_rect(fill = "azure"))+
    scale_fill_manual(values=biomes_color)
  
  ggsave(paste0("output/vin_proposal/mapp_",i,".pdf"),map)
  micro_tmp <- microeco::microtable$new(otu_table = as.data.frame(t(reads_tmp)),
                           sample_table = samples_tmp,
                           tax_table = motus_tmp)
  micro_tmp$tidy_dataset()
  
  micro_tmp$sample_table%<>%mutate(nb_asvs = colSums(micro_tmp$otu_table>0),
                                   nb_reads = colSums(micro_tmp$otu_table))
  micro_tmp$sample_table %<>% mutate(dumb_group = "sample")
  micro_tmp$otu_table %<>% filter(rowSums(.)>1)
  micro_tmp$tidy_dataset()
  
  micro_tmp$cal_abund()
  
  plots_comp <- NULL
  for( j in c("Order_conf_rdp","Genus_conf_rdp")){
   
      tmp16s1 <- trans_abund$new(micro_tmp,taxrank = j,ntaxa = 8,groupmean = "dumb_group")
      
      plot_data1 <- tmp16s1$data_abund
      use_taxanames <- tmp16s1$data_taxanames
      plot_data1$Taxonomy[!plot_data1$Taxonomy %in% c(use_taxanames,"unidentified")] <- "Others"
      plot_data1 %<>% 
        dplyr::group_by(!!!syms(c("Taxonomy","Sample"))) %>% 
        dplyr::summarise(Abundance = sum(Abundance)) %>% 
        as.data.frame(stringsAsFactors = FALSE)
      plot_data1$Taxonomy %<>% factor(., levels = c(use_taxanames[-9], "Others","unidentified"))
      plot_data1$label <- paste0(round(plot_data1$Abundance, 1),"%")
      donut_comp <- ggpubr::ggdonutchart(plot_data1,"Abundance",
                                 fill="Taxonomy",
                                 label = "label",
                                 color = "white",
                                 palette =  c(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(9),"lightgrey"),
                                 guide = guide_legend(reverse = TRUE) )+
        theme(legend.position = "right")
      
      tmp16s2 <- trans_abund$new(micro_tmp,taxrank = j,ntaxa = 8)
      plot_data2 <- tmp16s2$data_abund
      plot_data2$Taxonomy[!plot_data2$Taxonomy %in% c(use_taxanames,"unidentified")] <- "Others"
      plot_data2$Taxonomy %<>% factor(., levels = rev(c(use_taxanames[-9], "Others","unidentified")))
      plot_data2 <- left_join(plot_data2,site_coordinates[,c("Sample","BIOME_NAME")])
      plot_data2%<>%
        arrange(BIOME_NAME)%>%
        mutate(Sample=as.factor(Sample))%>%
        mutate(Sample=forcats::fct_relevel(Sample,unique(as.character(Sample))))
      bar_comp <- (plot_data2%>%
                     ggplot(aes(x=Sample,y=Abundance,fill=Taxonomy),color=NA)+
                     geom_bar(stat = "identity",width=1,position='fill')+
                     scale_fill_manual(values =  rev(c(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(9),"lightgrey")),
                                       guide = guide_legend(reverse = TRUE) )+
                     scale_y_continuous(expand = c(0,0))+
                     theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           plot.margin = unit(c(0,0,0,0),units="mm"),
                           axis.title.x = element_blank()))/(
                             ggplot(plot_data2,aes(x=Sample,y=1,fill=BIOME_NAME))+
                               geom_tile()+
                               scale_fill_manual(values=biomes_color)+
                               theme_void())+plot_layout(guides="collect",heights = c(9,1))+plot_annotation(title = i)
      
      
      
      plots_comp[[j]] <- donut_comp + bar_comp  
      }
  p_rich <- micro_tmp$sample_table%>%
  ggplot(aes(x=BIOME_NAME,y=nb_asvs,color=BIOME_NAME))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = .05)+
    scale_color_manual(values = biomes_color,name="Biomes")+
    ggpubr::theme_classic2()+
    xlab("")+
    ylab("Richness")+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
  
  ggsave(paste0("output/vin_proposal/p_rich_",i,".pdf"),p_rich)
  ggsave(paste0("output/vin_proposal/p_comp_order_",i,".pdf"),plots_comp[["Order_conf_rdp"]])
  ggsave(paste0("output/vin_proposal/p_comp_genus_",i,".pdf"),plots_comp[["Genus_conf_rdp"]])
## Taxa summary
  
    d <- micro_tmp$tax_table[,c("Order_conf_rdp","Class_conf_rdp","Phylum_conf_rdp")]
    d%<>%mutate(Order_conf_rdp=paste0("o__",Order_conf_rdp),
                Class_conf_rdp=paste0("c__",Class_conf_rdp),
                Phylum_conf_rdp=paste0("p__",Phylum_conf_rdp))
    tab <- micro_tmp$otu_table 
    
    tax.lvl <- "Order"
    order_tot <- agg.table.taxo(tab = tab, # aggregate my data at the targeted tax level
                               tax.lvl = tax.lvl,
                               tax.table = d) %>%
      apply(1,sum) %>% # sum the number of reads
      sort %>% # sort in ascending order
      identity
      tax.lvl <- "Class"
    class_tot <- agg.table.taxo(tab = tab, # aggregate my data at the targeted tax level
                               tax.lvl = tax.lvl,
                               tax.table = d) %>%
      apply(1,sum) %>% # sum the number of reads
      sort %>% # sort in ascending order
      identity
    tax.lvl <- "Phylum"
    phylum_tot <- agg.table.taxo(tab = tab, # aggregate my data at the targeted tax level
                                tax.lvl = tax.lvl,
                                tax.table = d) %>%
      apply(1,sum) %>% # sum the number of reads
      sort %>% # sort in ascending order
      identity
    top10_order <- names(order_tot[(length(order_tot)-ifelse(length(order_tot)<10,length(order_tot),10)):length(order_tot)])
    top10_order <- top10_order[which(top10_order!="o__")]
    top10_class <- names(class_tot[(length(class_tot)-ifelse(length(class_tot)<10,length(class_tot),10)):length(class_tot)])
    top10_class <- top10_class[which(top10_class!="c__")]
    top10_phylum <- names(phylum_tot[(length(phylum_tot)-ifelse(length(phylum_tot)<10,length(phylum_tot),10)):length(phylum_tot)])
    top10_phylum <- top10_phylum[which(top10_phylum!="p__")]
    
    d <- d %>% mutate(Phylum_conf_rdp=ifelse(!Phylum_conf_rdp%in%c("p__",top10_phylum),"p__others",Phylum_conf_rdp),
                      Class_conf_rdp=ifelse(!Class_conf_rdp%in%c("c__",top10_class),"c__others",Class_conf_rdp),
                      Order_conf_rdp=ifelse(!Order_conf_rdp%in%c("o__",top10_order),"o__others",Order_conf_rdp))
    d <- d %>%  mutate(Phylum_conf_rdp=ifelse(Phylum_conf_rdp=="p__","p__unknown",Phylum_conf_rdp),
                       Class_conf_rdp=ifelse(Class_conf_rdp=="c__","c__unknown",Class_conf_rdp),
                       Order_conf_rdp=ifelse(Order_conf_rdp=="o__","o__unknown",Order_conf_rdp))
 
    
    d$reads <- micro_tmp$taxa_sums()
    
    df <- d %>%
      ggsankey::make_long(1:3)
    
    TotalCount = nrow(d)
    # Step 2
    dagg <- df%>%
      dplyr::group_by(node)%>%
      tally()
    
    dagg <- dagg%>%
      dplyr::group_by(node)%>%
      dplyr::mutate(pct = n/TotalCount)%>%
      dplyr::mutate(tax_lvl=ifelse(grepl("c__",node),"class",
                                   ifelse(grepl("o__",node),"order","phylum")))
    
    
    # Table recap
    
    phylum <- filter(dagg,tax_lvl=="phylum")
    phylum %<>% mutate(node=gsub("p__","",node),
                       abun = paste0(n," (",round(pct*100, digits=2),"%)"))%>%
      dplyr::select(-c(2:4))
    
    class <- filter(dagg,tax_lvl=="class")
    class %<>% mutate(node=gsub("c__","",node),
                      abun = paste0(n," (",round(pct*100, digits=2),"%)"))%>%
      dplyr::select(-c(2:4))
    
    order <- filter(dagg,tax_lvl=="order")
    order %<>% mutate(node=gsub("o__","",node),
                      abun = paste0(n," (",round(pct*100, digits=2),"%)"))%>%
      dplyr::select(-c(2:4))
    
    if(nrow(phylum)<12){
      phylum <- rbind(phylum,data.frame(node=rep(NA,12-nrow(phylum)),abun=rep(NA,12-nrow(phylum))))
    }
    
    smry_tab <- cbind(phylum,class,order)
    smry_tab %<>% 
      dplyr::add_row(.after = 10)%>%
      flextable()%>%
      set_header_labels(node...1="Name",abun...2="Count",node...3="Name",abun...4="Count",node...5="Name",abun...6="Count")%>%
      add_header_row(top = TRUE, values = c("Phylum", "Class","Order"), colwidths = c(2,2,2))%>%
      align(align = "center", part = "all") %>%
      autofit()%>%
      bold(part='header')
    
    
    save_as_docx(smry_tab,path=paste0("output/vin_proposal/taxa_smry_",i,".docx"))

  
  
}





