library(readr)

count_16s_and_NTSI <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy179-[FROGSFUNC_1_placeseqs_and_copynumbers__frogsfunc_marker.tsv].tsv")
norm_cluster_by_site <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy183-[FROGSFUNC_2_functions__frogsfunc_functions_marker_norm.tsv].tsv")
weighed_NTSI <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy184-[FROGSFUNC_2_functions__frogsfunc_functions_weighted_nsti.tsv].tsv")
picrust2_excluded <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy185-[FROGSFUNC_2_functions__frogsfunc_functions_asv_excluded.tsv].tsv")
fct_kegg <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy189-[FROGSFUNC_2_functions___frogsfunc_functions_unstrat_KO.tsv].tsv")
path_kegg <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy191-[FROGSFUNC_3_pathways__frogsfunc_pathways_unstrat.tsv].tsv")
biom_test <- read_tsv("R_scripts/paper_16s/picrust2/Galaxy198-[FROGS_BIOM_to_TSV__abundance.tsv].tsv")

length(unique(fct_kegg$classification))


fct_kegg_list <- NULL
fct_kegg_list$metaD <- fct_kegg[,1:4]
fct_kegg_list$pseudo_com <- fct_kegg[,5:ncol(fct_kegg)]

library(magrittr)
library(dplyr)
library(ggplot2)
fct_kegg%>%
  select(-c(2:4))%>%
  reshape2::melt("classification")%>%
  filter(classification!="NA")%>%
  ggplot(aes(x=variable,y=value,fill=classification))+
  geom_bar(stat="identity")

path_kegg%>%
  select(-c(2:4))%>%
  reshape2::melt("classification")%>%
  filter(classification!="NA")%>%
  mutate(path=gsub(";.*$","",classification))%>%
  ggplot(aes(x=variable,y=value,fill=path))+
  geom_bar(stat="identity")

library(stringr)

path_kegg%>%
  select(-c(2:4))%>%
  reshape2::melt("classification")%>%
  filter(classification!="NA")%>%
  mutate(path=gsub(";.*$","",classification))%>%
  filter(path=="Metabolism")%>%
  mutate(path2=str_extract(classification,"(;.+?;)"))%>%
  group_by(variable)%>%
  mutate(val_norm=value/sum(value))%>%
  ggplot(aes(x=variable,y=val_norm,fill=path2))+
  geom_bar(stat="identity")




path_kegg%>%
  select(-c(2:4))%>%
  reshape2::melt("classification")%>%
  filter(classification!="NA")%>%
  mutate(path=gsub(";.*$","",classification))%>%
  filter(grepl("Carbohydrate",classification))%>%
  group_by(variable)%>%
  mutate(val_norm=value/sum(value))%>%
  ggplot(aes(x=variable,y=val_norm,fill=classification))+
  geom_bar(stat="identity")

path_kegg%>%
  select(-c(2:4))%>%
  reshape2::melt("classification")%>%
  filter(classification!="NA")%>%
  mutate(path=gsub(";.*$","",classification))%>%
  filter(grepl("Amino",classification))%>%
  group_by(variable)%>%
  mutate(val_norm=value/sum(value))%>%
  ggplot(aes(x=variable,y=val_norm,fill=classification))+
  geom_bar(stat="identity")

library(tidyr)

path_kegg%>%
  select(-c(2:4))%>%
  reshape2::melt("classification")%>%
  filter(classification!="NA")%>%
  separate_wider_delim(cols = classification,
                       delim = ";",
                       names = paste0("path_",1:4))%>%
  View()



meco_16s_singl$tidy_dataset()
meco_16s_singl$sample_table %<>% mutate(dumb_group = "sample")
meco_16s_singl$cal_abund()
meco_16s_singl <- trans_abund$new(meco_16s_singl,taxrank = "Order_conf_rdp",ntaxa = 8,groupmean = "dumb_group")

plot_data <- meco_16s_singl$data_abund
use_taxanames <- c(meco_16s_singl$data_taxanames,"unidentified")
plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
plot_data %<>% 
  dplyr::group_by(!!!syms(c("Taxonomy","Sample"))) %>% 
  dplyr::summarise(Abundance = sum(Abundance)) %>% 
  as.data.frame(stringsAsFactors = FALSE)
plot_data$Taxonomy %<>% factor(., levels = c(use_taxanames[-9], "Others", "unidentified"))
plot_data$label <- paste0(round(plot_data$Abundance, 1),"%")
donut_comp <- ggdonutchart(plot_data,"Abundance",
                           fill="Taxonomy",
                           label = "label",
                           color = "white",
                           palette = c(colorRampPalette(brewer.pal(8, "Dark2"))(9),"lightgrey"))+
  theme(legend.position = "right")
