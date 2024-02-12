data_sat <- read.csv("ressource/Site_edaphic_data/MappPtsSampled.csv")


trans_hellinger <- labdsv::hellinger(t(meco_23s$otu_table))
trans_hellinger %<>% as.data.frame()
# View(trans_hellinger)
trans_hellinger$siteCol <- rownames(trans_hellinger)
trans_hellinger %<>% mutate(siteCol = gsub("[A-Z]|-|(?=_).*","",siteCol,perl = T))
trans_hellinger$Sample <- trans_hellinger$siteCol


site_coordinates <- read.csv("ressource/Site_edaphic_data/site_coord_MAPP.csv") 
site_coordinates %<>% 
    mutate(Sample=ifelse(nchar(Sample)<3,paste0("0",Sample),Sample))%>%
    mutate(Sample=ifelse(nchar(Sample)<3,paste0("0",Sample),Sample))
site_coordinates %<>% 
    mutate(Sample=ifelse(grepl("^8",Sample),paste0("0",Sample),Sample))


test <- left_join(trans_hellinger,site_coordinates[,c(1:3)])
test %<>% relocate(Sample,siteCol,X,Y) %>% select(-1)


data_sat %<>% arrange(Sample)%>% 
    mutate(Sample=ifelse(nchar(Sample)<3,paste0("0",Sample),Sample))%>%
    mutate(Sample=ifelse(nchar(Sample)<3,paste0("0",Sample),Sample))
site_coordinates %<>% 
    mutate(Sample=ifelse(grepl("^8",Sample),paste0("0",Sample),Sample))
test %<>% arrange(siteCol)

test_data_sat <- filter(data_sat,Sample%in%test$siteCol)
test %<>% filter(!siteCol%in%c("116","119","166"))

test_data_sat %<>% rename(siteCol=Sample)
test_data_sat %<>% select(c("siteCol","X","Y","aet_lt","bulk","def_lt","evi_lt","gHM","gpp_lt","height",'lai_lt',"moist_lt","pdsi_lt","ph","pr_lt","ssm_lt","susm_lt","swe_lt","tmmn_lt","tmmx_lt","vpd_lt"))

argh <- gdm::formatsitepair(bioData = test,
                            bioFormat = 1,
                            dist="horn",
                            abundance = T,
                            siteColumn = "siteCol",
                            XColumn = "X",
                            YColumn = "Y",
                            predData = test_data_sat)


var_cor <- GGally::ggpairs(
    test_data_sat %>% select(c("X","Y","aet_lt","bulk","def_lt","evi_lt","gHM","gpp_lt","height",'lai_lt',"moist_lt","npp_lt","pdsi_lt","pet_lt","ph","pr_lt","ssm_lt","susm_lt","swe_lt","tmmn_lt","tmmx_lt","vap_lt","vpd_lt"))%>% mutate(across(.fns = ~replace(., . ==  -9999 , NA))),
    progress = T,
    cardinality_threshold = 105
)
ggsave("Figures/cov_cor.png",var_cor,height = 15,width = 15)

gdm_13s <- gdm::gdm(argh,
         geo=T)
summary(gdm_13s)
plot(gdm_13s)







test_data_sat2 <- filter(data_sat,Sample%in%test$siteCol)
test_data_sat2 %<>% rename(siteCol=Sample)
test_data_sat2 %<>% select(c("siteCol","X","Y","aet_03","bulk","def_03","evi_03","gHM","gpp_03","height",'lai_03',"moist_03","pdsi_03","ph","pr_03","ssm_03","susm_03","swe_03","tmmn_03","tmmx_03","vpd_03"))


argh2 <- gdm::formatsitepair(bioData = test,
                             bioFormat = 1,
                             dist="horn",
                             abundance = T,
                             siteColumn = "siteCol",
                             XColumn = "X",
                             YColumn = "Y",
                             predData = test_data_sat2) 

gdm_13s_march <- gdm::gdm(argh2,
                          geo=T)


plot(gdm_13s_march)


# CHECK the low dissimilarites