install.packages("tmap")
library(tmap)
site_coord <- read.csv("ressource/Site_edaphic_data/site_coord_MAPP.csv")

install.packages("sf")
library(sf)
str(site_coord)
site_coord%<>%mutate(Y=as.numeric(Y))
site_coord_no_na <- site_coord%>%drop_na(X,Y)
MAPP_site_sf <- st_as_sf(site_coord_no_na,coords = c("X","Y")) %>% st_set_crs(3857)


MAPP_site_sf%>%
    ggplot()+
    geom_sf()

data("World")
World %>%
    ggplot()+
    geom_sf()


world_map <- map_data("world")



#Creat a base plot with gpplot2
p <- ggplot() + coord_fixed() +
    xlab("") + ylab("")


p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                 colour="#9e7221", fill="#c9a563")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = 'lightgrey', colour = 'white'), 
          axis.line = element_line(colour = "white"), legend.position="none",
          axis.ticks=element_blank(), axis.text.x=element_blank(),
          axis.text.y=element_blank())+
    geom_point(data=site_coord,
               aes(x=Y,y=X),color="darkblue",fill="steelblue",pch=21,size=3,alpha=.7)+
    # Convert to polar coordinates
    coord_map("ortho", orientation = c(90, 0, 0))
