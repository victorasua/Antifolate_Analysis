
#library packages
library(tidyverse)
library(maptools) 
library(RColorBrewer) 
library(classInt) 
library(sf)
library(ggiraph)
library(classInt)
library(ggpubr)
library(cowplot)
library(scico) #scientific colors
library(rgdal)
library(ggmap)
library(ggspatial)


#### Uganda spatial data
ug_district_facet <- readRDS("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Documents/PhD_Files/thesis_data/data_raw_tes/MapPlotData.rds")


################### Study Site Names ##############
mi_study_districts <- ug_district_facet %>% 
  filter(district_name %in% c("Koboko", 
                              "Arua", 
                              "Lamwo", 
                              "Kaabong",
                              "Agago",
                              "Katakwi",
                              "Amolatar",
                              "Hoima",
                              "Kapchorwa",
                              "Tororo",
                              "Jinja",
                              "Mubende",
                              "Kasese",
                              "Kanungu",
                              "Kabale",
                              "Kole"))


################# Number genotype per site #########

##2 Read facility pgs locatiob
facility <- readxl::read_xls("/Users/victorasua/Desktop/prism_r/prx_all_data/raw_data/rawdata/map_data_2021.xls")
facility <- facility %>% data.frame() %>% dplyr::select(district, facility, N, lon, lat) %>% unique() %>% drop_na()
names(facility)[4] <- "lat"
names(facility)[5] <- "lon"

facility$N[facility$district=="Agago"]<-488
facility$N[facility$district=="Lamwo"]<-486
facility$N[facility$district=="Kole"]<-496
facility$N[facility$district=="Amolatar"]<-479
facility$N[facility$district=="Kaabong"]<-402
facility$N[facility$district=="Katakwi"]<-398
facility$N[facility$district=="Kapchorwa"]<-243
facility$N[facility$district=="Tororo"]<-444
facility$N[facility$district=="Jinja"]<-494
facility$N[facility$district=="Mubende"]<-491
facility$N[facility$district=="Kasese"]<-376
facility$N[facility$district=="Hoima"]<-396
facility$N[facility$district=="kanungu"]<-497
facility$N[facility$district=="Kabale"]<-309
facility$N[facility$district=="Arua"]<-492
facility$N[facility$district=="Koboko"]<-399

facility$district[facility$district=="kanungu"]<-"Kanungu"
### the plot
##site coordinates
coords <- read.csv("/Users/victorasua/Documents/PROJECTS/IMMRSE-U/2024_data_analysis/raw_data/survey_site_gps_cords.csv")
coord <- coords |> dplyr::select(District, lon, lat) %>% data.frame()
data_sf <- st_as_sf(coord, coords = c("lat", "lon"), crs = 4326)

########### Combine number genotyped and site

siteN <- facility |> select(district, N)
mi_study_districts2 <- mi_study_districts |> 
  left_join(siteN,by = c("district_name" = "district" )) |>
  mutate(siteN = paste(district_name, N))

# print words on separate lines
mi_study_districts2 <- mi_study_districts2 %>%
  dplyr::mutate(siteN = gsub(" ", "\n", siteN))


pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.1.pdf", 
    width = 6.8, height = 7)

ggplot ( ) + 
  geom_sf(data = ug_district_facet, lwd = .1, fill = "white", size = 0.01 )  +
  geom_sf(data = mi_study_districts, lwd = .1, fill = "lightblue4", size = 0.01, alpha = 0.05) +
  geom_sf_text(data = mi_study_districts2, aes(label = siteN), 
               size = 4, color = "orange4") +
  labs(fill = "Prevalence")+
  theme_bw() +
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis labels
    axis.ticks = element_blank()  # Remove axis ticks
  ) + 
  annotation_north_arrow(
    location = "tl",  # Location of the arrow (tl = top left, tr = top right, bl = bottom left, br = bottom right)
    which_north = "true",  # Use true north
    pad_x = unit(0.5, "in"),  # Padding from the plot edges
    pad_y = unit(0.5, "in"),
    height = unit(0.45, "in"),  # Height of the north arrow
    width = unit(0.2, "in")  # Width of the arrow
  )

dev.off()



