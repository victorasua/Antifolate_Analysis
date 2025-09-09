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
library(dplyr)
library(purrr)
library(broom)




# data
dr_ibdSEGT_genot8 <- read.csv( "/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/dr_ibdSEGT_genot8.csv")

propIBD <- dr_ibdSEGT_genot8 %>%
  filter(!is.na(Year), !is.na(trueGeno), site_grp == "withinSite") %>% 
  ungroup() %>% 
  filter(!trueGeno%in%c("PfDHFR I164L mixed", "PfDHPS A581G mixed"))  %>%
select(iid1,iid2,site,Year,trueGeno,GenPropShared,site_grp)  

# plotting
pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.5.pdf", 
    width = 7, height = 9) 

propIBD |>
  ggplot(aes(x = factor(Year), y = GenPropShared, color = trueGeno)) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 0.65) +
  scale_color_manual(values = c(
    "PfDHFR I164L" = "#E41A1C",   # red
    "PfDHPS A581G" = "#E69F00",   # blue
    "Reference"    = "#377EB8"), 
    labels = c(
      "PfDHFR I164L" = "PfDHFR 164L",
      "PfDHPS A581G" = "PfDHPS 581G",
      "Reference"    = "Wild-type")
  ) +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  facet_wrap(~ site, ncol = 3) +
  labs(
    #title = "Shared IBD Segments by Site and Year",
    x = "Year",
    y = "Proportion of genome IBD",
    color = "Genotype:"
  ) +
  #theme_pubr() +
  theme(strip.text = element_text(size = 9, face = "bold",colour = "black"),
        axis.title = element_text(size = 9, face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 270, hjust = 1,colour = "black"),
        axis.text = element_text(size = 8, face = "bold",colour = "black"),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 8, hjust = 0.5, vjust = 0.5, margin = margin(l = 1, unit = "pt"),colour = "black"),
        legend.title = element_text(size = 8, face = "bold",colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.margin = margin(t = -12, b = 0, unit = "pt"),  
        panel.spacing.x = unit(0.2, "line"),
        panel.spacing.y = unit(0.2, "line"),
        panel.background = element_blank(),
        #plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
        panel.grid.minor.y = element_blank(),   # turn off minor gridlines
        panel.grid.major.x = element_line(color = "grey80", linewidth = 0.1),   # no vertical grid lines
        panel.grid.minor.x = element_blank()
        
  )

dev.off()


