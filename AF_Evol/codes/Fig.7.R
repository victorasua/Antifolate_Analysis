
library(flexmix)
library(lattice)
library(moimix)
library(tidyverse)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)
library(readr)
library(SeqVarTools)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(vcfR)
library(ape)
library(isoRelate)
library(stringr)
library(forcats)
library(data.table)
library(stringi)
library(data.table)
library(ggraph)
library(cowplot)
library(igraph)


#gen_dat<-read.csv("/Users/victorasua/Documents/PROJECTS/Genomic/AF_analysis_2024.10.11/PCA/gen_dat.csv")
gen_dat<-read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/PCA/gen_dat.csv")
gen_dat_AF <- gen_dat|>filter(aa_change=="Ile164Leu"| aa_change=="Ala581Gly")
gen_dat_AF$iid <- str_remove_all(gen_dat_AF$id, pattern =  "-P.*")
names(gen_dat_AF)[1] <- "iid1"

dr_gen_dat <- gen_dat |> dplyr::select(id, aa_change, genotype)|>
  pivot_wider(names_from = aa_change, values_from = genotype) |> unique() |> 
  dplyr::select(id, Ile164Leu, Ala581Gly) 


dr_gen_dat$Ile164Leu[dr_gen_dat$Ile164Leu==0]<-"Reference"
dr_gen_dat$Ile164Leu[dr_gen_dat$Ile164Leu==1]<-"PfDHFR I164L"
dr_gen_dat$Ile164Leu[dr_gen_dat$Ile164Leu==2]<-"PfDHFR I164L"

dr_gen_dat$Ala581Gly[dr_gen_dat$Ala581Gly==0]<-"Reference"
dr_gen_dat$Ala581Gly[dr_gen_dat$Ala581Gly==1]<-"PfDHPS A581G"
dr_gen_dat$Ala581Gly[dr_gen_dat$Ala581Gly==2]<-"PfDHPS A581G"

dr_genotype2 <- dr_gen_dat |> mutate(trueGeno=case_when(Ile164Leu=="Reference"&Ala581Gly=="Reference" ~ "Reference",
                                                        Ile164Leu=="PfDHFR I164L"&Ala581Gly=="Reference"| Ile164Leu=="PfDHFR I164L"&is.na(Ala581Gly) ~ "PfDHFR I164L",
                                                        Ile164Leu=="Reference"&Ala581Gly=="PfDHPS A581G"| is.na(Ile164Leu)&Ala581Gly=="PfDHPS A581G" ~ "PfDHPS A581G",
)) |> dplyr::select(id, trueGeno)
names(dr_genotype2)[1] <- "iid" 


# my_p_clusters <- readRDS("/Users/victorasua/Documents/PROJECTS/Genomic/AF_analysis_2024.10.11/my_p_clusters_25.rds")
# my_p_clusters <- readRDS("/Users/victorasua/Documents/PROJECTS/Genomic/AF_analysis_2024.10.11/my_p_clusters_5.rds")
# my_p_clusters <- readRDS("/Users/victorasua/Documents/PROJECTS/Genomic/AF_analysis_2024.10.11/my_p_clusters_75.rds")

#my_p_clusters_50 <- readRDS("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/my_p_clusters_5.rds")
#my_p_clusters_75 <- readRDS("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/my_p_clusters_75.rds")
#my_p_clusters_95 <- readRDS("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/my_p_clusters_pruned_imputed.rds")
my_p_clusters_Dmnt_alleles_90b <- readRDS("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/my_p_clusters_clusters_Dmnt_alleles_90b.rds")


library(iplot)
# scico::scico_palette_show()

# plot the network of clusters
#iplot <- my_p_clusters_75$i.network
#iplot <- my_p_clusters_95$i.network
iplot <- my_p_clusters_Dmnt_alleles_90b$i.network

plot(iplot,
     vertex.label=NA,
     vertex.size=2,
     vertex.color="#67a9cf",
     vertex.frame.color=NA,
     layout=layout_with_kk)


#renamed
V(iplot)$name <- str_remove_all(V(iplot)$name, pattern =  "-P.*")

#fullnames
fullnames <- data.frame(V(iplot)$name)
names(fullnames)[1] <- "iid"
names(dr_genotype2)[1] <- "iid"


# #familynames (phenotypes)
# pheno_dhfr164_all$genotype[pheno_dhfr164_all$genotype==0] <- "Reference"
# pheno_dhfr164_all$genotype[pheno_dhfr164_all$genotype==1] <- "Mixed"
# pheno_dhfr164_all$genotype[pheno_dhfr164_all$genotype==2] <- "Mutant"

familynames <- left_join(fullnames, dr_genotype2, by="iid") |> distinct(iid, .keep_all = TRUE)
familynames$yr <- str_remove_all(str_sub(familynames$iid, 4,6), "-") 
familynames$site <- str_remove_all(str_sub(familynames$iid, 1,3), "-")
familynames$site[familynames$site=="LAM"]<-"LA"
familynames$site[familynames$site=="TR"]<-"TO"
familynames$site[familynames$site=="kS"]<-"KS"


familynames$yr[familynames$yr=="00"] <- "2016"
familynames$yr[familynames$yr=="02"] <- "2017"
familynames$yr[familynames$yr=="05"] <- "2020"
familynames$yr[familynames$yr=="06"] <- "2021"
familynames$yr[familynames$yr=="07"] <- "2022"


keySites <- c("AM", "AR", "KB", "KN", "HO", "JI", "MU", "KS", "LA")
keySites2 <- c("AM", "AR", "JI", "LA")

familynames$site[!familynames$site %in% keySites2] <- NA

familynames <- familynames |> mutate(trueGeno2 = case_when(trueGeno=="Reference" & is.na(site) |
                                                             trueGeno=="PfDHFR I164L" & is.na(site) |
                                                             trueGeno== "PfDHPS A581G" & is.na(site) |
                                                             is.na(trueGeno) & is.na(site)  ~ NA,
                                                           trueGeno=="Reference" & !is.na(site) ~ "Reference",
                                                           trueGeno=="PfDHFR I164L" & !is.na(site) ~ "PfDHFR I164L",
                                                           trueGeno== "PfDHPS A581G" & !is.na(site) ~ "PfDHPS A581G",
                                                           is.na(trueGeno) & !is.na(site)  ~ NA,
))

V(iplot)$year = familynames$yr
V(iplot)$site = familynames$site
V(iplot)$genotype = familynames$trueGeno2




# Ensure necessary libraries are loaded
## create a new iplot data (network data)
 iplot_90 <- iplot
 iplot2 <- iplot_90

 V(iplot2)$year

 # iplot_75 <- iplot
 # iplot3 <- iplot_75

 
# iplot_50 <- iplot
# iplot4 <- iplot_50
 


# remove values with missing data
iplot2_clean <- igraph::delete_vertices(iplot2, which(is.na(V(iplot2)$genotype)))
iplot3_clean <- igraph::delete_vertices(iplot3, which(is.na(V(iplot3)$genotype)))
iplot4_clean <- igraph::delete_vertices(iplot4, which(is.na(V(iplot4)$genotype)))

plot(iplot2_clean,
     vertex.label=NA,
     vertex.size=3,
     vertex.color="#67a9cf",
     vertex.frame.color=NA,
     layout=layout_with_kk)

# png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/IBDNetwork.95IBD_20250715_site.png", 
#    width = 25, height = 20, units = "cm", res = 300 ) 

pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.7a.pdf", 
    width = 3.4, height = 4.5 ) 

set.seed(44)
#genotype_95_4sites <-   
  ggraph(iplot2_clean, layout = 'fr') +  # Fruchterman-Reingold layout
  geom_edge_link(aes(), size = 0.01, color = "black") +  # Edge properties
  geom_node_point(aes(shape = genotype, color = site), size = 1.5) +  # Vertex shapes and color
  scale_shape_manual(
    values = c("PfDHPS A581G" = 18, "PfDHFR I164L" = 10, "Reference" = 20)
  ) +  # Map shapes to ggplot shape codes
  scale_color_manual(
    values = c("JI" = "goldenrod", "KB" = "gold4", "AM" = "#2bce48", "AR" = "#0067a5", 
               "KS" = "#a1caf1", "HO" = "#A52A2A", "KN" = "magenta3", "MU" = "#000000", 
               "LA" = "#FF0000"),
    labels = c("JI" = "Jinja", "KB" = "Kabale", "AM" = "Amolatar", "AR" = "Arua",  
               "KS" = "Kasese", "HO" = "Hoima", "KN" = "Kanungu", "MU" = "Mubende",  
               "LA" = "Lamwo")  # Rename sites in the legend
  ) +  # Assign color to regions
  theme(legend.position = "right") +  # Position legend on the right
  guides(
    shape = guide_legend(title = "Genotype", override.aes = list(size = 5), order = 1, direction = "vertical"),  # Rename and order genotypes in the legend
    color = guide_legend(title = "Sites", override.aes = list(size = 5), order = 2, direction = "vertical")  # Title for sites legend
  ) +  # Rename color legend for sites
  theme_void() +
    theme(
      legend.key.width = unit(0.1, "cm"),
      legend.key.height = unit(0.1, "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 9, color = "black", face = "bold"),
      legend.text = element_text(size = 8, color = "black", face = "bold"),
      legend.box        = "horizontal",   # arrange legends side by side
      legend.box.just   = "center"
    )

  dev.off()
  
  #png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/IBDNetwork.95IBD_20250715_year.png", 
   #   width = 25, height = 20, units = "cm", res = 300 ) 
  pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.7b.pdf", 
      width = 3.4, height = 4.5) 
  
  set.seed(44)
  #genotype_95_4sites <-   
  ggraph(iplot2_clean, layout = 'fr') +  # Fruchterman-Reingold layout
    geom_edge_link(aes(), size = 0.01, color = "black") +  # Edge properties
    geom_node_point(aes(shape = genotype, color = year), size = 1.5) +  # Vertex shapes and color
    scale_shape_manual(
      values = c("PfDHPS A581G" = 18, "PfDHFR I164L" = 10, "Reference" = 20)
    ) +  # Map shapes to ggplot shape codes
    scale_color_manual(
      values = c("2016" = "goldenrod", "HO"  = "gold4", "2020" = "#2bce48", "2021" = "#0067a5", 
                  "KN" = "#a1caf1", "2017" = "#A52A2A", "2022" = "magenta3", "MU" = "#000000", 
                 "LA" = "#FF0000")
      #labels = c("JI" = "Jinja", "KB" = "Kabale", "AM" = "Amolatar", "AR" = "Arua",  
      #           "KS" = "Kasese", "HO" = "Hoima", "KN" = "Kanungu", "MU" = "Mubende",  
      #           "LA" = "Lamwo"
      )   +  # Assign color to regions
    guides(
      shape = guide_legend(title = "Genotype", override.aes = list(size = 5), order = 1, direction = "vertical"),  # Rename and order genotypes in the legend
      color = guide_legend(title = "Year", override.aes = list(size = 5), order = 2, direction = "vertical")  # Title for sites legend
    ) +  # Rename color legend for sites
    theme_void() +
    theme(
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.position = "bottom",
      legend.title = element_text(size = 9, color = "black", face = "bold"),
      legend.text = element_text(size = 8, color = "black", face = "bold"),      
      legend.box        = "horizontal",   # arrange legends side by side
      legend.box.just   = "center"
    )
    
  
  dev.off()
  
  ############################################################

  
  
