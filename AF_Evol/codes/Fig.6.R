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
library(SeqArray)/
  library(vcfR)
library(ape)
library(isoRelate)
library(stringr)
library(forcats)
library(data.table)
library(stringi)
library(dplyr)
library(stringr)
library(tidyverse)

# data
dr_ibdSEGT_genot8 <- read.csv( "/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/dr_ibdSEGT_genot8.csv")

dr_ibdSEGT_genot9 <- dr_ibdSEGT_genot8 %>%
  filter(!is.na(Year), !is.na(trueGeno), site_grp == "withinSite"
  ) %>% 
  ungroup() %>% 
  filter(!trueGeno%in%c("PfDHFR I164L mixed", "PfDHPS A581G mixed"))  %>%
  select(iid1,iid2,site,Year,trueGeno,GenPropShared,site_grp)  

names(dr_ibdSEGT_genot9)[5]<-"Genotype"
dr_ibdSEGT_genot9$Genotype[dr_ibdSEGT_genot9$Genotype=="PfDHFR I164L"]<-"PfDHFR 164L"
dr_ibdSEGT_genot9$Genotype[dr_ibdSEGT_genot9$Genotype=="PfDHPS A581G"]<-"PfDHPS 581G"
dr_ibdSEGT_genot9$Genotype[dr_ibdSEGT_genot9$Genotype=="Reference"]<-"Wild type"

# Step 1: Create summary with all combinations
plot_data <- dr_ibdSEGT_genot9 %>%
  group_by(site, Genotype) %>%
  summarise(
    Proportion = mean(GenPropShared > 0.5, na.rm = TRUE),
    Count = sum(!is.na(GenPropShared)),
    .groups = "drop"
  ) %>% filter(Count > 4) |>
  # Complete missing site/genotype combos
  complete(site, Genotype, fill = list(Proportion = 0, Count = 0)) %>%
  # Calculate confidence intervals (0 if count is 0)
  mutate(
    CI_Lower = ifelse(Count > 0,
                      Proportion - 1.96 * sqrt((Proportion * (1 - Proportion)) / Count),
                      0),
    CI_Upper = ifelse(Count > 0,
                      Proportion + 1.96 * sqrt((Proportion * (1 - Proportion)) / Count),
                      0)
  ) %>%
  filter(!is.na(site), !is.na(Genotype))



pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.6.pdf", 
    width = 7, height = 3.5) 

# png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/ProportionIBD>0.5.png", 
#    width = 0.75, height = 3.5, units = "cm", res = 300 ) 

ggplot(plot_data, aes(x = site, y = Proportion, fill = Genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, col = "black") +
  scale_fill_manual( values = c("#E69F00",  "#938197FF", "#377EB8"),
                     na.translate = F,
                     name = "Genotype group",
                     guide = guide_legend(
                       direction = "horizontal",
                       label.position = "bottom",
                       title.position = "top",
                       title.hjust = 0.5,
                       nrow = 1 ) ) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.0, 0.0)) ) +
  geom_errorbar(
    aes(ymin = CI_Lower, ymax = CI_Upper),
    position = position_dodge(width = 0.8),
    width = 0.1
  ) +
  theme_minimal() +
  labs(
    y = "Proportion with IBD > 0.5",
    fill = "Genotype group"
  ) + 
  
  theme(#axis.ticks = element_blank(),
    legend.key.width = unit(0.35, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.text = element_text(size = 8, hjust = 0.5, vjust = 1, face = "bold",colour = "black"),
    legend.title = element_text(size = 8, face = "bold",colour = "black"),
    legend.position = "bottom",
    legend.box.margin = margin(t = -12, b = 0, unit = "pt"),
    legend.box = "horizontal", 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9, face = "bold",colour = "black"),
    axis.text = element_text(size = 8, face = "bold",colour = "black"),
    axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
    panel.grid.minor.y = element_blank(),   # turn off minor gridlines
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.1),   # no vertical grid lines
    panel.grid.minor.x = element_blank()
    
    )

dev.off()


########## Table summary: Proportion wiht IBD > 0.5
proportion_0.5CI_Genotype <- dr_ibdSEGT_genot9 %>%
  group_by(site, Genotype) %>%
  summarise(
    Proportion = mean(GenPropShared > 0.5, na.rm = TRUE),
    Count = sum(!is.na(GenPropShared)),
    CI_Lower = Proportion - 1.96 * sqrt((Proportion * (1 - Proportion)) / Count),
    CI_Upper = Proportion + 1.96 * sqrt((Proportion * (1 - Proportion)) / Count),
    .groups = "drop"
  )


proportion_0.5CI_Genotype_DF <- data.frame(proportion_0.5CI_Genotype) 
proportion_0.5CI_Genotype_DF2 <- proportion_0.5CI_Genotype_DF %>%
 # filter(Count > 4) |>
  mutate(
    Proportion = formatC(Proportion, format = "f", digits = 2),
    Count = round(Count),  # If Count is not already an integer, it will be rounded to nearest integer
    CI_Lower = formatC(CI_Lower, format = "f", digits = 2),
    CI_Upper = formatC(CI_Upper, format = "f", digits = 2),
    Prop_CI_N = paste0(Proportion, " (", CI_Lower, "-", CI_Upper, ")", ", ", Count) ) |> 
  select(site, Genotype, Prop_CI_N) |> 
  pivot_wider(names_from = Genotype, values_from = Prop_CI_N) 





