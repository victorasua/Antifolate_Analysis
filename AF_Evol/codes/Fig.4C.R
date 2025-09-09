


dr_ibdSEGT_genot8 <- read.csv( "/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/dr_ibdSEGT_genot8.csv")


dr_ibdSEGT_genot9 <- dr_ibdSEGT_genot8 %>%
  filter(!is.na(Year), !is.na(trueGeno)#, site_grp == "withinSite"
  ) %>% 
  ungroup() %>% 
  filter(!trueGeno%in%c("PfDHFR I164L mixed", "PfDHPS A581G mixed"))  %>%
  select(iid1,iid2,site,Year,trueGeno,GenPropShared,site_grp)  

names(dr_ibdSEGT_genot9)[5]<-"Genotype"


dr_ibdSEGT_genot9$Genotype[dr_ibdSEGT_genot9$Genotype=="PfDHFR I164L"]<-"PfDHFR 164L"
dr_ibdSEGT_genot9$Genotype[dr_ibdSEGT_genot9$Genotype=="PfDHPS A581G"]<-"PfDHPS 581G"
dr_ibdSEGT_genot9$Genotype[dr_ibdSEGT_genot9$Genotype=="Reference"]<-"Wild type"

####### Ploting non clonal samples
# png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/PropGenomeShared_Non_Clonal_Samples_3.png", 
#    width = 7, height = 12, units = "cm", res = 300 )
# Step 1: Calculate medians
means <- dr_ibdSEGT_genot9 %>% filter(Genotype%in%c( "Wild type","PfDHFR 164L","PfDHPS 581G")) %>%
  filter(GenPropShared < 0.85) %>%
  select(Genotype, GenPropShared) %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("Wild type", "PfDHFR 164L", "PfDHPS 581G"))) %>%
  filter(!is.na(Genotype)) %>%
  group_by(Genotype) %>%
  summarise(mean_val = mean(GenPropShared, na.rm = TRUE)) %>%
  mutate(label = sprintf("%.3f", mean_val))

# Step 2: Define manual p-value label
p_label_df <- data.frame(
  group1 = "PfDHFR 164L",
  group2 = "Wild type",
  y.position = 0.95,              # adjust as needed
  label = "p < 0.0001"
)



pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.4c.pdf",  
    width = 3, height = 4)

# Step 3: Create plot
dr_ibdSEGT_genot9 %>% filter(Genotype%in%c( "Wild type","PfDHFR 164L","PfDHPS 581G")) %>%
  filter(GenPropShared < 0.85) %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("Wild type", "PfDHFR 164L", "PfDHPS 581G"))) %>%
  filter(!is.na(Genotype)) %>%
  ggplot(aes(x = Genotype, y = GenPropShared)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.25, alpha = 0.25, size = 0.8) +
  geom_violin(fill = "lightblue", alpha = 0.75) +
  stat_summary(fun = mean, geom = "point", color = "brown4", size = 4,
               shape = 23, fill = "yellow", alpha = 0.75) +
  geom_text(data = means,
            aes(x = Genotype, y = mean_val - 0.05, label = label),
            inherit.aes = FALSE, color = "brown2", size = 2.5, fontface = "bold") +
  #stat_pvalue_manual(p_label_df, label = "label", tip.length = 0) +
  stat_compare_means(comparisons = list(c("PfDHFR 164L", "Wild type"), c("PfDHFR 164L", "PfDHPS 581G"), c("Wild type", "PfDHPS 581G")), 
                     label = "p.signif") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.015, 0.05)) ) +
  
  labs(y = "Proportion of genome IBD") +
  theme_linedraw() +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 9, face = "bold",colour = "black"),
    axis.text = element_text(size = 8, face = "bold",colour = "black"),
    axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0,colour = "black"),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
    panel.grid.minor.y = element_blank(),   # turn off minor gridlines
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.1),   # no vertical grid lines
    panel.grid.minor.x = element_blank()
    
  )


dev.off()

