
##### Plotting Clonal samples
# Plot proportion of zeros by trueGeno
dr_ibdSEGT_genot8 <- read.csv( "/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/dr_ibdSEGT_genot8.csv")

dr_ibdSEGT_genot9 <- dr_ibdSEGT_genot8 %>%
  filter(!is.na(Year), !is.na(trueGeno)#, site_grp == "withinSite"
         ) %>% 
  ungroup() %>% 
  filter(!trueGeno%in%c("PfDHFR I164L mixed", "PfDHPS A581G mixed"))  %>%
  select(iid1,iid2,site,Year,trueGeno,GenPropShared,site_grp)  

dr_ibdSEGT_genot9_Clonal <- dr_ibdSEGT_genot9 %>%
  mutate(is_clonal = GenPropShared > 0.85)

names(dr_ibdSEGT_genot9_Clonal)[5]<-"Genotype"

clonal_summary <- dr_ibdSEGT_genot9_Clonal %>%
  group_by(Genotype) %>%
  summarise(prop_clonal = mean(is_clonal)) # %>% filter(Genotype%in%c("Wild type","PfDHFR 164L","PfDHPS 581G"))

clonal_summary$Genotype[clonal_summary$Genotype=="PfDHFR I164L"]<-"PfDHFR 164L"
clonal_summary$Genotype[clonal_summary$Genotype=="PfDHPS A581G"]<-"PfDHPS 581G"
clonal_summary$Genotype[clonal_summary$Genotype=="Reference"]<-"Wild type"

pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.4a.pdf",  
    width = 2, height = 4)

ggplot(clonal_summary, aes(x = Genotype, y = prop_clonal)) +
  geom_col(fill = "orange4", col = "black", alpha = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.0))
  ) +
  labs( y = "% Clonal (IBD > 0.9), all samples") +
  theme_minimal() +
  #scale_y_continuous(breaks = seq(0, 100, 25)) +
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


# Percentage clonal from Jinja
JI_clonal_by_genotype <- dr_ibdSEGT_genot9_Clonal %>%
  filter(is_clonal == 1) %>% 
  group_by(Genotype) %>%
  summarise(
    n_clonal_total = n(),
    n_clonal_jinja = sum(site == "Jinja", na.rm = T)
  ) %>%
  mutate(pct_clonal_jinja = n_clonal_jinja / n_clonal_total)


JI_clonal_by_genotype$Genotype[JI_clonal_by_genotype$Genotype=="PfDHFR I164L"]<-"PfDHFR 164L"
JI_clonal_by_genotype$Genotype[JI_clonal_by_genotype$Genotype=="PfDHPS A581G"]<-"PfDHPS 581G"
JI_clonal_by_genotype$Genotype[JI_clonal_by_genotype$Genotype=="Reference"]<-"Wild type"


pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.4b.pdf", 
    width = 2, height = 4)

ggplot(JI_clonal_by_genotype, aes(x = Genotype, y = pct_clonal_jinja)) +
  geom_col(fill = "orange4", col = "black", alpha = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.0))
  ) +
  labs(y = "% Clonal (IBD > 0.9) from Jinja") +
  theme_minimal() +
  theme_linedraw() +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 9, face = "bold",colour = "black"),
    axis.text = element_text(size = 8, face = "bold",colour = "black"),
    axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0,colour = "black"),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.1),
    panel.grid.minor.y = element_blank(),   # turn off minor gridlines
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.1),   # no vertical grid lines
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4)
    
  )
    
dev.off()


