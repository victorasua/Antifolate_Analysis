#########################    load the following libraries.   #########################

# Install most recent released version
# install.packages(
# 'coiaf',
# repos=c("https://plasmogenepi.r-universe.dev", "https://cloud.r-project.org") )

library(stringr)
library(dplyr)
library(ggplot2)
library(ggplotify)
library(tidyr)
library(data.table)
library(purrr)
library(ggpubr)
library(coiaf)

##########################################################################################
#########################    Import required data files    ###############################


###.1 Data files for the extracted Allele depth (AD) for PfDHFR I164L position
# Data for calculating within sample allele frequency for position dhfr 164. Important to support drug resistance frequency data


prx00_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx00_dhfr_164_AD.csv")))
prx00_dhfr_164_AD<-prx00_dhfr_164_AD[-c(1:2),]%>%data.frame()

prx02_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx02_dhfr_164_AD.csv")))
prx02_dhfr_164_AD<-prx02_dhfr_164_AD[-c(1:2),]%>%data.frame()

prx03_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx03_dhfr_164_AD.csv")))
prx03_dhfr_164_AD<-prx03_dhfr_164_AD[-c(1:2),]%>%data.frame()

prx04_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx04_dhfr_164_AD.csv")))
prx04_dhfr_164_AD<-prx04_dhfr_164_AD[-c(1:2),]%>%data.frame()

prx05_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx05_dhfr_164_AD.csv")))
prx05_dhfr_164_AD<-prx05_dhfr_164_AD[-c(1:2),]%>%data.frame()

prx06_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx06_dhfr_164_AD.csv")))
prx06_dhfr_164_AD<-prx06_dhfr_164_AD[-c(1:2),]%>%data.frame()

prx07_dhfr_164_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx07_dhfr_164_AD.csv")))
prx07_dhfr_164_AD<-prx07_dhfr_164_AD[-c(1:2),]%>%data.frame()

# Combine AD data files
AD_data_dhfr<-rbind(prx00_dhfr_164_AD,
                    prx02_dhfr_164_AD,
                    prx03_dhfr_164_AD,
                    prx04_dhfr_164_AD,
                    prx05_dhfr_164_AD,
                    prx06_dhfr_164_AD,
                    prx07_dhfr_164_AD)

colnames(AD_data_dhfr)<-"id_AD"

## Split string to seprate sample IDs from read-depth
AD_data_dhfr_df <- cbind(str_remove_all(str_sub(AD_data_dhfr$id_AD, 1, 18), "="),  
                         str_sub(AD_data_dhfr$id_AD, 19) ) %>% data.frame()

prx_AD_dhfr_df2<-AD_data_dhfr_df%>%filter(X2!=".") %>% 
  filter(!str_detect(X1, ",0")) %>%
  filter(!str_detect(X1, "\\.$")) %>%
  filter(!str_detect(X2, "=."))
prx_AD_dhfr_df2$X2<-str_remove_all(prx_AD_dhfr_df2$X2, "=")

## Add Ref and Alt read-depths
prx_AD_dhfr_df2$Ref<-str_remove_all(string=prx_AD_dhfr_df2$X2, pattern = ",.*$") 
prx_AD_dhfr_df2$Alt<-str_remove_all(string=prx_AD_dhfr_df2$X2, pattern = "^.*,") 


## Turn Read_depths to numeric values
prx_AD_dhfr_df2$Ref1<-as.numeric(prx_AD_dhfr_df2$Ref)
prx_AD_dhfr_df2$Alt1<-as.numeric(prx_AD_dhfr_df2$Alt)


## GenotypeCalling for PfDHFR 164
prx_AD_WSAF_Genotype_dhfr <- prx_AD_dhfr_df2 %>%
  mutate(total_rd=prx_AD_dhfr_df2$Ref1+prx_AD_dhfr_df2$Alt1) %>%
  mutate(REF=prx_AD_dhfr_df2$Ref1 / total_rd,
         ALT=prx_AD_dhfr_df2$Alt1 / total_rd,
         Ref_allele=case_when(Ref1>4~"1",
                              Ref1<=4~"0"),
         Alt_allele=case_when(Alt1>4~"1",
                              Alt1<=4~"0"),
         PfDHFR_164_genotype=case_when(Ref_allele=="1"&Alt_allele=="0"~"0",
                                       Ref_allele=="1"&Alt_allele=="1"~"1",
                                       Ref_allele=="0"&Alt_allele=="1"~"2",)) %>%
  filter(total_rd>8)

prx_AD_WSAF_Genotype_dhfr$X1<-str_remove_all(prx_AD_WSAF_Genotype_dhfr$X1, "-P.*$")
names(prx_AD_WSAF_Genotype_dhfr)[1]<-"id"
names(prx_AD_WSAF_Genotype_dhfr)[2]<-"AD"



#####################################################################################
#### WSAF for PfDHFR 164L (mixed and mutant)
#wasf_164L
wsaf_pfdhfr164 <- prx_AD_WSAF_Genotype_dhfr %>%
  select(id,REF,Alt1,ALT,PfDHFR_164_genotype) %>%
  drop_na() %>%
  filter(PfDHFR_164_genotype==2|PfDHFR_164_genotype==1) %>%
  filter(Alt1>4) 

####plot 164 alone
pfdhfr164 <- wsaf_pfdhfr164 %>%
  select(id,REF,Alt1,ALT,PfDHFR_164_genotype) %>%
  drop_na() %>%
  filter(PfDHFR_164_genotype==2|PfDHFR_164_genotype==1) %>%
  filter(Alt1>4) %>%
  mutate(WSAF=case_when(ALT>0&ALT<0.05~"0-5",
                        ALT>0.05&ALT<0.2010~"5-20",
                        ALT>0.2011&ALT<0.4010~"21-40",
                        ALT>0.4011&ALT<0.6010~"41-60",
                        ALT>0.6011&ALT<0.8010~"61-80",
                        ALT>0.8011&ALT<1.00~"81-99",
                        ALT==1.00~"100"))%>%
  mutate(WSAF=factor(WSAF,level=c("0-5","5-20","21-40","41-60","61-80","81-99","100"))) %>%
  
  ggplot(aes(x = WSAF, fill = factor(PfDHFR_164_genotype)))+
  geom_bar(stat = "count", position = "stack") +
  labs(#title="PfDHFR 164L",
    x="WSAF(%)",
    y="Number of Samples")+
  scale_fill_manual(values = c("1" = "#1f77b4", "2" = "#d95f02"), 
                    name = "Allele:",
                    labels = c("Mixed", "Mutant")) +
  theme_pubclean() +
  theme(#panel.grid.major = element_line(size=0.1, colour = "gray"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  coord_flip()





##############@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PfDHPS 581G @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#############################
###.2 Data files for the extracted Allele depth (AD) for PfDHPS 518G position
# Data for calculating within sample allele frequency for position dhfr 164. Important to support drug resistance frequency data


prx00_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx00_dhps_518_AD.csv")))
prx00_dhps_518_AD<-prx00_dhps_518_AD[-c(1:2),]%>%data.frame()

prx02_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx02_dhps_518_AD.csv")))
prx02_dhps_518_AD<-prx02_dhps_518_AD[-c(1:2),]%>%data.frame()

prx03_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx03_dhps_518_AD.csv")))
prx03_dhps_518_AD<-prx03_dhps_518_AD[-c(1:2),]%>%data.frame()

prx04_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx04_dhps_518_AD.csv")))
prx04_dhps_518_AD<-prx04_dhps_518_AD[-c(1:2),]%>%data.frame()

prx05_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx05_dhps_518_AD.csv")))
prx05_dhps_518_AD<-prx05_dhps_518_AD[-c(1:2),]%>%data.frame()

prx06_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx06_dhps_518_AD.csv")))
prx06_dhps_518_AD<-prx06_dhps_518_AD[-c(1:2),]%>%data.frame()

prx07_dhps_518_AD<-data.frame(t(fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/AD_out/prx07_dhps_518_AD.csv")))
prx07_dhps_518_AD<-prx07_dhps_518_AD[-c(1:2),]%>%data.frame()

# Combine AD data files
AD_data_dhps<-rbind(prx00_dhps_518_AD,
                    prx02_dhps_518_AD,
                    prx03_dhps_518_AD,
                    prx04_dhps_518_AD,
                    prx05_dhps_518_AD,
                    prx06_dhps_518_AD,
                    prx07_dhps_518_AD)

colnames(AD_data_dhps)<-"id_AD"

## Split string to seprate sample IDs from read-depth
prx_AD_dhps_df<-cbind(str_remove_all(str_sub(AD_data_dhps$id_AD, 1, 18), "="),  
                      str_sub(AD_data_dhps$id_AD, 19) ) %>%data.frame()

prx_AD_dhsp_df2<-prx_AD_dhps_df%>%filter(X2!=".") %>% 
  filter(!str_detect(X1, ",0")) %>%
  filter(!str_detect(X1, "\\.$")) %>%
  filter(!str_detect(X2, "=."))
prx_AD_dhsp_df2$X2<-str_remove_all(prx_AD_dhsp_df2$X2, "=")

## Add Ref and Alt read-depths
prx_AD_dhsp_df2$Ref<-str_remove_all(string=prx_AD_dhsp_df2$X2, pattern = ",.*$") 
prx_AD_dhsp_df2$Alt<-str_remove_all(string=prx_AD_dhsp_df2$X2, pattern = "^.*,") 


## Turn Read_depths to numeric values
prx_AD_dhsp_df2$Ref1<-as.numeric(prx_AD_dhsp_df2$Ref)
prx_AD_dhsp_df2$Alt1<-as.numeric(prx_AD_dhsp_df2$Alt)


## GenotypeCalling for PfDHPS 581G
prx_AD_WSAF_Genotype_dhps <- prx_AD_dhsp_df2 %>%
  mutate(total_rd=prx_AD_dhsp_df2$Ref1+prx_AD_dhsp_df2$Alt1) %>%
  mutate( REF=prx_AD_dhsp_df2$Ref1 / total_rd,
          ALT=prx_AD_dhsp_df2$Alt1 / total_rd,
          Ref_allele=case_when(Ref1>4~"1",
                               Ref1<=4~"0"),
          Alt_allele=case_when(Alt1>4~"1",
                               Alt1<=4~"0"),
          PfDHPS_581_genotype=case_when(Ref_allele=="1"&Alt_allele=="0"~"0",
                                        Ref_allele=="1"&Alt_allele=="1"~"1",
                                        Ref_allele=="0"&Alt_allele=="1"~"2",)) %>%
  filter(total_rd>=10)


prx_AD_WSAF_Genotype_dhps$X1<-str_remove_all(prx_AD_WSAF_Genotype_dhps$X1, "-P.*$")
names(prx_AD_WSAF_Genotype_dhps)[1]<-"id"
names(prx_AD_WSAF_Genotype_dhps)[2]<-"AD"



#####################################################################################
#### WSAF for PfDHPS 581G (mixed and mutant)


pfdhps581<-prx_AD_WSAF_Genotype_dhps%>%
  select(id,REF,Alt1,ALT,PfDHPS_581_genotype)%>%
  drop_na()%>%
  filter(PfDHPS_581_genotype==2|PfDHPS_581_genotype==1) %>%
  filter(Alt1>4) %>%
  mutate(WSAF=case_when(ALT>0&ALT<0.05~"0-5",
                        ALT>0.05&ALT<0.2010~"5-20",
                        ALT>0.2011&ALT<0.4010~"21-40",
                        ALT>0.4011&ALT<0.6010~"41-60",
                        ALT>0.6011&ALT<0.8010~"61-80",
                        ALT>0.8011&ALT<1.00~"81-99",
                        ALT==1.00~"100"))%>%
  mutate(WSAF=factor(WSAF,level=c("0-5","5-20","21-40","41-60","61-80","81-99","100")))%>%
  
  ggplot(aes(x = WSAF, fill = factor(PfDHPS_581_genotype)))+
  geom_bar(stat = "count", position = "stack") +
  labs(#title="PfDHPS 518G",
    x="",
    y="Number of Samples")+
  scale_fill_manual(values = c("1" = "#1f77b4", "2" = "#d95f02"), 
                    name = "Allele:",
                    labels = c("Mixed", "Mutant")) +
  theme_pubclean() +
  theme(#panel.grid.major = element_line(size=0.1, colour = "gray"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"))+
  coord_flip()

library(cowplot)
####### Ploting combined pfdhfr and pfdhps read depths
png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/WSAF_164L_581G.png", 
    width = 15, height = 10, units = "cm", res = 300 )

ggdraw() +
  draw_plot(pfdhfr164, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(pfdhps581, x = 0.5, y = 0, width = 0.5, height = 1) + # Add a frame around Cluster_10
  annotate("text", x = 0.05, y = 0.95, label = "A", size = 3.5, color = "black", fontface = "bold") + # Pfdhfr 164
  annotate("text", x = 0.575, y = 0.95, label = "B", size = 3.5, color = "black", fontface = "bold") # PfDHPS 581


dev.off()


####### COmbine WSAF for 164L and 581G
#wasf_164L
wsaf_pfdhfr164 <- prx_AD_WSAF_Genotype_dhfr %>%
  mutate(total_rd=as.numeric(prx_AD_WSAF_Genotype_dhfr$Ref)+as.numeric(prx_AD_WSAF_Genotype_dhfr$Alt)) %>%
  select(id, total_rd, Alt1, ALT, PfDHFR_164_genotype) %>%
  drop_na() %>%
  filter(PfDHFR_164_genotype==2|PfDHFR_164_genotype==1) %>%
  #filter(Ref==0) %>%
  mutate(gene="PfDHFR 164L", genotype=PfDHFR_164_genotype)

#wasf_581G
wsaf_pfdhps581 <- prx_AD_WSAF_Genotype_dhps%>%
  mutate(total_rd=as.numeric(prx_AD_WSAF_Genotype_dhps$Ref)+as.numeric(prx_AD_WSAF_Genotype_dhps$Alt)) %>%
  select(id, total_rd, Alt1, ALT,PfDHPS_581_genotype) %>%
  drop_na()%>%
  filter(PfDHPS_581_genotype==2|PfDHPS_581_genotype==1) %>%
  mutate(gene="PfDHPS 581G", genotype=PfDHPS_581_genotype)

gene_labels <- c(
  "PfDHFR 164L" = "PfDHFR 164L",
  "PfDHPS 581G" = "PfDHPS 581G"
)

## combine
# pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/WSAF_164L_581G_boxplot_3.pdf", 
#    width = 12, height = 14, units = "cm", res = 300 )

pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/WSAF_164L_581G_boxplot_3.pdf", 
    width = 3, height = 5)

bind_rows(wsaf_pfdhfr164, wsaf_pfdhps581) %>%
  select( id, total_rd, Alt1, ALT, gene, genotype) %>%
  mutate(genotype = factor(genotype, levels = c(1, 2), labels = c("Mixed", "Mutant"))) %>%
  ggplot(aes(x = genotype, y = ALT)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.7, color = "black") +   # dot plot
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "lightblue2", width = 0.4) +  # boxplot
  facet_wrap(~gene,
             labeller = labeller(gene = gene_labels)
  ) +
  labs(x = "Genotype", y = "WSAF") +
  theme_minimal() +
  theme(#panel.grid.major = element_line(size=0.1, colour = "gray"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.title = element_text(face = "bold", size = 9, colour = "black"),
    axis.text = element_text(face = "bold", size = 8, colour = "black"),
    strip.text = element_text(face = "bold", size = 9, colour = "black"),
    panel.grid.minor.y = element_blank(),   # turn off minor gridlines
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.1),   # no vertical grid lines
    panel.grid.minor.x = element_blank()
  )


dev.off()

