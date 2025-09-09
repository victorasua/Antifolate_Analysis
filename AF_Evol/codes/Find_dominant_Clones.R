######################################################################################
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


######. IBD ANALYSIS with dominant clones
##########################################################################################
#########################    Import required data files    ###############################


###.1 Data files for the extracted Allele depth (AD) for PfDHFR I164L position
# Data for calculating within sample allele frequency for position dhfr 164. Important to support drug resistance frequency dat
data_AD<-fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/IBC20241010_biallelicFilt_AD.tsv", header = T)
View(data_AD)
#names(data_AD)<- as.character(str_remove_all(str_sub(colnames(data_AD), 1, 9), "=.*$"))
names(data_AD)<- as.character(str_remove_all(colnames(data_AD), "=.*$"))
colnames(data_AD)[1:2] <- c("CHROM", "POS")


# data format to long form

data_AD_long <- data_AD %>%
  pivot_longer(cols = c(3:5292), names_to = "sample", values_to = "AD") %>%
  mutate(AD2 = str_remove_all(AD, "^.*.=")) %>% 
  select(CHROM, POS, sample, AD2 ) %>%
  filter(AD2!=".,.")%>%
  separate(AD2, into = c("REF", "ALT"), sep = ",", convert = TRUE)
View(data_AD_long)



#########################    FINDING SAMPLES WITH DOMINANT CLONES    ###############################

# Step 1: Compute total depth and major allele frequency (MAF)
wsaf_classified <- data_AD_long %>%
  mutate(
    total_depth = REF + ALT,
    mj_af = pmax(REF, ALT) / total_depth
  ) %>%
  filter(total_depth >= 5)  # sites with enough coverage

View(wsaf_classified)

# Step 3: Classify each sample
dominant_clone_summary <- wsaf_classified %>%
  group_by(sample) %>%
  summarise(
    total_sites = n(),
    sites_with_maf_95 = sum(mj_af >= 0.95),
    proportion_high_maf = sites_with_maf_95 / total_sites
  ) %>%
  mutate(
    dominant_clone = proportion_high_maf >= 0.95  # Step 4: ≥95% of positions have mj_af ≥0.8
  )

#dominant_clone_summary

# View samples classified as having a dominant clone
dominant_samples <- dominant_clone_summary %>%
  filter(dominant_clone == TRUE) %>% filter(total_sites >1000) # 1000-1380 positions

# Optional: View summary table
dominant_samples


# List of samples with dominant clones
dominant_samples %>% select(sample) %>% write.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/sample_names95.txt")


