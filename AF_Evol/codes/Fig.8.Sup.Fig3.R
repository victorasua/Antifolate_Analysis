
# packages
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(ggplotify)
library(ggplot2)
library(data.table)
library(stringr)
library(stringi)
library(tidyr)
library(Cairo)
library(cluster)



setwd("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/haplotypeReanalysis")

# Meta data
metdat<-fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/PCA/gen_dat.csv")
metdat<-metdat %>% filter(!is.na(genotype)) %>% 
  filter(aa_change=="Ile164Leu"|aa_change=="Ala581Gly")%>%
  mutate(allele=case_when(aa_change=="Ile164Leu"&genotype==0~"PfDHFR I164",
                          aa_change=="Ile164Leu"&genotype==2~"PfDHFR 164L",
                          aa_change=="Ile164Leu"&genotype==1~"PfDHFR I164L",
                          aa_change=="Ala581Gly"&genotype==0~"PfDHPS A581",
                          aa_change=="Ala581Gly"&genotype==2~"PfDHPS 581G",
                          aa_change=="Ala581Gly"&genotype==1~"PfDHPS A581G"))
metdat<-data.frame(metdat)
metdat$year<-str_remove_all(str_sub(metdat$id, 4, 6), "-")

yr<-c("07"="2022",
      "06"="2021",
      "05"="2020",
      "04"="2019",
      "03"="2018",
      "02"="2017",
      "00"="2016")

sites<-c("TO"="Tororo","AR"="Arua","MU"="Mubende","JI"="Jinja","KB"="Kabale","AG"="Agago","LA"="Lamwo","KO"="Kole","AM"="Amolatar",
         "KN"="Kanungu","HO"="Hoima","KBG"="Kaabong","KBK"="Koboko","KS"="Kasese","KTK"="Katakwi","KAP"="Kapchorwa")


regions<-c("Agago"="North","Amolatar"="North","Arua"="N.West","Koboko"="N.West","Hoima"="West","Kapchorwa"="East","Kabale"="S.West",
           "Kabale"="S.West","Kanungu"="S.West","Kaabong"="North","Kole"="North","Kasese"="S.West","Katakwi"="East","Lamwo"="North",
           "Mubende"="Central","Tororo"="East","Jinja"="E.Central")


# Replace sites and year values
metadata <- metdat %>%
  mutate(year = recode(year, !!!yr)) %>% select(id, site, year, allele)
metadata <- metadata %>%
  mutate(site = recode(site, !!!sites)) 



# # Genotype data chr4
# gdat_chr4<-fread("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/haplotypeReanalysis/DRHAP.chr4_GenotypeTable_ReCoded.tsv", 
#                  sep = "\t")

gdat_chr4a <- read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/20250604_Reanalysis/haplotypeReanalysis/DRHAP.chr4_GenotypeTable_ReCoded.tsv",
                       sep = "\t", 
                       check.names = FALSE)

#gdat_chr8b <- read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/20250604_Reanalysis/haplotypeReanalysis/DRHAP.chr8_GenotypeTable_ReCoded.tsv",
#                    sep = "\t")


gdat_chr4 <- gdat_chr4a %>% 
  select(-c(REF,ALT,QUAL)) %>% t()
colnames(gdat_chr4) <- gdat_chr4[2, ]
gdat_chr4 <- data.frame(gdat_chr4[-c(1:3), ])


gdat_chr4$id <- str_remove_all(rownames(gdat_chr4), pattern =  "-P.*") 
gdat_chr4 <- gdat_chr4 %>%
  select(id, everything())
rownames(gdat_chr4)<-NULL


id_excld<-c("3D7","7G8","V1S","DD2","HB3","HC","LC","3d7","Dd2","NTC","PM2")
gdat_chr4_f<-gdat_chr4[!gdat_chr4$id%in%id_excld, ]

# combine Meta data and Genotype data
gdat_metdat<-left_join(gdat_chr4_f,metadata,by="id") %>% select(id, site, year, allele, everything()) %>% drop_na()


# Add region
gdat_metdat$region<-gdat_metdat$site
gdat_metdat <- gdat_metdat %>%
  mutate(region = recode(region, !!!regions)) %>% select(id, site, year, allele, region, everything())


############################# Drop polyclonal Samples #################################

DominatCloneLIst<-read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/sample_names95.txt")
DominatCloneLIst<-str_remove_all(DominatCloneLIst$sample, "-P.*$")
gdat_metdat_Domnt<-gdat_metdat[gdat_metdat$id%in%DominatCloneLIst, ]

############################# Drop polyclonal Samples #################################
# # Convert relevant columns (from 5th to last) to numeric
# gdat_metdat[, 6:ncol(gdat_metdat)] <- lapply(gdat_metdat[, 6:ncol(gdat_metdat)], as.numeric)
# 
# # Find IDs that have any value of 1 across any of the columns
# ids_with_1. <- unique(gdat_metdat$id[apply(gdat_metdat[, 6:ncol(gdat_metdat)], 1, function(x) any(x == 1))])
# 
# 
# # Drop all rows with those IDs
# gdat_metdat_filtered <- gdat_metdat[!gdat_metdat$id %in% ids_with_1., ]



# Position of 164 , 748576


################################################################################
#####@@@@@@@@@ FUNCTION TO GENERATE HEATMAPS FOR DHFR 164L @@@@@@@@@##########
# generate_dhfr_heatmap <- function(data, region_name = NULL, heatmap_name = NULL, col_annotation_colors = NULL, plot = TRUE) {
#   
#   if (is.null(heatmap_name)) {
#     heatmap_name <- paste0("result_", region_name)
#   }


############################  dhfr 164L #############################

  # Step 1: Filter data by allele and region
  data_filtered <- gdat_metdat_Domnt %>%
    filter(
      allele %in% c("PfDHFR 164L")#,
     # region == region_name
    )

# Sample data
set.seed(123)  # Optional: for reproducibility
#sampled_data <- data_filtered[sample(nrow(data_filtered), 200), ]
sampled_data <- data_filtered

  # Step 2: Split genotype and metadata
  genotype_data <- sampled_data[, 6:ncol(sampled_data)]
  genotype_data <- data_filtered
  rownames(genotype_data) <- genotype_data$id
  # Filter ranges 
  # 748576 + 50000
  # 748576 - 50000
  genotype_25kb <- genotype_data[, 144:418]  # 50kb region window
  genotype_matrix <- as.matrix(genotype_25kb)
  genotype_matrix_2 <- as.matrix (genotype_data)
  
  # Change to numeric
  genotype_matrix_num <- matrix(
    as.numeric(genotype_matrix_2),
    nrow = nrow(genotype_matrix_2),
    ncol = ncol(genotype_matrix_2),
    dimnames = dimnames(genotype_matrix_2)
  )[, 6:ncol(genotype_matrix_2)]
  
  
  # remove non variant columns
  # Threshold
  threshold <- 0.95
  
  # Function to check if column has >=80% same values
  is_low_variation <- function(x) {
    prop_max <- max(table(x)) / length(x)
    return(prop_max >= threshold)
  }
  
  # Apply across columns and keep those that do NOT meet condition
  variant_genotype_matrix <- genotype_matrix_num[, !apply(genotype_matrix_num, 2, is_low_variation)]
  
  # Check dimensions
  dim(genotype_matrix_num)
  dim(variant_genotype_matrix)
  
  
  # Step 3: Prepare metadata
  metadata <- data_filtered[, 1:5]
  metadata <- metadata[order(metadata$id, decreasing = FALSE), ]
  
  px <- data.frame(
    Allele = metadata$allele,
    Year = metadata$year,
    #Region = as.factor(metadata$region),
    Site = as.factor(metadata$site)
  )
  
  # Step 4: Default color scheme (can be overridden)
  col_annotation_colors <- NULL
  
  if (is.null(col_annotation_colors)) {
    col_annotation_colors <- list(
      # Region = c(
      #   "North" = "#A65628",
      #   "N.West" = "#F0E442",
      #   "West" = "#7570B3",
      #   "E.Central" = "#E7298A",
      #   "East" = "#66A61E",
      #   "S.West" = "#D95F02",
      #   "Central" = "#1F78B4"
      # ) ,

      # "light_peach"   = "#FFE5B4",
      # "apricot"       = "#FDB77E",
      # "light_orange"  = "#FFA64C",
      # "vivid_orange"  = "#FF7F00",
      # "pumpkin"       = "#FF7518",
      # "orange"        = "#FFA500",
      # "burnt_orange"  = "#CC5500",
      # "deep_orange"   = "#993300"
      # 
        Site = c(
          "Agago" = "#D79A6D",
          "Amolatar" = "#4D470F",
          "Arua" = "#9C931D",
          "Hoima" = "#A01F66",
          "Jinja" = "#FF7F00" ,
          "Kabale" = "#39355A",
          "Kaabong" = "#1F78B4",
          "Kanungu" = "#7570B3",
          "Kole"  =  "#FBF9C4" ,
          "Kasese" = "#D1CFEA",
          "Katakwi" = "#2F4F10",
          "Lamwo"  =  "#6F3B1A",
          "Mubende" = "#F7C2DA",
          "Tororo" = "#A7CA5D"
        ),
      Allele = c(
        "PfDHFR I164" = "#75B7D1",
        "PfDHFR 164L" =  "#FFE5B4",
        "PfDHPS 581G" = "#2D5765"
      ),
      Year = c(
        "2016" = "#D3DAD0",
        "2017" = "#B2BEB5",
        "2020" = "#9CA59F",
        "2021" = "#868F89",
        "2022" = "#5A635E"
      )
    )
  }
  
  # Step 5: Annotation setup
  hk <- rowAnnotation(
    df = px,
    col = col_annotation_colors,
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 10),
    annotation_name_rot = 90,
    annotation_legend_param = list(
      legend_height = unit(12, "cm"),
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    ),
    border = TRUE
  )
  
  # Step 6: Color palette

  pal <- colorRamp2(c(0, 1, 2), c("#4d194d", "#b980b9",  "#f1e5f1"))
  
  ### remove non-variant columns 
  # Remove columns where all values are the same
  #$ variant_genotype_matrix <- genotype_matrix[, apply(genotype_matrix, 2, function(x) length(unique(x)) > 1)]
  
  # Check dimensions before and after
  # dim(genotype_matrix)
  # dim(variant_genotype_matrix)
  
  
# Step 7: Generate heatmap
dhfr <- 
  Heatmap(
        #variant_genotype_matrix,
    genotype_matrix_num,
        col = pal,
        row_dend_side = "left",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        name = "heatmap_name",  # <-- dynamic name
        right_annotation = hk,
        row_names_gp = gpar(fontsize = 0),
        column_names_gp = gpar(fontsize = 0),
        column_names_rot = 45,
        width = 0.5,
        show_heatmap_legend = T,
        heatmap_legend_param = list(
          title = "Genotype", 
          at = c(0, 2),
          labels = c("REF", "ALT"),
          color_bar = "discrete",
          direction = "horizontal",               
          title_gp = gpar(fontsize= 9),
          labels_gp = gpar(fontsize = 8)
        ),
        row_dend_width = unit(14, "mm"),
        row_gap = unit(0.3, "mm"),
        border = TRUE
      )



# pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Haplotype_sharing_dhfr-20250730.pdf", 
#    width = 17.5, height = 10, units = "cm", res = 450 ) 

#pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Haplotype_sharing_dhfr-20250905.pdf", 
#    width = 7, height = 4) 

#pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.8a.pdf", 
 #   width = 7, height = 4.5) 
pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Sup.Fig.3a.pdf", 
    width = 7, height = 4.5) 

draw(
  dhfr,     
  ht_gap = unit(rep(1, 6), "mm"),
  annotation_legend_side = "right" ,
  heatmap_legend_side = "right"
)

col_index <- which(colnames( dhfr@matrix ) == "X748576")

decorate_heatmap_body("heatmap_name", {
  grid.lines(
    x = unit(c(col_index, col_index) / ncol(dhfr@matrix), "npc") ,
    y = unit(c(0, 1), "npc") ,
    gp = gpar(col = "red3", lwd = 1.5)
  )
})


dev.off() 





############################***2**** dhfr I164 #############################

# Step 1: Filter data by allele and region
data_filtered <- gdat_metdat_Domnt %>%
  filter(
    allele %in% c("PfDHFR I164")#,
    # region == region_name
  )

# Sample data
set.seed(123)  # Optional: for reproducibility
sampled_data <- data_filtered[sample(nrow(data_filtered), 200), ]
#sampled_data <- data_filtered

# Step 2: Split genotype and metadata
genotype_data <- sampled_data[, 6:ncol(sampled_data)]
#genotype_data <- data_filtered
rownames(genotype_data) <- genotype_data$id
# Filter ranges 
# 748576 + 50000
# 748576 - 50000
genotype_25kb <- genotype_data[, 144:418]  # 50kb region window
genotype_matrix <- as.matrix(genotype_25kb)
genotype_matrix_2 <- as.matrix (genotype_data)

# Change to numeric
genotype_matrix_num <- matrix(
  as.numeric(genotype_matrix_2),
  nrow = nrow(genotype_matrix_2),
  ncol = ncol(genotype_matrix_2),
  dimnames = dimnames(genotype_matrix_2)
)[, 6:ncol(genotype_matrix_2)]


# remove non variant columns
# Threshold
threshold <- 0.95

# Function to check if column has >=80% same values
is_low_variation <- function(x) {
  prop_max <- max(table(x)) / length(x)
  return(prop_max >= threshold)
}

# Apply across columns and keep those that do NOT meet condition
variant_genotype_matrix <- genotype_matrix_num[, !apply(genotype_matrix_num, 2, is_low_variation)]

# Check dimensions
dim(genotype_matrix_num)
dim(variant_genotype_matrix)


# Step 3: Prepare metadata
metadata <- sampled_data[, 1:5]
metadata <- metadata[order(metadata$id, decreasing = FALSE), ]

px <- data.frame(
  Allele = metadata$allele,
  Year = metadata$year,
  #Region = as.factor(metadata$region),
  Site = as.factor(metadata$site)
)

# Step 4: Default color scheme (can be overridden)
col_annotation_colors <- NULL

if (is.null(col_annotation_colors)) {
  col_annotation_colors <- list(
    # Region = c(
    #   "North" = "#A65628",
    #   "N.West" = "#F0E442",
    #   "West" = "#7570B3",
    #   "E.Central" = "#E7298A",
    #   "East" = "#66A61E",
    #   "S.West" = "#D95F02",
    #   "Central" = "#1F78B4"
    # ) ,
    
    # "light_peach"   = "#FFE5B4",
    # "apricot"       = "#FDB77E",
    # "light_orange"  = "#FFA64C",
    # "vivid_orange"  = "#FF7F00",
    # "pumpkin"       = "#FF7518",
    # "orange"        = "#FFA500",
    # "burnt_orange"  = "#CC5500",
    # "deep_orange"   = "#993300"
    # 
    Site = c(
      "Agago" = "#D79A6D",
      "Amolatar" = "#4D470F",
      "Arua" = "#9C931D",
      "Hoima" = "#A01F66",
      "Jinja" = "#FF7F00" ,
      "Kabale" = "#39355A",
      "Kaabong" = "#1F78B4",
      "Kanungu" = "#7570B3",
      "Kole"  =  "#FBF9C4" ,
      "Kasese" = "#D1CFEA",
      "Katakwi" = "#2F4F10",
      "Lamwo"  =  "#6F3B1A",
      "Mubende" = "#F7C2DA",
      "Tororo" = "#A7CA5D",
      "Kapchorwa" = "#FBF9C4",
      "Koboko" = "#9C931D"
        
    ),
    Allele = c(
      "PfDHFR I164" = "#75B7D1",
      "PfDHFR 164L" =  "#FFE5B4",
      "PfDHPS 581G" = "#2D5765"
    ),
    Year = c(
      "2016" = "#D3DAD0",
      "2017" = "#B2BEB5",
      "2020" = "#9CA59F",
      "2021" = "#868F89",
      "2022" = "#5A635E"
    )
  )
}

# Step 5: Annotation setup
hk <- rowAnnotation(
  df = px,
  col = col_annotation_colors,
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 10),
  annotation_name_rot = 90,
  annotation_legend_param = list(
    legend_height = unit(12, "cm"),
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  ),
  border = TRUE
)

# Step 6: Color palette

pal <- colorRamp2(c(0, 1, 2), c("#4d194d", "#b980b9",  "#f1e5f1"))

### remove non-variant columns 
# Remove columns where all values are the same
#$ variant_genotype_matrix <- genotype_matrix[, apply(genotype_matrix, 2, function(x) length(unique(x)) > 1)]

# Check dimensions before and after
# dim(genotype_matrix)
# dim(variant_genotype_matrix)


# Step 7: Generate heatmap
dhfr <- 
  Heatmap(
    #variant_genotype_matrix,
    genotype_matrix_num,
    col = pal,
    row_dend_side = "left",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    name = "heatmap_name",  # <-- dynamic name
    right_annotation = hk,
    row_names_gp = gpar(fontsize = 0),
    column_names_gp = gpar(fontsize = 0),
    column_names_rot = 45,
    width = 0.5,
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      title = "Genotype", 
      at = c(0, 2),
      labels = c("REF", "ALT"),
      color_bar = "discrete",
      direction = "horizontal",               
      title_gp = gpar(fontsize= 9),
      labels_gp = gpar(fontsize = 8)
    ),
    row_dend_width = unit(14, "mm"),
    row_gap = unit(0.3, "mm"),
    border = TRUE
  )



# pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Haplotype_sharing_dhfr-20250730.pdf", 
#    width = 17.5, height = 10, units = "cm", res = 450 ) 

#pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Haplotype_sharing_dhfr-20250905.pdf", 
#    width = 7, height = 4) 

#pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.8a.pdf", 
#   width = 7, height = 4.5) 
pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Suppl.Fig.3a.pdf", 
    width = 7, height = 4.5) 

draw(
  dhfr,     
  ht_gap = unit(rep(1, 6), "mm"),
  annotation_legend_side = "right" ,
  heatmap_legend_side = "right"
)

col_index <- which(colnames( dhfr@matrix ) == "X748576")

decorate_heatmap_body("heatmap_name", {
  grid.lines(
    x = unit(c(col_index, col_index) / ncol(dhfr@matrix), "npc") ,
    y = unit(c(0, 1), "npc") ,
    gp = gpar(col = "red3", lwd = 1.5)
  )
})


dev.off() 





####################### DHPS PLOT HEATMAPS #######################


#######################################################################

# Genotyping data chr8
gdat_chr8a <- read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/20250604_Reanalysis/haplotypeReanalysis/DRHAP.chr8_GenotypeTable_ReCoded.tsv",
                       sep = "\t", 
                       check.names = FALSE)

gdat_chr8<-gdat_chr8a %>% 
  select(-c(REF,ALT,QUAL)) %>% t()
colnames(gdat_chr8) <- gdat_chr8[2, ]
gdat_chr8 <- data.frame(gdat_chr8[-c(1:3), ])


gdat_chr8$id <- str_remove_all(rownames(gdat_chr8), pattern =  "-P.*") 
gdat_chr8 <- gdat_chr8 %>%
  select(id, everything())
rownames(gdat_chr8)<-NULL

id_excld<-c("3D7","7G8","V1S","DD2","HB3","HC","LC","3d7","Dd2","NTC","PM2")
gdat_chr8_f<-gdat_chr8[!gdat_chr8$id%in%id_excld, ]

# Step 3: Prepare metadata
metadata <- gdat_dhps1[, 1:5]
metadata <- metadata[order(metadata$id, decreasing = FALSE), ]

# combine Meta data and Genotype data
gdat_metdat2<-left_join(gdat_chr8_f, metadata,by="id") %>% select(id, site, year, allele, everything()) %>% drop_na()


############################# Drop polyclonal Samples #################################


gdata_metdata <- gdat_metdat2[gdat_metdat2$id%in%DominatCloneLIst,]


############################# Drop polyclonal Samples #################################

############## PfDHFR data and collected in 2022 ###########
# gdat_dhfr1<-gdat_metdat %>%
gdat_dhps1<-gdata_metdata %>% 
  filter(
    allele %in% c("PfDHPS 581G")  #& 
    #year == "2022"
  )


########### select only genotype data
gdat_dhps<-gdat_dhps1[, 5:ncol(gdat_dhps1)]


# 550117, 581G Position
#25kb region= 525117 - 575117 (198:418)

rownames(gdat_dhps)<-gdat_dhps1$id
gdat_dhps_25kb<-gdat_dhps[,198:541] ## To only look at 50kb region (-/+25kb)

gdat_dhps_25kb<-as.matrix(gdat_dhps_25kb)
gdat_dhps_25kb_num <- apply(gdat_dhps_25kb, 2, as.numeric)
rownames(gdat_dhps_25kb_num) <- rownames(gdat_dhps_25kb)
colnames(gdat_dhps_25kb_num) <- colnames(gdat_dhps_25kb)


pp_dhps<-as.matrix(gdat_dhps_25kb)

# rownames_old <- rownames(pp_dhps)
# colnames_old <- colnames(pp_dhps)



########## metadata dhps data ################
met_dhps<-gdat_dhps1[,1:4]
met_dhps<-met_dhps[order(met_dhps$id, decreasing = F),]

met_dhps$region<-met_dhps$site
met_dhps <- met_dhps %>%
  mutate(region = recode(region, !!!regions)) 


px2=data.frame(Allele=met_dhps$allele,
               Year=met_dhps$year,
               Site=as.factor(met_dhps$site) ) 

colk11 <- list(
  
  Allele = c(
    "PfDHFR 164L" = "#FFE5B4",
    "PfDHPS 581G" = "#2D5765"#,
    #"PfDHPS A581" = "#2D5765"
  ),
  
  Year = c(
    "2016" = "#D3DAD0",
    "2017" = "#B2BEB5",
    "2020" = "#9CA59F",
    "2021" = "#868F89",
    "2022" = "#5A635E"
  ),
  
  Site = c(
    "Agago" = "#D79A6D",
    "Amolatar" = "#4D470F",
    "Arua" = "#9C931D",
    "Hoima" = "#A01F66",
    "Jinja" = "#FF7F00",
    "Kabale" = "#39355A",
    "Kaabong" = "#1F78B4",
    "Kanungu" = "#7570B3",
    "Kole" = "#FBF9C4",
    "Kasese" = "#D1CFEA",
    "Katakwi" = "#2F4F10",
    "Lamwo" = "#6F3B1A",
    "Mubende" = "#F7C2DA",
    "Tororo" = "#A7CA5D",
    "Kapchorwa" = "#FBF9C4",
    "Koboko" = "#9C931D"
  )
)




hk1<-rowAnnotation(df=px2, 
                   col = colk11,
                   show_annotation_name=F,
                   annotation_name_gp = gpar(fontsize=10),
                   annotation_name_rot=90,
                   annotation_legend_param = list(
                     legend_height = unit(12, "cm"),
                     title_gp=gpar(fontsize=9),
                     labels_gp = gpar(fontsize = 8)
                   ),
                   border = T)




# "#733c73", "#c0a0c0","#f1e5f1", "#32CD32"

#"#703c70""#926192""#b585b5""#d8aad8"
#plot(hk)
pal1<-colorRamp2(c(0,1,2), 
                 c("grey82","#d8aad8","#4d194d"))

pal2 <- colorRamp2(c(0, 1, 2), c("#4d194d", "#b980b9",  "#f1e5f1"))


# Filter out non variant positions
variant_dhsp_genotype_matrix <- gdat_dhps_25kb_num[, apply(gdat_dhps_25kb_num, 2, function(x) length(unique(x)) > 1)]

# Check dimensions before and after
dim(variant_dhsp_genotype_matrix)
dim(gdat_dhps_25kb_num)

variant_dhsp_genotype_matrix <- apply(variant_dhsp_genotype_matrix, 2, as.numeric)



dhps <- Heatmap(variant_dhsp_genotype_matrix,
                col=pal2,
                row_dend_side = "left",
                cluster_rows = T,
                cluster_columns = F,
                name = "dhpsht",
                right_annotation =  hk1,
                row_names_gp = gpar(fontsize=0),
                column_names_gp  = gpar(fontsize=0),
                column_names_rot = 45,
                width = 0.01,
                show_heatmap_legend = T,
                heatmap_legend_param = list(
                  title = "Genotype", at = c(0, 2),
                  labels = c("REF", "ALT"),
                  color_bar = "discrete",
                  title_gp=gpar(fontsize=9),
                  labels_gp = gpar(fontsize = 8)
                ),
                row_dend_width  = unit(7.5, "mm"),
                row_gap = unit(0.5, "mm"), 
                border = T)



#png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Haplotype_sharing_dhps20250814.png", 
#   width = 15, height = 10, units = "cm", res = 450 ) 

pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Suppl.Fig.3b.pdf", 
    width = 7, height = 4.5) 

draw(
  dhps,     
  ht_gap = unit(rep(1, 6), "mm"),
  annotation_legend_side = "right" ,
  heatmap_legend_side = "right"
)


col_index <- which(colnames( dhps@matrix ) == "X550036")

decorate_heatmap_body("dhpsht", {
  grid.lines(
    x = unit(c(col_index, col_index) / ncol(dhps@matrix), "npc") ,
    y = unit(c(0, 1), "npc") ,
    gp = gpar(col = "red3", lwd = 1)
  )
})

dev.off()   











#############################################################
##### Finding number of clusters
# ---- 1. Example Data ----
genotype_matrix_num


# ---- 2. Create Heatmap ----
ht <- Heatmap(genotype_matrix_num, cluster_rows = TRUE, cluster_columns = F)
ht_drawn <- draw(ht)


# ---- 3. Optimal k using silhouette ----
find_optimal_k <- function(dist_mat, k_range = 2:20) {
  sil_width <- sapply(k_range, function(k) {
    pam_fit <- pam(dist_mat, diss = TRUE, k = k)
    pam_fit$silinfo$avg.width
  })
  best_k <- k_range[which.max(sil_width)]
  return(best_k)
}

# For rows
dist_rows <- dist(genotype_matrix_num)
optimal_k_rows <- find_optimal_k(dist_rows)
cat("Optimal number of row clusters:", optimal_k_rows, "\n")


# ---- 4. Cut dendrogram using optimal k ----
row_dend <- row_dend(ht_drawn)
row_clusters <- cutree(as.hclust(row_dend), k = 9)


# ---- 5. Output cluster sizes ----
cat("\nRow cluster sizes:\n")
print(table(row_clusters))


# ---- 6. Combine cluster assignments with names ----
row_cluster_df <- data.frame(Gene = names(row_clusters), RowCluster = row_clusters)

# View the results
print(row_cluster_df)

# ---- 6. Cluster analysis
row_cluster_df$site <- str_remove_all(str_sub(row_cluster_df$Gene, 1, 3), "-")
row_cluster_df$RowCluster <- as.factor(row_cluster_df$RowCluster)


## Total number members in each cluster
row_cluster_df %>%
  group_by(site, 
           RowCluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  #filter(n>3) %>%
  arrange(RowCluster)

## At how many sites do we see members of each cluster group
# row_cluster_df %>%
#   group_by(site, RowCluster) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(RowCluster) %>%
#   summarise(num_sites = n_distinct(site))%>%
#   arrange(num_sites)

cluster_summary <- row_cluster_df %>%
  group_by(site, RowCluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(RowCluster) %>%
  summarise(
    num_sites = n_distinct(site),
    cluster_size = sum(n),
    sites = paste(sort(unique(site)), collapse = ", "),  # list sites
    .groups = "drop"
  ) %>%
  arrange(RowCluster)

cluster_summary

write.csv(cluster_summary, "/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/cluster_summary.csv", row.names = FALSE)




  
      

  # Extract the dendrogram and calculate clusters
 
 
 row_dend <- row_dend(dhfr)
 row_clusters <- cutree(row_dend, k = 4)  # choose k manually
 table(row_clusters) 
   
####################### PLOT HEATMAPS #######################
######## all site

generate_dhfr_heatmap( gdat_metdat_filtered)
##### 1. Plotting heatmaps

      

##### 2. Join HeatMaps
# 1. Heatman plotting
ht_list = result_North%v%result_NWest%v%result_West%v%result_S.West%v%result_Central%v%result_E.Central%v%result_East

draw(ht_list)




#######################################################################

# Genotyping data chr8
gdat_chr8a <- read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/AF_analysis_2024.10.11/20250604_Reanalysis/haplotypeReanalysis/DRHAP.chr8_GenotypeTable_ReCoded.tsv",
                       sep = "\t", 
                       check.names = FALSE)

gdat_chr8<-gdat_chr8a %>% 
  select(-c(REF,ALT,QUAL)) %>% t()
colnames(gdat_chr8) <- gdat_chr8[2, ]
gdat_chr8 <- data.frame(gdat_chr8[-c(1:3), ])


gdat_chr8$id <- str_remove_all(rownames(gdat_chr8), pattern =  "-P.*") 
gdat_chr8 <- gdat_chr8 %>%
  select(id, everything())
rownames(gdat_chr8)<-NULL

id_excld<-c("3D7","7G8","V1S","DD2","HB3","HC","LC","3d7","Dd2","NTC","PM2")
gdat_chr8_f<-gdat_chr8[!gdat_chr8$id%in%id_excld, ]

# Step 3: Prepare metadata
metadata <- gdat_dhps1[, 1:5]
metadata <- metadata[order(metadata$id, decreasing = FALSE), ]

# combine Meta data and Genotype data
gdat_metdat2<-left_join(gdat_chr8_f, metadata,by="id") %>% select(id, site, year, allele, everything()) %>% drop_na()


############################# Drop polyclonal Samples #################################


gdata_metdata <- gdat_metdat2[gdat_metdat2$id%in%DominatCloneLIst,]


############################# Drop polyclonal Samples #################################

############## PfDHFR data and collected in 2022 ###########
# gdat_dhfr1<-gdat_metdat %>%
gdat_dhps1<-gdata_metdata %>% 
  filter(
    allele %in% c("PfDHPS 581G")  #& 
    #year == "2022"
  )


########### select only genotype data
gdat_dhps<-gdat_dhps1[, 5:ncol(gdat_dhps1)]


# 550117, 581G Position
#25kb region= 525117 - 575117 (198:418)

rownames(gdat_dhps)<-gdat_dhps1$id
gdat_dhps_25kb<-gdat_dhps[,198:541] ## To only look at 50kb region (-/+25kb)

gdat_dhps_25kb<-as.matrix(gdat_dhps_25kb)
gdat_dhps_25kb_num <- apply(gdat_dhps_25kb, 2, as.numeric)
rownames(gdat_dhps_25kb_num) <- rownames(gdat_dhps_25kb)
colnames(gdat_dhps_25kb_num) <- colnames(gdat_dhps_25kb)


pp_dhps<-as.matrix(gdat_dhps_25kb)

# rownames_old <- rownames(pp_dhps)
# colnames_old <- colnames(pp_dhps)



########## metadata dhps data ################
met_dhps<-gdat_dhps1[,1:4]
met_dhps<-met_dhps[order(met_dhps$id, decreasing = F),]

met_dhps$region<-met_dhps$site
met_dhps <- met_dhps %>%
  mutate(region = recode(region, !!!regions)) 


px2=data.frame(Site=as.factor(met_dhps$site), 
               Genotype=met_dhps$allele,
               Year=met_dhps$year ) 
colk11 <- list(
  
  Genotype = c(
    "PfDHFR 164L" = "#FFE5B4",
    "PfDHPS 581G" = "#2D5765"
  ),
  
  Year = c(
    "2016" = "#D3DAD0",
    "2017" = "#B2BEB5",
    "2020" = "#9CA59F",
    "2021" = "#868F89",
    "2022" = "#5A635E"
  ),
  
  Site = c(
    "Agago" = "#D79A6D",
    "Amolatar" = "#4D470F",
    "Arua" = "#9C931D",
    "Hoima" = "#A01F66",
    "Jinja" = "#FF7F00",
    "Kabale" = "#39355A",
    "Kaabong" = "#1F78B4",
    "Kanungu" = "#7570B3",
    "Kole" = "#FBF9C4",
    "Kasese" = "#D1CFEA",
    "Katakwi" = "#2F4F10",
    "Lamwo" = "#6F3B1A",
    "Mubende" = "#F7C2DA",
    "Tororo" = "#A7CA5D"
  )
)

  


hk1<-rowAnnotation(df=px2, 
                   col = colk11,
                   show_annotation_name=F,
                   annotation_name_gp = gpar(fontsize=12),
                   annotation_name_rot=90,
                   annotation_legend_param = list(
                     legend_height = unit(12, "cm"),
                     title_gp=gpar(fontsize=12),
                     labels_gp = gpar(fontsize = 10)
                   ),
                   border = T)


# "#733c73", "#c0a0c0","#f1e5f1", "#32CD32"

#"#703c70""#926192""#b585b5""#d8aad8"
#plot(hk)
pal1<-colorRamp2(c(0,1,2), 
                 c("grey82","#d8aad8","#4d194d"))

pal2 <- colorRamp2(c(0, 1, 2), c("#4d194d", "#b980b9",  "#f1e5f1"))


# Filter out non variant positions
variant_dhsp_genotype_matrix <- gdat_dhps_25kb_num[, apply(gdat_dhps_25kb_num, 2, function(x) length(unique(x)) > 1)]

# Check dimensions before and after
dim(variant_dhsp_genotype_matrix)
dim(gdat_dhps_25kb_num)

variant_dhsp_genotype_matrix <- apply(variant_dhsp_genotype_matrix, 2, as.numeric)



dhps <- Heatmap(variant_dhsp_genotype_matrix,
        col=pal2,
        row_dend_side = "left",
        cluster_rows = T,
        cluster_columns = F,
        name = "dhpsht",
        right_annotation =  hk1,
        row_names_gp = gpar(fontsize=0),
        column_names_gp  = gpar(fontsize=0),
        column_names_rot = 45,
        width = 0.01,
        show_heatmap_legend = T,
        heatmap_legend_param = list(
          title = "Genotype", at = c(0, 2),
          labels = c("REF", "ALT"),
          color_bar = "discrete",
          title_gp=gpar(fontsize=8),
          labels_gp = gpar(fontsize = 6)
        ),
        row_dend_width  = unit(7.5, "mm"),
        row_gap = unit(0.5, "mm"), 
        border = T)



png("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Haplotype_sharing_dhps20250814.png", 
    width = 15, height = 10, units = "cm", res = 450 ) 



draw(
  dhps,     
  ht_gap = unit(rep(1, 6), "mm"),
  annotation_legend_side = "right" ,
  heatmap_legend_side = "right"
)


col_index <- which(colnames( dhps@matrix ) == "X550036")

decorate_heatmap_body("dhpsht", {
  grid.lines(
    x = unit(c(col_index, col_index) / ncol(dhps@matrix), "npc") ,
    y = unit(c(0, 1), "npc") ,
    gp = gpar(col = "red3", lwd = 1.5)
  )
})



dev.off()   





#################################
## HeatMap Concatenation
ht_list = ht_dhfr %v% ht_dhps
draw(ht_list)

