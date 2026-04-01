############################################################
# DHFR 164L vs I164 Haplotype Heatmap (Reproducible Version)
# Author: <your name>
############################################################

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(grid)
})

############################################################
# MAIN FUNCTION
############################################################

plot_dhfr_heatmap <- function(
    input_file,
    output_pdf,
    core_snp_pos = 748576,
    sample_size = 300,
    seed = 123
) {
  
  set.seed(seed)
  
  ##########################################################
  # 1. LOAD DATA
  ##########################################################
  
  df <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # harmonise site names
  df$site[df$site == "Kabale"] <- "Rukiga"
  
  meta_cols <- c("id", "site", "year", "allele", "region")
  geno_cols <- setdiff(names(df), meta_cols)
  
  ##########################################################
  # 2. FILTER VARIANTS (5–95% ALT FREQUENCY)
  ##########################################################
  
  geno_mat <- apply(df[, geno_cols], 2, function(x) as.numeric(trimws(x)))
  
  alt_freq <- colMeans(geno_mat == 2, na.rm = TRUE)
  
  keep_cols <- names(alt_freq)[alt_freq >= 0.05 & alt_freq <= 0.95]
  
  df <- cbind(df[, meta_cols], df[, keep_cols, drop = FALSE])
  
  ##########################################################
  # 3. FILTER BY BIOLOGICAL CRITERIA
  ##########################################################
  
  df <- df %>%
    filter(
      allele %in% c("PfDHFR 164L", "PfDHFR I164"),
      !region %in% c("Kapchorwa", "Koboko")
    )
  
  ##########################################################
  # 4. BALANCED SAMPLING
  ##########################################################
  
  d1 <- df %>% filter(allele == "PfDHFR 164L")
  d2 <- df %>% filter(allele == "PfDHFR I164")
  
  n_each <- sample_size / 2
  
  if (nrow(d1) < n_each | nrow(d2) < n_each) {
    stop("Not enough samples in one allele group.")
  }
  
  sampled <- bind_rows(
    d1[sample(nrow(d1), n_each), ],
    d2[sample(nrow(d2), n_each), ]
  )
  
  ##########################################################
  # 5. PREP GENOTYPE MATRIX
  ##########################################################
  
  meta <- sampled[, meta_cols]
  geno <- apply(sampled[, setdiff(names(sampled), meta_cols)], 2, as.numeric)
  
  rownames(geno) <- sampled$id
  
  # remove invariant SNPs
  geno <- geno[, apply(geno, 2, function(x) length(unique(x)) > 1), drop = FALSE]
  
  ##########################################################
  # 6. ALIGN METADATA
  ##########################################################
  
  meta <- meta[match(rownames(geno), meta$id), ]
  
  if (any(is.na(meta$id))) {
    stop("Metadata alignment failed.")
  }
  
  annotation_df <- data.frame(
    Site = factor(meta$site),
    Allele = factor(meta$allele),
    Year = factor(meta$year)
  )
  
  ##########################################################
  # 7. ANNOTATION COLORS
  ##########################################################
  
  col_list <- list(
    Site = c(
      Agago = "#D79A6D", Amolatar = "#4D470F", Arua = "#9C931D",
      Hoima = "#A01F66", Jinja = "#FF7F00", Rukiga = "#39355A",
      Kaabong = "#1F78B4", Kanungu = "#7570B3", Kasese = "#D1CFEA",
      Katakwi = "#2F4F10", Lamwo = "#6F3B1A", Mubende = "#F7C2DA",
      Tororo = "#A7CA5D"
    ),
    Allele = c(
      "PfDHFR I164" = "#75B7D1",
      "PfDHFR 164L" = "#FFE5B4"
    ),
    Year = c(
      "2016" = "#D3DAD0",
      "2017" = "#B2BEB5",
      "2020" = "#9CA59F",
      "2021" = "#868F89",
      "2022" = "#5A635E"
    )
  )
  
  ha <- rowAnnotation(
    df = annotation_df,
    col = col_list,
    show_annotation_name = FALSE,
    simple_anno_size = unit(2.5, "mm")
  )
  
  ##########################################################
  # 8. SNP POSITION HANDLING
  ##########################################################
  
  col_pos <- as.numeric(gsub("X", "", colnames(geno)))
  
  if (!core_snp_pos %in% col_pos) {
    warning("Core SNP not found in dataset.")
  }
  
  # select 5 before and 5 after
  before <- col_pos[col_pos < core_snp_pos]
  after  <- col_pos[col_pos > core_snp_pos]
  
  pick_n <- function(x, n = 5) {
    if (length(x) <= n) return(x)
    x[round(seq(1, length(x), length.out = n))]
  }
  
  selected_pos <- sort(unique(c(
    pick_n(before, 5),
    core_snp_pos,
    pick_n(after, 5)
  )))
  
  ##########################################################
  # 9. HEATMAP MATRIX PREP
  ##########################################################
  
  colnames(geno) <- col_pos
  
  display_names <- ifelse(col_pos %in% selected_pos, col_pos, "")
  
  geno_display <- geno
  colnames(geno_display) <- display_names
  
  ##########################################################
  # 10. COLOUR SCALE
  ##########################################################
  
  pal <- colorRamp2(c(0, 1, 2), c("#6d194d", "#f1e5f1", "#b990b9"))
  
  ##########################################################
  # 11. HEATMAP
  ##########################################################
  
  ht <- Heatmap(
    geno_display,
    name = "Genotype",
    col = pal,
    
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    
    row_names_gp = gpar(fontsize = 0),
    
    right_annotation = ha,
    
    column_names_gp = gpar(fontsize = 6.5),
    
    heatmap_legend_param = list(
      title = "Genotype",
      at = c(0, 2),
      labels = c("REF", "ALT")
    )
  )
  
  ##########################################################
  # 12. OUTPUT
  ##########################################################
  
  pdf(output_pdf, width = 6, height = 6)
  draw(ht)
  
  # core SNP marker
  col_index <- which(colnames(geno_display) == core_snp_pos)
  
  if (length(col_index) == 1) {
    decorate_heatmap_body("Genotype", {
      grid.lines(
        x = unit(c(col_index, col_index) / ncol(geno_display), "npc"),
        y = unit(c(0, 1), "npc"),
        gp = gpar(col = "red3", lwd = 1.2)
      )
    })
  }
  
  dev.off()
  
  return(ht)
}

############################################################
# (EDIT PATHS ONLY HERE)
############################################################

# ht <- plot_dhfr_heatmap(
#   input_file = "data/SNP_data_4_HaplotypeStructure.csv",
#   output_pdf = "output/dhfr_heatmap.pdf"
# )

############################################################
# SESSION INFO (FOR PUBLICATION)
############################################################

sessionInfo()

