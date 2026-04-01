############################################################
# IBD NETWORK ANALYSIS (DHFR / DHPS) - PUBLICATION VERSION
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(stringr)
})

rm(list = ls())

############################################################
# 1. INPUT PARAMETERS 
############################################################
#"/Users/victorasua/Library/CloudStorage/.../PCA/gen_dat.csv"
#"/Users/victorasua/Library/CloudStorage/.../my_p_clusters_5.rds"

input_file <- "data/gen_dat.csv"
cluster_rds <- "data/my_p_clusters_5.rds"

output_dir <- "figures"

dir.create(output_dir, showWarnings = FALSE)

############################################################
# 2. LOAD DATA
############################################################

gen_dat <- read.csv(input_file, stringsAsFactors = FALSE)

############################################################
# 3. DEFINE GENOTYPE TABLE
############################################################

dr_gen_dat <- gen_dat |>
  select(id, aa_change, genotype) |>
  distinct() |>
  pivot_wider(names_from = aa_change, values_from = genotype) |>
  select(id, Ile164Leu, Ala581Gly)

############################################################
# 4. RECODE GENOTYPES
############################################################

recode_geno <- function(x, gene) {
  case_when(
    x == 0 ~ "Reference",
    x %in% c(1, 2) ~ gene,
    TRUE ~ NA_character_
  )
}

dr_gen_dat <- dr_gen_dat |>
  mutate(
    Ile164Leu  = recode_geno(Ile164Leu, "PfDHFR I164L"),
    Ala581Gly  = recode_geno(Ala581Gly, "PfDHPS A581G")
  )

############################################################
# 5. DEFINE COMBINED GENOTYPE
############################################################

dr_genotype2 <- dr_gen_dat |>
  mutate(
    trueGeno = case_when(
      Ile164Leu == "Reference" & Ala581Gly == "Reference" ~ "Reference",
      Ile164Leu == "PfDHFR I164L" & (Ala581Gly == "Reference" | is.na(Ala581Gly)) ~ "PfDHFR I164L",
      Ala581Gly == "PfDHPS A581G" & (Ile164Leu == "Reference" | is.na(Ile164Leu)) ~ "PfDHPS A581G",
      TRUE ~ NA_character_
    )
  ) |>
  select(iid = id, trueGeno)

############################################################
# 6. LOAD NETWORK
############################################################

my_p_clusters <- readRDS(cluster_rds)

iplot <- my_p_clusters$i.network

############################################################
# 7. CLEAN NODE IDS
############################################################

V(iplot)$name <- str_remove(V(iplot)$name, "-P.*")

node_df <- tibble(iid = V(iplot)$name) |>
  left_join(dr_genotype2, by = "iid")

############################################################
# 8. ADD METADATA (SITE + YEAR)
############################################################

node_df <- node_df |>
  mutate(
    site = str_sub(iid, 1, 2),
    yr_raw = str_sub(iid, 4, 6),
    year = recode(yr_raw,
                  "00" = "2016",
                  "02" = "2017",
                  "05" = "2020",
                  "06" = "2021",
                  "07" = "2022",
                  .default = NA_character_)
  ) |>
  mutate(
    site = recode(site,
                  "LAM" = "LA",
                  "TR"  = "TO",
                  "KS"  = "KS",
                  "KB"  = "RK",
                  "kS"  = "KS")
  )

############################################################
# 9. FILTER KEY SITES
############################################################

key_sites <- c("AM", "AR", "HO", "RK", "KN", "JI", "LA")

node_df <- node_df |>
  mutate(site = ifelse(site %in% key_sites, site, NA))

############################################################
# 10. ALIGN TO GRAPH
############################################################

V(iplot)$genotype <- node_df$trueGeno[match(V(iplot)$name, node_df$iid)]
V(iplot)$site     <- node_df$site[match(V(iplot)$name, node_df$iid)]
V(iplot)$year     <- node_df$year[match(V(iplot)$name, node_df$iid)]

############################################################
# 11. DROP MISSING
############################################################

iplot_clean <- delete_vertices(iplot, which(is.na(V(iplot)$genotype)))

############################################################
# 12. COLOR MAPS
############################################################

geno_cols <- c(
  "PfDHPS A581G" = "gold4",
  "PfDHFR I164L" = "#A52A2A",
  "Reference"    = "#67a9cf"
)

year_cols <- c(
  "2016" = "#1B9E77",
  "2017" = "#D95F02",
  "2020" = "#7570B3",
  "2021" = "#E7298A",
  "2022" = "#66A61E"
)

site_cols <- c(
  "JI" = "goldenrod", "RK" = "gold4", "AM" = "#2bce48",
  "AR" = "#0067a5", "KS" = "#a1caf1", "HO" = "#A52A2A",
  "KN" = "magenta3", "LA" = "#FF0000"
)

############################################################
# 13. APPLY COLORS
############################################################

set_vertex_color <- function(graph, column, palette) {
  V(graph)$color <- palette[V(graph)[[column]]]
  graph
}

iplot_geno <- set_vertex_color(iplot_clean, "genotype", geno_cols)
iplot_year <- set_vertex_color(iplot_clean, "year", year_cols)
iplot_site <- set_vertex_color(iplot_clean, "site", site_cols)

############################################################
# 14. COMMON PLOTTING FUNCTION
############################################################

plot_network <- function(graph, layout, legend_list, legend_title, file) {
  
  pdf(file, width = 9/2.54, height = 10/2.54, bg = "transparent")
  
  set.seed(44)
  
  layout[,1] <- layout[,1] - 0.2
  
  par(mar = c(3, 0, 0, 0))
  
  plot(
    graph,
    layout = layout,
    vertex.size = 4,
    vertex.color = V(graph)$color,
    vertex.frame.color = "black",
    vertex.frame.width = 0.5,
    edge.color = "black",
    edge.width = 0.5,
    vertex.label = NA
  )
  
  legend(
    x = -1.1, y = -1.05,
    legend = legend_list$labels,
    col = legend_list$cols,
    pch = 19,
    cex = 0.7,
    bty = "o",
    title = legend_title,
    horiz = TRUE
  )
  
  dev.off()
}

############################################################
# 15. LAYOUT
############################################################

layout_fr <- layout_with_fr(
  iplot_geno,
  niter = 100,
  area = vcount(iplot_geno)^3,
  repulserad = 1000
)

############################################################
# 16. EXPORT FIGURES
############################################################

plot_network(
  iplot_geno,
  layout_fr,
  list(
    labels = c("PfDHPS 581G", "PfDHFR 164L", "Reference"),
    cols   = c("gold4", "#A52A2A", "#67a9cf")
  ),
  "Genotype",
  file.path(output_dir, "ibd_genotype_50.pdf")
)

plot_network(
  iplot_year,
  layout_fr,
  list(
    labels = names(year_cols),
    cols   = year_cols
  ),
  "Year",
  file.path(output_dir, "ibd_year_50.pdf")
)

plot_network(
  iplot_site,
  layout_fr,
  list(
    labels = names(site_cols),
    cols   = site_cols
  ),
  "Site",
  file.path(output_dir, "ibd_site_50.pdf")
)

############################################################
# END
############################################################