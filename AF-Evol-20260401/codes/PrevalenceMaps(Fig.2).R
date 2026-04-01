############################################################
# Plasmodium falciparum DHFR/DHPS Prevalence Mapping
# Publication-ready spatial workflow
############################################################

library(tidyverse)
library(sf)
library(ggpubr)
library(scales)
library(stringr)

############################################################
# 1. DATA INPUT
############################################################

data2 <- read_csv(
  "/data/data_filtered/data2.csv"
)

uganda_adm2 <- read_sf(
  "/shapefiles/uga_admbnda_adm2_ubos_20200824.shp"
)

############################################################
# 2. CLEAN AND FILTER GENETIC DATA
############################################################

af_prev <- data2 %>%
  select(id, site, gene, aa_change, year, genotype) %>%
  drop_na() %>%
  filter(gene %in% c("dhfr-ts", "dhps"))

# harmonise genotype coding (collapse heterozygotes if needed)
af_prev <- af_prev %>%
  mutate(genotype = if_else(genotype == 2, 1, genotype))

############################################################
# 3. SAMPLE ID CLEANING
############################################################

af_prev_clean <- af_prev %>%
  mutate(
    sampID = id %>%
      str_remove_all("-PRX-06-1") %>%
      str_remove_all("-PRX-07-1") %>%
      str_replace("-\\d{2}$", "")
  ) %>%
  distinct()

############################################################
# 4. PREVALENCE ESTIMATION
############################################################

prevalence_df <- af_prev_clean %>%
  group_by(site, year, gene, aa_change) %>%
  summarise(
    total = n(),
    mutant = sum(genotype == 1),
    prevalence = mutant / total,
    .groups = "drop"
  )

############################################################
# 5. SITE STANDARDISATION
############################################################

site_map <- c(
  "TO" = "Tororo",
  "AR" = "Arua",
  "MU" = "Mubende",
  "JI" = "Jinja",
  "KB" = "Kabale",
  "AG" = "Agago",
  "LA" = "Lamwo",
  "KO" = "Kole",
  "AM" = "Amolatar",
  "KN" = "Kanungu",
  "HO" = "Hoima",
  "KBG" = "Kaabong",
  "KBK" = "Koboko",
  "KS" = "Kasese",
  "KTK" = "Katakwi",
  "KAP" = "Kapchorwa"
)

prevalence_df <- prevalence_df %>%
  mutate(site = recode(site, !!!site_map)) %>%
  mutate(site = if_else(site == "Kabale", "Rukiga", site))

############################################################
# 6. MERGE WITH SPATIAL DATA
############################################################

map_df <- uganda_adm2 %>%
  left_join(prevalence_df, by = c("ADM2_EN" = "site")) %>%
  select(ADM2_EN, geometry, year, aa_change, prevalence)

############################################################
# 7. DEFINE ANALYTICAL SUBSETS (NO REPEATED CODE BLOCKS)
############################################################

years <- c("2016-17", "2018-19", "2020", "2021", "2022")

sites_keep <- c(
  "Arua", "Lamwo", "Agago", "Amolatar", "Kole", "Jinja",
  "Tororo", "Kanungu", "Rukiga", "Mubende",
  "Kasese", "Hoima", "Kaabong"
)

map_filtered <- map_df %>%
  filter(year %in% years | is.na(year)) %>%
  mutate(
    prevalence = replace_na(prevalence, 0)
  )

############################################################
# 8. HANDLE MISSING YEARS (CROSS PRODUCT APPROACH)
############################################################

map_complete <- map_filtered %>%
  select(-year) %>%
  crossing(year = years) %>%
  bind_rows(map_filtered) %>%
  distinct()

############################################################
# 9. REBUILD SPATIAL OBJECT
############################################################

geom_ref <- map_complete %>%
  group_by(ADM2_EN) %>%
  slice(1) %>%
  ungroup() %>%
  select(ADM2_EN, geometry)

map_sf <- map_complete %>%
  left_join(geom_ref, by = "ADM2_EN") %>%
  st_as_sf()

############################################################
# 10. BIN PREVALENCE
############################################################

map_sf <- map_sf %>%
  filter(!is.na(prevalence)) %>%
  mutate(
    prev_bin = cut(
      prevalence,
      breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, 1),
      labels = c("Zero", ">0–0.2", ">0.2–0.4",
                 ">0.4–0.6", ">0.6–0.8", ">0.8–1.0")
    ),
    aa_change = recode(
      aa_change,
      "Ile164Leu" = "I164L",
      "Ala581Gly" = "A581G"
    )
  )

############################################################
# 11. COLOUR SCALE
############################################################

pal <- c(
  "Zero" = "#CACACA",
  ">0–0.2" = "#E9DFF0",
  ">0.2–0.4" = "#C3A9D6",
  ">0.4–0.6" = "#7B3FA0",
  ">0.6–0.8" = "#4B0082",
  ">0.8–1.0" = "#120019"
)

############################################################
# 12. MAP PLOT
############################################################

map_plot <- ggplot(map_sf) +
  geom_sf(aes(fill = prev_bin), color = NA) +
  facet_grid(aa_change ~ year) +
  scale_fill_manual(values = pal, na.translate = FALSE) +
  coord_sf(datum = NA) +
  theme_void() +
  theme_pubr(border = TRUE) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines")
  ) +
  labs(fill = "Mutant prevalence")

############################################################
# 13. SAVE OUTPUTS
############################################################

ggsave(
  filename = "figures/DHFR_DHPS_prevalence_map.pdf",
  plot = map_plot,
  width = 21/2.54,
  height = 12/2.54,
  device = cairo_pdf
)

