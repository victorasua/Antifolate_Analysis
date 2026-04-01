############################################################
# Uganda Study Site Map (ADM2 + Water + Sampling Sites)
# Publication-ready spatial visualization
############################################################

library(tidyverse)
library(sf)
library(ggpubr)
library(ggspatial)

############################################################
# 1. INPUT DATA
############################################################

uganda_adm2 <- read_sf(
  "~/data/uga_admbnda_adm2_ubos_20200824.shp"
)

uga_water <- read_sf(
  "~/data/Ug_Waterbodies.shp"
)

facility <- read_csv(
  "~/data/map_data_2021.csv"
) %>%
  select(district, facility, N, lon, lat) %>%
  distinct() %>%
  drop_na()

############################################################
# 2. CLEAN DISTRICT NAMES
############################################################

facility <- facility %>%
  mutate(
    district = recode(district,
                      "Kabale" = "Rukiga",
                      "kanungu" = "Kanungu")
  )

############################################################
# 3. DEFINE STUDY DISTRICTS
############################################################

study_districts <- c(
  "Koboko", "Arua", "Lamwo", "Kaabong",
  "Agago", "Katakwi", "Amolatar", "Hoima",
  "Kapchorwa", "Tororo", "Jinja", "Mubende",
  "Kasese", "Kanungu", "Rukiga", "Kole"
)

mi_districts <- uganda_adm2 %>%
  filter(ADM2_EN %in% study_districts)

############################################################
# 4. ADD SAMPLE SIZE ANNOTATIONS
############################################################

site_counts <- facility %>%
  group_by(district) %>%
  summarise(N = first(N), .groups = "drop")

mi_districts <- mi_districts %>%
  left_join(site_counts, by = c("ADM2_EN" = "district")) %>%
  mutate(label = paste0(ADM2_EN, "\n", N))

############################################################
# 5. CREATE LABEL COORDINATES (ROBUST METHOD)
############################################################

labels_df <- mi_districts %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  mutate(
    geometry = st_point_on_surface(geometry)
  ) %>%
  cbind(st_coordinates(.)) %>%
  as_tibble() %>%
  mutate(
    label = mi_districts$label
  )

# optional manual nudge for Kanungu
labels_df <- labels_df %>%
  mutate(
    X = if_else(label %in% c("Kanungu"), X + 0.2, X)
  )

############################################################
# 6. CREATE NATIONAL BOUNDARY (DISSOLVED ADM2)
############################################################

uganda_boundary <- uganda_adm2 %>%
  st_make_valid() %>%
  st_union() %>%
  st_as_sf()

############################################################
# 7. CLIP WATER BODIES
############################################################

uga_water_crop <- uga_water %>%
  st_make_valid() %>%
  st_intersection(st_transform(uganda_boundary, st_crs(uga_water)))

############################################################
# 8. PLOT
############################################################

map_uganda <- ggplot() +
  
  # National boundary
  geom_sf(
    data = uganda_boundary,
    fill = "white",
    color = "gray60",
    linewidth = 0.4
  ) +
  
  # Water bodies
  geom_sf(
    data = uga_water_crop,
    fill = "#2c7fb8",
    color = NA,
    alpha = 0.7
  ) +
  
  # Study districts outline
  geom_sf(
    data = mi_districts,
    fill = NA,
    color = "#fdae61",
    linewidth = 0.4
  ) +
  
  # Labels
  geom_text(
    data = labels_df,
    aes(X, Y, label = label),
    size = 3,
    fontface = "bold"
  ) +
  
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

############################################################
# 9. EXPORT FIGURE
############################################################

ggsave(
  filename = "figures/Uganda_study_sites_map.pdf",
  plot = map_uganda,
  width = 12,
  height = 12,
  units = "cm",
  device = cairo_pdf
)

ggsave(
  filename = "figures/Uganda_study_sites_map.png",
  plot = map_uganda,
  width = 12,
  height = 12,
  units = "cm",
  dpi = 300
)