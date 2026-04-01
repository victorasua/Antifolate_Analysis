############################################################
# IBD–Genotype Association Analysis (DHFR / DHPS)
# Publication-ready version
############################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tibble)
library(scales)

############################################################
# 1. INPUT
############################################################

df <- read.csv(
  "data/IBD-Genotype-Data.csv",
  sep = "\t",
  stringsAsFactors = FALSE
)

############################################################
# 2. GENOTYPE ENCODING FUNCTION
############################################################

encode_pair <- function(df, genotype_table, new_col) {
  
  df <- df %>%
    left_join(genotype_table, by = c("sample1" = "id")) %>%
    rename(g1 = genotype)
  
  df <- df %>%
    left_join(genotype_table, by = c("sample2" = "id")) %>%
    rename(g2 = genotype)
  
  df <- df %>%
    mutate(
      !!new_col := case_when(
        g1 == 0 & g2 == 0 ~ "WT-WT",
        g1 == 1 & g2 == 1 ~ "MT-MT",
        (g1 == 1 & g2 == 0) | (g1 == 0 & g2 == 1) ~ "WT-MT",
        TRUE ~ NA_character_
      )
    )
  
  df
}

############################################################
# 3. LOAD GENOTYPE TABLES (assumed precomputed)
############################################################

# expected format: id | genotype (0/1)
# only_164_genotype
# only_581_genotype
# both_164_581_genotype

############################################################
# 4. BUILD DATASETS
############################################################

df_164 <- encode_pair(df, only_164_genotype, "pair_164") %>%
  select(sample1, sample2, fract_sites_IBD, pair_164) %>%
  drop_na()

df_581 <- encode_pair(df, only_581_genotype, "pair_581") %>%
  select(sample1, sample2, fract_sites_IBD, pair_581) %>%
  drop_na()

df_both <- encode_pair(df, both_164_581_genotype, "pair_both") %>%
  select(sample1, sample2, fract_sites_IBD, pair_both) %>%
  drop_na()

############################################################
# 5. IBD CATEGORISATION FUNCTION
############################################################

add_ibd_category <- function(df) {
  df %>%
    mutate(
      ibd_category = case_when(
        fract_sites_IBD < 0.1 ~ "Unrelated (<0.1)",
        fract_sites_IBD > 0.5 ~ "High (>0.5)",
        TRUE ~ "Intermediate (0.1–0.5)"
      )
    )
}

############################################################
# 6. PREPARE DATA
############################################################

df_164 <- add_ibd_category(df_164)

############################################################
# 7. SITE FILTERING
############################################################

key_sites <- c("Amolatar", "Arua", "Hoima", "Jinja",
               "Kabale", "Kanungu", "Lamwo")

df_164 <- df_164 %>%
  mutate(siteStatus = Sample2Sites) %>%
  filter(siteStatus %in% key_sites) %>%
  mutate(siteStatus = recode(siteStatus, "Kabale" = "Rukiga"))

############################################################
# 8. AGGREGATE COUNTS
############################################################

df_counts <- df_164 %>%
  group_by(siteStatus, ibd_category, pair_164) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = pair_164, values_from = count, values_fill = 0)

############################################################
# 9. CHI-SQUARED TEST FUNCTION
############################################################

run_chi <- function(df_site) {
  
  df_site <- df_site %>%
    complete(
      ibd_category = c("Unrelated (<0.1)", "Intermediate (0.1–0.5)", "High (>0.5)"),
      fill = list(`MT-MT` = 0, `WT-MT` = 0, `WT-WT` = 0)
    )
  
  mat <- df_site %>%
    select(ibd_category, `MT-MT`, `WT-MT`, `WT-WT`) %>%
    column_to_rownames("ibd_category") %>%
    as.matrix()
  
  chi <- chisq.test(mat)
  
  resid <- as.data.frame(chi$stdres) %>%
    rownames_to_column("ibd_category") %>%
    pivot_longer(-ibd_category,
                 names_to = "pair_type",
                 values_to = "std_residual") %>%
    mutate(siteStatus = df_site$siteStatus[1])
  
  list(
    chi_test = tibble(
      siteStatus = df_site$siteStatus[1],
      chisq = unname(chi$statistic),
      df = unname(chi$parameter),
      p_value = chi$p.value
    ),
    residuals = resid
  )
}

############################################################
# 10. APPLY CHI-SQUARED ACROSS SITES
############################################################

results <- df_counts %>%
  group_by(siteStatus) %>%
  group_split() %>%
  map(run_chi)

chi_results <- map_dfr(results, "chi_test")
resid_all <- map_dfr(results, "residuals")

############################################################
# 11. PLOT: STANDARDISED RESIDUALS
############################################################

resid_all <- resid_all %>%
  mutate(
    siteStatus = factor(siteStatus),
    ibd_category = factor(ibd_category)
  )

pdf("figures/IBD_residuals.pdf", width = 7, height = 3)

ggplot(resid_all,
       aes(x = pair_type,
           y = siteStatus,
           fill = std_residual)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(std_residual, 1)), size = 2.5) +
  facet_wrap(~ ibd_category) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(face = "bold")
  )

dev.off()

############################################################
# 12. PLOT: PROPORTIONAL STACKED BARPLOTS
############################################################

df_long <- df_counts %>%
  pivot_longer(cols = c(`MT-MT`, `WT-MT`, `WT-WT`),
               names_to = "pair_type",
               values_to = "count") %>%
  group_by(siteStatus, pair_type) %>%
  mutate(prop = count / sum(count))

pdf("figures/IBD_proportions.pdf", width = 7, height = 4)

ggplot(df_long,
       aes(x = pair_type,
           y = prop,
           fill = ibd_category)) +
  geom_col() +
  facet_wrap(~ siteStatus) +
  scale_y_continuous(labels = percent) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90)
  )

dev.off()