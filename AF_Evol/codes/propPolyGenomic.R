
library(dplyr)
library(purrr)
library(broom)
library(dplyr)
library(purrr)
library(FSA)
library(trend)


#Statistical evaluation for proprtion of polygenomics infections
# Data source for proportion poly-genomic
Prop_Poly_Genomic_Overall<-readRDS("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Prop_Poly_Genomic_Overall.rds")

# 1. Chi-square test for trend results
trend_results <- Prop_Poly_Genomic_Overall %>%
  group_by(site) %>%
  filter(sum(count_coi_gt_2) > 0, n() > 1) %>%
  summarise(test = list(prop.trend.test(count_coi_gt_2, total)), .groups = "drop") %>%
  mutate(tidy_test = map(test, broom::tidy)) %>%
  unnest(tidy_test) %>%
  select(site, statistic, p.value)

# 2. Compute slope for direction
direction_results <- Prop_Poly_Genomic_Overall %>%
  group_by(site) %>%
  filter(sum(count_coi_gt_2) > 0, n() > 1) %>%
  summarise(
    slope = coef(lm(prop_coi_gt_2 ~ as.numeric(year)))[2],
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(slope > 0, "Increasing", "Decreasing"))

# 3. Merge
final_results <- left_join(trend_results, direction_results, by = "site") %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "<0.001",
      p.value < 0.01 ~ "<0.01",
      p.value < 0.05 ~ "<0.05",
      TRUE ~ "NS"
    )
  )

final_results
