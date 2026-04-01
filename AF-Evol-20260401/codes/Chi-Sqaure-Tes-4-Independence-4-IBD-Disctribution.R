############################################################
# Chi-squared analysis of pair-type frequencies by IBD category
# Publication-ready version
############################################################

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

############################################################
# FUNCTION: run_chi_by_site
############################################################

run_chi_by_site <- function(df_site) {
  
  #----------------------------------------------------------
  # 1. Ensure complete IBD category structure
  #----------------------------------------------------------
  
  df_site <- df_site %>%
    complete(
      ibd_category = c("Unrelated (<0.1)", "High (>0.5)"),
      fill = list(
        `MT-MT` = 0,
        `WT-MT` = 0,
        `WT-WT` = 0
      )
    )
  
  #----------------------------------------------------------
  # 2. Construct contingency table
  #----------------------------------------------------------
  
  mat <- df_site %>%
    select(ibd_category, `MT-MT`, `WT-MT`, `WT-WT`) %>%
    tibble::column_to_rownames("ibd_category") %>%
    as.matrix()
  
  # Safety check: ensure valid contingency structure
  if (any(mat < 0)) {
    stop("Negative values found in contingency table.")
  }
  
  if (nrow(mat) < 2) {
    stop("Not enough IBD categories for Chi-squared test.")
  }
  
  #----------------------------------------------------------
  # 3. Chi-squared test
  #----------------------------------------------------------
  
  chi_out <- suppressWarnings(chisq.test(mat))
  
  #----------------------------------------------------------
  # 4. Standardized residuals
  #----------------------------------------------------------
  
  resid_long <- as.data.frame(chi_out$stdres) %>%
    rownames_to_column("ibd_category") %>%
    pivot_longer(
      cols = c("MT-MT", "WT-MT", "WT-WT"),
      names_to = "pair_type",
      values_to = "std_residual"
    ) %>%
    mutate(siteStatus = df_site$siteStatus[1]) %>%
    select(siteStatus, ibd_category, pair_type, std_residual)
  
  #----------------------------------------------------------
  # 5. Return structured output
  #----------------------------------------------------------
  
  list(
    chi_test = tibble(
      siteStatus = unique(df_site$siteStatus),
      chisq_statistic = unname(chi_out$statistic),
      df = unname(chi_out$parameter),
      p_value = chi_out$p.value
    ),
    residuals = resid_long
  )
}

############################################################
# APPLY FUNCTION ACROSS SITES
############################################################

results_list <- df %>%
  group_split(siteStatus) %>%
  map(run_chi_by_site)

############################################################
# COMBINE RESULTS
############################################################

chi_results <- map_dfr(results_list, "chi_test")

resid_all <- map_dfr(results_list, "residuals")

############################################################
# OPTIONAL: INSPECT OUTPUT
############################################################

chi_results
resid_all