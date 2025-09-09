
library(dplyr)
library(purrr)
library(broom)
library(dplyr)
library(purrr)
library(FSA)
library(trend)

MeanTrendData <- meanIBD %>% select(site, Year, trueGeno, mean_GenPropShared)

mk_test_mean <- MeanTrendData %>% 
  group_by(site, trueGeno) %>%
  filter(sum(mean_GenPropShared) > 0, n() > 2) %>%
  summarise(test = list(mk.test(mean_GenPropShared)), .groups = "drop") %>%
  mutate(tidy_test = map(test, broom::tidy)) %>%
  unnest(tidy_test) %>%
  select(site, trueGeno, statistic, p.value)

