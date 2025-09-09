library(tidyverse)
#library(dplyr)
#library(ggplot2)
library(RColorBrewer)
#library(readr)
library(readxl)
#library(stringr)
#library(tidyr)
##library(vroom)
library(lessR)
#library(gplots)
library(geneHapR)
library(haplotypes)
library(vcfR)
library(purrr)

#data2 <- read.csv("/Users/victorasua/Documents/PROJECTS/PhD_Files/thesis_data/data_filtered/data2.csv")
data2 <- read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Documents/PhD_Files/thesis_data/data_filtered/data2.csv")

View(data.frame(unique(data2$id)) |> arrange(data2$id))
length(unique(data2$id)|>order())

## Summary of genes with available data
colnames(data2)

af_prev <-  data2 |> 
  dplyr::select(id, site, gene, aa_change, year, genotype) |> 
  drop_na() |>
  #select(gene, site, aa_change, genotype) |>
  filter(gene=="dhfr-ts"|gene=="dhps")

poly0 <- c( "Ala16Val",  "Asn51Ile" , "Cys59Arg" , "Ile164Leu", "Ser108Asn", "Ser108Thr", "Ala437Gly", "Ala581Gly", "Ala613Ser", "Ala613Thr", "Asn521fs", 
            "Ile431Val", "Lys540Asn", "Lys540Glu", "Ser436Ala", "Ser436Cys", "Ser436His", "Ser436Phe", "Lys27Glu",  "Gly406Ala", "Ile206Leu", "Phe217Cys",
            "Leu22Phe",  "Val20Ile",  "Glu108Lys", "Gly669Gly", "Ile104Leu", "Met648Leu", "Phe649Val", "Leu53Arg" )

poly <- c("Val20Ile", "Lys27Glu", "Asn51Ile", "Leu53Arg", "Cys59Arg", "Ser108Asn", "Ile164Leu", "Leu22Phe", 
          "Ile104Leu", "Glu108Lys", "Ile206Leu", "Phe217Cys", "Gly406Ala", "Ser436Ala", "Ser436Cys", 
          "Ser436His", "Ala437Gly", "Lys540Glu", "Ala581Gly", "Met648Leu", "Phe649Val", "Gly669Gly")

poly1 <- c("Asn51Ile", "Cys59Arg", "Ser108Asn", "Ile164Leu", "Ser436Ala", "Ser436Cys", 
           "Ser436His", "Ala437Gly", "Lys540Glu", "Ala581Gly")

polym2 <- c( "Ala675Val", "Arg561His", "Cys469Phe", "Cys469Tyr")


#af_prev$aa_change <- factor(af_prev$aa_change, levels = poly)

af_prev2 <- af_prev
af_prev2$genotype[af_prev2$genotype==2] <- 1


af_prev3 <- af_prev2 %>%
  mutate(samp_id = str_remove_all(str_remove_all(id, "-PRX-06-1"), "-PRX-07-1")) %>%
  mutate(samp_id_split = str_split(samp_id, "-")) %>%
  unnest_wider(samp_id_split, names_sep = "_") |> 
  mutate(id_no=str_sub(samp_id_split_3, -2, -1),
         sampID=paste0(samp_id_split_1, "-",
                       samp_id_split_2, "-",
                       id_no)) |> 
  dplyr::select(sampID, site, year, gene, aa_change, genotype) |> unique()

# Calculate prevalence and confidence intervals for each aa_change
df_grp <- af_prev3|>
  group_by(site, 
           year,
           gene, 
           aa_change, 
           genotype) |> 
  summarise(n=n())


prevalence_df <- df_grp %>% #df %>%
  group_by(site, 
           year,
           gene,
           aa_change) %>%
  summarise(
    total_cases = sum(n),
    positive_cases = sum(ifelse(genotype == 1, n, 0))
  ) %>%
  rowwise() %>%
  mutate(
    prevalence = positive_cases / total_cases #,
   # ci = list(prop.test(positive_cases, total_cases)$conf.int)
  ) 





#######################################################################################

# Updated shape file
uganda_adm2 <- read_sf("/Users/victorasua/Library/CloudStorage/Box-Box/IMMRSE-U Study/Melissa Drug Resistance Code/shape_files/uga_admbnda_ubos_20200824_shp/uga_admbnda_adm2_ubos_20200824.shp")


# IMMRS data summaries
immrsPrev <- read.csv("/Users/victorasua/Library/CloudStorage/Box-Box/Drug Resistance Paper/databases/1_simple_tabulations_DIS.csv")
immrsPrev <- immrsPrev %>% dplyr::select(District, collection, mutation_name, prev_mix_mut)

# filter key k13 mutations
dhfr_dhps_loci <- c( "164", "581" )


immrsPrev2 <- immrsPrev %>%
  filter(str_detect(mutation_name, str_c(dhfr_dhps_loci, collapse = "|")))
colnames(immrsPrev2)
names(immrsPrev2)[]<-c("site","year","aa_change","prevalence")

## renaming values
yr <- c("1"="2023_a",
        "2"="2023_b",
        "3"="2024_a",
        "4"="2024_b")

sites <- c("TO"="Tororo","AR"="Arua","MU"="Mubende","JI"="Jinja","KB"="Kabale","AG"="Agago","LA"="Lamwo","KO"="Kole","AM"="Amolatar",
           "KN"="Kanungu","HO"="Hoima","KBG"="Kaabong","KBK"="Koboko","KS"="Kasese","KTK"="Katakwi","KAP"="Kapchorwa")

immrsPrev3 <- immrsPrev2
immrsPrev3 <- immrsPrev3%>%
  mutate(year=recode(year, !!!yr))

# PRX3 data summaries
#data2 <- read.csv("/Users/victorasua/Documents/PROJECTS/PhD_Files/thesis_data/data_filtered/data2.csv")
prxData <- read.csv("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Documents/PhD_Files/thesis_data/data_filtered/data2.csv")


## Summary of genes with available data
dhfr_dhps <-  prxData |> select(id, site, gene, aa_change, year, genotype) |> 
  drop_na() |>
  #select(gene, site, aa_change, genotype) |>
  filter(gene%in%c("dhfr-ts", 
                   "dhps"))

dhfr_dhps$genotype[dhfr_dhps$genotype == 1] <- 2

## Calculate prevalences
dhfr_dhps_prev <- dhfr_dhps|>
  group_by(site,
           year,
           gene,
           aa_change,
           genotype) |>
  summarise(n=n()) |>
  mutate(prev=n/sum(n)*100) |>
  ungroup()   |>
  filter(genotype == 2 #|
         # genotype == 1 
  ) |>
  select(site, year, gene, aa_change, genotype, prev) |> 
  pivot_wider(names_from = year, values_from = prev) |>
  filter(str_detect(aa_change, str_c(c( #"72", "73", "74", "75", 
    "164", "581"), collapse = "|"))) |> unnest()

dhfr_dhps_prev1 <- dhfr_dhps_prev |> 
  mutate(site = recode(site, !!!sites)) |>
  pivot_longer(cols = c(5:9), names_to = "year", values_to = "prevalence") |>
  mutate(prevalence = replace_na(prevalence, 0)) |> select(site, year, aa_change, prevalence)


# ## COmbine immrs-u and prx data sets
# 
# tpt_map_data <- rbind(immrsPrev3, crt_mdr_prev1) |> unnest()
# 
# tpt_map_data$aa_change[tpt_map_data$aa_change == "crt_K76T"] <- "Lys76Thr"
# tpt_map_data$aa_change[tpt_map_data$aa_change == "mdr1_N86Y"] <- "Asn86Tyr"
# tpt_map_data$aa_change[tpt_map_data$aa_change == "mdr1_D1246Y"] <- "Asp1246Tyr"
# 
# tpt_map_data <- tpt_map_data |>
#   filter(str_detect(aa_change, str_c(c( #"72", "73", "74", "75", 
#     "Lys76Thr", "Asn86Tyr", "Asp1246Tyr"), collapse = "|"))) |> unnest()
# 
# tpt_map_data2 <- tpt_map_data |> pivot_wider(names_from = year, values_from = prevalence) |>
#   mutate("2023" = `2023_a` + `2023_b`/2,
#          "2024" = `2024_a` + `2024_b`/2,
#          "2020-21" = `2020` + `2021`/2) |>
#   select(site, aa_change, `2016-17`, `2018-19`, "2020-21", "2022", "2023", "2024") %>%
#   pivot_longer(cols = c(3:8), names_to = "year", values_to = "prevalence", )
# 

# 1. Split NA and non-NA year rows
SNP_data_map_AF <- uganda_adm2 %>%
  dplyr::left_join(dhfr_dhps_prev1, by = c("ADM2_EN" = "site")) %>% 
  dplyr::select(ADM2_EN, geometry, year, aa_change, prevalence)




######################@@@@@@@@@@@@@@@@ MAPPING Prevalence --2 of Resistance genes @@@@@@@@@@@@@@@@@##################

# 1. Split NA and non-NA year rows
ug_na <- SNP_data_map_AF %>%
  filter(is.na(year))

ug_non_na <- SNP_data_map_AF %>%
  filter(!is.na(year))

# 2. Replicate NA rows for each year
ug_na_filled <- ug_na %>%
  select(-year) %>%       # Remove the existing year column
  crossing(year = c("2016-17", "2018-19", "2020", "2021", "2022")) %>%
  select(ADM2_EN, geometry, year, aa_change, prevalence)
#dplyr::mutate(Prev = NA)       # Optionally reset 'mean' to NA for new rows


# 3. Combine back together
SNP_data_map_AF_filled <- bind_rows(ug_non_na, ug_na_filled) #%>%
# pivot_wider(names_from = "year", values_from = "prevalence")


# fill all values for all sites
geom_ref <- SNP_data_map_AF_filled %>%
  group_by(ADM2_EN) %>%           # Group by district
  slice_head(n = 1) %>%                 # Take first row per group (fast!)
  ungroup() %>%                         
  select(ADM2_EN, geometry)       # Keep only needed columns

SNP_data_map_AF_plot <- SNP_data_map_AF_filled %>%
  complete(ADM2_EN, aa_change, year, fill = list(prevalence = NA)) %>%
  left_join(geom_ref, by = "ADM2_EN") %>% 
  select(-c("geometry.x") ) %>%
  select(ADM2_EN, "geometry.y", everything()) %>%
  dplyr::rename(geometry=geometry.y) %>%  
  st_as_sf()


SNP_data_map_AF_plot_filtered <- SNP_data_map_AF_plot %>%
  group_by(aa_change, year) %>%
  filter(any(!is.na(prevalence))) %>%
  ungroup()


SNP_data_map_tpt_plot_filtered_binned <- SNP_data_map_AF_plot_filtered %>% 
  mutate(prev_bin = case_when(prevalence == 0  ~ "Not detected",
                              prevalence > 0 & prevalence < 20 ~ "01–20",
                              prevalence > 20 & prevalence < 40 ~ "20–40",
                              prevalence > 40 & prevalence < 60 ~ "40–60",
                              prevalence > 60 & prevalence < 80 ~ "60–80",
                              prevalence > 80 & prevalence < 100 ~ "80–100",
                              TRUE ~ as.character(prevalence)),
         prev_bin2 = factor(prev_bin, levels = (c("Not detected", "01–20", "20–40", "40–60", "60–80", "80–100" ))))

# color scales
color3b <- c("Not detected"="#B0B0B0",
             "01–20" = "#FCFDBFFF",
             "20–40" = "#FB8861FF",
             "40–60" = "#B63679FF", 
             "60–80" = "#51127CFF",
             "80–100" = "#000004FF")



SNP_Map_AF_plot_data <- SNP_data_map_tpt_plot_filtered_binned %>%
  filter(!st_is_empty(geometry),
         aa_change %in% c("Ala581Gly", "Ile164Leu" )) %>%
  mutate(aa_change = factor(aa_change,
                            levels = c("Ala581Gly", "Ile164Leu" )))



#@#############################################
# png("/Users/victorasua/Library/CloudStorage/Box-Box/Drug Resistance Paper/figures/PropMutantMap4_tpt.png", 
#     width = 22, height = 14.4, units = "cm", res = 300 ) 

pdf("/Users/victorasua/Library/CloudStorage/OneDrive-Personal/Genomic/mipAnalysis/Fig.3.pdf",  
    width = 7, height = 3.86)

ggplot(data = SNP_Map_AF_plot_data, aes(fill = prev_bin2)) +
  geom_sf() +
  coord_sf(datum = NA) +
  scale_fill_manual(
    values = color3b,
    na.translate = F,
    name = "% mutant",
    guide = guide_legend(
      direction = "horizontal",
      label.position = "bottom",
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  
  theme_void() +
  facet_grid(aa_change ~ year , switch = "y"
  ) +
  theme_linedraw() +
  theme_pubr(border = TRUE, margin = F) +
  theme(axis.ticks = element_blank(),
        legend.key.width = unit(0.75, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 8, hjust = 0, vjust = 1, colour = "black" ),
        legend.title = element_text(size = 9, vjust = 1, colour = "black" ),
        legend.position = "bottom",
        legend.spacing.x = unit(-1, "pt"),
        legend.box = "horizontal",
        legend.box.margin = margin(t = -12, b = 0, unit = "pt"),
        plot.background = element_rect(),
        panel.spacing.x = unit(0.2, "line"),
        panel.spacing.y = unit(0.2, "line"),
        strip.text = element_text(size = 9, colour = "black", face = "bold"),
        strip.background = element_rect(size = 0.2),
        panel.border = element_rect(color = "black", size = 0.2),  # Adjust the border thickness (size)
        panel.grid.major = element_blank(),  # Optionally remove grid lines if needed
        panel.grid.minor = element_blank()
        #plot.margin = unit(c(0., 0, 0, 0), "cm"
  )

dev.off()


