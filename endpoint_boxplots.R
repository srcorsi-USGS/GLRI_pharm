library(toxEval)
library(dplyr)
library(tidyr)
library(ggplot2)

source(file = "plot_tox_endpoints_manuscript.R")

path_to_file <- 'processed_data/pharm_data.xlsx' 
tox_list <- create_toxEval(path_to_file)
ACC <- get_ACC(tox_list$chem_info$CAS)
ACC <- remove_flags(ACC = ACC,
                    flagsShort = c('Borderline','OnlyHighest','GainAC50','Biochemical'))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c('ATG','NVS','OT','TOX21','CEETOX','APR','CLD','TANGUAY','NHEERL_PADILLA','NCCT_SIMMONS','ACEA'),
                             remove_groups = c('Background Measurement','Undefined'))

chemical_summary <- get_chemical_summary(tox_list, ACC, filtered_ep)

threshold <- 0.0001
siteThreshold <- 5

endpoints_sites_hits <- filter(chemical_summary,EAR > 0) %>%
  group_by(endPoint,site,date) %>%
  summarize(EARsum = sum(EAR)) %>%
  group_by(site,endPoint) %>%
  summarize(EARmax = max(EARsum)) %>%
  filter(EARmax >= threshold) %>%
  group_by(endPoint) %>%
  summarize(numSites = n_distinct(site)) %>%
  arrange(desc(numSites)) %>%
  filter(numSites >= siteThreshold)

priority_endpoints <- endpoints_sites_hits$endPoint

chemicalSummaryPriority <- filter(chemical_summary, endPoint %in% priority_endpoints)

AOP_crosswalk <- read.csv("raw_data/AOP_crosswalk.csv", stringsAsFactors = FALSE)

AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

eps_with_ids <- unique(AOP$endPoint)

chemicalSummaryPriority$has_AOP <- "AOP Undefined"
chemicalSummaryPriority$has_AOP[chemicalSummaryPriority$endPoint %in% eps_with_ids] <- "AOP Associated"

endpointPlot <- plot_tox_endpoints_manuscript(chemicalSummaryPriority, 
                                              category = "Chemical", 
                                              font_size = 8,title = " ",
                                              pallette = c("steelblue", "white"))

endpointPlot
ggsave(endpointPlot, filename = "endpoint.png")
