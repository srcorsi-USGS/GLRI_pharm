library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(toxEval)

df <- readxl::read_xlsx("raw_data/Initial Prospective Study Sites for WY18 Pharm and PFAS SamplingSRC.xlsx") %>%
  select(site=STAID, ww_percent = `Percent of streamflow attributable to WW`) %>%
  arrange(desc(ww_percent))

################################################################
# Concentration:
path_to_file <- 'processed_data/pharm_data_concentrations.xlsx' 
tox_list_conc <- create_toxEval(path_to_file)
chemical_summary_conc <- get_chemical_summary(tox_list_conc)

chem_sum_1 <- chemical_summary_conc %>%
  left_join(df, by="site")

chem_sum_2 <- chem_sum_1 %>%
  filter(EAR > 0) %>%
  group_by(shortName, ww_percent, date) %>%
  summarize(sum_conc = sum(EAR, na.rm = TRUE),
            n_chems = n()) %>%
  group_by(shortName, ww_percent) %>%
  summarize(sum_conc_max = max(sum_conc),
            max_chem = max(n_chems)) %>%
  filter(!is.na(ww_percent))

conc_plot <- ggplot() +
  geom_point(data = chem_sum_2,
             aes(x = ww_percent,
                 y = sum_conc_max)) +
  geom_text_repel(data = filter(chem_sum_2, 
                                sum_conc_max > 5 |
                                  ww_percent > 0.1),
                  aes(x = ww_percent,
                      y = sum_conc_max,
                      label = shortName)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Percent of streamflow attributable to WW") +
  ylab("Sum of concentration (of max sample)") +
  ggtitle("Wastewater % vs Sum of concentrations (max sample)")
conc_plot

dir.create("plots",showWarnings = FALSE)
ggsave(conc_plot, filename = "plots/conc_log_labels.png", width = 8, height = 6)

chem_sum_3 <- chem_sum_1 %>%
  filter(EAR > 0) %>%
  group_by(shortName, ww_percent) %>%
  summarize(max_chem = length(unique(CAS))) %>%
  filter(!is.na(ww_percent))

count_plot <- ggplot() +
  geom_point(data = chem_sum_3,
             aes(x = ww_percent,
                 y = max_chem)) +
  geom_text_repel(data = filter(chem_sum_3, 
                                max_chem > 20 |
                                  ww_percent > 20),
                  aes(x = ww_percent,
                      y = max_chem,
                      label = shortName)) +
  xlab("Percent of streamflow attributable to WW") +
  theme_minimal() +
  ylab("Number of Chemicals detected (of max sample)") +
  ggtitle("Wastewater % vs Number of chemicals")
count_plot
dir.create("plots",showWarnings = FALSE)
ggsave(count_plot, filename = "plots/count_label.png", width = 11, height = 6)

ww_order <- chem_sum_1 %>%
  select(shortName, ww_percent) %>%
  distinct() %>%
  arrange(desc(ww_percent)) %>%
  pull(shortName)

chem_sum_4 <- chem_sum_1 %>%
  filter(EAR > 0) %>%
  group_by(shortName, ww_percent, date) %>%
  summarize(sum_conc = sum(EAR, na.rm = TRUE),
            n_chems = n()) %>%
  mutate(date = as.numeric(as.factor(date)),
         sample = dense_rank(date)) %>%
  ungroup() %>%
  mutate(shortName = factor(shortName, levels = ww_order)) %>%
  filter(!is.na(ww_percent)) 

conc_sums <- ggplot() +
  geom_point(data = chem_sum_4,size = 3,
           aes(x = shortName, y = sum_conc, group = sample)) +
  ylab("Concentration sums per sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.title.x = element_blank()) +
  ggtitle("Sites ordered by decreasing WW %")

conc_sums
ggsave(conc_sums, filename = "plots/conc_site.png", width = 11, height = 5)

################################################################
# EAR
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

chemical_summary <- chemical_summary %>%
  left_join(df, by="site")

ww_order <- chemical_summary %>%
  select(shortName, ww_percent) %>%
  distinct() %>%
  arrange(desc(ww_percent)) %>%
  pull(shortName)

chem_sum_ear <- chemical_summary %>%
  filter(EAR > 0) %>%
  group_by(shortName, ww_percent, date) %>%
  summarize(sum_conc = sum(EAR, na.rm = TRUE),
            n_chems = n()) %>%
  mutate(date = as.numeric(as.factor(date)),
         sample = dense_rank(date)) %>%
  ungroup() %>%
  mutate(shortName = factor(shortName, levels = ww_order)) %>%
  filter(!is.na(ww_percent)) 

ear_sums <- ggplot() +
  geom_point(data = chem_sum_ear, size = 3,
             aes(x = shortName, y = sum_conc, group = sample)) +
  ylab("EAR sums per sample") +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.title.x = element_blank()) +
  ggtitle("Sites ordered by decreasing WW %")

ear_sums
ggsave(ear_sums, filename = "plots/ear_site_log.png", width = 11, height = 5)

chemical_summary <- chemical_summary %>%
  filter(EAR > 0)

gd <- graph_chem_data(chemical_summary, mean_logic = FALSE) %>%
  left_join(df, by="site") %>%
  left_join(select(tox_list$chem_site, site=SiteID, shortName=`Short Name`), by="site") %>%
  mutate(shortName = factor(shortName, levels = ww_order)) %>%
  filter(!is.na(ww_percent)) 

max_ear_chem <- ggplot() +
  geom_point(data = gd, size = 3,
             aes(x = shortName, y = meanEAR)) +
  ylab("max EAR for each chemical") +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.title.x = element_blank()) +
  ggtitle("Sites ordered by decreasing WW %")

max_ear_chem
ggsave(max_ear_chem, filename = "plots/max_ear_each_chem_site_log.png", width = 11, height = 5)

gd_sum <- gd %>%
  group_by(shortName, ww_percent) %>%
  summarize(sum_max = sum(meanEAR)) 
  

max_ear_sums <- ggplot() +
  geom_point(data = gd_sum, size = 2.5,
             aes(x = shortName, y = sum_max)) +
  ylab("Sum of max EAR") +
  scale_y_log10(limits = c(5e-4,0.3), labels=toxEval:::fancyNumbers) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 15,
                                   hjust = 1, vjust = 0.3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  ggtitle("Sites ordered by decreasing WW %")

max_ear_sums
ggsave(max_ear_sums, filename = "plots/sum_of_max_ear.png", width = 11, height = 5)

gd_max <- gd %>%
  group_by(shortName, ww_percent) %>%
  summarize(max_single = max(meanEAR)) %>%
  ungroup()

max_ear_combo <- ggplot() +
  geom_point(data = gd_sum, size = 3,
             aes(x = shortName, y = sum_max)) +
  geom_point(data = gd_max, size = 3, color="red",
             aes(x = shortName, y = max_single)) +
  ylab("Sum of max EAR") +
  scale_y_log10(limits = c(5e-4,0.3), labels=toxEval:::fancyNumbers) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 15,
                                   hjust = 1, vjust = 0.3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  ggtitle("Sites ordered by decreasing WW %")

max_ear_combo
ggsave(max_ear_combo, filename = "plots/combo_ear.png", width = 11, height = 5)


gd_percent_change <- gd_sum %>%
  left_join(gd_max, by=c("shortName","ww_percent")) %>%
  ungroup() %>%
  mutate(percent_of_max = 100*max_single/sum_max,
         other_percent = 100*(sum_max-max_single)/sum_max)

percent_ear_combo <- ggplot() +
  geom_point(data = gd_percent_change, size = 3,
             aes(x = shortName, y = other_percent)) +
  ylab("Percent of total sum by other chemical [%]") +
  # scale_y_log10(limits = c(5e-4,0.3), labels=toxEval:::fancyNumbers) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 15,
                                   hjust = 1, vjust = 0.3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  ggtitle("Sites ordered by decreasing WW %")

percent_ear_combo
ggsave(percent_ear_combo, filename = "plots/percentage_other.png", width = 11, height = 5)
