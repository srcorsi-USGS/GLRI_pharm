library(toxEval)
library(dplyr)
#NOTE: Add path to path_to_file!!!
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
chemical_summary <- filter(chemical_summary,EAR>0)

######################################
bio_plot <- plot_tox_boxplots(chemical_summary, 
                              category = 'Chemical',
                              mean_logic = FALSE,font_size = 14,
                              hit_threshold = NA,
                              title = 'Summing EARs for all chemicals of a sample, taking the max of each site',
                              plot_ND = TRUE)
bio_plot
png("plots/detected_chem_boxplots.png", width = 1000, height = 600, res = 100)
grid::grid.draw(bio_plot)
dev.off()
