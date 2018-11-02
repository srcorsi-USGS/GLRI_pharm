library(toxEval)
library(ggplot2)
# Special funtion:
source(file = "combo_graph_function.R")

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

# benchmarks <- tox_list$chem_info %>%
#   select(CAS, parameter_nm) %>%
#   mutate(endPoint = "Concentration",
#          Value = 1,
#          groupCol = "") %>%
#   left_join(distinct(select(toxEval::tox_chemicals,
#                             CAS = Substance_CASRN,
#                             Chemical = Substance_Name)), by="CAS")
# 
# benchmarks$Chemical[is.na(benchmarks$Chemical)] <- benchmarks$parameter_nm[is.na(benchmarks$Chemical)]
# benchmarks$Chemical <- gsub(", water, filtered (0.2 micron filter), recoverable, nanograms per liter", "", benchmarks$Chemical)
# benchmarks$Chemical <- gsub(", nanograms per liter", "", benchmarks$Chemical, fixed = TRUE)
# benchmarks$Chemical <- gsub(", water", "", benchmarks$Chemical, fixed = TRUE)
# benchmarks$Chemical <- gsub(", recoverable", "", benchmarks$Chemical, fixed = TRUE)
# benchmarks$Chemical <- gsub(", filtered", "", benchmarks$Chemical, fixed = TRUE)
# benchmarks$Chemical <- gsub("(0.2 micron filter)", "", benchmarks$Chemical, fixed = TRUE)
# benchmarks$Chemical <- gsub("(0.2 micron filter0", "", benchmarks$Chemical, fixed = TRUE)
# 
# 
# 
# library(openxlsx)
# 
# exclude <- read.csv("raw_data/exclude.csv", stringsAsFactors = FALSE)
# exclude <- select(exclude, CAS, endPoint, chnm=X)
# 
# list_of_datasets <- list("Data" = tox_list$chem_data,
#                          "Chemicals" = tox_list$chem_info,
#                          "Sites" = tox_list$chem_site,
#                          "Exclude" = exclude,
#                          "Benchmarks" = benchmarks)
# 
# write.xlsx(list_of_datasets, file = "processed_data/pharm_data_concentrations.xlsx", append=TRUE)

path_to_file <- 'processed_data/pharm_data_concentrations.xlsx' 
tox_list_conc <- create_toxEval(path_to_file)
chemical_summary_conc <- get_chemical_summary(tox_list_conc)

# Remove non-detects:
chemical_summary <- chemical_summary[chemical_summary$EAR > 0,]
chemical_summary_conc <- chemical_summary_conc[chemical_summary_conc$EAR > 0,]

remaining_classes <- unique(c(unique(as.character(chemical_summary$Class)),unique(as.character(chemical_summary_conc$Class))))
chemical_summary$Class <- droplevels(chemical_summary$Class)
chemical_summary_conc$Class <- droplevels(chemical_summary_conc$Class)
chemical_summary$chnm <- droplevels(chemical_summary$chnm)
chemical_summary_conc$chnm <- droplevels(chemical_summary_conc$chnm)
######################################
conc_plot <- plot_tox_boxplots(chemical_summary_conc, 
                              category = 'Chemical',
                              mean_logic = FALSE,
                              hit_threshold = NA, plot_ND = TRUE)
conc_plot
dir.create("plots", showWarnings = FALSE)
ggsave(conc_plot, filename = "plots/concentration.png", height = 12, width = 12)

graphData_tox <- graph_chem_data(chemical_summary)
graphData_tox$guide_side <- "ToxCast"#expression(atop(ToxCast,Maximum~EAR[Chem]~per~Site))

graphData_conc <- graph_chem_data(chemical_summary_conc, sum_logic = FALSE)
graphData_conc$guide_side <- "Concentration"#expression(Concentration)

graphData_tox$guide_up <- "A"
graphData_conc$guide_up <- "A"

toxPlot_conc <- combo_plot_matches(graphData_tox, graphData_conc, 
                                   thres_1 = NA, thres_2 = NA, 
                                   drop = FALSE, grid = FALSE)
toxPlot_conc
ggsave(toxPlot_conc, filename = "plots/comparison.pdf", width = 11, height = 9)
ggsave(toxPlot_conc, filename = "plots/comparison.png", width = 11, height = 9)


toxPlot_conc_drop <- combo_plot_matches(graphData_tox, graphData_conc, 
                                   thres_1 = NA, thres_2 = NA, 
                                   drop = TRUE, grid = FALSE)
toxPlot_conc_drop
ggsave(toxPlot_conc_drop, filename = "plots/comparison_drop.pdf", width = 11, height = 9)
ggsave(toxPlot_conc_drop, filename = "plots/comparison_drop.png", width = 11, height = 9)
