#########################
# Plot:
#########################
library(toxEval)
library(ggplot2)
library(dplyr)

path_to_file <- 'processed_data/pharm_data.xlsx' 
tox_list <- create_toxEval(path_to_file)

chem_data <- as.data.frame(tox_list$chem_data)


chem_detected <- chem_data %>% 
  group_by(CAS) %>% 
  summarize(conc_max = max(Value)),
            conc_med = median(Value),
            conc_min = median(Value),
            num_samples = sum(!is.na(Value)),
            num_detections = sum(Value>0)) %>%
  arrange(conc_max)
