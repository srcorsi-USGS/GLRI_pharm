#########################
# Plot:
#########################
library(toxEval)
library(ggplot2)
library(dplyr)

path_to_file <- 'processed_data/pharm_data.xlsx' 
tox_list <- create_toxEval(path_to_file)

chem_data <- as.data.frame(tox_list$chem_data)
chem_info <- as.data.frame(tox_list$chem_info)

chem_detected <- chem_data %>% 
  group_by(CAS) %>% 
  summarize(conc_max = max(Value),
            conc_med = median(Value),
            conc_min = median(Value),
            num_samples = sum(!is.na(Value)),
            num_detections = sum(Value>0))

chem_detected <- left_join(chem_info[,c("CAS","parameter_nm")],chem_detected) %>%
  arrange(desc(conc_max))

write.csv(chem_detected,file="./Analysis not pushed/Pharms_detected.csv",row.names = FALSE)
