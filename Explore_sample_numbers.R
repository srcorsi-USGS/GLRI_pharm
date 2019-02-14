# Determine how many chems show up in different situations
# How many in toxcast, detections, ...
#
library(toxEval)
library(USGSHydroTools)
library(dplyr)
library(readxl)


path_to_file <- 'processed_data/pharm_data.xlsx' 
tox_list <- create_toxEval(path_to_file)
chems <- tox_list$chem_info
dim(chems)

conc <- tox_list$chem_data
samples <- paste(conc$SiteID,conc$`Sample Date`)
unique(samples) #107 total samples in data set

unique(conc$remark_cd)

conc <- filter(conc,Value > 0)

unique(conc$CAS) #57 chemicals detected


ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep)

chemicalSummary <- get_chemical_summary(tox_list, ACClong, filtered_ep)
unique(chemicalSummary$chnm) #46 total chemicals in toxcast

chemicalSummary <- chemicalSummary %>%
  filter(EAR > 0)
unique(chemicalSummary$chnm) #22 detected chemicals in toxcast

length(unique(chemicalSummary$chnm)) #How many detected chems in toxcast = 22


# which site has max EAR for each chem
max_EARchem <- chemicalSummary %>%
  group_by(site,shortName,date,chnm) %>%
  summarize(EAR=sum(EAR)) %>%
  filter(EAR > 0) %>%
  group_by(chnm) %>%
  filter(EAR == max(EAR))

methenamine <- filter(chemicalSummary,chnm=="Methenamine")
max(methenamine$EAR)

#Metformin
metformin_CAS <- tox_list$chem_info$CAS[grep("Metformin",tox_list$chem_info$parameter_nm)]
metformin <- filter(tox_list$chem_data,CAS==metformin_CAS)
summary(metformin$Value)

ACCMetformin <- filter(ACClong,CAS == metformin_CAS)


