library(dplyr)
library(tidyr)

#######################
# Bringing in dataset #
#######################
#Import and give col names to RDB file
# headers<-read.csv("M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2018/Data analysis/Raw/NWISReportingToolDownloads/rdb", skip=104, header=FALSE, nrows=1, as.is=TRUE, sep="\t")
# headers<-unname(unlist(headers[1,]))
# df<-read.csv("M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2018/Data analysis/Raw/NWISReportingToolDownloads/rdb", header=FALSE, as.is=TRUE,skip=106, sep="\t", colClasses=c("V4"="character","V10"="character","V14"="character"))
# names(df)<-headers
# rm(headers)
# df$QW_ANL_SCHED<-substr(df$QW_ANL_SET_NO,1,4)
# df$QW_STATE_DB_RECNO<-paste(df$SITEFILE_NWIS_HOST_NM,df$QW_SAMPLE_DB_NO,df$QW_RECORD_NO,sep="_")
# 
# site_selection <- read.csv("M:/QW Monitoring Team/GLRI Toxics/GLRI II/WY2018/WY2018 planning/Site selection/WW influence Bar chart/Initial Prospective Study Sites for WY18 Pharm and PFAS SamplingSRC.csv", stringsAsFactors = FALSE)
# saveRDS(site_selection, "site_selection.rds")
# 
# #Drop samples not done with pharm schedule
# df <- df %>%
#   filter(QW_ANL_SCHED==2440,
#          QW_SAMP_TYPE_NM %in% c("Regular","Composite (time)"))
# 
# saveRDS(df, file = "raw_data/raw_data.rds")

#########################################
# Start from local copy:
#########################################
site_selection <- readRDS("raw_data/site_selection.rds")
df <- readRDS("raw_data/raw_data.rds")

#########################
# Get Site tab ready:
#########################
site_nos <- df %>%
  select(site_no = SITE_NO) %>%
  distinct() %>%
  mutate(site_no = gsub(" ","", site_no)) %>%
  pull(site_no)

siteInfo <- dataRetrieval::readNWISsite(siteNumbers = site_nos)

Sites <- siteInfo %>%
  select(SiteID = site_no, station_nm, 
         dec_lon = dec_long_va, dec_lat = dec_lat_va) %>%
  arrange(dec_lon) %>%
  left_join(select(site_selection, SiteID=STAID, `Short Name` = Short.Name, Watershed, Lake), by="SiteID")

Sites$`Short Name`[Sites$SiteID == "430128087533602"] <- "Milw Outer Harbor"
Sites$Watershed[Sites$SiteID == "430128087533602"] <- "Milwaukee"
Sites$Lake[Sites$SiteID == "430128087533602"] <- "Michigan"

Sites$`Short Name`[Sites$SiteID == "04093503"] <- "Burns at 20"
Sites$Watershed[Sites$SiteID == "04093503"] <- "Burns"
Sites$Lake[Sites$SiteID == "04093503"] <- "Michigan"

Sites$`Short Name`[Sites$SiteID == "413459087081301"] <- "Salt Creek"
Sites$Watershed[Sites$SiteID == "413459087081301"] <- "Burns"
Sites$Lake[Sites$SiteID == "413459087081301"] <- "Michigan"

Sites$`Short Name`[Sites$SiteID == "412227083344400"] <- "N Br Portage US"
Sites$Watershed[Sites$SiteID == "412227083344400"] <- "Portage"
Sites$Lake[Sites$SiteID == "412227083344400"] <- "Erie"

Sites$`Short Name`[Sites$SiteID == "04195000"] <- "N Br Portage DS"
Sites$Watershed[Sites$SiteID == "04195000"] <- "Portage"
Sites$Lake[Sites$SiteID == "04195000"] <- "Erie"

Sites$`Short Name`[Sites$SiteID == "04161908"] <- "Clinton Sterling"
Sites$Watershed[Sites$SiteID == "04161908"] <- "Clinton"
Sites$Lake[Sites$SiteID == "04161908"] <- "Lake St. Clair"

Sites$`Short Name`[Sites$SiteID == "04237500"] <- "Seneca Baldwinsville"
Sites$Watershed[Sites$SiteID == "04237500"] <- "Oswego"
Sites$Lake[Sites$SiteID == "04237500"] <- "Ontario"

rm(siteInfo, site_selection)
#########################
# Get Data tab ready:
#########################
chem_data <- df %>%
  select(SiteID = SITE_NO, pcode = QW_PARM_CD,
         `Sample Date` = QW_SAMPLE_START_LOCAL_DISP_FM,
         chem_nm = QW_PARM_NM, remark_cd = QW_REMARK_CD, Value = RESULT_VA_TX,
         report_lv_cd = QW_RPT_LEV_CD, reporting_level = QW_RPT_LEV_VA_TX) %>%
  mutate(`Sample Date` = as.POSIXct(strptime(as.character(`Sample Date`), format = "%Y%m%d%H%M")),
         pcode = as.character(pcode))

#########################
# Get Chemical tab ready:
#########################
chem_info <- dataRetrieval::readNWISpCode(unique(chem_data$pcode))

# Get Paul Bradley's classes:
# paul_class <- readxl::read_xlsx("SESQA MIX.TOXEVAL.NO DEGRADATES.4-COL.2017 06 23.xlsx", sheet = "Chemicals")
# 
# chem_info <- chem_info %>%
#   rename(CAS = casrn) %>%
#   left_join(paul_class, by="CAS")

# Original Fulong paper:
# Therapeutic <- readxl::read_excel("Furlong 2017 tables 1 and 2.xlsx", sheet = "Table 2")
# Class_table <- Therapeutic[-1,] %>%
#   select(CAS = CASRNa, Class = `Therapeutic use`) %>%
#   distinct() %>%
#   filter(!is.na(CAS),
#          !duplicated(CAS))

# Just for us:
Therapeutic <- readxl::read_excel("raw_data/2440 and other compounds_Use_Kow'.xlsx", 
                                  sheet = "T2-Phase II Pharmaceuticals",skip = 3)

Class_table <- Therapeutic %>%
  select(CAS = CASRNa, Class = `Broad Pharmacological Activity Category`) %>%
  filter(!is.na(CAS),
         !duplicated(CAS)) %>%
  distinct()

Class_table$CAS <- gsub("–","-",Class_table$CAS)
Class_table$CAS <- gsub("-","-",Class_table$CAS)
Class_table$CAS[Class_table$CAS == "66357-59-3"] <- "66357-35-5"
Class_table$CAS[Class_table$CAS == "56392-17-7"] <- "51384-51-1"
Class_table$CAS[Class_table$CAS == "33286-22-5"] <- "42399-41-7"
  
chem_info <- chem_info %>%
  rename(CAS = casrn) %>%
  left_join(Class_table, by="CAS")
#########################
# Cleanup:
#########################

find_surrogates <- grep("surrogate", chem_info$parameter_nm, ignore.case = TRUE)
chem_info <- chem_info[-find_surrogates,]

# Drop non-pharms:
chem_info <- chem_info %>%
  filter(!(CAS %in% c("1912-24-9","29385-43-1","51-03-6")))

# Deal with chems without CAS:
chem_info$CAS[chem_info$parameter_cd == "67999"] <- "67999"
chem_info$CAS[chem_info$parameter_cd == "67512"] <- "67512"
chem_info$CAS[chem_info$parameter_cd == "67460"] <- "67460"

# c("100-97-0","132-22-9","51481-61-9","141-83-3","846-50-4",
#   "60142-96-3","54910-89-3","79617-96-2","10540-29-1",
#   "65277-42-1","61869-08-7","60-87-7","50-48-6")

chem_info$Class[chem_info$CAS == "97825-25-7"] <- "Veterinary Growth Promoter"
chem_info$Class[chem_info$CAS == "100-97-0"] <- "Antibiotic"
chem_info$Class[chem_info$CAS == "132-22-9"] <- "OTC-chronic condition"
chem_info$Class[chem_info$CAS == "51481-61-9"] <- "OTC-Chronic Condition" #Antacid
chem_info$Class[chem_info$CAS == "141-83-3"] <- "Degradate" #Antacid
chem_info$Class[chem_info$CAS == "846-50-4"] <- "Antidepressant-Neurochemical Modulation" #Sleep medicine
chem_info$Class[chem_info$CAS == "60142-96-3"] <- "OTC-Chronic Condition" #Seizures
chem_info$Class[chem_info$CAS == "54910-89-3"] <- "Antidepressant-Neurochemical Modulation" #increasing the amount of serotonin
chem_info$Class[chem_info$CAS == "79617-96-2"] <- "Antidepressant-Neurochemical Modulation" #panic attacks
chem_info$Class[chem_info$CAS == "10540-29-1"] <- "Chemotherapeutic"
chem_info$Class[chem_info$CAS == "65277-42-1"] <- "Antifungal"
chem_info$Class[chem_info$CAS == "61869-08-7"] <- "Antidepressant-Neurochemical Modulation"
chem_info$Class[chem_info$CAS == "60-87-7"] <- "Antihistamine"
chem_info$Class[chem_info$CAS == "50-48-6"] <- "Antidepressant-Neurochemical Modulation"
chem_info$Class[chem_info$CAS == "4205-90-7"] <- "Cardiovascular Care"
chem_info$Class[chem_info$CAS == "186691-13-4"] <- "Asthma Relief" #COPD?
chem_info$Class[chem_info$CAS == "76963-41-2"] <- "OTC-chronic condition"
chem_info$Class[chem_info$CAS == "52-53-9"] <- "Cardiovascular Care"
chem_info$Class[chem_info$CAS == "300-62-9"] <- "Antidepressant-Neurochemical Modulation"
chem_info$Class[chem_info$CAS == "1159-82-6"] <- "Degradate"
chem_info$Class[chem_info$CAS == "130-95-0"] <- "Antimalarial, flavoring"
chem_info$Class[chem_info$CAS == "67512"] <- "OTC-Chronic Condition"
chem_info$Class[chem_info$CAS == "67460"] <- "OTC-Chronic Condition"
chem_info$Class[chem_info$CAS == "67999"] <- "Degradate"

# chem_info$Class[is.na(chem_info$Class)] <- "Unknown"

chem_data <- chem_data %>%
  filter(pcode %in% chem_info$parameter_cd) %>%
  left_join(select(chem_info, CAS, pcode=parameter_cd, units = parameter_units), by="pcode") %>%
  select(CAS, SiteID, `Sample Date`, Value, remark_cd, DL = reporting_level, units)

chem_data$Value[chem_data$units == "ng/l"] <- chem_data$Value[chem_data$units == "ng/l"]/1000
chem_data$units[chem_data$units == "ng/l"] <- "ug/l"
chem_data$SiteID <- gsub(" ","", chem_data$SiteID)
chem_info$parameter_units[chem_info$parameter_units == "ng/l"] <- "ug/l"

chem_data$Value[chem_data$remark_cd == "<"] <- 0 

Sites$site_grouping <- Sites$Watershed

all((unique(chem_data$SiteID) %in% Sites$SiteID))

#########################
# Save:
#########################
library(openxlsx)

exclude <- read.csv("raw_data/exclude.csv", stringsAsFactors = FALSE)
exclude <- select(exclude, CAS, endPoint, chnm=X)

list_of_datasets <- list("Data" = chem_data, 
                         "Chemicals" = chem_info,
                         "Sites" = Sites,
                         "Exclude" = exclude)
dir.create("processed_data", showWarnings = FALSE)
write.xlsx(list_of_datasets, file = "processed_data/pharm_data.xlsx", append=TRUE)

