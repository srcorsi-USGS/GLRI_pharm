library(dplyr)
library(toxEval)
library(tidyr)
library(ggplot2)

ear_thresh <- 0.001
siteThres <- 10
# ep_percent_thres <- 0.5

AOP_crosswalk <- read.csv("AOP_crosswalk.csv", stringsAsFactors = FALSE)
AOP <- AOP_crosswalk %>%
  select(endPoint=Component.Endpoint.Name, ID=AOP..) %>%
  distinct()

relevance <- read.csv("AOP_relevance.csv", stringsAsFactors = FALSE) %>% 
  mutate(Relevant = tools::toTitleCase(Relevant)) %>%
  rename(ID=AOP,
         endPoint = Endpoint.s.)

relevance$Relevant <- gsub("NO","No",relevance$Relevant)

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

boxData_max <- chemical_summary %>%
  filter(EAR > 0) %>%
  left_join(AOP, by="endPoint") %>%
  group_by(ID, chnm, CAS, site, date) %>%
  summarize(maxEAR = max(EAR, na.rm = TRUE),
            endPoint_used = endPoint[which.max(EAR)])

boxData <- boxData_max %>%
  group_by(ID,site,date) %>%
  summarize(total = sum(maxEAR))  %>%
  group_by(ID, site) %>%
  summarize(maxMaxEAR = max(total, na.rm = TRUE),
            date_picked = date[which.max(total)]) %>%
  ungroup() %>%
  filter(!is.na(ID)) %>%
  left_join(select(relevance, ID, Relevant), by="ID") %>%
  mutate(Relevant = replace_na(Relevant, "Unknown")) %>%
  mutate(ID = factor(ID),
         Relevant = factor(Relevant, levels = c("Yes","No","Maybe","Unknown"))) 

y_label <- expression(EAR[SiteAOP])
pretty_logs_new <- toxEval:::prettyLogs(boxData$maxMaxEAR)

boxplot_top <- ggplot(data = boxData) +
  geom_boxplot(aes(x=ID, y=maxMaxEAR, fill = Relevant)) +
  theme_bw() +
  xlab("AOP ID") +
  theme(#axis.ticks.x = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 12),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17, hjust = 0),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  scale_fill_manual(values=c("#E69F00", "#009E73", "white","#999999")) +
  scale_y_log10(y_label,
                labels=toxEval:::fancyNumbers,
                breaks=pretty_logs_new)+
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
boxplot_top
ggsave(boxplot_top, filename = "AOP_box.png", width = 11, height = 5)
