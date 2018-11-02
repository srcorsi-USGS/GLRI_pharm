#########################
# Plot:
#########################
library(toxEval)
library(ggplot2)
path_to_file <- 'processed_data/pharm_data.xlsx' 
tox_list <- create_toxEval(path_to_file)

Sites <- tox_list$chem_site

labels_watershed <- unique(Sites$Watershed)
unique_space_count <- 1
for(i in names(table(Sites$Watershed))[as.numeric(table(Sites$Watershed)) == 1]){
  labels_watershed[labels_watershed == i] <- paste(rep(" ",unique_space_count),collapse = "")
  unique_space_count <- unique_space_count + 1
}

site_occurances <- tox_list$chem_data %>%
  left_join(select(Sites, SiteID, `Short Name`, Watershed, Lake), by="SiteID") %>%
  group_by(`Short Name`, Watershed, Lake) %>%
  summarize(total_counts = length(unique(CAS[Value > 0]))) %>%
  ungroup() %>%
  mutate(Watershed = factor(Watershed, 
                            levels = unique(Sites$Watershed),
                            labels = labels_watershed))

site_summaries <- tox_list$chem_data %>%
  left_join(select(Sites, SiteID, `Short Name`, Watershed, Lake), by="SiteID") %>%
  group_by(`Short Name`, Watershed, Lake, `Sample Date`) %>%
  summarize(total_concentration = sum(Value)) %>%
  ungroup() %>%
  group_by(`Short Name`, Watershed, Lake) %>%
  summarize(total_concentration = max(total_concentration)) %>%
  ungroup() %>%
  mutate(Watershed = factor(Watershed, 
                            levels = unique(Sites$Watershed),
                            labels = labels_watershed))

occurance_count <- ggplot(data = site_occurances) +
  geom_col(aes(x = `Short Name`, y = total_counts), fill = "steelblue") +
  facet_grid(. ~ Watershed, scales="free", space="free") +
  geom_text(x=1, y=Inf, label  = "Watersheds:", vjust=-0.8,
            data = data.frame(Watershed = site_summaries$Watershed[1])) + 
  ylab("Total number of detected chemicals") +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "inside",
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour="black",fill = "transparent"),
        panel.spacing = unit(0.1, "lines"))
occurance_count

dir.create("plots",showWarnings = FALSE)
ggsave(occurance_count, filename = "plots/count.png", width = 11, height = 5)


concentration_watershed <- ggplot(data = site_summaries) +
  geom_col(aes(x = `Short Name`, y = total_concentration), fill = "steelblue") +
  facet_grid(. ~ Watershed, scales="free", space="free") +
  geom_text(x=1, y=Inf, label  = "Watersheds:", vjust=-0.8,
            data = data.frame(Watershed = site_summaries$Watershed[1])) + 
  ylab(expression(Sigma~"chemicals ["~mu*g/l~"]")) +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.placement = "inside",
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.border=element_rect(colour="black",fill = "transparent"),
        panel.spacing = unit(0.1, "lines"))
concentration_watershed

dir.create("plots",showWarnings = FALSE)
ggsave(concentration_watershed, filename = "plots/sum_conc_site.png", width = 11, height = 5)


chem_summaries <- tox_list$chem_data %>%
  group_by(CAS) %>%
  summarize(total_counts = length(unique(SiteID[Value > 0]))) %>%
  ungroup() %>%
  filter(total_counts > 0,
         !is.na(CAS)) %>%
  left_join(distinct(select(toxEval::tox_chemicals, 
                            CAS = Substance_CASRN, 
                            toxname = Substance_Name)), by="CAS") %>%
  left_join(distinct(select(tox_list$chem_info, CAS, srsname, parameter_nm, Class)), by="CAS") %>%
  arrange(desc(total_counts))

chem_summaries$toxname[is.na(chem_summaries$toxname)] <- chem_summaries$srsname[is.na(chem_summaries$toxname)]
chem_summaries$toxname[is.na(chem_summaries$toxname)] <- chem_summaries$parameter_nm[is.na(chem_summaries$toxname)]

chem_summaries$toxname <- gsub(", water, filtered (0.2 micron filter), recoverable, nanograms per liter", "", chem_summaries$toxname)
chem_summaries$toxname <- gsub(", nanograms per liter", "", chem_summaries$toxname, fixed = TRUE)
chem_summaries$toxname <- gsub(", water", "", chem_summaries$toxname, fixed = TRUE)
chem_summaries$toxname <- gsub(", recoverable", "", chem_summaries$toxname, fixed = TRUE)
chem_summaries$toxname <- gsub(", filtered", "", chem_summaries$toxname, fixed = TRUE)
chem_summaries$toxname <- gsub("(0.2 micron filter)", "", chem_summaries$toxname, fixed = TRUE)
chem_summaries$toxname <- gsub("(0.2 micron filter0", "", chem_summaries$toxname, fixed = TRUE)

chem_summaries$toxname <- factor(chem_summaries$toxname, levels = chem_summaries$toxname)

site_count <- ggplot(data = chem_summaries) +
  geom_col(aes(x = toxname, y = total_counts), fill = "steelblue") +
  ylab("Total number of sites with detections") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 12),
        plot.caption = element_text(size = 14),
        panel.background = element_blank(),
        axis.text.x = element_text(angle=90, size = 14,
                                   hjust = 1, vjust = 0.3)) +
  labs(caption = paste(nrow(Sites),"total sites"))

site_count
ggsave(site_count, filename = "plots/site_counts.png", width = 11, height = 5)

# #https://cida.usgs.gov/hbsl/apex/f?p=104:1:0::NO
# 
# hbsl <- data.table::fread("hbsl_data.csv") %>%
#   data.table::setDF() %>%
#   janitor::clean_names() %>%
#   rename(CAS = cas_registry_number_a,
#          Chemical = chemical_name,
#          pcode = usgs_parameter_code_s,
#          chronic_noncancer = chronic_noncancer_hhbp_µg_l_c_d,
#          mcl = mcl_µg_l_b,
#          noncancer_lde = noncancer_hbsl_µg_l_d_e,
#          carcinogenic_range = carcinogenic_hhbp_10_6_to_10_4_µg_l_c_d,
#          cancer_range = cancer_hbsl_10_6_to_10_4_µg_l_d_e,
#          Class = chemical_class) %>%
#   select(-benchmark_remarks) %>%
#   separate(carcinogenic_range, into = c("carcinogenic_min","carcinogenic_max"), sep = "-", convert = TRUE) %>%
#   separate(cancer_range, into = c("cancer_min","cancer_max"), sep = "-", convert = TRUE) %>%
#   filter(CAS != "") %>%
#   gather(endPoint, Value, -CAS, -Chemical, -pcode, -Class) %>%
#   filter(!is.na(Value)) %>%
#   filter(CAS %in% chem_info$CAS)
# 
# nrow(hbsl)
# 
# hbsl_2 <- readxl::read_xlsx("HBSL_supporting_toxicity_data.xlsx", skip = 3) %>%
#   rename(CAS = `CAS Number`, Chemical = `Chemical Name`) %>%
#   filter(CAS %in% chem_info$CAS,
#          !is.na(Value))

