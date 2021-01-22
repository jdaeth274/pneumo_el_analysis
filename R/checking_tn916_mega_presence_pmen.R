###############################################################################
### Double checking the presence of Mega and Tn916 in the pmen3 and pmen9 col #
###############################################################################

pmen_mega <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_mega_lib_updated_flanks/pmen_mega_merged_blast_file",
                      stringsAsFactors = FALSE)
pmen_tn916 <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_pub_run/pmen_tn916_merged_blast_file",
                       stringsAsFactors = FALSE)

pmen3_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_26_11_2020.csv",
                       stringsAsFactors = FALSE)

pmen9_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_14_12_2020.csv",
                       stringsAsFactors = FALSE)

## PMEN3 

pmen3_ids <- sub("#","_",pmen3_data$id)

pmen3_mega <- pmen_mega[pmen_mega$file_loc %in% pmen3_ids,]

pmen3_mega <- pmen3_mega[pmen3_mega$align > 2000,]
pmen3_mega_ids <- unique(pmen3_mega$file_loc)
pmen3_mega_ids <- str_split_fixed(pmen3_mega_ids, "_", 2)
pmen3_mega_ids <- paste(pmen3_mega_ids[,1], sub("_","#",pmen3_mega_ids[,2]), sep = "_")


pmen3_tn916 <- pmen_tn916 %>% filter(file_loc %in% pmen3_ids) %>% filter(align >= 7000)
pmen3_tn916_ids <- unique(pmen3_tn916$file_loc) 
pmen3_tn916_ids <- str_split_fixed(pmen3_tn916_ids, "_", 2)
pmen3_tn916_ids <- paste(pmen3_tn916_ids[,1], sub("_","#",pmen3_tn916_ids[,2]), sep = "_")

pmen3_orig_ids_mega <- pmen3_data[pmen3_data$MEGA__autocolour == "Yes_MEGA","id"]
pmen3_orig_ids_916 <- pmen3_data[pmen3_data$TN916__autocolour == "Yes_tn916","id"]

length(pmen3_orig_ids_916)
length(pmen3_tn916_ids)

length(pmen3_orig_ids_mega)
length(pmen3_mega_ids)


## PMEN9 ##

pmen9_ids <- sub("#","_",pmen9_data$id)

pmen9_mega <- pmen_mega[pmen_mega$file_loc %in% pmen9_ids,]

pmen9_mega <- pmen9_mega[pmen9_mega$align > 2000,]
pmen9_mega_ids <- unique(pmen9_mega$file_loc)
pmen9_mega_ids <- str_split_fixed(pmen9_mega_ids, "_", 2)
pmen9_mega_ids <- paste(pmen9_mega_ids[,1], sub("_","#",pmen9_mega_ids[,2]), sep = "_")


pmen9_tn916 <- pmen_tn916 %>% filter(file_loc %in% pmen9_ids) %>% filter(align >= 7000)
pmen9_tn916_ids <- unique(pmen9_tn916$file_loc) 
pmen9_tn916_ids <- str_split_fixed(pmen9_tn916_ids, "_", 2)
pmen9_tn916_ids <- paste(pmen9_tn916_ids[,1], sub("_","#",pmen9_tn916_ids[,2]), sep = "_")

pmen9_orig_ids_mega <- pmen9_data[pmen9_data$MEGA_presence__autocolour == "Yes_mega","id"]
pmen9_orig_ids_916 <- pmen9_data[pmen9_data$tn916_presence__autocolour == "Yes_916","id"]

length(pmen9_orig_ids_916)
length(pmen9_tn916_ids)

length(pmen9_orig_ids_mega)
length(pmen9_mega_ids)


pmen9_data <- pmen9_data %>% mutate(Tn916__autocolour = ifelse(id %in% pmen9_tn916_ids, "YES","NO")) %>%
  mutate(Tn1207__autocolour = ifelse(id %in% pmen9_mega_ids, "YES","NO"))
write.csv(pmen9_data, file = "~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_02_01_2021.csv",
          row.names = FALSE)







