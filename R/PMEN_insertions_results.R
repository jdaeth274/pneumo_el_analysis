###############################################################################
## Going through the hits for elements in the PMEN collections ################
###############################################################################

require(dplyr)

hits_combiner <- function(hits, misses){
  multi_hits <- which(misses$id %in% hits$id)
  if(length(multi_hits) > 0){
    misses <- misses[-multi_hits,]
  }
  
  hits_combo <- bind_rows(hits, misses) 
  
}


non_reccy_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_mega_lib_updated_flanks/pmen_mega_non_reccy_hits.csv", 
                           stringsAsFactors = FALSE) %>% rename(isolate_example = isolate_id)
reccy_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_mega_lib_updated_flanks/pmen_mega_reccy_hits.csv", 
                       stringsAsFactors = FALSE)
pmen_mega_hits_df <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_mega_pub_run/pmen_mega_hits_df.csv",
                              stringsAsFactors = FALSE)
pmen_mega_missing_df <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_mega_pub_run/pmen_mega_missing_df.csv",
                                 stringsAsFactors = FALSE)

pmen_mega <- hits_combiner(pmen_mega_hits_df, pmen_mega_missing_df)
total_hits_mega <- dplyr::bind_rows(reccy_hits, non_reccy_hits)

non_reccy_hits_916 <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_pub_run/pmen_tn916_non_reccy_hits.csv", 
                           stringsAsFactors = FALSE) %>% rename(isolate_example = isolate_id)
reccy_hits_916 <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_pub_run/pmen_tn916_reccy_hits.csv", 
                       stringsAsFactors = FALSE)

pmen_tn916_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_pub_run/pmen_tn916_hits_df.csv",
                            stringsAsFactors = FALSE)
pmen_tn916_missing <- read.csv("~/Dropbox/phd/insertion_site_analysis/pmen_tn916_pub_run/pmen_tn916_missing_df.csv",
                            stringsAsFactors = FALSE)

## seems to be duplicate ids in here 
which(pmen_tn916_hits$id %in% pmen_tn916_missing$id)


pmen_tn916 <- hits_combiner(pmen_tn916_hits, pmen_tn916_missing)


total_hits_tn916 <- dplyr::bind_rows(reccy_hits_916, non_reccy_hits_916)


## epi csvs

pmen3_epi_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_26_11_2020.csv",
                           stringsAsFactors = FALSE)

pmen9_epi_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_02_01_2021.csv",
                           stringsAsFactors = FALSE)


## PMEN3 mega 

pmen3_insertions <- total_hits_mega[total_hits_mega$cluster_name == "pmen3",]
pmen3_916_insertions <- total_hits_tn916[total_hits_tn916$cluster_name == "pmen3",]

pmen_mega_presence <- pmen_mega[,c("id"),drop = FALSE] %>% mutate(mega_isa__autocolour = "Yes")
pmen_tn916_presence <- pmen_tn916[,c("id"), drop = FALSE] %>% mutate(tn916_isa__autocolour = "Yes")
pmen3_epi_data <- pmen3_epi_data %>% left_join(y = pmen_mega_presence) %>% 
  mutate(mega_isa__autocolour = ifelse(is.na(mega_isa__autocolour),"No","Yes")) %>%
  left_join( y = pmen_tn916_presence) %>% 
  mutate(tn916_isa__autocolour = ifelse(is.na(tn916_isa__autocolour),"No","Yes"))
count(pmen3_epi_data, MEGA__autocolour)
count(pmen3_epi_data, mega_isa__autocolour)

count(pmen3_epi_data, TN916__autocolour)
count(pmen3_epi_data, tn916_isa__autocolour)

pmen9_epi_data <- pmen9_epi_data %>% left_join( y = pmen_mega_presence) %>% 
  mutate(mega_isa__autocolour = ifelse(is.na(mega_isa__autocolour),"No","Yes")) %>%
  left_join(y = pmen_tn916_presence) %>% 
  mutate(tn916_isa__autocolour = ifelse(is.na(tn916_isa__autocolour),"No","Yes"))
count(pmen9_epi_data, Tn1207__autocolour)
count(pmen9_epi_data, mega_isa__autocolour)
count(pmen9_epi_data, Tn916__autocolour)
count(pmen9_epi_data, tn916_isa__autocolour)


##################################################


write.csv(pmen3_epi_data, file = "~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_19_01_2021.csv",
          row.names = FALSE, quote = FALSE)

write.csv(pmen9_epi_data, file = "~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_19_01_2021.csv",
          row.names = FALSE, quote = FALSE)






