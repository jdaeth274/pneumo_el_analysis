###############################################################################
## Making S1 data table, need to combine GPS collection and PMEN collection ###
###############################################################################

require(dplyr)

## get the total gps accessions then remove the PMEN from this collection 

total_accessions_data <- read.csv("~/Dropbox/phd/elife_paper/data/all_gps_isolates_t2_gladstone.csv",
                                  stringsAsFactors = FALSE)

pmen_accession_data <- read.table("~/Dropbox/phd/elife_paper/data/pmen_accessions.tsv", 
                                  header = TRUE)


not_in_gps_sanger_id <- which(!(pmen_accession_data$Lane_name %in% total_accessions_data$Taxon))
not_in_gps_ERR_id <- which(!(pmen_accession_data$Lane_accession %in% total_accessions_data$ERR))

## Not in GPS published PMEN isolates 

pmen_non_gps_pub <- pmen_accession_data[not_in_gps_ERR_id,]

## Not in GPS at all 

pmen_non_gps_at_all <- pmen_non_gps_pub[-grep("STDY",pmen_non_gps_pub$Sample_name),]
