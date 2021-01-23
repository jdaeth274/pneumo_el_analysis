###############################################################################
## Rscript to collate cluster info into one csv ###############################
###############################################################################

require(stringr, quietly = TRUE)

## read in the cluster fasta_list 

cluster_fastas <- readLines("./cluster_fasta_list.txt")
gps_num_init <- str_split_fixed(cluster_fastas,"cluster_",2)[,2]
gps_nums <- str_split_fixed(gps_num_init, "_list", 2)[,1]

gpsc_name <- paste("gpsc.", gps_nums, sep = "")

lane_ids <- sub("\\.velvet\\.fasta","",basename(cluster_fastas))

out_csv <- cbind.data.frame(lane_ids, gpsc_name)
colnames(out_csv) <- c("Lane_name", "gpsc")

write.csv(out_csv, file = paste("./",gpsc_name[1],"_fastas_list.csv",sep = ""),
          row.names = FALSE)
