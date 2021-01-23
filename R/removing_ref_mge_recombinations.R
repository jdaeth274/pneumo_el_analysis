###############################################################################
## removing the mge recombinations from the overall r/m data ##################
###############################################################################
require(dplyr)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

tn916_refs <- readLines("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_realtered_act_refs.txt")
mega_refs <- readLines("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_realtered_act_refs.txt")



tn916_merged_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_merged_blast_file",
                              stringsAsFactors = FALSE) %>% mutate(subject = sub("_(?!_)[0-9]*$","",subject, perl = TRUE))
mega_merged_hits <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_merged_blast_file",
                             stringsAsFactors = FALSE) %>% mutate(subject = sub("_(?!_)[0-9]*$","",subject, perl = TRUE))

ref_iso_fasta <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/gps_reference_isolate_fasta.csv",
                          stringsAsFactors = FALSE) %>% mutate(reference = sub("\\.velvet\\.fasta","",basename(reference))) %>%
  select(reference, cluster_name) %>% aggregate(reference ~ cluster_name,data = .,getmode)


tn916_refs_db <- tn916_merged_hits[grep(paste(paste(tn916_refs, "$",sep = ""), collapse = "|"), tn916_merged_hits$subject),]
mega_refs_db <- mega_merged_hits[grep(paste(paste(mega_refs, "$",sep = ""), collapse = "|"), mega_merged_hits$subject),]


## so now we have the mges in the references, lets find what clusters the references represent and work from there 


recombination_lengths <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/gps_recombination_lengths.csv",
                                  stringsAsFactors = FALSE)

removing_mges <- function(recombination_lengths, ref_iso_fasta, tn916_merged_hits, mega_merged_hits){

  revamped_df <- NULL
  for(cluster in unique(recombination_lengths$cluster)){
    data_subset <- recombination_lengths[recombination_lengths$cluster == cluster,]
    reference_id <- ref_iso_fasta[ref_iso_fasta$cluster_name == cluster, "reference"]
    if(reference_id %in% tn916_merged_hits$subject){
      tn916_row <- tn916_merged_hits[tn916_merged_hits$subject == reference_id,]
      tn916_loc <- sort(as.numeric(tn916_row[1,6:7]))
      ## remove completely enveloping recombinations
      envelops <- which(data_subset$start <= tn916_loc[1] & data_subset$end >= tn916_loc[2])
      ## upstream overhang
      upstream <- which(data_subset$start >= tn916_loc[1] & data_subset$end >= tn916_loc[2] &
                          data_subset$start <= tn916_loc[2])
      ## downstream overhang
      downstream <- which(data_subset$start <= tn916_loc[1] & data_subset$end <= tn916_loc[2] &
                            data_subset$end >= tn916_loc[1])
      ## within tn916
      within <- which(data_subset$start >= tn916_loc[1] & data_subset$end <= tn916_loc[2])
      
      overlappers <- unique(c(envelops, upstream, downstream, within))
      
      if(length(overlappers) > 0){
        data_subset <- data_subset[-overlappers,]
      }
      
    }
    if(reference_id %in% mega_merged_hits$subject){
      mega_row <- mega_merged_hits[mega_merged_hits$subject == reference_id,]
      mega_loc <- sort(as.numeric(mega_row[1,6:7]))
      ## remove completely enveloping recombinations
      envelops <- which(data_subset$start <= mega_loc[1] & data_subset$end >= mega_loc[2])
      ## upstream overhang
      upstream <- which(data_subset$start >= mega_loc[1] & data_subset$end >= mega_loc[2] &
                          data_subset$start <= mega_loc[2])
      ## downstream overhang
      downstream <- which(data_subset$start <= mega_loc[1] & data_subset$end <= mega_loc[2] &
                            data_subset$end >= mega_loc[1])
      ## within mega
      within <- which(data_subset$start >= mega_loc[1] & data_subset$end <= mega_loc[2])
      
      overlappers <- unique(c(envelops, upstream, downstream, within))
      
      if(length(overlappers) > 0){
        data_subset <- data_subset[-overlappers,]
      }
      
    }
    revamped_df <- bind_rows(revamped_df, data_subset)
    
  }
  
  lost_rows <- nrow(recombination_lengths) - nrow(revamped_df)
  lost_MGEs <- nrow(recombination_lengths[recombination_lengths$MGE == "Yes",]) - nrow(revamped_df[revamped_df$MGE == "Yes",])
  cat(paste("Lost this many rows in total: ", lost_rows, "\n", sep = "\n"))
  cat(paste("Of which this many from MGE recombinations: ", lost_MGEs, "\n", sep = "\n"))
  
  return(revamped_df)
  
}

sans_mges <- removing_mges(recombination_lengths, ref_iso_fasta, tn916_merged_hits, mega_merged_hits)

mge_sans_mges <- sans_mges[sans_mges$MGE == "Yes",]
nonmge_sans_mges <- sans_mges[sans_mges$MGE == "No",]

median(mge_sans_mges$length)
median(nonmge_sans_mges$length)
wilcox.test(mge_sans_mges$length, nonmge_sans_mges$length)



histo_mges <- ggplot(data = sans_mges) + geom_histogram(aes(sans_mges$density,
                                                          fill = MGE), alpha  = 1,
                                                      colour = "black") +
  scale_x_log10(limits = c(0.0001, 0.1)) + scale_y_log10()+
  labs(y = "Count", x = "Density") + theme(legend.position = "none")

histo_mges_length <- ggplot(data = sans_mges) + geom_histogram(aes(sans_mges$length,
                                                                 fill = MGE), alpha = 1,
                                                             colour = "black") +
  scale_x_log10(limits = c(10, 1e6)) + scale_y_log10() +
  labs( y = "Count", x = "Bases") + scale_fill_discrete(breaks = c("Yes","No")) + theme(legend.position = "none")

compo_dot_plot_mges <- ggplot(data = sans_mges) + geom_point(aes(x = length, y = density, colour = MGE, alpha = MGE, fill = MGE)) + 
  scale_colour_discrete(breaks = c("Yes","No")) + scale_alpha_manual(values = c("Yes" = 1,"No" = 0.05), breaks = c("Yes","No")) +
  scale_fill_discrete(breaks = c("Yes","No"))


compo_dot_plot_log_mges <- compo_dot_plot_mges + scale_x_log10(limits = c(10,1e6)) + scale_y_log10(limits = c(0.0001, 0.1)) +
  labs( y = "Density", x = "Bases") + 
  stat_density2d(aes(x = length, y = density)) 

wilcox.test()
out_png_name <- "~/Dropbox/phd/elife_paper/figures/gps_reccy_lengths_sans_mge.png"
png(filename = out_png_name, width = 17, height = 15, units = "cm", res = 1000)
# ggdraw() +
#   draw_plot(histo_mges_length, x = 0, y = .5, width = .5, height = .5) +
#   draw_plot(histo_mges, x = .5, y = .5, width = .5, height = .5) +
#   draw_plot(compo_dot_plot_log_mges, x = 0, y = 0, width = 1, height = 0.5) +
#   draw_plot_label(label = c("A", "B", "C"), size = 15,
#                   x = c(0, 0.5, 0), y = c(1, 1, 0.5))
ggarrange(compo_dot_plot_log_mges,
          ggarrange(histo_mges, histo_mges_length, ncol = 2, labels = c("B", "C")),
          nrow = 2, labels = "A", common.legend = TRUE, legend = "bottom") 

dev.off()


