###############################################################################
## Plotting the top hits for each species through flank length for mega and ###
## Tn916 from the GPS collection ##############################################
###############################################################################

require(ggplot2)
require(dplyr)
require(stringr)

###############################################################################
## Lets take in all the results and plot them out #############################
###############################################################################

folder_to_res <- function(results_folder, graph_name, insert_name){
  ## Function to take in the output from the flanks_only_search and output the over flanks graphs 
  ## Results folder should be the main run folder, flanks veccy the main 
  
  folders <- list.files(results_folder, full.names = TRUE, include.dirs = TRUE)
  folders <- sub("//","/",folders)
  whole_blast_res <- folders[grep("_blast_results", folders)]
  before_blast_res <- folders[grep("_before_flank_blast_res", folders)]
  after_blast_res <- folders[grep("_after_flank_blast_res", folders)]
  
  ## whole blast res 
  tic("Running whole flanks sum up")
  whole_res <- NULL
  flanks_veccy <- as.integer(sub("_blast_results","",basename(whole_blast_res)))
  
  for(k in whole_blast_res){
    #next
    current_graphed_res <- files_load_getter(k, genome = "whole")
    whole_res <- bind_rows(whole_res, current_graphed_res)
    
  }
  
  toc()
  ## before blast res 
  tic("Before flanks sum up")
  before_res <- NULL
  flanks_veccy <- as.integer(sub("_before_flank_blast_res","",basename(before_blast_res)))
  for(k in before_blast_res){
    #next
    
    current_graphed_res <- files_load_getter(k, "upstream")
    before_res <- bind_rows(before_res, current_graphed_res)
    
  }
  
  
  toc()
  
  ## after blast res 
  tic("After flanks sum up")
  after_res <- NULL
  flanks_veccy <- as.integer(sub("_after_flank_blast_res","",basename(after_blast_res)))
  for(k in after_blast_res){
    
    current_graphed_res <- files_load_getter(k, "downstream")
    after_res <- bind_rows(after_res, current_graphed_res)
    
  }
  toc()
  
  total_res <- bind_rows(whole_res, before_res, after_res)
  total_graph_df <- reshape2::melt(total_res, id.vars = c("arm","flanks","position"))
  total_graph_df$grouping <- paste(total_graph_df$arm, total_graph_df$variable, sep = "-")
  ## Make a graph 
  
  total_graph <- ggplot(data = total_graph_df, aes(x = flanks, y = value, group = grouping, color = variable)) +
    geom_line(aes(linetype = arm)) + facet_wrap(~ position, ncol = 3) 
  
  
  
  return(list(graph_df = total_graph_df, graph = total_graph))
  
}



files_load_getter <- function(dir_to_files, genome){

  last_character <- base::substr(dir_to_files,
                                 nchar(dir_to_files),
                                 nchar(dir_to_files))
  if(last_character != "/"){
    dir_to_files <- paste(dir_to_files, "/", sep = "")
  }
  summo_csv_file <- list.files(dir_to_files, pattern = "*species_compo*")
  summary_csv <- paste(dir_to_files, summo_csv_file,sep = "")
  
  summo_csv <- read.csv(summary_csv, stringsAsFactors = FALSE)
  blast_results <- list.files(dir_to_files, pattern = "*list*.csv")
  
  counted_summo <- plyr::count(summo_csv$insertion_point)  

  ## Ok so lets use the summo csv with the top species as the top hit
  ## Then we'll create a function to assess what species is what in terms 
  ## of proportion. 

  current_flanks <- as.integer(str_split_fixed(basename(dir_to_files), "_", 2)[1])
  
  current_prop <- species_props(summo_csv, current_flanks, genome = genome)
  
  return(current_prop)
}

species_props <- function(summo_csv, flank_length, genome){

  ## Output table of proportions for the 4 species in the ref db for 
  ## reference and insert isolates
  
  ref_db <- summo_csv[summo_csv$insertion_point == "reference",]
  
  mitis_num_ref <- length(grep("mitis", ref_db$top_species))
  oralis_num_ref <- length(grep("oralis", ref_db$top_species))
  pseudo_num_ref <- length(grep("pseudo", ref_db$top_species))
  pneumo_num_ref <- nrow(ref_db) - (mitis_num_ref + oralis_num_ref + pseudo_num_ref)
  
  
  treatment_db <- summo_csv[summo_csv$insertion_point != "reference",]
  
  mitis_num_tre <- length(grep("mitis", treatment_db$top_species))
  oralis_num_tre <- length(grep("oralis", treatment_db$top_species))
  pseudo_num_tre <- length(grep("pseudo", treatment_db$top_species))
  pneumo_num_tre <- nrow(treatment_db) - (mitis_num_tre + oralis_num_tre + pseudo_num_tre)
  
  prop_df <- data.frame(matrix(ncol = 7, nrow = 2))
  colnames(prop_df) <- c("pneumo","mitis","oralis","pseudopneumoniae", "arm","flanks","position")
  
  props_ref <- c(pneumo_num_ref, mitis_num_ref, oralis_num_ref, pseudo_num_ref) / nrow(ref_db)
  props_tre <- c(pneumo_num_tre, mitis_num_tre, oralis_num_tre, pseudo_num_tre) / nrow(treatment_db)
  
  prop_df[1,1:4] <- props_ref
  prop_df[2,1:4] <- props_tre
  prop_df$arm <- c("reference","treatment")
  prop_df$flanks <- flank_length
  prop_df$position <- genome
  
  return(prop_df)
  
}


tn916_res <- folder_to_res("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/",
              graph_name = "test")
tn916_res$graph


mega_res <- folder_to_res("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/",
                          graph_name = "test")
mega_res$graph


mega_df <- mega_res$graph_df
tn916_df <- tn916_res$graph_df


mega_df$MGE <- "Tn1207.1"
tn916_df$MGE <- "Tn916"


total_mge_df <- bind_rows(mega_df, tn916_df)

head(total_mge_df)

## remove whole to be in line with no whole concat 
total_mge_df <- total_mge_df[total_mge_df$position != "whole",]

total_mge_df$faceting <- paste(total_mge_df$position, total_mge_df$MGE, sep = ".")
total_mge_df$faceting <- factor(total_mge_df$faceting, levels = c("upstream.Tn1207.1",
                                                             "downstream.Tn1207.1",
                                                             "whole.Tn1207.1",
                                                             "upstream.Tn916",
                                                             "downstream.Tn916",
                                                             "whole.Tn916"))

facet_labels <- c(downstream.Tn1207.1 = "Tn1207.1 Downstream",
                  upstream.Tn1207.1 = "Tn1207.1 Upstream",
                  whole.Tn916 = "Tn916 Whole",
                  downstream.Tn916 = "Tn916 Downstream",
                  upstream.Tn916 = "Tn916 Upstream",
                  whole.Tn1207.1 = "Tn1207.1 Whole")

total_graph <- ggplot(data = total_mge_df, aes(x = flanks, y = value, group = grouping, color = variable)) +
  geom_line(aes(linetype = arm)) + facet_wrap(~ faceting,labeller = labeller(faceting = facet_labels)) +
  scale_color_discrete(labels = c("S. pneumoniae","S. mitis","S. oralis", "S. pseudopneumoniae"),
                       name = "Species") +
  scale_linetype_discrete(labels = c("Reference","MGE"), name = "Flank isolate") +
  labs(x = "Flank length", y = "Proportion as top species")
total_graph

ggsave(filename = "~/Dropbox/phd/elife_paper/figures/gps_species_over_flank.png",
       width = 17, height = 12, dpi = 600, units = "cm")


