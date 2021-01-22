###############################################################################
## Simpsons diversity index on insertion types in GPS collection ##############
###############################################################################

require(vegan)


tn916_data <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/tn916_res/gps_tn916_run_free/gps_tn916_hits_df.csv",
                       stringsAsFactors = FALSE)
mega_data <- read.csv("~/Dropbox/phd/insertion_site_analysis/data/gps_run_data/mega_res/gps_mega_run4/gps_mega_hits_df.csv",
                      stringsAsFactors = FALSE)

colnames(tn916_data)
diversity(tn916_data$insert_name, index = "simpson")
diversity(mega_data$insert_name, index = "simpson")

## manual diversity ##
simpsons_function <- function(insert_data, by_cluster = FALSE){

  if(by_cluster){
    tn916_simpsons <- NULL
    clusters <- unique(insert_data$cluster_name)
    for(cluster in clusters){
      current_dat <- insert_data[insert_data$cluster_name == cluster,]
      tn916_counter <- current_dat %>% count(insert_name) %>% mutate(numerator = n*(n-1))
      tn916_denominator <- sum(tn916_counter[,2]) * (sum(tn916_counter[,2]) - 1) 
      tn916_numerator <- sum(tn916_counter[,3])
      
      current_simpsons <- 1 - (tn916_numerator/tn916_denominator)
      if(is.nan(current_simpsons))
        current_simpsons <- 0
      
      row_out <- cbind.data.frame(cluster, current_simpsons)
      colnames(row_out) <- c("cluster","index")
      tn916_simpsons <- bind_rows(tn916_simpsons, row_out)
    }
  }else{
  
    tn916_counter <- insert_data %>% count(insert_name) %>% mutate(numerator = n*(n-1))
    tn916_denominator <- sum(tn916_counter[,2]) * (sum(tn916_counter[,2]) - 1) 
    tn916_numerator <- sum(tn916_counter[,3])
    
    tn916_simpsons <- 1 - (tn916_numerator/tn916_denominator)
  }
  return(tn916_simpsons)
}

by_cluster <- simpsons_function(tn916_data, by_cluster = TRUE)
by_cluster$x_val <- "Yes"

ggplot(data = by_cluster, aes(x = x_val, y = index)) + geom_violin()

tn916_counter <- tn916_data %>% count(insert_name) %>% mutate(numerator = n*(n-1))
tn916_denominator <- sum(tn916_counter[,2]) * (sum(tn916_counter[,2]) - 1) 
tn916_numerator <- sum(tn916_counter[,3])

tn916_simpsons <- 1 - (tn916_numerator/tn916_denominator)


by_cluster_mega <- simpsons_function(mega_data, by_cluster = TRUE)
by_cluster_mega$x_val <- "Yes"

ggplot(data = by_cluster_mega, aes(x = x_val, y = index)) + geom_violin() + geom_point(position = position_jitter(width = 0.2))


mega_counter <- mega_data %>% count(insert_name) %>% mutate(numerator = n*(n-1))
mega_denominator <- sum(mega_counter[,2]) * (sum(mega_counter[,2]) - 1) 
mega_numerator <- sum(mega_counter[,3])

mega_simpsons <- 1 - (mega_numerator/mega_denominator)







