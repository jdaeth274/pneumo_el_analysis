###############################################################################
## Skygrowth with covariate on HPC cluster ######################################
###############################################################################

require(devtools)
require(skygrowth)

parse_args <- function(){
  input_args <- commandArgs(trailingOnly = TRUE)
  
  if(length(input_args) != 7){
    cat("Not enough arguments, need 6, you have:", length(input_args), "\n")
    cat("\n")
    cat("Usage: Rscript --vanilla ./skygrowth_covariate.R <dated_tree_obj> <dating_df> <iterations> <covariate_data> <covariate_column> <max_sample_time> <output_loc>")
    cat("\n")
    cat("\n")
    stop("Not enough input files")
  }
  
  return(input_args)
}

input_args <- parse_args()

tree_object <- input_args[1]
dating_csv <- input_args[2]
iterations <- as.integer(input_args[3])
covariate_data <- input_args[4]
covariate_column <- as.integer(input_args[5])
max_sample_time <- as.integer(input_args[6])
output <- input_args[7]

cat("Loading up data \n")
load(tree_object)
timed_tree <- central_res$tree
dec_data <- read.csv(dating_csv, stringsAsFactors = FALSE)
new_labs <- paste(timed_tree$tip.label, dec_data$dec_date, sep = "_") 
timed_tree$tip.label <- new_labs

mac_csv <- read.csv(covariate_data, stringsAsFactors = FALSE)

ecdc_reinert_data_total <- mac_csv[,c(1,covariate_column)] 
colnames(ecdc_reinert_data_total) <- c("time","var")
ecdc_reinert_data_total$var <- scale(ecdc_reinert_data_total$var)

cat("Running Skygrowth \n")


mcmc_covar_skygrowth <- skygrowth.mcmc.covar(timed_tree,~var,
                                             ecdc_reinert_data_total,
                                             maxSampleTime = 2008.5,
                                             res = 100,iter = iterations,
                                             control=list(thin=1e5))

mcmc_covar_skygrowth$time <- (max(dec_data$dec_date) + mcmc_covar_skygrowth$time)

plot(mcmc_covar_skygrowth) 
growth.plot(mcmc_covar_skygrowth)


save(mcmc_covar_skygrowth, 
     file = output)


