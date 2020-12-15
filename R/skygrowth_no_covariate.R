###############################################################################
## Skygrowth no covariate on HPC cluster ######################################
###############################################################################

require(devtools)
require(skygrowth)

parse_args <- function(){
  input_args <- commandArgs(trailingOnly = TRUE)
  
  if(length(input_args) != 4){
    cat("Not enough arguments, need 4, you have:", length(input_args), "\n")
    cat("\n")
    cat("Usage: Rscript --vanilla ./skygrowth_no_covariate.R <dated_tree_obj> <dating_df> <iterations> <output_loc>")
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
output <- input_args[4]

cat("Loading up data \n")
load(tree_object)
timed_tree <- central_res$tree
dec_data <- read.csv(dating_csv, stringsAsFactors = FALSE)
new_labs <- paste(timed_tree$tip.label, dec_data$dec_date, sep = "_") 
timed_tree$tip.label <- new_labs


mcmcfit <- skygrowth.mcmc(timed_tree,res = 100,mhsteps = iterations,
                          control=list(thin=1e5))
par(mfrow=c(1,2))
mcmcfit$time <- (max(dec_data$dec_date) + mcmcfit$time)

save(mcmcfit, file = output)
