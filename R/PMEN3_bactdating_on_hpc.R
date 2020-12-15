###############################################################################
## Script to run BactDating on tree with ordered dating csv ###################
###############################################################################

require(devtools)
require(ape)

## set up the input options ##

input_args <- parse_args()

run_install <- input_args[1]
gubbins_loc <- input_args[2]
dating_df <- input_args[3]
tree_model <- input_args[4]
iterations <- as.integer(input_args[5])
output_file <- input_args[6]

## install bactdating 
if(run_install == "Y"){
  devtools::install_github("xavierdidelot/BactDating")
}
require(BactDating)
## Run tree

re_ordered_central_df <- read.csv(dating_df, stringsAsFactors = FALSE)

## Get gubb tree 

total_gub_tree <- loadGubbins(gubbins_loc)

gub_labs <- total_gub_tree$tip.label
split_labs <- str_split_fixed(gub_labs, "_\\.",2)

total_gub_tree$tip.label <- split_labs[, 1]

central_gub_tree <- keep.tip(total_gub_tree, tip = re_ordered_central_df$id)

dec_datio <- re_ordered_central_df$dec_date

central_res <- BactDating::bactdate(ape::unroot(central_gub_tree),
                                                dec_datio, showProgress = T,
                                                nbIts = iterations, model = tree_model, useRec = T)

save(central_res, file =  output_file)

