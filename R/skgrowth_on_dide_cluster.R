###############################################################################
## Running skygrowth on cluster ###############################################
###############################################################################

require(ape)
require(phytools)
require(stringr)
require(didehpc)
require(devtools)
require(ggplot2)
require(skygrowth)
## Watch out this is the mrc-ide version of the package, there is a package 
## with the same name on CRAN that doesn't work
remotes::install_github("mrc-ide/buildr")
require(buildr)

## HAve to build the windows binary myself, use the altered source file from github!

#devtools::check_win_release(pkg = "/home/jd2117/Downloads/BactDating/", email = "jd2117@ic.ac.uk")
#devtools::build(pkg = "~/Documents/BactDating/", binary = TRUE)
#devtools::build(pkg = "~/Documents/BactDating/", binary = FALSE)

setwd("~/HOMES_drive2/")
options(didehpc.username = "jd2117",didehpc.home = "~/HOMES_drive2/",didehpc.cluster = "fi--didemrchnb")
didehpc::didehpc_config(cores = 3,parallel = FALSE)

context::context_log_start()

root <- "./skygrowth_cluster_no_vignettes_zip"


ctx<- context::context_save(root, packages = c("skygrowth", "phytools","ggplot2","ape"),
                            package_sources = provisionr::package_sources(#local = "~/Downloads/BactDating_1.0.12_knitr.zip"))#,
                              #github = "mrc-ide/skygrowth"))#, ## Just github 
                              local = "~/Downloads/skygrowth_0.3.1_no_vignettes.zip")) ## Windows binary))
                              #local_drat = "~/Downloads/BactDating/"))

## From this it looks like either the local knitr zip or no_vignettes zip is the best move forward
## Seems to get to the cross stage and then conks out. Lets try with skygrwoth quickly as well. 

config <- didehpc::didehpc_config(cores = 1, parallel = FALSE)

system.time(obj <- didehpc::queue_didehpc(ctx,config))

## simple check 
didehpc::web_login()
t <- obj$enqueue(sessionInfo())

