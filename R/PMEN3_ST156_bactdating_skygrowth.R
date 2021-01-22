###############################################################################
## ST156 skygrowth and bactdating #############################################
###############################################################################

require(BactDating)
require(skygrowth)
require(ggplot2)
require(dplyr)
require(ape)
require(cowplot)
require(ggpubr)

###############################################################################
## Load up the BactDating results #############################################
###############################################################################

load("~/Dropbox/phd/elife_paper/figures_data/ten_mil_PMEN3_ST156.RData", verbose = TRUE)
timed_tree <- central_res$tree
dec_data <- read.csv("~/Dropbox/phd/elife_paper/figures_data/ST156_ordered_dating.csv")

new_labs <- paste(timed_tree$tip.label, dec_data$dec_date, sep = "_") 

timed_tree$tip.label <- new_labs

plot(central_res, 'trace')

plot(central_res, 'treeCI')

## lets take the ST156 clade from the whole bactdating 
load("~/Dropbox/phd/PMEN3/whole_tree_relaxed_gamma_bactdating",
     verbose = TRUE)
timed_tree <- central_res$tree
ST156_tree <- keep.tip(timed_tree, dec_data$id)
plot(ST156_tree)

new_labs <- paste(ST156_tree$tip.label, dec_data$dec_date, sep = "_")

ST156_tree$tip.label <- new_labs

###############################################################################
## Confidence intervals are huge!. Lets try and run a skygrowth on this #######
## lineage ####################################################################
###############################################################################


initial_pop_test <- skygrowth.map(tre = ST156_tree)
plot(initial_pop_test)

###############################################################################
## So that looks like its a simple up and down in more recent times, will now #
## run an mcmc chain version to see how this looks ############################
###############################################################################

mcmcfit <- skygrowth.mcmc(ST156_tree,res = 100,mhsteps = 2e+07,
                          control=list(thin=1e5))
par(mfrow=c(1,2))
mcmcfit$time <- (max(dec_data$dec_date) + mcmcfit$time)

plot(mcmcfit) 
growth.plot(mcmcfit)
plot(mcmcfit$tau,type='l',ylab='tau',xlab='Iterations')
plot(mcmcfit,logy=FALSE)

mcmcfit$growthrate

save(mcmcfit, file = "~/Dropbox/phd/PMEN3/R_DATA_analysis/data/skygrowth_data/central_clade_20mill_fit_pub.RData")

sort(dec_data$dec_date)

load("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/german_clade_analysis/skygrowth_162_clade",
     verbose = TRUE)

###############################################################################
## Lets get the two plots together ############################################
###############################################################################

growth_rate_data <- as.data.frame(mcmcfit$growthrate_ci) 
colnames(growth_rate_data) <- c("lower","median","upper")
growth_rate_data$time <- mcmcfit$time

ne_data <- as.data.frame(mcmcfit$ne_ci)
colnames(ne_data) <- c("lower","median","upper")
ne_data$time <- mcmcfit$time

ne_plot <- ggplot(data = ne_data, aes(x = time, y = median)) +
  geom_line(colour = "midnightblue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "midnightblue", colour = "midnightblue", 
              alpha = 0.2) +
  scale_y_log10() + xlim(c(1980,2020)) + labs(x = NULL, y = "Ne")
ne_plot

growth_plot <- ggplot(data = growth_rate_data, aes( x = time, y = median)) +
  geom_line(colour = "midnightblue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), colour = "midnightblue", fill = "midnightblue",
              alpha = 0.2) +
  xlim(c(1980,2020)) + labs(x = "Time",y = "Growth rate")
growth_plot

combined_plot <- plot_grid(ne_plot, growth_plot, align = "v", ncol = 1, nrow = 2,
                           labels = c("B","C"))
combined_plot











