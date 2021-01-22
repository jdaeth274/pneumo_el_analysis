###############################################################################
## Plotting the ACE pen reconstruction ########################################
###############################################################################

require(ape)
require(phytools)
require(stringr)
require(didehpc)
require(devtools)
require(ggplot2)

ace_data_format <- function(epi_csv, tree){
  
  tip_labs <- tree$tip.label
  missing_ids <- which(!(epi_csv[,1] %in% tip_labs))
  if(length(missing_ids) > 0)
    epi_csv <- epi_csv[-missing_ids,]
  
  
  character_res <- setNames(as.factor(epi_csv[,2]), epi_csv[,1])
  character_mat <- to.matrix(character_res, levels(character_res))
  ids_no_data <- names(character_res[character_res == "ND"])
  
  if (length(ids_no_data) > 0){
    
    character_mat <- character_mat[,-which(colnames(character_mat) == "ND")]
  }
  ## Now we'll use the frequencies of the remaining values to estimate the ND from.
  
  freqs <- plyr::count(character_res)
  freqs <- freqs[-which(freqs$x == "ND"),]
  freqs_tot <- sum(freqs[,2])
  
  vals <- NULL
  
  for(jj in 1:nrow(freqs)){
    current_num <- freqs[jj,2]
    fraction <- current_num/freqs_tot
    
    vals <- c(vals,fraction)
    
  }
  
  if(length(ids_no_data) > 0){
    
    for(ll in 1:length(ids_no_data)){
      character_mat[ids_no_data[ll], ] <- vals
    }
    
  }  
  cols_character<-setNames(c(palette(rainbow(ncol(character_mat)))),colnames(character_mat))
  cols_character<-setNames(c(palette(rainbow(ncol(character_mat)))),colnames(character_mat))
  
  return(list(matrixio = character_mat, character_cols = cols_character))
}

cluster_run_set_up <- function(tree, epi_data, run_name, col_names_to_run_through){
  
  run_ids <- NULL
  
  for(k in 1:length(col_names_to_run_through)){
    
    current_col_to_use <- col_names_to_run_through[k]
    
    current_epi <- cbind(epi_data$id, epi_data[,which(colnames(epi_data) == current_col_to_use)])  
    orig_pmen3_ace_data <- ace_data_format(current_epi, tree)
    orig_pmen3_ace_mat <- orig_pmen3_ace_data$matrixio
    
    current_run_name <- paste(run_name, current_col_to_use)
    
    character_ace_pmen3_orig <- obj$enqueue(make.simmap(tree, orig_pmen3_ace_mat, model="ER", nsim=100, Q = "mcmc"),
                                            name = current_run_name)
    
    Sys.sleep(5)
    print(character_ace_pmen3_orig$status())
    print(character_ace_pmen3_orig$log())
    
    pmen3_three_cats_whole_tree <- character_ace_pmen3_orig$id
    run_ids <- append(run_ids, pmen3_three_cats_whole_tree)
    
  }
  return(run_ids)
  
}

run_locally <- function(tree, epi_data, run_name, col_names_to_run_through){
  
  run_ids <- NULL
  
  for(k in 1:length(col_names_to_run_through)){
    
    current_col_to_use <- col_names_to_run_through[k]
    
    current_epi <- cbind(epi_data$id, epi_data[,which(colnames(epi_data) == current_col_to_use)])  
    orig_pmen3_ace_data <- ace_data_format(current_epi, tree)
    orig_pmen3_ace_mat <- orig_pmen3_ace_data$matrixio
    
    current_run_name <- paste(run_name, current_col_to_use)
    
    character_ace_pmen3_orig <- make.simmap(tree, orig_pmen3_ace_mat, model="ER", nsim=100, Q = "mcmc")
    
    
    Sys.sleep(5)
    
    
    
    
  }
  return(character_ace_pmen3_orig)
  
}

ace_results_plotter <- function(character_ace, cols_character, input_tree, input_title ){
  
  system.time(character_pd <- summary(character_ace))
  
  plot(sample(character_ace,1)[[1]],cols_character,fsize=0.01,ftype="i",
       ylim=c(-2,Ntip(input_tree)))
  nodelabels(pie= character_pd$ace,piecol=cols_character,cex=0.2)
  add.simmap.legend(colors=cols_character,prompt=FALSE,x=0,y=-4,
                    vertical=FALSE)
  title(main = input_title, line = -1)
  
  character_recon <- list(ace = character_ace, pd = character_pd)
  
  
  
  return(character_recon)
  
}

plot_res <- function(cluster_res_ids, col_names, epi_data, tree, out_pdf){
  
  pdf(file = out_pdf, paper = "A4r", width = 10, height = 8)
  
  for(k in 1:length(cluster_res_ids)){
    cat("\r","On result number:",k)
    pmen3_mn_obj <- obj$task_get(cluster_res_ids[k])
    current_col_to_use <- col_names[k]
    if(pmen3_mn_obj$status() != "COMPLETE"){
      print(paste("Run not complete for:",col_names[k]))
      next
    }
    pmen3_mn_res <- pmen3_mn_obj$result()
    
    ace_data <- cbind(epi_data$id, epi_data[,which(colnames(epi_data) == current_col_to_use)])  
    formatted_ace_data <- ace_data_format(ace_data, tree)
    
    ace_results_plotter(character_ace = pmen3_mn_res, cols_character = formatted_ace_data$character_cols,
                        tree, input_title = current_col_to_use)
    
  }
  
  dev.off()  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


###############################################################################
## Load up the ace data, then choose the colours and plot #####################
###############################################################################
pmen3_micro_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_26_11_2020.csv",
                             stringsAsFactors = FALSE)
load("~/Dropbox/phd/PMEN3/whole_tree_relaxed_gamma_bactdating",
     verbose = TRUE)

pmen3_dated_tree <- central_res$tree

pmen3_dated_ace_res <- ace_data_format(pmen3_micro_data[,c(3,66)], tree = pmen3_dated_tree)

load("~/Dropbox/phd/PMEN3/pen_res/pmen3_pen_cdc_rf_21_10_2020.RData", verbose = TRUE)

pen_cols <- gg_color_hue(4)

pmen3_dated_ace_res$character_cols <- pen_cols[c(1,3,4)]
names(pmen3_dated_ace_res$character_cols) <- c("I","R","S")
pmen3_dated_whole <- ace_results_plotter(pmen3_cdc_rf_res,
                                         pmen3_dated_ace_res$character_cols,
                                         pmen3_dated_tree)#,
#input_title = "PMEN3 Dated Penicillin resistance")
max_h <-max(nodeHeights(pmen3_dated_tree))
tick.spacing<-10
min.tick<-min(nodeHeights(pmen3_dated_tree))
obj<-axis(1,pos=-3,at=seq(max_h,min.tick,by=-tick.spacing),cex.axis=0.5,
          labels=FALSE)
text(obj, rep(-8, length(obj)), max_h - obj, cex = 0.6)


svg(filename = "~/Dropbox/phd/elife_paper/figures/PMEN3_pen_ace.svg", width = 10, height = 8)

pmen3_dated_whole <- ace_results_plotter(pmen3_cdc_rf_res,
                                         pmen3_dated_ace_res$character_cols,
                                         pmen3_dated_tree)#,
#input_title = "PMEN3 Dated Penicillin resistance")

dev.off()

###############################################################################
## Try and get the scale bar and uncertainty around the nodes from the ########
## bactdating plots ###########################################################
###############################################################################

load("~/Dropbox/phd/PMEN3/whole_tree_relaxed_gamma_bactdating",
     verbose = TRUE)
dec_data <- read.csv("~/Dropbox/phd/elife_paper/figures_data/ST156_ordered_dating.csv")
timed_tree <- central_res$tree
ST156_tree <- keep.tip(timed_tree, dec_data$id)
plot(ST156_tree)

## Need to recalibrate the central res so:
##  - The root.time is the central clade time 
##  - The tree is the ST156 tree
##  - The CI is only for the nodes and tips in the ST156 tree 

ST156_tree$node.label
## So we have nodes 82 to 566 in the ST156 tree
which(central_res$tree$tip.label %in% ST156_tree$tip.label)
## And tip nodes 82 to 567!

## Lets change the CI frame now then this has tips then nodes in ordering 
central_res$CI <- central_res$CI[c(82:567,746:1230),]

central_res$tree <- ST156_tree

central_res$rootdate <- 1984.038

central_res$tree$root.time <- 1984.038

plot(central_res, 'treeCI', show.tip.label = FALSE)

bactdate_ci_recon_plot <- function(x, character_ace, cols_character){
  #browser()
  system.time(character_pd <- summary(character_ace))
  show.axis = TRUE
  xl=c(min(x$CI),max(x$CI))-x$tree$root.time
  #plot.phylo(x$tree, x.lim=xl, show.tip.label = FALSE)
  
  plot(sample(character_ace,1)[[1]],cols_character,fsize=0.01,ftype="i",
       ylim=c(-2,Ntip(x$tree)), xlim=xl, mar = c(3, 1, 1, 1))
  nodelabels(pie= character_pd$ace,piecol=cols_character,cex=0.2)
  pre=pretty(xl+x$tree$root.time)
  if (show.axis) axis(1,pre-x$tree$root.time,pre)
  #    axisPhylo(backward = F)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  transblue=grDevices::rgb(0,0,1,0.4)
  transred =grDevices::rgb(1,0,0,0.4)
  transblack = grDevices::rgb(0.5,0.5,0.5,0.75)
  for(i in 1:(Nnode(x$tree)+Ntip(x$tree)))
    if (x$CI[i,1]!=x$CI[i,2])
      lines(x=c(x$CI[i,1],x$CI[i,2])-x$tree$root.time,
            y=rep(obj$yy[i],2),lwd=2,lend=0,
            col=ifelse(i<=Ntip(x$tree),transblue,transred))
  #points(obj$xx[1:x$tree$Nnode+Ntip(x$tree)],
  #       obj$yy[1:x$tree$Nnode+Ntip(x$tree)],pch=19,col="blue",
  #       cex=1)
}

bactdate_ci_hybird_plot <- function(x){
  #browser()
  show.axis = TRUE
  xl=c(min(x$CI),max(x$CI))-x$tree$root.time 
  #xl[1] <- xl[1] - 20 ## add in a further 20 decline to space out tree in graph 
  #plot.phylo(x$tree, x.lim=xl, show.tip.label = FALSE)
  pre=pretty(xl+x$tree$root.time)
  plot.phylo(x$tree, show.tip.label = FALSE, x.lim = c(-5,35), edge.color = "midnightblue")
  if (show.axis) axis(1,pre-x$tree$root.time,pre)
  #    axisPhylo(backward = F)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  transblue=grDevices::rgb(0,0,1,0.4)
  transred =grDevices::rgb(1,0,0,0.75)
  transblack = grDevices::rgb(0.5,0.5,0.5,1)
  for(i in 1:(Nnode(x$tree)+Ntip(x$tree)))
    if (x$CI[i,1]!=x$CI[i,2])
      lines(x=c(x$CI[i,1],x$CI[i,2])-x$tree$root.time,
            y=rep(obj$yy[i],2),lwd=1,lend=0,
            col=ifelse(i<=Ntip(x$tree),transblue,transred))
  #points(obj$xx[1:x$tree$Nnode+Ntip(x$tree)],
  #       obj$yy[1:x$tree$Nnode+Ntip(x$tree)],pch=19,col="blue",
  #       cex=1)
}

bactdate_ci_hybird_plot(central_res)

plot_grid(ne_plot, ncol = 1, nrow = 1, labels = "A")


svg(filename = "~/Dropbox/phd/elife_paper/figures/PMEN3_ST156_bactdating.svg", width = 10, height = 8)
bactdate_ci_hybird_plot(central_res)
dev.off()



###############################################################################
## From the whole noder in the streamline pen change through time script, we ##
## can see that the first resistant node appeared in 1984, this was node ######
## number 746, so this is the number we want for the credible interval ########
###############################################################################

central_res$CI[746,]




