###############################################################################
## Trying out gg tree for the gubbins figures #################################
###############################################################################

require(dplyr)
require(ggplot2)
require(ggtree)

library(ggtree)
## remote_folder <- paste0("https://raw.githubusercontent.com/katholt/",
##                         "plotTree/master/tree_example_april2015/")
remote_folder <- "~/Dropbox/phd/ggtree_examples/plotTree/tree_example_april2015/" 

## read the phylogenetic tree
tree <- read.tree(paste0(remote_folder, "tree.nwk"))

## read the sampling information data set
info <- read.csv(paste0(remote_folder,"info.csv"))

## read and process the allele table
snps<-read.csv(paste0(remote_folder, "alleles.csv"), header = F,
               row.names = 1, stringsAsFactor = F)
snps_strainCols <- snps[1,] 
snps<-snps[-1,] # drop strain names
colnames(snps) <- snps_strainCols

gapChar <- "?"
snp <- t(snps)
lsnp <- apply(snp, 1, function(x) {
  x != snp[1,] & x != gapChar & snp[1,] != gapChar
})
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]

## read the trait data
bar_data <- read.csv(paste0(remote_folder, "bar.csv"))

## visualize the tree 
p <- ggtree(tree, ladderize = FALSE) 

## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% info + geom_tippoint(aes(color=location))

## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = location), shape = '|') +
  geom_facet(panel = "Trait", data = bar_data, geom = ggstance::geom_barh, 
             aes(x = dummy_bar_value, color = location, fill = location), 
             stat = "identity", width = .6) +
  theme_tree2(legend.position=c(.05, .85))


###############################################################################
## Lets try it with multiple pmen3 datasets ###################################
###############################################################################

pmen3_phandango_data <- read.csv("~/Dropbox/phd/elife_paper/figures_data/pmen3_microreact_pen_res.csv",
                                 stringsAsFactors = FALSE)
pmen3_tree <- read.tree("~/Dropbox/phd/PMEN3/gubbins_data/ungapped.PMEN3.node_labelled.final_tree.tre")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ace_groupInfo <- split(pmen3_phandango_data$id, pmen3_phandango_data$Clade)
lineages <- dplyr::count(pmen3_phandango_data, Clade)
colours_tree <- colorRampPalette(brewer.pal(12,"Paired"))(nrow(lineages) + 1)
colours_tree <- gg_color_hue(nrow(lineages) + 1)
names(colours_tree) <- c(lineages$Lineage,"0")

pmen3_ggtree <- pmen3_phandango_data[,c(1,4:10)]
pmen3_ggtree$bar_value <- 1

ace_tree <- groupOTU(pmen3_tree, ace_groupInfo)
ace_tree <- ggtree(ace_tree, aes(colour = group), layout = "rectangular", ladderize = FALSE) 
ace_tree <- ace_tree +  geom_facet(panel = "Trait", data = pmen3_ggtree, geom = ggstance::geom_barh, 
                       mapping = aes(y = id,x = bar_value, color = Country, fill = Country), 
                       stat = "identity", width = .6)

ace_tree
pmen3_ggtree <- pmen3_phandango_data[,c(1,4:10)]
rownames(pmen3_ggtree) <- pmen3_ggtree$id
pmen3_ggtree <- pmen3_ggtree[,-1]



genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
colnames(genotype) <- sub("\\.", "", colnames(genotype))
genotype[genotype == "trig"] <- "TRIG"
genotype[genotype == "pdm"] <- "Pdm/09"

p <- gheatmap(p, genotype, width=.4, offset=7, colnames=F) %>% 
  scale_x_ggtree








