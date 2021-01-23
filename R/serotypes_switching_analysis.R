###############################################################################
## Serotypes PMEN3

pmen3_changes <- read.csv("~/Dropbox/phd/elife_paper/data/pmen3_st_changes.csv",
                          stringsAsFactors = FALSE)
count(pmen3_changes, insertion_node)

pmen3_tree <- read.tree(file = "~/Dropbox/phd/PMEN3/gubbins_data/ungapped.PMEN3.node_labelled.final_tree.tre")
pmen3_states <- read.table("~/Dropbox/phd/elife_paper/data/states_res.txt", sep = "\t", comment.char = "",
                           header = TRUE)


pmen3_tree_nodelables <- as.data.frame(pmen3_tree$node.label)
colnames(pmen3_tree_nodelables) <- "node"
pmen3_tree_nodelables <- pmen3_tree_nodelables %>% left_join(pmen3_states)
pmen3_tree$node.label <- pmen3_tree_nodelables$mega

pmen3_tree_tiplables <- as.data.frame(pmen3_tree$tip.label)
colnames(pmen3_tree_tiplables) <- "node"
pmen3_tree_tiplables <- pmen3_tree_tiplables %>% mutate(node = sub("_\\..*$","",node)) %>% left_join(pmen3_states)

ace_groupInfo <- split(pmen3_tree_tiplables$node, pmen3_tree_tiplables$mega)

acer_tree <- groupOTU(pmen3_tree, ace_groupInfo)
ggtree(acer_tree, aes(colour = group),ladderize = FALSE) + geom_nodepoint(aes(colour = label))



anole.tree<-read.tree("http://www.phytools.org/eqg2015/data/anole.tre")
svl <- read.csv("http://www.phytools.org/eqg2015/data/svl.csv",
                row.names=1)
svl <- as.matrix(svl)[,1]
fit <- phytools::fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)

td <- data.frame(node = nodeid(anole.tree, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(anole.tree, d, by = 'node')

p1 <- ggtree(tree, aes(color=trait), layout = 'circular', 
             ladderize = FALSE, continuous = TRUE, size=2) +
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 

p2 <- ggtree(tree, layout='circular', ladderize = FALSE, size=2.8) + 
  geom_tree(aes(color=trait), continuous=T, size=2) +  
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(aes(color=trait), hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(.05, .85)) 

cowplot::plot_grid(p1, p2, ncol=2, labels=c("A", "B"))    



