###############################################################################
## Circular plots for penicillin resistance among other things ################
###############################################################################

require(ggtree)
pmen3_tree <- read.tree("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/microreact_total_tree_epi_data.nwk")
pmen3_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_26_11_2020.csv",
                                  stringsAsFactors = FALSE) %>% rename(Resistance = cdc_only_RF_SIR__autocolour)

pmen9_tree <- read.tree("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/altered_branch_length_575_tree.nwk")
pmen9_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_14_12_2020.csv",
                       stringsAsFactors = FALSE) %>% rename(Resistance = cdc_only_RF_SIR__autocolour)

pmen3_ace_data <- pmen3_data %>% rename(Resistance = cdc_only_RF_SIR__autocolour)


ace_data2 <- acinetobacter_lineage %>% mutate(Lineage = paste("GC",as.character(combined_Cluster__autocolour), sep = "")) %>%
  mutate(Lineage = ifelse(Lineage == "GC1", "GC-2", as.character(Lineage))) %>%
  mutate(Lineage = ifelse(Lineage == "GC2", "GC-1", Lineage)) %>% mutate(Lineage = ifelse(grepl("-",Lineage),
                                                                                          sub("-","",Lineage),
                                                                                          Lineage))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ace_groupInfo <- split(pmen3_ace_data$id, pmen3_ace_data$Resistance)
lineages <- dplyr::count(pmen3_ace_data, Resistance)
colours_tree <- gg_color_hue(nrow(lineages) + 1)
names(colours_tree) <- c(lineages$Resistance,"0")

pmen3_tree <- groupOTU(pmen3_tree, ace_groupInfo)
pmen3_tree_plot <- ggtree(pmen3_tree, aes(colour = group), layout = "rectangular") +
  scale_color_manual(values = colours_tree, name = "Penicillin resistance")


###############################################################################
## Simply tip labels only #####################################################
###############################################################################


p <- ggtree(ape::unroot(pmen3_tree))
pmen3_tree_dat <- pmen3_ace_data[,c("id","Resistance")]
p <- p %<+% pmen3_tree_dat + geom_tippoint(aes(color=Resistance))
p

p <- ggtree(pmen9_tree, layout = "circular")
pmen9_tree_dat <- pmen9_data[,c("id","Resistance")]
p <- p %<+% pmen9_tree_dat + geom_tippoint(aes(color=Resistance))
p



pmen3_tree_plot
png(filename = "~/Dropbox/phd/LSA/acba_tree.png", width = 21, height = 17,units = "cm", res = 1000)
print(ace_tree)
dev.off()


library(tibble)
library(tidyr)
library(Biostrings)
library(treeio)
library(ggplot2)
library(ggtree)

tree <- read.tree("data/HPV58.tree")

clade <- c(A3 = 92, A1 = 94, A2 = 108, B1 = 156, B2 = 159, C = 163, D1 = 173, D2 = 176)
tree <- groupClade(tree, clade)
cols <- c(A1 = "#EC762F", A2 = "#CA6629", A3 = "#894418", B1 = "#0923FA", 
          B2 = "#020D87", C = "#000000", D1 = "#9ACD32",D2 = "#08630A")

## visualize the tree with tip labels and tree scale
p <- ggtree(tree, aes(color = group), ladderize = FALSE) %>% rotate(rootnode(tree)) + 
  geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2.5) +
  geom_treescale(x = 0, y = 1, width = 0.002) + 
  scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Lineage",
                     breaks = c("A1", "A2", "A3", "B1", "B2", "C", "D1", "D2")) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme_tree2(legend.position = c(.1, .88))
## Optional
## add labels for monophyletic (A, C and D) and paraphyletic (B) groups  
p <- p + geom_cladelabel(94, "italic(A1)", color = cols[["A1"]], offset = .003, align = TRUE, 
                         offset.text = -.001, barsize = 1.2, extend = c(0, 0.5), parse = TRUE) +
  geom_cladelabel(108, "italic(A2)", color = cols[["A2"]], offset = .003, align = TRUE, 
                  offset.text = -.001, barsize = 1.2, extend = 0.5, parse = TRUE) +
  geom_cladelabel(131, "italic(A3)", color = cols[["A3"]], offset = .003, align = TRUE, 
                  offset.text = -.001, barsize = 1.2, extend = c(0.5, 0), parse = TRUE) +
  geom_cladelabel(92, "italic(A)", color = "darkgrey", offset = .00315, align = TRUE, 
                  offset.text = 0.0002, barsize = 2, fontsize = 5, parse = TRUE) +
  geom_cladelabel(156, "italic(B1)", color = cols[["B1"]], offset = .003, align = TRUE, 
                  offset.text = -.001, barsize = 1.2, extend = c(0, 0.5), parse = TRUE) +
  geom_cladelabel(159, "italic(B2)", color = cols[["B2"]], offset = .003, align = TRUE, 
                  offset.text = -.001, barsize = 1.2, extend = c(0.5, 0), parse = TRUE) +
  geom_strip(65, 71, "italic(B)", color = "darkgrey", offset = 0.00315, align = TRUE, 
             offset.text = 0.0002, barsize = 2, fontsize = 5, parse = TRUE) + 
  geom_cladelabel(163, "italic(C)", color = "darkgrey", offset = .0031, align = TRUE, 
                  offset.text = 0.0002, barsize = 3.2, fontsize = 5, parse = TRUE) +
  geom_cladelabel(173, "italic(D1)", color = cols[["D1"]], offset = .003, align = TRUE, 
                  offset.text = -.001, barsize = 1.2, extend = c(0, 0.5), parse = TRUE) +
  geom_cladelabel(176, "italic(D2)", color = cols[["D2"]], offset = .003, align = TRUE, 
                  offset.text = -.001, barsize = 1.2, extend = c(0.5, 0), parse = TRUE) +
  geom_cladelabel(172, "italic(D)", color = "darkgrey", offset = .00315, align = TRUE, 
                  offset.text = 0.0002, barsize = 2, fontsize = 5, parse = TRUE) 
## Optional
## display support values
p <- p + geom_nodelab(aes(subset = (node == 92), label = "*"), 
                      color = "black", nudge_x = -.001, nudge_y = 1) +
  geom_nodelab(aes(subset = (node == 155), label = "*"), 
               color = "black", nudge_x = -.0003, nudge_y = -1) +
  geom_nodelab(aes(subset = (node == 158), label = "95/92/1.00"), 
               color = "black", nudge_x = -0.0001, nudge_y = -1, hjust = 1) +
  geom_nodelab(aes(subset = (node == 162), label = "98/97/1.00"), 
               color = "black", nudge_x = -0.0001, nudge_y = -1, hjust = 1) +
  geom_nodelab(aes(subset = (node == 172), label = "*"), 
               color = "black", nudge_x = -.0003, nudge_y = -1) 


################################################################################
## Can't get the ggtree plots to look quite right, but I don't like the colour #
## scale for the microreact tree, so will add in a Resistance colour category ##
################################################################################

pmen3_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_26_11_2020.csv",
                       stringsAsFactors = FALSE) %>% rename(Resistance = cdc_only_RF_SIR__autocolour)
pmen9_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_14_12_2020.csv",
                       stringsAsFactors = FALSE) %>% rename(Resistance = cdc_only_RF_SIR__autocolour)


pmen3_cols <- gg_color_hue(4)
pmen3_resdf <- cbind.data.frame(count(pmen3_data, Resistance)[,1], pmen3_cols)
colnames(pmen3_resdf) <- c("Resistance", "Resistance__colour")


pmen3_data <- pmen3_data %>% left_join(pmen3_resdf, by = c("Resistance" = "Resistance"))


pmen9_cols <- gg_color_hue(4)
pmen9_resdf <- cbind.data.frame(count(pmen9_data, Resistance)[,1], pmen9_cols[c(1,3:4)])
colnames(pmen9_resdf) <- c("Resistance", "Resistance__colour")

pmen9_data <- pmen9_data %>% left_join(pmen9_resdf, by = c("Resistance" = "Resistance"))

write.csv(pmen3_data,"~/Dropbox/phd/elife_paper/figures_data/pmen3_microreact_pen_res.csv",
          row.names = FALSE)
write.csv(pmen9_data, file = "~/Dropbox/phd/elife_paper/figures_data/pmen9_microreact_pen_res.csv", 
          row.names = FALSE)



pmen3_resdf$y_val <- 1

ggplot(pmen3_resdf, aes(x = Resistance, y = y_val, fill = Resistance, colour = Resistance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pmen3_resdf$Resistance__colour[c(4,1,3,2)], breaks = c("S", "I", "R", "ND")) + 
  scale_color_manual(values = pmen3_resdf$Resistance__colour[c(4,1,3,2)], breaks = c("S", "I", "R", "ND")) 
