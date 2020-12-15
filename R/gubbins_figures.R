###############################################################################
## Tweaking the phanango plots to include the ST and the element presence/ ####
## absence ####################################################################
###############################################################################

require(dplyr)
require(stringr)
require(RColorBrewer)
require(ggplot2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


pmen3_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_26_11_2020.csv",
                       stringsAsFactors = FALSE)
pmen9_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_26_11_2020.csv",
                       stringsAsFactors = FALSE)


## For the moment lets just make the ST present/ identify key clades for the 
## PMEN9 lineage. We'll see if we want to add in further data on resistance and 
## element presence.

## Need to append _57195_E01.1 to the tree ids 

pmen3_st156 <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN3_st156_clade.csv",
                        stringsAsFactors = FALSE) %>% select(c("id","In.Silico.St__autocolour")) %>%
  rename(ST = In.Silico.St__autocolour)
pmen3_st143 <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN3_st143_clade.csv",
                        stringsAsFactors = FALSE) %>% select(c("id","In.Silico.St__autocolour")) %>%
  rename(ST = In.Silico.St__autocolour)

clade_colors <- gg_color_hue(3)

countries <- nrow(count(pmen3_data, Country__autocolour))
country_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(countries)

country_colours <- cbind.data.frame(count(pmen3_data, Country__autocolour)[,1], country_cols)
colnames(country_colours) <- c("country", "colour")

pmen3_phandango <- pmen3_data[,c("id","In.Silico.St__autocolour","Country__autocolour")] %>% 
  rename(ST = In.Silico.St__autocolour) %>% mutate(Clade = ifelse(id %in% pmen3_st156$id, "ST156","ST162")) %>%
  mutate(Clade = ifelse(id %in% pmen3_st143$id, "ST143",Clade)) %>% 
  mutate(id = paste(id, "57195_E01.1", sep = "_.")) %>%
  mutate("Clade:colour" = ifelse(Clade == "ST156",clade_colors[1], ifelse(Clade == "ST143", clade_colors[2], clade_colors[3]))) %>%
  rename(Country = Country__autocolour) %>% left_join(y = country_colours, by = c("Country" = "country")) %>%
  rename("Country:colour" = colour)

write.csv(pmen3_phandango, file = "~/Dropbox/phd/elife_paper/figures_data/PMEN3_gubbins_clade_only.csv",
          row.names = FALSE)


## PLot out the legends for these 

pmen3_phandango$x_val <- seq(1, nrow(pmen3_phandango))
pmen3_phandango$y_val <- seq(1, nrow(pmen3_phandango))
pmen3_phandango$density <- rnorm(nrow(pmen3_phandango), mean = 10)

ggplot(data = pmen3_phandango, aes(x = Country, y = y_val , colour = Country, fill = Country)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = country_colours$colour, breaks = country_colours$country[2:nrow(country_colours)]) +
  scale_color_manual(values = country_colours$colour, breaks = country_colours$country[2:nrow(country_colours)])

ggplot(data = pmen3_phandango, aes(x = Clade, y = y_val , colour = Clade, fill = Clade)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = clade_colors, breaks = country_colours$country[2:nrow(country_colours)]) +
  scale_color_manual(values = clade_colors)



###############################################################################
## PMEN9 gubbins now ##########################################################
###############################################################################

pmen9_data <- pmen9_data %>% 
  mutate(country__autocolour = ifelse(country__autocolour == "BRAZIL","Brazil",country__autocolour))

pmen9_german <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN9_german_clade.csv",
                         stringsAsFactors = FALSE)
pmen9_USA <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN9_USA_clade.csv",
                         stringsAsFactors = FALSE)
pmen9_SA <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN9_SA_clade.csv",
                         stringsAsFactors = FALSE)

clade_colors <- gg_color_hue(4)

countries <- nrow(count(pmen9_data, country__autocolour))
country_cols <- colorRampPalette(brewer.pal(12, "Paired"))(countries)

country_colours <- cbind.data.frame(count(pmen9_data, country__autocolour)[,1], country_cols)
colnames(country_colours) <- c("country", "colour")

pmen9_phandango <- pmen9_data[,c("id","country__autocolour")] %>% 
  rename(Country = country__autocolour) %>% mutate(Clade = ifelse(id %in% pmen9_german$id, "German","Other")) %>%
  mutate(Clade = ifelse(id %in% pmen9_USA$id, "USA",Clade)) %>% mutate(Clade = ifelse(id %in% pmen9_SA$id, "South African",Clade)) %>%
  mutate("Clade:colour" = ifelse(Clade == "German",clade_colors[1],
                                 ifelse(Clade == "USA", clade_colors[2],
                                        ifelse(Clade == "South African",
                                               clade_colors[3],clade_colors[4])))) %>% 
  left_join(y = country_colours, by = c("Country" = "country")) %>%
  rename("Country:colour" = colour)

write.csv(pmen9_phandango, file = "~/Dropbox/phd/elife_paper/figures_data/PMEN9_gubbins_clade_only.csv",
          row.names = FALSE)


## PLot out the legends for these 

pmen9_phandango$x_val <- seq(1, nrow(pmen9_phandango))
pmen9_phandango$y_val <- seq(1, nrow(pmen9_phandango))
pmen9_phandango$density <- rnorm(nrow(pmen9_phandango), mean = 10)

ggplot(data = pmen9_phandango, aes(x = Country, y = y_val , colour = Country, fill = Country)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = country_colours$colour, breaks = country_colours$country[2:nrow(country_colours)]) +
  scale_color_manual(values = country_colours$colour, breaks = country_colours$country[2:nrow(country_colours)])

ggplot(data = pmen9_phandango, aes(x = Clade, y = y_val , colour = Clade, fill = Clade)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = clade_colors) +
  scale_color_manual(values = clade_colors)






