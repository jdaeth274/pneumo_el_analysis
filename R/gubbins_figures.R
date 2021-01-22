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


pmen3_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_19_01_2021.csv",
                       stringsAsFactors = FALSE) %>% mutate(Country__autocolour = ifelse(Country__autocolour == "S. Africa","South Africa",Country__autocolour))
pmen9_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_19_01_2021.csv",
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
pmen3_19a <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN3_19A_clade.csv",
                      stringsAsFactors = FALSE)

clade_colors <- gg_color_hue(3)

countries <- nrow(count(pmen3_data, Country__autocolour))
country_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(countries)

country_colours <- cbind.data.frame(count(pmen3_data, Country__autocolour)[,1], country_cols)
colnames(country_colours) <- c("country", "colour")
country_colours <- country_colours %>%
  mutate(country = ifelse(grepl("RUSSIAN",country),"Russia",ifelse(grepl("INDIA", country),"India",country)))


serotypes <- pmen3_data %>% rename(serotype = serotype__autocolour) %>%
  mutate(serotype = ifelse(serotype == "coverage too low", "NT", serotype)) %>% count(serotype)
sero_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(serotypes))

sero_colours <- cbind.data.frame(serotypes[,1], sero_cols)
colnames(sero_colours) <- c("serotype", "colour")



pmen3_phandango <- pmen3_data[,c("id","In.Silico.St__autocolour","Country__autocolour", "serotype__autocolour")] %>% 
  rename(ST = In.Silico.St__autocolour) %>% mutate(Clade = ifelse(id %in% pmen3_st156$id, "ST156","ST162")) %>%
  mutate(Clade = ifelse(id %in% pmen3_st143$id, "ST143",Clade)) %>% 
  mutate(id = paste(id, "57195_E01.1", sep = "_.")) %>%
  mutate("Clade:colour" = ifelse(Clade == "ST156",clade_colors[1], ifelse(Clade == "ST143", clade_colors[2], clade_colors[3]))) %>%
  rename(Country = Country__autocolour) %>% left_join(y = country_colours, by = c("Country" = "country")) %>%
  rename("Country:colour" = colour) %>% rename(Serotype = serotype__autocolour) %>%
  mutate(Serotype = ifelse(Serotype == "coverage too low", "NT", Serotype)) %>% 
  left_join( sero_colours, by = c("Serotype" = "serotype")) %>% rename("Serotype:colour" = colour) %>%
  mutate(Country = ifelse(grepl("RUSSIAN",Country),"Russia",ifelse(grepl("INDIA", Country),"India",Country)))

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
  geom_bar(stat = "identity") + scale_fill_manual(values = clade_colours) +
  scale_color_manual(values = clade_colours)

  ggplot(data = pmen3_phandango, aes(x = Serotype, y = y_val , colour = Serotype, fill = Serotype)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = sero_colours$colour) +
  scale_color_manual(values = sero_colours$colour) + labs(fill = "Serotype", colour = "Serotype")


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
pmen9_china <- read.csv("~/Dropbox/phd/elife_paper/figures_data/pmen9_china_clade.csv",
                        stringsAsFactors = FALSE)

clade_colors <- gg_color_hue(4)

countries <- nrow(count(pmen9_data, country__autocolour))
country_cols <- colorRampPalette(brewer.pal(12, "Paired"))(countries)

country_colours <- cbind.data.frame(count(pmen9_data, country__autocolour)[,1], country_cols)
colnames(country_colours) <- c("country", "colour") 
country_colours <- country_colours %>%
  mutate(country = ifelse(grepl("RUSSIAN",country),"Russia",ifelse(grepl("INDIA", country),"India",country)))

serotypes <- pmen9_data %>% rename(serotype = serotype__autocolour) %>%
  mutate(serotype = ifelse(serotype == "coverage too low", "NT", serotype)) %>% count(serotype)
sero_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(serotypes))

sero_colours <- cbind.data.frame(serotypes[,1], sero_cols)
colnames(sero_colours) <- c("serotype", "colour")



pmen9_phandango <- pmen9_data[,c("id","country__autocolour", "serotype__autocolour")] %>% 
  rename(Country = country__autocolour) %>% mutate(Clade = ifelse(id %in% pmen9_german$id, "German","Other")) %>%
  mutate(Clade = ifelse(id %in% pmen9_USA$id, "USA",Clade)) %>% mutate(Clade = ifelse(id %in% pmen9_SA$id, "South African",Clade)) %>%
  mutate("Clade:colour" = ifelse(Clade == "German",clade_colors[1],
                                 ifelse(Clade == "USA", clade_colors[2],
                                        ifelse(Clade == "South African",
                                               clade_colors[3],clade_colors[4])))) %>% 
  left_join(y = country_colours, by = c("Country" = "country")) %>%
  rename("Country:colour" = colour) %>% rename(Serotype = serotype__autocolour) %>%
  mutate(Serotype = ifelse(Serotype == "coverage too low", "NT", Serotype)) %>% 
  left_join( sero_colours, by = c("Serotype" = "serotype")) %>% rename("Serotype:colour" = colour) %>%
  mutate(Country = ifelse(grepl("RUSSIAN",Country),"Russia",ifelse(grepl("INDIA", Country),"India",Country)))

write.csv(pmen9_phandango, file = "~/Dropbox/phd/elife_paper/figures_data/PMEN9_gubbins_clade_only.csv",
          row.names = FALSE)


## PLot out the legends for these 

pmen9_phandango$x_val <- seq(1, nrow(pmen9_phandango))
pmen9_phandango$y_val <- seq(1, nrow(pmen9_phandango))
pmen9_phandango$density <- rnorm(nrow(pmen9_phandango), mean = 10)

plot_cols <- country_colours[-c(1,12),]
ggplot(data = pmen9_phandango, aes(x = Country, y = y_val , colour = Country, fill = Country)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = country_colours$colour, breaks = plot_cols$country) +
  scale_color_manual(values = country_colours$colour, breaks = plot_cols$country)

ggplot(data = pmen9_phandango, aes(x = Clade, y = y_val , colour = Clade, fill = Clade)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = clade_colors,breaks = c("German","South African","USA","Other")) +
  scale_color_manual(values = clade_colors, breaks = c("German","South African","USA","Other"))

ggplot(data = pmen9_phandango, aes(x = Serotype, y = y_val , colour = Serotype, fill = Serotype)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = sero_colours$colour) +
  scale_color_manual(values = sero_colours$colour)




###############################################################################
## Make the Same figures for penicillin resistance, trimethoprim resistance ###
## and MGE presence absence ###################################################
###############################################################################

pmen3_data <- read.csv("~/Dropbox/phd/PMEN3/R_DATA_analysis/data/epi_data/epi_data_19_01_2021.csv",
                       stringsAsFactors = FALSE) %>% mutate(Country__autocolour = ifelse(Country__autocolour == "S. Africa","South Africa",Country__autocolour)) %>%
  mutate(Resistance = cdc_only_RF_SIR__autocolour)
pmen9_data <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/epi_data/epi_data_19_01_2021.csv",
                       stringsAsFactors = FALSE) %>%  mutate(Resistance = cdc_only_RF_SIR__autocolour) %>% 
  mutate(country__autocolour = ifelse(country__autocolour == "BRAZIL","Brazil",country__autocolour))

## Set up clade colour 

clade_colours <- gg_color_hue(4)

clades <- c("MLST156","MLST162","MLST143","19A")
clade_colours_df <- cbind.data.frame(clades, clade_colours)
colnames(clade_colours_df) <- c("Clade", "colour") 



## Set up country colours

countries <- nrow(count(pmen3_data, Country__autocolour))
country_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(countries)

country_colours <- cbind.data.frame(count(pmen3_data, Country__autocolour)[,1], country_cols)
colnames(country_colours) <- c("country", "colour")
country_colours <- country_colours %>%
  mutate(country = ifelse(grepl("RUSSIAN",country),"Russia",ifelse(grepl("INDIA", country),"India",country)))

## Set up Serotypes colours
serotypes <- pmen3_data %>% rename(serotype = serotype__autocolour) %>%
  mutate(serotype = ifelse(serotype == "coverage too low", "NT", serotype)) %>% count(serotype)
sero_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(serotypes))

sero_colours <- cbind.data.frame(serotypes[,1], sero_cols)
colnames(sero_colours) <- c("serotype", "colour")

## Set up penicillin resistance colours 

pen_cols <- toupper(c("#ffeda0","#808080","#feb24c","#f03b20"))[c(3,2,4,1)]#[c(4,2,1,3)]
#pen_cols <- sub("FF","",pen_cols)
pmen3_resdf <- cbind.data.frame(count(pmen3_data, Resistance)[,1], pen_cols)
colnames(pmen3_resdf) <- c("Penicillin", "Penicillin:colour")

trimethoprim_cols <- pen_cols[2:4]
pmen3_trim_df <- cbind.data.frame(count(pmen3_data, trimethoprim__autocolour)[,1], trimethoprim_cols)
colnames(pmen3_trim_df) <- c("Trimethoprim", "Trimethoprim:colour")

sulfa_cols <- pen_cols[2:4]
pmen3_sulfa_df <- cbind.data.frame(count(pmen3_data, sulfa__autocolour)[,1], sulfa_cols)
colnames(pmen3_sulfa_df) <- c("Sulfamethoxazole", "Sulfamethoxazole:colour")

mge_presence_cols <- toupper(c("#91bfdb","#fc8d59"))

pmen3_phandango <- pmen3_data[,c("id","In.Silico.St__autocolour","Country__autocolour", "serotype__autocolour", "cdc_only_RF_SIR",
                                 "trimethoprim__autocolour","sulfa__autocolour", "mega_isa__autocolour", "tn916_isa__autocolour")] %>% 
  rename(ST = In.Silico.St__autocolour) %>% mutate(Clade = ifelse(id %in% pmen3_st156$id, "MLST156","MLST162")) %>%
  mutate(Clade = ifelse(id %in% pmen3_st143$id, "MLST143",Clade)) %>% 
  mutate(Clade = ifelse(id %in% pmen3_19a$id, "19A", Clade)) %>%
  mutate(id = paste(id, "57195_E01.1", sep = "_.")) %>%
  left_join(y = clade_colours_df, by = c("Clade" = "Clade")) %>% rename("Clade:colour" = colour) %>%
  #mutate("Clade:colour" = ifelse(Clade == "MLST156",clade_colors[1], ifelse(Clade == "ST143", clade_colors[2], clade_colors[3]))) %>%
  rename(Country = Country__autocolour) %>% left_join(y = country_colours, by = c("Country" = "country")) %>%
  rename("Country:colour" = colour) %>% rename(Serotype = serotype__autocolour) %>%
  mutate(Serotype = ifelse(Serotype == "coverage too low", "NT", Serotype)) %>% 
  left_join( sero_colours, by = c("Serotype" = "serotype")) %>% rename("Serotype:colour" = colour) %>%
  mutate(Country = ifelse(grepl("RUSSIAN",Country),"Russia",ifelse(grepl("INDIA", Country),"India",Country))) %>%
  rename(Penicillin = cdc_only_RF_SIR) %>%  left_join(pmen3_resdf, by = c("Penicillin" = "Penicillin")) %>%
  rename(Trimethoprim = trimethoprim__autocolour) %>% left_join(pmen3_trim_df, by = c("Trimethoprim" = "Trimethoprim")) %>%
  rename(Sulfamethoxazole = sulfa__autocolour) %>% left_join(pmen3_sulfa_df, by = c("Sulfamethoxazole" = "Sulfamethoxazole")) %>%
  rename("Tn1207.1" = mega_isa__autocolour) %>% rename("Tn916" = tn916_isa__autocolour) %>%
  mutate("Tn1207.1:colour" = ifelse(Tn1207.1 == "Yes", mge_presence_cols[2], mge_presence_cols[1])) %>%
  mutate("Tn916:colour" = ifelse(Tn916 == "Yes", mge_presence_cols[2], mge_presence_cols[1])) 


colnames(pmen3_phandango)
pmen3_phandango <- pmen3_phandango[,c(1,10,2:4,5,7,6,8,9,11:18)]
write.csv(pmen3_phandango, file = "~/Dropbox/phd/elife_paper/figures_data/pmen3_microreact_pen_res.csv", 
          row.names = FALSE)

## PMEN9 cols now 


pmen9_german <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN9_german_clade.csv",
                         stringsAsFactors = FALSE)
pmen9_USA <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN9_USA_clade.csv",
                      stringsAsFactors = FALSE)
pmen9_SA <- read.csv("~/Dropbox/phd/elife_paper/figures_data/PMEN9_SA_clade.csv",
                     stringsAsFactors = FALSE)


clade_colours <- gg_color_hue(4)

clades <- c("Chinese", "German","South African","USA", "Other")
clade_colours_df <- cbind.data.frame(clades, c(clade_colours, "#000000"))
colnames(clade_colours_df) <- c("Clade", "colour") 

countries <- nrow(count(pmen9_data, country__autocolour))
country_cols <- colorRampPalette(brewer.pal(12, "Paired"))(countries)

country_colours <- cbind.data.frame(count(pmen9_data, country__autocolour)[,1], country_cols)
colnames(country_colours) <- c("country", "colour") 
country_colours <- country_colours %>%
  mutate(country = ifelse(grepl("RUSSIAN",country),"Russia",ifelse(grepl("INDIA", country),"India",country)))

serotypes <- pmen9_data %>% rename(serotype = serotype__autocolour) %>%
  mutate(serotype = ifelse(serotype == "coverage too low", "NT", serotype)) %>% count(serotype)
sero_cols <- mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nrow(serotypes))

sero_colours <- cbind.data.frame(serotypes[,1], sero_cols)
colnames(sero_colours) <- c("serotype", "colour")


pmen9_cols <- pen_cols
pmen9_resdf <- cbind.data.frame(count(pmen9_data, Resistance)[,1], pmen9_cols[c(1,3:4)])
colnames(pmen9_resdf) <- c("Penicillin", "Penicillin:colour")

trimethoprim_cols <- pen_cols[3:4]
pmen9_trim_df <- cbind.data.frame(count(pmen9_data, trimethoprim__autocolour)[,1], trimethoprim_cols)
colnames(pmen9_trim_df) <- c("Trimethoprim", "Trimethoprim:colour")

sulfa_cols <- pen_cols[3:4]
pmen9_sulfa_df <- cbind.data.frame(count(pmen9_data, sulfa__autocolour)[,1], sulfa_cols)
colnames(pmen9_sulfa_df) <- c("Sulfamethoxazole", "Sulfamethoxazole:colour")


pmen9_phandango <- pmen9_data[,c("id","country__autocolour", "serotype__autocolour", "tn916_isa__autocolour",
                                 "mega_isa__autocolour", "cdc_only_RF_SIR__autocolour","sulfa__autocolour",
                                 "trimethoprim__autocolour")] %>% 
  rename(Country = country__autocolour) %>% mutate(Clade = ifelse(id %in% pmen9_german$id, "German","Other")) %>%
  mutate(Clade = ifelse(id %in% pmen9_USA$id, "USA",Clade)) %>% mutate(Clade = ifelse(id %in% pmen9_SA$id, "South African",Clade)) %>%
  mutate(Clade = ifelse(id %in% pmen9_china$id, "Chinese",Clade)) %>%
  left_join(clade_colours_df, by = c("Clade" = "Clade")) %>% rename("Clade:colour" = colour) %>%
  #mutate("Clade:colour" = ifelse(Clade == "German",clade_colors[1],
  #                               ifelse(Clade == "USA", clade_colors[2],
  #                                      ifelse(Clade == "South African",
   #                                            clade_colors[3],clade_colors[4])))) %>% 
  left_join(y = country_colours, by = c("Country" = "country")) %>%
  rename("Country:colour" = colour) %>% rename(Serotype = serotype__autocolour) %>%
  mutate(Serotype = ifelse(Serotype == "coverage too low", "NT", Serotype)) %>% 
  left_join( sero_colours, by = c("Serotype" = "serotype")) %>% rename("Serotype:colour" = colour) %>%
  mutate(Country = ifelse(grepl("RUSSIAN",Country),"Russia",ifelse(grepl("INDIA", Country),"India",Country))) %>%
  rename(Penicillin = cdc_only_RF_SIR__autocolour) %>%  left_join(pmen9_resdf, by = c("Penicillin" = "Penicillin")) %>%
  rename(Trimethoprim = trimethoprim__autocolour) %>% left_join(pmen9_trim_df, by = c("Trimethoprim" = "Trimethoprim")) %>%
  rename(Sulfamethoxazole = sulfa__autocolour) %>% left_join(pmen9_sulfa_df, by = c("Sulfamethoxazole" = "Sulfamethoxazole")) %>%
  rename("Tn1207.1" = mega_isa__autocolour) %>% rename("Tn916" = tn916_isa__autocolour) %>%
  mutate("Tn1207.1:colour" = ifelse(Tn1207.1 == "Yes", mge_presence_cols[2], mge_presence_cols[1])) %>%
  mutate("Tn916:colour" = ifelse(Tn916 == "Yes", mge_presence_cols[2], mge_presence_cols[1])) 



pmen9_phandango <- pmen9_phandango[,c(1,2:3,9,6,7,8,5,4,10:17)]

write.csv(pmen9_phandango, file = "~/Dropbox/phd/elife_paper/figures_data/pmen9_microreact_pen_res.csv",
          row.names = FALSE)

## plot resistance and mge presence legends 

pmen3_phandango$x_val <- 1
pmen3_phandango$y_val <- 2
pmen3_phandango$MGE <- pmen3_phandango$`Tn916:colour`
pmen9_phandango$x_val <- 1
pmen9_phandango$y_val <- 2

ggplot(pmen3_phandango, aes(x = x_val, y = y_val, color = `Tn916:colour`, fill = `Tn916:colour`)) + 
  geom_bar(stat = "identity") + scale_fill_manual(labels = c("Absent", "Present"), values = mge_presence_cols) + 
  scale_color_manual(labels = c("Absent", "Present"), values = mge_presence_cols) + 
  labs(colour = "MGE presence", fill = "MGE presence")

ggplot(pmen3_phandango, aes(x = x_val, y = y_val, color = Penicillin, fill = Penicillin)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = pen_cols[c(4,1,3,2)] ,breaks = c("S","I","R","ND")) + 
  scale_color_manual(values = pen_cols[c(4,1,3,2)], breaks = c("S","I","R","ND")) + 
  labs(colour = "Resistance", fill = "Resistance")
ggplot(pmen9_phandango, aes(x = x_val, y = y_val, color = Clade, fill = Clade)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = clade_colours_df$colour[1:5], breaks = clade_colours_df$Clade[1:5]) +
  scale_color_manual(values = clade_colours_df$colour[1:5], breaks = clade_colours_df$Clade[1:5])

## co-trim stats

count(pmen9_phandango, Trimethoprim)
count(pmen9_phandango, Sulfamethoxazole)

SA_phandango <- pmen9_phandango[pmen9_phandango$Clade == "South African",]
count(SA_phandango, Trimethoprim)
count(SA_phandango, Sulfamethoxazole)

SA_country_phand <- pmen9_phandango[pmen9_phandango$Country == "South Africa",]


## pmen3

ST143 <- pmen3_phandango[pmen3_phandango$Clade == "ST143",]
ST_162 <- pmen3_phandango[pmen3_phandango$Clade == "ST162",]

count(ST143, Trimethoprim)
count(ST143, Sulfamethoxazole)

SA_pmen3 <- pmen3_phandango[pmen3_phandango$Country == "South Africa",]
count(SA_pmen3, Trimethoprim)
count(SA_pmen3, Sulfamethoxazole)

