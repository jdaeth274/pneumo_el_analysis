###############################################################################
## Going back over the reinert ecdc data merging ##############################
###############################################################################

reinert_data <- read.csv("~/Dropbox/phd/elife_paper/data/Reinert_paper_digitezed.csv",
                         stringsAsFactors = FALSE)
ecdc_2008 <- read.csv("~/Dropbox/phd/elife_paper/data/german_1997_2008_macrolide.csv",
                      stringsAsFactors = FALSE)
ecdc_europe_total <- read.csv("~/Dropbox/phd/elife_paper/data/ecdc_europe_wide_data_97_2016.csv",
                              stringsAsFactors = FALSE)
ecdc_combined <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/german_clade_analysis/germany_ecdc_reinert_australian_pen_and_mac_usage.csv",
                          stringsAsFactors = FALSE)
updated_ecdc <- read.csv("~/Dropbox/phd/PMEN9/PMEN9_R_project/data/german_clade_analysis/updated_ecdc_pen.csv",
                         stringsAsFactors = FALSE)
updated_german_ecdc <- updated_ecdc[,c(1:6)]
colnames(updated_german_ecdc) <- sub("Germany","update",colnames(updated_german_ecdc))

ecdc_combined[5:20,]
ecdc_2008
head(ecdc_europe_total)

colnames(ecdc_2008)
colnames(ecdc_europe_total)
colnames(reinert_data)
colnames(ecdc_combined)

ecdc_2008_short <- ecdc_2008 %>% select(Year, Macrolide_usage.DDD.1000.) %>% rename(short_08 = Macrolide_usage.DDD.1000.)
ecdc_europe_total_short <- ecdc_europe_total %>% select(Year, Germany_mac) %>% rename(short_eur = Germany_mac)
reinert_data_short <- reinert_data %>% select(Year, All_macrolide_use) %>% rename(short_rein = All_macrolide_use)

compo_short <- ecdc_combined %>% select(Year, Germany_mac, Germany_pen) %>% left_join(ecdc_2008_short) %>% left_join(ecdc_europe_total_short) %>%
  left_join(reinert_data_short) %>% mutate(scale = short_rein / Germany_mac) %>% mutate(scale_08 = short_rein / short_08) 

just_macs_ecdc <- c(compo_short$Germany_mac[1:5],compo_short$short_08[6:18])

just_macs <- cbind.data.frame(compo_short$Year[1:18], just_macs_ecdc)
colnames(just_macs) <- c("Year","macrolides")

write.csv(just_macs, file = "~/Dropbox/phd/PMEN9/PMEN9_R_project/data/german_clade_analysis/german_mac_only.csv",
          row.names = FALSE)

## actual scaling factor now for macs 

head(compo_short)
new_scale <- mean(compo_short$scale[6:9])
compo_short$short_rein[1:5] / new_scale

###############################################################################
## Switch back over to the updated pen and Mac values #########################
###############################################################################

updated_ecdc <- ecdc_combined %>% left_join(updated_german_ecdc)

updated_ecdc$update_mac[1:5] <- compo_short$short_rein[1:5] / new_scale

## Now for the pencillin, so we have 4.1 and 4.4 for 92 and 94 from the australian paper
## now we'll impute a linear trend from 94 to 97 for the rest of the values 

updated_ecdc$update_tot[1:3] <- updated_ecdc$Germany_pen[1:3]
pen_95_96 <- round(seq(4.4,4.39, length.out = 4)[2:3], digits = 2)

updated_ecdc$update_tot[4:5] <- pen_95_96
updated_ecdc$update_ratio[1:5] <- updated_ecdc$update_mac[1:5] / updated_ecdc$update_tot[1:5]

write.csv(updated_ecdc, file = "~/Dropbox/phd/PMEN9/PMEN9_R_project/data/german_clade_analysis/updated_combined_ecdc.csv",
          row.names = FALSE)

