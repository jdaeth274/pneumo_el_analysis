###############################################################################
## Serotypes PMEN3

pmen3_changes <- read.csv("~/Dropbox/phd/elife_paper/data/pmen3_st_changes.csv",
                          stringsAsFactors = FALSE)
count(pmen3_changes, insertion_node)
