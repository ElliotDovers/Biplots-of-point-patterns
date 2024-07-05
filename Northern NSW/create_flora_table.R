# load the species information
load("flora_species.RDATA")

# separate the possibly duplicate species with an A and B
sp.info$full_name[sp.info$spid %in% paste0("nsw", c(39, 40))] <- paste(sp.info$full_name[sp.info$spid %in% paste0("nsw", c(39, 40))], c("A", "B"))
# create nicer names for the species groups
sp.info$type <- factor(sp.info$group, levels = c("ot", "ou", "rt", "ru"),
                       labels = c("Open-forest Tree", "Open-forest Understorey", "Rainforest Tree", "Rainforest Understorey"))

library(xtable)
print(xtable(sp.info[,c("full_name", "type", "n_po")], caption = "Flora species included in the northern NSW dataset.", label = "tab:supp_nsw_flora", align = c("l", "l", "l", "r")))
