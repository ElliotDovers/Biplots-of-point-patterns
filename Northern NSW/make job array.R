# need to determine some appropriate numbers of basis functions to look at

library(glmmTMB)

# get the data
source("get_and_format_flora_data.R")
rm(dat_pa)

# load the species information
load("flora_species.RDATA")

# get the groups to which the flora species belongs to subset the data
flora_groups <- unique(sp.info$group)

# subset the presence-only data (in the case we are using restricted species)
pres <- pres[pres$group %in% flora_groups, ]

# want to use species with at least 10 presences
pres <- pres[pres$spid %in% names(table(pres$spid))[table(pres$spid) >= 10], ]


# this is an exhaustive list of realised k (which goes up to k = 681 post-pruning)
k <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 300, 325, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 925, 975, 1050, 1100, 1175, 1250, 1325, 1400, 1475, 1550, 1625, 1700, 1800, 1875, 1975, 2050)
radii <- NULL
actual_ks <- NULL

for (i in k) {
  bfs <- make_basis(i, domain, from.package = "scampr", coord.names = c("xm", "ym"))
  tmp <- bf_matrix(bfs, rbind(pres[,c("xm", "ym")], quad[,c("xm", "ym")]))
  bf.check <- apply(tmp, 2, function(x){sum(x > 0)})
  bfs <- bfs[bf.check != 0, ]
  actual_ks <- c(actual_ks, nrow(bfs))
  radii <- c(radii, bfs$scale[1])
}

d <- c(2, 4, 6)

# expand out all combinations
tab <- data.frame(expand.grid(k, d))
colnames(tab) <- c("k", "d")
tab$job <- 1:nrow(tab)
tab$projected_k <- actual_ks[match(tab$k, k)]

save(tab, file = "job_array.RDATA")
