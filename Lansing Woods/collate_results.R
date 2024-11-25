## Collate results of the convergence of the quadrature approximation ##########

home.wd <- getwd()

# initialise the result storage
dat <- NULL

res.list <- list()
res.objs <- list.files(paste0(home.wd, "/results"))[grepl("res_", list.files(paste0(home.wd, "/results")), fixed = T)]
for (job in 1:length(res.objs)) {
  load(paste0(home.wd, "/results/", res.objs[job]))
  res.list[[job]] <- res
}
dat <- do.call(rbind, res.list)
rm(res.list)
# re-order
dat <- dat[order(dat$q, dat$rep), ]

save(dat, file = "quad.conv.RDATA")
