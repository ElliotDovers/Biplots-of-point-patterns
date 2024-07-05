# make job array script

true_k <- seq(50, 650, by = 50)
# fit_k <- seq(50, 650, by = 50)
sim <- 1:100

# expand out all combinations
# tab <- data.frame(expand.grid(true_k, fit_k, sim))
tab <- data.frame(expand.grid(true_k, sim))
colnames(tab) <- c("true_k", "sim")
# colnames(tab) <- c("true_k", "fit_k", "sim")
tab$job <- 1:nrow(tab)

save(tab, file = "job_array.RDATA")
