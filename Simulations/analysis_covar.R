home.wd <- getwd()
library(dplyr)
library(ggplot2)
plot.res <- 500
# toggle to show a subset of the panels
sub.actual.k <- c(1,2,4,6,8,13)

# analyse the simulation results


# collate or load the covariate-included results
if (file.exists("collated_results_covar.RDATA")) {
  load("collated_results_covar.RDATA")
} else {
  # initialise the result storage
  dat1 <- NULL
  
  # change into appropriate result folder
  setwd(paste(home.wd, "results_covar", sep = "/"))
  # get a list of all the individual result files
  res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
  # inner loop through individual sim-by-model files
  for (job in res.list) {
    # load in the individual simulation results
    load(job)
    # add to the data
    dat1 <- rbind(dat1, res)
    rm(res)
  }
  
  # something stuffed up appending job information so correct here
  # dat1$true_k <- unlist(dat1$true_k)
  # dat1$sim <- unlist(dat1$sim)
  # dat1$job <- unlist(dat1$job)
  
  setwd(home.wd)
  save(list = "dat1", file = "collated_results_covar.RDATA")
}

# calculate the scenario averages for the known truth-based metrics
scen_summaries <- dat1 %>% group_by(true_k, set_k, fitted_k, actual_k) %>% summarise(mae_avg = mean(mae), mae_sd = sd(mae),
                                                                              rmse1 = sqrt(mean(se_beta1)),
                                                                              rmse2 = sqrt(mean(se_beta2)),
                                                                              rmse3 = sqrt(mean(se_beta3)),
                                                                              rmse4 = sqrt(mean(se_beta4)),
                                                                              rmse5 = sqrt(mean(se_beta5)),
                                                                              cp1 = mean(cover_beta1),
                                                                              cp2 = mean(cover_beta2),
                                                                              cp3 = mean(cover_beta3),
                                                                              cp4 = mean(cover_beta4),
                                                                              cp5 = mean(cover_beta5)
)

# calculate the mean and std. dev. for each metric path, for each simulation (so they can be standardised)
path_summaries <- dat1 %>% group_by(actual_k, sim) %>% summarise(ll_avg = mean(ll, na.rm = T), aic_avg = mean(aic, na.rm = T), bic_avg = mean(bic, na.rm = T),
                                                                 ll_sd = sd(ll, na.rm = T), aic_sd = sd(aic, na.rm = T), bic_sd = sd(bic, na.rm = T)
)

# loop through actual k scenarios and simulations to calculate the standardised metrics
dat1$ll_std <- NULL
dat1$aic_std <- NULL
dat1$bic_std <- NULL
for (k in sort(unique(dat1$actual_k))) {
  for (s in 1:max(dat1$sim)) {
    # scaled and centered or just centered
    # dat1$ll_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$ll[dat1$actual_k == k & dat1$sim == s] - path_summaries$ll_avg[path_summaries$actual_k == k & path_summaries$sim == s]) / path_summaries$ll_sd[path_summaries$actual_k == k & path_summaries$sim == s]
    dat1$ll_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$ll[dat1$actual_k == k & dat1$sim == s] - path_summaries$ll_avg[path_summaries$actual_k == k & path_summaries$sim == s])
    dat1$aic_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$aic[dat1$actual_k == k & dat1$sim == s] - path_summaries$aic_avg[path_summaries$actual_k == k & path_summaries$sim == s])
    dat1$bic_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$bic[dat1$actual_k == k & dat1$sim == s] - path_summaries$bic_avg[path_summaries$actual_k == k & path_summaries$sim == s])
  }
}

################################################################################
# Plot of MAE vs k and log-Likelihood ##########################################
################################################################################
png(filename = paste0(home.wd, "/figures/mae_vs_k_ll.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3,2)) # adjust according to sub.actual.k
for (i in sort(unique(dat1$actual_k))[sub.actual.k]) {
  # set the index of the loop control variable
  idx <- which(i == sort(unique(dat1$actual_k))[sub.actual.k])
  # determine margins based on the plot number
  if (idx == 1) {
    par(mar = c(2.5, 4.1, 1.5, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 2) {
    par(mar = c(2.5, 2.1, 1.5, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  } else if (idx == 3) {
    par(mar = c(3.25, 4.1, 0.75, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 4) {
    par(mar = c(3.25, 2.1, 0.75, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  } else if (idx == 5) {
    par(mar = c(4, 4.1, 0, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 6) {
    par(mar = c(4, 2.1, 0, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  }
  
  # subset the data
  tmp.dat <- dat1[dat1$actual_k == i, ]
  tmp.avg <- scen_summaries[scen_summaries$actual_k == i, ]
  
  # plot the selection metric
  with(tmp.dat, plot(fitted_k, ll_std, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = "", ylim = range(tmp.dat$ll_std, na.rm = T), xaxt = "n"))#if (idx %in% 5:6) {NULL} else {xaxt = "n"}))
  for (s in 1:max(dat1$sim)) {
    with(tmp.dat[tmp.dat$sim == s, ], lines(fitted_k, ll_std, col = ggplot2::alpha("grey10", 0.2)))
  }
  axis(4, ylim = range(tmp.dat$ll_std), col.axis = "black", col = "black")
  # only add label if the plot is on the right
  if (idx %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")}
  # add plot title
  mtext(text = c("A", "B", "C", "D", "E", "F")[idx], side = 3, cex = 1, line = 0.25, adj = 0)
  
  # plot the performance metric
  par(new = T)
  with(tmp.avg, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n", xaxt = "n"))
  if (idx %in% 5:6) {axis(1, ylim = range(tmp.avg$fitted_k))}
  axis(2, ylim = range(tmp.avg$mae_avg), col.axis = "royalblue")
  with(tmp.avg, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.avg, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  # only add label if the plot is on the left
  if (idx %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  # add x axis labels if the plots are those on the bottom
  if (idx == 5) {mtext(text = expression("#"~Basis~functions~"in"~each~of~the~(italic(r))~shared~and~(italic(m))~independent~fields), side = 1, srt = 90, cex = 1, xpd = T, outer = T, line = -1)}
  rm(tmp.dat, tmp.avg)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

################################################################################
# Plot of MAE vs k and AIC #####################################################
################################################################################
png(filename = paste0(home.wd, "/figures/mae_vs_k_aic.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3,2)) # adjust according to sub.actual.k
for (i in sort(unique(dat1$actual_k))[sub.actual.k]) {
  # set the index of the loop control variable
  idx <- which(i == sort(unique(dat1$actual_k))[sub.actual.k])
  # determine margins based on the plot number
  if (idx == 1) {
    par(mar = c(2.5, 4.1, 1.5, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 2) {
    par(mar = c(2.5, 2.1, 1.5, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  } else if (idx == 3) {
    par(mar = c(3.25, 4.1, 0.75, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 4) {
    par(mar = c(3.25, 2.1, 0.75, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  } else if (idx == 5) {
    par(mar = c(4, 4.1, 0, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 6) {
    par(mar = c(4, 2.1, 0, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  }
  
  # subset the data
  tmp.dat <- dat1[dat1$actual_k == i, ]
  tmp.avg <- scen_summaries[scen_summaries$actual_k == i, ]
  
  # plot the selection metric
  with(tmp.dat, plot(fitted_k, aic_std, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = "", ylim = range(tmp.dat$aic_std, na.rm = T), xaxt = "n"))#if (idx %in% 5:6) {NULL} else {xaxt = "n"}))
  for (s in 1:max(dat1$sim)) {
    with(tmp.dat[tmp.dat$sim == s, ], lines(fitted_k, aic_std, col = ggplot2::alpha("grey10", 0.2)))
  }
  axis(4, ylim = range(tmp.dat$aic_std), col.axis = "black", col = "black")
  # only add label if the plot is on the right
  if (idx %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Centered AIC", side = 4, line = 2.5, cex = 0.75, col = "black")}
  # add plot title
  mtext(text = c("A", "B", "C", "D", "E", "F")[idx], side = 3, cex = 1, line = 0.25, adj = 0)
  
  # plot the performance metric
  par(new = T)
  with(tmp.avg, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n", xaxt = "n"))
  if (idx %in% 5:6) {axis(1, ylim = range(tmp.avg$fitted_k))}
  axis(2, ylim = range(tmp.avg$mae_avg), col.axis = "royalblue")
  with(tmp.avg, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.avg, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  # only add label if the plot is on the left
  if (idx %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  # add x axis labels if the plots are those on the bottom
  if (idx == 5) {mtext(text = expression("#"~Basis~functions~"in"~each~of~the~(italic(r))~shared~and~(italic(m))~independent~fields), side = 1, srt = 90, cex = 1, xpd = T, outer = T, line = -1)}
  rm(tmp.dat, tmp.avg)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

################################################################################
# Plot of MAE vs k and EDF #####################################################
################################################################################
png(filename = paste0(home.wd, "/figures/mae_vs_k_edf.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3,2)) # adjust according to sub.actual.k
for (i in sort(unique(dat1$actual_k))[sub.actual.k]) {
  # set the index of the loop control variable
  idx <- which(i == sort(unique(dat1$actual_k))[sub.actual.k])
  # determine margins based on the plot number
  if (idx == 1) {
    par(mar = c(2.5, 4.1, 1.5, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 2) {
    par(mar = c(2.5, 2.1, 1.5, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  } else if (idx == 3) {
    par(mar = c(3.25, 4.1, 0.75, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 4) {
    par(mar = c(3.25, 2.1, 0.75, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  } else if (idx == 5) {
    par(mar = c(4, 4.1, 0, 2.1))
    # mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")
  } else if (idx == 6) {
    par(mar = c(4, 2.1, 0, 4.1))
    # mtext("Centered log-Likelihoods", side = 4, line = 2.5, cex = 0.75, col = "black")
  }
  
  # subset the data
  tmp.dat <- dat1[dat1$actual_k == i, ]
  tmp.avg <- scen_summaries[scen_summaries$actual_k == i, ]
  
  # plot the selection metric
  with(tmp.dat, plot(fitted_k, edf, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = "", ylim = range(tmp.dat$edf, na.rm = T), xaxt = "n"))#if (idx %in% 5:6) {NULL} else {xaxt = "n"}))
  for (s in 1:max(dat1$sim)) {
    with(tmp.dat[tmp.dat$sim == s, ], lines(fitted_k, edf, col = ggplot2::alpha("grey10", 0.2)))
  }
  axis(4, ylim = range(tmp.dat$edf), col.axis = "black", col = "black")
  # only add label if the plot is on the right
  if (idx %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Centered AIC", side = 4, line = 2.5, cex = 0.75, col = "black")}
  # add plot title
  mtext(text = c("A", "B", "C", "D", "E", "F")[idx], side = 3, cex = 1, line = 0.25, adj = 0)
  abline(a = -(1 + 5 + 5),b = 7, col = "red", lty = "dashed") # need to work out what this ratio should look like (account for the d+m fields as well as fixed effects etc. that don't count toward EDF)
  
  # plot the performance metric
  par(new = T)
  with(tmp.avg, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n", xaxt = "n"))
  if (idx %in% 5:6) {axis(1, ylim = range(tmp.avg$fitted_k))}
  axis(2, ylim = range(tmp.avg$mae_avg), col.axis = "royalblue")
  with(tmp.avg, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.avg, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  # only add label if the plot is on the left
  if (idx %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  # add x axis labels if the plots are those on the bottom
  if (idx == 5) {mtext(text = expression(Basis~Dimension~(italic(k))~of~each~of~the~italic(d)+italic(m)~Latent~Fields), side = 1, srt = 90, cex = 1, xpd = T, outer = T, line = -1)}
  rm(tmp.dat, tmp.avg)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()