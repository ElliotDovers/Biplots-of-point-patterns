home.wd <- getwd()

# collate or load the intercept-only results
if (file.exists("collated_results.RDATA")) {
  load("collated_results.RDATA")
} else {
  # initialise the result storage
  dat <- NULL
  
  # change into appropriate result folder
  setwd(paste(home.wd, "results", sep = "/"))
  # get a list of all the individual result files
  res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
  # inner loop through individual sim-by-model files
  for (job in res.list) {
    # load in the individual simulation results
    load(job)
    # add to the data
    dat <- rbind(dat, res)
    rm(res)
  }
  setwd(home.wd)
  save(list = "dat", file = "collated_results.RDATA")
}

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

library(dplyr)
library(ggplot2)
plot.res <- 500

tmp.dat <- dat0 %>% group_by(true_k, fitted_k) %>% summarise(mae_avg = mean(mae), mae_sd = sd(mae), bic_avg = mean(bic), bic_sd = sd(bic), edf_avg = mean(edf), edf_sd = sd(edf))

png(filename = paste0(home.wd, "/figures/mae_vs_k_bic.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2), mar = c(4.1, 3.1, 0, 3.1))
for (i in sort(unique(dat0$true_k))[1:6]) {
  with(tmp.dat[tmp.dat$true_k == i, ], plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                                            xlab = "", ylab = expression(paste(MAE,": |", mu, " - ", hat(mu), "|"))))
  if (which(i == sort(unique(dat0$true_k))[1:6]) %in% 5:6) {title(xlab = "Basis Dim. (k)")}
  with(tmp.dat[tmp.dat$true_k == i, ], arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5))
  abline(v = i, col = "green4", lty = "dashed")
  if (which(i == sort(unique(dat0$true_k))[1:6]) %in% c(1,3,5)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2, cex = 0.5)}
  par(new = T)
  with(tmp.dat[tmp.dat$true_k == i, ], plot(fitted_k, bic_avg, col = "darkgrey", yaxt = "n", pch = NA, ylab = "", xlab = ""))
  with(tmp.dat[tmp.dat$true_k == i, ], lines(fitted_k, bic_avg, col = "darkgrey"))
  axis(4, ylim = range(tmp.dat$bic_avg[tmp.dat$true_k == i]), col.axis = "darkgrey", col = "darkgrey")
  if (which(i == sort(unique(dat0$true_k))[1:6]) %in% c(2,4,6)) {mtext("BIC", side = 4, line = 2, cex = 0.75, col = "darkgrey")}
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png(filename = paste0(home.wd, "/figures/mae_vs_k_bic_boxplots.png"), res = plot.res, width = 5.5 * plot.res, height = 4.5 * plot.res)
par(mfrow = c(3, 2), mar = c(4.1, 3.1, 0, 3.1))
for (i in sort(unique(dat0$true_k))[1:6]) {
  boxplot(dat0$bic[dat0$true_k == i] ~ dat0$fitted_k[dat0$true_k == i], outline = F, ylab = "", xlab = "", staplelty = 0)
  if (which(i == sort(unique(dat0$true_k))[1:6]) %in% 5:6) {title(xlab = "Basis Dim. (k)")}
  if (which(i == sort(unique(dat0$true_k))[1:6]) %in% c(1,3,5)) {mtext("BIC", side = 2, line = 2, cex = 0.75)}
  par(new = T)
  with(tmp.dat[tmp.dat$true_k == i, ], plot(factor(fitted_k), mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                                            xlab = "", ylab = "", yaxt = "n", col = "white", border = alpha("white", 0)))
  with(tmp.dat[tmp.dat$true_k == i, ], points(factor(fitted_k), mae_avg, col = "royalblue"))
  with(tmp.dat[tmp.dat$true_k == i, ], lines(factor(fitted_k), mae_avg, col = "royalblue"))
  with(tmp.dat[tmp.dat$true_k == i, ], arrows(y0 = mae_avg - mae_sd, x0 = 1:length(unique(fitted_k)), y1 = mae_avg + mae_sd, x1 = 1:length(unique(fitted_k)), code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = which(i == sort(unique(dat0$true_k))), col = "red", lty = "dotted", lwd = 2)
  axis(4, ylim = range(tmp.dat$bic_avg[tmp.dat$true_k == i]), col.axis = "royalblue")
  if (which(i == sort(unique(dat0$true_k))[1:6]) %in% c(2,4,6)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 4, line = 2, cex = 0.5, col = "royalblue")}
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png(filename = paste0(home.wd, "/figures/mae_vs_k_edf.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2), mar = c(4.1, 3.1, 0, 3.1))
for (i in sort(unique(tmp.dat$true_k))[1:6]) {
  with(tmp.dat[tmp.dat$true_k == i, ], plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                                            xlab = "", ylab = expression(paste(MAE,": |", mu, " - ", hat(mu), "|"))))
  if (which(i == sort(unique(tmp.dat$true_k))[1:6]) %in% 5:6) {title(xlab = "Basis Dim. (k)")}
  with(tmp.dat[tmp.dat$true_k == i, ], arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5))
  abline(v = i, col = "green4", lty = "dashed")
  if (which(i == sort(unique(tmp.dat$true_k))[1:6]) == 1) {
    legend("bottomright", lty = c(rep("dotted", 3), "dashed"), col = c("red", "goldenrod", "goldenrod1", "green4"),
           legend = c("1:1", "2:3", "1:2", "True k"), bty = "n", text.col = c(rep("darkgrey", 3), "black"))
  }
  if (which(i == sort(unique(tmp.dat$true_k))[1:6]) %in% c(1,3,5)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2, cex = 0.5)}
  par(new = T)
  with(tmp.dat[tmp.dat$true_k == i, ], plot(fitted_k, edf_avg, col = "darkgrey", yaxt = "n", pch = NA, ylab = "", xlab = ""))
  with(tmp.dat[tmp.dat$true_k == i, ], lines(fitted_k, edf_avg, col = "darkgrey"))
  abline(a = 0, b = 1, col = "red", lty = "dotted")
  abline(a = 0, b = 2/3, col = "goldenrod", lty = "dotted")
  abline(a = 0, b = 1/2, col = "goldenrod1", lty = "dotted")
  axis(4, ylim = range(tmp.dat$bic_avg[tmp.dat$true_k == i]), col.axis = "darkgrey", col = "darkgrey")
  if (which(i == sort(unique(tmp.dat$true_k))[1:6]) %in% c(2,4,6)) {mtext("EDF", side = 4, line = 2, cex = 0.75, col = "darkgrey")}
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# With a covariate #############################################################

dat_avg <- dat1 %>% group_by(true_k, set_k, fitted_k, actual_k) %>% summarise(mae_avg = mean(mae), mae_sd = sd(mae), bic_avg = mean(bic), bic_sd = sd(bic),
                                                             edf_avg = mean(edf), edf_sd = sd(edf), ll_avg = mean(ll, na.rm = T), ll_sd = sd(ll, na.rm = T), aic_avg = mean(aic, na.rm = T), aic_sd = sd(aic, na.rm = T), ll_na = sum(is.na(ll)), aic_na = sum(is.na(ll)),
                                                             rmse1 = sqrt(mean(se_beta1)), rmse2 = sqrt(mean(se_beta2)), rmse3 = sqrt(mean(se_beta3)), rmse4 = sqrt(mean(se_beta4)), rmse5 = sqrt(mean(se_beta5)),
                                                             cp1 = mean(cover_beta1), cp2 = mean(cover_beta2), cp3 = mean(cover_beta3), cp4 = mean(cover_beta4), cp5 = mean(cover_beta5)
)

# calculate the edf to k ratio
dat$edf.ratio <- dat$edf_avg / (2 * dat$fitted_k)

# toggle to show a subset of the panels
sub.actual.k <- 1:6

png(filename = paste0(home.wd, "/figures/mae_vs_k_edf_covar.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2)) # adjust according to sub.actual.k
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n"))
  axis(2, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}
  # set up a smoother so we can be approximate where the line crosses ##########
  tmp.m <- mgcv::gam(edf_avg ~ s(fitted_k), data = tmp.dat)
  xs <- min(tmp.dat$fitted_k):max(tmp.dat$fitted_k)
  cross.idx <- min(which(predict(tmp.m, newdata = data.frame(fitted_k = xs)) / (xs*2) < 0.5))
  rect(xleft = xs[cross.idx],
       ybottom = min(tmp.dat$mae_avg) - 10, xright = max(tmp.dat$fitted_k) + 50,
       ytop = max(tmp.dat$mae_avg) * 1.2, border = NA, col = "grey")
  ##############################################################################
  with(tmp.dat, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) == 1) {
    legend("bottomright", lty = c("dashed", "dotted"), col = c("red", "royalblue"),
           legend = c("1:2 (EDF:k * d)", "True k"), bty = "n", text.col = c("black", "royalblue"))
  }
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  par(new = T)
  with(tmp.dat, plot(fitted_k, edf_avg, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = ""))
  with(tmp.dat, lines(fitted_k, edf_avg, col = "black"))
  abline(a = 0, b = 1, col = "red", lty = "dashed")
  axis(4, ylim = range(tmp.dat$edf_avg), col.axis = "black", col = "black")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext(expression(Avg.~EDF~of~italic(d)~Fields), side = 4, line = 2.5, cex = 0.75, col = "black")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# works because there is not much variability in EDF, as in:
png(filename = paste0(home.wd, "/figures/mae_vs_k_edf_covar_boxplots.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2), mar = c(4.1, 3.1, 0, 3.1))
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  boxplot(dat1$edf[dat1$actual_k == i] ~ dat1$fitted_k[dat1$actual_k == i], outline = F, ylab = "", xlab = "", staplelty = 0)
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {title(xlab = "Basis Dim. (k)")}
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext("EDF", side = 2, line = 2, cex = 0.75)}
  par(new = T)
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(factor(fitted_k), mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                                            xlab = "", ylab = "", yaxt = "n", col = "white", border = alpha("white", 0)))
  with(tmp.dat, points(factor(fitted_k), mae_avg, col = "royalblue"))
  with(tmp.dat, lines(factor(fitted_k), mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = 1:length(unique(fitted_k)), y1 = mae_avg + mae_sd, x1 = 1:length(unique(fitted_k)), code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = which(i == sort(unique(dat$actual_k))), col = "blue", lty = "dotted", lwd = 2)
  axis(4, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 4, line = 2, cex = 0.5, col = "royalblue")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# can also look for signal in log-likelihoods, however we don't see it since the scale over all sims obscures it:
png(filename = paste0(home.wd, "/figures/mae_vs_k_ll_covar_boxplots.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2), mar = c(4.1, 3.1, 0, 3.1))
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  boxplot(dat1$ll[dat1$actual_k == i] ~ dat1$fitted_k[dat1$actual_k == i], outline = F, ylab = "", xlab = "", staplelty = 0)
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {title(xlab = "Basis Dim. (k)")}
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext("log-Likelihood", side = 2, line = 2, cex = 0.75)}
  par(new = T)
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(factor(fitted_k), mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     xlab = "", ylab = "", yaxt = "n", col = "white", border = alpha("white", 0)))
  with(tmp.dat, points(factor(fitted_k), mae_avg, col = "royalblue"))
  with(tmp.dat, lines(factor(fitted_k), mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = 1:length(unique(fitted_k)), y1 = mae_avg + mae_sd, x1 = 1:length(unique(fitted_k)), code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = which(i == sort(unique(dat$actual_k))), col = "blue", lty = "dotted", lwd = 2)
  axis(4, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 4, line = 2, cex = 0.5, col = "royalblue")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

tmp <- dat1[dat1$true_k == 200, ]
plot(tmp$fitted_k, tmp$ll, type = "n", ylim = c(-3, 2))
for (s in 1:100) {
  lines(tmp$fitted_k[tmp$sim == s], (tmp$ll[tmp$sim == s] - mean(tmp$ll[tmp$sim == s], na.rm = T)) / sd(tmp$ll[tmp$sim == s], na.rm = T))
}

# calculate the means for each metric path, for each simulation (so they can be standardised)
path_summaries <- dat1 %>% group_by(actual_k, sim) %>% summarise(ll_avg = mean(ll, na.rm = T), aic_avg = mean(aic, na.rm = T), bic_avg = mean(bic, na.rm = T),
                                                                 ll_sd = sd(ll, na.rm = T), aic_sd = sd(aic, na.rm = T), bic_sd = sd(bic, na.rm = T)
)

# loop through actual k scenarios and simulations to calculate the standardised metrics
dat1$ll_std <- NULL
dat1$aic_std <- NULL
dat1$bic_std <- NULL
curve_df <- list()
for (k in sort(unique(dat1$actual_k))) {
  tmp.curves <- NULL
  for (s in 1:max(dat1$sim)) {
    dat1$ll_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$ll[dat1$actual_k == k & dat1$sim == s] - path_summaries$ll_avg[path_summaries$actual_k == k & path_summaries$sim == s]) / path_summaries$ll_sd[path_summaries$actual_k == k & path_summaries$sim == s]
    dat1$aic_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$aic[dat1$actual_k == k & dat1$sim == s] - path_summaries$aic_avg[path_summaries$actual_k == k & path_summaries$sim == s]) / path_summaries$aic_sd[path_summaries$actual_k == k & path_summaries$sim == s]
    dat1$bic_std[dat1$actual_k == k & dat1$sim == s] <- (dat1$bic[dat1$actual_k == k & dat1$sim == s] - path_summaries$bic_avg[path_summaries$actual_k == k & path_summaries$sim == s]) / path_summaries$bic_sd[path_summaries$actual_k == k & path_summaries$sim == s]
    tmp.curves <- cbind(tmp.curves, dat1$ll_std[dat1$actual_k == k & dat1$sim == s])
  }
}

# But on average there is a clear signal
png(filename = paste0(home.wd, "/figures/mae_vs_k_ll_covar.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(5, 3)) # adjust according to sub.actual.k
for (i in sort(unique(dat1$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat1[dat1$actual_k == i, ]
  tmp.dat <- dat_avg[dat_avg$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n"))
  axis(2, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}
  
  with(tmp.dat, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  par(new = T)
  
  with(tmp.dat, plot(fitted_k, ll_std, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = "", ylim = range(tmp.dat$ll_std, na.rm = T)))
  for (s in 1:max(dat1$sim)) {
    with(tmp.dat[tmp.dat$sim == s, ], lines(fitted_k, ll_std, col = ggplot2::alpha("grey10", 0.2)))
  }
  axis(4, ylim = range(tmp.dat$ll_std), col.axis = "black", col = "black")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Avg. log-Likelihood", side = 4, line = 2.5, cex = 0.75, col = "black")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# But on average there is a clear signal
png(filename = paste0(home.wd, "/figures/mae_vs_k_ll_covar.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2)) # adjust according to sub.actual.k
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n"))
  axis(2, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}

  with(tmp.dat, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  par(new = T)
  with(tmp.dat, plot(fitted_k, ll_avg, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = ""))
  with(tmp.dat, lines(fitted_k, ll_avg, col = "black"))
  axis(4, ylim = range(tmp.dat$ll_avg), col.axis = "black", col = "black")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Avg. log-Likelihood", side = 4, line = 2.5, cex = 0.75, col = "black")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# can also check with AIC and BIC - spoiler: they penalise too hard
png(filename = paste0(home.wd, "/figures/mae_vs_k_aic_covar.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2)) # adjust according to sub.actual.k
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n"))
  axis(2, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}
  
  with(tmp.dat, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  par(new = T)
  with(tmp.dat, plot(fitted_k, aic_avg, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = ""))
  with(tmp.dat, lines(fitted_k, aic_avg, col = "black"))
  axis(4, ylim = range(tmp.dat$aic_avg), col.axis = "black", col = "black")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Avg. AIC", side = 4, line = 2.5, cex = 0.75, col = "black")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png(filename = paste0(home.wd, "/figures/mae_vs_k_bic_covar.png"), res = plot.res, width = 6.2 * plot.res, height = 6.2 * plot.res)
par(mfrow = c(3, 2)) # adjust according to sub.actual.k
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, mae_avg, ylim = range(c(mae_avg + mae_sd, mae_avg - mae_sd)),
                     pch = NA, xlab = "", ylab = "", yaxt = "n"))
  axis(2, ylim = range(tmp.dat$mae_avg), col.axis = "royalblue")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}
  
  with(tmp.dat, points(fitted_k, mae_avg, col = "royalblue"))
  with(tmp.dat, arrows(y0 = mae_avg - mae_sd, x0 = fitted_k, y1 = mae_avg + mae_sd, x1 = fitted_k, code = 3, angle = 90, length = 0.025, lwd = 1.5, col = "royalblue"))
  abline(v = i, col = "royalblue", lty = "dotted")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(paste(MAE,": |", mu, " - ", hat(mu), "|")), side = 2, line = 2.5, cex = 0.75, col = "royalblue")}
  par(new = T)
  with(tmp.dat, plot(fitted_k, bic_avg, col = "black", yaxt = "n", pch = NA, ylab = "", xlab = ""))
  with(tmp.dat, lines(fitted_k, bic_avg, col = "black"))
  axis(4, ylim = range(tmp.dat$bic_avg), col.axis = "black", col = "black")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(2, length(sub.actual.k), by = 2)) {mtext("Avg. BIC", side = 4, line = 2.5, cex = 0.75, col = "black")}
  rm(tmp.dat)
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# Other metrics ################################################################

png(filename = paste0(home.wd, "/figures/rmse_covar.png"), res = plot.res, width = 5.5 * plot.res, height = 4.5 * plot.res)
par(mfrow = c(3, 2)) # adjust according to sub.actual.k
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, rmse1, ylim = range(c(rmse1, rmse2, rmse3, rmse4, rmse5)), xlab = "", ylab = "", pch = NA))
  with(tmp.dat, lines(fitted_k, rmse1, col = 2))
  with(tmp.dat, lines(fitted_k, rmse2, col = 3))
  with(tmp.dat, lines(fitted_k, rmse3, col = 4))
  with(tmp.dat, lines(fitted_k, rmse4, col = 5))
  with(tmp.dat, lines(fitted_k, rmse5, col = 6))
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}
  abline(v = i, lty = "dashed")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(RMSE~hat(beta)[j]), side = 2, line = 2.5, cex = 0.75)}
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png(filename = paste0(home.wd, "/figures/cover_prob_covar.png"), res = plot.res, width = 5.5 * plot.res, height = 4.5 * plot.res)
par(mfrow = c(3, 2)) # adjust according to sub.actual.k
for (i in sort(unique(dat$actual_k))[sub.actual.k]) {
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {par(mar = c(4.1, 4.1, 0, 2.1))} else {par(mar = c(4.1, 2.1, 0, 4.1))}
  tmp.dat <- dat[dat$actual_k == i, ]
  with(tmp.dat, plot(fitted_k, cp1, ylim = range(c(cp1, cp2, cp3, cp4, cp5)), xlab = "", ylab = "", pch = NA))
  with(tmp.dat, lines(fitted_k, cp1, col = 2))
  with(tmp.dat, lines(fitted_k, cp2, col = 3))
  with(tmp.dat, lines(fitted_k, cp3, col = 4))
  with(tmp.dat, lines(fitted_k, cp4, col = 5))
  with(tmp.dat, lines(fitted_k, cp5, col = 6))
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% (length(sub.actual.k)-1):length(sub.actual.k)) {mtext(expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), side = 1, line = 2.5 , cex = 0.75)}
  abline(v = i, lty = "dashed")
  if (which(i == sort(unique(dat$actual_k))[sub.actual.k]) %in% seq(1, length(sub.actual.k), by = 2)) {mtext(expression(C.I.~Coverage~beta[j]), side = 2, line = 2.5, cex = 0.75)}
}
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()
