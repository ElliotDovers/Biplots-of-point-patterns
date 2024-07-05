home.wd <- getwd()

library(glmmTMB)
library(spatstat)
library(gclus)
library(corrplot)
library(ggplot2)

# set the plot resolution
plot.res <- 500
# set the working dir
home.wd <- getwd()
dir.create(file.path(home.wd, "figures"))

# get a list of the model info.
file.list <- list.files("model_checks")
# initialise storage
res.tab <- NULL
# loop through and load the individual results
for (i in 1:length(file.list)) {
  load(paste0("model_checks/", file.list[i]))
  if (!is.null(res.tab)) {
    for (i_col in 1:ncol(res.tab)) {
      if (!colnames(res.tab)[i_col] %in% colnames(res)) {
        eval(parse(text = paste0("res$", colnames(res.tab)[i_col],  " <- NA")))
      }
    }
  }
  res.tab <- rbind(res.tab, res)
  rm(res)
}
# re-order
res.tab <- res.tab[order(res.tab$p, res.tab$d, res.tab$bfs), ]

# select out the comparisons
res_d2_xi1 <- res.tab[res.tab$tmp == 1 & res.tab$d == 2, ]
res_d2_xi2 <- res.tab[res.tab$tmp == 2 & res.tab$d == 2, ]# & res.tab$msg == "relative convergence (4)" & !is.na(res.tab$msg), ]
res_d4_xi1 <- res.tab[res.tab$tmp == 1 & res.tab$d == 4, ]
res_d4_xi2 <- res.tab[res.tab$tmp == 2 & res.tab$d == 4, ]# & res.tab$msg == "relative convergence (4)" & !is.na(res.tab$msg), ]
# res_d2_xi2$ll[res_d2_xi2$msg == "singular convergence (7)"] <- NA

png(filename = paste0(home.wd, "/figures/nsw_model_selection_r2.png"), width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
par(mar = c(4.1, 4.1, 4.1, 0.5))
plot(res_d2_xi2$bfs, res_d2_xi2$ll, ylim = range(res_d2_xi1$ll, res_d2_xi2$ll, res_d4_xi1$ll, res_d4_xi2$ll, na.rm = T), #pch = 4, col = "royalblue",
     xlab = expression("#"~Basis~functions~"in"~each~of~the~latent~fields), ylab = "log-Likelihood")
lines(res_d2_xi2$bfs, res_d2_xi2$ll)#, col = "royalblue")
points(res_d2_xi1$bfs, res_d2_xi1$ll, pch = 5)#, col = "goldenrod3")
lines(res_d2_xi1$bfs, res_d2_xi1$ll, lty = "dotted")#col = "goldenrod3")
abline(v = res_d2_xi1$bfs[23], lty = "dashed")
abline(v = res_d2_xi2$bfs[35], lty = "dashed")
abline(v = res_d2_xi2$bfs[13], lty = "dashed")
text(x = res_d2_xi2$bfs[13] + 20, y = -14325, labels = expression(k[a]))
text(x = res_d2_xi1$bfs[23] + 20, y = -14325, labels = expression(k[b]))
text(x = res_d2_xi2$bfs[35] + 20, y = -14325, labels = expression(k[c]))
legend(x = 0, y = -13990, legend = c(expression(b[i] == Lambda~u[i]~+~epsilon[i]), expression(b[i] == Lambda~u[i]), "candidate k explored"), xpd = T, pch = c(1, 5, NA), lty = c("solid", "dotted", "dashed"), bty = "n")#, col = c("royalblue", "goldenrod3"))
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

##

png(filename = paste0(home.wd, "/figures/nsw_model_selection.png"), width = 6.2 * plot.res, height = 3.5 * plot.res, res = plot.res)
par(mfrow = c(1 , 2), mar = c(3.5, 3, 1.5, 0))
plot(res_d2_xi2$bfs, res_d2_xi2$ll, ylim = range(res_d2_xi1$ll, res_d2_xi2$ll, res_d4_xi1$ll, res_d4_xi2$ll, na.rm = T), xlim = c(0, 400),#pch = 4, col = "royalblue",
     xlab = "", ylab = "")
mtext(expression(paste("A: ", r==2)), 3, line = 0, cex = 2, at = 75)
mtext("log-Likelihood", 2, line = 1.95)
mtext(expression("#"~Basis~functions~"in"~each~of~the~latent~fields), 1, line = 2.25, at = 450, xpd = T)
lines(res_d2_xi2$bfs, res_d2_xi2$ll)#, col = "royalblue")
points(res_d2_xi1$bfs, res_d2_xi1$ll, pch = 5)#, col = "goldenrod3")
lines(res_d2_xi1$bfs, res_d2_xi1$ll, lty = "dotted")#col = "goldenrod3")
abline(v = res_d2_xi1$bfs[8], lty = "dashed", col = "royalblue")
abline(v = res_d2_xi1$bfs[15], lty = "dashed", col = "royalblue")
# abline(v = res_d2_xi1$bfs[17], lty = "dashed", col = "royalblue")
text(x = res_d2_xi1$bfs[8] + 20, y = -13425, labels = expression(k[a]))
text(x = res_d2_xi1$bfs[15] + 20, y = -13425, labels = expression(k[b]))
# text(x = res_d2_xi1$bfs[17] + 25, y = -13425, labels = expression(k[c]))

par(mar = c(3.5, 2, 1.5, 1))
plot(res_d4_xi2$bfs, res_d4_xi2$ll, ylim = range(res_d2_xi1$ll, res_d2_xi2$ll, res_d4_xi1$ll, res_d4_xi2$ll, na.rm = T), xlim = c(0, 400),#pch = 4, col = "royalblue",
     xlab = "", yaxt = "n", ylab = "")
mtext(expression(paste("B: ", r==4)), 3, line = 0, cex = 2, at = 69)
lines(res_d4_xi2$bfs, res_d4_xi2$ll)#, col = "royalblue")
points(res_d4_xi1$bfs, res_d4_xi1$ll, pch = 5)#, col = "goldenrod3")
lines(res_d4_xi1$bfs, res_d4_xi1$ll, lty = "dotted")#col = "goldenrod3")
legend("bottomright", legend = c(expression(b[i] == Lambda~u[i]~+~epsilon[i]), expression(b[i] == Lambda~u[i]), "candidate k explored"), xpd = T, pch = c(1, 5, NA), lty = c("solid", "dotted", "dashed"), bty = "n", col = c("black","black", "royalblue"))
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png(filename = paste0(home.wd, "/figures/nsw_model_selection_radii.png"), width = 6.2 * plot.res, height = 3.5 * plot.res, res = plot.res)
par(mfrow = c(1 , 2), mar = c(3.5, 3, 1.5, 0))

plot(res_d2_xi2$bf_rad, res_d2_xi2$ll, ylim = range(res_d2_xi1$ll, res_d2_xi2$ll, res_d4_xi1$ll, res_d4_xi2$ll, na.rm = T), #pch = 4, col = "royalblue",
     xlab = "", ylab = "")
mtext(expression(paste("A: ", r==2)), 3, line = 0, cex = 2, at = 150)
mtext("log-Likelihood", 2, line = 1.95)
mtext(expression("#"~Basis~functions~"in"~each~of~the~latent~fields), 1, line = 2.25, at = 700, xpd = T)
lines(res_d2_xi2$bf_rad, res_d2_xi2$ll)#, col = "royalblue")
points(res_d2_xi1$bf_rad, res_d2_xi1$ll, pch = 5)#, col = "goldenrod3")
lines(res_d2_xi1$bf_rad, res_d2_xi1$ll, lty = "dotted")#col = "goldenrod3")
abline(v = res_d2_xi1$bf_rad[23], lty = "dashed", col = "royalblue")
abline(v = res_d2_xi2$bf_rad[35], lty = "dashed", col = "royalblue")
abline(v = res_d2_xi2$bf_rad[13], lty = "dashed", col = "royalblue")
text(x = res_d2_xi2$bf_rad[13] + 30, y = -14325, labels = expression(k[a]))
text(x = res_d2_xi1$bf_rad[23] + 30, y = -14325, labels = expression(k[b]))
text(x = res_d2_xi2$bf_rad[35] + 30, y = -14325, labels = expression(k[c]))

par(mar = c(3.5, 2, 1.5, 1))
plot(res_d4_xi2$bf_rad, res_d4_xi2$ll, ylim = range(res_d2_xi1$ll, res_d2_xi2$ll, res_d4_xi1$ll, res_d4_xi2$ll, na.rm = T), #pch = 4, col = "royalblue",
     xlab = "", yaxt = "n", ylab = "")
mtext(expression(paste("B: ", r==4)), 3, line = 0, cex = 2, at = 100)
lines(res_d4_xi2$bf_rad, res_d4_xi2$ll)#, col = "royalblue")
points(res_d4_xi1$bf_rad, res_d4_xi1$ll, pch = 5)#, col = "goldenrod3")
lines(res_d4_xi1$bf_rad, res_d4_xi1$ll, lty = "dotted")#col = "goldenrod3")
legend("bottomright", legend = c(expression(b[i] == Lambda~u[i]~+~epsilon[i]), expression(b[i] == Lambda~u[i]), "candidate k explored"), xpd = T, pch = c(1, 5, NA), lty = c("solid", "dotted", "dashed"), bty = "n", col = c("black","black", "royalblue"))
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()



# check the same plot for d = 4
tmp1 <- res.tab[res.tab$tmp == 1 & res.tab$d == 4, ]
tmp2 <- res.tab[res.tab$tmp == 2 & res.tab$d == 4, ]# & res.tab$msg == "relative convergence (4)" & !is.na(res.tab$msg), ]
tmp2$ll[tmp2$msg == "singular convergence (7)"] <- NA

png(filename = paste0(home.wd, "/figures/nsw_model_selection_four_latent.png"), width = 6.2 * plot.res, height = 4.5 * plot.res, res = plot.res)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(tmp1$bfs, tmp1$ll, ylim = range(tmp1$ll, tmp2$ll, na.rm = T), pch = 4, col = "royalblue",
     xlab = expression("#"~Basis~functions~"in"~each~of~the~latent~fields), ylab = "log-Likelihood")
lines(tmp1$bfs, tmp1$ll, col = "royalblue")
points(tmp2$bfs, tmp2$ll, pch = 1, col = "goldenrod3")
lines(tmp2$bfs, tmp2$ll, col = "goldenrod3")
abline(v = tmp1$bfs[18], lty = "dashed")
abline(v = tmp2$bfs[15], lty = "dashed")
text(x = tmp2$bfs[15] + 20, y = -14200, labels = expression(k[xi["2"]]))
text(x = tmp1$bfs[18] + 20, y = -14200, labels = expression(k[xi["1"]]))
legend(x = 0, y = -14025, legend = c(expression(paste("Correlated Fields Only: ", xi[1])), expression(paste("Correlated & Independent Fields: ", xi[2]))), xpd = T, pch = c(4, 1), col = c("royalblue", "goldenrod3"), bty = "n")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()


# set up for complete model selection plot
res.tab$mod <- factor(paste0(res.tab$p, res.tab$flds), levels = c("A1", "A2", "B1", "B2", "C1", "C2"),
                      labels = c("Intercepts. +\nCor Flds", "Intercepts +\nCor & Ind Flds", "Linear Env. +\nCor Flds", "Linear Env. +\nCor & Ind Flds", "Quad. Env. +\nCor Flds", "Quad. Env. +\nCor & Ind Flds"))

png(filename = paste0(home.wd, "/figures/nsw_model_selection_all.png"), width = 6.2 * plot.res, height = 4.2 * plot.res, res = plot.res)
ggplot(data = res.tab[res.tab$d != 6, ], aes(x = bfs, y = ll, color = mod, shape = mod)) +
  facet_wrap(~ d, labeller = labeller(d = c(`2` = "r = 2", `4` = "r = 4", `6` = "r = 6"))) +
  geom_point() + geom_path() +
  scale_y_continuous(name = "log-Likelihood") + 
  scale_x_continuous(name = "# Basis Functions", limits = c(0, max(res.tab$bfs, na.rm = T))) +
  scale_color_manual(name = "Model\nSpec.", values = c("goldenrod2", "goldenrod4", "royalblue", "royalblue4", "magenta", "magenta4")) + 
  scale_shape_manual(name = "Model\nSpec.", values = rep(c(4, 1), 3)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black")
  )
dev.off()

################################################################################

# explore the candidate models based on log-likelihood FOR INTERCEPTS ONLY MODEL

# for a simpler model where the ll appears to first peak
job1 = 8

# for the max. ll achieved by mA1
job2 = 17

# for the max. ll acheived by mA2 (or near enough without being > # presences)
job3 = 15

################################################################################

# load in the relevant models
files <- c("mA1_8.RDATA", "mA2_8.RDATA", "mA1_15.RDATA", "mA2_15.RDATA", "mA1_17.RDATA", "mA2_17.RDATA", "mA1_47.RDATA", "mA2_47.RDATA")
for (i in files) {
  load(paste0(home.wd, "/models/", i))
  mod <- sapply(strsplit(sapply(strsplit(i, ".", fixed = T), function(x){x[1]}), "_", fixed = T), function(y){y[1]})
  job <- as.numeric(sapply(strsplit(sapply(strsplit(i, ".", fixed = T), function(x){x[1]}), "_", fixed = T), function(y){y[2]}))
  fld <- if (grepl("1", sapply(strsplit(sapply(strsplit(i, ".", fixed = T), function(x){x[1]}), "_", fixed = T), function(y){y[1]}), fixed = T)) {
    1
  } else {
    2
  }
  ks <- if (job == 8) {
    "a"
  } else if (job == 15) {
    "b"
  } else if (job == 17) {
    "c"
  } else {
    "r4"
  }
  # remove the data frame stored in the model as this wastes space and (I don't think) it affects anything we want to do
  eval(parse(text = paste0(mod, "$frame <- NULL")))
  eval(parse(text = paste0("m", fld, ks, " <- ", mod)))
  eval(parse(text = paste0("rm(", mod, ")")))
}

# get the data #################################################################
source("get_and_format_flora_data.R")
domain$lon <- domain$x; domain$lat <- domain$y; domain$x <- domain$xm; domain$y <- domain$ym; domain$xm <- NULL; domain$ym <- NULL
rm(dat_pa)
# load the species information
load("flora_species.RDATA")

# get the groups to which the flora species belongs to subset the data
flora_groups <- unique(sp.info$group)

# subset the presence-only data (in the case we are using restricted species)
pres <- pres[pres$group %in% flora_groups, ]

# want to use species with at least 10 presences
pres <- pres[pres$spid %in% names(table(pres$spid))[table(pres$spid) >= 10], ]

# get the species
species <- unique(pres$spid)

# re-order the species information
sp.info <- sp.info[order(sp.info$spid), ]
tmp <- strsplit(sp.info$full_name, " ", fixed = T)
sp.info$nice.names <- paste(substr(sapply(tmp, function(x){x[1]}), start = 1, stop = 1), substr(sapply(tmp, function(x){x[2]}), start = 1, stop = 4), sep = ". ")
rm(tmp)
# need to fix the duplicate names of one species (likely an error from the paper)
sp.info$nice.names[sp.info$nice.names == "S. lueh"] <- paste0(sp.info$nice.names[sp.info$nice.names == "S. lueh"], c(" A", " B"))

# crop the sp.info
sp.info <- sp.info[sp.info$spid %in% unique(pres$spid), ]

# Correlation ##################################################################

# function for obtaining the correlation matrix
vcov2corr <- function(V) {
  s <- sqrt(diag(V))
  return(diag(1/s) %*% V %*% diag(1/s))
}

cor1a <- attr(VarCorr(m1a)$cond$basis.functions, "correlation")
cor2a <- vcov2corr(VarCorr(m2a)$cond$basis.functions + VarCorr(m2a)$cond$basis.functions.1)
cor1b <- attr(VarCorr(m1b)$cond$basis.functions, "correlation")
cor2b <- vcov2corr(VarCorr(m2b)$cond$basis.functions + VarCorr(m2b)$cond$basis.functions.1)
cor1c <- attr(VarCorr(m1c)$cond$basis.functions, "correlation")
cor2c <- vcov2corr(VarCorr(m2c)$cond$basis.functions + VarCorr(m2c)$cond$basis.functions.1)
cor1r4 <- attr(VarCorr(m1r4)$cond$basis.functions, "correlation")
cor2r4 <- vcov2corr(VarCorr(m2r4)$cond$basis.functions + VarCorr(m2r4)$cond$basis.functions.1)

# cov1a <- VarCorr(m1a)$cond$basis.functions
# cov2a <- VarCorr(m2a)$cond$basis.functions + VarCorr(m2a)$cond$basis.functions.1
# cov1b <- VarCorr(m1b)$cond$basis.functions
# cov2b <- VarCorr(m2b)$cond$basis.functions + VarCorr(m2b)$cond$basis.functions.1
# cov1c <- VarCorr(m1c)$cond$basis.functions
# cov2c <- VarCorr(m2c)$cond$basis.functions + VarCorr(m2c)$cond$basis.functions.1

dimnames(cor1a) <- list(sp.info$nice.names, sp.info$nice.names)
dimnames(cor2a) <- list(sp.info$nice.names, sp.info$nice.names)
dimnames(cor1b) <- list(sp.info$nice.names, sp.info$nice.names)
dimnames(cor2b) <- list(sp.info$nice.names, sp.info$nice.names)
dimnames(cor1c) <- list(sp.info$nice.names, sp.info$nice.names)
dimnames(cor2c) <- list(sp.info$nice.names, sp.info$nice.names)

png(filename = paste0(home.wd, "/figures/nsw_corplots.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
par(mfrow = c(2, 3), mar = c(0,2.1,2.1,0))
corrplot(cor1a[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
mtext(expression(b[i] == Lambda~u[i]), side = 2, cex = 1.5, line = 0)
# mtext(expression(k[a]==149), cex = 1.5, line = -1)
mtext(expression(k[a]), cex = 1.5, line = -1)
corrplot(cor1b[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
# mtext(expression(k[b]==318), cex = 1.5, line = -1)
mtext(expression(k[b]), cex = 1.5, line = -1)
corrplot(cor1c[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
# mtext(expression(k[c]==580), cex = 1.5, line = -1)
mtext(expression(k[c]), cex = 1.5, line = -1)
corrplot(cor2a[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
mtext(expression(b[i] == Lambda~u[i]~+~epsilon[i]), side = 2, cex = 1.5, line = -0)
corrplot(cor2b[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
corrplot(cor2c[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

png(filename = paste0(home.wd, "/figures/nsw_corplots.png"), width = 5 * plot.res, height = 5 * plot.res, res = plot.res)
par(mfrow = c(2, 2), mar = c(0,2.1,2.1,0))
corrplot(cor1a[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
mtext(expression(paste(bold(b)[i] ==bold("\u039b"~u)[i])), side = 2, cex = 1.5, line = 0)
mtext(expression(k[a]), cex = 1.5, line = -1)
corrplot(cor1b[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
mtext(expression(k[b]), cex = 1.5, line = -1)
corrplot(cor2a[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
mtext(expression(paste(bold(b)[i] ==bold("\u039b"~u)[i]~+~bold("\u03b5")[i])), side = 2, cex = 1.5, line = -0)
corrplot(cor2b[order.single(cor2a), order.single(cor2a)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n")#, cl.pos = "n")
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

# png(filename = paste0(home.wd, "/figures/nsw_covplots.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
# par(mfrow = c(2, 3), mar = c(0,2.1,2.1,0))
# corrplot(cov1a[order.single(cov1b), order.single(cov1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n", is.corr = F)
# mtext(expression(bold(xi[1])), side = 2, las = 1, cex = 2)
# mtext(expression(k[a]), cex = 2, line = -1)
# corrplot(cov1b[order.single(cov1b), order.single(cov1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n", is.corr = F)
# mtext(expression(k[b]), cex = 2, line = -1)
# corrplot(cov1c[order.single(cov1b), order.single(cov1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n", is.corr = F)
# mtext(expression(k[c]), cex = 2, line = -1)
# corrplot(cov2a[order.single(cov1b), order.single(cov1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n", is.corr = F)
# mtext(expression(bold(xi[2])), side = 2, las = 1, cex = 2)
# # mtext(expression(k[a]), cex = 2)
# corrplot(cov2b[order.single(cov1b), order.single(cov1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n", is.corr = F)
# # mtext(expression(k[b]), cex = 2)
# corrplot(cov2c[order.single(cov1b), order.single(cov1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,2.1,0), tl.pos = "n", is.corr = F)
# # mtext(expression(k[c]), cex = 2)
# par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
# dev.off()

# compare the estimates between correlation plots
r2_est <- cor2a[order.single(cor2a), order.single(cor2a)]
r4_est <- cor2r4[order.single(cor2a), order.single(cor2a)]
cor.comp <- r2_est
cor.comp[upper.tri(cor.comp)] <- r4_est[upper.tri(r4_est)]

png(filename = paste0(home.wd, "/figures/nsw_correlation_comparison.png"), width = 5.4 * plot.res, height = 5.2 * plot.res, res = plot.res)
corrplot(cor.comp, method = "square", tl.cex = 1, tl.srt = 45, diag = T, tl.col = "white", mar = c(0,2.1,2.1,0))
abline(a = 23, b = -1, xpd = T, lty = "dashed")
mtext("r = 2", side = 2, line = 2, at = 13, cex = 3)
mtext("r = 4", side = 3, line = 2, at = 18, cex = 3)
dev.off()

# select out some species
selected.sp <- c(2, 4, 5, 6, 7, 26)
tmp.cols <- rep("goldenrod2", nrow(cor2a))
tmp.cols[selected.sp] <- "black"
png(filename = paste0(home.wd, "/figures/nsw_corplot_david.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
# corrplot(cor2a[order.single(cor1b), order.single(cor1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,0,0), is.corr = F, tl.col = tmp.cols[order.single(cor1b)])
corrplot(cor2a[order.single(cor1b), order.single(cor1b)], type = "lower", method = "square", tl.srt = 45, diag = T, mar = c(0,2.1,0,0), is.corr = F, tl.col = "black")
dev.off()

# Latent Fields ################################################################

# extract the fields over the domain
if (file.exists("fld1a.RDATA")) {
  load("fld1a.RDATA")
} else {
  fld1a <- get_field(m1a, domain)
  save(fld1a, file = "fld1a.RDATA")
}
if (file.exists("fld2a.RDATA")) {
  load("fld2a.RDATA")
} else {
  fld2a <- get_field(m2a, domain)
  save(fld2a, file = "fld2a.RDATA")
}
if (file.exists("fld1b.RDATA")) {
  load("fld1b.RDATA")
} else {
  fld1b <- get_field(m1b, domain)
  save(fld1b, file = "fld1b.RDATA")
}
if (file.exists("fld2b.RDATA")) {
  load("fld2b.RDATA")
} else {
  fld2b <- get_field(m2b, domain)
  save(fld2b, file = "fld2b.RDATA")
}
if (file.exists("fld1c.RDATA")) {
  load("fld1c.RDATA")
} else {
  fld1c <- get_field(m1c, domain)
  save(fld1c, file = "fld1c.RDATA")
}
if (file.exists("fld2c.RDATA")) {
  load("fld2c.RDATA")
} else {
  fld2c <- get_field(m2c, domain)
  save(fld2c, file = "fld2c.RDATA")
}

zlims <- range(as.numeric(fld1a[,1:2]), as.numeric(fld2a[,1:2]), as.numeric(fld1b[,1:2]), as.numeric(fld2b[,1:2]), as.numeric(fld1c[,1:2]), as.numeric(fld2c[,1:2]))
# png(filename = paste0(home.wd, "/figures/nsw_latent_fields.png"), width = 6.2 * plot.res, height = 5.8 * plot.res, res = plot.res)
# layout(matrix(c(1, 2, 2, 3, 3, 4:18), nrow = 4, ncol = 5, byrow = T), widths = c(0.05, rep(0.95 / 4, 4)), heights = c(0.05, rep(0.95 / 3, 3)))
# par(mar = c(0,0,1,0))
# plot(1, type = "n", axes = F)
# plot(1, type = "n", axes = F)
# text(x = 1, y = 1.2, labels = expression(xi[1]), cex = 2, xpd = T)
# plot(1, type = "n", axes = F)
# text(x = 1, y = 1.2, labels = expression(xi[2]), cex = 2, xpd = T)
# plot(1, type = "n", axes = F)
# text(x = 1, y = 1, labels = expression(k[a]), cex = 2, xpd = T)
# plot(vec2im(fld1a[,1], domain$x, domain$y), main = expression(paste(bold(z),"(s)",bold(u[1]))), box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld1a[,2], domain$x, domain$y), main = expression(paste(bold(z),"(s)",bold(u[2]))), box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld2a[,1], domain$x, domain$y), main = expression(paste(bold(z),"(s)",bold(u[1]))), box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld2a[,2], domain$x, domain$y), main = expression(paste(bold(z),"(s)",bold(u[2]))), box = F, zlim = zlims, ribbon = F)
# plot(1, type = "n", axes = F)
# text(x = 1, y = 1, labels = expression(k[b]), cex = 2, xpd = T)
# plot(vec2im(fld1b[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld1b[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld2b[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld2b[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(1, type = "n", axes = F)
# text(x = 1, y = 1, labels = expression(k[c]), cex = 2, xpd = T)
# plot(vec2im(fld1c[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld1c[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld2c[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# plot(vec2im(fld2c[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
# par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
# dev.off()

png(filename = paste0(home.wd, "/figures/nsw_latent_fields.png"), width = 3.5 * plot.res, height = 5.8 * plot.res, res = plot.res)
# layout(matrix(c(1,2,3,4,5,6,7,8,5,9,10,11,12,13,14,15,12,16,17,18), nrow = 5, ncol = 4, byrow = T), widths = c(0.05, rep(0.95 / 3, 3)), heights = c(0.05, rep(0.95 / 4, 4)))
par(mfrow = c(4, 3), mar = c(0,0,1.5,0))
plot(vec2im(fld1a[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(bold(xi[1]),":")), side = 2, line = -2.1, las = 1, at = c(240, 6800), cex = 1.5)
mtext(expression(paste(bold(z),"(s)",bold(u[1]))), side = 2, line = -3.5, las = 1, at = c(240, 6700))
mtext(expression(k[a]), side = 3, line = -0.3, at = c(470, 6900), cex = 1.5)
plot(vec2im(fld1b[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(k[b]), side = 3, line = -0.3, at = c(470, 6900), cex = 1.5)
plot(vec2im(fld1c[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(k[c]), side = 3, line = -0.3, at = c(470, 6900), cex = 1.5)

plot(vec2im(fld2a[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(bold(xi[2]),":")), side = 2, line = -2.1, las = 1, at = c(240, 6800), cex = 1.5)
mtext(expression(paste(bold(z),"(s)",bold(u[1]))), side = 2, line = -3.5, las = 1, at = c(240, 6700))
plot(vec2im(fld2b[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
plot(vec2im(fld2c[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)

plot(vec2im(fld1a[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(bold(xi[1]),":")), side = 2, line = -2.1, las = 1, at = c(240, 6800), cex = 1.5)
mtext(expression(paste(bold(z),"(s)",bold(u[2]))), side = 2, line = -3.5, las = 1, at = c(240, 6700))
plot(vec2im(fld1b[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
plot(vec2im(fld1c[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)

plot(vec2im(fld2a[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(bold(xi[2]),":")), side = 2, line = -2.1, las = 1, at = c(240, 6800), cex = 1.5)
mtext(expression(paste(bold(z),"(s)",bold(u[2]))), side = 2, line = -3.5, las = 1, at = c(240, 6700))
plot(vec2im(fld2b[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
plot(vec2im(fld2c[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)

par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()


# Biplotting ###################################################################

# set up the color guide scheme #

# set up a matrix to rotate the coordinates
theta <- -pi/10
R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
rot_coords <- t(R %*% t(apply(domain[,c("x", "y")], 2, function(x){x - min(x)})))
# set up the continuous spatial color scheme
x_colfn <- colorRampPalette(c("white", "black"))
y_colfn <- colorRampPalette(c("blue", "gold"))
domain$xcol <- x_colfn(nrow(domain))[as.numeric(cut(rot_coords[,1],breaks = nrow(domain)))]
domain$ycol <- y_colfn(nrow(domain))[as.numeric(cut(rot_coords[,2],breaks = nrow(domain)))]
for (i in 1:nrow(domain)) {
  domain$xycol[i] <- colorRampPalette(c(domain$xcol[i], domain$ycol[i]))(3)[2]
}
# function to add cols to the basis function nodes
get_bfs_cols <- function(mod) {
  bfs <- mod$basis.functions
  ext.win <- owin(range(c(domain$x, bfs$x)), range(c(domain$y, bfs$y)))
  bfs.pp <- ppp(x = bfs$x, y = bfs$y, window = ext.win)
  domain.pp <- ppp(x = domain$x, y = domain$y, window = ext.win)
  nn.bfs <- nncross(bfs.pp, domain.pp, what = "which")
  domain$xycol[nn.bfs]
}
# set the basis function node colors
m1a$basis.functions$xycol <- get_bfs_cols(m1a)
m2a$basis.functions$xycol <- get_bfs_cols(m2a)
m1b$basis.functions$xycol <- get_bfs_cols(m1b)
m2b$basis.functions$xycol <- get_bfs_cols(m2b)
m1c$basis.functions$xycol <- get_bfs_cols(m1c)
m2c$basis.functions$xycol <- get_bfs_cols(m2c)

# select out some species
selected.sp <- c(2, 4, 5, 6, 7, 19) # with D. anta
selected.sp <- c(2, 4, 5, 6, 7, 14) # with E. lati

# also need smaller labels for species
sp.info$short_label <- ""
sp.info$short_label[selected.sp][c(1,3,6,4,5,2)] <- paste0("sp", 1:6)

png(filename = paste0(home.wd, "/figures/nsw_biplot_comparison.png"), width = 6.2 * plot.res, height = 5.8 * plot.res, res = plot.res)
layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,6,7,8,9,10), nrow = 3, byrow = T), widths = rep(1/6, 6), heights = c(2/5, 2/5, 1/5))
par(mar = c(1,4.1,2.1,0))
biplot(m1a, asp = 1, bty = "n", score.col = m1a$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$short_label, cex = 0.75, rotate.by.theta = 1.7, load.name.cex = 0.9)
mtext(expression(paste(bold(b)[i] ==bold("\u039b"~u)[i])), side = 2, las = 0, cex = 1.5, line = 2.2)
mtext(expression(k[a]), cex = 1.5, line = -0.45)
biplot(m1b, asp = 1, bty = "n", score.col = m1b$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$short_label, cex = 0.75, rotate.by.theta = 1.86, load.name.cex = 0.9)
mtext(expression(k[b]), cex = 1.5, line = -0.45)
# biplot(m1c, asp = 1, bty = "n", score.col = m1c$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$short_label, cex = 0.75, rotate.by.theta = 1.725, load.name.cex = 0.9, alpha = 0.5, xlim = c(-2.2,2))
# mtext(expression(k[c]), cex = 1.5, line = -0.45)
par(mar = c(2.1,4.1,0,0))
biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$short_label, cex = 0.75, rotate.by.theta = 1.6, load.name.cex = 0.9)
mtext(expression(paste(bold(b)[i] ==bold("\u039b"~u)[i]~+~bold("\u03b5")[i])), side = 2, las = 0, cex = 1.5, line = 2.2)
biplot(m2b, asp = 1, bty = "n", score.col = m2b$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$short_label, cex = 0.75, rotate.by.theta = 1.64, load.name.cex = 0.9)#, alpha = 0.5, xlim = c(-2,1.5))
# biplot(m2c, asp = 1, bty = "n", score.col = m2c$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$short_label, cex = 0.75, rotate.by.theta = 1.875, load.name.cex = 0.9, alpha = 0.5)
par(mar = c(0,1.5,2.1,0))
for (i in sp.info$spid[selected.sp][c(1,3,6,4,5,2)]) {
  plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
  points(pres[pres$spid == i, c("xm", "ym")], pch = 4, cex = 0.6, col = "black")
  mtext(paste0(sp.info$nice.names[sp.info$spid == i], " (sp", which(i == sp.info$spid[selected.sp][c(1,3,6,4,5,2)]), ")\n "), side = 3, line = -1, adj = 0, cex = 0.85, at = 200, xpd = T)
}
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

eg.sp <- c(4, 7)
zlims_1 <- range(as.numeric(fld2a[, paste0("m", eg.sp[1])]), as.numeric(fld2b[, paste0("m", eg.sp[1])]), as.numeric(fld2c[, paste0("m", eg.sp[1])])) 
zlims_2 <- range(as.numeric(fld2a[, paste0("m", eg.sp[2])]), as.numeric(fld2b[, paste0("m", eg.sp[2])]), as.numeric(fld2c[, paste0("m", eg.sp[2])])) 
png(filename = paste0(home.wd, "/figures/nsw_eg_overfitting.png"), width = 4.5 * plot.res, height = 5 * plot.res, res = plot.res)
par(mfrow = c(2, 3), mar = c(0,4.6,0,0))
plot(vec2im(fld2a[, paste0("m", eg.sp[1])], domain$x, domain$y), main = "", box = F, ribbon = F, zlim = zlims_1)
points(pres[pres$spid == sp.info$spid[eg.sp[1]], c("xm", "ym")], pch = 4, cex = 0.75)
# mtext(expression(paste("E. blak: ", bold(z), "(", bold(s), ")", bold(v[j]))), side = 2, cex = 2, line = 0.5)
mtext(expression(paste(bold(z), "(", bold(s), ")", bold(v[j]))), side = 2, cex = 2, line = 1.5, at = c(248, 6260))
mtext("E. blak", side = 2, cex = 2, line = 0)
mtext(expression(k[a]), line = -2.5, at = c(470, 6900), cex = 2)
par(mar = c(0,2.3,0,2.3))
plot(vec2im(fld2b[, paste0("m", eg.sp[1])], domain$x, domain$y), main = "", box = F, ribbon = F, zlim = zlims_1)
points(pres[pres$spid == sp.info$spid[eg.sp[1]], c("xm", "ym")], pch = 4, cex = 0.75)
mtext(expression(k[b]), line = -2.5, at = c(470, 6900), cex = 2)
par(mar = c(0,0,0,1.5))
plot(vec2im(fld2c[, paste0("m", eg.sp[1])], domain$x, domain$y), main = "", box = F, ribbon = T, zlim = zlims_1)
points(pres[pres$spid == sp.info$spid[eg.sp[1]], c("xm", "ym")], pch = 4, cex = 0.75)
mtext(expression(k[c]), line = -2.5, at = c(470, 6900), cex = 2)
par(mar = c(0,4.6,0,0))
plot(vec2im(fld2a[, paste0("m", eg.sp[2])], domain$x, domain$y), main = "", box = F, ribbon = F, zlim = zlims_2)
points(pres[pres$spid == sp.info$spid[eg.sp[2]], c("xm", "ym")], pch = 4, cex = 0.75)
# mtext(expression(paste("E. camp: ", bold(z), "(", bold(s), ")", bold(v[j]))), side = 2, cex = 2, line = 0.5)
mtext("E. camp", side = 2, cex = 2, line = 0)
par(mar = c(0,2.3,0,2.3))
plot(vec2im(fld2b[, paste0("m", eg.sp[2])], domain$x, domain$y), main = "", box = F, ribbon = F, zlim = zlims_2)
points(pres[pres$spid == sp.info$spid[eg.sp[2]], c("xm", "ym")], pch = 4, cex = 0.75)
par(mar = c(0,0,0,1.5))
plot(vec2im(fld2c[, paste0("m", eg.sp[2])], domain$x, domain$y), main = "", box = F, ribbon = T, zlim = zlims_2)
points(pres[pres$spid == sp.info$spid[eg.sp[2]], c("xm", "ym")], pch = 4, cex = 0.75)
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

# main paper plot ##############################################################

png(filename = paste0(home.wd, "/figures/nsw_biplot.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
# layout(matrix(c(1,1,1,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2/3, 1/3))
layout(matrix(c(1,1,2,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2/3, 1/3))
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
mtext("A: Biplot with Domain Color Guide", side = 3, line = 0.5, adj = 0, xpd = T)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(m2a$basis.functions$x, m2a$basis.functions$y, bg = m2a$basis.functions$xycol, pch = 21, cex = 1.2)
legend(x = 0.1, y = -0.1, legend = c("Basis function nodes", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(3.1,2.1,1,0), xpd = F)
# biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$nice.names, cex = 0.95, rotate.by.theta = 0.12, load.name.cex = 1.2, show.all.arrows = T)
biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$nice.names, cex = 0.95, rotate.by.theta = 1.6, load.name.cex = 1.2, show.all.arrows = T, ylim = c(-2.1,1.8), xlim = c(-1.7, 1.7), yaxt = "n")#, alpha = 0.525)
axis(side = 2, at = c(-2, -1, 0, 1))
par(mar = c(0,1.5,2.0,0))
for (i in sp.info$spid[selected.sp][c(1,3,6,4,5,2)]) {
  plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
  points(pres[pres$spid == i, c("xm", "ym")], pch = 4, cex = 0.6, col = "black")
  if (i == sp.info$spid[selected.sp][1]) {
    mtext("B: North-Eastern NSW", side = 3, line = 1.5, adj = 0)
  }
  # mtext(paste0(sp.info$nice.names[sp.info$spid == i], "\n(", sp.info$group[sp.info$spid == i], ")"), side = 3, line = -1.4, adj = 0, cex = 0.85, at = c(260, 7200), xpd = T)
  mtext(sp.info$nice.names[sp.info$spid == i], side = 3, line = -0.65, adj = 0, cex = 0.85, at = c(260, 7300), xpd = T)
}
dev.off()

## for talk:

png(filename = paste0(home.wd, "/figures/talk_plot4.png"), width = 6.2 * plot.res, height = 5 * plot.res, res = plot.res)
# layout(matrix(c(1,1,1,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2/3, 1/3))
layout(matrix(c(1,1,2,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2/3, 1/3))
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# mtext("A: Biplot with Domain Color Guide", side = 3, line = 0.5, adj = 0, xpd = T)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(m2a$basis.functions$x, m2a$basis.functions$y, bg = m2a$basis.functions$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function nodes", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(3.1,2.1,0.5,0), xpd = F)
biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, show.responses = selected.sp, xlab = "", ylab = "", load.names = sp.info$nice.names, cex = 0.85, rotate.by.theta = pi/1.392, load.name.cex = 1.2, show.all.arrows = T)
par(mar = c(0,1.5,2.1,0))
for (i in sp.info$spid[selected.sp][c(1,3,6,4,5,2)]) {
  plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
  points(pres[pres$spid == i, c("xm", "ym")], pch = 4, cex = 0.6, col = "black")
  if (i == sp.info$spid[selected.sp][1]) {
    mtext("North-Eastern NSW", side = 3, line = 1.5, adj = 0)
  }
  # mtext(paste0(sp.info$nice.names[sp.info$spid == i], "\n(", sp.info$group[sp.info$spid == i], ")"), side = 3, line = -1.4, adj = 0, cex = 0.85, at = c(260, 7200), xpd = T)
  mtext(sp.info$nice.names[sp.info$spid == i], side = 3, line = -0.65, adj = 0, cex = 0.85, at = c(260, 7300), xpd = T)
}
dev.off()


png(filename = paste0(home.wd, "/figures/nsw_biplot_all_species.png"), width = 6.2 * plot.res, height = 5.8 * plot.res, res = plot.res)
layout(matrix(1:2, nrow = 1, byrow = T), widths = c(0.35,0.65), heights = 1)
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
mtext("Biplot with Domain Color Guide", side = 3, line = 0.5, adj = 0, xpd = T)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(m2a$basis.functions$x, m2a$basis.functions$y, bg = m2a$basis.functions$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function nodes", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(3.1,2.1,0.5,0), xpd = F)
biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, xlab = "", ylab = "", load.names = sp.info$nice.names, cex = 0.85, rotate.by.theta = 1.6, load.name.cex = 1)
mtext(expression(paste(bold(u[1]))), side = 1, line = 2, at = 0)
mtext(expression(paste(bold(u[2]))), side = 2, line = 2, at = 0)
dev.off()

# working for residual covariance analysis ##########################################
tmp <- get_bfs_coeff(m2a)
zlims_m2a <- range(as.numeric(fld2a[,1:2]))

# the enite species biplot
png(filename = paste0(home.wd, "/figures/nsw_biplot_intercept_only_all_species.png"), width = 6.2 * plot.res, height = 5.8 * plot.res, res = plot.res)
layout(matrix(1:4, nrow = 2, byrow = T), widths = c(0.2,0.8), heights = c(0.8, 0.2))
plot(vec2im(fld2a[,2], domain$x, domain$y), main = "", box = F, zlim = zlims_m2a, ribbon = F)
biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, xlab = "", ylab = "", load.names = sp.info$nice.names, cex = 0.85, load.name.cex = 1)
plot(1, type = "n", axes = F)
plot(vec2im(fld2a[,1], -domain$y, domain$x), main = "", box = F, zlim = zlims_m2a, ribbon = F)


par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
mtext("A: Biplot with Domain Color Guide", side = 3, line = 0.5, adj = 0, xpd = T)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(m2a$basis.functions$x, m2a$basis.functions$y, bg = m2a$basis.functions$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function nodes", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(3.1,2.1,0.5,0), xpd = F)
biplot(m2a, asp = 1, bty = "n", score.col = m2a$basis.functions$xycol, xlab = "", ylab = "", load.names = sp.info$nice.names, cex = 0.85, load.name.cex = 1)
par(mar = c(0,1.5,2.1,0))
plot(vec2im(fld2a[,1], domain$x, domain$y), main = "", box = F, zlim = zlims_m2a, ribbon = F)
plot(1, type = "n", axes = F)
plot(vec2im(fld2a[,2], -domain$y, domain$x), main = "", box = F, zlim = zlims_m2a, ribbon = F)
for (i in sp.info$spid[selected.sp][c(1,3,6,4,5,2)]) {
  plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
  points(pres[pres$spid == i, c("xm", "ym")], pch = 4, cex = 0.6, col = "black")
  if (i == sp.info$spid[selected.sp][1]) {
    mtext("B: North-Eastern NSW", side = 3, line = 1.5, adj = 0)
  }
  # mtext(paste0(sp.info$nice.names[sp.info$spid == i], "\n(", sp.info$group[sp.info$spid == i], ")"), side = 3, line = -1.4, adj = 0, cex = 0.85, at = c(260, 7200), xpd = T)
  mtext(sp.info$nice.names[sp.info$spid == i], side = 3, line = -0.65, adj = 0, cex = 0.85, at = c(260, 7300), xpd = T)
}
dev.off()
################################################################################
