library(spatstat)
library(glmmTMB)
library(corrplot)

# set the plot resolution
plot.res <- 500
# set the working dir
home.wd <- getwd()
dir.create(file.path(home.wd, "figures"))

# extract the regular grid of quadrature included in the glmmTMB data
data("lansing", package = "glmmTMB")
domain <- lansing[lansing$tree == "blackoak" & lansing$pt == 0, ]

# convert the lansing woods data into the required format
data("lansing", package = "spatstat.data")
pp <- as.data.frame(lansing)
colnames(pp)[3] <- "tree"
pp$pt <- 1
pp$wt <- 1e-6
# pp$wt <- 1e-8
# set up the quadrature points
n_q <- 10000
set.seed(1) # settng seed for the random quad pts for reproducability
quad <- data.frame(x = runif(n_q), y = runif(n_q), pt = 0, wt = lansing$window$units$multiplier^2 / n_q) # wt is area / # quad pts
set.seed(NULL)
# create the required replicates of the quad points
tmp <- quad[rep(seq_len(n_q), length(levels(pp$tree))), ]
tmp$tree <- rep(levels(pp$tree), each = n_q) # assign the tree classes
dat <- rbind(pp, tmp)
# # get data without the misc. tree category
# dat_all <- rbind(pp, tmp)
# dat <- dat_all[dat_all$tree != "misc", ]
# dat$tree <- factor(as.character(dat$tree))
rm(tmp)

# fit models with increasing numbers of basis functions and record edf, ll, aic and bic

if (file.exists("bf.search.RDATA")) {
  load("bf.search.RDATA")
} else {
  # initialise storage
  actual_k <- NULL
  ll0 <- NULL
  ll <- NULL
  # edf <- NULL
  # aic <- NULL
  # bic <- NULL

  ks <- c(seq(25, 400, by = 25), 500, 600, 700, 800, 900, 1000, 1500)
  for (i in ks) {
    bfs <- make_basis(k = i, domain)
    actual_k <- c(actual_k, nrow(bfs))
    # fit the model without overdispersion fields
    tmp.m0 <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, mv.flds = "correlated")
    # grab out the starting parameters
    start.pars <- lapply(split(tmp.m0$fit$parfull, names(tmp.m0$fit$parfull)), unname)
    start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * length(levels(dat$tree)))) # add additional field coefs
    start.pars$theta <- c(start.pars$theta, rep(-3, length(levels(dat$tree)))) # add additional field variances
    
    # fit the model
    tmp.m <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, start = start.pars)
    # # calculate the EDF
    # tmp.inf <- influence(tmp.m)
    # tmp.edf <- sum(tmp.inf[(length(fixef(tmp.m)[[1]]) + length(levels(dat$tree)) + 1):length(tmp.inf)])
    # # extract the fitted values
    # tmp.fitted <- predict(tmp.m)
    
    # store results
    # edf <- c(edf, tmp.edf)
    ll0 <- c(ll0, logLik(tmp.m0))
    ll <- c(ll, logLik(tmp.m))
    # aic <- c(aic, -2*as.numeric(logLik(tmp.m)) + sum(tmp.inf) * 2)
    # bic <- c(bic, -2*as.numeric(logLik(tmp.m)) + sum(tmp.inf) * log(sum(dat$pt == 1)))
    
  }
  bf.search <- data.frame(
    set_k = ks, fitted_k = actual_k, ll0, ll#, edf, aic, bic
  )
  save(bf.search, file = "bf.search.RDATA")
}

# # calculate the ratio of edf to k
# bf.search$edf_to_k <- bf.search$edf/(bf.search$fitted_k * 2)

# png(filename = paste0(home.wd, "/figures/lansing_k_edf_and_ll.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
# par(mfrow = c(2, 1), mar = c(2.1, 5.1, 2.1, 0))
# with(bf.search, plot(fitted_k, edf, xaxt = "n", xlab = "", ylab = ""))
# mtext(expression(atop(Combined~EDF~of,the~Latent~Fields)), side = 2, line = 2.5)
# title("A", adj = 0)
# title("B", adj = 0, line = -9.5)
# ### set up a smoother so we can be approximate where the line crosses ##########
# tmp.m <- mgcv::gam(edf ~ s(fitted_k), data = bf.search)
# xs <- min(bf.search$fitted_k):max(bf.search$fitted_k)
# cross.idx <- min(which(predict(tmp.m, newdata = data.frame(fitted_k = xs)) / (xs*3) < 0.5))
# rect(xleft = xs[cross.idx],
#      ybottom = min(bf.search$edf) * 0.69, xright = max(bf.search$fitted_k) + 50,
#      ytop = max(bf.search$edf) * 1.03, border = NA, col = "grey")
# ##############################################################################
# with(bf.search, points(fitted_k, edf))
# with(bf.search, lines(fitted_k, edf))
# abline(a=0,b=1.5,lty="dashed", col = "red", lwd = 1.5)
# legend("bottomright", lty = "dashed", col = "red", legend = "1:2", bty = "n", title = "EDF : k * d")
# par(mar = c(4.1, 5.1, 0.1, 0))
# with(bf.search, plot(fitted_k, ll, xaxt = "n", xlab = expression(Basis~Dimension~of~each~of~the~italic(d)~Latent~Fields), ylab = "log-Likelihood"))
# axis(1, at = c(50, 100, 144, 169, 200, 250, 300, 350, 400), cex.axis = 0.75)
# # rect(xleft = 144, ybottom = min(bf.search$ll) - 3, xright = 169, ytop = max(bf.search$ll) + 3, border = NA, col = "grey")
# with(bf.search, points(fitted_k, ll))
# with(bf.search, lines(fitted_k, ll))
# abline(v = bf.search$fitted_k[which.max(bf.search$ll)], col = "royalblue", lty = "dashed")
# legend("bottomright", lty = "dashed", col = "royalblue", legend = "k for max. log-Lik.", bty = "n", title = "")
# par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()

# png(filename = paste0(home.wd, "/figures/lansing_k.png"), width = 6.2 * plot.res, height = 3.8 * plot.res, res = plot.res)
# par(mar = c(4.1, 5.1, 0.1, 0))
# with(bf.search, plot(fitted_k, ll, xaxt = "n", xlim = c(45,1000), xlab = expression("#"~Basis~functions~"in"~each~of~the~(italic(r))~shared~and~(italic(m))~independent~fields), ylab = "log-Likelihood"))
# axis(1, at = c(50, 150, 225, seq(300, 1000, 100)), cex.axis = 0.75)
# # rect(xleft = 144, ybottom = min(bf.search$ll) - 3, xright = 169, ytop = max(bf.search$ll) + 3, border = NA, col = "grey")
# with(bf.search, points(fitted_k, ll))
# with(bf.search, lines(fitted_k, ll))
# abline(v = bf.search$fitted_k[which.max(bf.search$ll)], col = "royalblue", lty = "dashed")
# legend("bottomright", lty = "dashed", col = "royalblue", legend = "k for max. log-Lik.", bty = "n", title = "")
# par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()

png(filename = paste0(home.wd, "/figures/lansing_k.png"), width = 6.2 * plot.res, height = 3.8 * plot.res, res = plot.res)
par(mar = c(4.1, 4.1, 0.1, 0))
with(bf.search, plot(fitted_k, ll, xaxt = "n", ylim = range(c(ll0, ll)), xlim = c(45,1000), xlab = expression("#"~Basis~functions~"in"~each~of~the~(italic(r))~shared~and~(italic(m))~independent~fields), ylab = "log-Likelihood"))
axis(1, at = c(50, 144, 200, seq(300, 1000, 100)), cex.axis = 0.65)
# rect(xleft = 144, ybottom = min(bf.search$ll) - 3, xright = 169, ytop = max(bf.search$ll) + 3, border = NA, col = "grey")
with(bf.search, points(fitted_k, ll))
with(bf.search, lines(fitted_k, ll))
with(bf.search, points(fitted_k, ll0, pch = 5))
with(bf.search, lines(fitted_k, ll0, lty = "dotted"))
abline(v = bf.search$fitted_k[which.max(bf.search$ll)], col = "royalblue", lty = "dashed")
legend(x = 685, y = -47710, pch = c(1, 5, NA), lty = c("solid", "dotted", "dashed"), col = c("black", "black", "royalblue"), legend = c(expression(b[i] == Lambda~u[i]~+~epsilon[i]), expression(b[i] == Lambda~u[i]),"k for max. log-Lik."), bty = "n", title = "")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# set the selected number of basis functions to use
bfs <- make_basis(k = bf.search$fitted_k[which.max(bf.search$ll)], domain)

# WORKING FOR BASIS FUNCTIONS
tmp.bfs <- bfs
class(tmp.bfs)[2] <- "bf.df"

bfs <- make_basis(k = 64, domain)
plot(vec2im(domain$pt, domain$x, domain$y), xlim = c(-0.25, 1.25), ylim = c(-0.25, 1.25))
symbols(bfs$x, bfs$y, circles = bfs$scale, inches = F, add = T)
points(domain[c(1985, 3200), c("x", "y")])
fields::rdist(domain[c(1985, 8000), c("x", "y")])

rad <- NULL
zz <- NULL
tmp.pts <- domain[c(1985, 3200), c("x", "y")]
tmp.dist <- fields::rdist(tmp.pts)
for (i in 4:400) {
  bfs <- make_basis(k = i, domain)
  z <- bf_matrix(bfs, tmp.pts)
  zz[i] <- sum(z[1, ] * z[2, ])
  rad[i] <- bfs$scale[1]
  rm(bfs, z)
}
plot(rad, zz)
abline(v = tmp.dist[1,2] / 1.5)

# fit the models
m0 <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, mv.flds = "correlated")

start.pars <- lapply(split(m0$fit$parfull, names(m0$fit$parfull)), unname)
start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * length(levels(pp$tree)))) # add additional field coefs
start.pars$theta <- c(start.pars$theta, rep(-3, length(levels(pp$tree)))) # add additional field variances
start.pars$theta <- rep(-3, length(start.pars$theta))
m <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, start = start.pars)

m <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree)
m_opt <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, n_factors = 4, control = glmmTMBControl(optimizer = nlminb, start_method = list(method = "res")))

# calculate K functions and sim envelopes for redoak and maple #################

# subset the data for just Black Oak
dat_bla <- dat[dat$tree == "blackoak", ]
# fit an homogeneous Poisson process
m_bla_ipp <- mvlgcp(pt ~ 1, data = dat_bla, weights = dat_bla$wt)
# create a pixel image of the spatially constant intensity
int.im_bla = vec2im(rep(exp(getME(m_bla_ipp, "beta")) * lansing$window$units$multiplier^2, nrow(domain)), domain$x, domain$y)
int.im_bla$xrange <- c(0,1)
int.im_bla$yrange <- c(0,1)
# set the point pattern as a ppp object for spatstat
blackoak.pp <- ppp(x = dat_bla$x[dat_bla$pt==1], y = dat_bla$y[dat_bla$pt==1], window = owin())
# calculate the observed K function
K_obs_bla <- Kinhom(blackoak.pp, lambda = int.im_bla, correction = "border")
# simulate the envelopes/bounds
K_env_bla <- envelope(blackoak.pp, fun = Kinhom, simulate = expression(rpoispp(lambda = int.im_bla)), verbose = FALSE)

# subset the data for just Hickory
dat_hic <- dat[dat$tree == "hickory", ]
# fit an homogeneous Poisson process
m_hic_ipp <- mvlgcp(pt ~ 1, data = dat_hic, weights = dat_hic$wt)
# create a pixel image of the spatially constant intensity
int.im_hic = vec2im(rep(exp(getME(m_hic_ipp, "beta")) * lansing$window$units$multiplier^2, nrow(domain), nrow(domain)), domain$x, domain$y)
int.im_hic$xrange <- c(0,1)
int.im_hic$yrange <- c(0,1)
# set the point pattern as a ppp object for spatstat
hickory.pp <- ppp(x = dat_hic$x[dat_hic$pt==1], y = dat_hic$y[dat_hic$pt==1], window = owin())
# calculate the observed K function
K_obs_hic <- Kinhom(hickory.pp, lambda = int.im_hic, correction = "border")
# simulate the envelopes/bounds
K_env_hic <- envelope(hickory.pp, fun = Kinhom, simulate = expression(rpoispp(lambda = int.im_hic)), verbose = FALSE)

# subset the data for just Maple
dat_map <- dat[dat$tree == "maple", ]
# fit an homogeneous Poisson process
m_map_ipp <- mvlgcp(pt ~ 1, data = dat_map, weights = dat_map$wt)
# create the pixel images (uses the window supplied in the original gorillas data)
int.im_map = vec2im(rep(exp(getME(m_map_ipp, "beta")) * lansing$window$units$multiplier^2, nrow(domain), nrow(domain)), domain$x, domain$y)
int.im_map$xrange <- c(0,1)
int.im_map$yrange <- c(0,1)
# set the point pattern as a ppp object for spatstat
maple.pp <- ppp(x = dat_map$x[dat_map$pt==1], y = dat_map$y[dat_map$pt==1], window = owin())
# calculate the observed K functions
K_obs_map <- Kinhom(maple.pp, lambda = int.im_map, correction = "border")
# simulate the envelopes/bounds
K_env_map <- envelope(maple.pp, fun = Kinhom, simulate = expression(rpoispp(lambda = int.im_map)), verbose = FALSE)

# subset the data for just Misc. class
dat_mis <- dat[dat$tree == "misc", ]
# fit an homogeneous Poisson process
m_mis_ipp <- mvlgcp(pt ~ 1, data = dat_mis, weights = dat_mis$wt)
# create a pixel image of the spatially constant intensity
int.im_mis = vec2im(rep(exp(getME(m_mis_ipp, "beta")) * lansing$window$units$multiplier^2, nrow(domain), nrow(domain)), domain$x, domain$y)
int.im_mis$xrange <- c(0,1)
int.im_mis$yrange <- c(0,1)
# set the point pattern as a ppp object for spatstat
misc.pp <- ppp(x = dat_mis$x[dat_mis$pt==1], y = dat_mis$y[dat_mis$pt==1], window = owin())
# calculate the observed K function
K_obs_mis <- Kinhom(misc.pp, lambda = int.im_mis, correction = "border")
# simulate the envelopes/bounds
K_env_mis <- envelope(misc.pp, fun = Kinhom, simulate = expression(rpoispp(lambda = int.im_mis)), verbose = FALSE)

# subset the data for just Red Oak
dat_red <- dat[dat$tree == "redoak", ]
# fit an homogeneous Poisson process
m_red_ipp <- mvlgcp(pt ~ 1, data = dat_red, weights = dat_red$wt)
# create a pixel image of the spatially constant intensity
int.im_red = vec2im(rep(exp(getME(m_red_ipp, "beta")) * lansing$window$units$multiplier^2, nrow(domain), nrow(domain)), domain$x, domain$y)
int.im_red$xrange <- c(0,1)
int.im_red$yrange <- c(0,1)
# set the point pattern as a ppp object for spatstat
redoak.pp <- ppp(x = dat_red$x[dat_red$pt==1], y = dat_red$y[dat_red$pt==1], window = owin())
# calculate the observed K function
K_obs_red <- Kinhom(redoak.pp, lambda = int.im_red, correction = "border")
# simulate the envelopes/bounds
K_env_red <- envelope(redoak.pp, fun = Kinhom, simulate = expression(rpoispp(lambda = int.im_red)), verbose = FALSE)

# subset the data for just White Oak
dat_whi <- dat[dat$tree == "whiteoak", ]
# fit an homogeneous Poisson process
m_whi_ipp <- mvlgcp(pt ~ 1, data = dat_whi, weights = dat_whi$wt)
# create a pixel image of the spatially constant intensity
int.im_whi = vec2im(rep(exp(getME(m_whi_ipp, "beta")) * lansing$window$units$multiplier^2, nrow(domain), nrow(domain)), domain$x, domain$y)
int.im_whi$xrange <- c(0,1)
int.im_whi$yrange <- c(0,1)
# set the point pattern as a ppp object for spatstat
whiteoak.pp <- ppp(x = dat_whi$x[dat_whi$pt==1], y = dat_whi$y[dat_whi$pt==1], window = owin())
# calculate the observed K function
K_obs_whi <- Kinhom(whiteoak.pp, lambda = int.im_whi, correction = "border")
# simulate the envelopes/bounds
K_env_whi <- envelope(whiteoak.pp, fun = Kinhom, simulate = expression(rpoispp(lambda = int.im_whi)), verbose = FALSE)

png(filename = paste0(home.wd, "/figures/lansing_inhomK.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
# par(mfrow = c(3,2))
layout(matrix(c(1, 1:7), nrow = 4, byrow = T), heights = c(0.05, rep(0.95/3, 3)), widths = c(0.5, 0.5))

par(mar = c(0,0,0,0))
plot(1:5, 1:5, type = "n", xaxt = "n", yaxt = "n", main = "", bty = "n")
legend(x = 1.75, y = 3.5, legend = c("Observed", "Theoretic", "Sim. Bounds"),
       col = c("red", "black", "grey80"), lty = c("solid", "dashed", "solid"), cex = 1,
       lwd = c(2, 2, 2), bty = "n", xpd = T, horiz = T)

par(mar = c(2.4, 4.1, 1, 0))
plot(K_env_bla$r, K_env_bla$mmean, type = "n", ylim = range(c(K_env_bla$obs, K_env_bla$hi, K_env_bla$lo)),
     ylab = "", xlab = "", main = "", xaxt = "n")
text(x = 0.05, y = 0.175, labels = "blackoak")
polygon(c(rev(K_env_bla$r), K_env_bla$r), c(rev(K_env_bla$hi), K_env_bla$lo), col = 'grey80', border = NA)
lines(K_env_bla$r, K_env_bla$mmean, lty = "dashed")
lines(K_obs_bla$r, K_obs_bla$border, col = "red")

par(mar = c(2.4, 4, 1, 0.1))
plot(K_env_hic$r, K_env_hic$mmean, type = "n", ylim = range(c(K_env_hic$obs, K_env_hic$hi, K_env_hic$lo)),
     ylab = "", xlab = "", main = "", xaxt = "n", yaxt = "n")
text(x = 0.05, y = 0.175, labels = "hickory")
polygon(c(rev(K_env_hic$r), K_env_hic$r), c(rev(K_env_hic$hi), K_env_hic$lo), col = 'grey80', border = NA)
lines(K_env_hic$r, K_env_hic$mmean, lty = "dashed")
lines(K_obs_hic$r, K_obs_hic$border, col = "red")

par(mar = c(2.9, 4.1, 0.5, 0))
plot(K_env_map$r, K_env_map$mmean, type = "n", ylim = range(c(K_env_map$obs, K_env_map$hi, K_env_map$lo)),
     ylab = "", xlab = "", main = "", xaxt = "n")
text(x = 0.05, y = 0.175, labels = "maple")
polygon(c(rev(K_env_map$r), K_env_map$r), c(rev(K_env_map$hi), K_env_map$lo), col = 'grey80', border = NA)
mtext("K function", side = 2, xpd = T,  line = 2.5)
lines(K_env_map$r, K_env_map$mmean, lty = "dashed")
lines(K_obs_map$r, K_obs_map$border, col = "red")

par(mar = c(2.9, 4, 0.5, 0.1))
plot(K_env_mis$r, K_env_mis$mmean, type = "n", ylim = range(c(K_env_mis$obs, K_env_mis$hi, K_env_mis$lo)),
     ylab = "", xlab = "", main = "", xaxt = "n", yaxt = "n")
text(x = 0.05, y = 0.175, labels = "misc")
polygon(c(rev(K_env_mis$r), K_env_mis$r), c(rev(K_env_mis$hi), K_env_mis$lo), col = 'grey80', border = NA)
lines(K_env_mis$r, K_env_mis$mmean, lty = "dashed")
lines(K_obs_mis$r, K_obs_mis$border, col = "red")

par(mar = c(3.4, 4.1, 0, 0))
plot(K_env_red$r, K_env_red$mmean, type = "n", ylim = range(c(K_env_red$obs, K_env_red$hi, K_env_red$lo)),
     ylab = "", xlab = "", main = "")
text(x = 0.05, y = 0.175, labels = "redoak")
polygon(c(rev(K_env_red$r), K_env_red$r), c(rev(K_env_red$hi), K_env_red$lo), col = 'grey80', border = NA)
lines(K_env_red$r, K_env_red$mmean, lty = "dashed")
lines(K_obs_red$r, K_obs_red$border, col = "red")
mtext("Distance", side = 1, xpd = T,  line = 2.4, at = 0.265)

par(mar = c(3.4, 4, 0, 0))
plot(K_env_whi$r, K_env_whi$mmean, type = "n", ylim = range(c(K_env_whi$obs, K_env_whi$hi, K_env_whi$lo)),
     ylab = "", xlab = "", main = "", yaxt = "n")
text(x = 0.05, y = 0.175, labels = "whiteoak")
polygon(c(rev(K_env_whi$r), K_env_whi$r), c(rev(K_env_whi$hi), K_env_whi$lo), col = 'grey80', border = NA)
lines(K_env_whi$r, K_env_whi$mmean, lty = "dashed")
lines(K_obs_whi$r, K_obs_whi$border, col = "red")

dev.off()

# checking K function tests... these are all significant

library(GET)
library(scampr)

res <- list()
for(i in levels(dat$tree)) {
  tmp.m <- scampr(pt ~ 1, data = dat[dat$tree == i, ], quad.weights.name = "wt", include.sre = F)
  k_env <- kfunc_envelopes(tmp.m, return.data = T, nsims = 1000)
  tmp.GET <- global_envelope_test(curve_set(obs = k_env$obs, sim = attr(k_env, "sim_curves")[, -1], r = k_env$r))
  res[[which(i == levels(dat$tree))]] <- tmp.GET
  plot(tmp.GET)
  rm(tmp.m, k_env, tmp.GET)
}

# # Testing complete spatial randomness (CSR)
# #==========================================
# data("lansing", package = "spatstat.data")
# pp.list <- split(lansing)
# res <- list()
# for (sp in 1:length(pp.list)) {
#   X <- pp.list[[sp]]
#   
#   nsim <- 1999 # Number of simulations
#   
#   
#   # Illustration of general workflow for simple hypotheses
#   #=======================================================
#   # First illustrate the general workflow for the test by this example
#   # of CSR test for a point pattern X using the empirical L-function.
#   # Define the argument values at which the functions are evaluated
#   obs.L <- Lest(X, correction="translate")
#   r <- obs.L[['r']]
#   # The test function for the data
#   obs <- obs.L[['trans']] - r
#   # Prepare simulations and calculate test functions for them at same r as 'obs'
#   sim <- matrix(nrow=length(r), ncol=nsim)
#   for(i in 1:nsim) {
#     sim.X <- runifpoint(ex=X) # simulation under CSR
#     sim[, i] <- Lest(sim.X, correction="translate", r=r)[['trans']] - r
#   }
#   # Create a curve_set containing argument values, observed and simulated functions
#   cset <- curve_set(r=r, obs=obs, sim=sim)
#   # Perform the test
#   tmp <- global_envelope_test(cset, type="erl")
#   eval(parse(text = paste0("res_", names(pp.list)[sp], " <- tmp")))
#   plot(tmp) + ggplot2::ylab(expression(italic(hat(L)(r)-r)))
#   rm(X, obs.L, r, obs, sim, sim.X, cset, tmp)
# }


# set up the continuous spatial color scheme
x_colfn <- colorRampPalette(c("white", "black"))
y_colfn <- colorRampPalette(c("blue", "gold"))
domain$xcol <- x_colfn(nrow(domain))[as.numeric(cut(domain$x,breaks = nrow(domain)))]
domain$ycol <- y_colfn(nrow(domain))[as.numeric(cut(domain$y,breaks = nrow(domain)))]
domain$xycol <- NULL
for (i in 1:nrow(domain)) {
  domain$xycol[i] <- colorRampPalette(c(domain$xcol[i], domain$ycol[i]))(3)[2]
}

# set the same biscale on the basis function nodes (using functions from spatstat here to interpolate the colors)
ext.win <- owin(range(c(domain$x, bfs$x)), range(c(domain$y, bfs$y)))
bfs.pp <- ppp(x = bfs$x, y = bfs$y, window = ext.win)
domain.pp <- ppp(x = domain$x, y = domain$y, window = ext.win)
nn.bfs <- nncross(bfs.pp, domain.pp, what = "which")
bfs$xycol <- domain$xycol[nn.bfs]

## OLD PLOT
# png(filename = paste0(home.wd, "/figures/lansing_biplot.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
# layout(matrix(c(1:3,3), nrow = 2), widths = c(1/3, 2/3), heights = c(0.5, 0.5))
# par(mar = rep(2.1, 4))
# plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# # title(main = "Lansing Woods\nPlot Domain\nSpatial Colour Guide", line = -3)
# mtext("Lansing Woods ", side = 1, line = 0.1)
# par(xpd = T)
# points(bfs$x, bfs$y, bg = bfs$xycol, pch = 21)
# plot(vec2im(rep("grey85", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree %in% c("blackoak", "maple") & dat$pt == 1, c("x", "y")], pch = 4,
#        col = c("black", "goldenrod")[(dat$tree[dat$tree %in% c("blackoak", "maple") & dat$pt == 1] == "maple") + 1])
# legend(x = 0, y = 1.45, legend = c("Basis fn. nodes", "Blackoaks", "Maples"), col = c("black", "black", "goldenrod"), pch = c(1, 4, 4), bty = "n", xpd = T)
# par(mar = c(2.1,2.1,2.1,0), xpd = F)
# biplot(m, score.col = bfs$xycol, alpha = 1.15, xlab = "", ylab = "", bty = "n", cex = 1.2)
# dev.off()

# # Full plot
# png(filename = paste0(home.wd, "/figures/lansing_biplot_full.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
# layout(matrix(c(1,2,2,1,2,2,3,4,5), nrow = 3, byrow = T), widths = c(1/3, 1/3, 1/3), heights = c(1/3, 1/3, 1/3))
# par(mar = rep(2.1, 4))
# plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# mtext("A: Biplot with\nDomain Color Guide", side = 3, line = -0.75, adj = 0)
# # mtext("Lansing Woods", side = 3, line = -1)
# par(xpd = T)
# points(bfs$x, bfs$y, bg = bfs$xycol, pch = 21)
# legend(x = 0.15, y = -0.2, legend = c("Basis function nodes", "Black Oaks", "Maples"), col = c("black", "black", "goldenrod"), pch = c(1, 4, 4), bty = "n", xpd = T)
# 
# par(mar = c(2.1,2.1,2.1,0), xpd = F)
# biplot(m, score.col = bfs$xycol, alpha = 1.15, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5)
# 
# # plot(1, type = "n", axes = F)
# 
# par(mar = c(0,2.1,2.1,2.1))
# plot(vec2im(rep("grey85", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree %in% c("blackoak", "maple") & dat$pt == 1, c("x", "y")], pch = 4,
#        col = c("black", "brown", "goldenrod", "red", "white")[as.numeric(dat$tree[dat$tree %in% c("blackoak", "maple") & dat$pt == 1])])
# mtext("B: Lansing Woods", side = 3, line = 0.1, adj = 0)
# 
# par(mar = c(3.1,2.1,2.1,0))
# plot(K_env_red$r, K_env_red$mmean, type = "n", ylim = range(c(K_env_red$obs, K_env_red$hi, K_env_red$lo)),
#      ylab = "", xlab = "", main = "")
# polygon(c(rev(K_env_red$r), K_env_red$r), c(rev(K_env_red$hi), K_env_red$lo), col = 'grey80', border = NA)
# lines(K_env_red$r, K_env_red$mmean, lty = "dashed")
# lines(K_obs_red$r, K_obs_red$border, col = "red")
# mtext("C: Ripley's K: Red Oaks and Maples", side = 3, line = 0.1, adj = 0)
# mtext("K functions", side = 2, line = 1.8, cex = 0.75)
# mtext("distance", side = 1, line = 2, cex = 0.75, at = 0.275)
# 
# par(mar = c(3.1,1,2.1,1.1))
# plot(K_env_map$r, K_env_map$mmean, type = "n", ylim = range(c(K_env_map$obs, K_env_map$hi, K_env_map$lo)),
#      ylab = "", xlab = "", main = "", yaxt = "n")
# polygon(c(rev(K_env_map$r), K_env_map$r), c(rev(K_env_map$hi), K_env_map$lo), col = 'grey80', border = NA)
# lines(K_env_map$r, K_env_map$mmean, lty = "dashed")
# lines(K_obs_map$r, K_obs_map$border, col = "red")
# 
# dev.off()

tmp.names <- c("", "", "", "", "redoak", "whiteoak")
tmp.names <- rep("", 6)

# Plot with all data
png(filename = paste0(home.wd, "/figures/lansing_biplot.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
layout(matrix(c(1,1,2,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2.1/3, 0.9/3))
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
mtext("A: Biplot with\nDomain Color Guide", side = 3, line = -0.75, adj = 0)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(bfs$x, bfs$y, bg = bfs$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function knots", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(2.1,2.1,2.1,0), xpd = F)
biplot(m, score.col = bfs$xycol, load.names = tmp.names, alpha = 1.35, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5, load.lab.offset = -0.05, arrow.head.length = 0.08)
text(x = c(-1.65, -0.6, 0.2, 1.2, -0.35, 0.35), y = c(0.05, -0.95, 0.675, 0.15, -0.1, -0.2), labels = c("misc", "maple", "hickory", "blackoak", "redoak", "whiteoak"), cex = 1.5)
# biplot(m, score.col = bfs$xycol, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5)

# plot(1, type = "n", axes = F)

par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "blackoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "black")
mtext("B: Lansing Woods", side = 3, line = 0.5, adj = 0)
mtext("black oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "hickory" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "brown")
mtext("hickory", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "maple" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "goldenrod")
mtext("maple", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "misc" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "grey25")
mtext("misc", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "redoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "red")
mtext("red oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "whiteoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "seashell1")
mtext("white oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

dev.off()

# for talk:

png(filename = paste0(home.wd, "/figures/talk_plot3A.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
layout(matrix(c(1,1,2,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2.1/3, 0.9/3))
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# mtext("A: Biplot with\nDomain Color Guide", side = 3, line = -0.75, adj = 0)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(bfs$x, bfs$y, bg = bfs$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function knots", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(2.1,2.1,0,0), xpd = F)
biplot(m, score.col = bfs$xycol, alpha = 1.35, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5, load.lab.offset = -0.05)
# biplot(m, score.col = bfs$xycol, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5)

# plot(1, type = "n", axes = F)

par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "blackoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "black")
mtext("Lansing Woods Presence Data:", side = 3, line = 0.5, adj = 0)
mtext("black oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "hickory" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "brown")
mtext("hickory", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "maple" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "goldenrod")
mtext("maple", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "misc" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "grey25")
mtext("misc", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "redoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "red")
mtext("red oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("grey70", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "whiteoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "seashell1")
mtext("white oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

dev.off()

png(filename = paste0(home.wd, "/figures/talk_plot_without_data.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
layout(matrix(c(1,1,2,2,2,2,3,4,5,6,7,8), nrow = 2, byrow = T), widths = rep(1/6, 6), heights = c(2.1/3, 0.9/3))
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# mtext("A: Biplot with\nDomain Color Guide", side = 3, line = -0.75, adj = 0)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(bfs$x, bfs$y, bg = bfs$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function knots", ""), pch = c(1, NA), bty = "n", xpd = T)

par(mar = c(2.1,2.1,0,0), xpd = F)
biplot(m, score.col = bfs$xycol, alpha = 1.35, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5, load.lab.offset = -0.05)
# biplot(m, score.col = bfs$xycol, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5)

# plot(1, type = "n", axes = F)

par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("white", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree == "blackoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "black")
mtext("Lansing Woods Presence Data:", side = 3, line = 0.5, adj = 0, col = "white")
mtext("black oak", side = 3, line = -0.85, adj = 0, cex = 0.9, col = "white")

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("white", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree == "hickory" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "brown")
mtext("hickory", side = 3, line = -0.85, adj = 0, cex = 0.9, col = "white")

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("white", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree == "maple" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "goldenrod")
mtext("maple", side = 3, line = -0.85, adj = 0, cex = 0.9, col = "white")

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("white", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree == "misc" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "grey25")
mtext("misc", side = 3, line = -0.85, adj = 0, cex = 0.9, col = "white")

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("white", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree == "redoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "red")
mtext("red oak", side = 3, line = -0.85, adj = 0, cex = 0.9, col = "white")

# par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("white", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
# points(dat[dat$tree == "whiteoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6)#, col = "seashell1")
mtext("white oak", side = 3, line = -0.85, adj = 0, cex = 0.9, col = "white")

dev.off()

# compute the correlation matrix
cov.bf <- VarCorr(m)[[1]]$basis.functions # the var-cov matrix from the shared fields
cov.bf1 <- VarCorr(m)[[1]]$basis.functions.1 # the var-cov matrix from the independent fields
cov.bf_opt <- VarCorr(m_opt)[[1]]$basis.functions # the var-cov matrix from the shared fields
cov.bf1_opt <- VarCorr(m_opt)[[1]]$basis.functions.1 # the var-cov matrix from the independent fields

tmp.z <- bf_matrix(m$basis.functions, domain[,c("x", "y")])

tmp.scale <- NULL
for (i in 1:nrow(tmp.z)) {
  tmp.scale[i] <- sum(tmp.z[i, ] * tmp.z[i, ])
}

cor.comp_mu <- cov2cor(cov.bf + cov.bf1)
cor.comp_mu_opt <- cov2cor(cov.bf_opt + cov.bf1_opt)
cor.comp_mu_opt2 <- cov2cor(cov.bf_opt + cov.bf1_opt + diag(rep(VarCorr(m_opt)[[1]]$tree, length(levels(dat$tree)))))
cor.comp_mu <- cov2cor(mean(tmp.scale) * (cov.bf + cov.bf1))
cor.comp_xi <- cov2cor(cov.bf)
dimnames(cor.comp_mu) <- list(levels(dat$tree), levels(dat$tree))
dimnames(cor.comp_xi) <- list(levels(dat$tree), levels(dat$tree))

waag_res_xi_shared <- c(0.365, -0.605, -0.820, -0.975, -0.320, 0.435, 0.03, 0.155, -0.02, 0.05, 0.015, -0.655, 0.105, 0.09, -0.42)
waag_res_mu <- c(0.305, -0.565, -0.680, -0.905, -0.265, 0.410, 0.015, 0.085, -0.01, 0.025, 0.01, -0.405, 0.07, 0.06, -0.19)

cor.comp_mu[upper.tri(cor.comp_mu)] <- waag_res_mu
cor.comp_xi[upper.tri(cor.comp_xi)] <- waag_res_xi_shared

png(filename = paste0(home.wd, "/figures/lansing_correlation_comparison_mu.png"), width = 5.4 * plot.res, height = 5.2 * plot.res, res = plot.res)
corrplot(cor.comp_mu, method = "square", tl.cex = 1, tl.srt = 45, diag = T, tl.col = "black", addCoef.col = "black", mar = c(0,2.1,2.1,0))
abline(a = 7, b = -1, xpd = T, lty = "dashed")
mtext("Estimates from the proposed approach", side = 2, line = 3, at = 3)
mtext("Estimates from Waagepetersen et al. (2016)", side = 3, line = 2.5, at = 4.25)
dev.off()

png(filename = paste0(home.wd, "/figures/lansing_correlation_comparison_xi.png"), width = 5.4 * plot.res, height = 5.2 * plot.res, res = plot.res)
corrplot(cor.comp_xi, method = "square", tl.cex = 1, tl.srt = 45, diag = T, tl.col = "black", addCoef.col = "black", mar = c(0,2.1,2.1,0))
abline(a = 7, b = -1, xpd = T, lty = "dashed")
mtext("Estimates from the proposed approach", side = 2, line = 3, at = 3)
mtext("Estimates from Waagepetersen et al. (2016)", side = 3, line = 2.5, at = 4.25)
dev.off()

# Plot with all data
png(filename = paste0(home.wd, "/figures/lansing_biplot.png"), width = 6.2 * plot.res, height = 4.8 * plot.res, res = plot.res)
layout(matrix(c(1,1,2,2,2,3,4,5,6,7), nrow = 2, byrow = T), widths = rep(1/5, 5), heights = c(2/3, 1/3))
par(mar = c(2.1, 1.5, 2.1, 2.1))
plot(vec2im(domain$xycol, domain$x, domain$y), valuesAreColors = T, box = F, main = "")
mtext("A: Biplot with\nDomain Color Guide", side = 3, line = -0.75, adj = 0)
# mtext("Lansing Woods", side = 3, line = -1)
par(xpd = T)
points(bfs$x, bfs$y, bg = bfs$xycol, pch = 21)
legend(x = 0.1, y = -0.1, legend = c("Basis function nodes", "Tree Locations"), pch = c(1, 4), bty = "n", xpd = T)

par(mar = c(2.1,2.1,2.1,0), xpd = F)
biplot(m, score.col = bfs$xycol, alpha = 1.8, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5)
# biplot(m, score.col = bfs$xycol, xlab = "", ylab = "", bty = "n", cex = 1, load.name.cex = 1.5)

# plot(1, type = "n", axes = F)

par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("grey50", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "blackoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6, col = "black")
mtext("B: Lansing Woods", side = 3, line = 0.5, adj = 0)
mtext("black oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey50", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "hickory" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6, col = "brown")
mtext("hickory", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey50", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "maple" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6, col = "goldenrod")
mtext("maple", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,2.1,2.1,0))
plot(vec2im(rep("grey50", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "redoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6, col = "red")
mtext("red oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

# par(mar = c(0,1.5,2.1,0))
plot(vec2im(rep("grey50", nrow(domain)), domain$x, domain$y), valuesAreColors = T, box = F, main = "")
points(dat[dat$tree == "whiteoak" & dat$pt == 1, c("x", "y")], pch = 4, cex = 0.6, col = "seashell1")
mtext("white oak", side = 3, line = -0.85, adj = 0, cex = 0.9)

dev.off()

# Set up the Waagerpetersen 2016 comparison

# include the misc. tree category
dat <- lansing

# need to again assess the adequacy of number of basis functions on the new data

if (file.exists("bf.search_incl_misc.RDATA")) {
  load("bf.search_incl_misc.RDATA")
} else {
  # initialise storage
  actual_k <- NULL
  ll <- NULL

  ks <- c(seq(25, 400, by = 25), 500, 600, 700, 800, 900, 1000, 1500)
  for (i in ks) {
    bfs <- make_basis(k = i, domain)
    actual_k <- c(actual_k, nrow(bfs))
    # fit the model
    tmp.m <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree)
    # store results
    ll <- c(ll, logLik(tmp.m))
    
  }
  bf.search_misc <- data.frame(
    set_k = ks, fitted_k = actual_k, ll
  )
  save(bf.search_misc, file = "bf.search_incl_misc.RDATA")
}

png(filename = paste0(home.wd, "/figures/lansing_k_incl_misc.png"), width = 6.2 * plot.res, height = 3.8 * plot.res, res = plot.res)
par(mar = c(4.1, 5.1, 0.1, 0))
with(bf.search_misc, plot(fitted_k, ll, xaxt = "n", xlab = expression("#"~Basis~functions~"in"~each~of~the~(italic(r))~shared~and~(italic(m))~independent~fields), ylab = "log-Likelihood"))
# axis(1, at = c(50, 100, 150, 200, 225, 250, 300, 350, 400), cex.axis = 0.75)
axis(1, at = c(50, 144, seq(250, 1000, 100)), cex.axis = 0.75)
# rect(xleft = 144, ybottom = min(bf.search$ll) - 3, xright = 169, ytop = max(bf.search$ll) + 3, border = NA, col = "grey")
with(bf.search_misc, points(fitted_k, ll))
with(bf.search_misc, lines(fitted_k, ll))
abline(v = bf.search_misc$fitted_k[which.max(bf.search_misc$ll)], col = "royalblue", lty = "dashed")
legend("topright", lty = "dashed", col = "royalblue", legend = "k for max. log-Lik.", bty = "n", title = "")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# set the selected number of basis functions to use
bfs_incl_misc <- make_basis(k = bf.search_misc$fitted_k[which.max(bf.search_misc$ll)], domain)

# fit the model that includes the misc. category
m_misc <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs_incl_misc, response.id = dat$tree)

# compute the correlation matrix
cov.bf <- VarCorr(m_misc)[[1]]$basis.functions # the var-cov matrix from the shared fields
cov.bf1 <- VarCorr(m_misc)[[1]]$basis.functions.1 # the var-cov matrix from the independent fields
cor.comp_mu <- cov2cor(cov.bf + cov.bf1)
cor.comp_xi <- cov2cor(cov.bf)
dimnames(cor.comp_mu) <- list(levels(dat$tree), levels(dat$tree))
dimnames(cor.comp_xi) <- list(levels(dat$tree), levels(dat$tree))

waag_res_xi_shared <- c(0.365, -0.605, -0.820, -0.975, -0.320, 0.435, 0.03, 0.155, -0.02, 0.05, 0.015, -0.655, 0.105, 0.09, -0.42)
waag_res_mu <- c(0.305, -0.565, -0.680, -0.905, -0.265, 0.410, 0.015, 0.085, -0.01, 0.025, 0.01, -0.405, 0.07, 0.06, -0.19)

cor.comp_mu[upper.tri(cor.comp_mu)] <- waag_res_mu
cor.comp_xi[upper.tri(cor.comp_xi)] <- waag_res_xi_shared

png(filename = paste0(home.wd, "/figures/lansing_correlation_comparison_mu.png"), width = 5.4 * plot.res, height = 5.2 * plot.res, res = plot.res)
corrplot(cor.comp_mu, method = "square", tl.cex = 1, tl.srt = 45, diag = T, tl.col = "black", addCoef.col = "black", mar = c(0,2.1,2.1,0))
abline(a = 7, b = -1, xpd = T, lty = "dashed")
mtext("Estimates from the proposed approach", side = 2, line = 3, at = 3)
mtext("Estimates from Waagepetersen et al. (2016)", side = 3, line = 2.5, at = 4.25)
dev.off()

png(filename = paste0(home.wd, "/figures/lansing_correlation_comparison_xi.png"), width = 5.4 * plot.res, height = 5.2 * plot.res, res = plot.res)
corrplot(cor.comp_xi, method = "square", tl.cex = 1, tl.srt = 45, diag = T, tl.col = "black", addCoef.col = "black", mar = c(0,2.1,2.1,0))
abline(a = 7, b = -1, xpd = T, lty = "dashed")
mtext("Estimates from the proposed approach", side = 2, line = 3, at = 3)
mtext("Estimates from Waagepetersen et al. (2016)", side = 3, line = 2.5, at = 4.25)
dev.off()

## Convergence of the quadrature approximation #################################

library(glmmTMB)

# extract the regular grid of quadrature included in the glmmTMB data
data("lansing", package = "glmmTMB")
domain <- lansing[lansing$tree == "blackoak" & lansing$pt == 0, ]

# convert the lansing woods data into the required format
data("lansing", package = "spatstat.data")
pp <- as.data.frame(lansing)
colnames(pp)[3] <- "tree"
pp$pt <- 1
pp$wt <- 1e-6
# set up the quadrature points
n_q <- 10000
set.seed(1) # settng seed for the random quad pts for reproducability
quad <- data.frame(x = runif(n_q), y = runif(n_q), pt = 0, wt = lansing$window$units$multiplier^2 / n_q) # wt is area / # quad pts
set.seed(NULL)
# create the required replicates of the quad points
tmp <- quad[rep(seq_len(n_q), length(levels(pp$tree))), ]
tmp$tree <- rep(levels(pp$tree), each = n_q) # assign the tree classes
dat <- rbind(pp, tmp)

# set the basis functions
bfs <- make_basis(k = 150, domain)

# set the scope for q
tmp.qs <- c(1000, 5000, seq(10000, 100000, by = 20000), 100000)

# set the repetitions for each value of q
reps <- 10

# set the loop variable for q (with repetitions)
n_q <- rep(tmp.qs, each = reps)
rep.idx <- rep(1:reps, times = length(tmp.qs))
rm(tmp.qs)

# initialise storage
ll0 <- NULL
ll <- NULL

set.seed(1)
if (file.exists("quad.conv.RDATA")) {
  load("quad.conv.RDATA")
} else {

  for (i in n_q) {

    # create the data
    quad <- data.frame(x = runif(i), y = runif(i), pt = 0, wt = lansing$window$units$multiplier^2 / i)
    tmp <- quad[rep(seq_len(i), length(levels(pp$tree))), ]
    tmp$tree <- rep(levels(pp$tree), each = i) # assign the tree classes
    dat <- rbind(pp, tmp)
    rm(quad, tmp)
    
    # fit the model without overdispersion fields
    tmp.m0 <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, mv.flds = "correlated")
    # grab out the starting parameters
    start.pars <- lapply(split(tmp.m0$fit$parfull, names(tmp.m0$fit$parfull)), unname)
    start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * length(levels(dat$tree)))) # add additional field coefs
    start.pars$theta <- c(start.pars$theta, rep(-3, length(levels(dat$tree)))) # add additional field variances
    
    # fit the full model
    tmp.m <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree, start = start.pars)
 
    # store results
    ll0 <- c(ll0, logLik(tmp.m0))
    ll <- c(ll, logLik(tmp.m))
    rm(dat, tmp.m0, tmp.m)
    
  }
  quad.conv <- data.frame(
    q = n_q, rep = rep.idx, ll0, ll
  )
  save(quad.conv, file = "quad.conv.RDATA")
}

library(ggplot2)
library(dplyr)
png(filename = paste0(home.wd, "/figures/quad_convergence.png"), width = 6.2 * plot.res, height = 4.2 * plot.res, res = plot.res)
ggplot(dat, aes(x = q, y = ll)) +
  geom_point() + 
  geom_line(data = dat %>% group_by(q) %>% summarise(y_mean = mean(ll)), aes(x = q, y = y_mean)) +
  geom_vline(aes(xintercept = 10000, color = "red"), lty = "dashed") +
  scale_x_continuous(name = "# of randomly sampled quadrature points (q)", breaks = c(1000, seq(10000, 100000, by = 10000)), labels = paste0(c(1000, seq(10000, 100000, by = 10000)) / 1000, "K")) +
  scale_y_continuous(name = "Approximate log-Likelihood") +
  scale_color_manual(values = "red", name = "", labels = expression(q==10000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"))
dev.off()

# Fitted Fields ################################################################

sp_flds <- get_field(m, domain, which.response = 1:6)
lat_flds <- get_field(m, domain)

zlims <- range(cbind(sp_flds, lat_flds))
Lambda <- get_loadings(m)

tmp <- VarCorr(m)
sigma_j <- diag(tmp$cond$basis.functions.1)

png(filename = paste0(home.wd, "/figures/lansing_fitted_fields.png"), width = 6.2 * plot.res, height = 3.5 * plot.res, res = plot.res)
layout(matrix(c(1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,13,14,15,16,17,18), nrow = 3, ncol = 7, byrow = T), widths = c(0.1, rep(0.9/6, 6)), heights = rep(1/3,3))
par(mar = c(0,0.1,1.5,0))

plot(1, axes = F, bty = "n", type = "n")

plot(vec2im(sp_flds[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(dat[dat$tree==levels(dat$tree)[1] & dat$pt == 1, c("x", "y")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = levels(dat$tree)[1])
mtext(expression(paste(ln~mu[j],"(s)")), side = 2, line = 0.5, las = 1, xpd = T)
mtext(bquote(lambda[1]==.(round(Lambda[1,1],1))~","~.(round(Lambda[1,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(dat[dat$tree==levels(dat$tree)[2] & dat$pt == 1, c("x", "y")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = levels(dat$tree)[2])
mtext(bquote(lambda[2]==.(round(Lambda[2,1],1))~","~.(round(Lambda[2,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)


plot(vec2im(sp_flds[,3], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(dat[dat$tree==levels(dat$tree)[3] & dat$pt == 1, c("x", "y")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = levels(dat$tree)[3])
mtext(bquote(lambda[3]==.(round(Lambda[3,1],1))~","~.(round(Lambda[3,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,4], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(dat[dat$tree==levels(dat$tree)[4] & dat$pt == 1, c("x", "y")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = levels(dat$tree)[4])
mtext(bquote(lambda[4]==.(round(Lambda[4,1],1))~","~.(round(Lambda[4,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,5], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(dat[dat$tree==levels(dat$tree)[5] & dat$pt == 1, c("x", "y")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = levels(dat$tree)[5])
mtext(bquote(lambda[5]==.(round(Lambda[5,1],1))~","~.(round(Lambda[5,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(vec2im(sp_flds[,6], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
points(dat[dat$tree==levels(dat$tree)[6] & dat$pt == 1, c("x", "y")], col = ggplot2::alpha("black", 0.2), cex = 0.5)
title(main = levels(dat$tree)[6])
mtext(bquote(lambda[6]==.(round(Lambda[6,1],1))~","~.(round(Lambda[6,2],1))), side = 1, line = 0, xpd = T, cex = 0.6)

plot(1, axes = F, bty = "n", type = "n")

plot(vec2im(lat_flds[,1], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(xi[1]^(u),"(s)")), side = 2, line = -3, las = 1, xpd = T)

plot(vec2im(lat_flds[,2], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(xi[2]^(u),"(s)")), side = 2, line = -3, las = 1, xpd = T)

par(mar = c(0,2,4,1))
legend_image <- as.raster(matrix(rev(viridisLite::plasma(20)), ncol=1))
plot(c(0,8),c(zlims[1],zlims[2]),type = 'n', axes = F, xlab = "", ylab = "", main = "")
mtext("field values", line = 0.5, cex = 0.75, adj = 0, at = 1.75)
text(x=1.9 + 3, y = seq(round(zlims[1]),round(zlims[2]),l=5), labels = seq(round(zlims[1]),round(zlims[2]),l=5), xpd = T)
lines(c(0,1.01,1.01,0,0) + 3, c(zlims[1],zlims[1],zlims[2],zlims[2],zlims[1]))
for (i in seq(round(zlims[1]),round(zlims[2]),l=5)) {
  lines(c(1.01, 1.21) + 3, rep(i, 2))
}
rasterImage(legend_image, 0 + 3, zlims[1], 1 + 3, zlims[2])

par(mar = c(1.5,0.1,1.5,0))
plot(1, axes = F, bty = "n", type = "n")

plot(vec2im(lat_flds[,3], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(expression(paste(xi[j]^(epsilon),"(s)")), side = 2, line = 0.5, las = 1, xpd = T)
mtext(bquote(sigma[1]^2==.(round(sigma_j[1],2))), side = 1, line = 0.5, xpd = T, cex = 0.6)
plot(vec2im(lat_flds[,4], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(bquote(sigma[2]^2==.(round(sigma_j[2],2))), side = 1, line = 0.5, xpd = T, cex = 0.6)
plot(vec2im(lat_flds[,5], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(bquote(sigma[3]^2==.(round(sigma_j[3],2))), side = 1, line = 0.5, xpd = T, cex = 0.6)
plot(vec2im(lat_flds[,6], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(bquote(sigma[4]^2==.(round(sigma_j[4],2))), side = 1, line = 0.5, xpd = T, cex = 0.6)
plot(vec2im(lat_flds[,7], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(bquote(sigma[5]^2==.(round(sigma_j[5],2))), side = 1, line = 0.5, xpd = T, cex = 0.6)
plot(vec2im(lat_flds[,8], domain$x, domain$y), main = "", box = F, zlim = zlims, ribbon = F)
mtext(bquote(sigma[6]^2==.(round(sigma_j[6],2))), side = 1, line = 0.5, xpd = T, cex = 0.6)

par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()


# WORKING ######################################################################

scores <- get_bfs_coeff(m)[,1:2]
loadings <- get_loadings(m)

t(scores) %*% scores
scores <- scores %*% solve(diag(sqrt(diag(t(scores) %*% scores))))

s <- apply(scores, 2, norm, type = "2")
s_l <- apply(loadings, 2, norm, type = "2")
new.scores <- scores %*% diag(1/s)
t(new.scores) %*% new.scores

scores <- scores %*% diag(1/s)
loadings <- diag(sqrt(diag(loadings  %*% t(loadings)))) %*% loadings

diag(sqrt(diag(loadings  %*% t(loadings)))) %*% loadings %*% t(loadings) %*% diag(sqrt(diag(loadings  %*% t(loadings))))

plot(scores, asp = 1)
arrows(x0 = rep(0, nrow(loadings)), y0 = rep(0, nrow(loadings)), x1 = loadings[,1], y1 = loadings[,2], length = 0.0, angle = 30, lwd = 0.5)
text(x = loadings[,1], y = loadings[,2], labels = c("blackoak", "hickory", "maple", "redoak", "whiteoak"))
