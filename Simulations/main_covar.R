library(glmmTMB)

# determine job number from pbs script
job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# load the job array
load("job_array.RDATA")

################################################################################
# set parameters of the simulation according to the job array ##################
n_knots <- tab$true_k[tab$job == job]
# fit_k <- tab$fit_k[tab$job == job] # TRYING THIS WITHIN A LOOP
################################################################################

dat <- lansing[lansing$tree != "misc", ]
dat$tree <- factor(as.character(dat$tree))
domain <- dat[dat$tree == "blackoak" & dat$pt == 0, ]

bfs <- make_basis(k = n_knots, domain)

gen_mod <- mvlgcp(pt ~ (1|tree), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$tree)

# extract the random intercepts
intercepts <- getME(gen_mod, "b")[1:length(levels(dat$tree))]

# set the true beta values for each species response to the covariate
true.beta <- c(getME(gen_mod, "beta"), (1:length(levels(dat$tree)) - 3) / 2)

# extract the design matrices from the generating model ########################

# get the call
call.list <- as.list(gen_mod$call)
# indicate not to fit the model
call.list$doFit <- FALSE
# get the model structure (by re-running the model without fitting)
mod_str <- eval(as.call(call.list)) # NOTE: this gets the corrected Z for smoothers since it evaluates mvlgcp()

# adjust the random effect design matrix to be multiplied by factor loadings ###

# set the index of the basis function random effects
i_re <- 2

# extract the factor loadings
Lambda <- gen_mod$obj$env$report(gen_mod$fit$parfull)$fact_load[[i_re]]

# expand the factor loadings by blockReps to align with Z matrix
IxLambda <- kronecker(diag(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps), Lambda) # NOTE: THIS WAY IMPLIES WE NEED Z TO BE k1:resp1, k1:resp2,...,k2:resp1, k2:resp2

# get the current random effect indices within Z
n_re <- unlist(lapply(gen_mod$modelInfo$reStruc$condReStruc, function(x){ x[["blockReps"]] * x[["blockSize"]] }))
if (i_re == 1) {
  # for the first random effects start at 1
  re_idx <- 1:n_re[[i_re]]
} else {
  # for other random effects start at the previous index + 1 and up to the cumulative sum of random effect indices
  re_idx <- (cumsum(n_re)[[i_re - 1]] + 1):cumsum(n_re)[[i_re]]
}
# compute the new Z component
newZ_comp <- mod_str$data.tmb$Z[ , re_idx] %*% IxLambda
# put the current Z matrix back together again 
newZ <- as.matrix(cbind(mod_str$data.tmb$Z[ , 1:ncol(mod_str$data.tmb$Z) < min(re_idx)],
                        newZ_comp,
                        mod_str$data.tmb$Z[ , 1:ncol(mod_str$data.tmb$Z) > max(re_idx)])
)
# add the covariate into the fixed effect design matrix ########################

# create the appropriate design matrix to include a covariate by tree interaction
tmp <- glmmTMB(pt ~ x:tree, data = dat, doFit = F)
newX <- tmp$data.tmb$X

# combine into a single design matrix
des.mat <- as.matrix(cbind(newX, newZ))

# generate new latent factors
set.seed(tab$sim[tab$job == job])
Tilde_u <- rnorm(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps * 2)
b_independent_fld <- rnorm(n_re[[3]])
set.seed(NULL)

# combine all parameters in order to multiply by the above
pars <- c(true.beta, intercepts, Tilde_u, b_independent_fld)

# compute the log-intensity over all the data (we will extract the domain from the quadrature for each species)
flds <- des.mat %*% pars

# simulate point patterns from the intensity fields ############################

# initialise some storage
pp.list <- list() # list of each response point pattern
gen_flds <- NULL # data frame of generated fields
# par(mfrow = c(3,2), mar = c(0,0,1,0))
for (i in 1:length(levels(dat$tree))) {
  # get the fitted field from tree-specific quadrature
  gen_flds <- cbind(gen_flds, exp(flds[dat$tree == levels(dat$tree)[i] & dat$pt == 0]))
  pp.list[[i]] <- spatstat.random::rpoispp(lambda = vec2im(gen_flds[ , i], domain$x, domain$y))
  # plot(vec2im(gen_flds[,i], domain$x, domain$y), main = levels(dat$tree)[i])
  # plot(pp.list[[i]], add = T)
}
# par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

# combine into a data frame of new point patterns
pres <- as.data.frame(do.call("rbind", lapply(pp.list, data.frame)))
# add in the tree identifier
pres$tree <- factor(rep(levels(dat$tree), times = unlist(lapply(pp.list, function(x){x$n}))))
# add in the remaining columns
pres$pt <- 1
pres$wt <- 1e-10

# add in the quadrature for each species (replicates of the supplied background points)
tmp <- domain[rep(seq_len(nrow(domain)), length(levels(dat$tree))), ]
tmp$tree <- rep(levels(dat$tree), each = nrow(domain)) # assign the trees
simdat <- rbind(pres, tmp)
rm(tmp)

# fit MVLGCP with a variety of basis dimensions and record edf, BIC and predictive performance

# combine the true intensity fields into a vector of quadrature by trees
true_flds <- as.vector(as.matrix(gen_flds))

# initialise some storage
actual_k <- NULL
edf <- NULL
ll <- NULL
aic <- NULL
bic <- NULL
mae <- NULL
sqerr <- NULL
cover <- NULL

ks <- sort(unique(tab$true_k)) # may want to alter this so that this also comes from the job array # ks <- c(25, 50)
for (i in ks) {
  # set the basis functions
  new_bfs <- make_basis(k = i, domain)
  actual_k <- c(actual_k, nrow(new_bfs))
  # fit the model
  tmp.m <- mvlgcp(pt ~ (1|tree) + x:tree, data = simdat, weights = simdat$wt, basis.functions = new_bfs, response.id = simdat$tree)
  # extract the fitted values
  tmp.fitted <- predict(tmp.m)
  # plot(vec2im(tmp.fitted[simdat$tree=="blackoak" & simdat$pt == 0], simdat$x[simdat$tree=="blackoak" & simdat$pt == 0], simdat$y[simdat$tree=="blackoak" & simdat$pt == 0]))
  # points(simdat[simdat$tree == "blackoak" & simdat$pt == 1, c("x", "y")])
  tmp.betas <- getME(tmp.m, "beta")
  
  # store the log-likelihood
  ll <- c(ll, logLik(tmp.m))
  
  # Attempt to calculate the EDF
  try(assign("tmp.inf", influence(tmp.m)))
  if (exists("tmp.inf")) {
    # calculate the EDF
    tmp.edf <- sum(tmp.inf)
    # store results
    edf <- c(edf, tmp.edf)
    aic <- c(aic, -2*as.numeric(logLik(tmp.m)) + tmp.edf * 2)
    bic <- c(bic, -2*as.numeric(logLik(tmp.m)) + tmp.edf * log(sum(simdat$pt == 1)))
    rm(tmp.inf)
  } else {
    # store results
    edf <- c(edf, NA)
    aic <- c(aic, NA)
    bic <- c(bic, NA)
  }
  
  # calculate mean absolute error between the fitted and true intensity fields, across all species
  fitted_flds <- NULL
  for (j in levels(simdat$tree)) {
    fitted_flds <- c(fitted_flds, exp(tmp.fitted)[simdat$pt == 0 & simdat$tree == j])
  }
  mae <- c(mae,  mean(abs(fitted_flds - true_flds)))
  # calculate the squared error in beta estimates
  sqerr <- cbind(sqerr, (tmp.betas - true.beta)^2)
  # calculate the coverage of beta estimates
  ci <- confint(tmp.m)[1:length(tmp.betas), 1:2]
  cover <- cbind(cover,
                 true.beta < ci[,2] & true.beta > ci[,1]
  )
  rm(tmp.m, tmp.fitted, tmp.betas)
}

# par(mfrow = c(5,2), mar = c(0,0,1,0))
# for (i in 1:length(levels(dat$tree))) {
#   plot(vec2im(true_flds[(1:10201) + 10201 * (i - 1)], domain$x, domain$y), main = "")
#   plot(vec2im(fitted_flds[(1:10201) + 10201 * (i - 1)], domain$x, domain$y), main = "")
# }
# par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

# collate the results for storage
res_fixed <- data.frame(t(sqerr), t(cover))
colnames(res_fixed) <- c(paste0("se_beta", 1:(length(levels(dat$tree)) + 1)), paste0("cover_beta", 1:(length(levels(dat$tree)) + 1)))
res <- data.frame(
  set_k = ks, fitted_k = actual_k, edf, ll, aic, bic, mae,
  res_fixed,
  data.frame(tab[tab$job == job, ], actual_k = nrow(bfs))
)

# save results
save(res, file = paste0("results_covar/res_", job, ".RDATA"))