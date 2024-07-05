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

# THIS APPROACH DOESN'T WORK - INSTEAD GET THE GENERATING FIELDS FROM FITTED VALUES
# flds <- get_field(gen_mod, domain, which.response = 1:5) # THIS IS BROKEN
intercepts <- getME(gen_mod, "b")[1:length(levels(dat$tree))]

call.list <- as.list(gen_mod$call)
# indicate not to fit the model
call.list$doFit <- FALSE
# get the model structure (by re-running the model without fitting)
mod_str <- eval(as.call(call.list)) # NOTE: this gets the corrected Z for smoothers since it evaluates mvlgcp()

# FOR MVLGCP: if the model is from mvlgcp() and has a rr() structure we need to re-order the Z matrix
if (call.list[[1]] == "mvlgcp" & !is.null(gen_mod$col.idx)) {
  # separate out the basis functions from the other random effects
  other_re <- mod_str$data.tmb$Z[ , 1:(ncol(mod_str$data.tmb$Z) - length(mod_str$col.idx))]
  newZ <- mod_str$data.tmb$Z[ , ((ncol(mod_str$data.tmb$Z) - length(mod_str$col.idx)) + 1):ncol(mod_str$data.tmb$Z)]
  # re-order the basis function component according to the column indexing
  mod_str$data.tmb$Z <- cbind(other_re, newZ[ , order(gen_mod$col.idx)])
}

# set the index of the basis function random effects
i_re <- 2

# extract the factor loadings
Lambda <- gen_mod$obj$env$report(gen_mod$fit$parfull)$fact_load[[i_re]]

# latent factors are independent, standard normal so use identity penalty matrix
Sigma <- diag(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps * ncol(Lambda))

# also need to adjust the random effect design matrix by multiplying by the factor loadings

# expand the factor loadings by blockReps to align with Z matrix
IxLambda <- kronecker(diag(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps), Lambda) # NOTE: THIS WAY IMPLIES WE NEED Z TO BE k1:resp1, k1:resp2,...,k2:resp1, k2:resp2

# get the current random effect indices within Z
n_re <- unlist(lapply(gen_mod$modelInfo$reStruc$condReStruc, function(x){ x[["blockReps"]] * x[["blockSize"]] }))
if (i_re == 1) {
  # for the first random effects start at 1
  re_idx <- 1:n_re[[i_re]]
} else {
  # for other random effects start at the previous index + 1 and up to the cumulative sum of random effect indices
  re_idx <- (n_re[[i_re - 1]] + 1):cumsum(n_re)[[i_re]]
}
# compute the new Z component
newZ_comp <- mod_str$data.tmb$Z[ , re_idx] %*% IxLambda
# put the current Z matrix back together again
des.mat <- as.matrix(cbind(mod_str$data.tmb$X,
                     mod_str$data.tmb$Z[ , 1:ncol(mod_str$data.tmb$Z) < min(re_idx)],
                     newZ_comp,
                     mod_str$data.tmb$Z[ , 1:ncol(mod_str$data.tmb$Z) > max(re_idx)])
)

pars <- c(getME(gen_mod, "beta"), intercepts, rnorm(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps * 2))
flds <- des.mat %*% pars

# newU <- c(rnorm(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps),  rnorm(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps))
# fld1 <- tmpZ[ , 1:gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps] %*% rnorm(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps)
# fld2 <- tmpZ[ , (gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps + 1):ncol(tmpZ)] %*% rnorm(gen_mod$modelInfo$reStruc$condReStruc[[i_re]]$blockReps)

par(mfrow = c(5,2), mar = c(0,0,1,0))
for (i in 1:length(levels(dat$tree))) {
  plot(vec2im(fld1[dat$tree == levels(dat$tree)[i] & dat$pt == 0, ], domain$x, domain$y), main = levels(dat$tree)[i])
  plot(vec2im(fld2[dat$tree == levels(dat$tree)[i] & dat$pt == 0, ], domain$x, domain$y), main = levels(dat$tree)[i])
}
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

# initialise some storage
pp.list <- list() # list of each response point pattern
gen_flds <- NULL # data frame of generated fields
par(mfrow = c(2,3))
for (i in 1:length(levels(dat$tree))) {
  # plot(vec2im(predict(gen_mod)[dat$tree == levels(dat$tree)[i] & dat$pt == 0], dat$x[dat$tree == levels(dat$tree)[i] & dat$pt == 0], dat$y[dat$tree == levels(dat$tree)[i] & dat$pt == 0]), main = levels(dat$tree)[i])
  # points(dat[dat$tree == levels(dat$tree)[i] & dat$pt == 1, c("x", "y")])
  # gen_flds <- cbind(gen_flds, exp(predict(gen_mod)[dat$tree == levels(dat$tree)[i] & dat$pt == 0]))
  gen_flds <- cbind(gen_flds, exp(flds[dat$tree == levels(dat$tree)[i] & dat$pt == 0]))
  # eval(parse(text = paste0("domain$fld_", i, " <- exp(getME(gen_mod, 'beta') + intercepts[which(i == levels(dat$tree))] + flds[ , which(i == levels(dat$tree))])")))
  pp.list[[i]] <- spatstat.random::rpoispp(lambda = vec2im(gen_flds[ , i], domain$x, domain$y))
  plot(vec2im(log(gen_flds[,i]), domain$x, domain$y))
  plot(pp.list[[i]], add = T)
}
par(mfrow = c(1,1))

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

# initialise some storage
edf <- NULL
bic <- NULL
mae <- NULL

ks <- sort(unique(tab$true_k)) # may want to alter this so that this also comes from the job array
for (i in ks) {
  # set the basis functions
  new_bfs <- make_basis(k = i, domain)
  # fit the model
  tmp.m <- mvlgcp(pt ~ (1|tree), data = simdat, weights = simdat$wt, basis.functions = new_bfs, response.id = simdat$tree)
  # calculate the EDF
  tmp.edf <- sum(influence(tmp.m))
  # extract the fitted values
  tmp.fitted <- predict(tmp.m)
  # plot(vec2im(tmp.fitted[simdat$tree=="blackoak" & simdat$pt == 0], simdat$x[simdat$tree=="blackoak" & simdat$pt == 0], simdat$y[simdat$tree=="blackoak" & simdat$pt == 0]))
  # points(simdat[simdat$tree == "blackoak" & simdat$pt == 1, c("x", "y")])
  
  # store results
  edf <- c(edf, tmp.edf)
  bic <- c(bic, -2*as.numeric(logLik(tmp.m)) + tmp.edf * log(sum(simdat$pt == 1)))
  # calculate mean absolute error between the fitted and true intensity fields, across all species
  mae <- c(mae,  mean(abs(exp(tmp.fitted)[simdat$pt == 0] - exp(predict(gen_mod))[dat$pt == 0]))) # NOTE: something like KL Div. would require being able to predict from the generating model to the simulated point patterns
}

# collate the results for storage
res <- data.frame(
  fitted_k = ks, edf, bic, mae,
  t(replicate(length(ks), tab[tab$job == job, ], simplify = "data.frame"))
)

# save results
save(res, file = paste0("results/res_", job, ".RDATA"))