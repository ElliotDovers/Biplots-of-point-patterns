# Get the job array
load("job_array.RDATA")

# determine job number from pbs script
job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

## MODEL SETTINGS ##############################################################

# set the number of latent factor we want to use
d = tab$d[tab$job == job]
# set approximate number of basis function nodes to use
basis.dim = tab$k[tab$job == job]
# get the predictor set to be modelled
# pred_set <- as.character(tab$p[tab$job == job])
# get the basis function type to be used
# bf.type <- as.character(tab$bf.type[tab$job == job])
################################################################################

library(glmmTMB)

# get the data
source("get_and_format_flora_data.R")
rm(dat_pa)

# get the relevant predictors according to pred_set (these are all basically derived from elevation, we are just looking at moisture and temp. covars and soil fertility score and disturbance index)
# preds <- switch(pred_set,
#                 A = "None",
#                 B = c("mi", "tempann"),
#                 C = c("mi", "tempann", "mi2", "tempann2"),
# )

# get a vector of the factor terms
facts <- c("disturb", "soilfert", "vegsis")

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

# assign some numbers
m <- length(species) # responses
q <- nrow(quad) # quadrature points

# add in the quadrature for each species (replicates of the supplied background points)
tmp <- quad[rep(seq_len(q), m), ]
tmp$spid <- rep(species, each = q) # assign the species id's
dat <- rbind(pres, tmp)
n <- nrow(dat)
rm(tmp)

# set the coordinates to x and y so we can just use default labels
dat$lon <- dat$x; dat$lat <- dat$y; dat$x <- dat$xm; dat$y <- dat$ym; dat$xm <- NULL; dat$ym <- NULL
domain$lon <- domain$x; domain$lat <- domain$y; domain$x <- domain$xm; domain$y <- domain$ym; domain$xm <- NULL; domain$ym <- NULL

# set the basis functions

bfs <- make_basis(k = basis.dim, domain, from.package = "scampr")

# need to prune in the case of scampr bfs
if (is(bfs, "scampr.bf")) {
  tmp <- bf_matrix(bfs, rbind(pres[,c("xm", "ym")], quad[,c("xm", "ym")])) # via the data?
  # tmp <- bf_matrix(bfs, pres[,c("xm", "ym")]) # via the data?
  bf.check <- apply(tmp, 2, function(x){sum(x > 0)})
  bfs <- bfs[bf.check != 0, ]
  rm(tmp, bf.check)
}
rm(pres, quad)
gc()

################################################################################
# Fit an intercepts-only model #################################################
################################################################################

# need a function to strip environments from the model objects so they can be saved without being massive
rm.env <- function(m) {
  attr(attr(m$frame, "terms"), ".Environment") <- NULL
  attr(m$modelInfo$terms$cond$fixed, ".Environment") <- NULL
  attr(m$modelInfo$reTrms$cond$terms$fixed, ".Environment") <- NULL
  attr(m$modelInfo$allForm$combForm, ".Environment") <- NULL
  attr(m$modelInfo$allForm$formula, ".Environment") <- NULL
  attr(m$modelInfo$allForm$ziformula, ".Environment") <- NULL
  attr(m$modelInfo$allForm$dispformula, ".Environment") <- NULL
  m
}

# single correlated fields #####################################################

try(assign("mA1", mvlgcp(occ ~ (1|spid), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
              control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))

# save if the model fitted
if (exists("mA1")) {
  print("mA1 converged without start parameters (attempt 1)")
  mA1 <- rm.env(mA1)
  mA1$start.pars <- "cold"
  save(mA1, file = paste0(getwd(), "/models/mA1_", job, ".RDATA"))
}

# for the dual fields correlated and species-specific ##########################

# first try to fit with cold starts for parameters
try(assign("mA2", mvlgcp(occ ~ (1|spid), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both")
))
if (exists("mA2")) {
  print(paste0("mA2 converged without start parameters (attempt 1)", if(is.na(logLik(mA2))) {" however, the likelihood was NA"} ))
  mA2$start.pars <- "cold"
  if (is.na(logLik(mA2))) {
    rm(mA2)
  }
}

# next try with mA1 as starts for parameters
if (exists("mA1") & !exists("mA2")) {
  start.pars <- lapply(split(mA1$fit$parfull, names(mA1$fit$parfull)), unname)
  start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * m)) # add additional field coefs
  start.pars$theta <- c(start.pars$theta, rep(-3, m)) # add additional field variances
  try(assign("mA2", mvlgcp(occ ~ (1|spid), data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both",
                start = start.pars)
  ))
  rm(start.pars)
  # check that the model has a log-likelihood
  if (exists("mA2")) {
    print(paste0("mA2 fitted with starting parameters from mA1 (attempt 2)", if(is.na(logLik(mA2))) {" however, the likelihood was NA"} ))
    mA2$start.pars <- "warm"
    if (is.na(logLik(mA2))) {
      rm(mA2)
    }
  }
}

# save if the model fitted
if (exists("mA2")) {
  mA2 <- rm.env(mA2)
  save(mA2, file = paste0(getwd(), "/models/mA2_", job, ".RDATA"))
}

################################################################################
# Fit the model including linear env terms #####################################
################################################################################

# with correlated fields only ##################################################

try(assign("mB1", mvlgcp(occ ~ (1 | spid) + mi + tempann + (0 + mi | spid) + (0 + tempann | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))

# save if the model fitted
if (exists("mB1")) {
  print("mB1 converged without start parameters (attempt 1)")
  mB1 <- rm.env(mB1)
  mB1$start.pars <- "cold"
  save(mB1, file = paste0(getwd(), "/models/mB1_", job, ".RDATA"))
}

# with dual fields #############################################################

# first try to fit with cold starts for parameters
try(assign("mB2", mvlgcp(occ ~ (1 | spid) + mi + tempann + (0 + mi | spid) + (0 + tempann | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both")
))
if (exists("mB2")) {
  print(paste0("mB2 converged without start parameters (attempt 1)", if(is.na(logLik(mB2))) {" however, the likelihood was NA"} ))
  mB2$start.pars <- "cold"
  if (is.na(logLik(mB2))) {
    rm(mB2)
  }
}

# next try with mB1 as starts for parameters
if (exists("mB1") & !exists("mB2")) {
  start.pars <- lapply(split(mB1$fit$parfull, names(mB1$fit$parfull)), unname)
  start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * m))
  start.pars$theta <- c(start.pars$theta, rep(-3, m))
  try(assign("mB2", mvlgcp(occ ~ (1 | spid) + mi + tempann + (0 + mi | spid) + (0 + tempann | spid),
                           data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                           control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both",
                           start = start.pars)
  ))
  rm(start.pars)
  # check that the model has a log-likelihood
  if (exists("mB2")) {
    print(paste0("mB2 fitted with starting parameters from mB1 (attempt 2)", if(is.na(logLik(mB2))) {" however, the likelihood was NA"} ))
    mB2$start.pars <- "warm"
    # if (is.na(logLik(mB2))) {
    #   rm(mB2)
    # }
  }
}

# save if the model fitted
if (exists("mB2")) {
  mB2 <- rm.env(mB2)
  save(mB2, file = paste0(getwd(), "/models/mB2_", job, ".RDATA"))
}

################################################################################
# Fit the model including quadratic env. terms #################################
################################################################################

# with correlated fields only ##################################################

try(assign("mC1", mvlgcp(occ ~ (1 | spid) + mi + tempann + mi2 + tempann2 + (0 + mi |  spid) +
                           (0 + tempann | spid) + (0 + mi2 | spid) + (0 + tempann2 | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "correlated")
))

# save if the model fitted
if (exists("mC1")) {
  mC1 <- rm.env(mC1)
  mC1$start.pars <- "cold"
  print("mC1 converged without start parameters (attempt 1)")
  save(mC1, file = paste0(getwd(), "/models/mC1_", job, ".RDATA"))
}

# with dual fields #############################################################

# first try to fit with cold starts for parameters
try(assign("mC2", mvlgcp(occ ~ (1 | spid) + mi + tempann + mi2 + tempann2 + (0 + mi |  spid) +
                           (0 + tempann | spid) + (0 + mi2 | spid) + (0 + tempann2 | spid),
                         data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                         control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both")
))
if (exists("mC2")) {
  print(paste0("mC2 converged without start parameters (attempt 1)", if(is.na(logLik(mC2))) {" however, the likelihood was NA"} ))
  mC2$start.pars <- "cold"
  if (is.na(logLik(mC2))) {
    rm(mC2)
  }
}

# next try with mC1 as starts for parameters
if (exists("mC1") & !exists("mC1")) {
  start.pars <- lapply(split(mC1$fit$parfull, names(mC1$fit$parfull)), unname)
  start.pars$b <- c(start.pars$b, rep(0, nrow(bfs) * m))
  start.pars$theta <- c(start.pars$theta, rep(-3, m))
  try(assign("mC2", mvlgcp(occ ~ (1 | spid) + mi + tempann + mi2 + tempann2 + (0 + mi |  spid) +
                             (0 + tempann | spid) + (0 + mi2 | spid) + (0 + tempann2 | spid),
                           data = dat, weights = dat$wt, basis.functions = bfs, response.id = dat$spid, n_factors = d,
                           control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)), mv.flds = "both",
                           start = start.pars)
  ))
  rm(start.pars)
  # check that the model has a log-likelihood
  if (exists("mC2")) {
    print(paste0("mC2 fitted with starting parameters from mC1 (attempt 2)", if(is.na(logLik(mC2))) {" however, the likelihood was NA"} ))
    mC2$start.pars <- "warm"
    if (is.na(logLik(mC2))) {
      rm(mC2)
    }
  }
}

# save if the model fitted
if (exists("mC2")) {
  mC2 <- rm.env(mC2)
  save(mC2, file = paste0(getwd(), "/models/mC2_", job, ".RDATA"))
}

