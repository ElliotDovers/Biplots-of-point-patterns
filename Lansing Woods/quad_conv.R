## Convergence of the quadrature approximation #################################

library(glmmTMB)

# get the job number from the array
# job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# extract the regular grid of quadrature included in the glmmTMB data
data("lansing", package = "glmmTMB")
domain <- lansing[lansing$tree == "blackoak" & lansing$pt == 0, ]

# convert the lansing woods data into the required format
data("lansing", package = "spatstat.data")
pp <- as.data.frame(lansing)
colnames(pp)[3] <- "tree"
pp$pt <- 1
pp$wt <- 1e-6

# set the basis functions
bfs <- make_basis(k = 150, domain)

# set the scope for q
tmp.qs <- c(1000, 5000, seq(10000, 100000, by = 30000))

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
  load("quad.conv..RDATA")
} else {
  
  for (i in 1:length(n_q)) {
    
    # create the data
    quad <- data.frame(x = runif(n_q[i]), y = runif(n_q[i]), pt = 0, wt = lansing$window$units$multiplier^2 / n_q[i])
    tmp <- quad[rep(seq_len(n_q[i]), length(levels(pp$tree))), ]
    tmp$tree <- rep(levels(pp$tree), each = n_q[i]) # assign the tree classes
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
    gc()
    print(paste0("Completed rep ", rep.idx[i], " of q = ", n_q[i]))
    quad.conv <- data.frame(
      q = n_q, rep = rep.idx, ll0, ll
    )
    save(quad.conv, file = "working.quad.conv.RDATA")
  }
  quad.conv <- data.frame(
    q = n_q, rep = rep.idx, ll0, ll
  )
  save(quad.conv, file = "quad.conv.RDATA")
}