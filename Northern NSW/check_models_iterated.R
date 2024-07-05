home.wd <- getwd()

library(glmmTMB)

# Get the job array
load("job_array.RDATA")

# determine job number from pbs script
job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# data and more needed to be able to re-eval the model calls ###################
# get the data
source("get_and_format_flora_data.R")
rm(dat_pa)

# load the species information
load("flora_species.RDATA")

# get a vector of the factor terms
facts <- c("disturb", "soilfert", "vegsis")

# get the groups to which the flora species belongs to subset the data
flora_groups <- unique(sp.info$group)

# subset the presence-only data (in the case we are using restricted species)
pres <- pres[pres$group %in% flora_groups, ]

# ## TEMPORARY ##
# # reduce the species to the most prevalent
# pres <- pres[pres$spid %in% names(sort(table(pres$spid), decreasing = T))[1:6], ]

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

################################################################################

# initialise some storage
res <- tab[tab$job == job, ]
res <- cbind(res, tmp = 1:6)
res$p <- rep(c("A", "B", "C"), each = 2)
res$flds <- rep(1:2, times = 3)
res$ll <- NA
res$bfs <- NA
res$bf_rad <- NA
# res$aic <- NA
# res$bic <- NA
# res$full.edf <- NA
# res$fld1.edf <- NA
# res$fld2.edf <- NA
res$msg <- NA
res$starts <- NA

# jump into the folder storing saved models and load and analyse
setwd(paste0(home.wd, "/models"))
print(paste0("checking for job ", job))
# check the fitted model exists
for (i in 1:6) {
  tmp.mask <- res$tmp == i
  if (file.exists(paste0("m", res$p[tmp.mask], res$flds[tmp.mask], "_", job, ".RDATA"))) {
    print(paste0("found and working on job ", job, " model m", res$p[tmp.mask], res$flds[tmp.mask]))
    
    # load the model object
    try(load(paste0("m", res$p[tmp.mask], res$flds[tmp.mask], "_", job, ".RDATA")))
    if (exists(paste0("m", res$p[tmp.mask], res$flds[tmp.mask]))) {
      eval(parse(text = paste0("mod <- m", res$p[tmp.mask], res$flds[tmp.mask])))
      
      # adjust the call so that various objects from the model are found
      mod$call$basis.functions <- mod$basis.functions
      mod$call$n_factors <- mod$n_factors
      
      # determine the number of random effects in the model
      n_re <- unlist(lapply(mod$modelInfo$reStruc$condReStruc, function(x){ x[["blockReps"]] * x[["blockSize"]] }))
      fld1_id <- which(grepl("basis.functions", names(mod$modelInfo$reStruc$condReStruc), fixed = T) & sapply(mod$modelInfo$reStruc$condReStruc, function(x){x[4]}) == 9)
      fld2_id <- which(grepl("basis.functions", names(mod$modelInfo$reStruc$condReStruc), fixed = T) & sapply(mod$modelInfo$reStruc$condReStruc, function(x){x[4]}) != 9)
      # determine the correlated fields indices
      if (fld1_id == 1) {
        # for the first random effects start at 1
        fld1_idx <- 1:n_re[[fld1_id]]
      } else {
        # for other random effects start at the previous index + 1 and up to the cumulative sum of random effect indices
        fld1_idx <- (cumsum(n_re)[[fld1_id - 1]] + 1):cumsum(n_re)[[fld1_id]]
      }
      if (length(fld2_id) > 0) {
        if (fld2_id == 1) {
          # for the first random effects start at 1
          fld2_idx <- 1:n_re[[fld2_id]]
        } else {
          # for other random effects start at the previous index + 1 and up to the cumulative sum of random effect indices
          fld2_idx <- (cumsum(n_re)[[fld2_id - 1]] + 1):cumsum(n_re)[[fld2_id]]
        }
      }
      
      # record model information #################################################
      res$starts <- mod$start.pars
      res$msg[tmp.mask] <- mod$fit$message
      res$ll[tmp.mask] <- as.numeric(logLik(mod))
      res$bfs[tmp.mask] <- nrow(mod$basis.functions)
      res$bf_rad[tmp.mask] <- mod$basis.functions$scale[1]
      # # influence matrix can fail so will wrap a try()
      # try(assign("tmp.inf", influence(mod)))
      # if (exists("tmp.inf")) {
      #   res$full.edf[tmp.mask] <- sum(tmp.inf)
      #   res$fld1.edf[tmp.mask] <- sum(tmp.inf[fld1_idx])
      #   if (length(fld2_id) > 0) {
      #     res$fld2.edf[tmp.mask] <- sum(tmp.inf[fld2_idx])
      #     rm(fld2_idx)
      #   }
      #   res$aic[tmp.mask] <- (-2*as.numeric(logLik(mod))) + (sum(tmp.inf) * 2)
      #   res$bic[tmp.mask] <- (-2*as.numeric(logLik(mod))) + (sum(tmp.inf) * log(sum(dat$occ == 1)))
      #   rm(tmp.inf)
      # }
      print(paste0("job ", job, " model m", res$p[tmp.mask], res$flds[tmp.mask], " completed. Loaded model info."))
      eval(parse(text = paste0("rm(mod, m", res$p[tmp.mask], res$flds[tmp.mask], ", n_re, fld1_id, fld2_id, fld1_idx)")))
      # rm(mod)
    } else {
      print(paste0("job ", job, " model m", res$p[tmp.mask], res$flds[tmp.mask], " had a load error"))
    }
  } else {
    print(paste0("job ", job, " model m", res$p[tmp.mask], res$flds[tmp.mask], " is missing"))
  }
  rm(tmp.mask)
}

setwd(paste0(home.wd, "/model_checks"))
save(res, file = paste0("tab_", job, ".RDATA"))