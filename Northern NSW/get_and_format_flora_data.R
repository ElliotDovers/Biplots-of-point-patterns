library(disdat)
library(sp)

# set the region required
region = "NSW"
region_size <- 76.18 * 1000 # in km^2 (from Elith, Jane, et al. "Presence-only and presence-absence data for comparing species distribution modeling methods." Biodiversity informatics 15.2 (2020): 69-80)
# these are the categorical variables included:
categoricalvars <- c("vegsys", "disturb", "soilfert")
# these are the flora groups
flora_groups <- c("ot", "ou", "rt", "ru") # open-forest trees, open-forest understorey plants, rainforest trees, rainforest understorey plants

# get the relevant predictors
preds <- disPredictors(region)
preds <- preds[-13] # remove vegsys as this factor is likely correlated with the response (as it broadly categorises the main plant types in the area)

# reading presence-only and background species data for this region, one file per region:
pres <- disPo(region)
quad <- disBg(region)

# there is one clear outlying record of solar radiation that I will remove
quad <- quad[quad$solrad != 255, ]

# obtain the presence/absence data in long format (to use as test data), note there are different sites depending on the animal group
dat_pa <- NULL
for (grp in flora_groups) {
  tmp.env <- disEnv(region, grp)
  tmp.occ <- disPa(region, grp)[-(1:4)]
  tmp <- cbind(tmp.env, tmp.occ)
  tmp.dat <- reshape::melt(tmp, id.vars = colnames(tmp.env), measure.vars = colnames(tmp.occ), variable_name = "spid")
  dat_pa <- rbind(dat_pa, tmp.dat)
  rm(tmp.env, tmp.occ, tmp, tmp.dat)
}
colnames(dat_pa)[colnames(dat_pa) == "value"] <- "occ"

# get the species
species <- unique(pres$spid)

# set the quadrature weights
pres$wt <- 1e-6
quad$wt <- region_size / nrow(quad)

# load the full domain data
load("domain.RDATA")

# scale and center covariates in the PO data and quad points (as well as setting factors for categorical variables)
for(i in preds){
  if(i %in% categoricalvars){
    fac_col <- i
    if (fac_col == "soilfert") { # combining the soil fertility ratings beyond 3 (since these are otherwise not well represented)
      pres[ ,fac_col][pres[ ,fac_col] %in% c(4,5)] <- 3
      quad[ ,fac_col][quad[ ,fac_col] %in% c(4,5)] <- 3
      dat_pa[ ,fac_col][dat_pa[ ,fac_col] %in% c(4,5)] <- 3
      domain[ ,fac_col][domain[ ,fac_col] %in% c(4,5)] <- 3
    }
    domain[ ,fac_col] <- as.factor(domain[ ,fac_col])
    pres[ ,fac_col] <- as.factor(pres[ ,fac_col])
    quad[ ,fac_col] <- as.factor(quad[ ,fac_col])
    dat_pa[ ,fac_col] <- as.factor(dat_pa[ ,fac_col])
    # expand and relevel incase the PO or PA datasets are missing levels from the full domain
    pres[ ,fac_col] <- forcats::fct_expand(pres[,fac_col], levels(domain[,fac_col]))
    pres[ ,fac_col] <- forcats::fct_relevel(pres[,fac_col], levels(domain[,fac_col]))
    quad[ ,fac_col] <- forcats::fct_expand(quad[,fac_col], levels(domain[,fac_col]))
    quad[ ,fac_col] <- forcats::fct_relevel(quad[,fac_col], levels(domain[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_expand(dat_pa[,fac_col], levels(domain[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_relevel(dat_pa[,fac_col], levels(domain[,fac_col]))
  } else {
    # scale according to full domain means/sds (on the original "og" covariates)
    pres[ , i] <- as.vector(scale(pres[ , i], center = mean(domain[,paste0(i, "_og")]), scale = sd(domain[,paste0(i, "_og")])))
    quad[ , i] <- as.vector(scale(quad[ , i], center = mean(domain[,paste0(i, "_og")]), scale = sd(domain[,paste0(i, "_og")])))
    dat_pa[ , i] <- as.vector(scale(dat_pa[ , i], center = mean(domain[,paste0(i, "_og")]), scale = sd(domain[,paste0(i, "_og")])))
    # add in orthogonal second order terms
    eval(parse(text = paste0("pres$", i, "1 <- poly(as.matrix(pres[,i]), 2)[,1]")))
    eval(parse(text = paste0("quad$", i, "1 <- poly(as.matrix(quad[,i]), 2)[,1]")))
    eval(parse(text = paste0("dat_pa$", i, "1 <- poly(as.matrix(dat_pa[,i]), 2)[,1]")))
    eval(parse(text = paste0("domain$", i, "1 <- poly(as.matrix(domain[,i]), 2)[,1]")))
    eval(parse(text = paste0("pres$", i, "2 <- poly(as.matrix(pres[,i]), 2)[,2]")))
    eval(parse(text = paste0("quad$", i, "2 <- poly(as.matrix(quad[,i]), 2)[,2]")))
    eval(parse(text = paste0("dat_pa$", i, "2 <- poly(as.matrix(dat_pa[,i]), 2)[,2]")))
    eval(parse(text = paste0("domain$", i, "2 <- poly(as.matrix(domain[,i]), 2)[,2]")))
  }
}

# convert coordinates to UTM ###################################################

# set the domain data as a spatial points dataframe
coordinates(pres) <- c("x", "y")
coordinates(quad) <- c("x", "y")
coordinates(dat_pa) <- c("x", "y")
# set the current projection in long/lat
proj4string(pres) <- "+proj=longlat +datum=WGS84 +no_defs"
proj4string(quad) <- "+proj=longlat +datum=WGS84 +no_defs"
proj4string(dat_pa) <- "+proj=longlat +datum=WGS84 +no_defs"
# project onto UTM
tmp_pres <- as.data.frame(spTransform(pres, CRS("+proj=utm +zone=56 +south ellps=WGS84"))) # note zone 56 (south) for northern NSW
tmp_quad <- as.data.frame(spTransform(quad, CRS("+proj=utm +zone=56 +south ellps=WGS84"))) # note zone 56 (south) for northern NSW
tmp_dat_pa <- as.data.frame(spTransform(dat_pa, CRS("+proj=utm +zone=56 +south ellps=WGS84"))) # note zone 56 (south) for northern NSW
# convert these to 1km resolution (to match quadrature weight scale and otherwise too fine a scale for images etc.)
if (is.null(tmp_pres$x)) { # adjust for different labelling that occurs between sp package versions
  tmp_pres$x <- round(tmp_pres$coords.x1 / 1000)
  tmp_pres$y <- round(tmp_pres$coords.x2 / 1000)
  tmp_quad$x <- round(tmp_quad$coords.x1 / 1000)
  tmp_quad$y <- round(tmp_quad$coords.x2 / 1000)
  tmp_dat_pa$x <- round(tmp_dat_pa$coords.x1 / 1000)
  tmp_dat_pa$y <- round(tmp_dat_pa$coords.x2 / 1000)
} else {
  tmp_pres$x <- round(tmp_pres$x / 1000)
  tmp_pres$y <- round(tmp_pres$y / 1000)
  tmp_quad$x <- round(tmp_quad$x / 1000)
  tmp_quad$y <- round(tmp_quad$y / 1000)
  tmp_dat_pa$x <- round(tmp_dat_pa$x / 1000)
  tmp_dat_pa$y <- round(tmp_dat_pa$y / 1000)
}

# convert the domain data back into a data frame
pres <- as.data.frame(pres)
quad <- as.data.frame(quad)
dat_pa <- as.data.frame(dat_pa)
# add in the UTM locations
pres$xm <- tmp_pres$x
pres$ym <- tmp_pres$y
quad$xm <- tmp_quad$x
quad$ym <- tmp_quad$y
dat_pa$xm <- tmp_dat_pa$x
dat_pa$ym <- tmp_dat_pa$y
rm(tmp_pres, tmp_quad, tmp_dat_pa)