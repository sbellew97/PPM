library(raster)
library(spatstat)
library(dismo)
library(ecospat)
library(terra)
library(maptools)
library(ggplot2)
library(colorspace)

# reproducibility
rm(list = ls())
set.seed(2)

# parameters in the script
# desaggregate_threshold <- 10 / 60
desaggregate_threshold <- 0 # larger values mean we remove presence points that are too close to one another
# resolutions <- c(1, 2, 5, 7, 10, 20, 30, 40, 50, 60) / 60 # different resolutions considered (in arcminutes)
npoints <- floor(exp(seq(from = 4, to = 10, by = 0.25))) # number of background points to consider

##########################################
# load presence data with initial analysis
##########################################

# load data in df format
dattot <- read.csv("GBIF_capeweed.csv") 
dattot <- data.frame(lat = dattot$decimalLatitude,
                     lon = dattot$decimalLongitude,
                     year = dattot$year)
# subset to recent years for fitting, earlier years for validation
cutoff_year <- 2012
dat <- dattot[dattot$year > cutoff_year, names(dattot) != "year"]
datval <- dattot[dattot$year <= cutoff_year & dattot$year > 2009, names(dattot) != "year"]
names(dat) <- c("y", "x")
names(datval) <- c("y", "x")
# remove duplicates in both sets
dat <- dat[!duplicated(dat),]
datval <- datval[!duplicated(datval),]
# desaggregate the data in both sets, removing points closer than a given threshold from one another
# in the current version, we do not do any desaggregation as this loses some dimension of the spatial data
dat.thin <- ecospat.occ.desaggregation(dat,
                                       min.dist = desaggregate_threshold)
dat.val.thin <- ecospat.occ.desaggregation(datval,
                                           min.dist = desaggregate_threshold)

#####################
# load covariate data
#####################

# load data
covraster.coarse <- stack("covariates.tif")
names(covraster.coarse) <- c('bio8',
                             'bio14',
                             'bio18',
                             'pH')
# extend covariate data, basically fixing xmax
e <- extent(130, 154, -44, -30) 
covraster.coarse_extended <- terra::extend(covraster.coarse, e)
# plot resulting covariates
plot(covraster.coarse_extended)

############
# set window
############

# contruct a window corresponding to land in SE Australia
pixel <- !is.na(rowSums(values(covraster.coarse_extended)))
pixel <- matrix(pixel,
                nrow = nrow(covraster.coarse_extended),
                ncol = ncol(covraster.coarse_extended),
                byrow = TRUE)
pixel <- pixel[c(nrow(pixel):1), ]
victoriapixel <- owin(xrange = c(130, 154),
                      yrange = c(-44, -30),
                      mask = pixel)

# convert covariates to spatstat format
bio8 <- maptools::as.im.RasterLayer(covraster.coarse_extended$bio8)
bio14 <- maptools::as.im.RasterLayer(covraster.coarse_extended$bio14)
bio18 <- maptools::as.im.RasterLayer(covraster.coarse_extended$bio18)
pH <- maptools::as.im.RasterLayer(covraster.coarse_extended$pH)
covlist <- list(bio8 = bio8,
                bio14 = bio14,
                bio18 = bio18,
                pH = pH)

# Remove presence points which are not inside the window
dat.thin <- dat.thin[inside.owin(dat.thin$x, dat.thin$y, victoriapixel), ]
dat.val.thin <- dat.val.thin[inside.owin(dat.val.thin$x, dat.val.thin$y, victoriapixel), ]

# convert (thinned) presence points to ppp format
configuration <- ppp(x = dat.thin$x, 
                     y = dat.thin$y, 
                     window = victoriapixel)
configuration_val <- ppp(x = dat.val.thin$x, 
                         y = dat.val.thin$y, 
                         window = victoriapixel)

result <- lapply(npoints, function(n) {
  ################################################################################################################
  # selecting grid points from total background raster (grid resolution = resolution(covraster.coarse)*gridfactor)
  ################################################################################################################
  windowarea <- volume(configuration$window)
  res <- sqrt(windowarea / n)
  cov <- values(covraster.coarse_extended)
  # gridfactor controls how many of the background points we keep
  gridfactor <- res / mean(res(covraster.coarse_extended))
  print(gridfactor)
  row <- ceiling(gridfactor*c((gridfactor + 1)/(2 * gridfactor) + 0:(nrow(covraster.coarse_extended) / gridfactor - 1)))
  col <- ceiling(gridfactor*c((gridfactor + 1)/(2 * gridfactor) + 0:(ncol(covraster.coarse_extended) / gridfactor - 1)))
  cellnum <- terra::cellFromRowColCombine(covraster.coarse_extended, row, col)
  # sgrid.loc contains subsetted background points
  sgrid.loc <- xyFromCell(covraster.coarse_extended,
                          cellnum)
  # sgrid.loctot contains all background points
  sgrid.loctot <- coordinates(covraster.coarse_extended)
  
  # plot the window
  # plot(victoriapixel)
  # compute the background points in ppp format, for both the subsetted set and the full one
  sgrid.loc.ipp <- ppp(sgrid.loc[, 'x'],
                       sgrid.loc[, 'y'],
                       window = victoriapixel)
  sgrid.loc.ipptot <- ppp(sgrid.loctot[, 'x'],
                          sgrid.loctot[, 'y'],
                          window = victoriapixel)
  # plot the subsetted data
  # plot(sgrid.loc.ipp)
  
  ###########################
  # prepare data for spatstat
  ###########################
  
  # convert the background ppp objects along with presence points to quadscheme format used with ppm
  datadummylogi <- quadscheme.logi(configuration,
                                   sgrid.loc.ipp)
  datadummyppm <- quadscheme(configuration,
                             sgrid.loc.ipp,
                             method = "grid")
  datadummylogitot <- quadscheme.logi(configuration,
                                      sgrid.loc.ipptot)
  datadummyppmtot <- quadscheme(configuration,
                                sgrid.loc.ipptot,
                                method = "grid")
  
  ###########
  # fit model
  ###########
  
  capeweedppm <- ppm(datadummyppm,
                     trend = ~ 1 + bio8 + bio14 + bio18 + pH + bio8 : pH,
                     covariates = covlist,
                     method = 'mpl')
  capeweedlogi <- ppm(datadummylogi,
                      trend = ~ 1 + bio8 + bio14 + bio18 + pH + bio8 : pH,
                      covariates = covlist,
                      method = 'logi')
  # plot predictions
  plot.ppm(capeweedppm,
           covariates = covlist,
           locations = victoriapixel,
           superimpose = TRUE,
           se = FALSE,
           col = rainbow)
  
  plot.ppm(capeweedlogi,
           covariates = covlist,
           locations = victoriapixel,
           se = FALSE,
           superimpose = TRUE,
           col = rainbow)
  
  #############################
  # log likelihood calculations
  #############################
  
  # logistic log-likelihood
  likelihoodblr <- function(param) {
    dum <- nrow(X.back1tot) / windowarea
    sum(log(exp(X.po %*% param) / (exp(X.po %*% param) + dum)), na.rm = TRUE) + 
      sum(log(dum / (exp(X.back1tot %*% param) + dum)), na.rm = TRUE)
    
  }
  
  # Poisson log-likelihood
  likelihoodIPPBTtot = function(param) {
    lambda = exp(X.backtot %*% param)
    mu = lambda * area.backtot
    logL.pp = sum(X.po %*% param, na.rm = TRUE) - sum(mu, na.rm = TRUE)
    
    (-1) * sum(logL.pp)
  }
  
  # add cross covariate to the raster stack
  bio8pH <- raster(covraster.coarse_extended)
  values(bio8pH) <- values(covraster.coarse_extended)[, 'bio8'] * values(covraster.coarse_extended)[, 'pH']
  covraster.coarse_extended_plus_cross <- stack(covraster.coarse_extended,
                                                bio8pH)
  names(covraster.coarse_extended_plus_cross)[names(covraster.coarse_extended_plus_cross) == "layer"] <- "bio8.pH"
  # compute the values of the covariates at both presences and presences + absences
  X.po <- data.frame(x = datadummyppmtot$data$x,
                     y = datadummyppmtot$data$y)
  X.backtot <- rbind(data.frame(x = datadummyppmtot$data$x,
                                y = datadummyppmtot$data$x),
                     data.frame(x = datadummyppmtot$dummy$x,
                                y = datadummyppmtot$dummy$y))
  X.po <- cbind(intercept = 1, extract(covraster.coarse_extended_plus_cross,
                                       X.po))
  X.backtot <- cbind(intercept = 1, extract(covraster.coarse_extended_plus_cross,
                                            X.backtot))
  area.backtot <- w.quad(datadummyppmtot)
  # Nb: I think commented code is incorrect, the actual window is smaller
  # windowarea <- (xmax(covraster.coarse_extended)-xmin(covraster.coarse_extended))*(ymax(covraster.coarse_extended)-ymin(covraster.coarse_extended))
  X.back1tot <- cbind(x = datadummylogitot$dummy$x,
                      y = datadummylogitot$dummy$y)
  X.back1tot <- cbind(intercept = 1, extract(covraster.coarse_extended_plus_cross,
                                             X.back1tot))
  
  maxllfitIPP <- -likelihoodIPPBTtot(as.vector(coefficients(capeweedppm))) - sum(log(1:configuration$n))
  maxllblr <- likelihoodblr(as.vector(coefficients(capeweedlogi)))
  
  ######################
  # validation using AUC
  ######################
  
  predicted_intensity <- predict.ppm(capeweedppm, type = "intensity", ngrid = dim(covraster.coarse_extended)[c(1,2)])
  resppm <- spatstat.core::auc(X = configuration_val, covariate = predicted_intensity)
  predicted_intensity <- predict.ppm(capeweedlogi, type = "intensity", ngrid = dim(covraster.coarse_extended)[c(1,2)])
  reslogi <- spatstat.core::auc(X = configuration_val, covariate = predicted_intensity)
  
  list(auclogi = reslogi,
       aucppm = resppm,
       lllogi = maxllblr,
       llpoisson = maxllfitIPP,
       fitppm = capeweedppm,
       fitlogi = capeweedlogi,
       npoints = length(cellnum),
       resolution = res)
})

df <- data.frame(minutes = sapply(result, function(x) x$resolution) * 60,
                 npoints = sapply(result, function(x) x$npoints),
                 auclogi = sapply(result, function(x) x$auclogi),
                 aucppm = sapply(result, function(x) x$aucppm),
                 lllogi = sapply(result, function(x) x$lllogi),
                 llpoisson = sapply(result, function(x) x$llpoisson))

ggplot(data = df) + 
  geom_line(aes(x = log(npoints), y = auclogi, colour = "BLR")) +
  geom_line(aes(x = log(npoints), y = aucppm, colour = "PPM-BT")) +
  geom_point(aes(x = log(npoints), y = auclogi, colour = "BLR")) +
  geom_point(aes(x = log(npoints), y = aucppm, colour = "PPM-BT")) +
  ylab("AUC") + 
  xlab("Logarithm of the number of background points") +
  scale_color_discrete(name = "") +
  theme_minimal()

ggplot(data = df) + 
  geom_line(aes(x = log(npoints), y = lllogi, colour = "BLR")) +
  geom_line(aes(x = log(npoints), y = llpoisson, colour = "PPM-BT")) +
  geom_point(aes(x = log(npoints), y = lllogi, colour = "BLR")) +
  geom_point(aes(x = log(npoints), y = llpoisson, colour = "PPM-BT")) +
  ylab("Log-likelihood") + 
  xlab("Logarithm of the number of background points") +
  scale_color_discrete(name = "") +
  theme_minimal()

# best fit
fit <- result[[which.max(sapply(result, function(x) x$npoints))]]$fitppm
# summary
summary(fit)
# number of dummy points
length(fit$Q$dummy$x)

# plot the fit
predicted_intensity <- as.data.frame(rasterToPoints(raster(log(predict.ppm(fit, type = "intensity", ngrid = dim(covraster.coarse_extended)[c(1,2)])))))
a <- "Training"
b <- "Validation"

names(predicted_intensity)[names(predicted_intensity) == 'layer'] <- "Log-intensity"

ggplot(data = predicted_intensity, aes_string(x = 'x', y = 'y')) +
  geom_tile(aes(fill = `Log-intensity`), alpha = 0.5) +
  scale_fill_continuous_sequential(palette = "Purple-Yellow") +
  geom_point(data = as.data.frame(configuration), aes(x = x, y = y, color = a, shape = a), size = 2, alpha = 0.1) +
  geom_point(data = as.data.frame(configuration_val), aes(x = x, y = y, color = b, shape = b), size = 2, alpha = 0.1) +
  scale_color_manual(values = c("black", "orange")) + # Obtained by wesanderson::wes_palette("Darjeeling1")
  scale_shape_manual(values = c(16, 17)) +
  xlab(NULL) + # Remove x labels
  ylab(NULL) + # Remove y labels
  coord_equal() +
  theme_minimal(base_size = 15) +
  guides(shape = guide_legend(override.aes = list(alpha = 1, size = 3), title = ""),
         color = guide_legend(alpha = 1, title = ""))
