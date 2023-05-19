library(raster)
library(spatstat)
library(dismo)
library(ecospat)
library(terra)
library(maptools)
library(ggplot2)
library(colorspace)
library(Metrics)

# reproducibility
rm(list = ls())
set.seed(2)

bioclim <- raster::getData('worldclim', 
                           var = 'bio',
                           res = 0.5,
                           lon = 115,
                           lat = -34)

e <- extent(116, 119, -34, -31)
bioclim <- crop(bioclim, e)
bioclim <- bioclim[[c(3, 8, 9, 13, 14)]]
bio3 <- subset(bioclim, 1, drop = TRUE)
bio8 <- subset(bioclim, 2, drop = TRUE)
bio9 <- subset(bioclim, 3, drop = TRUE)
bio13 <- subset(bioclim, 4, drop = TRUE)
bio14 <- subset(bioclim, 5, drop = TRUE)

s <- raster(ncol = ncol(bioclim), 
            nrow = nrow(bioclim), 
            xmn = xmin(bioclim), 
            xmx = xmax(bioclim), 
            ymn = ymin(bioclim), 
            ymx = ymax(bioclim))
r <- s
values(r) <- values(bio3)
names(r) <- 'bio3'
s <- addLayer(s, r)
values(r) <- values(bio8)
names(r) <- 'bio8'
s <- addLayer(s, r)
# values(r) <- values(bio9)
# names(r) <- 'bio9'
# s <- addLayer(s, r)
# values(r) <- values(bio13)
# names(r) <- 'bio13'
# s <- addLayer(s, r)
# values(r) <- values(bio14)
# names(r) <- 'bio14'
# s <- addLayer(s, r)
s <- scale(s, center = TRUE, scale = TRUE)

#create covariate list for use in spatstat from the covariate raster

bio3 <- scale(bio3, center = TRUE, scale = TRUE)
bio8 <- scale(bio8, center = TRUE,scale  =TRUE)
bio9 <- scale(bio9, center = TRUE, scale = TRUE)
bio13 <- scale(bio13, center = TRUE, scale = TRUE)
bio14 <- scale(bio14, center = TRUE, scale = TRUE)

covlist <- list(a = as.im.RasterLayer(bio3),
                b = as.im.RasterLayer(bio8))#,
                # c = as.im.RasterLayer(bio9),
                # d = as.im.RasterLayer(bio13),
                # e = as.im.RasterLayer(bio14))
w <- as.owin(covlist[[1]])

# How many presence points do we want to consider?
npresence <- readline(prompt = "How many presence points do you want to consider? Input a single character, l for ~ 100, m for ~ 1000, and h for ~10000 ")

X <- cbind(1,
           values(s)[, 'bio3'],
           values(s)[, 'bio8'])#,
           # values(s)[, 'bio9'],
           # values(s)[, 'bio13'],
           # values(s)[, 'bio14'])
gridsizeratio = sum(!is.na(X)) / length(X)
if(npresence == "l") {
  parameter <- c(beta1 = 5 - log(1.4) + 1 * log(0.1), beta2 = 0.2, beta3 = 0.3)#, beta4 = -0.2, beta5 = 0.26, beta6 = -0.2)
} else if(npresence == "m") {
  parameter <- c(beta1 = 5 - log(1.4) - 0 * log(0.1), beta2 = 0.2, beta3 = 0.3)#, beta4 = -0.2, beta5 = 0.26, beta6 = -0.2)
} else if(npresence == "h") {
  parameter <- c(beta1 = 5 - log(1.4) - 1 * log(0.1), beta2 = 0.2, beta3 = 0.3)#, beta4 = -0.2, beta5 = 0.26, beta6 = -0.2)
} else {
  stop("Incorrect number of presence points")
}

maxlambda <- max(exp(X %*% parameter), na.rm = TRUE)

# Vector, each entry containing possible number of background points that we wish to consider
npoints <- c(10^(0.5 * log10(100)), 100, 10^(mean(c(log10(100), log10(1000)))), 1000, 10^(mean(c(log10(1000), log10(10000)))), 10000, 10^(mean(c(log10(10000), log10(100000)))), 100000)
N <- 1000 # number of replications

windowarea <- (xmax(bioclim) - xmin(bioclim)) * (ymax(bioclim) - ymin(bioclim))

sgrid.loctot <- coordinates(s)
X.back1tot <- cbind(1, values(s))

result <- lapply(npoints, function(n) {
  # Get the grid resolution from the required number of points
  res <- sqrt(windowarea / n)
  s.loc <- xyFromCell(s, 1:ncell(s))
  # set up background grid (background resolution = resolution(s)*gridfactor)
  gridfactor <- res / mean(res(s))
  row <- ceiling(gridfactor * c((gridfactor + 1) / (2 * gridfactor) + 0:(nrow(s) / gridfactor - 1)))
  col <- ceiling(gridfactor * c((gridfactor + 1) / (2 * gridfactor) + 0:(ncol(s) / gridfactor - 1)))
  cellnum <- cellFromRowColCombine(s, row, col)
  print(length(cellnum))
  #cellnum <- 1:ncell(s)
  sgrid.loc <- xyFromCell(s, cellnum)
  sgrid.bio3 <- values(s)[, 'bio3'][cellnum]
  sgrid.bio8 <- values(s)[, 'bio8'][cellnum]
  # sgrid.bio9 <- values(s)[, 'bio9'][cellnum]
  # sgrid.bio13 <- values(s)[, 'bio13'][cellnum]
  # sgrid.bio14 <- values(s)[, 'bio14'][cellnum]
  X.back1 <- cbind(1,
                   sgrid.bio3,
                   sgrid.bio8)#,
                   # sgrid.bio9,
                   # sgrid.bio13,
                   # sgrid.bio14)
  
  # set initial values
  totlrbadd <- totwlrfit <- totfitIPP <- totfitIPPBT <- totppmIPP <- matrix(0, N, dim(X)[2])
  maxllfitIPP <- maxllfitIPPBT <- maxllwlrfit <- maxlllrbadd <- maxllppmIPP <- matrix(0, N)
  tmlrbadd <- tmwlrfit <- tmfitIPP <- tmppmIPP <- tmfitIPPBT <- 0
  
  # negative Poisson likelihood function (gridded only)
  
  likelihoodIPP <- function(param){
    lambda <- exp(X.back1 %*% param)
    mu <- lambda * windowarea / nrow(X.back1)
    logL.pp <- sum(X.po %*% param, na.rm = TRUE) - sum(mu, na.rm = TRUE)
    
    -sum(logL.pp)
  }
  
  likelihoodIPPtot <- function(param){
    lambda <- exp(X.back1tot %*% param)
    mu <- lambda * windowarea / nrow(X.back1tot)
    logL.pp <- sum(X.po %*% param, na.rm = TRUE) - sum(mu, na.rm = TRUE)
    
    -sum(logL.pp)
  }
  
  # negative Poisson likelihood function (X.back specifies Berman Turner, with area.back giving the weights)
  
  likelihoodIPPBT <- function(param){
    lambda <- exp(X.back %*% param)
    mu <- lambda * area.back
    logL.pp <- sum(X.po %*% param, na.rm = TRUE) - sum(mu, na.rm = TRUE)
    
    -sum(logL.pp)
  }
  
  likelihoodIPPBTtot <- function(param){
    lambda <- exp(X.backtot %*% param)
    mu <- lambda * area.backtot
    logL.pp <- sum(X.po %*% param, na.rm = TRUE) - sum(mu, na.rm = TRUE)
    
    -sum(logL.pp)
  }
  
  
  # negative weighted logistic regression likelihood function
  weight <- 1e3
  likelihoodwlr <- function(param) {
    alpha <- windowarea / (weight * dim(X.back1)[1])
    loglikewlr <- -(sum(X.po %*% param, na.rm = TRUE) - 
                      sum(log(1 + alpha * exp(X.po %*% param)), na.rm = TRUE) -
                      sum(weight * log(1 + alpha * exp(X.back1 %*% param)), na.rm = TRUE))
  }
  
  likelihoodwlrtot <- function(param) {
    alpha <- windowarea / (weight * dim(X.back1tot)[1])
    loglikewlr <- -(sum(X.po %*% param, na.rm = TRUE) -
                      sum(log(1 + alpha * exp(X.po %*% param)), na.rm = TRUE) - 
                      sum(weight * log(1 + alpha * exp(X.back1tot %*% param)), na.rm = TRUE))
  }
  
  # negative baddeley regression likelihood function
  
  likelihoodblr <- function(param){
    dum <- nrow(X.back1tot) / windowarea
    loglikblr <- sum(log(exp(X.po %*% param) / (exp(X.po %*% param) + dum)), na.rm = TRUE) + 
      sum(log(dum / (exp(X.back1tot %*% param) + dum)), na.rm = TRUE)
    loglikblr
    
  }
  
  set.seed(1)
  
  for(i in 1:N) {
    # sample presence points
    N.hpp <- rpois(1, maxlambda * windowarea)
    ind.hpp <- sample(1:ncell(s), size = N.hpp, replace = FALSE)
    loc.hpp <- s.loc[ind.hpp, ]
    lambda.hpp <- exp(X %*% parameter)[ind.hpp]
    
    ind.ipp <- runif(N.hpp, 0, 1) <= lambda.hpp / maxlambda
    N.ipp <- sum(ind.ipp, na.rm = TRUE)
    loc.ipp <- loc.hpp[ind.ipp, ]
    X.po <- X[ind.hpp[ind.ipp], ]
    X.poslope <- X.po[, 2:length(parameter)]
    
    #setting up Berman-Turner quadrature
    
    X.back <- rbind(X.po, X.back1)
    X.backtot <- rbind(X.po, X.back1tot)
    s.loc.ipp <- as.ppp(loc.ipp, W = w)
    sgrid.loc.ipp <- as.ppp(sgrid.loc, W = w)
    datadummy <- quadscheme(s.loc.ipp, sgrid.loc.ipp, method = "grid")
    area.back <- w.quad(datadummy)
    
    sgrid.loc.ipptot <- as.ppp(sgrid.loctot, W = w)
    datadummytot <- quadscheme(s.loc.ipp, sgrid.loc.ipptot, method = "grid")
    area.backtot <- w.quad(datadummytot)
    
    # fit Poisson likelihood
    
    tm <- Sys.time()
    fitIPP <- optim(par = array(0, length(parameter)), fn = likelihoodIPP, method = "BFGS")
    tmfitIPP <- tmfitIPP + Sys.time() - tm
    totfitIPP[i, ] <- fitIPP$par
    maxllfitIPP[i, ]<- -likelihoodIPPtot(totfitIPP[i,]) - sum(log(1:N.ipp))
    
    # fit Poisson likelihood BT
    tm <- Sys.time()
    fitIPPBT <- optim(par = array(0, length(parameter)), fn = likelihoodIPPBT, method = "BFGS")
    tmfitIPPBT <- tmfitIPPBT + Sys.time() - tm
    totfitIPPBT[i, ] <- fitIPPBT$par
    maxllfitIPPBT[i, ] <- -likelihoodIPPBTtot(totfitIPPBT[i,]) - sum(log(1:N.ipp))
    
    # put covariates into dataframe for use in glm
    
    tm <- Sys.time()
    colnames(X.poslope) <- colnames(X.back1)[2:dim(X.back1)[2]]
    X.polr <- cbind(X.poslope, X = array(1, dim(X.poslope)[1]))
    X.dpp <- cbind(X.back1[,2:dim(X.back1)[2]], X = array(0, dim(X.back1)[1]))
    mydata <- rbind(X.polr, X.dpp)
    offst <- array(log(dim(X.back1)[1] / windowarea), dim(mydata)[1])
    mydata <- as.data.frame(mydata)
    mydata$offst <- offst
    
    #fit Baddeley Regression
    
    lrbadd <- glm(X ~ offset(-offst) + sgrid.bio3 + sgrid.bio8# + sgrid.bio9 + sgrid.bio13 + sgrid.bio14 + 
                  , 
                  data = mydata,
                  family = "binomial")
    tmlrbadd <- tmlrbadd + Sys.time() - tm
    totlrbadd[i, ] <- coefficients(lrbadd)
    maxlllrbadd[i, ] <- likelihoodblr(totlrbadd[i, ])
    
    # fit Weighted logistic regression
    
    tm <- Sys.time()
    wlrfit <- optim(par = array(0, length(parameter)), fn = likelihoodwlr, method = "BFGS")
    tmwlrfit <- tmwlrfit + Sys.time() - tm
    totwlrfit[i, ] <- wlrfit$par
    maxllwlrfit[i, ] <- -likelihoodwlrtot(totwlrfit[i, ]) - sum(log(1:N.ipp))
    
    # fit using ppm
    
    tm <- Sys.time() 
    mydata2 <- as.data.frame(rbind(X.po, X.back1))
    mydata2 <- mydata2[, 2:dim(mydata2)[2]]
    ppmIPP <- ppm(datadummy, 
                  trend = ~ sgrid.bio3 + sgrid.bio8,# + sgrid.bio9 + sgrid.bio13 + sgrid.bio14, 
                  data = mydata2, 
                  method = 'mpl'
    )
    tmppmIPP <- tmppmIPP + Sys.time() - tm
    totppmIPP[i, ] <- coefficients(ppmIPP)
    maxllppmIPP[i, ] <- logLik(ppmIPP) - sum(log(1:N.ipp))
    
  }
  
  summary_results <- list()
  summary_results$IPP <- matrix(NA, length(summary(totfitIPP[, 1])), ncol(totfitIPP))
  for(i in seq_len(ncol(totfitIPP))) {
    summary_results$IPP[, i] <- summary(totfitIPP[, i])
  }
  
  summary_results$IPPBT <- matrix(NA, length(summary(totfitIPPBT[, 1])), ncol(totfitIPPBT))
  for(i in seq_len(ncol(totfitIPPBT))) {
    summary_results$IPPBT[, i] <- summary(totfitIPPBT[, i])
  }
  
  summary_results$wlr <- matrix(NA, length(summary(totwlrfit[, 1])), ncol(totwlrfit))
  for(i in seq_len(ncol(totwlrfit))) {
    summary_results$wlr[, i] <- summary(totwlrfit[, i])
  }
  
  summary_results$lr <- matrix(NA, length(summary(totlrbadd[, 1])), ncol(totlrbadd))
  for(i in seq_len(ncol(totlrbadd))) {
    summary_results$lr[, i] <- summary(totlrbadd[, i])
  }
  
  summary_results$ppm <- matrix(NA, length(summary(totppmIPP[, 1])), ncol(totppmIPP))
  for(i in seq_len(ncol(totppmIPP))) {
    summary_results$ppm[, i] <- summary(totppmIPP[, i])
  }
  
  for(i in seq_len(length(summary_results))) {
    colnames(summary_results[[i]]) <- names(parameter)
  }
  
  for(i in seq_len(length(summary_results))) {
    rownames(summary_results[[i]]) <- names(summary(totppmIPP[, 1]))
  }
  
  # We can also calculate the rmse for each fit
  
  rmse_results <- list()
  rmse_results$IPP <- matrix(NA, N)
  rmse_results$IPPBT <- matrix(NA, N)
  rmse_results$lr <- matrix(NA, N)
  rmse_results$lrbadd <- matrix(NA, N)
  rmse_results$ppmipp <- matrix(NA, N)
  for(i in 1:N) {
    rmse_results$IPP[i] <- rmse(parameter, totfitIPP[i, ])
    rmse_results$IPPBT[i] <- rmse(parameter, totfitIPPBT[i, ])
    rmse_results$lr[i] <- rmse(parameter, totwlrfit[i, ])
    rmse_results$lrbadd[i] <- rmse(parameter, totlrbadd[i, ])
    rmse_results$ppmipp[i] <- rmse(parameter, totppmIPP[i, ])
  }
  
  # and the log-likelihood
  
  ll_results <- list()
  ll_results$IPP <- matrix(NA, N)
  ll_results$IPPBT <- matrix(NA, N)
  ll_results$lr <- matrix(NA, N)
  ll_results$lrbadd <- matrix(NA, N)
  ll_results$ppmipp <- matrix(NA, N)
  for(i in 1:N) {
    ll_results$IPP[i] <- maxllfitIPP[i, ]
    ll_results$IPPBT[i] <- maxllfitIPPBT[i, ]
    ll_results$lr[i] <- maxllwlrfit[i, ]
    ll_results$lrbadd[i] <- maxlllrbadd[i, ]
    ll_results$ppmipp[i] <- maxllppmIPP[i, ]
  }
  
  # as well as the beta coefficients...
  
  beta <- vector(length = ncol(totfitIPP), mode = "list")
  for(j in seq_len(ncol(totfitIPP))) {
    beta[[j]]$IPP <- matrix(NA, N)
    beta[[j]]$IPPBT <- matrix(NA, N)
    beta[[j]]$lr <- matrix(NA, N)
    beta[[j]]$lrbadd <- matrix(NA, N)
    beta[[j]]$ppmipp <- matrix(NA, N)
    for(i in 1:N) {
      beta[[j]]$IPP[i] <- totfitIPP[i, j]
      beta[[j]]$IPPBT[i] <- totfitIPPBT[i, j]
      beta[[j]]$lr[i] <- totwlrfit[i, j]
      beta[[j]]$lrbadd[i] <- totlrbadd[i, j]
      beta[[j]]$ppmipp[i] <- totppmIPP[i, j]
    }
  }
  
  list(summary = summary_results,
       rmse = rmse_results,
       ll = ll_results,
       npoints = length(cellnum),
       resolution = res,
       beta = beta)
})

# save(result, file = "result_high_presence.RData")

df <- c()
for(i in seq_len(length(result))) {
  df <- rbind(df,
              cbind(result[[i]]$resolution * 60,
                    result[[i]]$npoints,
                    result[[i]]$rmse$IPP,
                    result[[i]]$ll$IPP,
                    result[[i]]$beta[[1]]$IPP,
                    result[[i]]$beta[[2]]$IPP,
                    result[[i]]$beta[[3]]$IPP,
                    # result[[i]]$beta[[4]]$IPP,
                    # result[[i]]$beta[[5]]$IPP,
                    # result[[i]]$beta[[6]]$IPP,
                    "PPM-GR"))
  df <- rbind(df,
              cbind(result[[i]]$resolution * 60,
                    result[[i]]$npoints,
                    result[[i]]$rmse$IPPBT,
                    result[[i]]$ll$IPPBT,
                    result[[i]]$beta[[1]]$IPPBT,
                    result[[i]]$beta[[2]]$IPPBT,
                    result[[i]]$beta[[3]]$IPPBT,
                    # result[[i]]$beta[[4]]$IPPBT,
                    # result[[i]]$beta[[5]]$IPPBT,
                    # result[[i]]$beta[[6]]$IPPBT,
                    "PPM-BT"))
  
  df <- rbind(df,
              cbind(result[[i]]$resolution * 60,
                    result[[i]]$npoints,
                    result[[i]]$rmse$lr,
                    result[[i]]$ll$lr,
                    result[[i]]$beta[[1]]$lr,
                    result[[i]]$beta[[2]]$lr,
                    result[[i]]$beta[[3]]$lr,
                    # result[[i]]$beta[[4]]$lr,
                    # result[[i]]$beta[[5]]$lr,
                    # result[[i]]$beta[[6]]$lr,
                    "IWLR"))
  
  df <- rbind(df,
              cbind(result[[i]]$resolution * 60,
                    result[[i]]$npoints,
                    result[[i]]$rmse$lrbadd,
                    result[[i]]$ll$lrbadd,
                    result[[i]]$beta[[1]]$lrbadd,
                    result[[i]]$beta[[2]]$lrbadd,
                    result[[i]]$beta[[3]]$lrbadd,
                    # result[[i]]$beta[[4]]$lrbadd,
                    # result[[i]]$beta[[5]]$lrbadd,
                    # result[[i]]$beta[[6]]$lrbadd,
                    "BLR"))
  
  df <- rbind(df,
              cbind(result[[i]]$resolution * 60,
                    result[[i]]$npoints,
                    result[[i]]$rmse$ppmipp,
                    result[[i]]$ll$ppmipp,
                    result[[i]]$beta[[1]]$ppmipp,
                    result[[i]]$beta[[2]]$ppmipp,
                    result[[i]]$beta[[3]]$ppmipp,
                    # result[[i]]$beta[[4]]$ppmipp,
                    # result[[i]]$beta[[5]]$ppmipp,
                    # result[[i]]$beta[[6]]$ppmipp,
                    "PPM-BT (spatstat)"))
  
}
df <- as.data.frame(df)
colnames(df) <- c("resolution", "npoints", "rmse", "ll", "intercept", "beta1", "beta2", #"beta3", "beta4", "beta5", 
                  "method")
df$rmse <- as.numeric(as.character(df$rmse))
df$ll <- as.numeric(as.character(df$ll))
df$resolution <- as.numeric(as.character(df$resolution))
df$npoints <- as.numeric(as.character(df$npoints))
df$intercept <- as.numeric(as.character(df$intercept))
df$beta1 <- as.numeric(as.character(df$beta1))
df$beta2 <- as.numeric(as.character(df$beta2))
# df$beta3 <- as.numeric(as.character(df$beta3))
# df$beta4 <- as.numeric(as.character(df$beta4))
# df$beta5 <- as.numeric(as.character(df$beta5))

# Fix offset for Baddeley logistic regression
a <- sd(df$ll[df$method != "BLR"]) / sd(df$ll[df$method == "BLR"])
b <- (mean(df$ll[df$method != "BLR" & df$npoints == max(df$npoints)]) - a * mean(df$ll[df$method == "BLR" & df$npoints == max(df$npoints)]))
df$ll[df$method == "BLR"] <- a * df$ll[df$method == "BLR"] + b

# Remove two versions of PPM-BT
df <- df[df$method != "PPM-BT (spatstat)", ]

df_medians <- df[!duplicated(df[, c("resolution", "npoints", "method")]), ]
for(i in seq_len(nrow(df_medians))) {
  df_medians$rmse[i] <- median(df$rmse[df$resolution == df_medians$resolution[i] & 
                                         df$npoints == df_medians$npoints[i] &
                                         df$method == df_medians$method[i]])
  df_medians$ll[i] <- median(df$ll[df$resolution == df_medians$resolution[i] & 
                                     df$npoints == df_medians$npoints[i] &
                                     df$method == df_medians$method[i]])
  df_medians$intercept[i] <- median(df$intercept[df$resolution == df_medians$resolution[i] & 
                                                   df$npoints == df_medians$npoints[i] &
                                                   df$method == df_medians$method[i]])
  df_medians$beta1[i] <- median(df$beta1[df$resolution == df_medians$resolution[i] & 
                                           df$npoints == df_medians$npoints[i] &
                                           df$method == df_medians$method[i]])
  df_medians$beta2[i] <- median(df$beta2[df$resolution == df_medians$resolution[i] & 
                                           df$npoints == df_medians$npoints[i] &
                                           df$method == df_medians$method[i]])
  # df_medians$beta3[i] <- median(df$beta3[df$resolution == df_medians$resolution[i] & 
  #                                          df$npoints == df_medians$npoints[i] &
  #                                          df$method == df_medians$method[i]])
  # df_medians$beta4[i] <- median(df$beta4[df$resolution == df_medians$resolution[i] & 
  #                                          df$npoints == df_medians$npoints[i] &
  #                                          df$method == df_medians$method[i]])
  # df_medians$beta5[i] <- median(df$beta5[df$resolution == df_medians$resolution[i] & 
  #                                          df$npoints == df_medians$npoints[i] &
  #                                          df$method == df_medians$method[i]])
}

# Code below computes bias and se
df$intercept_bias <- df$intercept - parameter[1]
df$beta1_bias <- df$beta1 - parameter[2]
df$beta2_bias <- df$beta2 - parameter[3]

bias <- expand.grid(npoints = unique(df$npoints), method = unique(df$method))
bias$mean_intercept_bias <- sapply(seq_len(nrow(bias)), function(i) mean(df$intercept_bias[df$npoints == bias$npoints[i] & 
                                                                                             df$method == bias$method[i]]))
bias$mean_beta1_bias <- sapply(seq_len(nrow(bias)), function(i) mean(df$beta1_bias[df$npoints == bias$npoints[i] & 
                                                                                     df$method == bias$method[i]]))
bias$mean_beta2_bias <- sapply(seq_len(nrow(bias)), function(i) mean(df$beta2_bias[df$npoints == bias$npoints[i] & 
                                                                                     df$method == bias$method[i]]))


bias$mean_intercept_sd <- sapply(seq_len(nrow(bias)), function(i) sd(df$intercept[df$npoints == bias$npoints[i] & 
                                                                                    df$method == bias$method[i]]))
bias$mean_beta1_sd <- sapply(seq_len(nrow(bias)), function(i) sd(df$beta1[df$npoints == bias$npoints[i] & 
                                                                            df$method == bias$method[i]]))
bias$mean_beta2_sd <- sapply(seq_len(nrow(bias)), function(i) sd(df$beta2[df$npoints == bias$npoints[i] & 
                                                                            df$method == bias$method[i]]))

bias








beta1 <- lapply(result, function(r) data.frame(IPP = r$beta[[2]]$IPP,
                                               IPPBT = r$beta[[2]]$IPPBT,
                                               lr = r$beta[[2]]$lr,
                                               lrbadd = r$beta[[2]]$lrbadd))

i <- 2
ggplot() +
  geom_point(aes(x = beta1[[i]]$lrbadd,
                 y = beta1[[7]]$IPPBT)) +
  theme_minimal(base_size = 20) +
  xlab(paste0("BLR estimate of BIO3 with ", result[[i]]$npoints, " dummy points")) +
  ylab("True PPM likelihood estimate of BIO3") +
  geom_hline(yintercept = parameter[2], colour = "red") +
  geom_vline(xintercept = parameter[2], colour = "red")

ggplot(data = df) + 
  geom_boxplot(aes(x = resolution, y = rmse, fill = method, colour = method, 
                   group = paste(method, resolution)),
               size = 0.2, alpha = 0.7, width = 0.1) +
  geom_line(data = df_medians, aes(x = resolution, y = rmse, colour = method)) +
  theme_minimal() +
  xlab("Resolution (arcminutes)") +
  ylab("RMSE") +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")


png(file = "rmse.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = rmse, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = rmse, colour = method), size = 1.5) +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("RMSE") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

png(file = "intercept.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = intercept, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = intercept, colour = method), size = 1.5) +
  geom_line(aes(x = npoints, y = parameter[1]), colour = "black") +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("Intercept") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

png(file = "bio3.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = beta1, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = beta1, colour = method), size = 1.5) +
  geom_line(aes(x = npoints, y = parameter[2]), colour = "black") +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("BIO3") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

png(file = "bio8.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = beta2, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = beta2, colour = method), size = 1.5) +
  geom_line(aes(x = npoints, y = parameter[3]), colour = "black") +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("BIO8") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

png(file = "bio9.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = beta3, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = beta3, colour = method), size = 1.5) +
  geom_line(aes(x = npoints, y = parameter[4]), colour = "black") +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("BIO9") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

png(file = "bio13.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = beta4, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = beta4, colour = method), size = 1.5) +
  geom_line(aes(x = npoints, y = parameter[5]), colour = "black") +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("BIO13") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

png(file = "bio14.png", bg = "white", width = 1100, height = 500)
ggplot(data = df[df$npoints > 50, ]) + 
  geom_boxplot(aes(x = npoints, y = beta5, fill = method, colour = method, 
                   group = paste(method, npoints)),
               size = 0.5, alpha = 0.7, width = 0.15, position = position_dodge(0.1)) +
  geom_line(data = df_medians[df_medians$npoints > 50, ], aes(x = npoints, y = beta5, colour = method), size = 1.5) +
  geom_line(aes(x = npoints, y = parameter[6]), colour = "black") +
  theme_minimal(base_size = 20) +
  xlab("Number of dummy points") +
  ylab("BIO14") +
  scale_x_continuous(trans = 'log10') +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
dev.off()

average_presence <- mean(exp(X %*% parameter)) * windowarea

ggplot(data = df) + 
  geom_boxplot(aes(x = log(average_presence) / log(npoints), y = rmse, fill = method, colour = method, 
                   group = paste(method, log(average_presence) / log(npoints))),
               size = 0.2, alpha = 0.7, width = 0.05) +
  theme_minimal() +
  xlab("Log(presence points) / log(dummy points)") +
  ylab("RMSE") +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")

ggplot(data = df) + 
  geom_boxplot(aes(x = log(npoints), y = ll, fill = method, colour = method, 
                   group = paste(method, log(npoints))),
               size = 0.5, alpha = 0.5, width = 0.1, position = position_dodge(0)) +
  geom_line(data = df_medians, aes(x = log(npoints), y = ll, colour = method), size = 1.5) +
  theme_minimal(base_size = 20) +
  xlab("Logarithm of the number of dummy points") +
  ylab("Log-likelihood") +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")

# min number of dummy points
min(sapply(result, function(x) x$npoints))
max(sapply(result, function(x) x$resolution)) * 60

# max number of dummy points
max(sapply(result, function(x) x$npoints))
min(sapply(result, function(x) x$resolution)) * 60
