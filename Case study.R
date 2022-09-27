library(raster)
library(spatstat)
library(dismo)
library(ecospat)

set.seed(2)

#load presence data with initial analysis

dattot <- read.csv("GBIF_capeweed.csv") 
dattot <- cbind(dattot$decimalLatitude,dattot$decimalLongitude,dattot$year)
dattot<- data.frame(dattot)
names(dattot)<-c("lat","lon","year")
dat <- dattot[dattot$year > 2012, ]
datval <- dattot[dattot$year < 2013 & dattot$year >2009, ]
dat <- data.frame(cbind(dat$lon,dat$lat))
names(dat)<-c("x","y")
datval <- data.frame(cbind(datval$lon,datval$lat))
names(datval)<-c("x","y")
dat <- dat[!duplicated(dat),]
datval <- datval[!duplicated(datval),]
dat.thin <- ecospat.occ.desaggregation(dat,min.dist=10/60)
dat.val.thin <- ecospat.occ.desaggregation(datval,min.dist=10/60)

#load covariate data

covraster.coarse <- stack("covariates.tif")
names(covraster.coarse) <- c('bio8','bio14','bio18','pH')
e <- extent(130,154,-44,-30) 
covraster.coarse1 <-extend(covraster.coarse,e)

#selecting grid points from total background raster (grid resolution = resolution(covraster.coarse)*gridfactor)
cov <- values(covraster.coarse)
gridfactor <- 2
row <- ceiling(gridfactor*c((gridfactor+1)/(2*gridfactor)+0:(nrow(covraster.coarse)/gridfactor-1)))
col <- ceiling(gridfactor*c((gridfactor+1)/(2*gridfactor)+0:(ncol(covraster.coarse)/gridfactor-1)))
cellnum <- cellFromRowColCombine(covraster.coarse,row,col)
#sgrid.loc <- raster::coordinates(covraster.coarse)
sgrid.loc <- xyFromCell(covraster.coarse,cellnum)
sgrid.loctot <- coordinates(covraster.coarse)

#set window

pixel <- !is.na(rowSums(cov))
pixel <- matrix(pixel,nrow=nrow(covraster.coarse),ncol=ncol(covraster.coarse),byrow=TRUE)
pixel <- pixel[c(nrow(pixel):1),]
victoriapixel <- owin(xrange=c(130,154),yrange=c(-44,-30),mask=pixel)
sgrid.loc.ipp <- ppp(sgrid.loc[,'x'],sgrid.loc[,'y'],window = victoriapixel)
sgrid.loc.ipptot <- ppp(sgrid.loctot[,'x'],sgrid.loctot[,'y'],window = victoriapixel)

#prepare data for spatstat

configuration <- ppp(x = dat.thin$x, y = dat.thin$y, window = victoriapixel)
configuration_val <- ppp(x = dat.val.thin$x, y = dat.val.thin$y, window = victoriapixel)
datadummylogi <- quadscheme.logi(configuration,sgrid.loc.ipp,method="grid")
datadummyppm <- quadscheme(configuration,sgrid.loc.ipp,method="grid")
datadummylogitot <- quadscheme.logi(configuration,sgrid.loc.ipptot,method="grid")
datadummyppmtot <- quadscheme(configuration,sgrid.loc.ipptot,method="grid")





bio8 <- as.im.RasterLayer(covraster.coarse[[1]])
bio14 <- as.im.RasterLayer(covraster.coarse[[2]])
bio18 <- as.im.RasterLayer(covraster.coarse[[3]])
pH <- as.im.RasterLayer(covraster.coarse[[4]])
covlist <- list(bio8=bio8,bio14=bio14,bio18=bio18,pH=pH)

#fit model

capeweedppm <- ppm(datadummyppm,trend=~bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='mpl')
capeweedlogi <- ppm(datadummylogi,trend=~bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='logi')
plot.ppm(capeweedppm,covariates=covlist,locations=victoriapixel,superimpose=TRUE,col=rainbow)
plot.ppm(capeweedlogi1,covariates=covlist,locations=outer_window,superimpose=TRUE)

#log likelihood calculations

likelihoodblr <- function(param){
  dum<-nrow(X.back1tot)/windowarea
  loglikblr <- sum(log(exp(X.po%*%param)/(exp(X.po%*%param)+dum)))+sum(log(dum/(exp(X.back1tot%*%param)+dum)))
  loglikblr
  
}

likelihoodIPPBTtot = function(param){
  lambda = exp(X.backtot %*% param)
  mu = lambda * area.backtot
  logL.pp = sum(X.po %*% param,na.rm=TRUE) - sum(mu,na.rm=TRUE)
  
  (-1)*sum(logL.pp)
}

bio8pH <- raster(covraster.coarse)
values(bio8pH) <- values(covraster.coarse)[,'bio8']*values(covraster.coarse)[,'pH']
covraster.coarse1 <- stack(covraster.coarse,bio8pH)
X.po <- cbind(datadummyppmtot$data$x,datadummyppmtot$data$y)
X.backtot <- rbind(cbind(datadummyppmtot$data$x,datadummyppmtot$data$x),cbind(datadummyppmtot$dummy$x,datadummyppmtot$dummy$y))
X.po <-cbind(1,extract(covraster.coarse1,X.po))
X.backtot <-cbind(1,extract(covraster.coarse1,X.backtot))
area.backtot <-w.quad(datadummyppmtot)
windowarea <- (xmax(covraster.coarse)-xmin(covraster.coarse))*(ymax(covraster.coarse)-ymin(covraster.coarse))
maxllfitIPP<- -likelihoodIPPBTtot(as.vector(coefficients(capeweedppm)))-sum(log(1:configuration$n))
X.back1tot <- cbind(datadummylogitot$dummy$x,datadummylogitot$dummy$y)
X.back1tot <-cbind(1,extract(covraster.coarse1,X.back1tot))

maxllblr <- likelihoodblr(as.vector(coefficients(capeweedlogi)))

#validation using AUC

predicted_intensity <- predict.ppm(capeweedppm1, type = "intensity",ngrid=dim(covraster.coarse)[c(1,2)])
resppm <- spatstat.core::auc(X = configuration_val, covariate = predicted_intensity)
predicted_intensity <- predict.ppm(capeweedlogi1, type = "intensity",ngrid=dim(covraster.coarse)[c(1,2)])
reslogi <- spatstat.core::auc(X = configuration_val, covariate = predicted_intensity)
