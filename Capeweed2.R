#capeweed modelling
library(raster)
library(mapview)
library(usdm)
library(spatstat)
library(maptools)
library(Metrics)
library(plotrix)
library(dismo)
library(slga)


#load presence data with initial analysis
setwd("/Users/jamesbellew/Desktop/Sean Honours")
dattot <- read.csv("GBIF_capeweed.csv")
dattot <- cbind(dattot$decimalLatitude,dattot$decimalLongitude,dattot$year)
dattot<- data.frame(dattot)
names(dattot)<-c("lat","lon","year")
#dattot <- read.csv("capeweed data.csv")
dat <- dattot[dattot$year > 2012, ]
datval <- dattot[dattot$year < 2013 & dattot$year >2009, ]
pattern <- ppp(dat$lon,dat$lat,c(130,154),c(-44,-30))
summary(pattern)
patternval <- ppp(datval$lon,datval$lat,c(130,154),c(-44,-30))
summary(pattern)
plot(Kest(pattern,correction=c('bord','bord.modif')))
plot(density(pattern))

dat <- data.frame(cbind(dat$lon,dat$lat))
names(dat)<-c("x","y")
datval <- data.frame(cbind(datval$lon,datval$lat))
names(datval)<-c("x","y")

dat <- dat[!duplicated(dat),]
datval <- datval[!duplicated(datval),]
dat.thin <- ecospat.occ.desaggregation(dat,min.dist=10/60)
dat.val.thin <- ecospat.occ.desaggregation(datval,min.dist=10/60)
#load covariate data

pH <- raster('pHc_000_005_EV_N_P_AU_NAT_C_20140801.tif')
sand <- raster('SND_000_005_EV_N_P_AU_TRN_N_20140801.tif')
P <- raster('PTO_000_005_EV_N_P_AU_TRN_N_20140801.tif')
bioclim1<-raster::getData('worldclim',var='bio',res=0.5,lon=120,lat=-60)
bioclim2<-raster::getData('worldclim',var='bio',res=0.5,lon=150,lat=-60)
bioclim3<-raster::getData('worldclim',var='bio',res=0.5,lon=120,lat=0)
bioclim4<-raster::getData('worldclim',var='bio',res=0.5,lon=150,lat=0)
bioclimtot <- merge(bioclim1,bioclim2,bioclim3,bioclim4)
prec1<-raster::getData('worldclim',var='prec',res=0.5,lon=120,lat=-60)
prec2<-raster::getData('worldclim',var='prec',res=0.5,lon=150,lat=-60)
prec3<-raster::getData('worldclim',var='prec',res=0.5,lon=120,lat=0)
prec4<-raster::getData('worldclim',var='prec',res=0.5,lon=150,lat=0)
prectot <- merge(prec1,prec2,prec3,prec4)
bioclimtot<-raster::getData('worldclim',var='bio',res=2.5)
bio8<-raster("wc2.1_2.5m_bio_8.tif")
bio14<-raster("wc2.1_2.5m_bio_14.tif")
bio18<-raster("wc2.1_2.5m_bio_18.tif")
bioclimtot<-stack(bio8,bio14,bio18)
#e <- extent(floor(min(dat$lon)),ceiling(max(dat$lon)),floor(min(dat$lat)),ceiling(max(dat$lat)))
e <- extent(130,154-res(bioclimtot)[[1]],-44,-30) 
bioclim.coarse <- crop(bioclimtot,e)
#prec.coarse<-crop(prectot,e)
#prec.coarse<-scale(prec.coarse,center=TRUE,scale=TRUE)
pH<-resample(pH,bioclim.coarse)
#P <- crop(P,e)
#sand <- crop(sand,e)
#pH.fine <- raster::shift(pH,dx=0.0004,dy=0.00042)
#pH.coarse <- aggregate(pH.fine,fact=10,na.rm=TRUE)
#P.fine <- raster::shift(P,dx=0.0004,dy=0.00042)
#P.coarse <- aggregate(P.fine,fact=10,na.rm=TRUE)
#sand.fine <- raster::shift(sand,dx=0.0004,dy=0.00042)
#sand.coarse <- aggregate(sand.fine,fact=10,na.rm=TRUE)
#bioclim.coarse <- exclude(bioclim.coarse)
#bioclim.fine <- disaggregate(bioclim,fact=10)
#bioclim.fine <- crop(bioclim.fine,extent(P))
#bioclim.coarse <- crop(bioclim.coarse,extent(P))
#covraster.fine <- raster::stack(bioclim.fine[[c(9,14)]],sand.fine,pH.fine,P.fine)
#covraster.fine <- scale(covraster.fine,center=TRUE,scale=TRUE)
covraster.coarse <- raster::stack(bioclim.coarse,pH)
covraster.coarse <- scale(covraster.coarse,center=TRUE,scale=TRUE)
#covraster<-aggregate(covraster,fact=4,na.rm=TRUE)
names(covraster.coarse) <- c('bio8','bio14','bio18','pH')
e <- extent(130,154,-44,-30) 
covraster.coarse1 <-extend(covraster.coarse,e)
#selecting grid points from total background raster
cov <- values(covraster.coarse)
gridfactor <- 8
row <- ceiling(gridfactor*c((gridfactor+1)/(2*gridfactor)+0:(nrow(covraster.coarse)/gridfactor-1)))
col <- ceiling(gridfactor*c((gridfactor+1)/(2*gridfactor)+0:(ncol(covraster.coarse)/gridfactor-1)))
cellnum <- cellFromRowColCombine(covraster.coarse,row,col)
#sgrid.loc <- raster::coordinates(covraster.coarse)
sgrid.loc <- xyFromCell(covraster.coarse,cellnum)
sgrid.loctot <- coordinates(covraster.coarse)

#outer_window <- owin(xrange=c(130,153.9583),yrange=c(-44,-30))
#m <- spatstat::pixellate(spatstat::ppp(dat$lon,dat$lat,window=outer_window),eps=c(1,1))>0
#m$v[!m$v] <- NA
#inner_im_window <- spatstat::as.owin(m)
#victoriam <- as.owin(m)
pixel <- !is.na(rowSums(cov))
pixel <- matrix(pixel,nrow=nrow(covraster.coarse),ncol=ncol(covraster.coarse),byrow=TRUE)
pixel <- pixel[c(nrow(pixel):1),]
victoriapixel <- owin(xrange=c(130,154),yrange=c(-44,-30),mask=pixel)
#inner_im_window <- intersect.owin(victoriam,victoriapixel)
sgrid.loc.ipp <- ppp(sgrid.loc[,'x'],sgrid.loc[,'y'],window = victoriapixel)
sgrid.loc.ipptot <- ppp(sgrid.loctot[,'x'],sgrid.loctot[,'y'],window = victoriapixel)


configuration <- ppp(x = dat.thin$x, y = dat.thin$y, window = victoriapixel)
configuration_val <- ppp(x = dat.val.thin$x, y = dat.val.thin$y, window = victoriapixel)
datadummylogi <- quadscheme.logi(configuration,sgrid.loc.ipp,method="grid")
datadummyppm <- quadscheme(configuration,sgrid.loc.ipp,method="grid")
datadummylogitot <- quadscheme.logi(configuration,sgrid.loc.ipptot,method="grid")
datadummyppmtot <- quadscheme(configuration,sgrid.loc.ipptot,method="grid")
#densityim<-density(configuration)
#density.back<-densityim[inner_im_window]
#density.values<-densityim[configuration]
#X.values <- extract(covraster.coarse,coords(configuration))
#X.values <- cbind(X.values,density.values)
#X.back <- extract(covraster.coarse,coords(sgrid.loc.ipp))
#X.back <- cbind(X.back,density.back)
#cov <- rbind(X.values,X.back)
#covdf <- as.data.frame(cov)
#X <-as.data.frame(X.back)
#colnames(covdf) <- c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11','bio12','bio13'
          #           ,'bio14','bio15','bio16','bio17','bio18','bio19','sand','pH','P','density')
#capeweedppm <- ppm(datadummyppm,trend=~bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='mpl')

#fit and plot results

#capeweedppm <- ppm(datadummyppm,trend=~bio14+bio15+bio18+bio3+pH+sand,data=covdf,method='mpl')
#capeweedlogi <- ppm(datadummylogi,trend=~bio14+bio15+bio18+bio3+pH+sand,data=covdf,method='logi')
#capeweedppm2 <- ppm(datadummyppm,trend=~density+pH,data=covdf,method='mpl')
#capeweedlogi2 <- ppm(datadummylogi,trend=~bio1+bio4+bio14+bio15+sand+pH,data=covdf,method='logi')



#using ppmlasso instead

#data<-cbind(coords(sgrid.loc.ipp),X.back)
#lassofit<-ppmlasso(~+bio14+bio15+bio18+bio3+pH+sand,sp.xy=coords(configuration),coord=c('x','y'),env.grid=data,sp.scale=res(covraster.coarse),criterion='nlgcv')

#plot (across E Aus)

#bio1 <- as.im.RasterLayer(covraster.coarse[[1]])
#bio2 <- as.im.RasterLayer(covraster.coarse[[2]])
#bio3 <- as.im.RasterLayer(covraster.coarse[[3]])
#bio4 <- as.im.RasterLayer(covraster.coarse[[4]])
#bio5 <- as.im.RasterLayer(covraster.coarse[[5]])
#bio6 <- as.im.RasterLayer(covraster.coarse[[6]])
#bio7 <- as.im.RasterLayer(covraster.coarse[[7]])
bio8 <- as.im.RasterLayer(covraster.coarse[[1]])
#bio9 <- as.im.RasterLayer(covraster.coarse[[9]])
#bio10 <- as.im.RasterLayer(covraster.coarse[[10]])
#bio11 <- as.im.RasterLayer(covraster.coarse[[11]])
#bio12 <- as.im.RasterLayer(covraster.coarse[[12]])
#bio13 <- as.im.RasterLayer(covraster.coarse[[13]])
bio14 <- as.im.RasterLayer(covraster.coarse[[2]])
#bio15 <- as.im.RasterLayer(covraster.coarse[[15]])
#bio16 <- as.im.RasterLayer(covraster.coarse[[16]])
#bio17 <- as.im.RasterLayer(covraster.coarse[[17]])
bio18 <- as.im.RasterLayer(covraster.coarse[[3]])
#bio19 <- as.im.RasterLayer(covraster.coarse[[19]])
#sand <- as.im.RasterLayer(covraster.coarse[[20]])
pH <- as.im.RasterLayer(covraster.coarse[[4]])
#P <- as.im.RasterLayer(covraster.coarse[[22]])
#covlist <- list(bio1=bio1,bio2=bio2,bio3=bio3,bio4=bio4,bio5=bio5,bio6=bio6,bio7=bio7,bio8=bio8,bio9=bio9,
                #bio10=bio10,bio11=bio11,bio12=bio12,bio13=bio13,bio14=bio14,bio15=bio15,bio16=bio16,
                #bio17=bio17,bio18=bio18,bio19=bio19,sand=sand,P=P,pH=pH)
covlist <- list(bio8=bio8,bio14=bio14,bio18=bio18,pH=pH)
capeweedppm <- ppm(datadummyppm,trend=~bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='mpl')
capeweedlogi <- ppm(datadummylogi,trend=~bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='logi')
plot.ppm(capeweedppm,covariates=covlist,locations=victoriapixel,superimpose=TRUE,col=rainbow)
plot.ppm(capeweedlogi1,covariates=covlist,locations=outer_window,superimpose=TRUE)

#log likelihood

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



#plot (study window only)

#plot.ppm(capeweedppm1,covariates=covlist,locations=inner_im_window,superimpose=FALSE)
#plot.ppm(capeweedlogi,covariates=covlist,locations=inner_im_window,superimpose=TRUE)

#validation using AUC

predicted_intensity <- predict.ppm(capeweedppm1, type = "intensity",ngrid=dim(covraster.coarse)[c(1,2)])
resppm <- spatstat.core::auc(X = configuration_val, covariate = predicted_intensity)
predicted_intensity <- predict.ppm(capeweedlogi1, type = "intensity",ngrid=dim(covraster.coarse)[c(1,2)])
reslogi <- spatstat.core::auc(X = configuration_val, covariate = predicted_intensity)

resppm2 <- resppm
reslogi2 <- reslogi
maxllfitIPP2 <- maxllfitIPP
maxllblr2 <- maxllblr
capeweedppm2 <- capeweedppm
capeweedlogi2 <-capeweedlogi

resppm3 <- resppm
reslogi3 <- reslogi
maxllfitIPP3 <- maxllfitIPP
maxllblr3 <- maxllblr
capeweedppm3 <- capeweedppm
capeweedlogi3 <-capeweedlogi

resppm4 <- resppm
reslogi4 <- reslogi
maxllfitIPP4 <- maxllfitIPP
maxllblr4 <- maxllblr
capeweedppm4 <- capeweedppm
capeweedlogi4 <-capeweedlogi

resppm6 <- resppm
reslogi6 <- reslogi
maxllfitIPP6 <- maxllfitIPP
maxllblr6 <- maxllblr
capeweedppm6 <- capeweedppm
capeweedlogi6 <-capeweedlogi

resppm8 <- resppm
reslogi8 <- reslogi
maxllfitIPP8 <- maxllfitIPP
maxllblr8 <- maxllblr
capeweedppm8 <- capeweedppm
capeweedlogi8 <-capeweedlogi

resppm12 <- resppm
reslogi12 <- reslogi
maxllfitIPP12 <- maxllfitIPP
maxllblr12 <- maxllblr
capeweedppm12 <- capeweedppm
capeweedlogi12 <-capeweedlogi

resppm16 <- resppm
reslogi16 <- reslogi
maxllfitIPP16 <- maxllfitIPP
maxllblr16 <- maxllblr
capeweedppm16 <- capeweedppm
capeweedlogi16 <-capeweedlogi

resppm24 <- resppm
reslogi24 <- reslogi
maxllfitIPP24 <- maxllfitIPP
maxllblr24 <- maxllblr
capeweedppm24 <- capeweedppm
capeweedlogi24 <-capeweedlogi

resppm48 <- resppm
reslogi48 <- reslogi
maxllfitIPP48 <- maxllfitIPP
maxllblr48 <- maxllblr
capeweedppm48 <- capeweedppm
capeweedlogi48 <-capeweedlogi



plot(2.5*c(1,2,3,4,6,8,12,16,24),883.619+c(maxllblr1,maxllblr2,maxllblr3,maxllblr4,maxllblr6,maxllblr8,maxllblr12,maxllblr16,maxllblr24),type='o',col='red',xlab= 'resolution (arcminutes)',ylab='log-likelihood',main='log-likelihood')
lines(2.5*c(1,2,3,4,6,8,12,16,24),c(maxllfitIPP1,maxllfitIPP2,maxllfitIPP3,maxllfitIPP4,maxllfitIPP6,maxllfitIPP8,maxllfitIPP12,maxllfitIPP16,maxllfitIPP24),type='o')
legend("bottomleft",c("IPP-BT","BLR"), fill=c("black","red") ,bty="n")
plot(2.5*c(1,2,3,4,6,8,12,16,24),c(resppm1,resppm2,resppm3,resppm4,resppm6,resppm8,resppm12,resppm16,resppm24),type='o',ylim=c(0.77,0.82),xlab='resolution (arcminutes)',ylab='AUC',main='AUC')
lines(2.5*c(1,2,3,4,6,8,12,16,24),c(reslogi1,reslogi2,reslogi3,reslogi4,reslogi6,reslogi8,reslogi12,reslogi16,reslogi24),type='o',col='red')
legend("bottomleft",c("IPP-BT","BLR"), fill=c("black","red") ,bty="n")

plot(raster(predicted_intensity),col=rainbow(99,start=0,end=0.65,rev=TRUE))
#further analysis

plot(Kinhom(X=configuration,lambda=capeweedppm7))
plot(Kinhom(X=configuration,lambda=capeweedlogi))
plot(envelope(configuration,Kinhom,nsim=99,simulate=expression(rpoispp(predict.ppm(capeweedppm1,covariates=covlist,locations=victoriapixel)))))

anova.ppm(capeweedppm,test='Chisq')
anova.ppm(capeweedlogi,test='Chisq')

#using ppmlasso

data<-cbind(coords(sgrid.loc.ipp),X.back)

lassofit<-ppmlasso(~+bio14+bio15+bio18+bio3+pH+sand,sp.xy=coords(configuration),coord=c('x','y'),env.grid=data,sp.scale=1,criterion='nlgcv')


i=1
plot(res,ylim=c(2.047,2.09),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='intercept')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")
i=2
plot(res,ylim=c(0.03,0.24),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='bio14')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")
i=3
plot(res,ylim=c(0.31,0.38),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='bio15')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")
i=
plot(res,ylim=c(-0.53,-0.39),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='bio18')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")
i=5
plot(res,ylim=c(0.2,0.26),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='bio3')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")
i=6
plot(res,ylim=c(-0.524,-0.35),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='pH')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")
i=7
plot(res,ylim=c(0.29,0.57),c(param[2,i],param[3,i],param[4,i],param[5,i],param[6,i],param[7,i],param[8,i]),type='o',xlab='n0',ylab='parameter',main='sand')
lines(res,c(param1[2,i],param1[3,i],param1[4,i],param1[5,i],param1[6,i],param1[7,i],param1[8,i]),col='red',type='o')
legend("topright",c("IPP","logi"), fill=c("black","red") ,bty="n")