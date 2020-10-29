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
setwd("/Volumes/MUSIC\ BACK\ UP/Uni\ files")

#load presence data with initial analysis

dattot <- read.csv("capeweed data.csv")
dat <- dattot[dattot$year > 2014 & dattot$country == "E_AUS", ]
datval <- dattot[dattot$year < 2015 & dattot$year >2011 & dattot$country == "E_AUS", ]
pattern <- ppp(dat$lon,dat$lat,c(min(dat$lon),max(dat$lon)),c(min(dat$lat),max(dat$lat)))
summary(pattern)
patternval <- ppp(datval$lon,datval$lat,c(min(dat$lon),max(dat$lon)),c(min(dat$lat),max(dat$lat)))
summary(pattern)
plot(Kest(pattern,correction=c('bord','bord.modif')))
plot(density(pattern))

#load covariate data

pH <- raster('PHC_000_005_EV_N_P_AU_TRN_N_20140801.tif')
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
e <- extent(floor(min(dat$lon)),ceiling(max(dat$lon)),floor(min(dat$lat)),ceiling(max(dat$lat)))
e <- extent(130,154,-44,-30) 
bioclim.coarse <- crop(bioclimtot,e)
prec.coarse<-crop(prectot,e)
prec.coarse<-scale(prec.coarse,center=TRUE,scale=TRUE)
pH <- crop(pH,e)
P <- crop(P,e)
sand <- crop(sand,e)
pH.fine <- raster::shift(pH,dx=0.0004,dy=0.00042)
pH.coarse <- aggregate(pH.fine,fact=10,na.rm=TRUE)
P.fine <- raster::shift(P,dx=0.0004,dy=0.00042)
P.coarse <- aggregate(P.fine,fact=10,na.rm=TRUE)
sand.fine <- raster::shift(sand,dx=0.0004,dy=0.00042)
sand.coarse <- aggregate(sand.fine,fact=10,na.rm=TRUE)
#bioclim.coarse <- exclude(bioclim.coarse)
#bioclim.fine <- disaggregate(bioclim,fact=10)
#bioclim.fine <- crop(bioclim.fine,extent(P))
#bioclim.coarse <- crop(bioclim.coarse,extent(P))
#covraster.fine <- raster::stack(bioclim.fine[[c(9,14)]],sand.fine,pH.fine,P.fine)
#covraster.fine <- scale(covraster.fine,center=TRUE,scale=TRUE)
covraster.coarse <- raster::stack(bioclim.coarse,sand.coarse,pH.coarse,P.coarse,r)
covraster.coarse <- scale(covraster.coarse,center=TRUE,scale=TRUE)
#covraster<-aggregate(covraster,fact=4,na.rm=TRUE)
names(covraster.coarse) <- c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11','bio12','bio13'
                             ,'bio14','bio15','bio16','bio17','bio18','bio19','sand','pH','P','bias')

#selecting grid points from total background raster


configuration <- spatstat::ppp(x = dat$lon, y = dat$lat, window = victoriapixel)
configuration_val <- spatstat::ppp(x = datval$lon, y = datval$lat, window = victoriapixel)
datadummylogi <- quadscheme.logi(configuration,sgrid.loc.ipp,method="grid")
datadummyppm <- quadscheme(configuration,sgrid.loc.ipp,method="grid")
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
capeweedppm160 <- ppm(datadummyppm,trend=~bias+bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='mpl')
capeweedlogi160 <- ppm(datadummylogi,trend=~bias+bio8+bio14+bio18+pH+bio8*pH,covariates = covlist,method='logi')
#fit and plot results

#capeweedppm <- ppm(datadummyppm,trend=~bio14+bio15+bio18+bio3+pH+sand,data=covdf,method='mpl')
#capeweedlogi <- ppm(datadummylogi,trend=~bio14+bio15+bio18+bio3+pH+sand,data=covdf,method='logi')
#capeweedppm2 <- ppm(datadummyppm,trend=~density+pH,data=covdf,method='mpl')
#capeweedlogi2 <- ppm(datadummylogi,trend=~bio1+bio4+bio14+bio15+sand+pH,data=covdf,method='logi')



#using ppmlasso instead

#data<-cbind(coords(sgrid.loc.ipp),X.back)
#lassofit<-ppmlasso(~+bio14+bio15+bio18+bio3+pH+sand,sp.xy=coords(configuration),coord=c('x','y'),env.grid=data,sp.scale=res(covraster.coarse),criterion='nlgcv')

#plot (across E Aus)

bio1 <- as.im.RasterLayer(covraster.coarse[[1]])
bio2 <- as.im.RasterLayer(covraster.coarse[[2]])
bio3 <- as.im.RasterLayer(covraster.coarse[[3]])
bio4 <- as.im.RasterLayer(covraster.coarse[[4]])
bio5 <- as.im.RasterLayer(covraster.coarse[[5]])
bio6 <- as.im.RasterLayer(covraster.coarse[[6]])
bio7 <- as.im.RasterLayer(covraster.coarse[[7]])
bio8 <- as.im.RasterLayer(covraster.coarse[[8]])
bio9 <- as.im.RasterLayer(covraster.coarse[[9]])
bio10 <- as.im.RasterLayer(covraster.coarse[[10]])
bio11 <- as.im.RasterLayer(covraster.coarse[[11]])
bio12 <- as.im.RasterLayer(covraster.coarse[[12]])
bio13 <- as.im.RasterLayer(covraster.coarse[[13]])
bio14 <- as.im.RasterLayer(covraster.coarse[[14]])
bio15 <- as.im.RasterLayer(covraster.coarse[[15]])
bio16 <- as.im.RasterLayer(covraster.coarse[[16]])
bio17 <- as.im.RasterLayer(covraster.coarse[[17]])
bio18 <- as.im.RasterLayer(covraster.coarse[[18]])
bio19 <- as.im.RasterLayer(covraster.coarse[[19]])
sand <- as.im.RasterLayer(covraster.coarse[[20]])
pH <- as.im.RasterLayer(covraster.coarse[[21]])
P <- as.im.RasterLayer(covraster.coarse[[22]])
bias <- as.im.RasterLayer(covraster.coarse[[23]])
covlist <- list(bio1=bio1,bio2=bio2,bio3=bio3,bio4=bio4,bio5=bio5,bio6=bio6,bio7=bio7,bio8=bio8,bio9=bio9,
                bio10=bio10,bio11=bio11,bio12=bio12,bio13=bio13,bio14=bio14,bio15=bio15,bio16=bio16,
                bio17=bio17,bio18=bio18,bio19=bio19,sand=sand,P=P,pH=pH,bias=bias)
capeweedppm <- ppm(datadummyppm,trend=~bias+bio8+bio14+bio18+pH+bio8*pH+bio18*pH,covariates = covlist,method='mpl')
plot.ppm(capeweedppm80,covariates=covlist,locations=victoriapixel,superimpose=FALSE,col=rainbow)
plot.ppm(capeweedlogi1,covariates=covlist,locations=outer_window,superimpose=TRUE)

#plot (study window only)

plot.ppm(capeweedppm1,covariates=covlist,locations=inner_im_window,superimpose=FALSE)
plot.ppm(capeweedlogi,covariates=covlist,locations=inner_im_window,superimpose=TRUE)

#validation using AUC
capeweedppm <- ppm(datadummyppm,trend=~bias+bio8+autprec+bio18+pH,covariates = covlist,method='mpl')
predicted_intensity <- predict(capeweedppm, type = "intensity")
res <- spatstat::auc(X = configuration_val, covariate = predicted_intensity)

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