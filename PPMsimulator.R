remove(list=ls())
library(raster)
library(mvtnorm)
library(spatstat)
"function simulates intensity function of model lambda=parameter1*w+parameter2*u, 
where u and w are covariates that are inbuilt within the function"
PPMsimulator <- function(parameter1,parameter2,minx,maxx,miny,maxy) {
  #setting up window
  window<-as.owin(c(minx,maxx,miny,maxy))
  windowarea <- (maxx-minx)*(maxy-miny)
  npixels.x <- (maxx-minx)/0.002
  npixels.y <- (maxy-miny)/0.002
  s <- raster(ncol=npixels.x, nrow=npixels.y, xmn=minx, xmx=maxx, ymn=miny, ymx=maxy)
  s.loc <- xyFromCell(s, 1:ncell(s))
  #covariate u
  ucov <- dmvnorm(s.loc,mean=c(2,3),sigma=matrix(c(1,0.2,0.2,1),ncol=2))
  ucov <- (ucov - mean(ucov))/sd(ucov)
  values(s) <- ucov
  names(s) <- 'u'
  #covariate w
  wcov<-0.5*dmvnorm(s.loc,mean=c(1,2),sigma=matrix(c(2,0.5,0.5,1),ncol=2))
  wcov <- (wcov - mean(wcov))/sd(wcov)
  r<-raster(s)
  values(r) <- wcov
  names(r) <- 'w'
  s <- addLayer(s, r)
  t<-values(s)
  #assigning lambda values across the grid
  B <- c(parameter1,parameter2)
  X <- cbind(values(s)[,'w'], values(s)[,'u'])
  values(s) <- exp(X%*%B)
  maxlambda <- max(values(s))
  #simulating homogeneous poisson
  N.hpp <- rpois(1,maxlambda*windowarea)
  ind.hpp <- sample(1:ncell(s), size=N.hpp, replace=FALSE)
  loc.hpp <- s.loc[ind.hpp, ]
  lambda.hpp <- values(s)[ind.hpp]
  #thinning to produce inhomogeneous poisson
  ind.ipp <- runif(N.hpp, 0,1) <= lambda.hpp/maxlambda
  N.ipp <- sum(ind.ipp)
  loc.ipp <- loc.hpp[ind.ipp, ]
  X.po <- X[ind.hpp[ind.ipp], ]
  plot(loc.ipp)
  #fitting IPP using Maximum Likelihood function and numerical methods
  likelihoodIPP=function(param){
    l<- -sum(X.po%*%c(param[1],param[2]))+windowarea/(npixels.x*npixels.y-N.ipp)*(sum(exp(X%*%c(param[1],param[2])))-sum(exp(X.po%*%c(param[1],param[2]))))
  }
  fitIPP <- optim(par=c(parameter1+0.5,parameter2-0.5),fn=likelihoodIPP)
  #logistic regression
  dummyintensity<-4*npixels.x/4*(length(X.po)/2)/windowarea
  N.dpp <- rpois(1,dummyintensity*windowarea)
  ind.dpp <- sample((1:ncell(s))[-ind.hpp[ind.ipp]], size=N.dpp, replace=FALSE)
  X.dpp <- X[ind.dpp,]
  X.dpp <-cbind(X.dpp,array(0,length(X.dpp/2)))
  X.polr <-cbind(X.po,array(1,length(X.po/2)))
  mydata<-rbind(X.polr,X.dpp)
  mydata<-as.data.frame(mydata)
  lrfit <- glm(V3~V1+V2,data=mydata,family="binomial")
  #poisson regression
  IPPglm <- glm(V3~V1+V2,data=mydata,family="poisson")
  #Weighted logistic regression
  weight <- 10000
  "X.dpp <- X[ind.dpp,]
  X.dpp <-cbind(X.dpp,array(-weight,length(X.dpp/2)))
  mydata<-rbind(X.polr,X.dpp)
  mydata<-as.data.frame(mydata)
  wlrfit <- glm(V3~V1+V2,data=mydata,family=)"
  likelihoodwlr <- function(param){
    weight = 1000
    loglikewlr <- -(sum(log(windowarea/(weight*N.dpp))+X.po%*%c(param[1],param[2]))+sum(weight*log(1+windowarea/(weight*N.dpp)*exp(X.po%*%c(param[1],param[2]))))-sum(weight*log(1+windowarea/(weight*N.dpp)*exp(X%*%c(param[1],param[2]))))-sum(log(1+windowarea/(weight*N.dpp)*exp(X.po%*%c(param[1],param[2])))))
  }
  wlrfit <- optim(par=c(parameter1+0.5,parameter2-0.5),fn=likelihoodwlr)
  #Summary
  mylist<-list(lrfit,fitIPP,IPPglm,wlrfit)
  return(mylist)
  
}


