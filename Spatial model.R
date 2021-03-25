library(spBayes)
library(coda)
library(scatterplot3d)
par(mar=c(5,6,4,4)+.1)

###### Spatial Model (Gaussian case)

#####Make the data frame

##### Must download the data from: https://www.kaggle.com/sudalairajkumar/daily-temperature-of-major-cities
# and save it to the data frame city_temperature
city_temperature <- read.csv("~/Documents/Lecture Notes/Year 4/Project IV/Final Code/city_temperature.csv")
#####

citytemp2000=city_temperature[which(city_temperature$Year==2000),]
ctjun2000=citytemp2000[which(citytemp2000$Month==6),]
ct1jun2000=ctjun2000[which(ctjun2000$Day == 1),]
ct1jun2000=ct1jun2000[,-c(3,5,6,7)]
ctafrica=ct1jun2000[which(ct1jun2000$Region=='Africa'),]

citiesv = c(
  'Libreville',
  'Windhoek',
  'Capetown',
  'Kampala',
  'Lusaka',
  'Maputo',
  'Lilongwe',
  'Brazzaville',
  'Dar Es Salaam',
  'Bangui',
  'Addis Ababa',
  'Cairo'
)

afcities = data.frame(
  'City'=citiesv,
  'Easting'=c(550563,714444,261791,453545,637707,458823,584013,532397,519706,228629,474018,329653)+625000*c(0,1,2,4,3,4,4,1,5,2,5,4),
  'Northing'=c(43119,7503256,6243850,38421,8295178,7129376,8476034,9528267,9251009,482496,995458,3323794)+87500000*c(6,2,0,7,4,2,4,5,5,7,8,10)
)
# Googled these and shifted them by the UTM zones to make them relatively correct

plot(afcities$Easting,afcities$Northing,col=c(1:8,1:4),pch=c(rep(19,8),rep(20,4)))
legend('topleft',legend=citiesv,col=c(1:8,1:4),pch=c(rep(19,8),rep(20,4)))

spmoddat=data.frame(
  'City'=citiesv,
  'Eastkm'=afcities$Easting/1000,
  'Northkm'=afcities$Northing/1000,
  'Temp'=c(79.8,61.4,57.7,72.4,67.5,68.3,64.5,77.8,74.8,78.2,65.9,79.8)
)
#Note: values for Dar Es Salaam and Lilongwe are from the 3rd and 5th of June respectively, as they had no value on the 1st.

#####Fit the model

set.seed(17)

coords=as.matrix(spmoddat[,2:3])

n.samples <- 2000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors.1 <- list("beta.Norm"=list(70, diag(20,1)),
                 "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

cov.model='exponential'

n.report=500
verbose=TRUE

m.1 <- spLM(Temp~1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report,data=spmoddat)

burn.in <- 0.5*n.samples

#####recover beta and spatial random effects

m.1 <- spRecover(m.1, verbose=FALSE)
par(mfrow=c(1,2))
plot(as.vector(m.1$p.beta.recover.samples[1:1000]),col=6,type='l',xlim=c(0,2000),ylab=expression(beta))
lines(x=1001:2000,y=as.vector(m.1$p.beta.recover.samples[1001:2000]),col=1,type='l')
plot(as.vector(m.1$p.theta.recover.samples[1:1000,1]),col=6,type='l',xlim=c(0,2000),ylab=expression(sigma^2))
lines(x=1001:2000,y=as.vector(m.1$p.theta.recover.samples[1001:2000,1]),col=1,type='l')
par(mfrow=c(1,1))

m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.1$p.beta.recover.samples)$quantiles[c(3,1,5)],2)

m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]


##### And prediction with surface plot

cogrid=function(coords){
  n=length(coords[,1])
  final=matrix(0,nrow=n^2,ncol=2)
  final[,1]=rep(coords[,1],n)
  final[,2]=rep(coords[,2],each=n)
  return(final)
}

pred.coords=cogrid(coords)

pred.covars=matrix(1)

set.seed(5)
ypred=matrix(0,nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    coordsmat=matrix(pred.coords[12*(j-1)+i,],ncol=2)
    m.1.pred=spPredict(m.1,pred.coords = coordsmat,pred.covars=pred.covars,start=burn.in,verbose=FALSE)
    ypred[i,j]=mean(m.1.pred$p.y.predictive.samples[1:5])
  }
}
xord=order(coords[,1])
yord=order(coords[,2])
for(i in 1:12){
  for(j in 1:12){
    ypred[i,j]=ypred[xord[i],yord[j]]
  }
}

z=ypred
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=sort(coords[,1]),y=sort(coords[,2]),z=ypred,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)


##### Non-Gaussian Case

#Our GLM will be a binary variable on whether or not the temperature is above 70.
set.seed(1)
spglmdat=spmoddat
spglmdat$Temp=rep(0,12)
spglmdat$Temp[which(spmoddat$Temp>70)]=1

weights=rep(1,12)
fit <- glm((Temp/weights)~1, weights=weights, family="binomial",data=spglmdat)
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length
starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0)
tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5)
priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1))
amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43)

m.2=spGLM(Temp~1,weights=weights,data=spglmdat,coords=coords,starting=starting,
          tuning=tuning,priors=priors,cov.model=cov.model,amcmc=amcmc,n.samples=n.samples,n.report=10)
burn.in <- 0.5*n.samples
sub.samps <- burn.in:n.samples

round(summary(window(m.2$p.beta.theta.samples, start=burn.in))$quantiles[,c(1,3,5)],3)

beta.hat <- m.2$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.2$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(as.matrix(rep(1,12))%*%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(12, size=weights, prob=p.hat)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##### Prediction for the GLM version

set.seed(1876)
ypred2=matrix(0,nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    coordsmat=matrix(pred.coords[12*(j-1)+i,],ncol=2)
    m.2.pred=spPredict(m.2,pred.coords = coordsmat,pred.covars=pred.covars,start=burn.in,verbose=FALSE)
    ypred2[i,j]=mean(m.2.pred$p.y.predictive.samples[1:500])
  }
}
for(i in 1:12){
  for(j in 1:12){
    ypred2[i,j]=ypred2[xord[i],yord[j]]
  }
}
z=ypred2
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=sort(coords[,1]),y=sort(coords[,2]),z=ypred2,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)
































