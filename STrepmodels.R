library(spBayes)
library(coda)
library(scatterplot3d)
par(mar=c(5,6,4,4)+.1)


##### Sort the data out!
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

##### Must download the data from: https://www.kaggle.com/sudalairajkumar/daily-temperature-of-major-cities
# and save it to the data frame city_temperature
city_temperature <- read.csv("~/Documents/Lecture Notes/Year 4/Project IV/Final Code/city_temperature.csv")
#####

citytemp2000=city_temperature[which(city_temperature$Year==2000),]
ctaf=citytemp2000[which(citytemp2000$City %in% citiesv),]
ctaf=ctaf[which(ctaf$Month %in% 6:9),]
ctaf=ctaf[which(ctaf$Day %in% 1:5),]

afcities = data.frame(
  'City'=citiesv,
  'Easting'=c(550563,714444,261791,453545,637707,458823,584013,532397,519706,228629,474018,329653)+625000*c(0,1,2,4,3,4,4,1,5,2,5,4),
  'Northing'=c(43119,7503256,6243850,38421,8295178,7129376,8476034,9528267,9251009,482496,995458,3323794)+87500000*c(6,2,0,7,4,2,4,5,5,7,8,10)
) #These are googled values and combined into the same UTC zone

stdat=data.frame(
  'City'=rep(citiesv,each=4),
  'Eastkm'=rep(afcities$Easting/1000,each=4),
  'Northkm'=rep(afcities$Northing/1000,each=4),
  'Month'=rep(6:9,12),
  'Temp'=c(79.8,77.0,75.7,77.5,61.4,51.8,66.1,71.0,57.7,56.3,58.1,58.2,72.4,70.0,71.5,73.2,67.5,59.8,67.1,75.2,68.3,65.8,67.9,73.8,
           64.5,52.8,61.1,68.6,77.8,72.4,71.8,78.2,74.8,72.9,65.0,75.2,78.2,77.3,77.5,76.7,65.9,58.4,59.2,61.1,79.8,79.8,82.7,80.9)
) #A temperature is the first available temperature for that city/month combo

##### Ok nice! Now make the continuous time model.

coords = stdat[,2:4]
starting <- list("beta"= 50, "phi"=3/0.5, "sigma.sq"=2, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
cov.model <- "exponential"
priors <- list("phi.Unif"=list(3/1, 3/0.1),
                "sigma.sq.IG"=list(2, 2),
                "tau.sq.IG"=c(2, 1),
                "beta.norm"=list(50,10))

set.seed(87)

m.3=spSVC(Temp~1,coords=coords,starting=starting,
           tuning=tuning,priors=priors,cov.model=cov.model,n.samples=2000,data=stdat)
plot(m.3$p.theta.samples)

par(mfrow=c(1,2))
plot(as.vector(m.3$p.theta.samples[1:1000,1]),type='l',xlim=c(0,2000),ylab=expression(sigma^2),col=6)
lines(x=1001:2000,y=as.vector(m.3$p.theta.samples[1001:2000,1]),type='l')

m.3=spRecover(m.3,start=1,thin=1,n.omp.threads=1)
plot(as.vector(m.3$p.beta.recover.samples[1:1000]),type='l',xlim=c(0,2000),ylab=expression(beta),col=6)
lines(x=1001:2000,y=as.vector(m.3$p.beta.recover.samples[1001:2000]),type='l')
par(mfrow=c(1,1))

m.3=spRecover(m.3,start=1000,thin=1,n.omp.threads=1)
round(summary(m.3$p.beta.recover.samples)$quant[c(1,3,5)],2)
round(summary(m.3$p.theta.recover.samples)$quant[,c(1,3,5)],2) #Summary of the posterior

##### Looks sweet! Set to work on prediction when you're ready.

cogrid2=function(coords,t){
  n=length(coords[,1])
  final=matrix(0,nrow=n^2,ncol=3)
  final[,1]=rep(coords[,1],n)
  final[,2]=rep(coords[,2],each=n)
  final[,3]=rep(t,n^2)
  return(final)
}

cogrid3 = function(coords3d){
  final=list(matrix(0,nrow=12,ncol=12),matrix(0,nrow=12,ncol=12),matrix(0,nrow=12,ncol=12),matrix(0,nrow=12,ncol=12))
  for(t in 6:9){
    coordst=coords3d[which(coords3d$Month==t),]
    final[[t-5]]=cogrid2(coordst,t)
  }
  mat=rbind(final[[1]],final[[2]],final[[3]],final[[4]])
  return(mat)
}

pred.coords=cogrid3(coords)

pred.covars=matrix(1)

### June surface

par(mfrow=c(2,2))
set.seed(5)
ypredjun=matrix(0,nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    if(j!=i){
      coordsmat=matrix(pred.coords[12*(j-1)+i,],ncol=3)
      m.3.pred=spPredict(m.3,pred.coords = coordsmat,pred.covars=pred.covars,start=1000,verbose=FALSE)
      ypredjun[i,j]=mean(m.3.pred$p.y.predictive.samples)
    }
    else{
      ypredjun[i,i]=mean(m.3$p.y.samples[4*(i-1)+1,])
    }
  }
}
xord=order(coords[which(coords$Month==6),1])
yord=order(coords[which(coords$Month==6),2])
for(i in 1:12){
  for(j in 1:12){
    ypredjun[i,j]=ypredjun[xord[i],yord[j]]
  }
}

z=ypredjun
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=sort(coords[which(coords$Month==6),1]),y=sort(coords[which(coords$Month==6),2]),z=ypredjun,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)

###July surface

set.seed(5)
ypredjul=matrix(0,nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    if(j!=i){
      coordsmat=matrix(pred.coords[12*(j-1)+i+144,],ncol=3)
      m.3.pred=spPredict(m.3,pred.coords = coordsmat,pred.covars=pred.covars,start=1000,verbose=FALSE)
      ypredjul[i,j]=mean(m.3.pred$p.y.predictive.samples)
    }
    else{
      ypredjul[i,i]=mean(m.3$p.y.samples[4*(i-1)+2,])
    }
  }
}
xord=order(coords[which(coords$Month==7),1])
yord=order(coords[which(coords$Month==7),2])
for(i in 1:12){
  for(j in 1:12){
    ypredjul[i,j]=ypredjul[xord[i],yord[j]]
  }
}

z=ypredjul
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=sort(coords[which(coords$Month==7),1]),y=sort(coords[which(coords$Month==7),2]),z=ypredjul,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)

### Aug surface

set.seed(5)
ypredaug=matrix(0,nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    if(j!=i){
      coordsmat=matrix(pred.coords[12*(j-1)+i+144*2,],ncol=3)
      m.3.pred=spPredict(m.3,pred.coords = coordsmat,pred.covars=pred.covars,start=1000,verbose=FALSE)
      ypredaug[i,j]=mean(m.3.pred$p.y.predictive.samples)
    }
    else{
      ypredaug[i,i]=mean(m.3$p.y.samples[4*(i-1)+3,])
    }
  }
}
xord=order(coords[which(coords$Month==8),1])
yord=order(coords[which(coords$Month==8),2])
for(i in 1:12){
  for(j in 1:12){
    ypredaug[i,j]=ypredaug[xord[i],yord[j]]
  }
}

z=ypredaug
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=sort(coords[which(coords$Month==8),1]),y=sort(coords[which(coords$Month==8),2]),z=ypredaug,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)

### Sep surface

set.seed(5)
ypredsep=matrix(0,nrow=12,ncol=12)
for(i in 1:12){
  for(j in 1:12){
    if(j!=i){
      coordsmat=matrix(pred.coords[12*(j-1)+i+144*3,],ncol=3)
      m.3.pred=spPredict(m.3,pred.coords = coordsmat,pred.covars=pred.covars,start=1000,verbose=FALSE)
      ypredsep[i,j]=mean(m.3.pred$p.y.predictive.samples)
    }
    else{
      ypredsep[i,i]=mean(m.3$p.y.samples[4*i,])
    }
  }
}
xord=order(coords[which(coords$Month==8),1])
yord=order(coords[which(coords$Month==8),2])
for(i in 1:12){
  for(j in 1:12){
    ypredsep[i,j]=ypredsep[xord[i],yord[j]]
  }
}

z=ypredsep
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=sort(coords[which(coords$Month==9),1]),y=sort(coords[which(coords$Month==9),2]),z=ypredsep,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)
par(mfrow=c(1,1))


######## Discrete time model

y.t=matrix(nrow=144,ncol=4)
for(j in 1:4){
  y.t[seq(1,144,13),j]=as.numeric(stdat[which(stdat$Month==(j+5)),5])
}
y.t=as.data.frame(y.t)

p=1
N.t=ncol(y.t)
n=nrow(y.t)
starting <- list("beta"=rep(50,N.t*p), "phi"=rep(3/0.5, N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=as.matrix(1))

tuning <- list("phi"=rep(5, N.t)) 

priors <- list("beta.0.Norm"=list(rep(50,p), diag(100,p)),
               "phi.Unif"=list(rep(3/1, N.t), rep(3/0.1, N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
               "sigma.eta.IW"=list(1, diag(0.1,1)))

mods <- lapply(paste(colnames(y.t),1,sep='~'), as.formula)

n.samples <- 2000

set.seed(764800)

m.1 <- spDynLM(mods, data=y.t, coords=pred.coords[1:144,1:2],
               starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
               cov.model="exponential", n.samples=n.samples, n.report=25)

betvec=as.vector(m.1$p.beta.samples[,1])
par(mfrow=c(1,2))
plot(betvec[1:1000],type='l',col=6,xlim=c(0,2000),ylab=expression(beta[1]))
lines(x=1001:2000,y=betvec[1001:2000])
sigvec=as.vector(m.1$p.theta.samples[,1])
plot(sigvec[1:1000],type='l',col=6,xlim=c(0,2000),ylim=c(0,350),ylab=expression(sigma[1]^2))
lines(x=1001:2000,y=sigvec[1001:2000])
par(mfrow=c(1,1))

burn.in=0.5*n.samples

quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}

beta <- apply(m.1$p.beta.samples[burn.in:n.samples,], 2, quant)
beta.0 <- beta[,grep("Intercept", colnames(beta))]
plot(1:N.t, beta.0[1,], pch=19, cex=0.5, xlab="months", ylab="beta.0", ylim=range(beta.0))
arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[3,], length=0.02, angle=90)
arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[2,], length=0.02, angle=90)

theta <- apply(m.1$p.theta.samples[burn.in:n.samples,], 2, quant)
sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
tau.sq <- theta[,grep("tau.sq", colnames(theta))]
phi <- theta[,grep("phi", colnames(theta))]

plot(1:N.t, sigma.sq[1,], pch=19, cex=0.5, xlab="months", ylab="sigma.sq", ylim=range(sigma.sq))
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[3,], length=0.02, angle=90)
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[2,], length=0.02, angle=90)

plot(1:N.t, tau.sq[1,], pch=19, cex=0.5, xlab="months", ylab="tau.sq", ylim=range(tau.sq))
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[3,], length=0.02, angle=90)
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[2,], length=0.02, angle=90) 
#These were not used in the report, but are nice to see anyway

y.hat <- apply(m.1$p.y.samples[,burn.in:n.samples], 1, quant)
y.hat.med <- matrix(y.hat[1,], ncol=N.t)
y.hat.up <- matrix(y.hat[3,], ncol=N.t)
y.hat.low <- matrix(y.hat[2,], ncol=N.t)
y.hat.mean = matrix(apply(m.1$p.y.samples[,burn.in:n.samples], 1, mean),ncol=N.t)

##### Make the surfaces, lad!

coords2d=stdat[4*(1:12),c('Eastkm','Northkm')]

### Jun surface
par(mfrow=c(2,2))
ypredjun=matrix(0,nrow=12,ncol=12)
for (i in 1:12){
  for(j in 1:12){
    ypredjun[i,j] = y.hat.mean[12*(i-1)+j,1]
  }
}
for(i in 1:12){
  for(j in 1:12){
    ypredjun[i,j]=ypredjun[xord[i],yord[j]]
  }
}

z=ypredjun
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=coords2d[xord,1],y=coords2d[yord,2],z=ypredjun,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)

### Jul surface

ypredjul=matrix(0,nrow=12,ncol=12)
for (i in 1:12){
  for(j in 1:12){
    ypredjul[i,j] = y.hat.mean[12*(i-1)+j,2]
  }
}
for(i in 1:12){
  for(j in 1:12){
    ypredjul[i,j]=ypredjul[xord[i],yord[j]]
  }
}

z=ypredjul
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=coords2d[xord,1],y=coords2d[yord,2],z=ypredjul,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)

### Aug surface

ypredaug=matrix(0,nrow=12,ncol=12)
for (i in 1:12){
  for(j in 1:12){
    ypredaug[i,j] = y.hat.mean[12*(i-1)+j,3]
  }
}
for(i in 1:12){
  for(j in 1:12){
    ypredaug[i,j]=ypredaug[xord[i],yord[j]]
  }
}

z=ypredaug
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=coords2d[xord,1],y=coords2d[yord,2],z=ypredaug,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)

### Sep surface

ypredsep=matrix(0,nrow=12,ncol=12)
for (i in 1:12){
  for(j in 1:12){
    ypredsep[i,j] = y.hat.mean[12*(i-1)+j,4]
  }
}
for(i in 1:12){
  for(j in 1:12){
    ypredsep[i,j]=ypredsep[xord[i],yord[j]]
  }
}

z=ypredsep
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(x=coords2d[xord,1],y=coords2d[yord,2],z=ypredsep,
      theta=0,phi=90,col=colors[z.facet.range],expand=0.19,box=FALSE,
      shade=NA,border = 7)
par(mfrow=c(1,1))
































