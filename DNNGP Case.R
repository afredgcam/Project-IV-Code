library(spNNGP)
library(spBayes)
library(coda)
par(mar=c(5,6,4,4)+.1)

##### NNGP with straight line cities

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

afcities = data.frame(
  'City'=citiesv,
  'Easting'=c(550563,714444,261791,453545,637707,458823,584013,532397,519706,228629,474018,329653)+625000*c(0,1,2,4,3,4,4,1,5,2,5,4),
  'Northing'=c(43119,7503256,6243850,38421,8295178,7129376,8476034,9528267,9251009,482496,995458,3323794)+87500000*c(6,2,0,7,4,2,4,5,5,7,8,10)
) #This was Googled for each city and the UTM zones have been merged

stdat=data.frame(
  'City'=rep(citiesv,each=4),
  'Eastkm'=rep(afcities$Easting/1000,each=4),
  'Northkm'=rep(afcities$Northing/1000,each=4),
  'Month'=rep(6:9,12),
  'Temp'=c(79.8,77.0,75.7,77.5,61.4,51.8,66.1,71.0,57.7,56.3,58.1,58.2,72.4,70.0,71.5,73.2,67.5,59.8,67.1,75.2,68.3,65.8,67.9,73.8,
           64.5,52.8,61.1,68.6,77.8,72.4,71.8,78.2,74.8,72.9,65.0,75.2,78.2,77.3,77.5,76.7,65.9,58.4,59.2,61.1,79.8,79.8,82.7,80.9)
) #The Temp part is what has come from city_temperature, according to the month and city. Values are the first
# available temperature in that month.

nngpdf=stdat[which(stdat$City%in%citiesv[c(4,6,7,9,11,12)]),]
nngpdf=nngpdf[,c(1,4,3,5)]
northun=nngpdf$Northkm[4*(0:5)+1]

##### Make the model
set.seed(2000)
n.samples=2000
starting=list("phi"=3/0.5, "sigma.sq"=5, "tau.sq"=1)
tuning=list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
priors=list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
cov.model='exponential'
coords=as.matrix(nngpdf[,c(2,3)])
coords[,1]=coords[,1]*100000
m.s=spNNGP(Temp~1, coords=coords, starting=starting, method="latent", n.neighbors=5,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=1, data=nngpdf)

summary(m.s)[,c(1,3,5)]
par(mfrow=c(1,2))
bet=as.vector(m.s$p.beta.samples)
plot(bet[1:1000],type='l',col=6,ylab=expression(beta),xlim=c(0,2000))
lines(x=1001:2000,bet[1001:2000])

sig=as.vector(m.s$p.theta.samples[,1])
plot(sig[1:1000],type='l',col=6,ylab=expression(sigma^2),xlim=c(0,2000))
lines(x=1001:2000,sig[1001:2000])
par(mfrow=c(1,1))

wint=apply(m.s$p.w.samples[,1001:2000],1,median)
bint=rep(median(m.s$p.beta.samples[1001:2000]),24)
yint=wint+bint
yintjun=yint[4*(0:5)+1]

##### Make a surface plot from this

yint=yint[order(-coords[,2])]

surf=matrix(0,nrow=6,ncol=4)
for (i in 1:6){
  surf[i,]=yint[(4*(i-1)+1):(4*i)]
}

z=surf
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

persp(y=(6:9),x=sort(northun),z=surf,
      theta=90,phi=90,col=colors[z.facet.range],expand=0.19,box=F,
      shade=NA,border = 7) #This has time increasing from left to right, and moving North from bottom to top. 



























