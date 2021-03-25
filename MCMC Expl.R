set.seed(1)

#Make the data 
n=100
mutrue=0.5
dat=rnorm(n,0.5,1)

#Make the target function and the met.-Hast. random walk function
target=function(mu){
  return(exp(-((n+1)/2)*(mu-((n+1)/n)*mean(dat))^2))
}
MHRW1D=function(x,func,samsize,burn,cova){
  samps=x
  for(i in 2:(samsize+burn)){
    prop = rnorm(n=1,x,cova)
    if(runif(1)<func(prop)/func(x)){
      x=prop
    }
    samps=rbind(samps,x)
  }
  return(samps)
}

#Set the inputs for it
samsize=1000
burn=1000
cova=matrix(0.1)
mu0=c(3)

#Get the MCMC Samples, with and without burn-in
samples=MHRW1D(mu0,target,samsize,burn,cova)
gsamps=sort(samples[1001:2000])

#Check for good MCMC properties
plot(samples,pch=19,col=c(rep(6,1000),rep(1,1000)),cex=0.01,main='Trace Plot of MCMC Samples',xlab='',ylab='')
lines(1:1000,samples[1:1000],col=6,type='l')
lines(1001:2000,samples[1001:2000],col=1,type='l')
legend(1000,2.8,legend=c('Burn-In','Samples'),col=c(6,1),lty=1,cex=1.25)

#Plot empirical density vs true density
xpts=seq(0,1,0.01)
plot(density(gsamps),xlab=expression(mu),main=expression(paste('Empirical and True Densities of Posterior Distribution for ', mu)))
lines(xpts,dnorm(xpts,((n+1)/n)*mean(dat),(n+1)^-0.5),col=6)
legend(0.3,3.5,legend=c('Empirical','True'),col=c(1,6),cex=1,lty=1)

#95% HPD credible set for mu
c(gsamps[.025*samsize],gsamps[.975*samsize])
