library(mvtnorm)
library(Boom)
library(scatterplot3d)
set.seed(28)
#Make the parent covariance function
parcov = function(sigsq,phi,tausq,x1,y1,t1,x2,y2,t2){
  if((x1==x2)&&(y1==y2)&&(t1==t2)){
    return(sigsq+tausq)
  }
  d=sqrt((x1-x2)^2+(y1-y2)^2)
  u=abs(t1-t2)
  return(sigsq*exp(-phi*(d+u)))
}

#Set the true values of the 'unknown' parameters
phitrue=0.5
sigmasqtrue=2
tausqtrue = 1

##Generate the synthetic data

#(a) make the frame with space-time locations (just error process)
df=data.frame(
  'xco'=c(0,1,0,1,0,1,0,1),
  'yco'=c(0,0,2,2,0,0,2,2),
  'time'=c(0,0,0,0,1,1,1,1),
  'temperature'=c(0,0,0,0,0,0,0,0)
)

#(b) make the cov matrix for the sample from the parent cov function
Rcov=matrix(nrow=8,ncol=8)
for(i in 1:8){
  for(j in 1:8){
    Rcov[i,j]=parcov(sigsq = sigmasqtrue, phi=phitrue, tausq=tausqtrue,
                     x1=df$xco[i],y1=df$yco[i],t1=df$time[i],
                     x2=df$xco[j],y2=df$yco[j],t2=df$time[j])
  }
}

#(c) sample from the desired distribution, noting that the random error is 
#already included.
temp=rmvnorm(n=1,mean=rep(0,8),sigma=Rcov)
df$temperature=as.vector(temp)

#(d) Make the true data for plots and such
dftrue = df
dftrue$temperature = df$temperature + 20

#Define the neighbour sets
Nl=array(dim=c(2,3,8))
Nl[1,,2]=c(df$xco[1],df$yco[1],df$time[1])
Nl[1,,3]=c(df$xco[1],df$yco[1],df$time[1])
Nl[2,,3]=c(df$xco[2],df$yco[2],df$time[2])
Nl[1,,4]=c(df$xco[2],df$yco[2],df$time[2])
Nl[2,,4]=c(df$xco[3],df$yco[3],df$time[3])
Nl[1,,5]=c(df$xco[1],df$yco[1],df$time[1])
Nl[2,,5]=c(df$xco[2],df$yco[2],df$time[2])
Nl[1,,6]=c(df$xco[2],df$yco[2],df$time[2])
Nl[2,,6]=c(df$xco[5],df$yco[5],df$time[5])
Nl[1,,7]=c(df$xco[3],df$yco[3],df$time[3])
Nl[2,,7]=c(df$xco[5],df$yco[5],df$time[5])
Nl[1,,8]=c(df$xco[4],df$yco[4],df$time[4])
Nl[2,,8]=c(df$xco[7],df$yco[7],df$time[7])

##For DNNGP model fitting, we need to find the matrix K, and hence the vectors a_N(l_i)

#Find the vectors/matrices of parcov applied to the neighbour sets with themselves
#or the original locations. Note that we have different things going on for locations
#1 and 2, as they have a different number of neighbours from the rest.
CNNinv=function(sigsq,phi,tausq,loc){
  if(loc==1){
    return(0)
  }
  if(loc==2){
    return(1/(sigsq+tausq))
  }
  offdiag=parcov(sigsq=sigsq,phi=phi,tausq=tausq,
                 x1=Nl[1,1,loc],y1=Nl[1,2,loc],t1=Nl[1,3,loc],
                 x2=Nl[2,1,loc],y2=Nl[2,2,loc],t2=Nl[2,3,loc])
  m=matrix(data=c(sigsq+tausq,offdiag,offdiag,sigsq+tausq),nrow=2,ncol=2)
  return(solve(m))
}
CNl=function(sigsq,phi,tausq,loc){
  if(loc==1){
    return(0)
  }
  if(loc==2){
    return(parcov(sigsq=sigsq,phi=phi,tausq=tausq,
                  x1=Nl[1,1,2],y1=Nl[1,2,2],t1=Nl[1,3,2],
                  x2=df$xco[2],y2=df$yco[2],t2=df$time[2]))
  }
  p=parcov(sigsq=sigsq,phi=phi,tausq=tausq,
           x1=Nl[,1,loc],y1=Nl[,2,loc],t1=Nl[,3,loc],
           x2=rep(df$xco[loc],times=2),y2=rep(df$yco[loc],times=2),
           t2=rep(df$time[loc],times=2))
  return(p)
}

#Create the F matrix
Fmat=function(sigsq,phi,tausq){
 vec=rep(0,times=8)
 vec[1]=sigsq+tausq
 vec[2]=sigsq+tausq - (1/(sigsq+tausq))*(CNl(sigsq=sigsq,phi=phi,tausq=tausq,loc=2))^2
 for (i in 3:8){
   vec[i] = sigsq+tausq - 
     CNl(sigsq=sigsq,phi=phi,tausq=tausq,loc=i)%*%
     CNNinv(sigsq=sigsq,phi=phi,tausq=tausq,loc=i)%*%
     CNl(sigsq=sigsq,phi=phi,tausq=tausq,loc=i)
 }
 return(diag(vec))
}

#Create the V matrix
Vmat=function(sigsq,phi,tausq){
  mat=diag(rep(1,times=8))
  mat[1,2]=-CNNinv(sigsq,phi,tausq,2)*CNl(sigsq,phi,tausq,2)
  mat[c(1,2),3]=-CNNinv(sigsq,phi,tausq,3)%*%CNl(sigsq,phi,tausq,3)
  mat[c(2,3),4]=-CNNinv(sigsq,phi,tausq,4)%*%CNl(sigsq,phi,tausq,4)
  mat[c(1,2),5]=-CNNinv(sigsq,phi,tausq,5)%*%CNl(sigsq,phi,tausq,5)
  mat[c(2,5),6]=-CNNinv(sigsq,phi,tausq,6)%*%CNl(sigsq,phi,tausq,6)
  mat[c(3,5),7]=-CNNinv(sigsq,phi,tausq,7)%*%CNl(sigsq,phi,tausq,7)
  mat[c(4,7),8]=-CNNinv(sigsq,phi,tausq,8)%*%CNl(sigsq,phi,tausq,8)
  return(mat)
}

#Create the K matrix
Kmat=function(sigsq,phi,tausq){
  inv=t(Vmat(sigsq,phi,tausq))%*%solve(Fmat(sigsq,phi,tausq))%*%Vmat(sigsq,phi,tausq)
  return(solve(inv))
}

##Construct the unnormalised posteriors with the following priors: 
#sigsq: IG(7,12) so mode = 12/(7+1)=3/2 and var = 12^2/(7-1)^2(7-2)=4/5
#phi: U(0,5) as Zhang 2004/BAN14 says to learn about sigsq
#tausq: IG(7,12), same reason as sigsq

#Here, x is the ordered vector of all unknown/unseen stuff (sigsq,phi,tausq,w(l_i))
unnormpost.dnngp=function(x){
  if((x[2]>5)|(x[2]<=0)){
    return(0)
  }
  f=1
  for(i in 1:8){
    f=f*dnorm(df$temperature[i],mean=x[3+i],sd=sqrt(x[3]))
  }
  f=f*dmvnorm(x[4:11],mean=rep(0,times=8),sigma = Kmat(x[1],x[2],x[3]))
  f=f*dinvgamma(x[1],7,12)*dinvgamma(x[3],7,12)
  return(f)
}

#We must make the 'parent' K matrix (C_RR)
CRR=function(sigsq,phi,tausq){
  mat=matrix(nrow=8,ncol=8)
  mat[1,1]=sigsq+tausq
  for(i in 2:8){
    mat[i,i]=sigsq+tausq
    for(j in 1:(i-1)){
      mat[i,j]=parcov(sigsq = sigsq, phi=phi,tausq=tausq,
                      x1=df$xco[i],y1=df$yco[i],t1=df$time[i],
                      x2=df$xco[j],y2=df$yco[j],t2=df$time[j])
      mat[j,i]=parcov(sigsq = sigsq, phi=phi,tausq=tausq,
                      x1=df$xco[i],y1=df$yco[i],t1=df$time[i],
                      x2=df$xco[j],y2=df$yco[j],t2=df$time[j])
    }
  }
  return(mat)
}

#The same as above, but with for the parent GP
unnormpost.parent=function(x){
  if((x[2]>5)|(x[2]<=0)){
    return(0)
  }
  f=1
  for(i in 1:8){
    f=f*dnorm(df$temperature[i],mean=x[3+i],sd=sqrt(x[3]))
  }
  f=f*dmvnorm(x[4:11],mean=rep(0,times=8),sigma = CRR(x[1],x[2],x[3]))
  f=f*dinvgamma(x[1],7,12)*dinvgamma(x[3],7,12)
  return(f)
}

##Now to do a Random Walk Metropolis-Hastings algorithm to sample from the above.
#x should be the starting values for all the unknowns, func the target function,
#samsize the desired sample size, burn the burn-in period, and cova the variance
#for the random walk.

MHRW=function(x,func,samsize,burn,cova){
  samps=x
  for(i in 2:(samsize+burn)){
    prop = rmvnorm(n=1,x,cova)
    if(runif(1)<func(prop)/func(x)){
      x=prop
    }
    samps=rbind(samps,x)
  }
  return(samps)
}

#Choose a variance matrix for the random walk (choose them to be indep.).
#We don't want it to move too quickly relative to how close to the edges it will be.
#For example, the variance of phi is smaller than 0.1 as the true value is 0.5, and
#it mustn't go below zero. I am concerened that 0.1 is too small, especially for ws,
#which I understand less in general.

covvec=rep(0.1,times=11)
covvec[3]=0.05
covar=diag(covvec)

#####The sampling.
#Starting values
start=c(2,2,2,rep(1,times=8))

samps=MHRW(x=start,func=unnormpost.dnngp,samsize=10000,burn=1000,cova=covar)
#This takes 50 seconds (53.09)

samps.parent=MHRW(x=start,func=unnormpost.parent,samsize=10000,burn=1000,cova=covar)
#This takes 30 seconds (32.25), but has less of my subpar coding

#Gives idea of posterior density of tau^2 for each model (similar)
hist(samps.parent[1001:11000,3])
hist(samps[1001:11000,3],50)

#Reveals poor mixing of the MCMC algorithm for phi
plot(samps.parent[,1],type='l')
plot(samps[,1],type='l')

#~"sigmasq*phi is identifyable, but not them individually" (Zhang, 2004). This
#shows that, despite this, we still don't have particularly good convergence, but
#it is as good as tausq, and better than the two before (sigsq, phi), individually.
sigphi=rep(0,length(samps[,1]))
for (i in 1:length(samps[,1])){
  sigphi[i]=samps[i,1]*samps[i,2]
}
plot(sigphi,type='l')

#Some playing around with starting values and covar, and even the true values
#reveals that both models get decent convergence for the ws but that the actual
#parameters are harder to pin down. 

#Importantly, both give similar performance when all else is kept the same.
#Annoyingly, the parent model is faster, but the DNNGP model has to go through 
#many more functions that I have written, which are slightly suboptimal.

###################### Where Next? ######################

#To continue this, I should try to predict the value of the process at three points:
#(a) p: (x,y,t) = (0,1,0)
#(b) q: (x,y,t) = (1,1,0)
#(c) r: (x,y,t) = (0,1,1)
#(d) s: (x,y,t) = (1,1,1)

## Note that 1 is the half way point for the y-coordinate 

#New points
ptnew=data.frame(
  'xco'=c(0,1,0,1),
  'yco'=c(1,1,1,1),
  'time'=c(0,0,1,1)
)

#New point neighbour sets
Nn=array(dim=c(2,3,4))
Nn[1,,1]=as.numeric(df[1,1:3])
Nn[2,,1]=as.numeric(df[3,1:3])
Nn[1,,2]=as.numeric(df[2,1:3])
Nn[2,,2]=as.numeric(df[4,1:3])
Nn[1,,3]=as.numeric(df[5,1:3])
Nn[2,,3]=as.numeric(df[7,1:3])
Nn[1,,4]=as.numeric(df[6,1:3])
Nn[2,,4]=as.numeric(df[8,1:3])

##For kriging, I need to make the a_k vector and the diag(f_k) matrix from my photo

#make a_k vector wise first
Cnninv=function(sigsq,phi,tausq,ptnum){
  offdiag=parcov(sigsq=sigsq,phi=phi,tausq=tausq,
                 x1=Nn[1,1,ptnum],y1=Nn[1,2,ptnum],t1=Nn[1,3,ptnum],
                 x2=Nn[2,1,ptnum],y2=Nn[2,2,ptnum],t2=Nn[2,3,ptnum])
  m=matrix(data=c(sigsq+tausq,offdiag,offdiag,sigsq+tausq),nrow=2,ncol=2)
  return(solve(m))
}

Cnp=function(sigsq,phi,tausq,ptnum){
  p=parcov(sigsq=sigsq,phi=phi,tausq=tausq,
           x1=Nn[,1,ptnum],y1=Nn[,2,ptnum],t1=Nn[,3,ptnum],
           x2=rep(ptnew$xco[ptnum],times=2),y2=rep(ptnew$yco[ptnum],times=2),
           t2=rep(ptnew$time[ptnum],times=2))
  return(p)
}

ashort=function(sigsq,phi,tausq,ptnum){
  return(Cnninv(sigsq,phi,tausq,ptnum)%*%Cnp(sigsq,phi,tausq,ptnum))
}

along=function(sigsq,phi,tausq,ptnum){
  a=rep(0,8)
  if(ptnum==1){
    a[c(1,3)]=ashort(sigsq,phi,tausq,ptnum)
  }
  if(ptnum==2){
    a[c(2,4)]=ashort(sigsq,phi,tausq,ptnum)
  }
  if(ptnum==3){
    a[c(5,7)]=ashort(sigsq,phi,tausq,ptnum)
  }
  if(ptnum==4){
    a[c(6,8)]=ashort(sigsq,phi,tausq,ptnum)
  }
  return(a)
}

akrig = function(sigsq,phi,tausq){
  m=matrix(data=c(along(sigsq,phi,tausq,1),
             along(sigsq,phi,tausq,2),
             along(sigsq,phi,tausq,3),
             along(sigsq,phi,tausq,4)),
           nrow=4,ncol=8)
  return(m)
}

#Then the variances
flist = function(sigsq,phi,tausq){
  f=c(0,0,0,0)
  f[1]=sigsq+tausq - 
    t(Cnp(sigsq,phi,tausq,1))%*%Cnninv(sigsq,phi,tausq,1)%*%Cnp(sigsq,phi,tausq,1)
  f[2]=sigsq+tausq - 
    t(Cnp(sigsq,phi,tausq,2))%*%Cnninv(sigsq,phi,tausq,2)%*%Cnp(sigsq,phi,tausq,2)
  f[3]=sigsq+tausq - 
    t(Cnp(sigsq,phi,tausq,3))%*%Cnninv(sigsq,phi,tausq,3)%*%Cnp(sigsq,phi,tausq,3)
  f[4]=sigsq+tausq - 
    t(Cnp(sigsq,phi,tausq,4))%*%Cnninv(sigsq,phi,tausq,4)%*%Cnp(sigsq,phi,tausq,4)
  return(f)
}

diagfmat= function(sigsq,phi,tausq){
  m=diag(flist(sigsq,phi,tausq))
  return(m)
}

#Then the parts of the matrix from the photo that no one has a record of...
Amat = function(sigsq,phi,tausq){
  m=akrig(sigsq,phi,tausq)%*%Kmat(sigsq,phi,tausq)%*%t(akrig(sigsq,phi,tausq))+diagfmat(sigsq,phi,tausq)
  return(m)
}

Bmat = function(sigsq,phi,tausq){
  m=akrig(sigsq,phi,tausq)%*%Kmat(sigsq,phi,tausq)
  return(m)
}

#Then the parts of the mvnorn
condmnvec=function(sigsq,phi,tausq,wr){
  v=akrig(sigsq,phi,tausq)%*%wr
  return(v)
}

condvarmat = function(sigsq,phi,tausq){
  m=Amat(sigsq,phi,tausq)-akrig(sigsq,phi,tausq)%*%t(Bmat(sigsq,phi,tausq))
  return(m)
}

##Now do the kriging
krigs = function(samples){
  newpreds = matrix(
    as.vector(rmvnorm(1,
            mean=condmnvec(samples[1001,1],samples[1001,2],samples[1001,3],
                           samples[1001,4:11]),
            sigma = condvarmat(samples[1001,1],samples[1001,2],samples[1001,3])
            )
  ), nrow=1,ncol=4)
  
  for(i in 1002:10000){
    next1= rmvnorm(1,
                   mean=condmnvec(samples[i,1],samples[i,2],samples[i,3],
                                  samples[i,4:11]),
                   sigma = condvarmat(samples[i,1],samples[i,2],samples[i,3])
    )
    newpreds=rbind(newpreds,next1)
  }
  return(newpreds)
}

newptdist = krigs(samps)

#And investigate: all nice histograms, but very centred at zero. 
hist(newptdist[,1])
hist(newptdist[,2])
hist(newptdist[,3])
hist(newptdist[,4])
plot(density(newptdist[,1]))

#Design the kriging for the parent function
ATmat=function(sigsq,phi,tausq){
  m=matrix(nrow=4,ncol=4)
  for(i in 1:4){
    for(j in 1:4){
      m[i,j] = parcov(sigsq,phi,tausq,
                      x1=ptnew$xco[i],y1=ptnew$yco[i],t1=ptnew$time[i],
                      x2=ptnew$xco[j],y2=ptnew$yco[j],t2=ptnew$time[j])
    }
  }
  return(m)
}
BTmat=function(sigsq,phi,tausq){
  m=matrix(nrow=4,ncol=8)
  for(i in 1:4){
    for(j in 1:8){
      m[i,j] = parcov(sigsq,phi,tausq,
                      x1=ptnew$xco[i],y1=ptnew$yco[i],t1=ptnew$time[i],
                      x2=df$xco[j],y2=df$yco[j],t2=df$time[j])
    }
  }
  return(m)
}

condmnvec.p=function(sigsq,phi,tausq,wr){
  v=BTmat(sigsq,phi,tausq)%*%solve(CRR(sigsq,phi,tausq))%*%wr
  return(v)
}

condvarmat.p=function(sigsq,phi,tausq){
  m=ATmat(sigsq,phi,tausq)-BTmat(sigsq,phi,tausq)%*%solve(CRR(sigsq,phi,tausq))%*%t(BTmat(sigsq,phi,tausq))
  return(m)
}

krigs.parent=function(samples){
  newpreds = matrix(
    as.vector(rmvnorm(1,
                      mean=condmnvec.p(samples[1001,1],samples[1001,2],samples[1001,3],
                                     samples[1001,4:11]),
                      sigma = condvarmat.p(samples[1001,1],samples[1001,2],samples[1001,3])
    )
    ), nrow=1,ncol=4)
  
  for(i in 1002:10000){
    next1= rmvnorm(1,
                   mean=condmnvec.p(samples[i,1],samples[i,2],samples[i,3],
                                  samples[i,4:11]),
                   sigma = condvarmat.p(samples[i,1],samples[i,2],samples[i,3])
    )
    newpreds=rbind(newpreds,next1)
  }
  return(newpreds)
}

npd.parent = krigs.parent(samps.parent)
20+mean(newptdist[,4])
20+mean(newptdist[,3])
20+mean(newptdist[,2])
20+mean(newptdist[,1])
20+mean(npd.parent[,1])
20+mean(npd.parent[,2])
20+mean(npd.parent[,3])
20+mean(npd.parent[,4]) 
#These means are interesting to compare...they show two neighbours is not enough here!

#Put everything into one data frame
dffull=data.frame(
  'xco'=c(0,1,0,1,0,1,0,1,0,1,0,1),
  'yco'=c(0,0,2,2,0,0,2,2,1,1,1,1),
  'time'=c(0,0,0,0,1,1,1,1,0,0,1,1),
  'temperature'= c(dftrue$temperature,20+mean(newptdist[,1]),20+mean(newptdist[,2]),20+mean(newptdist[,3]),20+mean(newptdist[,4]))
)

#Plot the locations: blue is t=0, and magenta is t=1
colkrig=c(rep(6,4),rep(5,4),6,6,5,5)
col2=c(rep(6,8),rep(5,4))
scatterplot3d(dffull[,1:3],color=col2,pch=c(1,1,1,1,19,19,19,19,1,1,19,19))

dffull.parent=data.frame(
  'xco'=c(0,1,0,1,0,1,0,1,0,1,0,1),
  'yco'=c(0,0,2,2,0,0,2,2,1,1,1,1),
  'time'=c(0,0,0,0,1,1,1,1,0,0,1,1),
  'temperature'= c(dftrue$temperature,20+mean(npd.parent[,1]),20+mean(npd.parent[,2]),20+mean(npd.parent[,3]),20+mean(npd.parent[,4]))
)

par(mfrow=c(1,2))
scatterplot3d(dffull[,c(1,2,4)],color=colkrig+2,type='h',pch=c(1,1,1,1,19,19,19,19,1,1,19,19))
scatterplot3d(dffull.parent[,c(1,2,4)],color=colkrig+2,type='h',pch=c(1,1,1,1,19,19,19,19,1,1,19,19))
par(mfrow=c(1,1)) #Compares the predictions of the two cases

#Prepare the standard persp plots
tempmat=array(dim=c(2,3,2))
tempmat[1,,1]=as.numeric(dffull$temperature[c(1,9,3)])
tempmat[2,,1]=as.numeric(dffull$temperature[c(2,10,4)])
tempmat[1,,2]=as.numeric(dffull$temperature[c(5,11,7)])
tempmat[2,,2]=as.numeric(dffull$temperature[c(6,12,8)])

z=tempmat[,,1]
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

par(mfrow=c(1,4))
for(i in 0:3){
  persp(x=c(0,1),y=c(0,1,2),z=tempmat[,,1],
        theta=0,phi=30*i,col=colors[z.facet.range],expand=0.19,box=FALSE,
        shade=NA,border = 7)
}
par(mfrow=c(1,1)) #This is the t=0 for DNNGP (sequentially tilted)

z=tempmat[,,2]
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

par(mfrow=c(1,4))
for(i in 0:3){
  persp(x=c(0,1),y=c(0,1,2),z=tempmat[,,2],
        theta=0,phi=30*i,col=colors[z.facet.range],expand=0.19,box=FALSE,
        shade=NA,border = 7)
}
par(mfrow=c(1,1)) #This is t=1 for DNNGP, also sequentially tilted.

tempmat2=array(dim=c(2,3,2))
tempmat2[1,,1]=as.numeric(dffull.parent$temperature[c(1,9,3)])
tempmat2[2,,1]=as.numeric(dffull.parent$temperature[c(2,10,4)])
tempmat2[1,,2]=as.numeric(dffull.parent$temperature[c(5,11,7)])
tempmat2[2,,2]=as.numeric(dffull.parent$temperature[c(6,12,8)])

z=tempmat2[,,1]
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

par(mfrow=c(1,4))
for(i in 0:3){
  persp(x=c(0,1),y=c(0,1,2),z=tempmat2[,,1],
        theta=0,phi=30*i,col=colors[z.facet.range],expand=0.19,box=FALSE,
        shade=NA,border = 7)
}
par(mfrow=c(1,1)) #This is t=0 for the full GP method (with tilt)

z=tempmat2[,,2]
col.pal<-colorRampPalette(c(5,6))
colors<-col.pal(100)
# height of facets
z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
# Range of the facet center on a 100-scale (number of colors)
z.facet.range<-cut(z.facet.center, 100)

par(mfrow=c(1,4))
for(i in 0:3){
  persp(x=c(0,1),y=c(0,1,2),z=tempmat2[,,2],
        theta=0,phi=30*i,col=colors[z.facet.range],expand=0.19,box=FALSE,
        shade=NA,border = 7)
}
par(mfrow=c(1,1)) #This is t=1 for the GP method with tilt.

#In the report I use:   plots to show the locations of the data and kriging points;
#                       plots to show the kriging surfaces;
#                       plots to compare the predictions;
#                       means to demonstrate this.




















