library(ggplot2)
library('infotheo')
source("zeroinflated.R")
n.set <- 100
set.size <- seq(10,109,1) #gene quantities
samplesize <- 100# cell quantities
n.latent<- 10 #potential variables(genes)
a.latent<- sample(1:7,replace=T,size=n.latent+1)
a.individual <- 2
b.latent <-  rep(1,n.latent+1)
b.individual<- 1
pp <- 0.3 #how to choose
p0=0.4
pw=0.3
pro0=0.1
k=15
set.seed(1)
pro=sample(c(0.3,0.2),replace=T,size=k)
set.seed(112)

#1. calculate the mean correlation
COR<-0
for(i in 1:20)
{
  NB = MZINB.set(set.size[i],samplesize,n.latent,a.latent,a.individual,b.latent,b.individual,pp,k,p0,pw,pro0,pro,pro1)
  real<-!sapply(1:set.size[i],function(x) all(t(NB)[,x]==0))
  NB1 <- NB[real,]
  COR[i] <- sum(cor(t(NB1))-diag(dim(NB1)[1]))/(dim(NB1)[1]*(dim(NB1)[1]-1))
}
hist(COR,freq=FALSE)
lines(density(COR),col="red",lty=2)

for(i in 1:100)
{
  k=1
  NB = MZINB.set(set.size[i],samplesize,n.latent,a.latent,a.individual,b.latent,b.individual,pp,k,p0,pw,pro0,pro)
  for(g in 1:(set.size[i]-1))
  {
    for(j in (g+1):set.size[i])
    {
      mii[k]=mi(t(NB)[,g],t(NB)[,j])
      k = k+1
    }
  }
  MI[i]=sum(mii)/(k-1)
}
#what if when one column of NB are all zero?

#plot of 2 genes
set.seed(223)
corev <- NULL

#1. one plot
output<-MZINB.set(set.size[1],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro,pro1)
out<-as.data.frame(t(output))
xx<-function(V4,V3){
  if(V4==0&V3==0)"red"else if((V3==0)&(V4!=0))"blue"else if((V4==0)&(V3!=0))"yellow"else "black"
}
labelcol<-sapply(1:100,function(i) xx(out$V1[i],out$V2[i]))
ggplot(data=out,aes(x=out$V1,y=out$V2,color=labelcol))+geom_point(size=2)+geom_jitter(alpha=0.3)
apply(NB,1,function(x){mean(x==0)})

#Then we change different parameters in functon to find the best expression
#for k: larger k leads to more 0；
NB.k<-function(k)
{
  return(mean(MZINB.set(set.size[1],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,rep(0.35,k),pro1)==0))
}
kk=5:100
k1=sapply(5:100,function(x){NB.k(x)})
ggplot(data.frame(kk,k1),aes(x=kk,y=k1))+geom_point()+geom_smooth()

#for p0
NB.p0<-function(pw)
{
  return(mean(MZINB.set(set.size[1],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)==0))
}
y=0
xx=seq(0.1,0.9,by=0.01)
y=sapply(xx,function(h){NB.p0(h)})
ggplot(data.frame(xx,y),aes(x=xx,y=y))+geom_point()+geom_smooth()

#for pro
NB.pro<- function(p)
{
  return(mean(MZINB.set(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,rep(p,k),pro1)==0))
}
#pro can have many forms, so we 
#for p1, larger p1 leads to larger amount 0
NB.p1<- function(pro)
{
  return(mean(MZINB.set(set.size[i],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,rep(pro,k),pro1)==0))
}
xx=seq(0,1,by=0.01)
y=sapply(xx,function(h){NB.p1(h)})
ggplot(data.frame(xx,y),aes(x=xx,y=y))+geom_point()+geom_smooth()

#correlation matrix(MI,)
corNB.x <-function(pro)#pro0 increase and the meancor increase; pw decrease and meancor increase
{
  for(i in 1:n.set)
  {
    NB = MZINB.set(set.size[i],samplesize,n.latent,a.latent,a.individual,b.latent,b.individual,pp,k,p0,pw,pw1,pro0,rep(pro,k),pro1)
    real<-!sapply(1:set.size[i],function(x) all(t(NB)[,x]==0))
    NB1 <- NB[real,]
    COR[i] <- sum(cor(t(NB1))-diag(dim(NB1)[1]))/(dim(NB1)[1]*(dim(NB1)[1]-1))
  }
  return(mean(COR))
}
xx=seq(0.1,0.9,by=0.02)
y=sapply(xx,function(h){corNB.x(h)})
ggplot(data.frame(xx,y),aes(x=xx,y=y))+geom_point()+geom_smooth()
corNB.k <-function(k)
{
  for(i in 1:n.set)
  {
    NB = MZINB.set(set.size[i],samplesize,n.latent,a.latent,a.individual,b.latent,b.individual,pp,k,p0,pw,pw1,pro0,rep(0.35,k),pro1)
    real<-!sapply(1:set.size[i],function(x) all(t(NB)[,x]==0))
    NB1 <- NB[real,]
    COR[i] <- sum(cor(t(NB1))-diag(dim(NB1)[1]))/(dim(NB1)[1]*(dim(NB1)[1]-1))
  }
  return(mean(COR))
}
kk=5:30
y=sapply(kk,function(h){corNB.k(h)})
ggplot(data.frame(kk,y),aes(x=kk,y=y))+geom_point()+geom_smooth()


truecorNB.k <-function(k)
{
  for(i in 1:n.set)
  {
    NB = MZINB.set(set.size[i],samplesize,n.latent,a.latent,a.individual,b.latent,b.individual,pp,k,p0,pw,pro0,rep(0.3,k),pro1)
    real<-!sapply(1:set.size[i],function(x) all(t(NB)[,x]==0))
    NB1 <- NB[real,]
    COR[i] <- sum(cor(t(NB1))-diag(dim(NB1)[1]))/(dim(NB1)[1]*(dim(NB1)[1]-1))
  }
  return(mean(COR))
}

