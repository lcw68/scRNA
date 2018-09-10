source("zeroinflated.R")
source("common.R")
source("BVZINB4.R")
source("numerical.R")
source("bvnb.R")
library(ggplot2)
library(rootSolve)
n.set=15
sample.p=50
set.size <- rep(10,n.set)#sample(10:25,replace=TRUE,size=n.set)#seq(10,39,1) #gene quantities
samplesize <- 100# cell quantities
n.latent<- 10 #potential variables(genes)
a.latent<- c(0.6,rep(1,n.latent))
a.individual <- 2
b.latent <-  rep(1,n.latent+1)
b.individual<- 1
pp <- 0.35 #how to choose
p0=0.3
pw= 0.3
pw1=0.13
pro0=0.05
pro1=0.8
k=25
pro=sample(c(0.33,0.3),replace=TRUE,size = k)
alph=0.98
MZINB.all <- function(n.set, set.size, samplesize, n.latent, a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)
{
  result <- list(NULL)
  for(i in 1:n.set)
  {
    dat1 <- MZINB.set(set.size[i],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)
    result[i] <- list(dat1)
  }
  return(result)
}


truecor.matrix <- function(NB)
{
  lr <- dim(NB)[1]
  sumcor = 0
  for(i in 1:(lr-1))
  {
    for(j in (i + 1):lr)
    {
      sumcor <-sumcor + ML.BvZINB4.2b(t(NB)[,i],t(NB)[,j],showFlag = FALSE, tol = 1e-8,boosting = FALSE)[12]
    }
  }
  d <- lr*(lr-1)/2
  return(sumcor/d)
}
cor.compare <- function(DAT,n.set,set.size,sample.p,alph)
{
  #using pearson correlation
  #screening
  set.size1 = set.size
  for(i in 1:n.set)
  {
    DAT[[i]]=DAT[[i]][apply(DAT[[i]],1,function(x){mean(x==0)})<alph,]
    set.size1[i]=c(dim(DAT[[i]])[1])
  }
  p.set1 <- NULL
  p.set2 <- NULL
  p.set3 <- NULL
  cortest1 <- NA
  cortest2 <- NA
  cortest3 <- NA
  zi <- NULL
  for(i in 1:n.set)
  {
    zi <- rbind(zi,cbind(rep(i,set.size1[i]),seq(1,set.size1[i],1)))
  }
  for(i in 1:n.set)
  {
    dat0 <- DAT[[i]]
    kk = dim(dat0)[1]
    cortest1[i] <- sum(cor(t(dat0))-diag(kk))/(kk*(kk-1))
    cortest2[i]<- truecor.matrix(DAT[[i]])
    
    meancor1 <- 0
    meancor2 <- 0
    for(j in 1:sample.p)
    {
      background.gene <- NULL
      index <- sample(1:sum(set.size1),replace = FALSE, size = set.size1[i])
      for(t in 1:set.size1[i])
      {
        background.gene <- rbind(background.gene,DAT[[zi[index,1][t]]][zi[index,2][t],])
      }
      kk <- dim(background.gene)[1]
      meancor1[j] <- sum(cor(t(background.gene))-diag(kk))/(kk*(kk-1))
      meancor2[j] <- truecor.matrix(background.gene)
    }
    p.set1[i] <- sum(meancor1>cortest1[i])/length(meancor1)
    p.set2[i] <- sum(meancor2>cortest2[i])/length(meancor2)
  }
  return(cbind(p.set1,p.set2))
}

DAT<- MZINB.all(n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)

ppset <- cor.compare(DAT,n.set,set.size,sample.p,alph)
ggplot(data=data.frame(1:n.set,ppset[,1]),aes(x=1:n.set,y=ppset[,1],col=ppset[,1]))+geom_point()+geom_hline(aes(yintercept=0.05))

ggplot(data=data.frame(1:n.set,ppset[,2]),aes(x=1:n.set,y=ppset[,2],col=ppset[,2]))+geom_point()+geom_hline(aes(yintercept=0.05))
