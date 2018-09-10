source("zeroinflated.R")
source("common.R")
source("BVZINB4.R")
source("numerical.R")
source("bvnb.R")
library(ggplot2)
library('infotheo')
n.set <- 30
set.size <- rep(15,n.set)#sample(10:25,replace=TRUE,size=n.set)#seq(10,39,1) #gene quantities
samplesize <- 100# cell quantities
n.latent<- 10 #potential variables(genes)
a.latent<- c(8,rep(2,n.latent))
a.individual <- 2
b.latent <-  rep(1,n.latent+1)
a.latent1 <- 0.5
b.individual<- 1
pp <- 0.3 #how to choose
p0=0.3
pw= 0.35
pw1=0.13
pro0=0.05
pro1=0.8
sample.p=500
k=25
pro=sample(c(0.3,0.35),replace=TRUE,size = k)
alph=0.99
MZINB.all <- function(n.set, set.size, samplesize, n.latent, a.latent,b.latent, a.individual,b.individual, pp, k, p0, pw, pw1, pro0, pro,pro1)
{
  result <- list(NULL)
  for(i in 1:n.set)
  {
    dat1 <- MZINB.set(set.size[i],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)
    result[i] <- list(dat1)
  }
  return(result)
}
MZINB.randall <- function(n.set,set.size,samplesize,a.latent1,b,latent,a.individual,b.individual,k,p0,pw,pro0,pro)
{
  result <- list(NULL)
  for(i in 1:n.set)
  {
    dat1 <- MZINB.rand(set.size[i],samplesize,a.latent1,b,latent,a.individual,b.individual,k,p0,pw,pro0,pro)
    result[i] <- list(dat1)
  }
  return(result)
}

MZINB.noncorall <- function(n.set,set.size,samplesize,a.individual,b.individual,pw)
{
  result <- list(NULL)
  for(i in 1:n.set)
  {
    dat1 <- MZINB.noncor(set.size[i],samplesize,a.individual,b.individual,pw)
    result[i] <- list(dat1)
  }
  return(result)
}

cormatrix.all<- function(n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
{
  Result <- list(NULL)
  for(i in 1:n.set)
  {
    A <- MZINB.set(set.size[i],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
    Result[i] <- list(COR.vector(A))
  }
  return(Result)
}

COR.all <- function(n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
{
  meancor <- 0
  meanmi <- 0
  for(i in 1:n.set)
  {
    #set.seed(100*i)
    dat1 <- MZINB.set(set.size[i],samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
    cordat <- COR.vector(dat1)
    meancor[i] <- mean(cordat$cor)
    meanmi[i] <- mean(cordat$mi)
  }
  return(cbind(meancor,meanmi))
}

DAT<- MZINB.all(n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)

cor.p <- function(DAT,n.set,set.size,sample.p,alph)
{
  #using pearson correlation
  #screening
  set.size1 = set.size
  for(i in 1:n.set)
  {
    DAT[[i]]=DAT[[i]][apply(DAT[[i]],1,function(x){mean(x==0)})<alph,]
    set.size1[i]=c(dim(DAT[[i]])[1])
  }
  B = NULL
  p.set <- NULL
  cortest <- NA
  zi <- NULL
  for(i in 1:n.set)
  {
    zi <- rbind(zi,cbind(rep(i,set.size1[i]),seq(1,set.size1[i],1)))
  }
  for(i in 1:n.set)
  {
    dat0 <- DAT[[i]]
    kk = dim(dat0)[1]
    cortest[i] <- sum(cor(t(dat0))-diag(kk))/(kk*(kk-1))
    meancor <- 0
    for(j in 1:sample.p)
    {
      background.gene <- NULL
      index <- sample(1:sum(set.size1),replace = FALSE, size = set.size1[i])
      for(t in 1:set.size1[i])
      {
        background.gene <- rbind(background.gene,DAT[[zi[index,1][t]]][zi[index,2][t],])
      }
      kk <- dim(background.gene)[1]
      meancor[j] <- sum(cor(t(background.gene))-diag(kk))/(kk*(kk-1))
    }
  
    #H0: cor <= 0.25
    p.set[i] <- sum(meancor[!is.na(meancor)]>cortest[i])/length(meancor[!is.na(meancor)])
    A=data.frame(rep(i,sample.p),meancor)
    B = rbind(B,A)
  }
  return(p.set)
}

cor.p.true <- function(DAT,n.set,set.size,sample.p,alph)
{
  #using pearson correlation
  #screening
  set.size1= set.size
  for(i in 1:n.set)
  {
    DAT[[i]]=DAT[[i]][apply(DAT[[i]],1,function(x){mean(x==0)})<alph,]
    set.size1[i]=c(dim(DAT[[i]])[1])
  }
  #B = NULL
  p.set <- NULL
  cortest <- NA
  zi <- NULL
  for(i in 1:n.set)
  {
    zi <- rbind(zi,cbind(rep(i,set.size1[i]),seq(1,set.size1[i],1)))
  }
  for(i in 1:n.set)
  {
    cortest[i] <- truecor.matrix(DAT[[i]])
    meancor <- 0
    for(j in 1:sample.p)
    {
      background.gene <- NULL
      index <- sample(1:sum(set.size1),replace = FALSE, size = set.size1[i])
      for(t in 1:set.size1[i])
      {
        background.gene <- rbind(background.gene,DAT[[zi[index,1][t]]][zi[index,2][t],])
      }
      meancor[j] <- truecor.matrix(background.gene)
    }
    #H0: cor <= 0.25
    p.set[i] <- sum(meancor[!is.na(meancor)]>cortest[i])/length(meancor[!is.na(meancor)])
    #A=data.frame(rep(set.size1[i],sample.p),meancor)
    #B = rbind(B,A)
  }
  return(p.set)
}





for(i in 1:n.set)
{
  print(apply(DAT[[i]],1,function(x){mean(x==0)}))
}


p.set <- cor.p(sample.p,n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro,pro1,alph)[[1]]
p.set 
q.set <- p.adjust(p.set,method="BH",n=length(p.set))
B <- cor.p(sample.p,n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro,pro1,alph)[[2]]
cortest<- cor.p(sample.p,n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro,pro1,alph)[[3]]
colnames(B)<-c("index","meancor")
B[,1]=as.factor(B[,1])
p <- ggplot(data=B,aes(x=index,y=meancor))+geom_boxplot(col="red")
p <- p + geom_point(data=data.frame(levels(B[,1]),cortest),aes(x=levels(B[,1]),y=cortest),col="blue")
p
ggplot(data=data.frame(levels(B[,1]),p.set),aes(x=levels(B[,1]),y=p.set,col=p.set))+geom_point()+geom_hline(aes(yintercept=0.05))

#ggplot(data=data.frame(set.size,p.set),aes(x=set.size,y=p.set,col=p.set))+geom_point()+geom_hline(aes(yintercept=0.05))
