#3. replicate all above functions rep times:
source("allset.R")
library(ggplot2)
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
  return(list(p.set,B,cortest))
}

Rep.set <- function(rep.time,sample.p,n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1,alph)
{
  p.value<- matrix(0,nrow=n.set,ncol=rep.time)
  for(i in 1:rep.time)
  {
    DAT<- MZINB.all(n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)
    p.value[,i] <- cor.p(DAT,n.set,set.size,sample.p,alph)
  }
  return(p.value)
}
q.value <- matrix(0,nrow=n.set,ncol=rep.time)
for(i in 1:rep.time)
{
  q.value[,i] <- p.adjust(p.value[,i],method="BH",n=length(p.value[,i]))
}

rep.time=30
p.value <-Rep.set(rep.time,sample.p,n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1,alph)
p.value
Rep.cor <- function(rep.time,n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
{
  overall <- list(NULL)
  for(i in 1:rep)
  {
    overall[i] <- list(cormatrix.all(n.set,set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro))
  }
  return(overall)
}
xlabel=1:rep
alpha = 0.05
accept.H0 = sapply(1:rep,function(i) sum(p.value[!is.na(p.value[,i]),i]>alpha)/length(!is.na(p.value[,i])))
ggplot(data.frame(xlabel,accept.H0),aes(x=xlabel,y=accept.H0))+geom_point()

qvalue = sapply(1:rep,function(i) sum(q.value[,i]<alpha)/length(q.value[,i]))
ggplot(data.frame(xlabel,qvalue),aes(x=xlabel,y=qvalue,col=qvalue))+geom_point()

