power <- function(replicate,cc,sample.p,n.set,set.size,samplesize,a.individual,b.individual,k,p0,pw,pro0,pro,alph)
{
  type1 = 0
  #H0 is true, so we generate the data without correlation(NB part)
  DAT<- MZINB.randall(n.set,set.size,samplesize,a.latent1,b,latent,a.individual,b.individual,k,p0,pw,pro0,pro)
  #using pearson correlation
  #screening
  set.size1 = set.size
  for(i in 1:n.set)
  {
    DAT[[i]]=DAT[[i]][apply(DAT[[i]],1,function(x){mean(x==0)})<alph,]
    set.size1[i]=dim(DAT[[i]])[1]
  }
  zi <- NULL
  for(i in 1:n.set)
  {
    zi <- rbind(zi,cbind(rep(i,set.size1[i]),seq(1,set.size1[i],1)))
  }
  cortest <- NULL
  for(i in 1:n.set)
  {
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
    judge = 0
    for(s in 1:replicate)
    {
      datest <- MZINB.rand(set.size1[i],samplesize,a.latent1,b,latent,a.individual,b.individual,k,p0,pw,pro0,pro)
      datest = datest[apply(datest,1,function(x){mean(x==0)})<alph,]
      kk = dim(datest)[1]
      cortest[s] <- sum(cor(t(datest))-diag(kk))/(kk*(kk-1))
      judge[s]= ifelse(cortest[s] > sort(meancor)[sample.p*(1-cc)],1,0)
    }
    type1[i] <- sum(judge)/replicate
  }
  return(type1)
}


powernoncor <- function(replicate,cc,sample.p,n.set,set.size,samplesize,a.individual,b.individual,pw,alph)
{
  type1 = 0
  #H0 is true, so we generate the data without correlation(NB part)
  DAT<- MZINB.noncorall(n.set,set.size,samplesize,a.individual,b.individual,pw)
  #using pearson correlation
  #screening
  set.size1 = set.size
  for(i in 1:n.set)
  {
    DAT[[i]]=DAT[[i]][apply(DAT[[i]],1,function(x){mean(x==0)})<alph,]
    set.size1[i]=c(dim(DAT[[i]])[1])
  }
  zi <- NULL
  for(i in 1:n.set)
  {
    zi <- rbind(zi,cbind(rep(i,set.size1[i]),seq(1,set.size1[i],1)))
  }
  cortest <- NULL
  for(i in 1:n.set)
  {
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
    judge = 0
    for(s in 1:replicate)
    {
      datest <- MZINB.noncor(set.size1[i],samplesize,a.individual,b.individual,pw)
      datest = datest[apply(datest,1,function(x){mean(x==0)})<alph,]
      kk = dim(datest)[1]
      cortest[s] <- sum(cor(t(datest))-diag(kk))/(kk*(kk-1))
      judge[s]= ifelse(cortest[s] > sort(meancor)[sample.p*(1-cc)],1,0)
    }
    type1[i] <- sum(judge)/replicate
  }
  return(type1)
}
type1error <- power(replicate,cc,sample.p,n.set,set.size,samplesize,a.individual,b.individual,pw,alph)
ggplot(data.frame(1:n.set,sort(type1)),aes(x=1:n.set,y=sort(type1)))+geom_point()+geom_hline(aes(yintercept=0.1))
power1 <- function(cc)
{
  return(mean(power(replicate,cc,sample.p,n.set,set.size,samplesize,a.individual,b.individual,pw,alph)))
}

ccc=seq(0,0.3,0.005)
y=sapply(ccc,function(x){power1(x)})
ggplot(data.frame(ccc,y),aes(x=ccc,y=y,col=y))+geom_point()


#realdata part

xv <- sample(1:23425,size=10)
yv <- sample(2:801,size=100)
data.select <- data[xv,yv]
colnames(data.select)<-NULL
rownames(data.select)<-NULL
data.select <- as.matrix(data.select)
apply(data.select,1,function(x){mean(x==0)})
data.select <- data.select[apply(data.select,1,function(x){mean(x==0)})<alph,]
d=dim(data.select)[1]
sum(cor(t(data.select))-diag(d))/(d*(d-1))
truecor.matrix(data.select)
