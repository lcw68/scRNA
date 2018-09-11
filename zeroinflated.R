#1.Negative Binomial
source("common.R")
source("BVZINB4.R")
source("numerical.R")
source("bvnb.R")

MVNB <- function(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp)
{
  r.gene <- rgamma((n.latent+1)*samplesize,a.latent,1/b.latent)
  r.gene <- matrix(r.gene,ncol=samplesize,nrow=n.latent+1,byrow=F)
  s.gene <- rgamma(samplesize*set.size,a.individual,1/b.individual)
  s.gene <- matrix(s.gene,ncol=samplesize,nrow=set.size,byrow=F)
  a <- matrix(0,ncol=n.latent,nrow=set.size)
  for(i in 1:set.size)
  {
    a[i,]= rbinom(n.latent,1,pp)
  }
  a <- cbind(rep(1,set.size),a)
  genesum <- s.gene + a %*% r.gene
  NB <- matrix(rpois(set.size*samplesize,genesum),nrow=set.size,ncol=samplesize,byrow=F)
  #every row is negative binomial
  return(NB)
}#this is a function which generates the multivariate negative binomial.

MVNB.noncor <- function(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,pro1)
{
  r.gene <- rgamma((n.latent+1)*samplesize,a.latent,1/b.latent)
  r.gene <- matrix(r.gene,ncol=samplesize,nrow=n.latent+1,byrow=F)
  s.gene <- rgamma(samplesize*set.size,a.individual,1/b.individual)
  s.gene <- matrix(s.gene,ncol=samplesize,nrow=set.size,byrow=F)
  select <- sample(1:set.size,replace=FALSE,size=set.size*pro1)
  s.gene.cor <- s.gene[select,]
  di <- dim(s.gene.cor)[1]
  indexa <- matrix(0,ncol=n.latent,nrow=di)
  for(i in 1:di)
  {
    indexa[i,]= rbinom(n.latent,1,pp)
  }
  indexa <- cbind(rep(1,di),indexa)
  genesum <- s.gene.cor + indexa %*% r.gene
  NB <- matrix(rpois(di*samplesize,genesum),nrow=di,ncol=samplesize,byrow=F)
  NB.noncor <- matrix(rpois((set.size-di)*samplesize,s.gene[-select,]),nrow=set.size-di,ncol=samplesize,byrow=F)
  NB1=rbind(NB,NB.noncor)
  return(NB1)
}
NB1=MVNB.noncor(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,pro1)
#2. Multivariate zeroinflated negative binomial
#we add some non-association part to every gene set
MZINB.set <- function(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pw1,pro0,pro,pro1)
{
  #set.size represents the amounts of genes
  #samplesize represents amounts of cells
  r.gene <- rgamma((n.latent+1)*samplesize,a.latent,1/b.latent)
  r.gene <- matrix(r.gene,ncol=samplesize,nrow=n.latent+1,byrow=F)
  s.gene <- rgamma(samplesize*set.size,a.individual,1/b.individual)
  s.gene <- matrix(s.gene,ncol=samplesize,nrow=set.size,byrow=F)
  select <- sample(1:set.size,replace=FALSE,size=set.size*pro1)
  s.gene.cor <- s.gene[select,]
  di <- dim(s.gene.cor)[1]
  indexa <- matrix(0,ncol=n.latent,nrow=di)
  for(i in 1:di)
  {
    indexa[i,]= rbinom(n.latent,1,pp)
  }
  indexa <- cbind(rep(1,di),indexa)
  genesum <- s.gene.cor + indexa %*% r.gene
  NB <- matrix(rpois(di*samplesize,genesum),nrow=di,ncol=samplesize,byrow=F)
  NB.noncor <- matrix(rpois((set.size-di)*samplesize,s.gene[-select,]),nrow=set.size-di,ncol=samplesize,byrow=F)
  #we first generate the negative binomial matrix
  q<-matrix(0,nrow=di,ncol=k)
  for(i in 1:k)
  {
    q[,i]<- rbinom(di,1,p0)
  }
  q <- cbind(rep(1,di),q)
  #q is the matrix that decide how to choose the binary variables
  z<-matrix(0,ncol=k,nrow=samplesize)
  w <-matrix(0,ncol=samplesize,nrow=set.size)
  for(i in 1:di)
  {
    w[i,]<-rbinom(samplesize,1,pw)
  }
  for(i in (di+1):set.size)
  {
    w[i,]<-rbinom(samplesize,1,pw1)
  }
  #w has the similar function as s.gene, it's like basic choice for every gene in case of the accidental abnormity
  for(i in 1:k)
  {
    z[,i] <- rbinom(samplesize,1,pro[i])
  }
  z0<-rbinom(samplesize,1,pro0)
  z <- cbind(z0,z)
  #z is the basic binary matrix which will combine the rows to generate the actual 0,1 sequence.
  z.new<-pmin(q%*%t(z)+w[1:di,],1)
  z.new1 <- 1-z.new
  #z.new1 is the real matrix that in every row, 1 means the gene is founded in the cell and 0 means not.
  for(i in 1:di)
  {
    NB[i,]= z.new1[i,]*NB[i,]
  }
  for(i in 1:(set.size[1]-di))
  {
    NB.noncor[i,]=w[i+di,]*NB.noncor[i,]
  }
  NB=rbind(NB,NB.noncor)
  return(NB)
}

# NB part cor = 0; zI part keep the same
MZINB.rand <- function(set.size,samplesize,a.latent1,b,latent,a.individual,b.individual,k,p0,pw,pro0,pro)
{
  #set.size represents the amounts of genes
  #samplesize represents amounts of cells
  r.gene <- rgamma(samplesize,a.latent1,1/b.latent)
  s.gene <- rgamma(samplesize*set.size,a.individual,1/b.individual)
  s.gene <- matrix(s.gene,ncol=samplesize,nrow=set.size,byrow=F)
  s.gene <- s.gene + rep(1,set.size[1])%*%t(r.gene)
  NB <- matrix(rpois(set.size*samplesize,s.gene),nrow=set.size,ncol=samplesize,byrow=F)
  #we first generate the negative binomial matrix
  w <-matrix(0,ncol=samplesize,nrow=set.size)
  for(i in 1:set.size)
  {
    w[i,]<-rbinom(samplesize,1,pw)
  }
  #w has the similar function as s.gene, it's like basic choice for every gene in case of the accidental abnormity
  z0<-rbinom(samplesize,1,pro0)
  #z is the basic binary matrix which will combine the rows to generate the actual 0,1 sequence.
  z.new<-pmin(z0+w,1)
  z.new1 <- 1-z.new
  #z.new1 is the real matrix that in every row, 1 means the gene is founded in the cell and 0 means not.
  for(i in 1:set.size)
  {
    NB[i,]= z.new1[i,]*NB[i,]
  }
  return(NB)
}

MZINB.noncor <- function(set.size,samplesize,a.individual,b.individual,pw)
{
  #set.size represents the amounts of genes
  #samplesize represents amounts of cells
  s.gene <- rgamma(samplesize*set.size,a.individual,1/b.individual)
  s.gene <- matrix(s.gene,ncol=samplesize,nrow=set.size,byrow=F)

  NB <- matrix(rpois(set.size*samplesize,s.gene),nrow=set.size,ncol=samplesize,byrow=F)
  #we first generate the negative binomial matrix
  w <-matrix(0,ncol=samplesize,nrow=set.size)
  for(i in 1:set.size)
  {
    w[i,]<-rbinom(samplesize,1,pw)
  }
  #w has the similar function as s.gene, it's like basic choice for every gene in case of the accidental abnormity
  z.new1 <- 1-w
  #z.new1 is the real matrix that in every row, 1 means the gene is founded in the cell and 0 means not.
  for(i in 1:set.size)
  {
    NB[i,]= z.new1[i,]*NB[i,]
  }
  return(NB)
}


#mutual information
mi <- function(r1,r2)
{
  H<-function(t){#求向量t的熵H(t)
    frequency<-as.data.frame(table(t)) #统计每个元素出现的频数
    f<-frequency$Freq/length(t) #f为频率向量
    return(-sum(f*log(f)))
  }
  return(H(r1)+H(r2)-H(paste(r1,r2)))
}

#give a list including pairwise correlation and mutual information
COR.vector <- function(dat1)
{
  lr <- dim(dat1)[1]
  numcol <- lr*(lr-1)/2
  output <- data.frame(i = rep(NA,numcol),j =rep(NA,numcol), cor = NA,mi = NA, ZINB = NA)
  k <- 1
  for(i in 1:(dim(dat1)[1]-1))
  {
    for(j in (i+1):dim(dat1)[1])
    {
      output[k,] <- c(i,j,cor(t(dat1)[,i],t(dat1)[,j]),mi(t(dat1)[,i],t(dat1)[,j]),NA)
      k = k+1
    }
  }
  return(output)
}

# We use TRUE correlation instead of pearson
purecor.vector <- function(dat1)
{
  lr <- dim(dat1)[1]
  numcol <- lr*(lr-1)/2
  output <- data.frame(i = rep(NA,numcol),j = rep(NA,numcol), purecor = NA)
  for(i in 1:(lr-1))
  {
    for(j in (i + 1):lr)
    {
      output[k,] <- c(i,j,ML.BvZINB4.2b(t(dat1[,i]),t(dat1[,j]),showFlag = FALSE, tol = 1e-8)[12])
      k = k + 1
    }
  }
  return(output)
}

if(FALSE){
  library(ggplot2)
  set.size <- 10 #gene quantities
  samplesize <- 100# cell quantities
  n.latent<-10 #potential variables(genes)
  a.latent<- sample(1:7,replace=T,size=n.latent+1)
  a.individual <- 2
  b.latent <-  rep(1,n.latent+1)
  b.individual<- 1
  pp <- 0.3 #how to choose
  p0=0.3
  pro0=0.2
  pw=0.3
  k=15
  pro=sample(c(0.3,0.2),replace=T,size=k)
  output<-MZINB.set(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
}

apply(NB1,1,function(x){mean(x==0)})
apply(NB2,1,function(x){mean(x==0)})


