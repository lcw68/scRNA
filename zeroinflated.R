BVNB <- function(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp)
{
  set.seed(1)
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

ZIFL <- function(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
{
  #set.size represents the amounts of genes
  #samplesize represents amounts of cells
  r.gene <- rgamma((n.latent+1)*samplesize,a.latent,1/b.latent)
  r.gene <- matrix(r.gene,ncol=samplesize,nrow=n.latent+1,byrow=F)
  s.gene <- rgamma(samplesize*set.size,a.individual,1/b.individual)
  s.gene <- matrix(s.gene,ncol=samplesize,nrow=set.size,byrow=F)
  indexa <- matrix(0,ncol=n.latent,nrow=set.size)
  for(i in 1:set.size)
  {
    indexa[i,]= rbinom(n.latent,1,pp)
  }
  indexa <- cbind(rep(1,set.size),indexa)
  genesum <- s.gene + indexa %*% r.gene
  NB <- matrix(rpois(set.size*samplesize,genesum),nrow=set.size,ncol=samplesize,byrow=F)
  #we first generate the negative binomial matrix
  q<-matrix(0,nrow=set.size,ncol=k)
  for(i in 1:k)
  {
    q[,i]<- rbinom(set.size,1,p0)
  }
  q <- cbind(rep(1,set.size),q)
  #q is the matrix that decide how to choose the binary variables
  z<-matrix(0,ncol=k,nrow=samplesize)
  z0<-rbinom(samplesize,1,pro0)
  w <-matrix(0,ncol=samplesize,nrow=set.size)
  for(i in 1:set.size)
  {
    w[i,]<-rbinom(samplesize,1,pw)
  }
  #w has the similar function as s.gene, it's like basic choice for every gene in case of the accidental abnormity
  for(i in 1:k)
  {
    z[,i] <- rbinom(samplesize,1,pro[i])
  }
  z <- cbind(z0,z)
  #z is the basic binary matrix which will combine the rows to generate the actual 0,1 sequence.
  z.new<-pmin(q%*%t(z)+w,1)
  z.new1 <- 1-z.new
  #z.new1 is the real matrix that in every row, 1 means the gene is founded in the cell and 0 means not.
  for(i in 1:set.size)
  {
    NB[i,]= z.new1[i,]*NB[i,]
  }
  return(NB)
}

if(FALSE){
  library(ggplot2)
  set.size <- 10 #gene quantities
  samplesize <- 100# cell quantities
  n.latent<-10 #potential variables(genes)
  a.latent<- rep(2,n.latent+1)#sample(1:5,replace=T,size=n.latent+1)
  a.individual <- 2
  b.latent <-  rep(1,n.latent+1)
  b.individual<- 1
  pp <- 0.3 #how to choose
  p0=0.3
  pro0=0.1
  k=10
  pro=sample(c(0.3,0.2),replace=T,size=k)
  output<-ZIFL(set.size,samplesize,n.latent,a.latent,b.latent,a.individual,b.individual,pp,k,p0,pw,pro0,pro)
}




