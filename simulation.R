R=0
#negbin
n.sim<-1000
n.set=1000
set.size=rep(20,1000)
samplesize=30
n.latent=10
a.latent <- sample(1:5,replace=T,size=n.latent+1)
a.individual <- rep(2,n.set)
pp <- 0.3
b.latent<- rep(1,n.latent+1)
b.individual<-rep(1,n.set)
#initial value
DATA <- function(set.size,samplesize,n.latent,a.latent,a.individual,b.latent,b.individual,pp)
{
  r.generate=matrix(rgamma(samplesize*(n.latent+1),a.latent,1/b.latent),nrow=n.latent+1,ncol=samplesize,byrow=F)
  tmp <- rgamma(set.size*samplesize,a.individual,1/b.individual)
  y.individual <- matrix(tmp,ncol=samplesize,nrow=set.size,byrow=F)
  a <- matrix(0,ncol=n.latent,nrow=set.size)
  for(i in 1:set.size)
  {
    a[i,]= rbinom(n.latent,1,pp)
  }
  a <- cbind(rep(1,set.size),a)
  y.sum <- y.individual + a %*% r.generate
  #if we want to use different b, then we adjust it now by multipling different scale.
  
  NB <- matrix(rpois(set.size*samplesize,y.sum),nrow=set.size,ncol=samplesize,byrow=F)
  k=set.size
  S<- sum(cor(t(NB))-diag(k))/(k*(k-1))
  return(list(NB,S,cor(t(NB))))
}
COR <- sapply(1:n.set,function(x) DATA(set.size[x],samplesize,n.latent,a.latent,a.individual[x],b.latent,b.individual[x],pp)[[2]])
hh=hist(COR,freq=FALSE)
lines(density(COR),col="red",lty=2)
Result <- lapply(1:n.set,function(x) DATA(set.size[x],samplesize,n.latent,a.latent,a.individual[x],b.latent,b.individual[x],pp)[[1]])

output<- do.call(rbind,Result)




#every genesets' pairwise correlation:

pairwise <- lapply(1:n.set,function(x) DATA(set.size[x],samplesize,n.latent,a.latent,a.individual[x],b.latent,b.individual[x],pp)[[3]])

