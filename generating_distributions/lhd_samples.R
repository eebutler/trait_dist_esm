load("PFT_mean_cov.RData")

library(lhs)
library(mvtnorm)

## generates maximin LHS sample from multivariate normal distribution
## n=number of samples (input points), seed=seed, mu=mean vector, c=covariance matrix

lhd_mult_norm=function(seed,n,mu,c){
  d=length(mu)
  set.seed(seed)
  lhs_unif=maximinLHS(n,d,eps=1e-8)
  lhs_norm=qnorm(lhs_unif)
  cholc=chol(c)
  lhs_mvnorm=t(t(cholc)%*%t(lhs_norm)+mu)
  lhs_mvnorm
}

seed=1
n=100
lhd_samples=lapply(subPFT_summ, function(x) lhd_mult_norm(seed,n,x$mean,x$cov))
lhd_samples=lapply(lhd_samples,exp)

for (i in seq(1,14)){
  lhd_samples[[i]] <- rbind(lhd_samples[[i]],colMeans(lhd_samples[[i]]))
}

write.csv(lhd_samples,'lhd_100.csv',row.names=F)
