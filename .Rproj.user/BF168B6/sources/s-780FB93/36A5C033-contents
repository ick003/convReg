#### Functions used for EDA and CPP ####

distribution.list <- function()
{
   # return the list of supported distribution
   return ( c('Pois','Nbinom','CoMPoisson','ZIP','HP','Binom','Multinom','Gauss','Lnorm', 'Gamma') )
}

##### Useful functions ####

#' Zero-truncated Poisson distribution
#' @param N sample size
#' @param lambda mean Poisson parameter
#' @keywords truncated Poisson
#' @export
rtpois <- function(N, lambda){qpois(runif(N, dpois(0, lambda), 1), lambda)}

pkolmogorov1x <- function(x, n) {
   if (x <= 0){return(0)}
   if (x >= 1){return(1)}
   j <- seq.int(from = 0, to = floor(n * (1 - x)))
   xj = x + j/n
   xj[which.min(xj)] <- max(min(xj), 0.0001)
   xj[which.max(xj)] <- min(max(xj), 0.9999)
   1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - (xj)) + (j - 1) * log(xj)))
}


#' Kolmogorov-Smirnov distribution
#' @param y sample values
#' @param n degree of freedom
#' @export
dks = function(y,n){
   pks1 = sapply(y,function(x) pkolmogorov1x(x,n))
   pks2 = pks1[2:length(pks1)]
   return(pks2 - pks1[1:(length(pks1)-1)])
}

dZ.lam = function(lambda,nu){
   dZ = colSums(matrix(unlist(lapply(1:99,function(k){exp(log(k)+ log(lambda)*k - nu*sum(log(1:k)))})), nrow=99, byrow=T))
   dZ[lambda==0]=0
   dZ[which(is.nan(dZ))] = 10^100
   dZ[which(dZ > 10^100)] = 10^100
   dZ[which(dZ < -10^100)] = -10^100
   return(dZ)

}

Zcomp = function(lambda,nu){
   nu[nu==Inf] = 10^6
   lambda[lambda == Inf] = 10^6
   Z = colSums(matrix(unlist(lapply(0:99,function(k){exp( log(lambda)*k - nu*log(factorial(k)))})), nrow=100, byrow=T))
   Z[lambda==0]=0
   Z[which(Z > 10^100)] = 10^100
   Z[which(is.nan(Z))] = 10^100
   Z[which(Z < -10^100)] = -10^100
   return(Z)

}

Zcomp.mu = function(lambda,nu){
   nu[nu==Inf] = 10^6
   lambda[lambda == Inf] = 10^6
   Z.mu = colSums(matrix(unlist(lapply(1:99,function(k){exp(log(k) + log(lambda)*k - nu*sum(log(1:k)))})), nrow=99, byrow=T))/Zcomp(lambda,nu)
   Z.mu[lambda==0]=0
   Z.mu[which(is.nan(Z.mu))] = 10^100

   return(Z.mu)
}

dZ.nu = function(lambda,nu){
   dZ= colSums(matrix(unlist(lapply(0:150,function(k){log(factorial(k))* lambda^k * ((factorial(k))^-nu)})), nrow=151, byrow=T))
   dZ[which(dZ > 10^100)] = 10^100
   dZ[which(is.nan(dZ))] = 10^100
   dZ[which(dZ < -10^100)] = -10^100
   return(dZ)

}

dcom2 = function(x,lambda,nu){

   #x.t = matrix(unlist(lapply(0:100,function(k){lambda^k * ((factorial(k))^-nu)})), nrow=101, byrow=T)

   return(lambda^x * (factorial(x))^-nu)

}

.dCoMPoissonLnorm <- function(x,mu,size,mug,sigmag,log=F){

   if(length(mu)==1){
      res = dlnorm(matrix(kronecker(x,0:100,FUN="-"), ncol=101, byrow=T),mug,sigmag) %*% dcmp(0:100,nu=size,lambda=mu)
   }else{
      rd.p = 0:100
      yx = cbind(x,mu,size,mug,sigmag)
      res = apply(yx,1,function(x) sum(dlnorm(x[1]-rd.p,x[4], x[5])*dcmp(rd.p,nu=x[3],lambda = x[2])))
      res = sapply(res,function(x) max(x,10^(-300)))
   }
   res[is.nan(res)] = 10^(-300)
   if(log==T){res=log(res)}
   return(res)

}

dTCoMPoissonGauss <- function(x,mu,size,mug,sigmag,k.max,f.k,j){
   rd.p = 0:k.max
   yx = cbind(x,mu,size,mug,sigmag)
   if(j == 1){res = apply(yx,1,function(x) sum(dnorm(x[1]-rd.p,x[4], x[5])*dcmp(rd.p,nu=x[3],lambda = x[2])*rd.p))}
   if(j == 2){res = apply(yx,1,function(x) sum(dnorm(x[1]-rd.p,x[4], x[5])*dcmp(rd.p,nu=x[3],lambda = x[2])*log(factorial(rd.p))))}
   if(j ==3){res = apply(yx,1,function(x) sum(dnorm(x[1]-rd.p,x[4], x[5])*dcmp(rd.p,nu=x[3],lambda = x[2])*(x[1]-rd.p)))}
   if(j ==4){res = apply(yx,1,function(x) sum(dnorm(x[1]-rd.p,x[4], x[5])*dcmp(rd.p,nu=x[3],lambda = x[2])*(x[1]-rd.p)^2))}

   res1 = sapply(res,function(x) max(x,10^(-300)))
   res1[which(is.nan(res1))]=10^(100)

   return(res1)
}

.dkCoMPoissonGauss <- function(x,mu,size,mug,sigmag,k){
   if(length(k)==1){
      res=dnorm(x-k,mug, sigmag)*dcomp(k,nu=size,lam = mu)
   }
   if(length(k)>1){
      res = matrix(NA,nrow=length(x), ncol = length(k))
      for(j in 1:length(x)){
         dd = dcom2(0:99,mu[j],size[j])
         res[j,] = dnorm(x[j]-k,mug[j], sigmag[j])*dd[k+1]/sum(dd)
      }
   }
   res[is.nan(res)] = 0
   return(res)
}

.dkCoMPoissonLnorm <- function(x,mu,size,mug,sigmag,k){
   if(length(k)==1){
      res=dlnorm(x-k,mug, sigmag)*dcomp(k,nu=size,lam = mu)
   }
   if(length(k)>1){
      res = matrix(NA,nrow=length(x), ncol = length(k))
      for(j in 1:length(k)){
         res[,j] = dlnorm(x-k[j],mug, sigmag)*dcomp(k[j],nu=size,lam = mu)
      }
   }
   res[which(is.nan(res))] = 10^(-300)
   return(res)
}

.dkBinomLnorm <- function(x,size,prob,mug,sigmag,k){
   if(length(k)==1){
      res=dlnorm(x-k,mug, sigmag)*dbinom(k,size=size,prob = prob)
   }
   if(length(k)>1){
      res = matrix(NA,nrow=length(x), ncol = length(k))
      for(j in 1:length(k)){
         res[,j] = dlnorm(x-k[j],mug, sigmag)*dbinom(k[j],size=size,prob = prob)
      }
   }
   res[which(is.nan(res))] = 10^(-300)
   return(res)
}

.dCoMPoissonGauss <- function(x,mu,size,mug,sigmag,log=F){
   if(length(mu)==1){
      res = dnorm(matrix(kronecker(x,0:100,FUN="-"), ncol=101, byrow=T),mug,sigmag) %*% dcmp(0:100,nu=size,lambda=mu)
   }else{
      res = rep(0,length(x))
      for(i in 1:length(x)){
         dd = dcom2(0:99,mu[i],size[i])
         res[i] = sum(dnorm(x[i]-0:99-mug[i],0,sigmag[i])*dd/sum(dd))
      }
      res[is.nan(res)] = 0
   }
   res = sapply(res,function(x) max(x,10^(-300)))
   if(log==T){res=log(res)}
   return(res)

}

.dkZIPGauss <- function(x,mu,size,mug,sigmag,k){
   if(length(k)==1){
      res=dnorm(x-k,mug, sigmag)*dZIP(k,pi=mu,lam = size)
   }
   if(length(k)>1){
      res = matrix(NA,nrow=length(x), ncol = length(k))
      for(j in 1:length(x)){
         dd = dZIP(k,pi=mu[j],lam = size[j])
         res[j,] = dnorm(x[j]-k,mug[j], sigmag[j])*dd[k+1]/sum(dd)
      }
   }
   res[is.nan(res)] = 0
   return(res)
}

############################### Gauss ###############################

#### Pois / Gauss ####

#' Poisson and Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rPoisGauss <- function(N,mu,size,mug,sigmag){
   rnb.sim = rpois(N,mu)
   LOS.sim = rnb.sim + rnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' Poisson and Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pPoisGauss <- function(q,mu,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pnorm(q-rd.p,mug, sigmag)*dpois(rd.p, mu))
   }
   if(length(q)>1){
      x=q
      res = sapply(x,function(x) sum(pnorm(x-rd.p,mug,sigmag)*dpois(rd.p, mu)))
   }
   return(res)
}

#' Poisson and Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qPoisGauss <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rPoisGauss(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rPoisGauss(100000,mu[i],sigma[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

##### Nbinom / Gauss #####

#' Negative binomial and Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rNbinomGauss <- function(N,mu,size,mug,sigmag){
   rnb.sim = rnbinom(N,size=size,mu=mu)
   LOS.sim = rnb.sim + rnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' Negative binomial and Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pNbinomGauss <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pnorm(q-rd.p,mug, sigmag)*dnbinom(rd.p,size=size,mu = mu))
   }
   if(length(q)>1){
         res = sapply(q,function(x) sum(pnorm(x-rd.p,mug,sigmag)*dnbinom(rd.p,size=size,mu = mu)))
   }
   return(res)
}

#' Negative binomial and Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qNbinomGauss <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
   sample = rNbinomGauss(100000,mu,size,mug,sigmag)
   q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rNbinomGauss(100000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### CoMP / Gauss ####

#' CoMP and Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rCoMPoissonGauss <- function(N,mu,size,mug,sigmag){
   if(length(mu)==1){
      rnb.sim = rcmp(N,nu=size,lambda=mu)
   }else{
      rnb.sim = apply(cbind(mu,size), 1, function(x) rcmp(N/length(mu),lambda = x[1],nu=x[2]))
   }
   LOS.sim = c(t(rnb.sim)) + rnorm(N,mug,sigmag)


   return(LOS.sim)

}

#' CoMP and Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pCoMPoissonGauss <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pnorm(q-rd.p,mug, sigmag)*dcomp(rd.p,nu=size,lam = mu,sumTo = 50), na.rm=T)
   }
   if(length(q)>1){
      x=q
      res = sapply(x,function(x) sum(pnorm(x-rd.p,mug,sigmag)*dcomp(rd.p,nu=size,lam = mu,sumTo = 50), na.rm=T))
   }
   return(res)

}

#' CoMP and Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qCoMPoissonGauss <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rCoMPoissonGauss(10000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rCoMPoissonGauss(10000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### ZIP / Gauss ####

#' ZIP and Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rZIPGauss <- function(N,mu,size,mug,sigmag){
   zip = rbinom(N,1,size)
   rnb.sim = rpois(N,mu)
   rnb.sim[which(zip == 1)] = 0
   LOS.sim = rnb.sim + rnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' ZIP and Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pZIPGauss <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pnorm(q-rd.p,mug, sigmag)*dZIP(rd.p, mu, size))
   }
   if(length(q)>1){
      x=q
      res = sapply(x,function(x) sum(pnorm(x-rd.p,mug,sigmag)*dZIP(rd.p, mu, size)))
   }
   return(res)
}

#' ZIP and Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qZIPGauss <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rZIPGauss(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rZIPGauss(100000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### HP / Gauss ####

#' HP and Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param pi parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rHPGauss <- function(N,mu,pi,mug,sigmag){
   hp = rbinom(N,1,pi)
   rnb.sim = rtpois(N,mu)
   rnb.sim[which(hp == 1)] = 0
   LOS.sim = rnb.sim + rnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' HP and Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pHPGauss <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pnorm(q-rd.p,mug, sigmag)*dHP(rd.p, mu, pi = size))
   }
   if(length(q)>1){
      x=q
      res = sapply(x,function(x) sum(pnorm(x-rd.p,mug,sigmag)*dHP(rd.p, mu, pi = size)))
   }
   return(res)
}

#' HP and Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qHPGauss <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rHPGauss(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rHPGauss(100000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### Binom / Gauss ####

#' Binomial and Gaussian convolution
#' @param N sample size
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rBinomGauss <- function(N,prob,size,mug,sigmag){
   rb.sim = rbinom(N,size=size,prob=prob)
   LOS.sim = rb.sim + rnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' Binomial and Gaussian convolution
#' @param q quantiles
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pBinomGauss <- function(q,prob,size,mug,sigmag){

   rd.p = 0:max(20,size)
   if(length(q)==1){
      res = sum(pnorm(q-rd.p,mug, sigmag)*dbinom(rd.p,size=size,prob = prob))
   }

   if(length(q)>1){
     res = sapply(q,function(x) sum(pnorm(x-rd.p,mug,sigmag)*dbinom(rd.p,size=size,prob = prob)))
   }
   return(res)
}

#' Binomial and Gaussian convolution
#' @param p probabilities
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qBinomGauss <- function(p,prob,size,mug,sigmag){

   if(length(size)==1){
      sample = rBinomGauss(100000,size=size,prob=prob,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(size)>1){
      q = NULL
      for(i in 1:length(size)){
         sample = rBinomGauss(100000,size=size[i],prob=prob[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }

   return(q)

}

############################## Lnorm ################################

#### Pois / Lnorm ####

#' Poisson and Log-Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rPoisLnorm <- function(N,mu,size,mug,sigmag){
   rnb.sim = rpois(N,mu)
   LOS.sim = rnb.sim + rlnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' Poisson and Log-Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pPoisLnorm <- function(q,mu,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(plnorm(q-rd.p,mug, sigmag)*dpois(rd.p, mu))
   }
   if(length(q)>1){
     x=q
     res = sapply(x,function(x) sum(plnorm(x-rd.p,mug,sigmag)*dpois(rd.p, mu)))
   }
   return(res)

}

#' Poisson and Log-Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qPoisLnorm <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rPoisLnorm(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rPoisLnorm(100000,mu[i],sigma[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### Nbinom / Lnorm ####

#' Negative binomial and Log-Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rNbinomLnorm <- function(N,mu,size,mug,sigmag){
   rnb.sim = rnbinom(N,size=size,mu=mu)
   LOS.sim = rnb.sim + rlnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' Negative binomial and Log-Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pNbinomLnorm <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(plnorm(q-rd.p,mug, sigmag)*dnbinom(rd.p,size=size,mu = mu))
   }
   if(length(q)>1){
      res = sapply(q,function(x) sum(plnorm(x-rd.p,mug,sigmag)*dnbinom(rd.p,size=size,mu = mu)))
   }
   return(res)
}

#' Negative binomial and Log-Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qNbinomLnorm <- function(p,mu,size,mug,sigmag){
   if(length(mu)==1){
      sample = rNbinomLnorm(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p)
   }
   if(length(mu)>1){
      q=NULL
      for(i in 1:length(mu)){
         sample = rNbinomLnorm(100000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p))
      }

   }
   return(q)

}

#### CoMP / Lnorm ####

#' CoMP and Log-Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rCoMPoissonLnorm <- function(N,mu,size,mug,sigmag){
   if(length(mu)==1){
      rnb.sim = rcmp(N,nu=size,lambda=mu)
   }else{
      rnb.sim = apply(cbind(mu,size), 1, function(x) rcmp(N/length(mu),lambda = x[1],nu=x[2]))
   }
   LOS.sim = c(t(rnb.sim)) + rlnorm(N,mug,sigmag)


   return(LOS.sim)

}

#' CoMP and Log-Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pCoMPoissonLnorm <- function(q,mu,size,mug,sigmag){
  l=1.5
  rd.p = 0:200
  if(length(q)==1){
    res = sum(plnorm(q-rd.p,mug, sigmag)*dcomp(rd.p,nu=size,lam = mu,sumTo = 50), na.rm=T)
  }
  if(length(q)>1){
    x=q
    res = sapply(x,function(x) sum(plnorm(x-rd.p,mug,sigmag)*dcomp(rd.p,nu=size,lam = mu,sumTo = 50), na.rm=T))
  }
  return(res)
}

#' CoMP and Log-Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qCoMPoissonLnorm <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rCoMPoissonLnorm(10000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rCoMPoissonLnorm(10000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### ZIP / Lnorm ####

#' ZIP and Log-Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rZIPLnorm <- function(N,mu,size,mug,sigmag){
   zip = rbinom(N,1,size)
   rnb.sim = rpois(N,mu)
   rnb.sim[which(zip == 1)] = 0
   LOS.sim = rnb.sim + rlnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' ZIP and Log-Gaussian convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pZIPLnorm <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(plnorm(q-rd.p,mug, sigmag)*dZIP(rd.p, mu, size))
   }
   if(length(q)>1){
      x=q
      res = sapply(x,function(x) sum(plnorm(x-rd.p,mug,sigmag)*dZIP(rd.p, mu, size)))
   }
   return(res)
}

#' ZIP and Log-Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qZIPLnorm <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rZIPLnorm(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rZIPLnorm(100000,mu[i],sigma[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### Binom / Lnorm ####

#' Binomial and Log-Gaussian convolution
#' @param N sample size
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rBinomLnorm <- function(N,prob,size,mug,sigmag){
   rb.sim = rbinom(N,size=size,prob=prob)
   LOS.sim = rb.sim + rlnorm(N,mug,sigmag)
   return(LOS.sim)
}

#' Binomial and Log-Gaussian convolution
#' @param q quantiles
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pBinomLnorm <- function(q,prob,size,mug,sigmag){

  rd.p = 0:max(20,size)
  if(length(q)==1){
    res = sum(plnorm(q-rd.p,mug, sigmag)*dbinom(rd.p,size=size,prob = prob))
  }

  if(length(q)>1){
    res = sapply(q,function(x) sum(plnorm(x-rd.p,mug,sigmag)*dbinom(rd.p,size=size,prob = prob)))
  }
  return(res)
}

#' Binomial and Log-Gaussian convolution
#' @param p probabilities
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qBinomLnorm <- function(p,prob,size,mug,sigmag){

   if(length(size)==1){
      sample = rBinomLnorm(100000,size,prob,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(size)>1){
      q = NULL
      for(i in 1:length(size)){
         sample = rBinomLnorm(100000,size[i],prob[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }

   return(q)

}

############################## Gamma ################################

#### Pois / Gamma ####

#' Poisson and Gamma convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rPoisGamma <- function(N,mu,size,mug,sigmag){
   rnb.sim = rpois(N,mu)
   LOS.sim = rnb.sim + rgamma(N,mug,sigmag)
   return(LOS.sim)
}

#' Poisson and Gamma convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pPoisGamma <- function(q,mu,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      x=seq(0,q,0.1)
      res = lapply(x,function(x) sum(dgamma(x-rd.p,mug, sigmag)*dpois(rd.p,mu)))
      #res = sum(unlist(lapply(res,function(x) max(x,10^(-300))))*0.1)
      res = sum(sapply(res,function(x) max(x,10^(-300)))*0.1)
   }
   if(length(q)>1){
      res = NULL
      for(j in 1:length(q)){
         x=seq(0,q[j],0.1)
         res.i = lapply(x,function(x) sum(dgamma(x-rd.p,mug,sigmag)*dpois(rd.p,mu)))
         res.i = sum(unlist(lapply(res.i,function(x) max(x,10^(-300))))*0.1)
         res.i = sum(sapply(res.i,function(x) max(x,10^(-300)))*0.1)
         res = c(res,res.i)
      }
   }
   return(res)

}

#' Poisson and Gamma convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qPoisGamma <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rPoisGamma(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rPoisGamma(100000,mu[i],sigma[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### Nbinom / Gamma ####

#' Negative binomial and Gamma convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rNbinomGamma <- function(N,mu,size,mug,sigmag){
   rnb.sim = rnbinom(N,size=size,mu=mu)
   LOS.sim = rnb.sim + rgamma(N,mug,sigmag)
   return(LOS.sim)
}

#' Negative binomial and Gamma convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pNbinomGamma <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pgamma(q-rd.p,mug, sigmag)*dnbinom(rd.p,size=size,mu = mu))
   }
   if(length(q)>1){
      res = sapply(q,function(x) sum(pgamma(x-rd.p,mug,sigmag)*dnbinom(rd.p,size=size,mu = mu)))
   }
   return(res)
}

#' Negative binomial and Gamma convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qNbinomGamma <- function(p,mu,size,mug,sigmag){
   if(length(mu)==1){
      sample = rNbinomGamma(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p)
   }
   if(length(mu)>1){
      q=NULL
      for(i in 1:length(mu)){
         sample = rNbinomGamma(100000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p))
      }

   }
   return(q)

}

#### CoMP / Gamma ####

#' CoMP and Gamma convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rCoMPoissonGamma <- function(N,mu,size,mug,sigmag){
   if(length(mu)==1){
      rnb.sim = rcmp(N,nu=size,lambda=mu)
   }else{
      rnb.sim = apply(cbind(mu,size), 1, function(x) rcmp(N/length(mu),lambda = x[1],nu=x[2]))
   }
   LOS.sim = c(t(rnb.sim)) + rgamma(N,mug,sigmag)


   return(LOS.sim)

}

#' CoMP and Gamma convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pCoMPoissonGamma <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      x=seq(0,q,0.1)
      res = lapply(x,function(x) sum(dgamma(x-rd.p,mug, sigmag)*dcmp(rd.p,nu=size,lambda = mu)))
      #res = sum(unlist(lapply(res,function(x) max(x,10^(-300))))*0.1)
      res = sum(sapply(res,function(x) max(x,10^(-300)))*0.1)
   }
   if(length(q)>1){
      res = NULL
      for(j in 1:length(q)){
         x=seq(0,q[j],0.1)
         res.i = lapply(x,function(x) sum(dgamma(x-rd.p,mug,sigmag)*dcmp(rd.p,nu=size,lambda = mu)))
         res.i = sum(unlist(lapply(res.i,function(x) max(x,10^(-300))))*0.1)
         res.i = sum(sapply(res.i,function(x) max(x,10^(-300)))*0.1)
         res = c(res,res.i)
      }
   }
   return(res)

}

#' CoMP and Gamma convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qCoMPoissonGamma <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rCoMPoissonGamma(10000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rCoMPoissonGamma(10000,mu[i],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### ZIP / Gamma ####

#' ZIP and Gamma convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rZIPGamma <- function(N,mu,size,mug,sigmag){
   zip = rbinom(N,1,size)
   rnb.sim = rpois(N,mu)
   rnb.sim[which(zip == 1)] = 0
   LOS.sim = rnb.sim + rgamma(N,mug,sigmag)
   return(LOS.sim)
}

#' ZIP and Gamma convolution
#' @param q quantiles
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pZIPGamma <- function(q,mu,size,mug,sigmag){
   l=1.5
   rd.p = 0:200
   if(length(q)==1){
      res = sum(pgamma(q-rd.p,mug, sigmag)*dZIP(rd.p, mu, size))
   }
   if(length(q)>1){
      x=q
      res = sapply(x,function(x) sum(pgamma(x-rd.p,mug,sigmag)*dZIP(rd.p, mu, size)))
   }
   return(res)
}

#' ZIP and Gamma convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qZIPGamma <- function(p,mu,size,mug,sigmag){

   if(length(mu)==1){
      sample = rZIPGamma(100000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(mu)>1){
      q = NULL
      for(i in 1:length(mu)){
         sample = rZIPGamma(100000,mu[i],sigma[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }
   return(q)

}

#### Binom / Gamma ####

#' Binomial and Gamma convolution
#' @param N sample size
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rBinomGamma <- function(N,prob,size,mug,sigmag){
   rb.sim = rbinom(N,size=size,prob=prob)
   LOS.sim = rb.sim + rgamma(N,mug,sigmag)
   return(LOS.sim)
}

#' Binomial and Gamma convolution
#' @param q quantiles
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
pBinomGamma <- function(q,prob,size,mug,sigmag){

   rd.p = 0:max(20,size)
   if(length(q)==1){
      x=seq(0,q,0.1)
      res = lapply(x,function(x) sum(dgamma(x-rd.p,mug,sigmag)*dnbinom(rd.p,size=size,prob = prob)))
      #res = sum(unlist(lapply(res,function(x) max(x,10^(-300))))*0.1)
      res = sum(sapply(res,function(x) max(x,10^(-300)))*0.1)
   }

   if(length(q)>1){
      res = NULL
      for(j in 1:length(q)){
         x=seq(0,q[j],0.1)
         res.i = lapply(x,function(x) sum(dgamma(x-rd.p,mug,sigmag)*dnbinom(rd.p,size=size,prob =prob)))
         #res.i = sum(unlist(lapply(res.i,function(x) max(x,10^(-300))))*0.1)
         res.i = sum(sapply(res.i,function(x) max(x,10^(-300)))*0.1)
         res = c(res,res.i)
      }
   }
   return(res)

}

#' Binomial and Gamma convolution
#' @param p probablities
#' @param prob parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qBinomGamma <- function(p,prob,size,mug,sigmag){

   if(length(size)==1){
      sample = rBinomGamma(100000,size,prob,mug,sigmag)
      q = quantile(sample,probs = p, na.rm=TRUE)
   }
   if(length(size)>1){
      q = NULL
      for(i in 1:length(size)){
         sample = rBinomGamma(100000,size[i],prob[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p, na.rm=TRUE))
      }



   }

   return(q)

}

############################ Multinom ##################################

#' Mulinomial and Gaussian convolution
#' @param x value
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @param log output on log-scale?
#' @keywords convolution
#' @export
dMultinomGauss <- function(x,mu,size,mug,sigmag,log=F){

   if(!is.matrix(mu)){
     mu.u = matrix(mu,nrow=1)
     res = dnorm(matrix(kronecker(x,0:(max(size)),"-"),ncol=max(size+1),byrow=T),mug,sigmag)%*%t(mu.u)}
   if(is.matrix(mu)){res = rowSums(dnorm(matrix(kronecker(x,0:(max(size)),"-"),ncol=max(size+1),byrow=T),mug,sigmag)*mu)}
   res = sapply(res,function(x) max(x,10^(-300)))
   res[is.nan(res)] = 10^(-300)
   if(log==T){res=log(res)}
   return(res)

}

#' Mulinomial and Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rMultinomGauss <- function(N,mu,size,mug,sigmag){

   if(!is.matrix(mu)){
      rnb.sim = sample(0:max(size),N,prob=mu, replace=T)
   }else{
      rnb.sim = apply(mu, 1, function(x) sample(0:max(size),1,prob=x))
   }
   LOS.sim = c(t(rnb.sim)) + rnorm(N,mug,sigmag)

   return(LOS.sim)
}

#' Mulinomial and Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qMultinomGauss <- function(p,mu,size,mug,sigmag){
   if(!is.matrix(mu)){
      sample = rMultinomGauss(10000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p)
   }
   if(is.matrix(mu)){
      q=NULL
      for(i in 1:nrow(mu)){
         sample = rMultinomGauss(10000,mu[i,],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p))
      }

   }
   return(q)

}

#' Mulinomial and Log-Gaussian convolution
#' @param N sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rMultinomLnorm <- function(N,mu,size,mug,sigmag){

   if(!is.matrix(mu)){
      rnb.sim = sample(0:max(size),N,prob=mu, replace=T)
   }else{
      rnb.sim = apply(mu, 1, function(x) sample(0:max(size),1,prob=x))
   }
   LOS.sim = c(t(rnb.sim)) + rlnorm(N,mug,sigmag)

   return(LOS.sim)
}

#' Mulinomial and Log-Gaussian convolution
#' @param x value
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @param log output on log-scale?
#' @keywords convolution
#' @export
dMultinomLnorm <- function(x,mu,size,mug,sigmag,log=F){

   if(!is.matrix(mu)){
     mu.u = matrix(mu,nrow=1)
     res = dlnorm(matrix(kronecker(x,0:(max(size)),"-"),ncol=max(size+1),byrow=T),mug,sigmag)%*%t(mu.u)}
   if(is.matrix(mu)){res = rowSums(dlnorm(matrix(kronecker(x,0:(max(size)),"-"),ncol=max(size+1),byrow=T),mug,sigmag)*mu)}
   res = sapply(res,function(x) max(x,10^(-300)))
   res[is.nan(res)] = 10^(-300)
   if(log==T){res=log(res)}
   return(res)

}

#' Mulinomial and Log-Gaussian convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qMultinomLnorm <- function(p,mu,size,mug,sigmag){
   if(!is.matrix(mu)){
      sample = rMultinomLnorm(10000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p)
   }
   if(is.matrix(mu)){
      q=NULL
      for(i in 1:nrow(mu)){
         sample = rMultinomLnorm(10000,mu[i,],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p))
      }

   }
   return(q)

}


#' Mulinomial and Gamma convolution
#' @param N Sample size
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
rMultinomGamma <- function(N,mu,size,mug,sigmag){

   if(!is.matrix(mu)){
      rnb.sim = sample(0:max(size),N,prob=mu, replace=T)
   }else{
      rnb.sim = apply(mu, 1, function(x) sample(0:max(size),1,prob=x))
   }
   LOS.sim = c(t(rnb.sim)) + rgamma(N,mug,sigmag)

   return(LOS.sim)
}

#' Mulinomial and Gamma convolution
#' @param x value
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @param log output on log-scale?
#' @keywords convolution
#' @export
dMultinomGamma <- function(x,mu,size,mug,sigmag,log=F){

   if(!is.matrix(mu)){
      mu.u = matrix(mu,nrow=1)
      res = dlnorm(matrix(kronecker(x,0:(max(size)),"-"),ncol=max(size+1),byrow=T),mug,sigmag)%*%t(mu.u)}
   if(is.matrix(mu)){res = rowSums(dgamma(matrix(kronecker(x,0:(max(size)),"-"),ncol=max(size+1),byrow=T),mug,sigmag)*mu)}
   res = sapply(res,function(x) max(x,10^(-300)))
   res[is.nan(res)] = 10^(-300)
   if(log==T){res=log(res)}
   return(res)

}

#' Mulinomial and Gamma convolution
#' @param p probabilities
#' @param mu parameter1-1
#' @param size parameter1-2
#' @param mug parameter2-1
#' @param sigmag parameter2-2
#' @keywords convolution
#' @export
qMultinomGamma <- function(p,mu,size,mug,sigmag){
   if(!is.matrix(mu)){
      sample = rMultinomGamma(10000,mu,size,mug,sigmag)
      q = quantile(sample,probs = p)
   }
   if(is.matrix(mu)){
      q=NULL
      for(i in 1:nrow(mu)){
         sample = rMultinomGamma(10000,mu[i,],size[i],mug[i],sigmag[i])
         q = rbind(q,quantile(sample,probs = p))
      }

   }
   return(q)

}



#mod: chriso: byte-code optimization of distribution function
# dNbinomLnorm <- cmpfun(.dNbinomLnorm)
dCoMPoissonGauss <- cmpfun(.dCoMPoissonGauss)
dkCoMPoissonGauss <- cmpfun(.dkCoMPoissonGauss)
dCoMPoissonLnorm <- cmpfun(.dCoMPoissonLnorm)
dkCoMPoissonLnorm <- cmpfun(.dkCoMPoissonLnorm)
# dkZIPGauss <- cmpfun(.dkZIPGauss)
