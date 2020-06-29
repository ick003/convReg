dloglikConvreg <- function(theta.var,
                           fixed,
                           idx.var,
                           idx.fixed,
                           y,
                           X.mu1,
                           X.sigma1,
                           X.mu2,
                           X.sigma2,
                           N.THETA, N.SIGMA, N.THETA.2,
                           fun1, fun2, fun3, fun4,
                           dist1,
                           dist2,
                           weights
){

   theta = matrix(NA,nrow = length(theta.var) + length(fixed$value))
   theta[idx.var]  = theta.var
   if(!(length(fixed$value) == 0)){theta[idx.fixed] = fixed$value}

   mu = fun1(X.mu1 %*% theta[1:N.THETA]) # formula.mu1
   size = fun2(X.sigma1 %*% theta[(N.THETA+1):(N.THETA+N.SIGMA)]) # formula.sigma1

   mu.g = fun3(X.mu2 %*% theta[(N.THETA+N.SIGMA + 1):(N.THETA+N.SIGMA + N.THETA.2)]) # formula.mu2
   sigma = fun4(X.sigma2 %*% theta[(N.THETA+N.SIGMA+N.THETA.2 + 1):length(theta)]) # formula.sigma2

   mu[which(mu>10^10)] = 10^10
   size[which(size>10^10)] = 10^10
   sigma[which(sigma>10^10)] = 10^10

   size[which(size<10^-10)] = 10^-10
   mu[which(mu<10^-10)] = 10^-10
   sigma[which(sigma<10^-10)] = 10^-10

   if(dist1=="Binom"){mu[which(is.nan(mu))] = 1}

   if(paste('d',dist1,dist2,sep="") == "dBinomGauss"){

      k.max = max(size)

      Y = matrix(y,nrow = length(y), ncol=k.max+1, byrow = F)
      K = matrix(0:k.max,nrow = length(y), ncol = k.max+1, byrow=T)

      num = dkBinomGauss(y,size,mu,mu.g,sigma)

      denom = rowSums(num)

      prob.num = (K - matrix(size*mu,ncol=k.max+1, nrow=length(y),byrow=F) )*(num)

      prob.grad = t(X.mu1 * weights) %*%  matrix(rowSums(prob.num)/denom,ncol=1)

      #size.num = (size*digamma(size + K) - K *(size)/(size+mu) - size * digamma(size) + size *log(size/(size+mu)) +size - size^2/(size+mu)) * (num)

      #size.grad = t(X.sigma1) %*%  matrix(rowSums(size.num)/denom,ncol=1)

      mug.num = ((Y-K)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mug.num)/denom, ncol=1)

      sigma.num = (-2*(Y-K)*c(mu.g)/sigma^2+ (Y-K)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(prob.grad,0,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") == "dPoisGauss"){

      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)

      num = dkPoisGauss(y,mu,mu.g,sigma)

      denom = rowSums(num)

      lam.num = (K - mu )*(num)

      lam.grad = t(X.mu1 * weights) %*%  matrix(rowSums(lam.num)/denom,ncol=1)

      #size.num = (size*digamma(size + K) - K *(size)/(size+mu) - size * digamma(size) + size *log(size/(size+mu)) +size - size^2/(size+mu)) * (num)

      #size.grad = t(X.sigma1) %*%  matrix(rowSums(size.num)/denom,ncol=1)

      mug.num = ((Y-K)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mug.num)/denom, ncol=1)

      sigma.num = (-2*(Y-K)*c(mu.g)/sigma^2+ (Y-K)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(lam.grad,0,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") == "dZIPGauss"){

      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)
      pi.h = fun1(X.sigma1 %*% theta[(N.THETA+1):(N.THETA+N.SIGMA)])

      num = dkPoisGauss(y,mu,mu.g,sigma)

      denom = rowSums(num)*(1-size) + size * dnorm(y,mu.g,sigma)

      lam.num = (K - mu )*(num)

      lam.grad = t(X.mu1 * weights) %*%  (matrix(rowSums(lam.num)/denom,ncol=1) * (1-size) + size * dnorm(y,mu.g,sigma) / denom)

      size.num = (1 - rowSums(num))*size * pi.h

      #size.num = (size*digamma(size + K) - K *(size)/(size+mu) - size * digamma(size) + size *log(size/(size+mu)) +size - size^2/(size+mu)) * (num)

      size.grad = t(X.sigma1 * weights) %*%  matrix(size.num/denom,ncol=1)

      mug.num = ((Y-K)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  (matrix(rowSums(mug.num)/denom, ncol=1)* (1 - size) + (y / sigma^2 - mu.g/sigma^2)*size * dnorm(y,mu.g,sigma))

      sigma.num = (-2*(Y-K)*c(mu.g)/sigma^2+ (Y-K)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  (matrix(rowSums(sigma.num)/denom,ncol=1)* (1 - size) + (-2*(y)*c(mu.g)/sigma^2+ (y)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (dnorm(y,mu.g,sigma)) * size)

      return(c(lam.grad,size.grad,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") == "dNbinomGauss"){

      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)

      #          dkN = function(k){dkNbinomGauss(y,mu,size,mu.g,sigma,k)}
      #          dkn = Vectorize(dkN)
      #          #num = dkn(0:99)
      num = sapply(0:99,function(k) dkNbinomGauss(y,mu,size,mu.g,sigma,k))
      denom = rowSums(num, na.rm=T)

      mu.num = (K - K*(mu/(mu+size)) - mu*size / (mu + size))*(num)

      mu.grad = t(X.mu1 * weights) %*%  matrix(rowSums(mu.num)/denom,ncol=1)

      size.num = (size*digamma(size + K) - K *(size)/(size+mu) - size * digamma(size) + size *log(size/(size+mu)) +size - size^2/(size+mu)) * (num)

      size.grad = t(X.sigma1 * weights) %*%  matrix(rowSums(size.num)/denom,ncol=1)

      mug.num = ((Y-K)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mug.num)/denom, ncol=1)

      sigma.num = (-2*(Y-K)*c(mu.g)/sigma^2+ (Y-K)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(mu.grad,size.grad,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") == "dCoMPoissonGauss"){

      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)
      K.fac = matrix(log(factorial(0:99)),nrow = length(y), ncol = 100, byrow=T)

      num =  dkCoMPoissonGauss(y,mu,size,mu.g,sigma,0:99)

      denom = rowSums(num)

      denom[denom==0] = 10^-300

      lambda.num = (K - dZ.lam(mu,size)/Zcomp(mu,size))*(num)

      lambda.grad = t(X.mu1 * weights) %*%  matrix(rowSums(lambda.num)/denom,ncol=1)

      nu.num = (-K.fac*size + size* dZ.nu(mu,size)/Zcomp(mu,size)) * (num)

      nu.grad = t(X.sigma1 * weights) %*%  matrix(rowSums(nu.num)/denom,ncol=1)

      mu.num = ((Y-K)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mu.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mu.num)/denom, ncol=1)

      sigma.num = (-2*(Y-K)*c(mu.g)/sigma^2+ (Y-K)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(lambda.grad,nu.grad,mu.grad,sigma.grad)[idx.var])

   }

   if(paste('d',dist1,dist2,sep="") == "dBinomLnorm"){

      k.max = max(size)

      Y = matrix(y,nrow = length(y), ncol=k.max+1, byrow = F)
      K = matrix(0:k.max,nrow = length(y), ncol = k.max+1, byrow=T)
      YK = Y-K
      idYK = which(YK <= 0, arr.ind=T)
      YK[idYK] = 0.1
      lYK = log(YK)
      lYK[idYK] = 0

      num = dkBinomLnorm(y,size,mu,mu.g,sigma)

      denom = rowSums(num)

      prob.num = (K - matrix(size*mu,ncol=k.max+1, nrow=length(y),byrow=F) )*(num)

      prob.grad = t(X.mu1 * weights) %*%  matrix(rowSums(prob.num)/denom,ncol=1)

      #size.num = (size*digamma(size + K) - K *(size)/(size+mu) - size * digamma(size) + size *log(size/(size+mu)) +size - size^2/(size+mu)) * (num)

      #size.grad = t(X.sigma1) %*%  matrix(rowSums(size.num)/denom,ncol=1)

      mug.num = ((lYK)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mug.num)/denom, ncol=1)

      sigma.num = (-2*(lYK)*c(mu.g)/sigma^2+ (lYK)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(prob.grad,0,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") == "dPoisLnorm"){

      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)
      YK = Y-K
      idYK = which(YK <= 0, arr.ind=T)
      YK[idYK] = 0.1
      lYK = log(YK)
      lYK[idYK] = 0

      num = dkPoisLnorm(y,mu,mu.g,sigma)

      denom = rowSums(num)

      lam.num = (K - mu)*(num)

      lam.grad = t(X.mu1 * weights) %*%  matrix(rowSums(lam.num)/denom,ncol=1)

      mug.num = ((lYK)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mug.num)/denom, ncol=1)

      sigma.num = (-2*(lYK)*c(mu.g)/sigma^2+ (lYK)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(lam.grad,0,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") ==  "dNbinomLnorm"){
      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)
      YK = Y-K
      idYK = which(YK <= 0, arr.ind=T)
      YK[idYK] = 0.1
      lYK = log(YK)
      lYK[idYK] = 0

      num = sapply(0:99,function(k) dkNbinomLnorm(y,mu,size,mu.g,sigma,k))
      denom = rowSums(num)

      mu.num = (K - K*(mu/(mu+size)) - mu*size / (mu + size))*(num)

      mu.grad = t(X.mu1 * weights) %*%  matrix(rowSums(mu.num)/denom,ncol=1)

      size.num = (size*digamma(size + K) - K *(size)/(size+mu) - size * digamma(size) + size *log(size/(size+mu)) +size - size^2/(size+mu)) * (num)

      size.grad = t(X.sigma1 * weights) %*%  matrix(rowSums(size.num)/denom,ncol=1)

      mug.num = ((lYK)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mug.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mug.num)/denom, ncol=1)

      sigma.num = (-2*(lYK)*c(mu.g)/sigma^2+ (lYK)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(mu.grad,size.grad,mug.grad,sigma.grad)[idx.var])
   }
   if(paste('d',dist1,dist2,sep="") ==  "dCoMPoissonLnorm"){
      Y = matrix(y,nrow = length(y), ncol=100, byrow = F)
      K = matrix(0:99,nrow = length(y), ncol = 100, byrow=T)
      K.fac = matrix(log(factorial(0:99)),nrow = length(y), ncol = 100, byrow=T)

      YK = Y-K
      idYK = which(YK <= 0, arr.ind=T)
      YK[idYK] = 0.1
      lYK = log(YK)
      lYK[idYK] = 0

      num = sapply(0:99,function(k) dkCoMPoissonLnorm(y,mu,size,mu.g,sigma,k))

      denom = rowSums(num)

      lambda.num = (K - dZ.lam(mu,size)/Zcomp(mu,size))*(num)

      lambda.grad = t(X.mu1 * weights) %*%  matrix(rowSums(lambda.num)/denom,ncol=1)

      nu.num = (-K.fac*size + size*dZ.nu(mu,size)/Zcomp(mu,size)) * (num)

      nu.grad = t(X.sigma1 * weights) %*%  matrix(rowSums(nu.num)/denom,ncol=1)

      mu.num = ((lYK)/sigma^2 -c(mu.g)/sigma^2)*(num)

      mu.grad = t(X.mu2 * weights) %*%  matrix(rowSums(mu.num)/denom, ncol=1)

      sigma.num = (-2*(lYK)*c(mu.g)/sigma^2+ (lYK)^2/sigma^2 + c(mu.g)^2/sigma^2-1) * (num)

      sigma.grad = t(X.sigma2 * weights) %*%  matrix(rowSums(sigma.num)/denom,ncol=1)

      return(c(lambda.grad,nu.grad,mu.grad,sigma.grad)[idx.var])
   }

 ##  return(as.numeric(dl))

}
