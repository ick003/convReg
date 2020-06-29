loglikConvreg <- function(theta.var,
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

   sigma[which(sigma<10^-10)] = 10^-10
   size[which(size<10^-10)] = 10^-10
   mu[which(mu<10^-10)] = 10^-10

   if(dist1=="Binom"){mu[which(is.nan(mu))] = 1}

      # mod: chriso, customizing name
   switch(paste('d',dist1,dist2,sep=""),

          "dBinomGauss"={ return(
             sum(
                weights*dBinomGauss(y,
                            mu,
                            size,
                            mu.g,
                            sigma,
                            log=TRUE),
                na.rm=T)
          )
          },
          "dPoisGauss"={ return(
             sum(
                weights*dPoisGauss(y,
                           mu,
                           mu.g,
                           sigma,
                           log=TRUE),
                na.rm=T)
          )
          },
          "dZIPGauss"={ return(
             sum(
                weights*dZIPGauss(y,
                           mu,
                           size,
                           mu.g,
                           sigma,
                           log=TRUE),
                na.rm=T)
          )
          },
          "dHPGauss"={ return(
             sum(
                weights*dHPGauss(y,
                          mu,
                          size,
                          mu.g,
                          sigma,
                          log=TRUE),
                na.rm=T)
          )
          },
          "dNbinomGauss"={
             return(
                sum(
                   weights*dNbinomGauss(y,
                                mu,
                                size,
                                mu.g,
                                sigma,
                                log=TRUE),
                   na.rm=T)
             )
          },
          "dCoMPoissonGauss"={
             return(
                sum(
                   weights*dCoMPoissonGauss(y,
                                    mu,
                                    size,
                                    mu.g,
                                    sigma,
                                    log=TRUE),
                   na.rm=T)
             )
          },



          "dBinomLnorm"={ return(
             sum(
                weights*dBinomLnorm(y,
                            mu,
                            size,
                            mu.g,
                            sigma,
                            log=TRUE),
                na.rm=T)
          )
          },
          "dPoisLnorm"={ return(
             sum(
                weights*dPoisLnorm(y,
                           mu,
                           mu.g,
                           sigma,
                           log=TRUE),
                na.rm=T)
          )
          },
          "dNbinomLnorm"={
             return(
                sum(
                   weights*dNbinomLnorm(y,
                                mu,
                                size,
                                mu.g,
                                sigma,
                                log=TRUE),
                   na.rm=T)
             )
          },
          "dCoMPoissonLnorm"={
             return(
                sum(
                   weights*dCoMPoissonLnorm(y,
                                    mu,
                                    size,
                                    mu.g,
                                    sigma,
                                    log=TRUE),
                   na.rm=T)
             )
          }



   )
}
