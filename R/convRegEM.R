convRegEM = function(formula.resp,formula.mu1,formula.sigma1, formula.mu2,formula.sigma2,
                     dist1 = "Nbinom", dist2 = "Gauss", df, theta.0,
                     X.mu1,X.sigma1,X.mu2,X.sigma2,fun1,fun2,fun3,fun4, weights,
                     quiet = FALSE, debug = FALSE, scale = FALSE, scaleInit = 1){


   if(debug){browser()}

   if(dist1 %in% c("Pois", "Nbinom", "ZIP", "HP", "CoMPoisson")){maxK = 49}
   if(dist1 %in% c("Binom")){maxK = theta.0$init[theta.0$indices$idx1b]}

   data = as.data.frame(df[rep(seq_len(nrow(df)), each=length(0:maxK)),])
   names(data) = names(df)
   data$k = 0:maxK

   wPar = weights

   theta = matrix(c(theta.0$init),ncol=1)

   theta.hist = matrix(theta, nrow = 1)
   i = 1
   diff=  1

   idx.m1 = theta.0$indices$idx1
   idx.s1 = theta.0$indices$idx1b
   idx.m2 = theta.0$indices$idx2
   idx.s2 = theta.0$indices$idx2b

   R.Cmd <- sprintf("data$resp=data$%s", all.vars(formula.resp))
   eval(parse(text=R.Cmd))


   if(dist2 %in% "Gauss"){wdist2 = "norm"}
   if(dist2 %in% "Lnorm"){wdist2 = "lnorm"}
   if(dist2 %in% "Gamma"){wdist2 = "gamma"}

   if(dist1 %in% "Pois"){wdist1 = "poisson"}
   if(dist1 %in% "Nbinom"){wdist1 = "nbinom"}
   if(dist1 %in% "ZIP"){wdist1 = "ZIP"}
   if(dist1 %in% "HP"){wdist1 = "HP"}
   if(dist1 %in% "CoMPoisson"){wdist1 = "comp"}
   if(dist1 %in% "Binom"){wdist1 = "binom"}


if(!scale){

   while(diff > 10e-8 & i < 500){
      i = i + 1
      if(!quiet){if((i / 5) == round(i/5)){print(i)}}

      mu = rep(fun1(X.mu1 %*% theta[idx.m1]), each = maxK+1) # formula.mu1
      if(!(dist1 %in% c("Pois"))){size = rep(fun2(X.sigma1 %*% theta[idx.s1]), each = maxK+1)} # formula.sigma1

      mu.g = rep(fun3(X.mu2 %*% theta[idx.m2]), each = maxK+1) # formula.mu2
      sigma = rep(fun4(X.sigma2 %*% theta[idx.s2]), each = maxK+1) # formula.sigma2

      mu[which(mu>10^10)] = 10^10
      if(!(dist1 %in% "Pois")){size[which(size>10^10)] = 10^10}
      sigma[which(sigma>10^10)] = 10^10

      sigma[which(sigma<10^-10)] = 10^-10
      if(!(dist1 %in% "Pois")){size[which(size<10^-10)] = 10^-10}
      mu[which(mu<10^-10)] = 10^-10

      if(dist1=="Pois"){
         R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + data$k, sigma)*dpois(data$k, mu) * wPar", wdist2)
      }else if(dist1=="Nbinom"){
         R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + data$k, sigma)*dnbinom(data$k, mu= mu, size = size)* wPar", wdist2)
      }else if(dist1=="CoMPoisson"){
         R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + data$k, sigma)*dcomp(data$k, lam= mu, nu = size, sumTo=50)* wPar", wdist2, wdist1)
      }else if(dist1=="ZIP"){
         R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + data$k, sigma)*dZIP(data$k, lam= mu, pi = size)* wPar", wdist2, wdist1)
      }else if(dist1=="Binom"){
         R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + data$k, sigma)*dbinom(data$k, prob= mu, size = size)* wPar", wdist2, wdist1)
      }else if(dist1=="HP"){
         R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + data$k, sigma)*dHP(data$k, lam= mu, pi = 1-size)* wPar", wdist2, wdist1)
      }
      eval(parse(text=R.Cmd))

      dw = data.frame(weights = weights, ID = rep(seq_len(nrow(df)), each=maxK+1))
      #ans = setDT(dw)[,.(A = sum(weights)), by = 'ID']

      #dw = data.frame(w = weights, ID = rep(seq_len(nrow(df)), each=50))

      agg.dw = aggregate(weights~ID, data = dw, sum)

      #weights = weights / rep(ans$A, each = maxK+1)

      weights = weights / rep(agg.dw$weights, each = maxK+1)

      data.t = data[which(weights > 1/(maxK+1)),]

      weights = weights[which(weights>1/(maxK+1))]

      if(dist2 %in% "Gauss"){
         formu = eval(parse(text=paste("resp~",formula.mu2[2],"",sep="")))
         res.2 = glm(formu, data = data.t, weights = weights, offset=data.t$k, family = gaussian)
         c21 = coef(res.2)
         c22 = log(sqrt(deviance(res.2)/df.residual(res.2)))
      }
      if(dist2 %in% "Lnorm"){
         formu = eval(parse(text=paste("log(resp)~",formula.mu2[2],"",sep="")))
         res.2 = glm(formu, data = data.t, weights = weights, offset=data.t$k, family = gaussian)
         c21 = coef(res.2)
         c22 = log(sqrt(deviance(res.2)/df.residual(res.2)))
      }
      if(dist2 %in% "Gamma"){
         formu = eval(parse(text=paste("resp~",formula.mu2[2],"",sep="")))
         res.2 = glm(formu, data = data.t, weights = weights, offset=data.t$k, family = Gamma)
         c21 = coef(res.2)
         c22 = log(sqrt(deviance(res.2)/df.residual(res.2)))
      }

      if(dist1 %in% "Pois"){
         formu = eval(parse(text=paste("k~",formula.mu1[2],sep="")))
         res.1 = glm(formu, data = data.t, weights = weights, family = poisson)
         c11 = coef(res.1)
      }

      if(dist1 %in% "Nbinom"){
         formu = eval(parse(text=paste("k~",formula.mu1[2],sep="")))
         res.1 = tryCatch({glm.nb(formu, data = data.t, weights = weights, init.theta = 1)},
                          error = function(err){stop("Error: NegBin regression error")})
         c11 = res.1$coefficients
         c12 = log(res.1$theta)
      }
      if(dist1 %in% "ZIP"){
         if(diff(range(data.t$k))>0){
         formu = eval(parse(text=paste("k~",formula.mu1[2], "|",formula.sigma1[2] ,sep="")))
         res.1 = zeroinfl(formu, data = data.t, weights = weights)
         c11 = res.1$coefficients$count
         c12 = res.1$coefficients$zero
         }else{
            c11=c11
            c12 = c12
         }

      }
      if(dist1 %in% "HP"){
         if(diff(range(data.t$k))>0){
         formu = eval(parse(text=paste("k~",formula.mu1[2], "|",formula.sigma1[2] ,sep="")))
         res.1 = hurdle(formu, data = data.t, weights = weights, dist = "poisson")
         c11 = res.1$coefficients$count
         c12 = res.1$coefficients$zero
         }else{
            c11=c11
            c12 = c12
         }
      }
      if(dist1 %in% "CoMPoisson"){

         data.tt = data.t
         data.tt$os = 0
         formu1 = eval(parse(text=paste("k~",formula.mu1[2],sep="")))
         nuFormula =eval(parse(text=paste("k~",formula.sigma1[2],"+offset(os)",sep="")))

         res.1 = glmComp(lamFormula = formu1, nuFormula = nuFormula, data = data.tt, weights = weights, sumTo = 50)
         c11 = coef(res.1)$beta
         c12 = coef(res.1)$zeta
      }
      if(dist1 %in% "Binom"){
         formu1 = eval(parse(text=paste("k/",theta[idx.s1],"~",formula.mu1[2],sep="")))
         #res.1 = multinom(formu1,data=data.t, weights = weights)
         res.1 = glm(formu1,data=data.t, weights = weights, family = "quasibinomial")
         c11 = coef(res.1)
         c12 = theta[idx.s1]
      }

      if(dist1 %in% "Pois"){
         theta[c(idx.m1, idx.m2, idx.s2)] = c(c11, c21, c22)
      }else{
         theta = c(c11, c12,c21, c22)
      }

      diff = max((theta.hist[nrow(theta.hist),] - theta)^2 / abs(theta), na.rm=T)
      theta.hist = rbind(theta.hist,c(theta))

   }

   mu = fun1(X.mu1 %*% theta[idx.m1]) # formula.mu1

   size = NULL
if(!(dist1 %in% "Pois")){
   size = fun2(X.sigma1 %*% theta[idx.s1])# formula.sigma1
   size[which(size>10^10)] = 10^10
   size[which(size<10^-10)] = 10^-10
}

   mu.g = fun3(X.mu2 %*% theta[idx.m2]) # formula.mu2
   sigma = fun4(X.sigma2 %*% theta[idx.s2]) # formula.sigma2

   mu[which(mu>10^10)] = 10^10

   sigma[which(sigma>10^10)] = 10^10

   sigma[which(sigma<10^-10)] = 10^-10

   mu[which(mu<10^-10)] = 10^-10

   weights.final = weights

   theta.est = theta

   if(dist1 %in% c("Pois" ,"Binom")){theta.est = theta[c(idx.m1, idx.m2, idx.s2)]}

   RET = list(par = theta.est, par.hist = theta.hist, m1 = mu, s1 = size,m2 = mu.g, s2 = sigma, weights = weights.final,
              opt1 = res.1, opt2 = res.2, wPar = wPar)

}

if(scale){

      scaleEst = scaleInit
      scaleHist = scaleEst

      if(debug){browser()}

      while(diff > 10e-8 & i < 500){

         i = i + 1
         if(!quiet){if((i / 20) == round(i/20)){print(i)}}

         mu = rep(fun1(X.mu1 %*% theta[idx.m1]), each = maxK+1) # formula.mu1
         if(!(dist1 %in% c("Pois"))){size = rep(fun2(X.sigma1 %*% theta[idx.s1]), each = maxK+1)} # formula.sigma1

         mu.g = rep(fun3(X.mu2 %*% theta[idx.m2]), each = maxK+1) # formula.mu2
         sigma = rep(fun4(X.sigma2 %*% theta[idx.s2]), each = maxK+1) # formula.sigma2

         mu[which(mu>10^10)] = 10^10
         if(!(dist1 %in% "Pois")){size[which(size>10^10)] = 10^10}
         sigma[which(sigma>10^10)] = 10^10

         sigma[which(sigma<10^-10)] = 10^-10
         if(!(dist1 %in% "Pois")){size[which(size<10^-10)] = 10^-10}
         mu[which(mu<10^-10)] = 10^-10

         if(dist1=="Pois"){
            R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + scaleEst*data$k, sigma)*dpois(data$k, mu) * wPar", wdist2)
         }else if(dist1=="Nbinom"){
            R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + scaleEst*data$k, sigma)*dnbinom(data$k, mu= mu, size = size)* wPar", wdist2)
         }else if(dist1=="CoMPoisson"){
            R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + scaleEst*data$k, sigma)*dcomp(data$k, lam= mu, nu = size, sumTo=50)* wPar", wdist2, wdist1)
         }else if(dist1=="ZIP"){
            R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + scaleEst*data$k, sigma)*dZIP(data$k, lam= mu, pi = size)* wPar", wdist2, wdist1)
         }else if(dist1=="Binom"){
            R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + scaleEst*data$k, sigma)*dbinom(data$k, prob= mu, size = size)* wPar", wdist2, wdist1)
         }else if(dist1=="HP"){
            R.Cmd <- sprintf("weights = d%s(data$resp,mu.g + scaleEst*data$k, sigma)*dHP(data$k, lam= mu, pi = 1-size)* wPar", wdist2, wdist1)
         }
         eval(parse(text=R.Cmd))
         #browser()

         dw = data.frame(weights = weights, ID = rep(seq_len(nrow(df)), each=maxK+1))
         #ans = setDT(dw)[,.(A = sum(weights)), by = 'ID']

         #dw = data.frame(w = weights, ID = rep(seq_len(nrow(df)), each=50))

         agg.dw = aggregate(weights~ID, data = dw, sum)

         #weights = weights / rep(ans$A, each = maxK+1)
         weights = weights / rep(agg.dw$weights, each = maxK+1)

         data.t = data[which(weights > 1/(maxK+1)),]

         weights = weights[which(weights>1/(maxK+1))]

         if(dist2 %in% "Gauss"){
            formu = eval(parse(text=paste("resp~k+",formula.mu2[2],"",sep="")))
            res.2 = glm(formu, data = data.t, weights = weights, family = gaussian)
            idxK = match("k",names(coef(res.2)))
            c21 = coef(res.2)[-idxK]
            scaleEst = coef(res.2)[idxK]
            c22 = log(sqrt(deviance(res.2)/df.residual(res.2)))
         }
         if(dist2 %in% "Lnorm"){
            formu = eval(parse(text=paste("log(resp)~k+",formula.mu2[2],"",sep="")))
            res.2 = glm(formu, data = data.t, weights = weights, family = gaussian)
            idxK = match("k",names(coef(res.2)))
            c21 = coef(res.2)[-idxK]
            scaleEst = coef(res.2)[idxK]
            c22 = log(sqrt(deviance(res.2)/df.residual(res.2)))
         }
         if(dist2 %in% "Gamma"){
            formu = eval(parse(text=paste("resp~k+",formula.mu2[2],"",sep="")))
            res.2 = glm(formu, data = data.t, weights = weights, family = Gamma)
            idxK = match("k",names(coef(res.2)))
            c21 = coef(res.2)[-idxK]
            scaleEst = coef(res.2)[idxK]
            c22 = log(sqrt(deviance(res.2)/df.residual(res.2)))
         }

         if(dist1 %in% "Pois"){
            formu = eval(parse(text=paste("k~",formula.mu1[2],sep="")))
            res.1 = glm(formu, data = data.t, weights = weights, family = poisson)
            c11 = coef(res.1)
         }

         if(dist1 %in% "Nbinom"){
            formu = eval(parse(text=paste("k~",formula.mu1[2],sep="")))
            res.1 = tryCatch({glm.nb(formu, data = data.t, weights = weights, init.theta = 1)},
                             error = function(err){stop("Error: NegBin regression error")})
            c11 = res.1$coefficients
            c12 = log(res.1$theta)
         }
         if(dist1 %in% "ZIP"){
            if(diff(range(data.t$k))>0){
               formu = eval(parse(text=paste("k~",formula.mu1[2], "|",formula.sigma1[2] ,sep="")))
               res.1 = zeroinfl(formu, data = data.t, weights = weights)
               c11 = res.1$coefficients$count
               c12 = res.1$coefficients$zero
            }else{
               c11=c11
               c12 = c12
            }

         }
         if(dist1 %in% "HP"){
            if(diff(range(data.t$k))>0){
               formu = eval(parse(text=paste("k~",formula.mu1[2], "|",formula.sigma1[2] ,sep="")))
               res.1 = hurdle(formu, data = data.t, weights = weights, dist = "poisson")
               c11 = res.1$coefficients$count
               c12 = res.1$coefficients$zero
            }else{
               c11=c11
               c12 = c12
            }
         }
         if(dist1 %in% "CoMPoisson"){

            data.tt = data.t
            data.tt$os = 0
            formu1 = eval(parse(text=paste("k~",formula.mu1[2],sep="")))
            nuFormula =eval(parse(text=paste("k~",formula.sigma1[2],"+offset(os)",sep="")))

            res.1 = glmComp(lamFormula = formu1, nuFormula = nuFormula, data = data.tt, weights = weights, sumTo = 50)
            c11 = coef(res.1)$beta
            c12 = coef(res.1)$zeta
         }
         if(dist1 %in% "Binom"){
            formu1 = eval(parse(text=paste("k/",theta[idx.s1],"~",formula.mu1[2],sep="")))
            #res.1 = multinom(formu1,data=data.t, weights = weights)
            res.1 = glm(formu1,data=data.t, weights = weights, family = "quasibinomial")
            c11 = coef(res.1)
            c12 = theta[idx.s1]
         }

         if(dist1 %in% "Pois"){
            theta[c(idx.m1, idx.m2, idx.s2)] = c(c11, c21, c22)
         }else{
            theta = c(c11, c12,c21, c22)
         }

         if(is.na(scaleEst)){scaleEst = scaleHist[i-1]}

         diff = max((theta.hist[nrow(theta.hist),] - theta)^2 / abs(theta), na.rm=T)
         diff = max(diff, (scaleHist[length(scaleHist)] - scaleEst)^2 / abs(scaleEst))
         theta.hist = rbind(theta.hist,c(theta))
         scaleHist = c(scaleHist,scaleEst)
      }

      mu = fun1(X.mu1 %*% theta[idx.m1]) # formula.mu1

      size = NULL
      if(!(dist1 %in% "Pois")){
         size = fun2(X.sigma1 %*% theta[idx.s1])# formula.sigma1
         size[which(size>10^10)] = 10^10
         size[which(size<10^-10)] = 10^-10
      }

      mu.g = fun3(X.mu2 %*% theta[idx.m2]) # formula.mu2
      sigma = fun4(X.sigma2 %*% theta[idx.s2]) # formula.sigma2

      mu[which(mu>10^10)] = 10^10

      sigma[which(sigma>10^10)] = 10^10

      sigma[which(sigma<10^-10)] = 10^-10

      mu[which(mu<10^-10)] = 10^-10

      weights.final = weights

      theta.est = theta

      if(dist1 %in% c("Pois" ,"Binom")){theta.est = theta[c(idx.m1, idx.m2, idx.s2)]}



      RET = list(par = theta.est, par.hist = theta.hist, m1 = mu, s1 = size,m2 = mu.g, s2 = sigma, weights = weights.final,
                 opt1 = res.1, opt2 = res.2, wPar = wPar, scale = scaleEst, scale.hist = scaleHist)

   }

   return(RET)
}
