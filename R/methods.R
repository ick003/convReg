fitConvreg <- function(paramList, method){

   theta.est = rep(0,ncol(paramList$X.mu1)+ncol(paramList$X.sigma1)+ncol(paramList$X.mu2)+ncol(paramList$X.sigma2))

   theta.name = c(paste("mu 1: ",colnames(paramList$X.mu1),sep=""),
                  paste("sigma 1: ",colnames(paramList$X.sigma1),sep=""),
                  paste("mu 2: ",colnames(paramList$X.mu2),sep=""),
                  paste("sigma 2: ",colnames(paramList$X.sigma2),sep=""))

   theta.est[paramList$init.theta$idx1] = paramList$init.theta$theta.mu1
   theta.est[paramList$init.theta$idx1b] = paramList$init.theta$theta.sigma1
   theta.est[paramList$init.theta$idx2] = paramList$init.theta$theta.mu2
   theta.est[paramList$init.theta$idx2b] = paramList$init.theta$theta.sigma2

   theta.est[paramList$idx.fixed] = paramList$fixed$value

   names(theta.est) = theta.name

   idx.na = which(is.na(theta.est))

   if(length(idx.na) > 0){
      if(!paramList$quiet){message(sprintf('NA coefficients estimated, dropping variables: %s', theta.name[idx.na]))}

      paramList$idx.fixed = c(paramList$idx.fixed, idx.na)
      paramList$idx.var   =  which(!(theta.name %in% c(paramList$fixed$name, theta.name[idx.na])))

      paramList$fixed$name = c(paramList$fixed$name,theta.name[idx.na])
      paramList$fixed$value = c(paramList$fixed$value, rep(0, length(idx.na)))
   }

   N.THETA   <- ncol(paramList$X.mu1)
   N.THETA.2 <- ncol(paramList$X.mu2)
   N.SIGMA   <- ncol(paramList$X.sigma1)

   if(method == "mle"){

      if(paramList$dist1 != "Multinom"){


         # optimization of parameter

         status = 1
         num = 0

         while( status != 0  & num < 50){
            num = num + 1
            if(!paramList$quiet){message(sprintf(' ** optimization attempt num %s : [please wait]', num ))}

            theta.0 = theta.est[paramList$idx.var]

            #
            ll.cv = function(theta){ loglikConvreg(theta,paramList$fixed, paramList$idx.var, paramList$idx.fixed,
                                                      paramList$y,paramList$X.mu1,paramList$X.sigma1, paramList$X.mu2,paramList$X.sigma2,
                                                      N.THETA, N.SIGMA, N.THETA.2,
                                                      paramList$fun1, paramList$fun2, paramList$fun3, paramList$fun4,
                                                      paramList$dist1,paramList$dist2,paramList$weights)}

            dll.cv = function(theta){ dloglikConvreg(theta,paramList$fixed, paramList$idx.var, paramList$idx.fixed,
                                                        paramList$y,paramList$X.mu1,paramList$X.sigma1, paramList$X.mu2,paramList$X.sigma2,
                                                        N.THETA, N.SIGMA, N.THETA.2,
                                                        paramList$fun1, paramList$fun2, paramList$fun3, paramList$fun4,
                                                        paramList$dist1,paramList$dist2,paramList$weights)}

            if(paramList$dist1 %in% c("ZIP","HP")){
               opt.reg = optim(theta.0, fn = ll.cv,gr = NULL, method = "BFGS", control = list(fnscale = -1), hessian=T)
            }else{
               opt.reg = optim(theta.0, fn = ll.cv, gr = dll.cv, method = "BFGS", control = list(fnscale = -1), hessian=T)
            }

            status = opt.reg$convergence
            theta.est[paramList$idx.var] = opt.reg$par
         }

      }
      scaleP = NULL
   }
   if(method == "em"){


      theta0 = list(init = c(paramList$init.theta[[1]],paramList$init.theta[[2]],paramList$init.theta[[3]],paramList$init.theta[[4]]),
                    indices = list(idx1 = paramList$init.theta$idx1,
                                   idx1b = paramList$init.theta$idx1b,
                                   idx2 = paramList$init.theta$idx2,
                                   idx2b = paramList$init.theta$idx2b))



      res.em = convRegEM(paramList$formula.resp,paramList$formula.mu1, paramList$formula.sigma1, paramList$formula.mu2, paramList$formula.sigma2,
                         dist1 = paramList$dist1, dist2 = paramList$dist2,df = paramList$data, theta.0 = theta0,
                         X.mu1 = paramList$X.mu1,X.sigma1 = paramList$X.sigma1,X.mu2 = paramList$X.mu2,X.sigma2 = paramList$X.sigma2,
                         fun1 = paramList$fun1,fun2 = paramList$fun2,fun3 = paramList$fun3,fun4 = paramList$fun4,
                         weights = paramList$weights, quiet = paramList$quiet, debug = paramList$debug,
                         scale = paramList$scale, scaleInit = paramList$scaleInit)

      scaleP = res.em$scale
      status = max(1-res.em$opt1$converged,1-res.em$opt2$converged)
      theta.est[paramList$idx.var] = res.em$par


      ll = function(x0){
         loglikConvreg(x0, paramList$fixed, paramList$idx.var, paramList$idx.fixed,
                        paramList$y,paramList$X.mu1,paramList$X.sigma1, paramList$X.mu2,paramList$X.sigma2,
                        N.THETA, N.SIGMA, N.THETA.2,
                        paramList$fun1, paramList$fun2, paramList$fun3, paramList$fun4,
                        paramList$dist1,paramList$dist2, paramList$weights)
      }

      Hessian = numDeriv::hessian(func = ll, x = theta.est[paramList$idx.var])

      opt.reg = list(value = ll(theta.est[paramList$idx.var]), counts = rep(nrow(res.em$par.hist),2),
                     message = "", convergence = status, hessian = Hessian)

   }

   profile.lik = list(
      param = (theta.est[paramList$idx.var]),
      other = list(paramList$fixed, paramList$idx.var, paramList$idx.fixed,
                   paramList$y,paramList$X.mu1,paramList$X.sigma1, paramList$X.mu2,paramList$X.sigma2,
                   N.THETA, N.SIGMA, N.THETA.2,
                   paramList$fun1, paramList$fun2, paramList$fun3, paramList$fun4,
                   paramList$dist1,paramList$dist2, paramList$weights)
   )

   return(list(opt = opt.reg, prof = profile.lik, theta = theta.est, scalePar = scaleP))
}

#' @export
vcov.convreg <- function(object, ...) {
   vcov <- -solve(object$hessian)
   return(vcov)
}

#' @export
coef.convreg <- function(object, ...) {
   return(list(theta_mu1 = object$estimation$Estimate[object$idx[[1]]],
               theta_sigma1 = object$estimation$Estimate[object$idx[[2]]],
               theta_mu2 = object$estimation$Estimate[object$idx[[3]]],
               theta_sigma2 = object$estimation$Estimate[object$idx[[4]]]))
}

#' @export
predict.convreg <- function(object, newdata = NULL, ...){

  obj = object

  if(is.null(newdata)){newdata = obj$data}

  idx.p1 = index.convreg(obj)[[1]]
  idx.p2 = index.convreg(obj)[[2]]
  idx.p3 = index.convreg(obj)[[3]]
  idx.p4 = index.convreg(obj)[[4]]

  dist1= unlist(strsplit(obj$distname,"/"))[1]
  dist2= unlist(strsplit(obj$distname,"/"))[2]

  formula.mu1 = as.formula(obj$formulas$mu1)
  formula.sigma1 = as.formula(obj$formulas$sigma1)
  formula.mu2 = as.formula(obj$formulas$mu2)
  formula.sigma2 = as.formula(obj$formulas$sigma2)

  X.mu1.model.frame = model.frame(formula.mu1,data = newdata,na.action=NULL, drop.unused.levels = TRUE)
  X.mu1.model.terms = attr(X.mu1.model.frame,"terms")
  X.mu1 = model.matrix(X.mu1.model.terms, X.mu1.model.frame)

  X.sigma1.model.frame = model.frame(formula.sigma1,data = newdata,na.action=NULL, drop.unused.levels = TRUE)
  X.sigma1.model.terms = attr(X.sigma1.model.frame,"terms")
  X.sigma1 = model.matrix(X.sigma1.model.terms, X.sigma1.model.frame)

  X.mu2.model.frame = model.frame(formula.mu2,data = newdata,na.action=NULL, drop.unused.levels = TRUE)
  X.mu2.model.terms = attr(X.mu2.model.frame,"terms")
  X.mu2 = model.matrix(X.mu2.model.terms, X.mu2.model.frame)

  X.sigma2.model.frame = model.frame(formula.sigma2,data = newdata,na.action=NULL, drop.unused.levels = TRUE)
  X.sigma2.model.terms = attr(X.sigma2.model.frame,"terms")
  X.sigma2 = model.matrix(X.sigma2.model.terms, X.sigma2.model.frame)


  if(dist1 != "Multinom"){mu1 = obj$transforms[[1]](X.mu1 %*% obj$estimation$Estimate[idx.p1])}
  if(dist1 == "Multinom"){mu1 = exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.p1], nrow = obj$estimation$Estimate[idx.p2])))) / rowSums( exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.p1], nrow = obj$estimation$Estimate[idx.p2])))))}

  sigma1 = obj$transforms[[2]](X.sigma1 %*% obj$estimation$Estimate[idx.p2])
  mu2 = obj$transforms[[3]](X.mu2 %*% obj$estimation$Estimate[idx.p3])
  sigma2 = obj$transforms[[4]](X.sigma2 %*% obj$estimation$Estimate[idx.p4])

  if(dist1 == "Nbinom" | dist1 == "Pois"){exp1 = mu1}
  if(dist1 == "Binom"){exp1 = mu1*sigma1}
  if(dist1 == "CoMPoisson"){exp1 = Zcomp.mu(mu1,sigma1)}
  if(dist1 == "ZIP"){exp1 = mu1 * (1 - sigma1)}
  if(dist1 == "HP"){exp1 = mu1 / (1 - exp(-mu1)) * (1-sigma1)}

  if(dist2 == "Gauss"){exp2 = mu2}
  if(dist2 == "Lnorm"){exp2 = exp(mu2 + sigma2^2/2)}

  y.fitted = exp1 + exp2

  RET = data.frame(y = obj$xdata$y ,pred = y.fitted, k = exp1, mu1 = mu1, s1 = sigma1, e = exp2, mu2 = mu2, s2= sigma2)

  return(RET)
}

#' @export
summary.convreg <- function(object,... ) {
  if(!inherits(object, "convreg"))
    stop("'summary.convreg' called on a non-'convreg' object")
  ## Here we should actually coerce the object to a 'convreg' object, dropping all the subclasses...
  ## Instead, we force the program to use convreg-related methods
  result <- object$maxim
  nParam <- length(object$estimation$Estimate)
  results <- object$estimation

  summary <- list(maximType=object$method,
                  distr=object$distname,
                  returnCode=object$code,
                  returnMessage=object$message,
                  loglik=object$loglik,
                  estimate=results,
                  NActivePar=nParam,
                  constraints=object$constraints,
                  indices = object$idx)
  class(summary) <- "summary.convreg"
  summary
}

#' @export
print.convreg <- function(object, ...) {
  obj = object
  cat("\nCall:\n")
  print(obj$call)
  cat("\nTheta_mu1:\n")
  print(obj$estimation$Estimate[obj$idx[[1]]])
  cat("\nTheta_mu2:\n")
  print(obj$estimation$Estimate[obj$idx[[3]]])
  cat("\nDegrees of Freedom:", obj$df.null,"Total (i.e. Null);", obj$df.residual , "Residual\n")
  cat("AIC:", obj$aic,"\n")
  cat("Log-Likelihood:", obj$loglik, "\n")
  invisible(obj)
}

#' @export
print.summary.convreg <- function(object, ... ) {
  cat("--------------------------------------------\n")
  cat("Convolution Regression Results","\n")
  cat(toupper(object$maximType), ", ", object$distr,"\n", sep="")
  cat("Return code ",object$returnCode, ": ", object$returnMessage, "\n", sep="")
  if(!is.null(object$estimate)) {
    cat("Log-Likelihood:", object$loglik, "\n")
    cat(object$NActivePar, " free parameters\n")
    cat("\n")
    cat(paste("Estimates: ",strsplit(object$distr,"/")[[1]][1] ,"\n"))
    printCoefmat(object$estimate[c(object$indices[[1]],object$indices[[2]]),])
    cat("\n")
    cat(paste("Estimates: ",strsplit(object$distr,"/")[[1]][2] ,"\n"))
    printCoefmat(object$estimate[c(object$indices[[3]],object$indices[[4]]),])

  }
  if(!is.null(object$constraints)) {
    cat("\nWarning: constrained likelihood estimation.",
        "Inference is probably wrong\n")
    cat("Constrained optimization based on", object$constraints$type,
        "\n")
    if(!is.null(object$constraints$code))
      cat("Return code:", object$constraints$code, "\n")
    # note: this is missing for 'constrOptim'
    if(!is.null(object$constraints$message))
      cat(object$constraints$message, "\n")
    # note: this is missing for 'constrOptim'
    cat(object$constraints$outer.iterations,
        " outer iterations, barrier value",
        object$constraints$barrier.value, "\n")
  }
  cat("--------------------------------------------\n")
}



#' Calculate quantile residuals
#' @param obj convreg object
#' @param ... additional parameters
#' @keywords convolution
#' @export
quantRes <- function(obj, ...) {

   p = predict(obj)

   dist1= unlist(strsplit(obj$distname,"/"))[1]
   dist2= unlist(strsplit(obj$distname,"/"))[2]

   if(dist1 == "Pois"){
     E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$mu2,p$s2),1,function(x) p%s%s(x[1],x[2], x[3],x[4]))", dist1, dist2)))
   }else{
     E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x) p%s%s(x[1],x[2], x[3],x[4], x[5]))", dist1, dist2)))
   }


   if(obj$scale){
      if(dist2 %in% "Gauss"){wdist2 = "norm"}
      if(dist2 %in% "Lnorm"){wdist2 = "lnorm"}
      if(dist2 %in% "Gamma"){wdist2 = "gamma"}

      if(dist1 %in% "Pois"){wdist1 = "poisson"
      E =  eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x){
               rd.p = 0:200
               res = sum(p%s(x[1]-obj$scalePar*rd.p,x[4], x[5])*dpois(rd.p,x[2]))}
            )", wdist2)))
            }
      if(dist1 %in% "Nbinom"){wdist1 = "nbinom"
      E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x){
               rd.p = 0:200
               res = sum(p%s(x[1]-obj$scalePar*rd.p,x[4], x[5])*dnbinom(rd.p,mu = x[2], size = x[3]))}
            )", wdist2)))
      }
      if(dist1 %in% "ZIP"){wdist1 = "ZIP"
      E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x){
               rd.p = 0:200
               res = sum(p%s(x[1]-obj$scalePar*rd.p,x[4], x[5])*dZIP(rd.p,x[2], x[3]))}
            )", wdist2)))
      }
      if(dist1 %in% "HP"){wdist1 = "HP"
      E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x){
               rd.p = 0:200
               res = sum(p%s(x[1]-obj$scalePar*rd.p,x[4], x[5])*dHP(rd.p,x[2], x[3]))}
            )", wdist2)))
      }
      if(dist1 %in% "CoMPoisson"){wdist1 = "comp"
      E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x){
               rd.p = 0:200
                                res = sum(p%s(x[1]-obj$scalePar*rd.p,x[4], x[5])*dcomp(rd.p,x[2], x[3]))}
      )", wdist2)))
            }
      if(dist1 %in% "Binom"){wdist1 = "binom"
      E = eval(parse(text = sprintf("apply(cbind(p$y,p$mu1,p$s1,p$mu2,p$s2),1,function(x){
               rd.p = 0:200
               res = sum(p%s(x[1]-obj$scalePar*rd.p,x[4], x[5])*dbinom(rd.p,prob=x[2], size = x[3]))}
            )", wdist2)))
      }


   }

   qres = qnorm(E)

   return(qres)
}

#' @export
meanVarianceConvreg <- function(mu1, sigma1,mu2, sigma2,dist1, dist2) {

   if(dist2 %in% "Gauss"){variance2 = sigma2^2}
   if(dist2 %in% "Lnorm"){variance2 = exp(2*mu2+sigma2)*(exp(sigma2)-1)}
   if(dist2 %in% c("Gamma")){variance2 = mu2 * sigma2^2}

   if(dist2 %in% "Gauss"){exp2 = mu2}
   if(dist2 %in% "Lnorm"){exp2 = exp(mu2 + sigma2^2/2)}
   if(dist2 %in% c("Gamma")){exp2 = mu2*sigma2}

   if(dist1 == "Nbinom"){exp1 = mu1}
   if(dist1 == "Pois"){exp1 = mu1}
   if(dist1 == "ZIP"){exp1 = mu1*(1-sigma1)}
   if(dist1 == "HP"){exp1 = mu1/(1-exp(-mu1)) * (1-sigma1)}
   if(dist1 == "CoMPoisson"){exp1 = Zcomp.mu(mu1,sigma1)}
   if(dist1 == "Binom"){exp1 = mu1*sigma1}
   if(dist1 == "Multinom"){exp1 = mu1*sigma1}

   if(dist1 == "Nbinom"){variance1 = mu1^2/sigma1 + mu1}
   if(dist1 == "Pois"){variance1 = mu1}
   if(dist1 == "ZIP"){variance1 = mu1*(1-sigma1)*(1 + sigma1 * mu1)}
   if(dist1 == "HP"){variance1 = 2*exp(2) / (exp(2)-1) * (1 - (2 / (exp(2)-1))) * (sigma1 * (1 - sigma1))}
   if(dist1 == "CoMPoisson"){variance1 = YY(mu1, sigma1,100) / Z(mu1, sigma1,100) - exp1^2}
   if(dist1 == "Binom"){variance1 = mu1*(1-(mu1/sigma1))}
   if(dist1 == "Multinom"){variance1 = mu1*(1-(mu1/sigma1))}


   mean = exp1 + exp2
   variance = variance1 + variance2


   return(list(mean = mean, variance = variance, e1 = exp1, e2=exp2, s1=variance1, s2=variance2))
}

index.convreg <- function(obj) {
   idx.mu1 = obj$idx[[1]]
   idx.sigma1 = obj$idx[[2]]
   idx.mu2 = obj$idx[[3]]
   idx.sigma2 = obj$idx[[4]]
   return(list(idx.mu1, idx.sigma1, idx.mu2, idx.sigma2))
}

#' Calculate marginal regression
#' @param obj convreg object
#' @param name.reg variable on which calculate the regression
#' @param bp boolean. If factor.
#' @param ... additional parameters
#' @keywords convolution
#' @export
marginalReg <- function(obj, name.reg,bp, ...) {

   idx.var = grep(name.reg ,rownames(obj$estimation))

   idx.data = grep(name.reg ,colnames(obj$data))

   idx.bs = grep("bs", rownames(obj$estimation))

   idx.mu1 = obj$idx[[1]]
   idx.sigma1 =obj$idx[[2]]
   idx.mu2 = obj$idx[[3]]
   idx.sigma2 = obj$idx[[4]]

   theta.mod = matrix(obj$estimation$Estimate, ncol=1)
   theta.mod[idx.var] = 0

   mu1.fixed = colMeans(obj$profile.lik$other[[5]] %*% theta.mod[idx.mu1])
   sigma1.fixed = colMeans(obj$profile.lik$other[[6]] %*% theta.mod[idx.sigma1])
   mu2.fixed = colMeans(obj$profile.lik$other[[7]] %*% theta.mod[idx.mu2])
   sigma2.fixed = colMeans(obj$profile.lik$other[[8]] %*% theta.mod[idx.sigma2])


   theta.var = matrix(0,nrow = length(obj$estimation$Estimate), ncol=1)
   theta.var[idx.var] = obj$estimation$Estimate[idx.var]

   full.X = cbind(obj$profile.lik$other[[5]],obj$profile.lik$other[[6]],obj$profile.lik$other[[7]],obj$profile.lik$other[[8]])

   if(!bp){range.var = matrix(seq(min(full.X[,idx.var]),max(full.X[,idx.var]),by=diff(range(full.X[,idx.var]))/100),ncol=1)}
   if(bp){range.var = diag(length(idx.var))}

   #var.X = matrix(full.X[1,], nrow = nrow(range.var), ncol = ncol(full.X), byrow=T)

   if(!is.null(idx.bs)){
     var.X = matrix(0, nrow = nrow(range.var), ncol = ncol(full.X), byrow=T)
     var.X[,idx.var] = range.var
   }else{
   if(min(idx.var == idx.bs)==1){
      newX = seq(min(obj$data[,idx.data]),max(obj$data[,idx.data]),by=diff(range(obj$data[,idx.data]))/100)
      range.var = matrix(newX, ncol=1)
      bbs = splines::bs(obj$data[,idx.data], degree = length(idx.bs))
      rrange.var = predict(bbs, newX)
      var.X = matrix(0, nrow = length(newX), ncol = ncol(full.X), byrow=T)
      var.X[,idx.var] = rrange.var
   }else{
     var.X = matrix(0, nrow = nrow(range.var), ncol = ncol(full.X), byrow=T)
     var.X[,idx.var] = range.var
   }
}
   X.mu1 = matrix(var.X[,idx.mu1], nrow = nrow(range.var))
   X.sigma1 = matrix(var.X[,idx.sigma1], nrow = nrow(range.var))
   X.mu2 = matrix(var.X[,idx.mu2], nrow = nrow(range.var))
   X.sigma2 = matrix(var.X[,idx.sigma2], nrow = nrow(range.var))

   mu1.var = X.mu1 %*% theta.var[idx.mu1]
   sigma1.var = X.sigma1 %*% theta.var[idx.sigma1]
   mu2.var = X.mu2 %*% theta.var[idx.mu2]
   sigma2.var = X.sigma2 %*% theta.var[idx.sigma2]

   if(bp){mu1.var = c(0, mu1.var)
   sigma1.var = c(0, sigma1.var)
   mu2.var = c(0, mu2.var)
   sigma2.var = c(0, sigma2.var)
   range.var = rbind(0, range.var)}

   mu.1 = obj$transforms[[1]](mu1.fixed+mu1.var)
   sigma.1 = obj$transforms[[2]](sigma1.fixed+sigma1.var)

   mu.2 =  obj$transforms[[3]](mu2.fixed+mu2.var)
   sigma.2 = obj$transforms[[4]](sigma2.fixed+sigma2.var)



   return(list(values = data.frame(x = range.var,mu1= mu.1, s1=sigma.1, mu2=mu.2, s2=sigma.2), names = rownames(obj$estimation)[idx.var]))
}

#' @export
print.testdist <- function(object, ...) {
   cat("--------------------------------------------\n")
   cat("Convolution Regression Distribution Tests","\n")
   cat("Assumed distribution: ", object$dist1," and ",object$dist2,"\n", sep="")
   cat("Tests results \n")
   cat("\n")
   printCoefmat(object$df, digits = 4, P.values = TRUE, has.Pvalue = TRUE)
   cat("\n")
}


