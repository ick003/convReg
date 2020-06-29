#' Convolutive regression function
#'
#' Convolutive regression function
#' @param formula.resp response
#' @param formula.mu1 formula for variable1 expectation
#' @param formula.sigma1 formula for variable1 variance
#' @param formula.mu2 formula for variable2 expectation
#' @param formula.sigma2 formula for variable2 variance
#' @param data dataset
#' @param dist1 chain of character to idnetify the distribution of variable1
#' @param dist2 chain of character to idnetify the distribution of variable2
#' @param fixed set the parameters that are fixed to some specified values
#' @param na.action action for NA values
#' @param method estimation method
#' @param debug getting into the function
#' @param quiet display comments.
#' @param weights weights.
#' @param scale Should the observations be scales
#' @param scaleInit Initial value of the scaling.
#' @param ... additional parameters
#' @keywords convreg
#' @export
#' @examples
#' set.seed(123)
#' e=0.75+rnorm(n=500,mean=0,sd=0.25)
#' x1 = rnorm(500,0,0.5)
#' x2 = runif(500,-2,1)
#' k= rnbinom(500,mu=exp(x2),size=1.5)
#' y= data.frame(obs=(k + e) , f1 = x1, f2 = x2)
#' par(mfrow = c(1,1), mar = c(5,4,4,2))
#' hist(y$obs, breaks = seq(min(y$obs)-1.5, max(y$obs)+3, 3), xlim = c(0,14*24))
#' res.reg.em =convreg( ~obs,
#'           formula.mu1 =~ f2,
#'           formula.mu2 =~ f1,
#'           data=y,
#'           dist1 = "Nbinom",
#'           method = "mle")
#' res.reg.em
#' summary(res.reg.em)
#' distplot(res.reg.em, delta=3, xlim=c(0,200))
convreg = function(formula.resp,
                   formula.mu1= ~ 1,
                   formula.sigma1= ~ 1,
                   formula.mu2 = ~ 1,
                   formula.sigma2 = ~ 1,
                   data,
                   dist1="Nbinom",                        # nass density
                   dist2="Gauss",
                   fixed = NULL,
                   na.action,
                   method = "mle",
                   debug=FALSE,
                   quiet = TRUE,
                   weights = NULL,
                   scale = FALSE,
                   scaleInit = 1,
                   ...){


   # case of mix of dist1 and dist2
   # case: dist1 == dist2 , stop with error messag e
   # case: dist1 != dist2 , but not in the dist.

   # method = "mle"
   #          "pen"
   #require(maxLik, CompGLM, compoisson)

   call <- match.call()

   if (missing(data)) {
      #data <- environment(formula.mu1)
      data <- environment(formula.resp)
   }

   transforms = rep(NA,4)

   # automatically asign to transformation type for specific distribution
   if(dist1 %in% c("Nbinom","CoMPoisson","Pois")){transforms[1:2] = c("pos","pos")}
   if(dist1 %in% c("Binom")){transforms[1:2] = c("01","id")}
   if(dist1 %in% c("Multinom")){transforms[1:2] = c("pos","id")}
   if(dist1 %in% c("ZIP")){transforms[1:2] = c("pos","01")}
   if(dist1 %in% c("HP")){transforms[1:2] = c("pos","01")}

   if(dist2 %in% c("Gauss","Lnorm")){transforms[3:4] = c("id","pos")}
   if(dist2 %in% c("Gamma")){transforms[3:4] = c("pos","pos")}

   data[which(is.na(data),arr.ind=T)] = 0

   keep.vars = c(all.vars(formula.resp),all.vars(formula.mu1),all.vars(formula.sigma1),all.vars(formula.mu2),all.vars(formula.sigma2))

   if(length(keep.vars)>1){data = data[,which(colnames(data) %in% keep.vars)]
   }else{data=as.data.frame(data.matrix(model.frame(formula.resp, data = data, na.action=NULL)))}

   if(is.null(weights)){weights = rep(1, nrow(data))}

   if(debug){browser()}

   y = data.matrix(model.frame(formula.resp, data = data, na.action=NULL))

   if(dist1 == "Binom" & length(which(fixed$name == "sigma 1: (Intercept)")) == 0){
      fixed$name = c(fixed$name, "sigma 1: (Intercept)")
      fixed$value = c(fixed$value, max(1,ceiling(y)))
   }

   # parameter estimation of 4

   X.mu1.model.frame = model.frame(formula.mu1,data = data,na.action=NULL, drop.unused.levels = TRUE)
   X.mu1.model.terms = attr(X.mu1.model.frame,"terms")
   X.mu1 = model.matrix(X.mu1.model.terms, X.mu1.model.frame)

   X.sigma1.model.frame = model.frame(formula.sigma1,data = data,na.action=NULL, drop.unused.levels = TRUE)
   X.sigma1.model.terms = attr(X.sigma1.model.frame,"terms")
   X.sigma1 = model.matrix(X.sigma1.model.terms, X.sigma1.model.frame)

   X.mu2.model.frame = model.frame(formula.mu2,data = data,na.action=NULL, drop.unused.levels = TRUE)
   X.mu2.model.terms = attr(X.mu2.model.frame,"terms")
   X.mu2 = model.matrix(X.mu2.model.terms, X.mu2.model.frame)

   X.sigma2.model.frame = model.frame(formula.sigma2,data = data,na.action=NULL, drop.unused.levels = TRUE)
   X.sigma2.model.terms = attr(X.sigma2.model.frame,"terms")
   X.sigma2 = model.matrix(X.sigma2.model.terms, X.sigma2.model.frame)


      # link to the parameter
   theta.mu1 = rep(0.1,ncol(X.mu1))
   theta.sigma1 = rep(0.1,ncol(X.sigma1))
   theta.mu2 = rep(0.1,ncol(X.mu2))
   theta.sigma2 = rep(0.1,ncol(X.sigma2))


      # Setup constraints on parameters
   fun1 = eval(parse(text=paste("f",transforms[1],sep="")))
   fun2 = eval(parse(text=paste("f",transforms[2],sep="")))
   fun3 = eval(parse(text=paste("f",transforms[3],sep="")))
   fun4 = eval(parse(text=paste("f",transforms[4],sep="")))
   if(!quiet){message(sprintf('links: [formula1: %s sigma1: %s, formula2 : %s sigma2 : %s]',
                   transforms[1],transforms[2],transforms[3],transforms[4]))}

   # setup initialization for the optimization algo
   # in future, we gonna do more comb of distribbutio.


   if(dist1 == "Pois"){fixed$name = c(fixed$name,"sigma 1: (Intercept)")
   fixed$value = c(fixed$value,0)}

   theta.name = c(paste("mu 1: ",colnames(X.mu1),sep=""),
                  paste("sigma 1: ",colnames(X.sigma1),sep=""),
                  paste("mu 2: ",colnames(X.mu2),sep=""),
                  paste("sigma 2: ",colnames(X.sigma2),sep=""))

   init.theta = initthetaConvreg(y,dist1,dist2,formula.mu1, formula.sigma1, formula.mu2, formula.sigma2,
                                   fixed,data,theta.name, scale, scaleInit, debug)

   idx.fixed = as.numeric(na.omit(match(fixed$name,theta.name)))
   idx.var   = which(!(theta.name %in% fixed$name))

   listParam = list(init.theta=init.theta,data = data,
                    formula.resp = formula.resp, formula.mu1 = formula.mu1, formula.sigma1 = formula.sigma1, formula.mu2 = formula.mu2, formula.sigma2 = formula.sigma2,
                    dist1 = dist1, dist2 = dist2,weights = weights, quiet = quiet, debug = debug,
                    X.mu1 = X.mu1, X.sigma1 = X.sigma1, X.mu2 = X.mu2, X.sigma2 = X.sigma2,y=y,
                    fun1 = fun1, fun2 = fun2, fun3=fun3,fun4=fun4,
                    idx.fixed = idx.fixed, idx.var = idx.var,fixed = fixed,
                    scale = scale, scaleInit = scaleInit)

   res = fitConvreg(listParam,method)

   opt.fin = res$opt

   profile.lik = res$prof

   Hessian = opt.fin$hessian

   theta.est = matrix(res$theta, ncol=1)

   attr(theta.est, "names") <- names(res$theta)

   if(dist1 != "Multinom"){

      rownames(theta.est) = attr(theta.est,"names")

      idx.mu1 = seq_len(ncol(X.mu1));
      idx.sigma1 = seq_len(ncol(X.sigma1)) + max(idx.mu1);
      idx.mu2 = seq_len(ncol(X.mu2)) + max(idx.sigma1);
      idx.sigma2 = seq_len(ncol(X.sigma2)) + max(idx.mu2);

      # estimation of distution parameter

      mu1 = as.matrix(fun1(X.mu1 %*% theta.est[idx.mu1]))

      sigma1 = as.matrix(fun2(X.sigma1 %*% theta.est[idx.sigma1]))


      if(dist1 == "HP"){sigma1 = 1-sigma1}

      mu2 = as.matrix(fun3(X.mu2 %*% theta.est[idx.mu2]))

      sigma2 = as.matrix(fun4(X.sigma2 %*% theta.est[idx.sigma2]))

      if(dist1 == "Nbinom" | dist1 == "Pois"){exp1 = mu1}
      if(dist1 == "Binom"){exp1 = mu1*sigma1}
      if(dist1 == "CoMPoisson"){exp1 = Zcomp.mu(mu1,sigma1)}
      if(dist1 == "ZIP"){exp1 = mu1 * (1 - sigma1)}
      if(dist1 == "HP"){exp1 = mu1 / (1 - exp(-mu1)) * (1-sigma1)}

      if(dist2 == "Gauss"){exp2 = mu2}
      if(dist2 == "Lnorm"){exp2 = exp(mu2 + sigma2^2/2)}

      y.fitted = exp1 + exp2

      if(scale){y.fitted = res$scalePar*exp1 + exp2}

      fisher_info <- ginv(-Hessian)
      s.e = rep(0,length(theta.est))
      s.e[idx.var]<-sqrt(diag(fisher_info))

      upper<-theta.est+1.96*s.e
      lower<-theta.est-1.96*s.e

      df.residuals = length(y) - length(idx.var)
      df.null = length(y) -1
      niter <- c(f.evals = opt.fin$counts[1L], g.evals = opt.fin$counts[2L])

      par.est = data.frame(Estimate = theta.est, Std.Error = s.e,
                           Low = lower, Upp = upper,
                           t.value = theta.est / s.e,
                           p.value = 2.0 * pt(-abs(theta.est / s.e), df = df.residuals),
                           row.names = rownames(theta.est))

   }

      xdata <- data.frame(
         time=1:length(y),
         y=y,
         y.fitted,
         X.mu1,
         X.sigma1,
         X.mu2,
         X.sigma2)

      names(xdata)[2] <- 'y'

      RET<-list(estimation=par.est, fitted = y.fitted, xdata=xdata, data = data, theta.i = c(init.theta$mu1, init.theta$sigma1,init.theta$mu2, init.theta$sigma2),
                idx = list(init.theta$idx1,init.theta$idx1b,init.theta$idx2,init.theta$idx2b),
                #residuals = residuals,
                hessian=Hessian, loglik = opt.fin$value,
                transforms = list(fun1,fun2,fun3,fun4), df.residuals = df.residuals,
                niter = niter,distname = paste(dist1, dist2, sep="/"), df.null = df.null,
                method = method, formulas = list(mu1 = formula.mu1, sigma1 = formula.sigma1,mu2 = formula.mu2,sigma2 = formula.sigma2),
                aic = 2*length(theta.est) - 2*opt.fin$value,code = opt.fin$convergence, profile.lik = profile.lik,
                message = opt.fin$message, call = call,scale = scale,
                scalePar = res$scalePar, model.terms = list(X.mu1.model.terms, X.sigma1.model.terms, X.mu2.model.terms, X.sigma2.model.terms))


   class(RET) <- 'convreg'
   if(!quiet){message (" ** done")}
   return(RET)

}
