#' Convolutive regression function testing
#'
#' Convolutive regression function testing
#' @param y.obs observations vector
#' @param dist1 chain of character to idnetify the distribution of variable1
#' @param dist2 chain of character to idnetify the distribution of variable2
#' @param ... additional parameters for convreg().
#' @keywords convreg test
#' @export
#' @examples
#' set.seed(123)
#' e=rnorm(n=150,mean=0,sd=0.05)
#' x = rbinom(150,1,0.5)
#' k=5 + x
#' y= data.frame(obs=k + e , f1 = as.factor(x))
#' resTest=testdist(y$obs)
testdist <- function(y.obs, dist1 = "Nbinom", dist2 = "Gauss",...) {

   n = length(y.obs)
   if(n > 500){n = 500}
   dat.temp = data.frame(obs = y.obs[sample(1:length(y.obs), n)])
   cv = convreg(formula.resp = ~obs, dist1 = dist1, dist2 = dist2, data = dat.temp, ...)

   #### Testing the chosen distribution ####

   # Chi-2 continuous

   res = y.obs
   #res = round(res)
   by.o = diff(range(res)) / (length(c(res))/20)
   seq.o = seq(min(res)-by.o, max(res)+by.o, by=by.o)
   O = sapply(seq.o, function(x) sum(res <= x))
   O[2:length(seq.o)] = O[2:length(seq.o)] - O[1:(length(seq.o)-1)]
   idx.O = which(O > 2)
   O = O[idx.O]
   x.o = seq.o[idx.O]
   #E = ppois(seq.o,mean(res)) * length(res)
   mu1 = cv$transform[[1]](coef(cv)$theta_mu1)
   sigma1 = cv$transform[[2]](coef(cv)$theta_sigma1)
   mu2 = cv$transform[[3]](coef(cv)$theta_mu2)
   sigma2 = cv$transform[[4]](coef(cv)$theta_sigma2)
   if(dist1 == "Pois"){
     E =  eval(parse(text = sprintf("p%s%s(seq.o, mu1, mu2, sigma2)", dist1, dist2)))
   }else{
     E =  eval(parse(text = sprintf("p%s%s(seq.o, mu1, sigma1,mu2, sigma2)", dist1, dist2)))
   }

   E = E * length(res)
   E[2:length(seq.o)] = E[2:length(seq.o)] - E[1:(length(seq.o)-1)]
   E = E[idx.O]

   chi.stat = sum((O[E>0]-E[E>0])^2 / E[E>0])
   df.chi = max(length(E[E>0])-4,1)
   p.val.chi <- 1-pchisq(chi.stat,df = df.chi)

   # G - testing

   g.stat = sum( 2 * O * log(O / E))
   p.val.g <- 1-pchisq(g.stat,df = df.chi)

   # K-S

   if(dist1 == "Pois"){
     ks.stat =  eval(parse(text = sprintf("ks.test(res, p%s%s, mu1, mu2, sigma2)", dist1, dist2)))
   }else{
     ks.stat =  eval(parse(text = sprintf("ks.test(res, p%s%s, mu1, sigma1,mu2, sigma2)", dist1, dist2)))
      }


   test.df = data.frame("statistic" = c(chi.stat, g.stat, ks.stat$statistic), "df" = c(df.chi, df.chi, NA), "p-value" = c(p.val.chi, p.val.g, ks.stat$p.value))

   rownames(test.df) <- c("Chi-2 test", "G-test", "KS test")

   RET = list(df = test.df,E = E,O = O,x.O = x.o, dist1 = dist1, dist2 = dist2, ks.list = list(res,mu1,sigma1,mu2,sigma2))

   class(RET) <- 'testdist'
   return(RET)
}
