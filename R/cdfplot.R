#' Convolutive regression function cdf plot
#'
#' Convolutive regression function cdf plot
#' @param obj convreg object
#' @param delta discretization parameter
#' @param lwd.pt line width for observations
#' @param cex.pt point size for observations
#' @param cex.text legend size.
#' @param ... additional plotting parameters
#' @keywords convreg
#' @export
#' @examples
#' set.seed(123)
#' e=rnorm(n=150,mean=0,sd=0.05)
#' x = rbinom(150,1,0.5)
#' k=5 + x
#' y= data.frame(obs=k + e , f1 = as.factor(x))
#' res=convreg( ~obs,
#'           formula.mu1 =~ 1,
#'           formula.mu2 =~ f1,
#'           data=y)
#' cdfplot(res)
cdfplot <- function(obj,delta = 0.5,lwd.pt = 2, cex.pt = 1,cex.text = 1, ...)
{

   idx.p1 = index.convreg(obj)[[1]]
   idx.p2 = index.convreg(obj)[[2]]
   idx.p3 = index.convreg(obj)[[3]]
   idx.p4 = index.convreg(obj)[[4]]

   dist1= unlist(strsplit(obj$distname,"/"))[1]
   dist2= unlist(strsplit(obj$distname,"/"))[2]

   y = obj$profile.lik$other[[4]]
   X.mu1 = obj$profile.lik$other[[5]]
   X.sigma1 = obj$profile.lik$other[[6]]
   X.mu2 = obj$profile.lik$other[[7]]
   X.sigma2 = obj$profile.lik$other[[8]]


   if(dist1 != "Multinom"){mu1 = obj$transforms[[1]](X.mu1 %*% obj$estimation$Estimate[idx.p1])}
   if(dist1 == "Multinom"){mu1 = exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.p1], nrow = obj$estimation$Estimate[idx.p2])))) / rowSums( exp(X.mu1 %*% cbind(rep(0,ncol(X.mu1)),t(matrix(obj$estimation$Estimate[idx.p1], nrow = obj$estimation$Estimate[idx.p2])))))}

   sigma1 = obj$transforms[[2]](X.sigma1 %*% obj$estimation$Estimate[idx.p2])
   mu2 = obj$transforms[[3]](X.mu2 %*% obj$estimation$Estimate[idx.p3])
   sigma2 = obj$transforms[[4]](X.sigma2 %*% obj$estimation$Estimate[idx.p4])

   idx.boot = sample(1:nrow(X.mu1), nrow(X.mu1), replace=T)

   if(dist1 != "Multinom"){R.Cmd <- sprintf("sample =matrix(r%s%s(N=nrow(X.mu1)*1000, mu1[idx.boot],sigma1[idx.boot], mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)", dist1, dist2)}
   if(dist1 == "Multinom"){R.Cmd <- sprintf("sample =matrix(r%s%s(N=nrow(X.mu1)*1000, mu1[idx.boot,],sigma1[idx.boot], mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)", dist1, dist2)}

   if(obj$scale){
     if(dist2 %in% "Gauss"){wdist2 = "norm"}
     if(dist2 %in% "Lnorm"){wdist2 = "lnorm"}
     if(dist2 %in% "Gamma"){wdist2 = "gamma"}

     if(dist1 %in% "Pois"){
       wdist1 = "pois"
       txt = "sample =matrix(obj$scalePar*r%s(nrow(X.mu1)*1000, mu1[idx.boot])+ r%s(nrow(X.mu1)*1000,mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)"
       R.Cmd <- sprintf(txt, wdist1, wdist2)
     }
     if(dist1 %in% "Nbinom"){
       wdist1 = "nbinom"
       txt = "sample =matrix(obj$scalePar*r%s(nrow(X.mu1)*1000, mu=mu1[idx.boot],size=sigma1[idx.boot])+ r%s(nrow(X.mu1)*1000,mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)"
       R.Cmd <- sprintf(txt, wdist1, wdist2)
     }
     if(dist1 %in% "ZIP"){
       wdist1 = "ZIP"
       txt = "sample =matrix(obj$scalePar*r%s(nrow(X.mu1)*1000, mu1[idx.boot],sigma1[idx.boot])+ r%s(nrow(X.mu1)*1000,mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)"
       R.Cmd <- sprintf(txt, wdist1, wdist2)
     }
     if(dist1 %in% "HP"){
       wdist1 = "HP"
       txt = "sample =matrix(obj$scalePar*r%s(nrow(X.mu1)*1000, mu1[idx.boot],sigma1[idx.boot])+ r%s(nrow(X.mu1)*1000,mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)"
       R.Cmd <- sprintf(txt, wdist1, wdist2)}
     if(dist1 %in% "CoMPoisson"){
       wdist1 = "comp"
       txt = "sample =matrix(obj$scalePar*r%s(nrow(X.mu1)*1000, mu1[idx.boot],sigma1[idx.boot])+ r%s(nrow(X.mu1)*1000,mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)"
       R.Cmd <- sprintf(txt, wdist1, wdist2)
     }
     if(dist1 %in% "Binom"){
       wdist1 = "binom"
       txt = "sample =matrix(obj$scalePar*r%s(nrow(X.mu1)*1000, prob = mu1[idx.boot],size=sigma1[idx.boot])+ r%s(nrow(X.mu1)*1000,mu2[idx.boot],sigma2[idx.boot]),ncol = 1000, byrow=F)"
       R.Cmd <- sprintf(txt, wdist1, wdist2)
     }

   }

   eval(parse(text=R.Cmd))
   par(mfrow = c(1,1))
   col.band = hcl(seq(0,360,length.out=15),alpha=0.2,c=95,l=30)
   col.est = hcl(seq(0,360,length.out=15),alpha=0.5,c=95,l=60)

   min.seq = min(min(sample), min(y))-delta
   max.seq = max(max(sample), max(y))+delta

   c.h = NULL
   for(j in 1:ncol(sample)){
      c.h = cbind(c.h,cumsum(hist(sample[,j],seq(min.seq,max.seq,delta),plot=F)$density)*delta)
   }
   c.h.f = apply(c.h,1,function(x) quantile(x,probs=c(0.025,0.975)))

   x.coords = hist(sample[,j],seq(min.seq,max.seq,delta),plot=F)$mids

   #dev.off()

   layout(rbind(1,2), heights=c(1,9))

   par(mar=c(0, 0, 0, 0))
   # c(bottom, left, top, right)
   plot.new()
   legend("center","groups", c("actual","population", "predicted"),
          lty = c(1, 1, 1), pch = c(20, NA, 20),col = c("black",col.est[1], col.est[5]),
          lwd = c(1,3,1), bty= "n", ncol = 3, cex = cex.text)
   par(mar=c(5.1, 4.1, 0.1, 2.1))
   plot(ecdf(y),main = "",xlab = "ecdf(x)",... )
   polygon(c(x.coords,rev(x.coords)),
           c(c.h.f[1,],rev(c.h.f[2,])),col = col.band[1],border=NA)
   h.dist=hist(sample,seq(min.seq,max.seq,delta), plot=F)
   points(h.dist$mids, cumsum(h.dist$density)*delta, type="l", col = col.est[1],lwd=lwd.pt, cex=cex.pt)
   plot(ecdf(obj$fitted), add=T, col = col.est[5], cex = cex.pt)

}
