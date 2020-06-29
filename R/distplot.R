#' Convolutive regression function distribution plot
#'
#' Convolutive regression function distribution plot
#' @param obj convreg object
#' @param delta discretization parameter
#' @param ylim plotting parameter
#' @param lwd.pt line width for observations
#' @param cex.pt point size for observations
#' @param cex.text legend size.
#' @param ... additional plotting parameters
#' @keywords distribution plot
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
#' distplot(res, cex.lab = 0.5, cex.axis = 0.5, lwd.pt = 0.5, cex.pt = 0.5, cex.text = 0.5)
distplot <- function(obj,delta = 0.5,ylim=NULL,lwd.pt = 2, cex.pt = 1, cex.text = 1, ...)
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
   if(dist1 == "HP"){sigma1 = 1 - sigma1}
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

   min.seq = min(min(sample), min(y), min(obj$fitted))-delta
   max.seq = max(max(sample), max(y), max(obj$fitted))+delta

   par(mfrow = c(1,1))
   col.est = hcl(seq(0,360,length.out=15),alpha=0.8,c=95,l=50)
   h = hist(y, breaks = seq(min.seq,max.seq,delta), plot=F)
   h.fit=hist(obj$fitted,seq(min.seq,max.seq,delta),plot=F)
   h.dist=hist(sample,seq(min.seq,max.seq,delta), plot=F)

   #dev.off()

   layout(rbind(1,2), heights=c(1,9))

   par(mar=c(0, 0, 0, 0))
   # c(bottom, left, top, right)
   plot.new()
   legend("center","groups", c("actual","population", "predicted"),
          lty = c(NA, 1, NA), pch = c(NA, 1, 1),col = c("black",col.est[1], col.est[5]),
          bty= "n", ncol = 3, lwd = c(NA,lwd.pt,NA),
          cex = cex.text,
          fill = c("lightgrey", NA, NA), border = c("black",NA,NA))
   par(mar=c(5.1, 4.1, 0.1, 2.1))
   if(is.null(ylim)){ylim = c(0,max(h$density,h.fit$density,h.dist$density))}
   plot(h,ylim = ylim,freq=F,col="lightgrey",
        main = "", xlab = "Values", ...)

   points(h.fit$mids, h.fit$density, type="p", col = col.est[5],lwd = lwd.pt, cex = cex.pt)
   points(h.dist$mids, h.dist$density, type="l", col = col.est[1],lwd = lwd.pt, cex = cex.pt)
   points(h.dist$mids, h.dist$density, type="p", col = col.est[1],lwd = lwd.pt, cex = cex.pt)
}
