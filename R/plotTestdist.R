#' Convolutive regression function testing plot
#'
#' Convolutive regression function testing plot
#' @param testdist testdist object
#' @keywords convreg
#' @export
#' @examples
#' set.seed(123)
#' e=rlnorm(n=150,mean=0,sd=0.05)
#' x = rnbinom(150,1,0.5)
#' k=5 + x
#' y= data.frame(obs=k + e , f1 = as.factor(x))
#' resTest=testdist(y$obs)
#' plotTestdist(resTest)
plotTestdist <- function(testdist){

   par(mfrow = c(2,2), mar = c(5.1,4.1,2.1,2.1))

   del = diff(range(testdist$x.O))/length(testdist$x.O)/2

   plot(testdist$x.O - del/8, testdist$O, ylab = "counts", xlab = "categories", cex = 1,
        ylim = c(-10, max(testdist$O,testdist$E) + 5), type = "h", lwd=10)
   points(testdist$x.O+del/8, testdist$E, cex=1, col = rgb(0,1,.5), type = "h", lwd = 10)

   weights.chi = (testdist$O-testdist$E)^2 / (testdist$E + 0.001)
   weights.g = log(testdist$O/(testdist$E+0.001)) * testdist$E
   norm = max(diff(range(weights.g)), max(weights.chi))
   weights.g = weights.g / norm
   weights.chi = weights.chi / norm


   abline(h = 0, col = "grey")

   arrows(testdist$x.O - del/8,-10 + weights.chi*10,testdist$x.O-del/8,-10, code = 0, col = "blue", lwd = 2)
   arrows(testdist$x.O + del/8,-5 + weights.g*10,testdist$x.O+del/8,-5, code = 0, col = "red", lwd = 2)


   x.chi = seq(0, min(max(testdist$df$statistic),3*testdist$df$df[1], na.rm=T) + testdist$df$df[1], 0.1)
   y.chi = dchisq(x.chi, df=testdist$df$df[1])
   plot(x.chi,y.chi , type = "l", ylab = "", xlab = "statistic")
   abline(v = testdist$df$statistic[1:2], col = c("blue", "red"), lwd = 2)

   for(i in 1:2){
      coli = c("blue", "red")
      if(min(testdist$df$statistic[i]) > 3*testdist$df$df[1]){
         arrows(.9*max(x.chi),i*max(y.chi)/3,max(x.chi),i*max(y.chi)/3, col = coli[i], length=.1)
      }
   }

   x.sort = sort(testdist$ks.list[[1]])
   plot(x.sort, (1:length(x.sort)) /length(x.sort ), type = "l", ylab = "cdf", xlab = "obs", lwd=2)

   if(testdist$dist1 == "Pois"){
     th.cdf = eval(parse(text = sprintf("p%s%s(x.sort, testdist$ks.list[[2]],testdist$ks.list[[4]], testdist$ks.list[[5]])", testdist$dist1, testdist$dist2)))
   }else{th.cdf = eval(parse(text = sprintf("p%s%s(x.sort, testdist$ks.list[[2]], testdist$ks.list[[3]], testdist$ks.list[[4]], testdist$ks.list[[5]])", testdist$dist1, testdist$dist2)))}

   points(x.sort, th.cdf, type = "l", col = rgb(0,1,.5), lwd = 2)

   delta = 1 / length(x.sort)
   x.ks = seq(0,testdist$df$statistic[3] + 10*delta,by = delta)

   plot(x.ks[2:length(x.ks)], dks(x.ks, length(x.sort)), type = "l", ylab = "pdf", xlab = "statistic")
   abline(v = testdist$df$statistic[3], col = "blue", lwd = 2)

}
