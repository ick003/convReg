qq_ci <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
   q.function <- eval(parse(text = paste0("q", distribution)))
   d.function <- eval(parse(text = paste0("d", distribution)))
   x <- na.omit(x)
   ord <- order(x)
   n <- length(x)
   P <- ppoints(length(x))
   df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
   
   if(is.null(line.estimate)){
      Q.x <- quantile(df$ord.x, c(0.25, 0.75))
      Q.z <- q.function(c(0.25, 0.75), ...)
      b <- diff(Q.x)/diff(Q.z)
      coef <- c(Q.x[1] - b * Q.z[1], b)
   } else {
      coef <- coef(line.estimate(ord.x ~ z))
   }
   
   zz <- qnorm(1 - (1 - conf)/2)
   SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
   fit.value <- coef[1] + coef[2] * df$z
   df$upper <- fit.value + zz * SE
   df$lower <- fit.value - zz * SE
   
   if(!is.null(labels)){ 
      df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
   }
   
   plot(df$z, df$ord.x, type="p", col="black",cex=0.5, main ="QQ-plot", xlab = "Theoretical Quantiles", ylab = "Sample Qauntiles")
   abline(coef[1],coef[2],col=rgb(1,0,0,0.8))
   polygon(c(df$z,rev(df$z)),c(df$upper,rev(df$lower)), col=rgb(1,0,0,0.2), border=NA)
}