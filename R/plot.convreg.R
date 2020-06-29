#' @export
plot.convreg <- function(x,...)
{
   obj = x
   res = quantRes(obj)
   par(mfrow = c(2,2))
   qq_ci(res)
   hist(res, xlab = "Randomized Quantile Residuals", ...)
   plot(obj$fitted, res, main="Residuals vs Fitted",
        xlab = "Fitted", ylab = "Randomized Quantile Residuals")
   plot(obj$data[,1],res, main="Residuals vs Observed",
        xlab = "Observed", ylab = "Randomized Quantile Residuals")
}
