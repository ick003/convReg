thetaML <- function (y, mu, n = sum(weights), weights, limit = 10, eps = .Machine$double.eps^0.25, 
          trace = FALSE) 
{
   score <- function(n, th, mu, y, w) sum(w * (digamma(th + 
                                                          y) - digamma(th) + log(th) + 1 - log(th + mu) - (y + 
                                                                                                              th)/(mu + th)))
   info <- function(n, th, mu, y, w) sum(w * (-trigamma(th + 
                                                           y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu + 
                                                                                                                 th)^2))
   if (inherits(y, "lm")) {
      mu <- y$fitted.values
      y <- if (is.null(y$y)) 
         mu + residuals(y)
      else y$y
   }
   if (missing(weights)) 
      weights <- rep(1, length(y))
   t0 <- n/sum(weights * (y/mu - 1)^2)
   it <- 0
   del <- 1
   if (trace) 
      message(sprintf("theta.ml: iter %d 'theta = %f'", it, 
                      signif(t0)), domain = NA)
   while ((it <- it + 1) < limit && abs(del) > eps) {
      t0 <- abs(t0)
      del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, 
                                                     mu, y, weights))
      t0 <- t0 + del
      if (trace) 
         message("theta.ml: iter", it, " theta =", signif(t0))
   }
   if (t0 < 0) {
      t0 <- 0
      warning("estimate truncated at zero")
      attr(t0, "warn") <- gettext("estimate truncated at zero")
   }
   if (it == limit) {
      warning("iteration limit reached")
      attr(t0, "warn") <- gettext("iteration limit reached")
   }
   attr(t0, "SE") <- sqrt(1/i)
   t0
}