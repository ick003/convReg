#' Conway-Maxwell-Poisson GLM
#'
#' Conway-Maxwell-Poisson GLM function
#' @param lamFormula formula for lambda
#' @param nuFormula formula for nu
#' @param data dataset
#' @param weights weights
#' @param lamStart initial values for lambda
#' @param nuStart initial values for nu
#' @param sumTo integer to sum to to calculate the normalization constant
#' @param method optimisation method
#' @param ... additional parameters
#' @keywords comp
#' @export
#' @examples
#' # Not run
glmComp <- function (lamFormula, nuFormula = NULL, data, weights = NULL,lamStart = NULL, 
          nuStart = NULL, sumTo = 100L, method = "BFGS", ...) 
{
   call <- match.call()
   if (missing(data)) {
      data <- environment(lamFormula)
   }
   lamModelFrame <- model.frame(lamFormula, data)
   lamModelTerms <- attr(lamModelFrame, "terms")
   y <- model.response(lamModelFrame)
   nobs <- NROW(y)
   if (any(y + 10 > sumTo)) {
      stop("Response within 10 of 'sumTo' argument")
   }
   xLam <- model.matrix(lamModelTerms, lamModelFrame)
   nLamVariables <- NCOL(xLam)
   if (is.null(nuFormula)) {
      xNu <- matrix(1, nobs, 1L)
      colnames(xNu) <- "(Intercept)"
      nuModelTerms <- NULL
      nuOffset <- rep.int(0L, nobs)
   }
   else {
      nuModelFrame <- model.frame(nuFormula, data)
      nuModelTerms <- attr(nuModelFrame, "terms")
      xNu <- model.matrix(nuModelTerms, nuModelFrame)
      nuOffset <- as.vector(model.offset(nuModelFrame))
   }
   nNuVariables <- NCOL(xNu)
   lamOffset <- as.vector(model.offset(lamModelFrame))
   if (is.null(lamOffset)) {
      lamOffset <- rep.int(0L, nobs)
   }
   poissonModel <- glm.fit(xLam, y, rep.int(1L, nobs), family = poisson(), 
                           offset = lamOffset)
   coefs <- poissonModel$coefficients
   if (!poissonModel$converged) {
      stop("attempt to fit poisson glm with glm.fit failed")
   }
   if (any(is.na(coefs))) {
      badVariables = paste(names(coefs)[is.na(coefs)], collapse = ", ")
      warning(sprintf("initial parameter estimates return NA from glm.fit, dropping variables: %s", 
                      badVariables))
      xLam <- xLam[, !is.na(coefs), drop = FALSE]
   }
   if (is.null(lamStart)) {
      lamStart <- coefs[!is.na(coefs)]
   }
   if (is.null(nuStart)) {
      nuStart <- rep.int(0L, nNuVariables)
      names(nuStart) <- colnames(xNu)
   }
   start <- c(nuStart, lamStart)
   nuIndexes <- seq_len(nNuVariables)
   lamIndexes <- seq_len(nLamVariables) + max(nuIndexes)
   fmin <- function(par) {
      zeta <- par[nuIndexes]
      nu <- exp(xNu %*% zeta + nuOffset)
      beta <- par[lamIndexes]
      lam <- exp(xLam %*% beta + lamOffset)
      pr <- dcomp(y, lam, nu, sumTo)
      return(-sum(weights*log(pr)))
   }
   gmin <- function(par) {
      out <- numeric(length(par))
      zeta <- par[nuIndexes]
      nu <- exp(xNu %*% zeta + nuOffset)
      beta <- par[lamIndexes]
      lam <- exp(xLam %*% beta + lamOffset)
      betaGrad <- t(xLam) %*% ((y - Y(lam, nu, sumTo)/Z(lam, 
                                                       nu, sumTo)) * weights)
      zetaGrad <- t(xNu) %*% (((-log(factorial(y)) + W(lam, nu, sumTo)/Z(lam, nu, sumTo)) * nu) * weights)
      out[lamIndexes] <- betaGrad
      out[nuIndexes] <- zetaGrad
      return(-out)
   }
   res <- optim(start, fmin, gmin, method = method, hessian = TRUE, 
                ...)
   beta <- res$par[lamIndexes]
   zeta <- res$par[nuIndexes]
   niter <- c(f.evals = res$counts[1L], g.evals = res$counts[2L])
   df.residual <- nobs - (nNuVariables + nLamVariables)
   fit <- list(call = call, beta = beta, zeta = zeta, logLik = -res$value, 
               terms = lamModelTerms, nuModelTerms = nuModelTerms, data = data, 
               lamOffset = lamOffset, nuOffset = nuOffset, convergence = res$convergence, 
               niter = niter, hessian = res$hessian, df.residual = df.residual, 
               df.null = nobs - 1L, nobs = nobs)
   class(fit) <- "Comp"
   return(fit)
}