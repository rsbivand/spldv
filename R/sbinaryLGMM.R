##### Functions for sbinaryLGMM #####

#' @title Estimation of SAR for binary models using Linearized GMM.
#' @description Estimation of SAR model for binary dependent variables (either Probit or Logit), using Linearized GMM estimator suggested by Klier and McMillen (2008). The model is: 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon,
#' }
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise; \eqn{\epsilon \sim N(0, 1)} if \code{link = "probit"} or \eqn{\epsilon \sim L(0, \pi^2/3)} \code{link = "logit"}.
#' 
#' @name sbinaryLGMM
#' @param formula A symbolic description of the model of the form \code{y ~ x | wx} where \code{y} is the binary dependent variable, \code{x} are the independent variables. The variables after \code{|} are those variables that enter spatially lagged: \eqn{WX}. 
#' @param data A \code{data.frame} containing the variables in the model.
#' @param listw A spatial weights matrix of class \code{listw}, \code{matrix}, or \code{Matrix}.
#' @param nins Integer. Order of spatial instruments to generate; default is \code{nins = 2}, such that \eqn{H = (Z, WZ, W^2Z)} are used as instruments.   
#' @param link Character string; distributional assumption on the error term. Either \code{"probit"} (default) or \code{"logit"}. 
#' @param x,object, an object of class \code{binlgmm}.
#' @param digits Integer. The number of digits for \code{summary} methods.
#' @param ... additional arguments.
#' @details 
#' 
#' The steps for the linearized spatial Probit/Logit model are the following:
#'
#' 1. Estimate the model by standard Probit/Logit model, in which spatial autocorrelation and heteroskedasticity are ignored. The estimated values are \eqn{\beta_0}. Calculate the generalized residuals assuming that \eqn{\lambda = 0} and the gradient terms \eqn{G_{\beta}} and \eqn{G_{\lambda}}.
#' 
#' 2. The second step is a two-stage least squares estimator of the linearized model. Thus regress \eqn{G_{\beta}} and \eqn{G_{\lambda}} on \eqn{H = (Z, WZ, W^2Z, ...., W^qZ)} and obtain the predicted values \eqn{\hat{G}}. Then regress \eqn{u_0 + G_{\beta}'\hat{\beta}_0} on \eqn{\hat{G}}. The coefficients are the estimated values of \eqn{\beta} and \eqn{\lambda}.  
#' 
#' The variance-covariance matrix can be computed using the traditional White-corrected coefficient covariance matrix from the last two-stage least squares estimator of the linearlized model. 
#' @examples
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # LGMM for probit using q = 3 for instruments
#' lgmm <- sbinaryLGMM(CRIMED ~ INC + HOVAL | INC,
#'                 link  = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"),
#'                 nins  = 3, 
#'                 data  = COL.OLD)
#' summary(lgmm)
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @keywords models
#' @return An object of class ``\code{bingmm}'', a list with elements:
#' \item{coefficients}{the estimated coefficients,}
#' \item{call}{the matched call,}
#' \item{X}{the X matrix, which contains also WX if the second part of the \code{formula} is used, }
#' \item{H}{the H matrix of instruments used,}
#' \item{y}{the dependent variable,}
#' \item{listw}{the spatial weight matrix,}
#' \item{link}{the string indicating the distribution of the error term,}
#' \item{fit}{an object of \code{lm} representing the T2SLS,}
#' \item{formula}{the formula,}
#' \item{data}{the data,}
#' \item{contrastsX}{the contrasts used in the first part of the formula,}
#' \item{contrastsD}{the contrasts used in the second part of the formula,}
#' \item{Xlevels}{a record of the levels of the factors used in fitting.}
#' @references 
#
#' Klier, T., & McMillen, D. P. (2008). Clustering of auto supplier plants in the United States: generalized method of moments spatial logit for large samples. Journal of Business & Economic Statistics, 26(4), 460-471.
#' 
#' Piras, G., & Sarrias, M. (2023). One or Two-Step? Evaluating GMM Efficiency for Spatial Binary Probit Models. Journal of choice modelling, 48, 100432. 
#' 
#' Piras, G,. & Sarrias, M. (2023). GMM Estimators for Binary Spatial Models in R. Journal of Statistical Software, 107(8), 1-33. 
#' @seealso \code{\link[spldv]{sbinaryGMM}}, \code{\link[spldv]{impacts.bingmm}}.
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update))
#' @import stats methods Formula
#' @importFrom sphet listw2dgCMatrix
#' @export 
sbinaryLGMM <- function(formula, data, 
                        listw = NULL, 
                        nins  = 2, 
                        link  = c("logit", "probit"), 
                        ...)
{
  # Based on McMillen's code 
  link  <- match.arg(link)
  
  # Spatial weight matrix (W): as CsparseMatrix
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("Neighbourhood list or listw format unknown")
  if(inherits(listw,"listw"))   W    <- sphet::listw2dgCMatrix(listw)	
  if(inherits(listw,"matrix"))  W    <- Matrix(listw)	
  if(inherits(listw,"Matrix"))  W    <- listw	
  
  # Model frame
  callT    <- match.call(expand.dots = TRUE)
  callF    <- match.call(expand.dots = FALSE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  f1       <- Formula(formula)
  if (length(f1)[2L] == 2L) Durbin <- TRUE else Durbin <- FALSE 
  mf$formula <- f1 
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  nframe   <- length(sys.calls())
  
  # Get variables and run some checks
  y  <- model.response(mf)
  if (anyNA(y)) stop("NAs in dependent variable")
  if (!all(y %in% c(0, 1, TRUE, FALSE))) stop("All dependent variables must be either 0, 1, TRUE or FALSE")
  if (!is.numeric(y)) y <- as.numeric(y)
  X  <- model.matrix(f1, data = mf, rhs = 1)
  # Added for prediction
  contrastsX <- attr(X, "contrasts")
  Xlevels    <- .getXlevels(attr(mf, "terms"), mf)
  if (Durbin){
    x.for.w    <- model.matrix(f1, data = mf, rhs = 2)
    contrastsD <- attr(x.for.w, "contrasts")
    name.wx    <- setdiff(colnames(x.for.w), "(Intercept)")
    
    if (!all(name.wx %in% colnames(X))) 
      warning("Some variables in WX do not appear in X. Check the formula if this is not intended.")
    
    WX           <- W %*% x.for.w[, name.wx, drop = FALSE]
    colnames(WX) <- paste0("lag_", name.wx)
    if (anyNA(WX)) stop("NAs in WX variable")
    X <- cbind(X, WX)
  } else {
    contrastsD <- NULL
  }
  if (anyNA(X)) stop("NAs in independent variables")
  N  <- nrow(X)
  K  <- ncol(X)
  sn <- nrow(W)
  if (N != sn) stop("Number of spatial units in W is different to the number of data")
  
  ## Generate initial instruments
  Z <- make.instruments(W, x = X, q = nins)
  H <- cbind(X, Z)
  # Let linearly independent columns
  H <- H[, qr(H)$pivot[seq_len(qr(H)$rank)]]
  P <- ncol(H)
  if (P < K) stop("Underspecified model")
  if (any(is.na(H))) stop("NAs in the instruments")
  
  ## Starting values
  sbinary <- glm(y ~ as.matrix(X) - 1, family = binomial(link = link), data = mf)
  b_init  <- sbinary$coef # Initial values of beta from a standard binary model
  
  ## First step is standard binary of y on X and compute the generalized residuals 
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  ddfun <- switch(link,
                  "logit"  = function(x) (1 - 2 * pfun(x)) * pfun(x) * (1 - pfun(x)),
                  "probit" = function(x) -x * dnorm(x))  
  ai   <- as.vector(X %*% b_init)
  q    <- 2*y - 1
  fa   <- pmax(dfun(q*ai), .Machine$double.eps)
  Fa   <- pmax(pfun(q*ai), .Machine$double.eps)
  ffa  <- ddfun(q*ai) 
  u0   <- q * (fa/Fa)
  grad <- -1 * as.vector(q^2 * ((ffa * Fa - fa^2) / (Fa^2))) # Common vector of the derivative
  Gb   <- grad * X
  Glam <- grad * (W %*% ai)
  G    <- cbind(Gb, Glam)
  
  # Second step
  Ghat        <- H %*% solve(crossprod(H)) %*% crossprod(H, G) #Predicted values of G
  epsilon     <- u0 + as.vector(Gb %*% b_init)
  fit         <- lm(epsilon ~ as.matrix(Ghat) + 0)
  bhat        <- coef(fit)
  names(bhat) <- c(colnames(X), "lambda") 
  
  ## Saving results
  out <- structure(
    list(
      coefficients = bhat, 
      call         = callF,
      X            = X, 
      H            = H, 
      y            = y, 
      listw        = W,
      link         = link, 
      fit          = fit, 
      formula     = f1,
      mf          = mf, 
      data        = data,
      contrastsX  = contrastsX,
      contrastsD  = contrastsD,
      Xlevels     = Xlevels
    ), 
    class = "binlgmm"
  )
  out
}

#### S3 methods ----
#' @rdname sbinaryLGMM
#' @export
coef.binlgmm <- function(object, ...){
  object$coefficients
}


#' @rdname sbinaryLGMM
#' @method vcov binlgmm
#' @importFrom car hccm
#' @export 
vcov.binlgmm <- function(object, ...){
  V <- car::hccm(object$fit)
  colnames(V) <- rownames(V) <- names(coef(object))
  return(V)
}

#' Get Model Summaries for use with "mtable" for objects of class binlgmm
#' 
#' A generic function to collect coefficients and summary statistics from a \code{binlgmm} object. It is used in \code{mtable}
#' 
#' @param obj a \code{binlgmm} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @return A list with an array with coefficient estimates and a vector containing the model summary statistics. 
#' @importFrom memisc getSummary
#' @method getSummary binlgmm
#' @export 
getSummary.binlgmm <- function(obj, alpha = 0.05, ...){
  if (inherits(obj, c("summary.binlgmm"))){
    coef <- obj$CoefTable
  } else {
    smry <- summary(obj)
    coef <- smry$Coef
  }
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  nrow(obj$X)
  sumstat <- c(N = N)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

#' @rdname sbinaryLGMM
#' @method print binlgmm
#' @export 
print.binlgmm <- function(x, 
                          digits = max(3, getOption("digits") - 3),
                          ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print.default(format(drop(x$coefficients), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname sbinaryLGMM
#' @method summary binlgmm
#' @export
summary.binlgmm <- function(object, ...){
  b <- coef(object$fit)
  names(b) <- names(object$coefficients)
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.binlgmm", "binlgmm")
  return(object)
}

#' @rdname sbinaryLGMM
#' @method print summary.binlgmm
#' @export
print.summary.binlgmm <- function(x,
                                  digits = max(3, getOption("digits") - 2),
                                  ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      SLM Binary Model by Linearized GMM \n")
  cat("        ------------------------------------------------------------\n")
  
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  
  cat(paste("\nSample size:", signif(nrow(x$X), digits)), "\n")
  invisible(x)
}


#' Predictions for Spatial Binary LGMM Models
#'
#' Computes predicted probabilities for spatial binary response models estimated via LGMM. 
#' Supports both probit and logit links, accounts for spatial heteroskedasticity, and optionally 
#' returns standard errors using the Delta method.
#'
#' @param object An object of class \code{binlgmm} returned by a spatial binary GMM estimation function.
#' @param newdata An optional data frame in which to look for variables with which to predict. 
#' If omitted, the original data used to fit the model is used.
#' @param Sinv Optional user-supplied spatial multiplier matrix \eqn{S = (I - \lambda W)^{-1}}. 
#' If \code{NULL}, it is computed using the spatial weight matrix.
#' @param het Logical. If \code{TRUE}, assumes a heteroskedastic error structure with spatially varying variances.
#' @param approximation Logical. If \code{TRUE}, uses power-series approximation to compute the inverse spatial matrix.
#' @param pw Integer. Power-order to use when \code{approximation = TRUE}.
#' @param ses Logical. If \code{TRUE}, standard errors of the predictions are computed using the Delta method.
#' @param theta Optional parameter vector (including \code{lambda}) to use for prediction instead of the estimated one.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The function computes predicted probabilities \eqn{\hat{p}_i = F(a_i)} where \eqn{a_i} is a spatially filtered linear index. 
#' In the presence of heteroskedasticity (\code{het = TRUE}), the normalization involves the row-wise standard deviation 
#' of the spatial multiplier. When \code{ses = TRUE}, standard errors are computed using the analytical Jacobian of the 
#' prediction function with respect to the parameters and the estimated variance-covariance matrix.
#'
#' @return A numeric vector of predicted probabilities if \code{ses = FALSE}. If \code{ses = TRUE}, returns a matrix with:
#' \describe{
#'   \item{\code{p_hat}}{Predicted probabilities.}
#'   \item{\code{Std. error}}{Standard errors of the predictions.}
#'   \item{\code{z value}}{Z-statistics.}
#'   \item{\code{Pr(> z)}}{Two-sided p-values.}
#' }
#'
#' @seealso \code{\link{sbinaryLGMM}}, \code{\link{vcov.binlgmm}}
#'
#' @examples
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # Estimate the model
#' lgmm <- sbinaryLGMM(CRIMED ~ INC + HOVAL | INC,
#'                     link  = "probit",
#'                     listw = spdep::nb2listw(COL.nb, style = "W"),
#'                     nins  = 3,
#'                     data  = COL.OLD)
#' 
#' # Predicted probabilities with SES
#' out <- predict(lgmm, ses = TRUE)
#' head(out, 5)
#' 
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @keywords prediction
#' @export 
#' @method predict binlgmm
predict.binlgmm <- function(object, 
                           newdata, 
                           Sinv = NULL,
                           het  = TRUE,
                           approximation = FALSE,
                           pw  = 5,
                           ses = FALSE,
                           theta = NULL, ...){
  
  if (!inherits(object, "binlgmm")) warning("calling predict.bingmm(<fake-bingmm-object>) ...")
  # Obtain data from formula
  
  # Extract formula and spatial weight matrix
  f1     <- object$formula
  Durbin <- (length(f1)[2L] == 2L)
  W      <- object$listw
  
  # Generate model frame and matrices
  if (missing(newdata) || is.null(newdata)){
    mf <- model.frame(f1,  data = object$data)
    X  <- model.matrix(f1, data = mf, rhs = 1)
  } else {
    # Generate model frame with new data
    mf <- model.frame(f1, newdata, xlev = object$Xlevels)
    X  <- model.matrix(f1, data = mf, rhs = 1, contrasts.arg = object$contrastsX)
  }
  if (Durbin){
    x.for.w      <- model.matrix(f1, data = mf, rhs = 2, contrasts.arg = object$contrastsD)
    name.wx      <- colnames(x.for.w)
    WX           <- crossprod(t(W), x.for.w)
    name.wx      <- name.wx[which(name.wx != "(Intercept)")]
    WX           <- WX[ , name.wx, drop = FALSE] # Delete the constant from the WX
    colnames(WX) <- paste0("lag_", name.wx)
    if (any(is.na(WX))) stop("NAs in WX variable")
    X <- cbind(X, WX)
  }
  
  # Get parameters
  n         <- nrow(X)
  theta.hat <- if(is.null(theta)) coef(object) else theta
  lambda    <- theta.hat["lambda"]
  betas     <- theta.hat[which(names(theta.hat) != "lambda")]
  
  # Generate link
  link <- object$link
  pfun <- switch(link, "probit" = pnorm, "logit"  = plogis)
  if (ses) dfun <- switch(link, "probit" = dnorm, "logit" = dlogis)
  
  
  # Generate S matrix or use S provided by user
  if (is.null(Sinv)) {
    if (approximation) {
      Sinv <- app_W(W, lambda, pw)
    } else {
      A    <- Matrix::Diagonal(n) - lambda * W
      Sinv <- Matrix::solve(A)
    }
  }
  
  # Generate linear index
  Xb  <- drop(X %*% betas)
  Sxb <- drop(Sinv %*% Xb)
  if (het){
    rownorms  <- sqrt(Matrix::rowSums(Sinv ^ 2))
    sigma_inv <- 1 / rownorms
    a         <- sigma_inv * Sxb
  } else {
    a  <- Sxb
  }
  
  # Predictions
  pred <- pfun(a)
  if (!ses) return(pred)
  
  # Standard errors via Delta method
  dfa <- dfun(a)
  
  # Derivative w.r.t beta
  if (het){
    SinvX        <- as.matrix(Sinv %*% X)
    SinvX_scaled <- SinvX * sigma_inv
    der_beta     <- dfa * SinvX_scaled  
  } else {
    SinvX        <- as.matrix(Sinv %*% X)
    der_beta     <- dfa * SinvX
  }
  
  
  # Derivative w.r.t lambda
  if (het){
    BW        <- Sinv %*% W
    BBt       <- Matrix::tcrossprod(Sinv)
    diag_term <- 2 * Matrix::rowSums(BW * BBt)
    Drho      <- - 0.5 * sigma_inv^3 * diag_term
    term1     <- Drho * Sxb
    
    term2           <- (Sinv %*% W %*% Sinv %*% Xb) * sigma_inv
    der_lambda      <- dfa * (term1 + term2)
  } else {
    der_lambda      <- dfa * drop(Sinv %*% W %*% Sinv %*% Xb)
  }
  
  
  # Combine Jacobian
  Jac        <- cbind(der_beta, der_lambda)
  
  # Compute VCOV and SE
  V   <- vcov(object)
  se  <- sqrt(rowSums((Jac %*% V) * Jac)) # Efficient diag(JVJ')
  
  # Return full prediction table
  z    <- pred / se
  pval <- 2 * pnorm(-abs(z))
  
  return(cbind(`p_hat` = pred, `Std. error` = se, `z value` = z, `Pr(> z)` = pval))
}

