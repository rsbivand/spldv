##### Functions for sbinaryGMM ####


#' @title GMM Estimation for Binary Spatial Autoregressive Models.
#' 
#' @description Estimates a Spatial Autoregressive (SAR) model for binary outcomes using the Generalized Method of Moments (GMM). The model supports both Probit and Logit links, with one-step or two-step GMM estimation. The structural equation is:
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon,
#' }
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise.
#' The error term \eqn{\varepsilon} follows a standard normal distribution (\code{link = "probit"}) or logistic distribution (\code{link = "logit"}).
#' 
#' @name sbinaryGMM
#' @param formula A symbolic description of the model of the form \code{y ~ x | wx} where \code{y} is the binary dependent variable, \code{x} are the independent variables. The variables after \code{|} are those variables that enter spatially lagged: \eqn{WX}. 
#' @param data A \code{data.frame} containing the variables in the model.
#' @param listw A spatial weights matrix of class \code{listw}, \code{matrix}, or \code{Matrix}.
#' @param nins Integer. Order of spatial instruments to generate; default is \code{nins = 2}, such that \eqn{H = (Z, WZ, W^2Z)} are used as instruments.  
#' @param link Character string; distributional assumption on the error term. Either \code{"probit"} (default) or \code{"logit"}. 
#' @param winitial Character string. A string indicating the initial moment-weighting matrix \eqn{\Psi}; it can be either \code{winitial = "optimal"} (the default, \eqn{(H'H/n)^{-1}}) or \code{winitial = "identity"}. 
#' @param s.matrix Character string. Only valid of \code{type = "twostep"} is used. This is a string indicating the type of variance-covariance matrix \eqn{\hat{S}} to be used in the second-step procedure; it can be \code{s.matrix = "robust"} (the default) or \code{s.matrix = "iid"}.
#' @param type Character string. A string indicating whether the one-step (\code{type = "onestep"}), or two-step GMM (\code{type = "twostep"}) should be computed.
#' @param gradient Logical. Only for testing procedures. Should the analytic gradient be used in the GMM optimization procedure? \code{TRUE} as default. If \code{FALSE}, then the numerical gradient is used. 
#' @param start Vector of optional starting values. If not \code{NULL}, the user must provide a vector of initial parameters for the optimization procedure. When \code{start = NULL}, \code{sbinaryGMM} uses the traditional Probit or Logit estimates as initial values for the parameters, and the correlation between \eqn{y} and \eqn{Wy} as initial value for \eqn{\lambda}.
#' @param cons.opt Logical. Should a constrained optimization procedure for \eqn{\lambda} be used? \code{FALSE} as default.  
#' @param approximation Logical. If \code{TRUE}, then \eqn{(I - \lambda W)^{-1}} is approximated as \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is \code{FALSE}.
#' @param pw Integer. The power used for the approximation \eqn{I + \lambda W + \lambda^2 W^2 + \lambda^3 W^3 + ... +\lambda^q W^q}. The default is 5.
#' @param tol.solve Tolerance for \code{solve()}.
#' @param fastmom Logical. Use faster approximate gradient for \eqn{\lambda}? Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, the code reports messages and some values during optimization. 
#' @param print.init Logical. If \code{TRUE} the initial parameters used in the optimization of the first step are printed. 
#' @param ... additional arguments passed to \code{maxLik}.
#' @param x,object,  an object of class \code{bingmm}.
#' @param vce Character string. A string indicating what kind of standard errors should be computed when using \code{summary}. For the one-step GMM estimator, the options are \code{"robust"} and \code{"ml"}. For the two-step GMM estimator, the options are \code{"robust"}, \code{"efficient"} and \code{"ml"}. The option \code{"vce = ml"} is an exploratory method that evaluates the VC of the RIS estimator using the GMM estimates.
#' @param R Integer. Only valid if \code{vce = "ml"}. It indicates the number of draws used to compute the simulated probability in the RIS estimator.  
#' @param method Character string. Only valid if \code{vce = "ml"}. It indicates the algorithm used to compute the Hessian matrix of the RIS estimator. The default is \code{"bhhh"}.  
#' @param digits Integer. The number of digits for \code{summary} methods.  
#' @details 
#' 
#' The data generating process is:
#' 
#' \deqn{
#' y^*= X\beta + WX\gamma + \lambda W y^* + \epsilon = Z\delta + \lambda Wy^{*} + \epsilon
#' }
#' where  \eqn{y = 1} if \eqn{y^*>0} and 0 otherwise; \eqn{\epsilon \sim N(0, 1)} if \code{link = "probit"} or \eqn{\epsilon \sim L(0, \pi^2/3)} if \code{link = "logit"}. The general GMM
#'   estimator minimizes: 
#' \deqn{
#'  J(\theta) = g'(\theta)\hat{\Psi} g(\theta),
#' }
#' where \eqn{\theta = (\beta, \gamma, \lambda)} and 
#' \deqn{
#' g = n^{-1}H'v,
#' }
#' where \eqn{v} is the generalized residuals. Let \eqn{Z = (X, WX)}, then the instrument matrix \eqn{H} contains the linearly independent
#' columns of \eqn{H = (Z, WZ, ..., W^qZ)}. The one-step GMM estimator minimizes \eqn{J(\theta)} setting either 
#' \eqn{\hat{\Psi} = I_p} if \code{winitial = "identity"} or \eqn{\hat{\Psi} = (H'H/n)^{-1}} if \code{winitial = "optimal"}. The two-step GMM estimator
#' uses an additional step to achieve higher efficiency by computing the variance-covariance matrix of the moments \eqn{\hat{S}} to weight the sample moments.
#' This matrix is computed using the residuals or generalized residuals from the first-step, which are consistent. This matrix is computed as
#'  \eqn{\hat{S} = n^{-1}\sum_{i = 1}^n h_i(f^2/(F(1 - F)))h_i'} if \code{s.matrix = "robust"} or 
#'   \eqn{\hat{S} = n^{-1}\sum_{i = 1}^n \hat{v}_ih_ih_i'}, where \eqn{\hat{v}} are the first-step generalized residuals. 
#'   
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @return An object of class ``\code{bingmm}'', a list with elements:
#' \item{coefficients}{the estimated coefficients,}
#' \item{call}{the matched call,}
#' \item{callF}{the full matched call,}  
#' \item{X}{the X matrix, which contains also WX if the second part of the \code{formula} is used, }
#' \item{H}{the H matrix of instruments used,}
#' \item{y}{the dependent variable,}
#' \item{listw}{the spatial weight matrix,}
#' \item{link}{the string indicating the distribution of the error term,}
#' \item{Psi}{the moment-weighting matrix used in the last round,}
#' \item{type}{type of model that was fitted,}
#' \item{s.matrix}{the type of S matrix used in the second round,}
#' \item{winitial}{the moment-weighting matrix used for the first step procedure,}
#' \item{opt}{object of class \code{maxLik},}
#' \item{approximation}{a logical value indicating whether approximation was used to compute the inverse matrix,}
#' \item{pw}{the powers for the approximation,}
#' \item{formula}{the formula,}
#' \item{data}{the data, }
#' \item{contrastsX}{the contrasts used in the first part of the formula,}
#' \item{contrastsD}{the contrasts used in the second part of the formula,}
#' \item{Xlevels}{a record of the levels of the factors used in fitting,}
#' \item{fastmom}{whether the model was fitted using the faster gradient for \eqn{\lambda}.}
#' 
#' @examples
#' \donttest{
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # Two-step (Probit) GMM estimator
#' ts <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "twostep",
#'                 verbose = TRUE)
#'                 
#'# Robust standard errors
#'summary(ts)
#'# Efficient standard errors
#'summary(ts, vce = "efficient")
#'
#'# One-step (Probit) GMM estimator 
#'os <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "onestep",
#'                 verbose = TRUE)
#'summary(os)
#'
#'# One-step (Logit) GMM estimator with identity matrix as initial weight matrix
#'os_l <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                   link = "logit", 
#'                   listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                   data = COL.OLD, 
#'                   type = "onestep",
#'                   winitial = "identity", 
#'                   verbose = TRUE)
#'summary(os_l)
#'
#'# Two-step (Probit) GMM estimator with WX
#'ts_wx <- sbinaryGMM(CRIMED ~ INC + HOVAL| INC + HOVAL,
#'                    link = "probit", 
#'                    listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                    data = COL.OLD, 
#'                    type = "twostep",
#'                    verbose = FALSE)
#'summary(ts_wx)
#'
#'# Constrained two-step (Probit) GMM estimator 
#'ts_c <- sbinaryGMM(CRIMED ~ INC + HOVAL,
#'                   link = "probit", 
#'                   listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                   data = COL.OLD, 
#'                   type = "twostep",
#'                   verbose = TRUE, 
#'                   cons.opt = TRUE)
#'summary(ts_c)
#'}
#' @references 
#' 
#' Pinkse, J., & Slade, M. E. (1998). Contracting in space: An application of spatial statistics to discrete-choice models. Journal of Econometrics, 85(1), 125-154.
#' 
#' Fleming, M. M. (2004). Techniques for estimating spatially dependent discrete choice models. In Advances in spatial econometrics (pp. 145-168). Springer, Berlin, Heidelberg.
#' 
#' Klier, T., & McMillen, D. P. (2008). Clustering of auto supplier plants in the United States: generalized method of moments spatial logit for large samples. Journal of Business & Economic Statistics, 26(4), 460-471.
#' 
#' LeSage, J. P., Kelley Pace, R., Lam, N., Campanella, R., & Liu, X. (2011). New Orleans business recovery in the aftermath of Hurricane Katrina. Journal of the Royal Statistical Society: Series A (Statistics in Society), 174(4), 1007-1027.
#' 
#' Piras, G., & Sarrias, M. (2023). One or Two-Step? Evaluating GMM Efficiency for Spatial Binary Probit Models. Journal of choice modelling, 48, 100432. 
#' 
#' Piras, G,. & Sarrias, M. (2023). GMM Estimators for Binary Spatial Models in R. Journal of Statistical Software, 107(8), 1-33. 
#' @seealso \code{\link[spldv]{sbinaryLGMM}}, \code{\link[spldv]{impacts.bingmm}}.
#' @keywords models
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update)) 
#' @import stats methods Formula maxLik
#' @importFrom sphet listw2dgCMatrix
#' @export 
sbinaryGMM <- function(formula, 
                       data, 
                       listw = NULL, 
                       nins = 2,                                        # number of instruments
                       link = c("probit", "logit"),                     # probit or logit?
                       winitial = c("optimal", "identity"),             # initial moment-weighing matrix
                       s.matrix = c("robust", "iid"),                   # estimate for second round moment-weighing matrix
                       type = c("onestep", "twostep"),                  # one- and two-step procedure
                       gradient = TRUE,                                 # use the gradient of J for the minimization?
                       start    = NULL,                                 # vector of starting values 
                       cons.opt = FALSE,                                # use constrained optimization for lambda?
                       approximation = FALSE,                           # use inverse approximation?
                       verbose = TRUE,                                  # print messages?
                       print.init = FALSE,                              # print initial values?
                       pw = 5,                                          # default powers for the approximation 
                       tol.solve = .Machine$double.eps,                 # tolerance for solve
                       fastmom   = FALSE,                               # faster gradient
                       ...){
  # winitial: Weight matrix of moment for first step. 
  #      1. If optimal, then  Psi=((1/n)Z'Z)^{-1} Eq(13)
  #      2. If identity, then Psi = I Eq(16)
  # s.matrix is the estimated variance matrix for the moments
  #     1. If iid then S = n^1sum_i u_i^2h_ih_i'
  #     2. If robust then as in the paper
  # In this version, we use the notation in paper. 
  
  # --- Argument validation ---
  winitial    <- match.arg(winitial)
  type        <- match.arg(type)
  s.matrix    <- match.arg(s.matrix)
  link        <- match.arg(link)
  
  # --- Spatial weights matrix ---
  if (inherits(listw, "listw")) {
    W <- sphet::listw2dgCMatrix(listw)
  } else if (inherits(listw, "matrix")) {
    W <- Matrix::Matrix(listw)
  } else if (inherits(listw, "Matrix")) {
    W <- listw
  } else {
    stop("Neighbourhood list or listw format unknown")
  }

  
  # --- Model frame setup ---
  callT      <- match.call(expand.dots = TRUE)
  callF      <- match.call(expand.dots = FALSE)
  mf         <- callT
  m          <- match(c("formula", "data"), names(mf), 0L)
  mf         <- mf[c(1L, m)]
  f1         <- Formula(formula)
  Durbin     <- (length(f1)[2L] == 2L)
  mf$formula <- f1 
  mf[[1L]]   <- as.name("model.frame")
  mf         <- eval(mf, parent.frame())
  nframe     <- length(sys.calls())
  
  # Set optimization defaults
  if (is.null(callT$method)) callT$method   <- 'bfgs'
  if (is.null(callT$iterlim)) callT$iterlim <- 10000
  if (is.null(callT$reltol)) callT$reltol   <- 1e-6
  callT$finalHessian <- FALSE # We do not require the Hessian. This speeds the optimization procedure. 
  
  # Extract and validate response variable
  y  <- model.response(mf)
  if (anyNA(y)) stop("NAs in dependent variable")
  if (!all(y %in% c(0, 1, TRUE, FALSE))) stop("All dependent variables must be either 0, 1, TRUE or FALSE")
  if (!is.numeric(y)) y <- as.numeric(y)
  
  # Construct design matrix
  X          <- model.matrix(f1, data = mf, rhs = 1)
  contrastsX <- attr(X, "contrasts")
  Xlevels    <- .getXlevels(attr(mf, "terms"), mf)
  
  # Hadle Durbin model case
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
  if (N != nrow(W)) stop("Number of spatial units in W is different to the number of data")
  
  ## Generate initial instruments 
  Z <- make.instruments(W, x = X, q = nins)
  H <- cbind(X, Z)
  
  # Drop collinear instruments using QR
  qr_H <- qr(H)
  H <- H[, qr_H$pivot[seq_len(qr_H$rank)]]
  P <- ncol(H)
  if (P < K) stop("Underspecified model (insufficient instruments)")
  if (anyNA(H)) stop("NAs in the instruments")
  
  ## Starting values for optimization of GMM
  if (is.null(start)){
    # Try initial values from probit/logit
    glm_attempt <- tryCatch(
      glm.fit(as.matrix(X), y, family = binomial(link = link)),
      warning = function(w) {
        warning("Warning during glm.fit(): ", conditionMessage(w))
        return(NULL)
      },
      error = function(e) {
        warning("Error during glm.fit(): ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (!is.null(glm_attempt) && !anyNA(glm_attempt$coefficients)) {
      b_init <- glm_attempt$coefficients
    } else {
      warning("Falling back to constant starting values for beta.")
      b_init <- rep(0.1, K)
    }
    
    # Handle potential NA or unstable lambda.init
    Wy <- as.numeric(crossprod(t(W), y))
    lambda.init <- suppressWarnings(cor(y, Wy))
    if (is.na(lambda.init)) {
      warning("Initial lambda correlation is NA; setting to 0.1.")
      lambda.init <- 0.1
    }
    
    theta <- c(b_init, lambda.init)
    names(theta) <- c(colnames(X), "lambda")
  } else {
    theta <- start
    if (length(start) != (K + 1)) stop("Incorrect number of intial parameters")
    names(theta) <- c(colnames(X), "lambda")
  }
  if (print.init) {
    cat("\nStarting Values:\n")
    print(theta)
  } 
  
  # --- Constraints for lambda ---
  if (cons.opt){
    sym          <- all(W == t(W))
    omega        <- eigen(W, only.values = TRUE, symmetric = sym)
    lambda_space <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
    # A %*% theta + B >= 0 
    A <- rbind(c(rep(0, K), 1), 
               c(rep(0, K), -1))  # Matrix of restrictions. See help(maxLik)
    B <- cbind(c(-1* (lambda_space[1] + .Machine$double.eps), lambda_space[2] - .Machine$double.eps))
    callT$constraints <- list(ineqA = A, ineqB = B)
  }
  
  # --- Intial Weight matrix Psi ---
  Psi <- if (winitial == "optimal") {
    Matrix::solve(crossprod(H)/N, tol = tol.solve)
    #chol2inv(chol(Matrix::crossprod(H)/N))  # More stable than solve() for symmetric matrix
  } else {
    Matrix::Diagonal(ncol(H))
  }
  
  # --- First-step GMM estimation ---
  if (verbose) cat("\nFirst-step GMM optimization based on", winitial, "initial weight matrix \n")
  opt <- callT
  opt$start <- theta
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHessian', 'reltol'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]]      <- as.name('maxLik')
  opt$logLik    <- as.name('J_minML')
  opt$gradient  <- as.name('gradient')
  opt$link      <- as.name('link')
  opt$tol.solve <- as.name('tol.solve')
  opt$listw     <- as.name('W')
  opt$approximation    <- as.name('approximation')
  opt$pw <- as.name('pw')
  opt[c('y', 'X', 'H', 'Psi')] <- list(as.name('y'), 
                                       as.name('X'), 
                                       as.name('H'), 
                                       as.name('Psi'))
  opt$fastmom     <- as.name('fastmom')
  x     <- eval(opt, sys.frame(which = nframe))
  b_hat <- coef(x)

  # --- Second-step GMM estimation ---
  if (type == "twostep") {
    # Make S matrix
    S   <- makeS(b_hat         = b_hat, 
                 y             = y, 
                 X             = X, 
                 H             = H, 
                 listw         = W, 
                 link          = link, 
                 wmatrix       = s.matrix, 
                 approximation = approximation, 
                 pw            = pw, 
                 fastmom       = fastmom)
    Psi <- Matrix::solve(S, tol = tol.solve)
    #Psi <- chol2inv(chol(S))  # More stable than solve() for symmetric matrix
    
    if (verbose) cat("\nSecond-step GMM optimization using S moment-weighing matrix\n")
    opt$start     <- b_hat
    opt[c('Psi')] <- list(as.name('Psi'))
    x             <- eval(opt, sys.frame(which = nframe))
    b_hat         <- coef(x)
  }
  
  # --- Output ---
  out <- structure(
    list(
      coefficients  = b_hat, 
      call          = callT,
      callF         = callF,
      X             = X, 
      H             = H, 
      y             = y, 
      listw         = W,
      link          = link, 
      Psi           = Psi, 
      type          = type, 
      s.matrix      = s.matrix, 
      winitial      = winitial, 
      opt           = x, 
      approximation = approximation,
      pw            = pw, 
      formula       = f1, 
      tol.solve     = tol.solve,
      mf            = mf,
      data          = data,
      contrastsX    = contrastsX,
      contrastsD    = contrastsD,
      Xlevels       = Xlevels,
      fastmom       = fastmom
    ), 
    class = "bingmm"
  )
  out
}

##### Other functions ----
# Approximate (I-\rho W)^-1
app_W <- function(listw, lambda, pw){
  n           <- dim(listw)[1]
  wi          <- listw
  lambda_loop <- lambda 
  Dn          <- Matrix::Diagonal(n) 
  pW          <- Matrix::Matrix(0, n, n)
  i <- 1
  while (i < pw) {
    wi          <- listw %*% wi 
    lambda_loop <- lambda * lambda_loop
    pW          <- pW + (lambda_loop * wi)
    i           <- i + 1
  }
  pW <- Dn + lambda*listw +  pW
}

# Moment function
momB_slm <- function(start, y, X, H, listw, link, approximation, pw, tol.solve = .Machine$double.eps){
  # This function generates: 
  #  (1) generalized residuals
  #  (2) the moment conditions 
  #  (3) the gradient dv/dtheta for SLM binary models
  K      <- ncol(X)
  N      <- nrow(X)
  beta   <- start[1:K]
  lambda <- start[K + 1]
  
  W <- listw
  A <- Matrix::Diagonal(N) - lambda * W
  B <- if (approximation) app_W(W, lambda, pw) else Matrix::solve(A, tol = tol.solve)
  
  # Avoid full Sigma_u if only diag is needed
  sigma     <- sqrt(Matrix::rowSums(B^2))  # More efficient than  diag(tcrossprod(B))
  sigma_inv <- 1 / sigma
  BX        <- B %*% X
  BXbeta    <- BX %*% beta
  ai        <- as.vector(BXbeta) * sigma_inv
  
  # Choose link functions efficiently
  if (link == "probit") {
    pfun <- pnorm
    dfun <- dnorm
    ddfun <- function(x) -x * dnorm(x)
  } else if (link == "logit") {
    pfun <- plogis
    dfun <- dlogis
    ddfun <- function(x) (1 - 2 * plogis(x)) * plogis(x) * (1 - plogis(x))
  }
  
  # Residuals calculation
  q    <- 2 * y - 1
  ai_q <- q * ai
  fa   <- pmax(dfun(ai_q), .Machine$double.eps)
  Fa   <- pmax(pfun(ai_q), .Machine$double.eps)
  ffa  <- ddfun(ai_q)
  
  # Generalized residuals and moment conditions
  v    <- q * (fa / Fa)                  
  g    <- Matrix::crossprod(H, v) / N
  
  # Make gradient
  der <- (ffa * Fa - fa^2) / (Fa^2)
  Gb  <- BX * der * sigma_inv 
  
  # NOTE: previous version of D_lambda incorrectly assumes that W is always symmetric
  part1  <- (B %*% (W %*% BXbeta)) * sigma_inv
  
  # Compute diag(A^{-1}W Sigma_u + Sigma_u W' (A^{-1})')
  # Let A = BW(BB') so M = A + A', and diag(M) = diag(A) + diag(A') = 2 * diag(A)
  # For diag(A) = rowSums(BW * BBt)
  #diag_term <- 2 * Matrix::diag(B %*% W %*% Matrix::tcrossprod(B))
  
  BW        <- B %*% W
  BBt       <- Matrix::tcrossprod(B)
  diag_term <- 2 * Matrix::rowSums(BW * BBt)
  Drho      <- - 0.5 * sigma_inv^3 * diag_term
  
  part2  <- Drho * BXbeta
  G_rho  <- der * (part1 + part2)
  
  G <- cbind(Gb, G_rho)
  out <- list(
    g = g,
    G = G,
    v = v,
    u = y - pfun(ai),
    vu = (fa^2) / (Fa * (1 - Fa))
  )
  return(out)
}

# Faster moment function
FmomB_slm <- function(start, y, X, H, listw, link, approximation, pw, tol.solve = .Machine$double.eps) {
  # This function generates: 
  #  (1) generalized residuals
  #  (2) the moment conditions 
  #  (3) the gradient dv/dtheta for SLM binary models
  K      <- ncol(X)
  N      <- nrow(X)
  beta   <- start[1:K]
  lambda <- start[K + 1]
  
  W <- listw
  A <- Matrix::Diagonal(N) - lambda * W
  B <- if (approximation) app_W(W, lambda, pw) else Matrix::solve(A, tol = tol.solve)
  
  # Avoid full Sigma_u if only diag is needed
  sigma     <- sqrt(Matrix::rowSums(B^2))  # More efficient than  diag(tcrossprod(B))
  sigma_inv <- 1 / sigma
  BX        <- B %*% X
  BXbeta    <- BX %*% beta
  ai        <- as.vector(BXbeta) * sigma_inv
  
  # Choose link functions efficiently
  if (link == "probit") {
    pfun <- pnorm
    dfun <- dnorm
    ddfun <- function(x) -x * dnorm(x)
  } else if (link == "logit") {
    pfun <- plogis
    dfun <- dlogis
    ddfun <- function(x) (1 - 2 * plogis(x)) * plogis(x) * (1 - plogis(x))
  }
  
  # Residuals calculation
  q    <- 2 * y - 1
  ai_q <- q * ai
  fa   <- pmax(dfun(ai_q), .Machine$double.eps)
  Fa   <- pmax(pfun(ai_q), .Machine$double.eps)
  ffa  <- ddfun(ai_q)
  
  # Generalized residuals and moment conditions
  v    <- q * (fa / Fa)                  
  g    <- Matrix::crossprod(H, v) / N
  
  # Make gradient
  der <- (ffa * Fa - fa^2) / (Fa^2)
  Gb  <- BX * der * sigma_inv 
  
  part1  <- (B %*% (W %*% BXbeta)) * sigma_inv
  
  # Approximation

  WB        <- W %*% B
  BWB       <- B %*% WB
  diag_term <- 2 * Matrix::rowSums(B * BWB)
  Drho      <- - 0.5 * sigma_inv^3 * diag_term
  
  part2  <- Drho * BXbeta
  G_rho  <- der * (part1 + part2)
  
  G <- cbind(Gb, G_rho)
  out <- list(
    g = g,
    G = G,
    v = v,
    u = y - pfun(ai),
    vu = (fa^2) / (Fa * (1 - Fa))
  )
  return(out)
}

# Objective function to be minimized
J_minML <- function(start, y, X, H, listw, link, Psi, gradient = FALSE, approximation = FALSE, pw = NULL, tol.solve = 1e-8, fastmom = FALSE) {
  # Retrieve moment conditions and derivatives
  getR <- if (fastmom) {
    FmomB_slm(start, y, X, H, listw, link, approximation, pw, tol.solve)
  } else {
    momB_slm(start, y, X, H, listw, link, approximation, pw, tol.solve)
  }
  
  g <- getR$g
  Psi_t <- t(Psi)  # pre-transpose once
  gP <- Matrix::crossprod(Psi_t, g)  # Psi' * g
  J <- -1 * Matrix::crossprod(g, gP)  # g' * Psi * g
  
  if (gradient) {
    G  <- getR$G
    GdivN <- G / nrow(X)  # scale G once
    HG <- Matrix::crossprod(H, GdivN)  # H' * (G / N)
    Gr <- -2 * crossprod(HG, gP)
    attr(J, "gradient") <- as.vector(Gr)
  }
  
  return(J)
}
# #Make S matrix: var-cov of moments
# makeS <- function(b_hat, y, X, H, listw, link, wmatrix, approximation, pw, fastmom){
#   N <- nrow(X)
#   evm <- if (fastmom) {
#     FmomB_slm(start = b_hat, y = y, X = X, H = H, listw = listw, link = link, approximation = approximation, pw = pw)
#   } else {
#     momB_slm(start = b_hat, y = y, X = X, H = X, listw = listw, link = link, approximation = approximation, pw = pw)
#   }
#   if (wmatrix == "iid") {
#     u_hat <- evm$v # predicted generalized residuals
#     Shat <- 0
#     for (i in 1:N) {
#       Shat <- Shat + (u_hat[i] ^ 2 * tcrossprod(H[i,])) #klier and McMillen (2008)
#     }
#     Shat <- Shat / N
#   }
#   if (wmatrix == "robust"){
#     vu       <- evm$vu # This is f^2 / ((1 - F)*F)
#     Shat <- 0
#     for (i in 1:N) {
#       Shat <- Shat + (H[i, ] %*% t(H[i, ] * vu[i]))    ## Equation(15)
#     }
#     #Shat <- Shat / (N^2)
#     Shat <- Shat / N
#   }
#   return(Shat)
# }

# Make S matrix: var-cov of moments
makeS <- function(b_hat, y, X, H, listw, link, wmatrix, approximation, pw, fastmom) {
  N <- nrow(X)

  evm <- if (fastmom) {
    FmomB_slm(start = b_hat, y = y, X = X, H = H, listw = listw, link = link, approximation = approximation, pw = pw)
  } else {
    momB_slm(start = b_hat, y = y, X = X, H = X, listw = listw, link = link, approximation = approximation, pw = pw)
  }

  if (wmatrix == "iid") {
    u_hat <- evm$v  # residuals
    # Vectorized version of: sum_i (u_hat[i]^2 * H[i,] %*% t(H[i,]))
    weights    <- u_hat^2
    H_weighted <- H * weights  # broadcasting weights over rows
    Shat       <- Matrix::crossprod(H, H_weighted) / N

  } else if (wmatrix == "robust") {
    vu         <- evm$vu  # variance scaling factor
    H_weighted <- H * vu  # broadcasting vu over rows
    Shat       <- Matrix::crossprod(H, H_weighted) / N
  } else {
    stop("Unsupported weighting matrix option")
  }
  return(Shat)
}

### Overtest
#over.test <- function(object, ...){
#  H <- object$H
#  K <- length(coef(object))
#  P <- ncol(H)
#  N <- nrow(H)
#  J <- -N * object$opt$maximum 
#  out <- list(statistic =  J, df = P - K, p.value =  pchisq(J, P - K, lower.tail =  FALSE))
#  return(out)
#}


#' @title Make instruments for spatial models
#' @name make.instruments
#' @usage make.instruments(listw, x, q)
#' @param listw object. An object of class \code{matrix}, or \code{Matrix}.  
#' @param x variable(s) to be lagged
#' @param q number of lags
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @keywords instruments
#' @export 
make.instruments <- function(listw, x, q = 3){
  # This function creates the instruments (WX, ...,W^qX) as in K&P
  W <- listw
  # Drop constant (if any)
  names.x <- colnames(x)
  if (names.x[1] == "(Intercept)") x <- matrix(x[,-1], dim(x)[1], dim(x)[2] - 1)
  names.x <- names.x[which(names.x != "(Intercept)")]
  sq1 <- seq(1, ncol(x) * q, ncol(x))
  sq2 <- seq(ncol(x), ncol(x) * q, ncol(x))
  Zmat <- matrix(NA, nrow = nrow(x), ncol = ncol(x) * q)
  names.ins <- c()
  for (i in 1:q) {
    Zmat[, sq1[i]:sq2[i]] <- as.matrix(W %*% x)
    x <- Zmat[, sq1[i]:sq2[i]]
    names.ins <- c(names.ins, paste(paste(replicate(i, "W"), collapse = ""), names.x, sep = "*"))
  }
  colnames(Zmat) <- names.ins
  return(Zmat)
}


### S3 methods ----

#' @rdname sbinaryGMM
#' @export
coef.bingmm <- function(object, ...){
  object$coefficients
}


#' @rdname sbinaryGMM
#' @method vcov bingmm
#' @export 
vcov.bingmm <- function(object, vce = c("robust", "efficient", "ml"), method = "bhhh", R = 1000, tol.solve = .Machine$double.eps, ...){
  # vce: indicates the estimator for the variance-covariance matrix when using two-step estimator
  #       1- if efficient, then the lowest bound of the vc is used. 
  vce           <- match.arg(vce)
  R             <- R
  method        <- method
  link          <- object$link
  listw         <- object$listw
  type          <- object$type 
  theta         <- object$coefficients
  s.matrix      <- object$s.matrix
  winitial      <- object$winitial
  approximation <- object$approximation
  pw            <- object$pw
  K             <- length(theta)
  y             <- object$y
  H             <- object$H
  X             <- object$X
  Psi           <- object$Psi
  N             <- nrow(X)
  fastmom       <- object$fastmom
  # FIXME:  Si traes el gradiente desde el optimizer no necessitas calcularlo otra vez aca...
  # Compute G only once at the optimum using exact gradient
  G      <- momB_slm(start = theta, y = y, X = X, H = H, listw = listw, link = link, approximation = approximation, pw = pw)$G
  G_bar  <- Matrix::crossprod(H, G)
  
  if (vce == "efficient"){
    if (type == "onestep") stop("Efficient VCE for one-step estimator not yet implemented")
    # Similar to Stata, we use the the Psi that was used to compute the final-round estimate 
    A <- Matrix::crossprod(G_bar, Psi %*% G_bar)
    V <- N * solve(A, tol = tol.solve)
  }
  if (vce == "robust"){
    # Based on Stata: Psi is the weight matrix requested with wmatrix() and it is calculated based on the residuals obtained after the first estimation step. 
    S       <- makeS(b_hat = theta, y = y, X = X, H = H, listw = listw, link = link, wmatrix = s.matrix, approximation = approximation, pw = pw, fastmom = fastmom)
    GPG     <- Matrix::crossprod(G_bar, Psi %*% G_bar)            # G'Psi G
    GPSPG   <- Matrix::crossprod(G_bar, Psi%*% S %*% (Psi %*% G_bar)) # G'Psi S Psi G
    inv_GPG <- Matrix::solve(GPG, tol = tol.solve)
    V       <- N * inv_GPG  %*% GPSPG %*% inv_GPG 
  }
  if (vce == "ml"){
    opt <- object$call
    opt$start <- theta
    m <- match(c('start'),
               names(opt), 0L)
    opt$iterlim <- 0
    opt <- opt[c(1L, m)]
    opt[[1]]     <- as.name('maxLik')
    opt$logLik   <- as.name('ll_ris')
    W <- listw
    opt$method <- as.name("method")
    opt$W    <- as.name('W')
    opt$R    <- as.name('R')
    opt[c('y', 'X')] <- list(as.name('y'), as.name('X'))
    x <- eval(opt)
    return(vcov(x))
  }
  colnames(V) <- rownames(V) <- names(theta)
  return(V)
}


#' @rdname sbinaryGMM
#' @method print bingmm
#' @export 
print.bingmm <- function(x, 
                         digits = max(3, getOption("digits") - 3),
                         ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print.default(format(drop(coef(x)), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}


#' @rdname sbinaryGMM
#' @method summary bingmm
#' @export
summary.bingmm <- function(object, vce = c("robust", "efficient", "ml"), method = "bhhh", R = 1000, tol.solve = .Machine$double.eps, ...){
  vce                 <- match.arg(vce)
  b                   <- object$coefficients
  fastmom             <- object$fastmom
  std.err             <- sqrt(diag(vcov(object, vce = vce, method = method, R = R, tol.solve = tol.solve)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.bingmm", "bingmm")
  return(object)
}


#' @rdname sbinaryGMM
#' @method print summary.bingmm
#' @export
print.summary.bingmm <- function(x,
                                 digits = max(5, getOption("digits") - 3),
                                 ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      SLM Binary Model by GMM \n")
  cat("        ------------------------------------------------------------\n")
  
  
  cat("\nCall:\n")
  cat(paste(deparse(x$callF), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  
  cat(paste("\nSample size:", signif(nrow(x$X), digits)), "\n")
  #cat(paste("\nJ at the optimum:", signif(-1 * x$opt$maximum, digits)))
  #cat(paste("\nInstruments:", paste(colnames(x$H))))
  invisible(x)
}

#' Get Model Summaries for use with "mtable" for objects of class bingmm
#' 
#' A generic function to collect coefficients and summary statistics from a \code{bingmm} object. It is used in \code{mtable}
#' 
#' @param obj a \code{bingmm} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @return A list with an array with coefficient estimates and a vector containing the model summary statistics. 
#' @importFrom memisc getSummary
#' @method getSummary bingmm
#' @export 
getSummary.bingmm <- function(obj, alpha = 0.05, ...){
  if (inherits(obj, c("summary.bingmm"))){
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
  #sumstat <- c(logLik = NA, deviance = NA, AIC = NA, BIC = NA, N = N, 
  #             LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
  #             Nagelkerke = NA)
  sumstat <- c(N = N)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

#' Predictions for Spatial Binary GMM Models
#'
#' Computes predicted probabilities for spatial binary response models estimated via GMM. 
#' Supports both probit and logit links, accounts for spatial heteroskedasticity, and optionally 
#' returns standard errors using the Delta method with robust or efficient variance estimators.
#'
#' @param object An object of class \code{bingmm} returned by a spatial binary GMM estimation function.
#' @param newdata An optional data frame in which to look for variables with which to predict. 
#' If omitted, the original data used to fit the model is used.
#' @param Sinv Optional user-supplied spatial multiplier matrix \eqn{S = (I - \lambda W)^{-1}}. 
#' If \code{NULL}, it is computed using the spatial weight matrix.
#' @param het Logical. If \code{TRUE}, assumes a heteroskedastic error structure with spatially varying variances.
#' @param approximation Logical. If \code{TRUE}, uses power-series approximation to compute the inverse spatial matrix.
#' @param pw Integer. Power-order to use when \code{approximation = TRUE}.
#' @param ses Logical. If \code{TRUE}, standard errors of the predictions are computed using the Delta method.
#' @param vce Type of variance-covariance estimator: \code{"robust"} (default), \code{"efficient"}, or \code{"ml"}.
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
#' @seealso \code{\link{sbinaryGMM}}, \code{\link{vcov.bingmm}}
#'
#' @examples
#' # Data set
#' data(oldcol, package = "spdep")
#' 
#' # Create dependent (dummy) variable
#' COL.OLD$CRIMED <- as.numeric(COL.OLD$CRIME > 35)
#' 
#' # Two-step (Probit) GMM estimator
#' ts <- sbinaryGMM(CRIMED ~ INC*HOVAL + factor(CP),
#'                 link = "probit", 
#'                 listw = spdep::nb2listw(COL.nb, style = "W"), 
#'                 data = COL.OLD, 
#'                 type = "twostep",
#'                 verbose = FALSE)
#'                 
#' # Compute just the predicted probability
#' p_hat <- predict(ts)
#' 
#' # Compute predicted probability and standard errors
#' res <- predict(ts, ses = TRUE, het = TRUE, vce = "efficient")
#' head(res)
#' @author Mauricio Sarrias and Gianfranco Piras. 
#' @keywords prediction
#' @export 
#' @method predict bingmm
predict.bingmm <- function(object, 
                           newdata, 
                           Sinv = NULL,
                           het  = TRUE,
                           approximation = FALSE,
                           pw  = 5,
                           ses = FALSE,
                           vce = c("robust", "efficient", "ml"),
                           theta = NULL, ...){
  
  if (!inherits(object, "bingmm")) warning("calling predict.bingmm(<fake-bingmm-object>) ...")
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
  vce <- match.arg(vce)
  V   <- vcov(object, vce = vce)
  se  <- sqrt(rowSums((Jac %*% V) * Jac)) # Efficient diag(JVJ')
  
  # Return full prediction table
  z    <- pred / se
  pval <- 2 * pnorm(-abs(z))
  
  return(cbind(`p_hat` = pred, `Std. error` = se, `z value` = z, `Pr(> z)` = pval))
}

