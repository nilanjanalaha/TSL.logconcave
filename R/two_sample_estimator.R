# Function for calculaing delta_hat by pooling samples and smoothed density
# Function for calculating the location estimator for two sample data
#' Logconcave estimator of the location-shift in the two sample model.
#'
#' Suppose m univariate observations \eqn{X_1, \ldots, X_m} are sampled from a density \eqn{g(x-\mu)},
#' and n univariate observations \eqn{Y_1, \ldots, Y_n} are sampled from the density \eqn{g(x-\mu-\Delta)}, where
#' \eqn{g} is an unknown log-concave density. This function computes
#' a one step estimator to estimte \eqn{\Delta}. This estimator relies on
#' the smoothed log-concave MLE estimator from the package \code{\link{logcondens}} to estimate \eqn{g}, and  is root-n consistent
#' for \eqn{\Delta} provided \eqn{g} is log-concave. 
#' 
#' 
#' @param dat       A list with two components: x and y, each being vector of possibly different lengths; represents the data.
#' @param eta A fraction between 0 and 1/2. Corresponds to the truncation level of the one step estimator.
#'              The default is 0.0001.
#' 
#' @details \code{eta:} If eta is zero, the function computes the one step estimator  without any
#'                      truncation. See Saha et al. (2023) for more details.
#'                        
#'@return   A vector of length two.   
#'\itemize{
#'\item \code{estimate:}  The estimated value of \eqn{\Delta}.
#'\item \code{FI:} The estimated Fisher information for estimating \eqn{\Delta}.
#'}
#'  
#'@references Saha R., Dey P., Laha N. (2023). \emph{Revisiting the two-sample location model 
#' with log-concavity assumption}. submitted.
#'@author \href{https://www.nilanjana91.de/}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@tamu.edu}.\cr
#'Ridhiman Saha, \email{riddhimansaha@@fas.harvard.edu}
#' 
#' @examples
#' x <- rlogis(100); y <- rlogis(150) + 0.1; 
#' pooled_smoothed(list(x=x, y=y), eta = 0.0001)
#' @export
pooled_smoothed <- function(dat, eta = 0.0001) {
  x <- dat$x; y <- dat$y
  n_1 <- length(x)
  n_2 <- length(y)
  n <- n_1 + n_2
  # Center and merge data
  mu_bar <- mean(x)
  delta_bar <- mean(y) - mean(x)
  
  # center data and estimate densities
  x_centered <- x - mu_bar
  y_centered <- y - mu_bar - delta_bar
  
  # Pool samples
  z <- c(x - mu_bar, y - mu_bar - delta_bar)
  z <- sort(z, index.return = TRUE)
  # Keep track of source in the merged vector
  x_indices <- which(z$ix <= n_1)
  y_indices <- which(z$ix > n_1)
  z <- z$x # This is our pooled sample, sorted.
  
  # Fit log-concave smoothed density on z
  res <- logcondens::logConDens(z, smoothed = TRUE)
  knots <- res$knots
  slopes <- diff(res$phi[res$IsKnot == 1]) / diff(knots)
  intercepts <- res$phi[res$IsKnot == 1][-1] - slopes * knots[-1]
  VarFn <- logcondens::LocalVariance(x = res$x, w = res$w, phi = res$phi)
  # std. dev. of gaussian kernel so that var of est smoothed density
  # equals sample s.d.
  gam <- sqrt(res$sig^2 - VarFn)
  
  # Function to find (CDF - eta) based on fitted density
  findCDF <- function(x, res, eta=0) {
      -eta + as.numeric(logcondens::evaluateLogConDens(x, res, which=5)[, "smooth.CDF"])
  }
  
  # Function for finding quantiles
  findQuantile <- function(eta, res) {
    # if eta = 0 or eta = 1, return some extreme value
    if(eta == 0) {
      return(min(res$xs))
      # return(2 * min(res$xs) - max(res$xs))
    }
    if(eta == 1) {
      return(max(res$xs))
      # return(2 * max(res$xs) - min(res$xs))
    }
    # For other values, find by uniroot() function
    # eta^th quantile. res is "smoothed" estimator
    if(is.unsorted(res$F.smoothed)) {
      res$F.smoothed <- sort(res$F.smoothed)
    }
    ind <- findInterval(eta, res$F.smoothed, rightmost.closed=TRUE)
    if(ind == 0) {
      # Then CDF at res$xs[1] is bigger then eta
      j <- 0
      upper <- res$xs[1]
      while(TRUE) {
        j <- 2^j
        if(findCDF(upper - j, res, eta=eta) < 0) {
          lower <- upper - j
          xi_1 <- uniroot(findCDF, lower=lower, upper=upper,
            res=res, eta=eta)$root
          break
        }
      }
    } else if (ind == length(res$xs)) {
      # Then CDF at res$xs[500] is smaller then eta
      # print("ewjkhg")
      j <- 0
      lower <- res$xs[length(res$xs)]
      while(TRUE) {
        j <- 2^j
        if(findCDF(lower + j, res, eta=eta) > 0) {
          upper <- lower + j
          xi_1 <- uniroot(findCDF, lower=lower, upper=upper,
            res=res, eta=eta)$root
          break
        }
      }
    } else if(abs(res$F.smoothed[ind] - eta) < 1e-6) {
        xi_1 <- res$xs[ind]
    } else if(abs(res$F.smoothed[ind+1] - eta) < 1e-6) {
        xi_1 <- res$xs[ind+1]
    } else {
        xi_1 <- uniroot(findCDF, lower=res$xs[ind], upper=res$xs[ind+1],
            res=res, eta=eta)$root
    }
    xi_1
  }
  
  xi_1 <- findQuantile(eta, res)
  xi_2 <- findQuantile(1-eta, res)
  
  # Evaluate \psi' at sample points
  # First find g'
  g_slope <- function(x, knots, slopes, intercepts, gam) {
    # x: vector on which dg/dx to be calculated
    # knots: vector of knots K_1, ..., k_L from logConDens
    # slopes: b_1, ..., b_{L-1} : slope of linear pieces b/w knots
    # intercepts: a_1, ... , a_{L-1]}
    # gam: Variance used in logConDens smoothing (or gaussian kernel variance)
    n <- length(x)
    L <- length(knots)
    part_2 <- rep(0, n)
    for(l in 1:(L-1)) {
      part_2 <- part_2 + slopes[l] *
        exp(intercepts[l] + slopes[l]*x + 0.5*slopes[l]^2 *gam^2) * (
          pnorm(knots[l+1], mean = x + slopes[l]*gam^2, sd = gam) -
            pnorm(knots[l], mean = x + slopes[l]*gam^2, sd = gam)
        )
    }
    return(part_2)
  }
  psi_slopes <- g_slope(z, knots, slopes, intercepts, gam) /
    logcondens::evaluateLogConDens(z, res, which=4)[, "smooth.density"]
  
  # Find truncated information
  integrand <- function(x, res, gam) {
    slopes <- diff(res$phi[res$IsKnot == 1]) / diff(knots)
    intercepts <- res$phi[res$IsKnot == 1][-1] - slopes * knots[-1]
    g_slopes <- g_slope(x, res$knots, slopes, intercepts, gam)
    g <- logcondens::evaluateLogConDens(x, res, which=4)[, "smooth.density"]
    g_slopes^2 / g
  }
  information <- integrate(integrand, lower=xi_1, upper=xi_2,
                           res=res, gam=gam)$value
  
  # One-step estimator
  include <- (z >= xi_1) & (z <= xi_2)
  # Now this is the integral w.r.t. empirical cdf
  delta_hat <- delta_bar +
    (mean(psi_slopes[x_indices] * include[x_indices]) / information) -
    (mean(psi_slopes[y_indices] * include[y_indices]) / information)
  
  return(c(delta_hat, information))
}
