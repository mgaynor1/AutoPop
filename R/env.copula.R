#' Environmental Copula -  Calculates correlation probability.
#'
#' @description Calculates the probability coefficient which allow cytotype and stage correlations
#'
#'
#' @param rho Correlation coefficient. Must be a single integer between 0 and 1.
#' @inheritParams gen.iter.f.choosy
#'
#' @returns A list of two matrix: X and U. Here, X is the matrix which was simulated from
#' a multivariate normal distribution. The U matrix is simply the X matrix transformed into a uniform distribution.
#' Each matrix has 3 columns, one for each cytotype. The length of the matrix matches the number of generations provided.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm
#'


env.copula <- function(rho, generations){
  # Create covariance matrix
  n <- 3 # Three cytotypes
  mean <- rep(0,n)
  # covariance matrix for the cytotypes
  SIGMA <- ((1-rho)*diag(n) + rho*matrix(1, nrow = n, ncol = n))
  # Simulate from a multivariate normal distribution
  X <- MASS::mvrnorm(n = generations, mu = mean, Sigma = SIGMA)
  # Transform to uniform using pnorm, which is the cdf for a normal distribution
  U <- pnorm(X)
  out <- list(X, U)
  return(out)
}


