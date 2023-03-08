#' MAR(1) Conditional least squares parameter estimates.
#'
#' @description  As described in  Ives et al. 2003, the Multivariate Autoregressive Model, also known as the MAR(1) model,
#'    is a discrete-time model for a multispecies stochastic community subject to environmental noise.  It is a Markov process.
#'    Given multispecies time series (log) abundance data without observation error, parameter estimation of this model parameters
#'    can be quickly done via Conditional Least Squares.
#'
#' @param comm.mat Data frame with log-scale abundance for each time step following burn-in.
#'   We expect one column per species and one row per time-step.
#' @param covariate Matrix of covariates (ex. precipitation per time step) with each column representing a variable and
#'   each row representing a time step. Defaults to "NULL" when no matrix is supplied.
#'
#' @returns Resulting data frame includes:
#'
#' * \strong{A}, vector of length p (p = number of species) containing the estimates of the intrinsic rate of natural increase of each species.
#' * \strong{B},  a matrix of p x p estimates where the diagonal elements represent the intra-specific, density-dependent effects.
#'      The elements \eqn{b_{ij}} gives the effect of the abundance of species j on per capita growth rate of species i.
#' * \strong{sigma}, a matrix of p x p, represents environmental noise variance-covariance matrix.
#' * \strong{C}, a vector of estimated coefficients for every covariate representing the effects of every covariate on every species.
#' * \strong{E}, vector representing stochastic environmental variability.
#' * \strong{Yhat}, \eqn{\hat{Y}} estimator of predicted abundance.
#' * \strong{R2}, \eqn{R^2}, proportion of explained variation in the log scale population abundance for each species.
#' * \strong{R2_D}, conditional \eqn{R^2}, proportion of variation in the change in one unit of time of the log scale population
#'      abundance explained by the model for each species.
#' * \strong{AIC}, Akaike information criterion.
#' * \strong{BIC}, Bayesian information criterion.
#' * \strong{lnlike}, maximized log likelihood of the MAR(1) model.
#' @md
#'
#' @importFrom MASS ginv
#'
#'
#' @references Ives, A. R., Dennis, B., Cottingham, K. L., Carpenter, S. R. (2003). Estimating community stability and ecological interactions from time-series data.
#' \emph{Ecological Monographs}, 72(2): 301 - 330
#'
#'
#'
#'

mars.cls	<-	function(comm.mat, covariate = NULL){
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package \"MASS\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Format input
  statevar <- comm.mat[2:nrow(comm.mat), ] # Y
  lagstate <- comm.mat[1:(nrow(comm.mat)-1), ] # X
  Q	<-	length(statevar[ ,1])	# length of the time series
  P	<-	length(statevar[1, ])	# number of species

  if (is.null(covariate) == FALSE) {
    R	<-	length(covariate[1, ])
  } else {
    R <- 0
  }

  # Initiate variates
  E	<-	matrix(0, nrow = Q, ncol = P)	# residuals for least squares
  A	<-	rep(0, P)				# intercepts fo the variate
  B	<-	matrix(0,nrow = P, ncol = P)	# parameters for the variates

  ## Parameters for the covariates
  if (is.null(covariate) == FALSE) {
    C	<-	matrix(0, nrow = P, ncol = R)
  } else{
    C	<-	NULL
  }

  Yhat <-	matrix(0, ncol = P, nrow = Q)	# estimates
  varY <-	matrix(0, ncol = P, nrow = Q)	# variates for total variance
  varDY <-	matrix(0, ncol = P, nrow = Q)	# variates for per capita variance

  # Calculate variates for each species
  for (dv in 1:P){
    Y	<-	statevar[ ,dv]

    # Define predictors
    if (is.null(covariate) == FALSE) {
      X <- as.matrix(cbind(rep(1, Q), lagstate, covariate[-nrow(covariate), ]))
    } else {
      X	<-	as.matrix(cbind(rep(1, Q), lagstate))
    }

    # CLS estimates
    beta	<-	MASS::ginv((t(X) %*% X)) %*% (t(X) %*% Y)

    # Calculate residuals for that dependent variate
    Yhat[,dv]	<-	X %*% beta
    E[,dv] <-	Y - Yhat[ ,dv]

    # Parameter estimates
    A[dv]	<-	beta[1,1]
    B[dv,]	<-	beta[2:(P+1),1]

    if (is.null(covariate) == FALSE) {
      C[dv,]	<-	t(beta[(P+2):nrow(beta), ])
    } else {
      C	<- NULL
    }

    # Accumulate Y variates for calculation of R^2 for observed vs. predicted
    varY[,dv]	<-	Y - mean(Y)

    # Accumulate Y variates for calculation of R^2 for observed change in x vs. predicted
    varDY[,dv]	<-	Y - lagstate[,dv]
  }

  # Assess the CLS fits

  # Calculate CLS log-likelihood
  sigma	<-	t(E) %*% E/Q
  lnlike 	<-	-Q*(P/2)*log(2*pi) - Q/2*log(det(sigma) + 0.001) - Q*P/2

  # Calculate CLS explained variance
  varMatrix	<-	t(varY) %*% varY/Q
  R2	<-	1 - (diag(sigma)/diag(varMatrix))

  # Calculate CLS explained variance in change in x
  varMatrix_D	<-	t(varDY) %*% varDY/Q
  R2_D <- 1 - diag(sigma)/diag(varMatrix_D)

  # Correlation matrix for process error
  d	<- diag(1/sqrt(diag(sigma)))
  corrmatrix	<-	d*sigma*d

  # Define the parameter count
  ## parameter count = # non-zero parameters in A,B,C and some of sigma
  par	<- sum(A!=0) + sum(sum(B!=0)) + P*(P+1)/2 + sum(sum(C!=0))

  # AIC and BIC
  lsAIC	<-	-2*lnlike + 2*par
  lsBIC	<-	-2*lnlike + par*log(Q)

  spp.names <- colnames(comm.mat)
  colnames(B) <- spp.names
  row.names(B) <- spp.names

  results	<-	list(A = A, # A is a vector of constants, maximum growth rate
                   B = B, # matrix of p x p where b_ij gives the effect of the abundance of species j on per capita growth rate of species i
                   sigma = sigma, # env noise variance covariance matrix
                   C = C, # covariant related
                   E = E, # vector representing stochastic env factors that are independent through time
                   Yhat = Yhat, # Estimator of predicted abundance
                   R2 = R2,
                   R2_D = R2_D,
                   AIC = lsAIC,
                   BIC = lsBIC,
                   lnlike = lnlike)
  return(results)
}



