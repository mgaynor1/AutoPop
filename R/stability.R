#' Stability metrics calculator following Ives et al. 2003
#'
#' @description  Based on Ives et al. 2003, we calculate four stability  metrics based on the estimates of the matrix of ecological interactions
#' B and the estimated variance-covariance matrix of environmental noise 'Sigma'.  These two quantities can be read from the output of the function `mars.cls()`.
#'
#' @param B Matrix of p x p (p = number of species), where \eqn{b_{ij}} gives the effect of the abundance of species j on per capita growth rate of species.
#'   This can be calculated with the `mars.cls()` function.
#' @param sigma Matrix of p x p, environmental noise variance-covariance matrix. This can be calculated with the `mars.cls()` function.
#'
#' @returns Resulting list which includes:
#' * \strong{var.prop}, at stationarity the var.prop is the variance proportion attributable to environmental noise. The smaller the values, the more stable the dynamics.
#' * \strong{mean.return.time} & \strong{var.return.time}, rate at which the transition distribution converges back to the stationary distribution. The less time it takes to return to the stationary distribution, the more stable the population.
#' * \strong{reactivity}, measures reaction to perturbations or the distance away from stationary a system moves in response to a disturbance. Again, smaller is better in terms of stability.
#' * \strong{sp.contribs}, squared eigenvalue representing the characteristic return rate of the variance of the transition distribution of the estimated MARS(1) Markov Process.
#'
#'
#' @references Ives, A. R., Dennis, B., Cottingham, K. L., Carpenter, S. R. (2003). Estimating community stability and ecological interactions frm time-series data. \emph{Ecological Monographs}, 72(2): 301 - 330.
#'
#'


stability	<- function(B, sigma){

  B	<-	as.matrix(B)
  nspecies <- nrow(B)
  stab1 <- det(B)^(2) # defines var.prop
  lams.vec <- Re(eigen(B)$values) # characteristic return rate of the mean of the transition distribution to the mean of the stationary distribution
  mean.return <- max(lams.vec) # return rate for the mean
  var.return <- max(Re(eigen(kronecker(B,B))$values)) # return rate for the variance
  Vinf.true  <- matrix(solve(diag(nspecies*nspecies) - kronecker(B,B))%*%as.vector(sigma) , nrow=nspecies,ncol=nspecies)
  stab3 <- -sum(diag(B))/sum(diag(Vinf.true))
  spp.contribs <- (1/lams.vec^2)/sum(1/lams.vec^2)
  names(spp.contribs) <- colnames(B)

  if(nrow(B) == 2) {
    return(data.frame(var.prop=stab1,
                mean.return.time=mean.return,
                var.return.time=var.return,
                reactivity=stab3,
                sp1.contribs = spp.contribs[1],
                sp2.contribs = spp.contribs[2]))
  } else {
  return(list(var.prop=stab1,
              mean.return.time=mean.return,
              var.return.time=var.return,
              reactivity=stab3,
              sp.contribs = spp.contribs))
  }
}

