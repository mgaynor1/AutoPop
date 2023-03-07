#' Mature survival calculation.
#'
#' @description Calculates the probability of mature survival based on a beta binomial distribution.
#'
#'
#' @param as.msurv Mean probability of mature survival (list).
#' @param  env.ci Proportion of environmental variance used to define mature survival rate per generation
#'   with a beta distribution. This number must be in between 0 and 1, but cannot be equal to 0 or 1.
#'
#' @returns Probability of mature survival.
#'
#' @importFrom stats na.omit rbeta



mature.surv.calc <- function(env.ci, as.msurv){
  var <- var.option(env.ci = env.ci, mu = as.msurv)
  alpha.beta <- alphabeta.calc(mu = as.msurv, var = var)
  absamp <- rbeta(n = 1, alpha.beta[1], alpha.beta[2])
  return(absamp)
}


