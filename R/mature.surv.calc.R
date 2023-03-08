#' Mature survival -  Calculates probability of mature survival.
#'
#' @description Calculates the probability of mature survival based
#' on a beta binomial distribution.
#'
#'
#' @param as.msurv Mean probability of mature survival. Must be a single integer between 0 and 1.
#' @param  env.ci Proportion of environmental variance used to define mature survival rate per generation.
#'   Must be an integer greater than or equal to 0 and less than 1.
#'
#' @returns Probability of mature survival.
#'
#' @importFrom stats na.omit rbeta
#'


mature.surv.calc <- function(env.ci, as.msurv){
  var <- var.option(env.ci = env.ci, mu = as.msurv)
  alpha.beta <- alphabeta.calc(mu = as.msurv, var = var)
  absamp <- rbeta(n = 1, alpha.beta[1], alpha.beta[2])
  return(absamp)
}


