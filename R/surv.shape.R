#' Survival shape -  Calculates the alpha and beta value for survival beta distributions
#'
#' @description Calculates the shape parameters needed to sample a beta distribution
#' for both immature survival and mature survival
#'
#'
#' @param raw.means Mean probability of mature survival or immature survival. Must be a list of three integer between 0 and 1.
#' @param  env.ci Proportion of environmental variance used to define mature survival rate per generation.
#'   Must be a single integer greater than or equal to 0 and less than 1.
#'
#' @returns Matrix of alpha and beta parameters where matrix[1,i] is the alpha value for the ith mean, and matrix[2,i]
#' is the beta value for the ith mean.
#'
#' @importFrom stats na.omit rbeta
#'


surv.shape <- function(env.ci, raw.means){
  var <- var.option(env.ci = env.ci, mu = raw.means)
  alpha.beta <- sapply(1:length(raw.means), function(item) alphabeta.calc(raw.means[item], var[item]))
  return(alpha.beta)
}


