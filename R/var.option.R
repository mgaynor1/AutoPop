#' Mature survival - Calculate variance from env.ci and mu.
#'
#' Calculates variance based on env.ci, the proportion of env. variance, and mu, the mean probability of mature survival.
#'
#'
#' @param env.ci  Proportion of environmental variance used to define mature survival rate per generation
#'   with a beta distribution. This number must be in between 0 and 1, but cannot be equal to 0 or 1.
#' @param mu Mean probability of mature survival (list).
#'
#' @return Calculated variance.
#'
#' @examples
#' variance <- var.option(env.ci = 0.9, mu = 0.6)

var.option <- function(env.ci, mu){
  var <- (env.ci*mu*(1-mu))
  return(var)
}
