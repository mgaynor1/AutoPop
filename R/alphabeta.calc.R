#' Mature survival - Calculates alpha and beta from mean and variance.
#'
#' @description Calculates the alpha and beta parameters based on mu, here the mean probability of mature survival,
#' and s calculated variance.
#'
#' @param mu Sample mean. Must be a single integer between 0 and 1.
#' @param var Derived variance based on var.option. Must be a single integer between 0 and 1.
#'
#' @returns List with the parameters alpha and beta.
#'


alphabeta.calc <- function(mu, var){
  alpha <- (mu * (((mu*(1-mu))/var) - 1))
  beta <- ((1-mu) * (((mu*(1-mu))/var) - 1))
  parm.list <- c(alpha, beta)
  return(parm.list)
}
