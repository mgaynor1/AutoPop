#' Mature survival alpha and beta calculation.
#'
#' @description Calculates the alpha and beta parameters based on mu, the mean probability of mature survival, and the calculated variance.
#'
#' @param mu Mean probability of mature survival (list).
#' @param var Derived variance based on var.option.
#'
#' @returns List with the parameters alpha and beta.

alphabeta.calc <- function(mu, var){
  alpha <- (mu * (((mu*(1-mu))/var) - 1))
  beta <- ((1-mu) * (((mu*(1-mu))/var) - 1))
  parm.list <- c(alpha, beta)
  return(parm.list)
}
