#' Mature survival - Calculates mean and variance from alpha and beta.
#'
#' @description This function can be used to convert alpha and beta back to mu and var.
#'
#'
#' @param alpha See alphabeta.calc.
#' @param beta See alphabeta.calc.
#'
#' @returns List with the parameters mu and var.
#'

muvar.calc <- function(alpha, beta){
  mu <- (alpha/(alpha+beta))
  var <- ((alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1)))
  parm.list <- c(mu, var)
  return(parm.list)
}
