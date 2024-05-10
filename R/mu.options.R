#' Max mean allowed
#'
#' Given an env.ci value, we identified allowed means between 0 and 1.
#'
#' @param env.ci  Proportion of environmental variance used to define mature survival rate per generation
#'   with a beta distribution. This number must be in between 0 and 1, but cannot be equal to 0 or 1.
#'
#' @return All values allowed and not allowed.
#'

mu.options <- function(env.ci){
  p.mu <- seq(0.001,1, by = 0.0001)
  notallowed <- c()
  allowed <- c()
  for(i in 1:length(p.mu)){

    mu <- p.mu[i]
    var <- (env.ci*(mu*(1-mu)))
    less <- (mu*(mu-1))

    if((var > less) == FALSE){
      notallowed <- c(notallowed, mu)
    }else{
      allowed <- c(allowed, mu)
    }
  }
 out <-  list(not_allowed = notallowed, allowed = allowed)
 return(out)
}

