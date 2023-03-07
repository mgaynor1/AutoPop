#' Simulate multiple generations
#'
#' @description This is a function that loops through a set number of generations
#'   to look at population dynamics overtime.
#'
#'
#' @param generations Number of generations to iterate.
#' @param init.pop Initial population size, ie. number of mature diploids in founding population.
#' @param  env.ci Proportion of environmental variance used to define mature survival rate per generation
#'   with a normal distribution where the mean is the previous generations maturation rate.
#' @param as.msurv probability of mature survival (list).
#' @inheritParams one.iter.f.choosy
#' @inheritParams cytotype_repro_mate
#' @inheritParams format.iter
#'
#' @returns A single data frame as defined by `format.iter()`.
#'



#source("R/one-iter-f-choosy.R")
#source("R/format-iter-choosy.R")
#source("R/mature.surv.calc.R")

gen.iter.f.choosy <- function(generations, init.pop, env.ci, aii.vec,
                              as.matur, as.msurv, d, gnum.base,
                              b, cc, s, mc, density.type = "all",
                              mate.lazy = FALSE){

  ## Set up first generation
  popvect <- list()
  first.vec <- c(0, 0, 0, init.pop, 0, 0, 0, 0, 0, 0, 0, 0)
  popvect[[1]] <- first.vec
  as.msurv.e <- list()
  as.msurv.e[[1]] <- as.msurv


  ## Loop for set number of generations #########################################
  for(i in 1:generations){
    # Environmental. Stochasticity
    m1e <- mature.surv.calc(env.ci = env.ci, as.msurv = as.msurv.e[[1]][1])
    m2e <- mature.surv.calc(env.ci = env.ci, as.msurv = as.msurv.e[[1]][2])
    m3e <- mature.surv.calc(env.ci = env.ci, as.msurv = as.msurv.e[[1]][3])
    as.msurv.e[[i + 1]] <- c(m1e, m2e, m3e)
    ct.vec <- popvect[[i]]
    popvect[[i + 1]] <- one.iter.f.choosy(ct.vec, aii.vec = aii.vec , as.matur = as.matur,
                                          as.msurv.e.set = as.msurv.e[[i + 1]], d = d, gnum.base = gnum.base,
                                          b = b, cc = cc, s = s, mc = mc, density.type, mate.lazy = mate.lazy)
  }


  onedf <- format.iter(popvect, generations)
  return(onedf)
}
