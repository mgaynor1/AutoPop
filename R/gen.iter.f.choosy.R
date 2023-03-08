#' Simulate multiple generations.
#'
#' @description Defined in detail in [Gaynor et al. 2023](). To summarize,
#'   this function simulates a stochastic stage-structured
#'   matrix population dynamics model for diploid, triploid, and autotetraploid
#'   perennial plants with overlapping generations and two life-stages (reproductively immature and reproductively mature).
#'   Population composition at time t + 1 is defined by reproduction,
#'   survival, and maturation. This is a function that loops through a set number of generations.
#'   Note, lists of three should always have values representing `c(diploids, triploid, autotetraploids)`.
#'
#'
#'
#' @param generations Number of generations to simulate. Must be a numeric value.
#' @param init.pop The number of mature diploids in the initial founding population. Must be a numeric value greater than 0.
#' @param  env.ci Proportion of environmental variance used to define mature survival rate per generation.
#'   Must be an integer greater than or equal to 0 and less than 1.
#' @param as.msurv The mean survival probability of a mature individual for each cytotype.
#'    Must be a list of three integers between 0 and 1. For example, `as.msurv = c(0.5, 0.3, 0.5)`.
#' @inheritParams one.iter.f.choosy
#' @inheritParams cytotype_repro_mate
#' @inheritParams form.autopop
#'
#' @returns A single data frame as defined by `form.autopop()`. Each row is a generation. The columns are as follows,
#' * V1: number of immature diploids.
#' * V2: number of immature triploids.
#' * V3: number of immature tetraploids.
#' * V4: number of mature diploids.
#' * V5: number of mature triploids.
#' * V6: number of mature tetraploids.
#' * V7: number of diploid offspring produced during t - 1.
#' * V8: number of triploid offspring produced during t - 1.
#' * V9: number of tetraploid offspring produced during t - 1.
#' * V10: number of gametes per diploid individual at t - 1.
#' * V11: number of gametes per triploid individual at t - 1.
#' * V12: number of gametes per tetraploid individual at t - 1.
#' * gen: generation.
#' * sum: total number of individuals.
#' * sum2x-sum4x: total number of each cytotype at that generation.
#' * V1a:V6a: relative abundance of V1:V6 in that generation.
#' * C2: relative abundance of all diploids (ie. sum2x/sum).
#' * C3: relative abundance of all triploids (ie. sum3x/sum).
#' * C4: relative abundance of all tetraploids (ie. sum4x/sum)
#'
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


  onedf <- form.autopop(popvect, generations)
  return(onedf)
}
