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
#' @inheritParams env.copula
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
#' * sum2x: sum of immature and mature diploids.
#' * sum3x: sum of immature and mature triploids.
#' * sum4x: sum of immature and mature tetraploids.
#' * V1a: relative abundance of immature diploids (V1).
#' * V2a: relative abundance of immature triploids (V2).
#' * V3a: relative abundance of immature tetraploids (V3).
#' * V4a: relative abundance of mature diploids (V4).
#' * V5a: relative abundance of mature triploids (V5).
#' * V6a: relative abundance of mature tetraploids (V6).
#' * C2: relative abundance of all diploids (sum2x/sum).
#' * C3: relative abundance of all triploids (sum3x/sum).
#' * C4: relative abundance of all tetraploids (sum4x/sum)
#' * i1e: diploid immature survival probability during t - 1.
#' * i2e: triploid immature survival probability during t - 1.
#' * i3e: tetraploid immature survival probability during t - 1.
#' * m1e: diploid mature survival probability during t - 1.
#' * m2e: triploid mature survival probability during t - 1.
#' * m3e: tetraploid mature survival probability during t - 1.


#source("R/one-iter-f-choosy.R")
#source("R/format-iter-choosy.R")
#source("R/mature.surv.calc.R")

gen.iter.f.choosy <- function(generations, init.pop, env.ci, env.sigma, aii.vec,
                              as.matur, as.msurv, d, gnum.base,
                              b, cc, s, mc,
                              gam.density.type = "all", is.density.type = "all",  rho, mate.lazy = FALSE,
                              rj = c(1,1),  yj = c(1,1),  env.immature.survival = FALSE){

  ## Set up first generation
  popvect <- list()
  first.vec <- c(0, 0, 0, init.pop, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  popvect[[1]] <- first.vec

  as.msurv.e <- list()
  as.msurv.e[[1]] <- c(0,0,0)
  Xi.Ui.vec <- env.copula(rho = rho, generations = generations)
  Xi.vec <- Xi.Ui.vec[[1]]
  Ui.vec <- Xi.Ui.vec[[2]]
  if(sum(as.msurv) != 0){
  msurv.alpha.beta <- surv.shape(env.ci = env.ci, raw.means = as.msurv)
  }
  ## Loop for set number of generations #########################################
  for(i in 1:generations){
    if(sum(as.msurv) != 0){
    # Mature Survival - Environmental. Stochasticity
    m1e <- qbeta(Ui.vec[i,1], msurv.alpha.beta[1,1], msurv.alpha.beta[2,1])
    m2e <- qbeta(Ui.vec[i,2], msurv.alpha.beta[1,2], msurv.alpha.beta[2,2])
    m3e <- qbeta(Ui.vec[i,3], msurv.alpha.beta[1,3], msurv.alpha.beta[2,3])
    as.msurv.e[[i + 1]] <- c(m1e, m2e, m3e)
    } else{
      as.msurv.e[[i + 1]] <- c(0, 0, 0)
    }
    ct.vec <- popvect[[i]]
    popvect[[i + 1]] <- one.iter.f.choosy(ct.vec, aii.vec = aii.vec, as.matur = as.matur,
                                          as.msurv.e.set = as.msurv.e[[i + 1]], d = d, gnum.base = gnum.base,
                                          b = b, cc = cc, s = s, mc = mc, mate.lazy = mate.lazy,
                                          env.sigma = env.sigma, Xi.is = Xi.vec[i,],  gam.density.type = gam.density.type,
                                          is.density.type = is.density.type,
                                          rj = rj,
                                          yj = yj, env.immature.survival = env.immature.survival)
  }


  onedf <- form.autopop(popvect, generations)
  survival_rates <- as.data.frame(do.call(rbind, popvect))[1:(generations+1),13:15]
  for(i in 1:generations){
    onedf$i1e[i] <- survival_rates$V13[i]
    onedf$i2e[i] <- survival_rates$V14[i]
    onedf$i3e[i] <- survival_rates$V15[i]
    onedf$m1e[i] <- as.msurv.e[[i]][1]
    onedf$m2e[i] <- as.msurv.e[[i]][2]
    onedf$m3e[i] <- as.msurv.e[[i]][3]
  }

  return(onedf)
}
