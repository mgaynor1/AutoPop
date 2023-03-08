#' Simulate a single generation.
#'
#' @description Defined in detail in [Gaynor et al. 2023](). To summarize,
#'   this function is a single time step in a stochastic stage-structured
#'   matrix population dynamics model for diploid, triploid, and autotetraploid
#'   perennial plants with two life-stages (reproductively immature and reproductively mature).
#'   Population composition at time t + 1 is defined by reproduction,
#'   survival, and maturation. Note, lists of three should always have values representing `c(diploids, triploid, autotetraploids)`.
#'
#' @param ct.vec Population composition at time t.
#' Indicates the sum of each type at time t, `ct.vec = c(2x_immature, 3x_immature, 4x_immature, 2x_mature, 3x_mature, 4x_mature)`
#' @param aii.vec The survival probability of an immature individual for each cytotype. Must be a list of three integers between 0 and 1. For example, `aii.vec = c(0.5, 0.3, 0.5)`.
#' @param as.matur The probability of maturation from an immature stage to the mature stage for each cytotype. Must be a list of three integers between 0 and 1. For example, `as.matur = c(0.5, 0.3, 0.5)`.
#' @param as.msurv.e.set The survival probability of a mature individual for each cytotype. Must be a list of three integers between 0 and 1. For example, `as.msurv.e.set = c(0.5, 0.3, 0.5)`.
#'  This list is defined within the `gen.iter.f.choosy()` function.
#' @param d Strength of density dependency on gamete production for each cytotype. Must be a list of three integers between 0 and 1. For example, `d = c(0.001, 0.009, 0.001)`.
#' @param gnum.base  Mean number of gametes per individual per cytotype. Must be a list of three numeric values. For example, `gnum.base = c(100, 100, 100)`.
#' @param mate.lazy Default = FALSE, this prevents selfing from occurring during outcrossing. However, this increases the computational time by 31x!
#' @inheritParams cytotype_repro_mate
#' @param density.type Default = "all", this sets the density at time t as all individuals at time t.
#' Alternatively, "like-cytotype" sets the density at time t for each cytotype based on only the total immature and mature individuals
#' of that cytotype at time t.
#' @returns List of 9, with 1:6 representing the number of individuals of each cytotype
#'   and at both stages at time t + 1. Items 7:9 are the number of gametes sampled for
#'   2x, 3x, and 4x individuals at time t.
#
#'
#' @importFrom stats na.omit rbinom
#'
#'
#'

# One iteration function
one.iter.f.choosy <- function(ct.vec, aii.vec,
                              as.matur, as.msurv.e.set,
                              d, gnum.base,
                              s, b,
                              cc, mc,
                              density.type = "all", mate.lazy = FALSE){

    # Empty vector for next generation
    ctp1.vec <- rep(0, 6)

    # Calculate current generations gnum.vec
    ## gnum.vec = gnum.base times the exponential of the product of strength of density dependency
    ## times the sum of all individuals at time t
    if(density.type == "all"){
      gnum.vec <- c((gnum.base[1]*exp(-d[1]*(sum(ct.vec[1:6])))),
                    (gnum.base[2]*exp(-d[2]*(sum(ct.vec[1:6])))),
                    (gnum.base[3]*exp(-d[3]*(sum(ct.vec[1:6])))))
      gnum.vec[is.na(gnum.vec)] <- 0
    } else if(density.type == "like-cytotype"){
      gnum.vec <- c((gnum.base[1]*exp(-d[1]*(sum(c(ct.vec[1], ct.vec[4]))))),
                    (gnum.base[2]*exp(-d[2]*(sum(c(ct.vec[2], ct.vec[5] ))))),
                    (gnum.base[3]*exp(-d[3]*(sum(c(ct.vec[3], ct.vec[6]))))))
      gnum.vec[is.na(gnum.vec)] <- 0
    } else print("Density type unavailable. This option only allows density.type to be set as all or like-cytotype.")

    # Define offspring for all cytotypes
    if(mate.lazy == FALSE){
      ## When NAs are produced by dhyper, it prints a warning. Our function deals with these NAs, so here we suppress the NA message
      off <- suppressWarnings(cytotype_repro_mate(cgen = ct.vec[4:6], b = b, cc = cc, gnum.vec = gnum.vec, s = s, mc = mc))
    } else if(mate.lazy == TRUE){
      off <- matrix_cytotype_repro_mate(cgen = ct.vec[4:6], b = b, cc = cc, gnum.vec = gnum.vec, s = s, mc = mc)
    }

    dip.gamps <- off[1]; if(is.na(dip.gamps)){dip.gamps <- 0}
    trip.gamps <- off[2]; if(is.na(trip.gamps)){trip.gamps <- 0}
    tetr.gamps <- off[3]; if(is.na(tetr.gamps)){tetr.gamps <- 0}

    # Diploids - two steps: 1. Maturation+Survival and 2. Gamete production (defined above)
    ## Step 1:
    ### Apply survival and maturation
    Ndip.imm <-  rbinom(n = 1, size = ct.vec[1], prob = exp(-aii.vec[1]*sum(ct.vec[1], ct.vec[4]))) #density
    Ndip.mat <-  rbinom(n = 1, size = ct.vec[1], prob = as.matur[1])
    Ndip.surv <- rbinom(n = 1, size = ct.vec[4], prob = as.msurv.e.set[1])

    ## Step 2:
    ### Sum each stage
    ctp1.vec[4] <- Ndip.mat + Ndip.surv
    ctp1.vec[1] <-  dip.gamps + Ndip.imm


    # Triploids
    ## Step 1:
    ### Apply survival and maturation
    Ntrip.imm <-  rbinom(n = 1, size = ct.vec[2], prob = exp(-aii.vec[2]*sum(ct.vec[2], ct.vec[5])))
    Ntrip.mat <-  rbinom(n = 1, size = ct.vec[2], prob = as.matur[2])
    Ntrip.surv <- rbinom(n = 1, size = ct.vec[5], prob = as.msurv.e.set[2])

    ## Step 2:
    ### Sum each stage
    ctp1.vec[5] <-  Ntrip.mat + Ntrip.surv
    ctp1.vec[2] <-  trip.gamps + Ntrip.imm


    # Tetraploids
    ## Step 1:
    ### Apply survival and maturation
    Ntetr.imm <-  rbinom(n = 1, size = ct.vec[3], prob = exp(-aii.vec[3]*sum(ct.vec[3], ct.vec[6])))
    Ntetr.mat <-  rbinom(n = 1, size = ct.vec[3], prob = as.matur[3])
    Ntetr.surv <- rbinom(n = 1, size = ct.vec[6], prob = as.msurv.e.set[3])

    ## Step 2:
    ### Sum each stage
    ctp1.vec[6] <-  Ntetr.mat + Ntetr.surv
    ctp1.vec[3] <-  tetr.gamps + Ntetr.imm

    # Combined and Return
    ctp1.vec <- c(ctp1.vec, off[1:6])
    return(ctp1.vec)
}

