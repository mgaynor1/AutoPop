#' Simulate a single generation.
#'
#' @description Defined in detail in [Gaynor et al. 2025](https://doi.org/10.1086/734411). To summarize,
#'   this function is a single time step in a stochastic stage-structured
#'   matrix population dynamics model for diploid, triploid, and autotetraploid
#'   perennial plants with two life-stages (reproductively immature and reproductively mature).
#'   Population composition at time t + 1 is defined by reproduction,
#'   survival, and maturation. Note, lists of three should always have values representing `c(diploids, triploid, autotetraploids)`.
#'
#' @param ct.vec Population composition at time t.
#' Indicates the sum of each type at time t, `ct.vec = c(2x_immature, 3x_immature, 4x_immature, 2x_mature, 3x_mature, 4x_mature)`
#' @param aii.vec The survival probability of an immature individual for each cytotype. Must be a list of three integers between 0 and 1. For example, `aii.vec = c(0.0005, 0.005, 0.0005)`.
#' @param as.matur The probability of maturation from an immature stage to the mature stage for each cytotype. Must be a list of three integers between 0 and 1. For example, `as.matur = c(0.60, 0.06, 0.60)`.
#' @param as.msurv.e.set The survival probability of a mature individual for each cytotype. Must be a list of three integers between 0 and 1. These are defined based on user supplied
#'  mean probability of mature survival (`as.msurv`), proportion of environmental variance (`env.ci`),
#'  This list is defined within the `gen.iter.f.choosy()` function.
#' @param d Strength of density dependency on gamete production for each cytotype. Must be a list of three integers between 0 and 1. For example, `d = c(0.001, 0.009, 0.001)`.
#' @param gnum.base  Mean number of gametes per individual per cytotype. Must be a list of three numeric values. For example, `gnum.base = c(100, 100, 100)`.
#' @param mate.lazy Default = FALSE, this prevents selfing from occurring during outcrossing. However, this increases the computational time by 31x!
#' @param gam.density.type Default = "all", this sets the density dependence type for number of gametes produced at time t as all individuals at time t.
#' Alternatively, "like-cytotype" sets the density at time t for each cytotype based on only the total immature and mature individuals
#' of that cytotype at time t.
#' @param rj List of two indicates the stage dependent density dependent impact of immature and mature individuals on  number of gametes produced. Default is c(1,1).
#' @param is.density.type Default = "all", this sets the density dependence type for immature surivival at time t as all individuals at time t.
#' Alternatively, "like-cytotype" sets the density at time t for each cytotype based on only the total immature and mature individuals
#' of that cytotype at time t.
#' @param yj List of two indicates the stage dependent density dependent impact of immature and mature individuals on the probability of immature survival. Default is c(1,1).
#' @param env.sigma Sigma value to define environmental variance corresponding to immature survival rate per generation.
#' @param Xi.is Sample from a multivariate normal distribution. Allows immature and mature survival probabilities be correlated,
#' as well as allowing for a correlation among cytotypes based on a supplied `rho` value.
#' @param env.immature.survival Default = FALSE. When FALSE, the mean immature survival rate is equal to the exponential of the inverse
#' immature survival rate (`aii`) times the sum of immature individuals and mature individuals of the cytotype indicated by `is.density.type`,
#' which may be modified by `yj`. When this variable equals TRUE, the previous mean, or the immature determinate survival rate, is transformed
#' by log(mu/(1-mu)) and the Johnson SB distribution is sampled given the `env.sigma` and `Xi.is` (defined in `gen.iter.f.choosy()`).
#' @inheritParams cytotype_repro_mate
#' @returns List of 15, with 1:6 representing the number of individuals of each cytotype
#'   and at both stages at time t + 1. Items 7:9 are the number of offspring generated for
#'   2x, 3x, and 4x individuals at time t. Items 10:12 are the number of gametes generated from
#'   2x, 3x, and 4x individuals at time t. Items 12:15 are the immature survival probabilities
#'   at time t.
#'
#'
#' @importFrom stats na.omit rbinom
#'
#'
#'

# One iteration function
one.iter.f.choosy <- function(ct.vec, aii.vec, as.matur, as.msurv.e.set,
                              d, gnum.base, s, b,
                              cc, mc, gam.density.type, is.density.type,
                              mate.lazy = FALSE, rj = c(1,1), yj = c(1,1), env.sigma,
                              Xi.is, env.immature.survival){

    # Empty vector for next generation
    ctp1.vec <- rep(0, 6)

    # Calculate current generations gnum.vec
    ## gnum.vec = gnum.base times the exponential of the product of strength of density dependency
    ## times the sum of all individuals at time t
    if(gam.density.type == "all"){
      gam.density <- ((rj[1]*sum(ct.vec[1:3]))+(rj[2]*sum(ct.vec[4:6])))
      gnum.vec <- c((gnum.base[1]*exp(-d[1]*gam.density)),
                    (gnum.base[2]*exp(-d[2]*gam.density)),
                    (gnum.base[3]*exp(-d[3]*gam.density)))
      gnum.vec[is.na(gnum.vec)] <- 0
    } else if(gam.density.type == "like-cytotype"){
      gnum.vec <- c((gnum.base[1]*exp(-d[1]*((rj[1]*ct.vec[1]) + (rj[2]*ct.vec[4])))),
                    (gnum.base[2]*exp(-d[2]*((rj[1]*ct.vec[2]) + (rj[2]*ct.vec[5])))),
                    (gnum.base[3]*exp(-d[3]*((rj[1]*ct.vec[3]) + (rj[2]*ct.vec[6])))))
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

    # Setup Immature Survival Probability
    if(is.density.type == "all"){
      is.density <- ((yj[1]*sum(ct.vec[1:3]))+(yj[2]*sum(ct.vec[4:6])))
      immature.det.surv <- c((exp(-aii.vec[1]*is.density)),
                                (exp(-aii.vec[2]*is.density)),
                                (exp(-aii.vec[3]*is.density)))
    } else if(is.density.type == "like-cytotype"){
      immature.det.surv <- c((exp(-aii.vec[1]*((yj[1]*ct.vec[1]) + (yj[2]*ct.vec[4])))),
                                (exp(-aii.vec[2]*((yj[1]*ct.vec[2]) + (yj[2]*ct.vec[5])))),
                                (exp(-aii.vec[3]*((yj[1]*ct.vec[3]) + (yj[2]*ct.vec[6])))))
    } else print("Density type unavailable. This option only allows density.type to be set as all or like-cytotype.")

    immature.det.surv[is.na(immature.det.surv)] <- 0
    if(env.immature.survival == TRUE){
    immature.mu<- c()
    as.isurv <- c()
    for(z in 1:3){
      immature.mu[z] <- log(immature.det.surv[z]/(1-immature.det.surv[z]))
      as.isurv[z] <- (exp(immature.mu[z] + (env.sigma*Xi.is[z]))/(1 + exp(immature.mu[z] + (env.sigma*Xi.is[z]))))
    }
    as.isurv[is.na(as.isurv)] <- 0
    } else{
      as.isurv <- immature.det.surv
    }



    # Diploids - two steps: 1. Maturation+Survival and 2. Gamete production (defined above)
    ## Step 1:
    ### Apply survival and maturation
    Ndip.imm <-  rbinom(n = 1, size = ct.vec[1], prob = as.isurv[1]) #density
    Ndip.mat <-  rbinom(n = 1, size = ct.vec[1], prob = as.matur[1])
    Ndip.surv <- rbinom(n = 1, size = ct.vec[4], prob = as.msurv.e.set[1])

    ## Step 2:
    ### Sum each stage
    ctp1.vec[4] <- Ndip.mat + Ndip.surv
    ctp1.vec[1] <-  dip.gamps + Ndip.imm


    # Triploids
    ## Step 1:
    ### Apply survival and maturation
    Ntrip.imm <-  rbinom(n = 1, size = ct.vec[2], prob = as.isurv[2])
    Ntrip.mat <-  rbinom(n = 1, size = ct.vec[2], prob = as.matur[2])
    Ntrip.surv <- rbinom(n = 1, size = ct.vec[5], prob = as.msurv.e.set[2])

    ## Step 2:
    ### Sum each stage
    ctp1.vec[5] <-  Ntrip.mat + Ntrip.surv
    ctp1.vec[2] <-  trip.gamps + Ntrip.imm


    # Tetraploids
    ## Step 1:
    ### Apply survival and maturation
    Ntetr.imm <-  rbinom(n = 1, size = ct.vec[3], prob = as.isurv[3])
    Ntetr.mat <-  rbinom(n = 1, size = ct.vec[3], prob = as.matur[3])
    Ntetr.surv <- rbinom(n = 1, size = ct.vec[6], prob = as.msurv.e.set[3])

    ## Step 2:
    ### Sum each stage
    ctp1.vec[6] <-  Ntetr.mat + Ntetr.surv
    ctp1.vec[3] <-  tetr.gamps + Ntetr.imm

    ctp1.vec[is.na(ctp1.vec)] <- 0
    # Combined and Return
    ctp1.vec <- c(ctp1.vec, off[1:6], as.isurv)
    return(ctp1.vec)
}

