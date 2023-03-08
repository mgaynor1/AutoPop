#' Lazy reproduction function.
#'
#' @description This is a function to generate the number of offspring produced per generation.
#'   Assumes a mixed mating system where both selfing and outcrossing occur.
#'   Also allows group-based assortative mating (see Otto et al. 2008).
#'
#' @param cgen Number of mature individuals of each cytotype in the current generation.
#'   Must be a list of three numeric values. For example, `cgen = c(100, 100, 100)`.
#' @param b Proportion of unreduced gamete formed by each diploid and tetraploid individual.
#'    Must be a single integer between 0 and 1.
#' @param cc Proportion of 3n gamete formation from each triploid individual.
#'   Must be a single integer between 0 and 1.
#' @param gnum.vec The mean number of gametes produced by each individual for each cytotype.
#'   Must be a list of three numeric values. For example, `gnum.vec = c(100, 100, 100)`.
#' @param s Selfing rate. Must be a single integer between 0 and 1.
#' @param mc Strength of mating choice. Must be a single integer between 0 and 1.
#'
#' @return Output is a generation log including number of 2x offspring, 3x offspring, 4x offspring, gametes sampled per 2x individual, gametes sampled per 3x individual, and gametes sampled per 4x individual.
#'
#'
#' @importFrom stats na.omit rpois rmultinom
#'
#'
#' @references Otto, S. P., Servedio, M. R., Nuismer, S. L. (2008). Frequency-dependent selection and the evolution of assortative mating. \emph{Genetics}, 179(4) 2091 - 2112. doi: 10.1534/genetics.107.084418
#'
#'
#'
#'

matrix_cytotype_repro_mate <- function(cgen, b, cc, gnum.vec, s, mc) {

  # Set number of individuals
  ndi <- as.numeric(cgen[1])
  nti <- as.numeric(cgen[2])
  nai <- as.numeric(cgen[3])

  # Set number of gametes per individual based on cytotype
  ngd <- as.numeric(rpois(n = 1, lambda = gnum.vec[1]))
  ngt <- as.numeric(rpois(n = 1, lambda = gnum.vec[2]))
  nga <- as.numeric(rpois(n = 1, lambda = gnum.vec[3]))

  # If no adults exist or no gametes are made, return 0s
  if (ngd == 0 & ngt == 0 & nga == 0) {
    genlog <- c(0, 0, 0, ngd, ngt, nga)
  } else if (ndi == 0 & nti == 0 & nai == 0) {
    genlog <- c(0, 0, 0, ngd, ngt, nga)
  } else{

    # Assign gametes a copy number
    ## Diploids produce 1n and 2n, Triploids produce 3n, and Autotetraploids produce 2n gametes.
    D1n <-  floor((ndi*ngd)*(1 - b)); if (is.na(D1n)) {D1n <- 0} # 1n from diploid
    D2n <-  floor((ndi*ngd)*b); if (is.na(D2n)) {D2n <- 0} # 2n from diploid
    T3n <-  floor((nti*ngt)*cc); if (is.na(T3n)) {T3n <- 0} # 3n from triploid
    A2n <-  floor((nai*nga)*(1 - b)); if (is.na(A2n)) {A2n <- 0} # 2n from tetraploid

    # Define gamete pools
    ## Diploid
    prob1nself <- (s*(D1n/(D1n + D2n))); if (is.na(prob1nself)) {prob1nself <- 0}
    prob2nself <- (s*(D2n/(D1n + D2n))); if (is.na(prob2nself)) {prob2nself <- 0}

    d1nself <- floor(D1n*prob1nself); if (is.na(d1nself)) {d1nself <- 0}
    d2nself <- floor(D2n*prob2nself); if (is.na(d2nself)) {d2nself <- 0}

    d1nout <- (D1n - d1nself); if (is.na(d1nout)) {d1nout <- 0}
    d2nout <- (D2n - d2nself); if (is.na(d2nout)) {d2nout <- 0}

    ## Triploid
    t3self <- floor(T3n*s); if (is.na(t3self)) {t3self <- 0 }
    t3nout <- (T3n - t3self); if (is.na(t3nout)) {t3nout <- 0 }

    ## Autotetraploid
    a2nself <- floor(A2n*s);  if (is.na(a2nself)) {a2nself <- 0 }
    a2nout <- (A2n - a2nself);  if (is.na(a2nout)) {a2nout <- 0 }

    # Pairing of gametes
    # Selfing
    ## prior-selfing rather than delayed selfing
    selfed2x <- as.integer(floor(d1nself/2));  if (is.na(selfed2x)) {selfed2x <- 0}
    selfed4x <- as.integer(sum(na.omit(c(floor(a2nself/2), floor(d2nself/2)))));  if (is.na(selfed4x)) {selfed4x <- 0}

    # Outcrossing
    ## 2 mating types are allowed
    ## (1) Outcrossing
    if (mc == 0) {
      ## Summarize
      (totaloutcrossing <- sum(na.omit(c(d1nout, d2nout, t3nout, a2nout))));  if (is.na(totaloutcrossing)) {totaloutcrossing <- 0}
      if (totaloutcrossing != 0) {
          p1nd <- (d1nout/totaloutcrossing);  if (is.na(p1nd)) {p1nd <- 0}
          p2nd <- (d2nout/totaloutcrossing);  if (is.na(p2nd)) {p2nd <- 0}
          p3nt <- (t3nout/totaloutcrossing);  if (is.na(p3nt)) {p3nt <- 0}
          p2na <- (a2nout/totaloutcrossing);  if (is.na(p2na)) {p2na <- 0}

          ## Random mating
          rowcol <- c("1n", "2n", "3n", "2n")
          mates <- matrix(nrow = 4, ncol = 4)
          row.names(mates) <- rowcol
          colnames(mates) <- rowcol

          mates[1,1] <- p1nd * p1nd # 2x #Like
          mates[1,2] <- p2nd * p1nd # 3x #Like
          mates[1,3] <- p3nt * p1nd # 4x
          mates[1,4] <- p2na * p1nd # 3x

          mates[2,1] <- p1nd * p2nd # 3x #Like
          mates[2,2] <- p2nd * p2nd # 4x #Like
          mates[2,3] <- p3nt * p2nd # 0
          mates[2,4] <- p2na * p2nd # 4x

          mates[3,1] <- p1nd * p3nt # 4x
          mates[3,2] <- p2nd * p3nt # 0
          mates[3,3] <- p3nt * p3nt # 0 #Like
          mates[3,4] <- p2na * p3nt # 0

          mates[4,1] <- p1nd * p2na # 3x
          mates[4,2] <- p2nd * p2na # 4x
          mates[4,3] <- p3nt * p2na # 0
          mates[4,4] <- p2na * p2na # 4x #Like

          # Multinomial samples
          mates.vec <- as.vector(mates)
          mates.vec[is.na(mates.vec)] <- 0

            if (sum(mates.vec) != 0) {
              out.t <- rmultinom(n = 1, size = abs(floor(totaloutcrossing/2)), prob = mates.vec)

              outcross2x <- out.t[1] ; if (is.na(outcross2x)) {outcross2x <- 0}
              outcross3x <- sum(na.omit(c(out.t[2], out.t[4], out.t[5], out.t[13]))); if (is.na(outcross3x)) {outcross3x <- 0}
              outcross4x <- sum(na.omit(c(out.t[3], out.t[6], out.t[8], out.t[9], out.t[14], out.t[16]))); if (is.na(outcross4x)) {outcross4x <- 0}

            } else{
              outcross2x <- 0
              outcross3x <- 0
              outcross4x <- 0
            }
     } else{
       outcross2x <- 0
       outcross3x <- 0
       outcross4x <- 0
      }
    ## (2) Choosy mating. ie. assortative mating
    } else if (mc > 0 & mc != 0) {
      ## Define pools
      ## 2x only pool
      tot2x <- sum(na.omit(c(d1nout, d2nout)));  if (is.na(tot2x)) {tot2x <- 0}
      if (tot2x != 0) {
        pr1nd_prob <- mc*(d1nout/(tot2x));  if (is.na(pr1nd_prob)) {pr1nd_prob <- 0}
        pr2nd_prob <- mc*(d2nout/(tot2x));  if (is.na(pr2nd_prob)) {pr2nd_prob <- 0}

        pool2x_1n <- floor(pr1nd_prob*(d1nout));  if (is.na(pool2x_1n)) {pool2x_1n <- 0}
        pool2x_2n <- floor(pr2nd_prob*(d2nout));  if (is.na(pool2x_2n)) {pool2x_2n <- 0}

        tot2x_pool <- sum(na.omit(c(pool2x_1n, pool2x_2n)));  if (is.na(tot2x_pool)) {tot2x_pool <- 0}
        p2x_prob1n <- (pool2x_1n/tot2x_pool); if (is.na(p2x_prob1n)) {p2x_prob1n <- 0}
        p2x_prob2n <- (pool2x_2n/tot2x_pool); if (is.na(p2x_prob2n)) {p2x_prob2n <- 0 }

        p2x_mates <- abs(as.vector(c((p2x_prob1n*p2x_prob1n), (p2x_prob1n*p2x_prob2n), (p2x_prob1n*p2x_prob2n), (p2x_prob2n*p2x_prob2n)))) #2x, 3x, 3x, 4x
        p2x_mates[is.na(p2x_mates)] <- 0

        if (sum(p2x_mates) != 0) {
            p2x_out.t <- rmultinom(n = 1, size = abs(floor(tot2x_pool/2)), prob = p2x_mates)
            p2x_2x <- p2x_out.t[1]; if (is.na(p2x_2x)) {p2x_2x <- 0}
            p2x_3x <- sum(na.omit(c(p2x_out.t[2], p2x_out.t[3]))); if (is.na(p2x_3x)) {p2x_3x <- 0}
            p2x_4x <- p2x_out.t[4]; if (is.na(p2x_4x)) {p2x_4x <- 0}
        } else{
            p2x_2x <- 0
            p2x_3x <- 0
            p2x_4x <- 0
        }
     } else{
        pr1nd_prob <- 0
        pr2nd_prob <- 0
        p2x_2x <- 0
        p2x_3x <- 0
        p2x_4x <- 0
      }

      ## 3x only
      #pool3x <- floor(mc*(t3nout)); if (is.na(pool3x)) { pool3x <- 0 }

      ## 4x only
      pool4x <- floor(mc*(a2nout)); if (is.na(pool4x)) {pool4x <- 0 }
      p4x_4x <- as.integer(pool4x/2); if (is.na(p4x_4x)) {p4x_4x <- 0}

      ## Random Mating
      poolrm_1nd <- floor((1 - pr1nd_prob)*(d1nout));  if (is.na(poolrm_1nd)) {poolrm_1nd <- 0}
      poolrm_2nd <- floor((1 - pr2nd_prob)*(d2nout));  if (is.na(poolrm_2nd)) {poolrm_2nd <- 0}
      poolrm_3nt <-  floor((1 - mc)*(t3nout));  if (is.na(poolrm_3nt)) {poolrm_3nt <- 0}
      poolrm_2na <- floor((1 - mc)*(a2nout));  if (is.na(poolrm_2na)) {poolrm_2na <- 0}

      total_poolrm <- sum(na.omit(c(poolrm_1nd, poolrm_2nd, poolrm_3nt, poolrm_2na))) ; if (is.na(total_poolrm)) {total_poolrm <- 0}
       if (total_poolrm != 0) {
        p1nd <- (poolrm_1nd/total_poolrm);  if (is.na(p1nd)) {p1nd <- 0}
        p2nd <- (poolrm_2nd/total_poolrm);  if (is.na(p2nd)) {p2nd <- 0}
        p3nt <- (poolrm_3nt/total_poolrm);  if (is.na(p3nt)) {p3nt <- 0}
        p2na <- (poolrm_2na/total_poolrm);  if (is.na(p2na)) {p2na <- 0}

        ## Random mating
        rowcol <- c("1n", "2n", "3n", "2n")
        mates <- matrix(nrow = 4, ncol = 4)
        row.names(mates) <- rowcol
        colnames(mates) <- rowcol

        mates[1,1] <- p1nd * p1nd # 2x #Like
        mates[1,2] <- p2nd * p1nd # 3x #Like
        mates[1,3] <- p3nt * p1nd # 4x
        mates[1,4] <- p2na * p1nd # 3x

        mates[2,1] <- p1nd * p2nd # 3x #Like
        mates[2,2] <- p2nd * p2nd # 4x #Like
        mates[2,3] <- p3nt * p2nd # 0
        mates[2,4] <- p2na * p2nd # 4x

        mates[3,1] <- p1nd * p3nt # 4x
        mates[3,2] <- p2nd * p3nt # 0
        mates[3,3] <- p3nt * p3nt # 0 #Like
        mates[3,4] <- p2na * p3nt # 0

        mates[4,1] <- p1nd * p2na # 3x
        mates[4,2] <- p2nd * p2na # 4x
        mates[4,3] <- p3nt * p2na # 0
        mates[4,4] <- p2na * p2na # 4x #Like

        # Multinomial samples
        mates.vec <- abs(as.vector(mates))
        mates.vec[is.na(mates.vec)] <- 0

          if (sum(mates.vec) != 0) {
            out.t <- rmultinom(n = 1, size = abs(floor(total_poolrm/2)), prob = mates.vec)

            prm_2x <- out.t[1] ; if (is.na(prm_2x)) {prm_2x <- 0}
            prm_3x <- sum(na.omit(c(out.t[2], out.t[4], out.t[5], out.t[13]))); if (is.na(prm_3x)) {prm_3x <- 0}
            prm_4x <- sum(na.omit(c(out.t[3], out.t[6], out.t[8], out.t[9], out.t[14], out.t[16]))); if (is.na(prm_4x)) {prm_4x <- 0}

          } else{
            prm_2x <- 0
            prm_3x <- 0
            prm_4x <- 0
          }
       } else{
        prm_2x <- 0
        prm_3x <- 0
        prm_4x <- 0
      }
      outcross2x <- sum(na.omit(c(p2x_2x, prm_2x))); if (is.na(outcross2x)) {outcross2x <- 0}
      outcross3x <- sum(na.omit(c(p2x_3x, prm_3x))); if (is.na(outcross3x)) {outcross3x <- 0}
      outcross4x <- sum(na.omit(c(p2x_4x, p4x_4x, prm_4x))); if (is.na(outcross4x)) {outcross4x <- 0}
    }
    # Sum offspring
    out2x <- sum(na.omit(c(outcross2x, selfed2x))); if (is.na(out2x)) {out2x <- 0}
    out3x <- outcross3x; if (is.na(out3x)) {out3x <- 0}
    out4x <- sum(na.omit(c(outcross4x, selfed4x))); if (is.na(out4x)) {out4x <- 0}

    genlog <- c(out2x, out3x, out4x, ngd, ngt, nga)

    }
  return(genlog)
}

