#' Format multiple generation simulation.
#'
#' @description This function formats the output of gen-iter-f-choosy.
#'
#'
#' @param popvect Output from gen-iter-f-choosy function.
#' @param generations Number of generations included.
#'
#' @returns Output file is a data frame where each row is a generation. The columns are as follows,
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

form.autopop <- function(popvect, generations){
    plot.pop <- as.data.frame(do.call(rbind, popvect))
    plot.pop$gen <- 1:(generations + 1)
    plot.pop[is.na(plot.pop)] = 0

    ## Total pop size
    for(f in 1:(generations + 1)){
      plot.pop$sum[f] <- sum(plot.pop[f, 1:6])
      plot.pop$sum2x[f] <- (sum(plot.pop$V1[f], plot.pop$V4[f]))
      plot.pop$sum3x[f] <- (sum(plot.pop$V2[f], plot.pop$V5[f]))
      plot.pop$sum4x[f] <- (sum(plot.pop$V3[f], plot.pop$V6[f]))
    }

    for(i in 1:(generations + 1)){

      plot.pop$V1a[i] <- (plot.pop$V1[i]/plot.pop$sum[i])
      plot.pop$V2a[i] <- (plot.pop$V2[i]/plot.pop$sum[i])
      plot.pop$V3a[i] <- (plot.pop$V3[i]/plot.pop$sum[i])
      plot.pop$V4a[i] <- (plot.pop$V4[i]/plot.pop$sum[i])
      plot.pop$V5a[i] <- (plot.pop$V5[i]/plot.pop$sum[i])
      plot.pop$V6a[i] <- (plot.pop$V6[i]/plot.pop$sum[i])
      plot.pop$C2[i] <- (sum(plot.pop$V1[i], plot.pop$V4[i])/plot.pop$sum[i])
      plot.pop$C3[i] <- (sum(plot.pop$V2[i], plot.pop$V5[i])/plot.pop$sum[i])
      plot.pop$C4[i] <- (sum(plot.pop$V3[i], plot.pop$V6[i])/plot.pop$sum[i])
    }
    return(plot.pop)

}

