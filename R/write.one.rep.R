#' Processing function to subset a single replicate from a model set.
#'
#' @description Subset a single replicate from every model set.
#'
#' @param set List containing paths to simulations wth replicates. This can be made with the function `list.files(path = "", full.names = TRUE)`.
#' @param setname Name for output folder.
#'
#' @returns `write.one.rep()` creates a folder "output/one_rep/{setname}" where a csv file for each replicate is made.


write.one.rep <- function(set, setname) {

   # Check if output directory exists
  outputdirectory <- paste0("output/one_rep/", setname)
  if(file.exists(outputdirectory) == FALSE) {
    system(paste0("mkdir ", outputdirectory))
    print(paste0(outputdirectory, "created"))
  } else {
    print(paste0(outputdirectory, "already exists :)"))
  }

  # Subset one_rep per model and write as csv
  setout <- c()
  for(j in 1:length(set)){
    name <- set[j]
    output_name <- paste0("output/one_rep/", setname, "/", name, ".csv")
    plot.pop.list <- get(set[j])
    plot.pop.one <- plot.pop.list[[1]]
    plot.pop.one$C2[is.na(plot.pop.one$C2)] <- 0
    plot.pop.one$C3[is.na(df$C3)] <- 0
    df$C4[is.na(df$C4)] <- 0
    write.csv(plot.pop.one, output_name)

  }
}
