#' CTS method for 2 sources
#'
#' Compute the CTS for the selected pairs or triplets of tracers.
#'
#' @param source Data frame containing the sediment sources from a dataset
#' @param mixture Data frame containing one of the dataset mixtures
#' @param sol Selected solution from the pairs/triplet functions
#' 
#' @return Data frame containing the CR score for each tracer.
#'
#' @export
#'
cts_2s <- function(source, mixture, sol) {
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  cols <- (ncol(source) - 1) / 2
  
  # Normalize each column for two sources
  for (col in 1:cols) {
    mx <- max(source[, col] + source[, cols + col])
    mn <- min(source[, col] - source[, cols + col])
    source[, col] <- (source[, col] - mn) / (mx - mn)
    source[, cols + col] <- source[, cols + col] / (mx - mn)
    mixture[, col] <- (mixture[, col] - mn) / (mx - mn)
  }
  
  tracer <- colnames(source)[1:cols]
  err <- numeric(cols)
  
  # Calculate errors for two sources
  for (col in 1:cols) {
    err[col] <- abs(source[1, col] * sol[1] + source[2, col] * sol[2] - mixture[1, col])
  }
  
  return(data.frame(tracer, err))
}
