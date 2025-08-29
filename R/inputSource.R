#' Input sediment sources
#'
#' The function select and extract the source samples of the dataset.
#'
#' @param data Data frame containing source and mixtures data
#' @param na.omit Boolean to omit or not NA values when computing the mean and SD
#'
#' @export
#' 
inputSource <- function(data, na.omit = T) {
  # reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  sources <- data[!data[,2] == levels(data[,2])[nlevels(data[,2])],]
  
  sources2 <- sources[ order(sources[,2]),]
  
  s_groups <- droplevels(sources2[,2])

  mean <- aggregate(sources2[,3:ncol(sources2)], list(s_groups), mean, na.rm = na.omit)
  mean[[1]] <- NULL
  #mean <- round(mean, 2)
  
  sd <- aggregate(sources2[,3:ncol(sources2)], list(s_groups), sd, na.rm = na.omit)
  
  colsd <- colnames(sd[,2:ncol(sd)])
  paste('D',colsd, sep='')
  sd.data <- as.data.frame(sd[,-1])
  colnames(sd.data) <- paste('D',colsd, sep='')
  #sd.data <- round(sd.data, 2)
  
  n <- data.frame(table(sources2[,2]))
  n[[1]] <- NULL
  n <- as.data.frame(n[-nrow(n),])
  colnames(n) <- 'n'
  
  id_sources <- c(1:nrow(mean))
  id <- paste('S',id_sources, sep='')
  
  sources <- cbind(id, mean, sd.data, n)
  
  return(sources)
}
