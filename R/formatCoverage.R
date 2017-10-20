#' Scale a vector
#'
#' @param x
#' @param len
#'
#' @return
#' @export FALSE
#'
#' @examples
scale_exon <- function(x, len=100){
  get_clusters <- function( clusters ){
    # takes a vector produced by which(x == 1)
    cluster_start <- c()
    cluster_length <- c()
    for( i in 1:length(clusters) ){
      if(i==1){
        n_cluster <- 1
        cluster_start[n_cluster] <- clusters[i]
        cluster_length[n_cluster] <- 1
        clu <- clusters[i]
      }
      if( clusters[i] == clu + 1){
        cluster_length[n_cluster] <- cluster_length[n_cluster] + 1
        clu <- clu + 1
        next
      }
      if( clusters[i] > clu + 1){
        n_cluster <- n_cluster + 1
        cluster_start[n_cluster] <- clusters[i]
        cluster_length[n_cluster] <- 1
        clu <- clusters[i]
      }
    }
    return( list(n_clusters = n_cluster, start = cluster_start, length = cluster_length))
  }
  # initialise new empty vector
  y <- rep(0, len)
  # x is a vector of 0 and 1 encoding the positions of clusters along a sequence
  # if x is all zero then return resized zero vector
  if( length( which(x==1) ) == 0 ){
    return(y)
  }
  L <- length(x)
  clusters <- get_clusters( which(x>0) )
  # scale each cluster

  for( i in 1:clusters$n_clusters){
    s <- clusters$start[i]
    e <- clusters$start[i] + ( clusters$length[i] - 1 )
    l <- clusters$length[i]
    pos <- (s-1) / L
    S <- (pos * len) + 1
    E <- ( e / L) * len
    y[S:E] <- 1
  }

  return(y)
}

#' Perform scaling operation on a single vector
#'
#' @param cov
#' @param left
#' @param right
#' @param centre
#' @param centre_length
#'
#' @return
#' @export
#'
#' @examples
scale_coverage <- function(cov, left,right,centre, centre_length){
  # fill left and right if wanted
  leftBin <- c()
  rightBin <- c()

  if( left > 0){
    leftBin <- cov[1:left]
  }

  if( right > 0){
    rightBin <- cov[ (length(cov) - right + 1 ):length(cov)]
  }

  centreCov <- cov[ (left+1):(length(cov) - right) ]

  if( centre == "scaled"){
    centreCov <- scale_exon( centreCov, centre_length)
  }
  if( centre == "proportion" | centre == "hidden"){
    centreCov <- rep(NA, centre_length)
  }

  # TODO: add proportion mode

  outCov <- c(leftBin, centreCov, rightBin)
  return(outCov)
}


#' Create scaled and smoothed vector for plotting
#'
#' @param coverage
#' @param smoothing
#' @param left
#' @param centre
#' @param right
#' @param centre_length
#'
#' @return
#' @export TRUE
#'
#' @examples
formatCoverage <- function(coverage,
                           smoothing = 30,
                           left = 0,
                           centre = NULL,
                           right = 0,
                           centre_length=1000){
  # given a format string, format the coverage
  # example:
  # a set of introns between 100-1000bp with 50bp flank either side
  # I want to see non-scaled coverage of the interior 50bp too
  # and the remainder should be scaled to a uniform length.
  #left = 100
  #right = 100
  # centre can be one of "scaled", "proportion", NULL
  #centre_length = 200
  #centre = "hidden"
  #cov <- coverage$cov[[1]]
  scaled <- lapply(coverage$cov,
                   FUN = function(cov){
                    scale_coverage(cov, left, right, centre, centre_length)
                     })
  names(scaled) <- "placeholder"
  scaled <- as.data.frame( do.call( what = rbind, args = scaled) )

  # TODO: this should be actually be the number of rows in A
  ## the python script could output this
  interval_count <- length(coverage$cov)
  print(interval_count)
  scaled.count <- colSums(scaled,na.rm = TRUE)
  scaled.norm <- scaled.count / interval_count

  smoothing <- smoothing
  if( smoothing > 0){
    scaled.smooth <- smoother::smth(scaled.norm, window = smoothing, method = "gaussian", tails = TRUE)
    return(scaled.smooth)
  }else{
    return(scaled.norm)
  }
}
