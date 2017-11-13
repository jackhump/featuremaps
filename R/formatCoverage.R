#' Scale a vector
#'
#' @param x
#' @param len
#'
#' @return
#' @export
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
scaleCoverage <- function(cov, left,right,centre, centre_length){
  # fill left and right if wanted
  leftBin <- c()
  rightBin <- c()

  if( left > 0){
    leftBin <- cov[1:left]
  }else{
    leftBin <- c()
  }

  if( right > 0){
    rightBin <- cov[ (length(cov) - right + 1 ):length(cov)]
  }else{
    rightBin <- c()
  }

  centreCov <- cov[ (left+1):(length(cov) - right) ]

  if( centre == "scaled"){
    centreCov <- scale_exon( centreCov, centre_length)
  }
  if( centre == "proportion" | centre == "hidden"){
    centreCov <- rep(NA, centre_length)
  }

  outCov <- c(leftBin, centreCov, rightBin)
  return(outCov)
}


#' For a single coverage vector report number of overlaps within centre
#'
#' @param cov
#' @param left
#' @param right
#' @param centre
#'
#' @return
#'
#' @examples
nOverlap <- function(cov, left,right){
  # for a single coverage vector
  # report number of overlaps within centre
  if( left > 0){
    leftBin <- cov[1:left]
  }else{
    leftBin <- c()
  }

  if( right > 0){
    rightBin <- cov[ (length(cov) - right + 1 ):length(cov)]
  }else{
    rightBin <- c()
  }

  centreCov <- cov[ (left+1):(length(cov) - right) ]

  nOverlap <- sum(centreCov)

  return(nOverlap)
}


#' Decide on how to smooth the coverage based on the left and right flanks and how the centre section is to be presented
#'
#' @param scaled.norm
#' @param smoothing
#' @param centre
#'
#' @return
#' @export
#'
#' @examples
#'
smoothCoverage <- function(scaled.norm, smoothing, left, right, centre){
  if( centre == "scaled"){
    scaled.smooth <- smoother::smth(scaled.norm, window = smoothing, method = "gaussian", tails = TRUE)
  }
  # if hidden mode smooth left and right flanks separately
  if( centre == "hidden" | centre == "proportion"){
    # flanked on both sides
    if( left > 0 & right > 0 ){
      scaled.smooth <- c(
        smoother::smth(scaled.norm[1:left], window = smoothing, method = "gaussian", tails = TRUE),
        scaled.norm[(left+1):(length(scaled.norm) - right )],
        smoother::smth(scaled.norm[ ( length(scaled.norm) - (right - 1) ):(length(scaled.norm))], window = smoothing, method = "gaussian", tails = TRUE)
      )
    }
    # flanked on right only
    if( left == 0 & right > 0){
      scaled.smooth <- c(
        scaled.norm[(left+1):(length(scaled.norm) - right )],
        smoother::smth(scaled.norm[ ( length(scaled.norm) - (right - 1) ):(length(scaled.norm))], window = smoothing, method = "gaussian", tails = TRUE)
      )
    }
    if( left > 0 & right == 0){
      scaled.smooth <- c(
        smoother::smth(scaled.norm[1:left], window = smoothing, method = "gaussian", tails = TRUE),
        scaled.norm[(left+1):(length(scaled.norm) - right )]
      )
    }
  }
  return(scaled.smooth)
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
#' @export
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

  # Illegal operations
  # if centre = "hidden" then left and/or right flank must be >0
  if( centre == "hidden" & left == 0 & right == 0 ){
    message("Hidden centre no operation to perform")
    return(NULL)
  }

  # smoothing must be smaller than the shortest nonzero flank or centre length if included (scaled)
  lengths <- c(right,left)
  if( centre == "scaled"){
    lengths <- c(lengths, centre_length)
  }
  if( centre != "scaled" & any(lengths != 0 ) ){
    if( smoothing > min(lengths[ lengths != 0]) ){
      message("Smoothing variable must be smaller than shortest nonzero segment")
      message(paste("Try setting smoothing to", ceiling(0.75 * min(lengths[lengths != 0]) ) ) )
      return(NULL)
    }
  }else{
    # if all lengths are 0
  }

  # add formatting choices to metadata
  format_metadata <- list(
    left_seg = left,
    right_seg = right,
    centre_seg = centre_length,
    centre_mode = centre,
    smoothing = smoothing
  )

  coverage$metadata <- c(coverage$metadata, format_metadata)


  scaled <- lapply(coverage$cov,
                   FUN = function(cov){
                    scaleCoverage(cov, left, right, centre, centre_length)
                     }
                   )
  names(scaled) <- "placeholder"
  scaled <- as.data.frame( do.call( what = rbind, args = scaled) )

  interval_count <- as.numeric(coverage$metadata$nA)
  scaled.count <- colSums(scaled,na.rm = TRUE)
  scaled.norm <- scaled.count / interval_count

  if( centre == "proportion"){

    overlap <- unlist(lapply( coverage$cov, FUN = function(cov){
        nOverlap(cov, left,right)
      }))
    overlap_sum <- sum(overlap > 0)
    overlap_prop <- overlap_sum / interval_count

    coverage$overlap <- overlap_prop
  }


  if( smoothing > 0){
    scaled.smooth <- smoothCoverage( scaled.norm, smoothing, left, right, centre )
    coverage$scaled <- scaled.smooth
    return(coverage)
  }else{
    coverage$scaled <- scaled.norm
    return(coverage)
  }
}
