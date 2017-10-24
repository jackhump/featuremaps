#' Create coverage object
#'
#' @param A
#' @param B
#' @param A_flank
#' @param stranded
#'
#' @return
#' @export FALSE
#'
#' @examples
createCoverage <- function(A, B, A_flank, stranded){

  path <- paste(system.file(package="RNAmaps"), "intersectRNA.py", sep="/")

  if( stranded ){
    strandFlag <- "--stranded"
  }else{
    strandFlag <- ""
  }
  command <- "python"
  # create temporary file
  outFile <- "RNAmaps_coverage_tmp"

  args <- c(A,B,outFile,A_flank, strandFlag)

  allArgs <- c(path, args)
  # this returns a vector of raw python output
  # internalise stdout of script
  metadata <- system2(command, args=allArgs,stdout = TRUE, stderr = FALSE)
  # pass as message to console
  message(paste(metadata, collapse = "\n") )

  return( list( outFile = outFile, metadata = metadata) )
}


#' Use pybedtools to intersect two bedfiles and output numerical vectors
#'
#' @param A
#' @param B
#' @param A_flank
#' @param stranded
#'
#' @return
#' @export
#'
#' @examples
intersectRNA <- function(A, B, A_flank=0, stranded=TRUE){
  # tests
  # test A and B exist
  #coverage_file <- "RNAmaps_coverage_tmp"
  coverage_object <- createCoverage(A, B, A_flank, stranded)
  coverage_raw <- readLines(coverage_object$outFile)
  # data is CSV - split
  coverage_raw <- sapply(coverage_raw, FUN = function(x) stringr::str_split(x, ","))
  # separate coordinate ids from binary coverage vector
  coverage <- list()
  # get out metadata from object
  # metadata is : delimited - LHS is key, RHS is value
  coverage_object$metadata <- coverage_object$metadata[ grepl( ":", coverage_object$metadata, fixed = TRUE)]
  coverage$metadata <- stringr::str_split_fixed( coverage_object$metadata, ":", 2)[,2]
  coverage$metadata <- stringr::str_trim(coverage$metadata, side = "left")
  names(coverage$metadata) <- stringr::str_split_fixed( coverage_object$metadata, ":", 2)[,1]

  # for each coverage line the first element is the coordinate of the A interval
  coverage$coords <- unlist(lapply(coverage_raw, FUN = function(x)
      return(x[1]) ))
  # the rest is the binary coverage vector
  coverage$cov <- lapply( coverage_raw, FUN = function(x)
    as.numeric(x[2:length(x)]
               ))
  # remove names
  names(coverage$coords) <- NULL
  names(coverage$cov) <- NULL

  return(coverage)
}







