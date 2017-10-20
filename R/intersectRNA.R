#' Create coverage object
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
  system2(command, args=allArgs)

  return(outFile)
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
  coverage_file <- createCoverage(A, B, A_flank, stranded)
  coverage_raw <- readLines(coverage_file)
  # split out values
  coverage_raw <- sapply(coverage_raw, FUN = function(x) stringr::str_split(x, ","))
  # separate coordinate ids from binary coverage vector
  coverage <- list()
  coverage$coords <- unlist(lapply(coverage_raw, FUN = function(x) return(x[1]) ))
  coverage$cov <- lapply( coverage_raw, FUN = function(x) as.numeric(x[2:length(x)]))

  return(coverage)
}







