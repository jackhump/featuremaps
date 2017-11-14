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
createCoverage <- function(A,
                           B,
                           A_flank,
                           stranded,
                           outFile=NA
                           ){

  path <- paste(system.file(package="RNAmaps"), "intersectRNA.py", sep="/")

  if( stranded ){
    strandFlag <- "--stranded"
  }else{
    strandFlag <- ""
  }
  command <- "python"
  # create temporary file
  if( is.na(outFile)){
    # generate random temp file name
    outFile <-"RNAmaps_coverage_tmp"
  }else{
    outFile <- outFile
  }
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
intersectRNA_python <- function(A, B, A_flank=0, stranded=TRUE){
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
  # convert to list
  coverage$metadata <- as.list(coverage$metadata)
  coverage$metadata$nA <- as.numeric(coverage$metadata$nA)
  coverage$metadata$nB <- as.numeric(coverage$metadata$nB)
  coverage$metadata$flank <- as.numeric(coverage$metadata$flank)

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



#' Intersect two bed files - using GenomicRanges in R
#'
#' @param A
#' @param B
#' @param stranded
#' @param A_flank
#' @param B_flank
#'
#' @return
#' @export
#'
#' @examples
intersectRNA <- function(A,B, stranded = TRUE, A_flank = NA, B_flank = NA){

  if(stranded){
    ignoreStrand <- FALSE
  }else{
    ignoreStrand <- TRUE
  }

  # read in bed files to genomic ranges
  a <- rtracklayer::import.bed(A)
  b <- rtracklayer::import.bed(B)

  # if flank requested then flank datasets
  # assume equal flank at each end
  if( !is.na(A_flank)){
    GenomicRanges::start(a) <-  GenomicRanges::start(a) - A_flank
    GenomicRanges::end(a) <-  GenomicRanges::end(a) + A_flank
  }
  if( !is.na(B_flank)){
    GenomicRanges::start(b) <-  GenomicRanges::start(b) - B_flank
    GenomicRanges::end(b) <-  GenomicRanges::end(b) - B_flank
  }

  # set seqlevels
  all_chrs <- intersect( seqlevels(a), seqlevels(b) )
  GenomeInfoDb::seqlevels(a, force = TRUE) <- all_chrs
  GenomeInfoDb::seqlevels(b, force = TRUE) <- all_chrs

  # remove duplicates from a
  a <- IRanges::reduce(a)

  # perform intersect
  intersect <- GenomicRanges::findOverlaps(a,b, ignore.strand = ignoreStrand)

  # just elements in A that overlap B
  a_list <- a[ unique(S4Vectors::queryHits(intersect)) ]

  # create list where each element is a
  # genomic ranges object of the elements in B that overlap one entry in A
  b_in_a <- GenomicRanges::intersect(a,b, ignore.strand = ignoreStrand)
  # convert to list - this is very time consuming!
  b_list <- split(b_in_a, S4Vectors::queryHits(intersect) )

  # for each element in a_list
  fillVector <- function(i){

    a_range <- a_list[i] # time consuming!
    a_start <- GenomicRanges::start(a_range)
    a_width <- GenomicRanges::width(a_range)


    strand <- as.character( GenomicRanges::strand( a_range ) )
    # create zero vector - this is time consuming for large vectors
    vector <- rep.int(0, times = a_width )

    b_range <- GenomicRanges::ranges(b_list[[i]]) # time consuming!
    b_start <- GenomicRanges::start(b_range)
    b_end <- GenomicRanges::end(b_range)
    # subtract b coordinates by a to get relative coordinates
    b_start <- b_start - a_start
    b_end <- b_end - a_start

    # walk through b_list and convert 0 to 1 at those coordinates
    # create list of start and end relative coordinates
    coords_in_b <- unlist(purrr::map2( b_start, b_end, seq) )

    # remove any that overhang
    overhangs <- which( coords_in_b > a_width | coords_in_b <  0  )
    if( length(overhangs) > 0){
      coords_in_b <- coords_in_b[ -overhangs ]
    }
    # set all elements within those ranges to 1
    vector[ coords_in_b ] <- 1
    # flip vector if on negative strand

    if( strand == "-"){
      vector <- rev.default(vector)
    }

    #message(i)
    return(vector)
  }
  #library(microbenchmark)
  #microbenchmark(fillVector(200), times =50)

  vectors <- purrr::map( seq_along(a_list), ~fillVector(.) )

  return(vectors)
}





