# OLD CODE

# Entire intron coverage created by whole_intron_coverage.py
# new method - using coverage across the entire intron
library("smoother")
library(ggplot2)
library(stringr)
library(optparse)
library(reshape2)
# Jack Humphrey
# parameters
smoothing = 30
flank = 300
exon_flank = 15
intron_length = 100
exon_length = 100 # has to be this due to flanking at each end
exon_height = 0.1
bar_width = 0.7 * ( exon_length - (2*exon_flank) )
# This needs to be done three times:
# Skipped exons
# Included exons
# Control exons
options(echo=TRUE)

opt <- parse_args(
  OptionParser(option_list=list(
    make_option( "--included", type="character", default=NULL, "coverage over the introns encompassing included exons"),
    make_option( "--skipped",type="character", default=NULL, "coverage over the introns encompassing skipped exons"),
    make_option( "--control",type="character", default=NULL, "coverage over a set of introns containing non-regulated exons"),
    make_option( "--skipped_exons",type="character", default=NULL, "a bed file of the skipped exons"),
    make_option( "--included_exons",type="character", default=NULL, "a bed file of the included exons"),
    make_option( "--control_exons",type="character", default=NULL, "a bed file of the control exons"),
    make_option( "--code",type="character", default=NULL, help = "the same dataset-specific code used throughout the pipeline"),
    make_option( "--no_scaled_intron", action="store_true",type = "logical", dest=NULL, default=FALSE, "whether to plot the scaled coverage"),
    make_option( "--no_scaled_exon", action="store_true", type = "logical", dest= NULL, default =FALSE, " where to included scaled exon coverage"),
    make_option( "--mode", type="character", default=NULL, "whether coverage represents motifs or iCLIP clusters"),
    make_option( "--input", type="character", "the clusters or motif used" ),
    #make_option( "--use_strand", action="store_true", type = "logical", default = FALSE, help = "whether to use the strand from the exon list to flip the coverage of negative strand introns. Not necessary for motifs as already flipped"),
    make_option( "--outFolder",type="character", default=NULL, help = "where you want the plot to go"),
    make_option( "--annotations", type="character", "biomart annotations") )
  )
)



outFolder <- opt$outFolder
included <- opt$included
skipped <- opt$skipped
control <- opt$control
skipped_exons <- opt$skipped_exons
included_exons <- opt$included_exons
control_exons <- opt$control_exons
code <- opt$code
no_scaled_intron <- opt$no_scaled_intron
no_scaled_exon <- opt$no_scaled_exon
mode <- opt$mode
input <- opt$input
biomart <- opt$annotations
#use_strand <- opt$use_strand

if( is.null(no_scaled_intron)){
  no_scaled_intron <- FALSE
}
if( is.null(no_scaled_exon)){
  no_scaled_exon <- FALSE
}
#if( is.null(use_strand) ){
#  use_strand <- FALSE
#}

#cat("use strand? ", use_strand)


exists <- 0
for( file in c(skipped,included,control,skipped_exons,included_exons,control_exons)){
  #print(file)
  if(!file.exists(file)){
    print(paste0(file," does not exist"))
    exists <- exists+ 1
  }
}

stopifnot(exists == 0)


# FUNCTIONS

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
  #print(clusters)
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

create_intron_df <- function(coverage){
  # create data frame of introns in the coverage list
  introns <- unlist(lapply(coverage, FUN = function(x) return(x[[1]]) ))
  names(introns) <- NULL
  split <- str_split_fixed(introns, ":", 2)
  intron_df <- data.frame(
    chr = split[,1],
    start = as.numeric( str_split_fixed(split[,2], "-", 2 )[,1] ),
    end = as.numeric( str_split_fixed(split[,2], "-", 2)[,2] ),
    stringsAsFactors = FALSE
  )
  return(intron_df)
}

create_matched_exons <- function(exons, intron_df){
  # match exons with introns - pick one exon per intron
  matched_exons <- apply( intron_df, MAR = 1, FUN = function(x){
    exon <-  exons[exons$V1 == x[[1]] & exons$V2 > as.numeric( x[[2]]) & exons$V3 < as.numeric(x[[3]]) ,]
    if( nrow(exon) > 1){
      exon <- exon[1,]
    }
    if( nrow(exon) == 0 ){
      exon <- exons[1,]
      exon[,1:ncol(exon)] <- NA
    }
    return(exon)
  } )
  matched_exons <- do.call(rbind, args = matched_exons)
  names(matched_exons)[1:3] <- c("chr", "start", "end")

  return(matched_exons)
}

smooth_scaled_coverage <- function( coverage, intron_df, matched_exons, num.exons ){
  coverage <- lapply( coverage, FUN = function(x) as.numeric(x[2:length(x)]))

  scaled_coverage <- lapply(1:length(coverage), FUN = function(i) {
    #print(i)
    cov <- coverage[[i]]
    #print(sum(cov))
    exon <- matched_exons[i,]
    intron <- intron_df[i,]
    # check for sanity
    #print("exon length:")
    #print(exon$end - exon$start)
    total_length <- exon_length + flank + intron_length + flank + exon_length + flank + intron_length + flank + exon_length
    if( is.na(exon$chr) ){
      print("exon is NA")
      return( rep(0, total_length ))
    }
    #stopifnot( length(cov) == intron$end - intron$start) # sanity check - motif sequence 1nt shorter than intron interval
    stopifnot( exon$chr == intron$chr & intron$start < exon$start & intron$end > exon$end) # double check the matching
    l <- length(cov)
    S <- intron$start
    E <- intron$end
    #if( use_strand == TRUE ){
    strand <- exon$V6
    # }else{
    #   strand <- "+"
    # }
    #print(strand)
    if( strand == "+"){
      s <- exon$start - S
      e <- exon$end - S
    }else{
      s <- E - exon$end
      e <- E - exon$start
    }
    #print(cov[(e-5):(e+5)])
    # extract raw coverage
    upstream_exon_with_flank <- cov[1:(exon_length+flank)]
    # danger of flank being larger than the distance between the start of the exon and the start of the intron
    if( s > flank){
      #print("threeSS looks good")
      threeSS_flank <- cov[(s-flank):(s-1)]
    }else{
      print("dodgy threeSS")
      threeSS_flank <- rep(0,flank)
      threeSS_flank[(flank-(s-2) ):flank] <- cov[1:(s-1)]
    }
    # exon 62 is partly NA for this
    if( e+flank < l){

      fiveSS_flank <- cov[(e):(e+flank)]
      #print(fiveSS_flank)
    }else{
      print("dodgy fiveSS")
      fiveSS_flank <- rep(0,flank)
      #fiveSS_flank
    }
    if( e + flank > l & s < flank ){
      downstream_exon_with_flank <- rep(0, flank + exon_length)
    }else{
      downstream_exon_with_flank <- cov[ (l - (flank + exon_length - 1) ):l]
    }
    # extract scaled coverage
    # introns may not exist due to overlapping flanking areas so in that case add fake coverage
    if( (1+exon_length+flank) < (s-flank-1) ){
      upstream_intron <- scale_exon( cov[(1+exon_length+flank):(s-flank-1)], intron_length)
    }else{
      upstream_intron <- rep(0,intron_length)
      # what to put in here?
    }
    if( (e+flank) < (l-flank)){
      downstream_intron <- scale_exon( cov[(e+flank):(l-flank)], intron_length)
    }else{
      downstream_intron <- rep(0,intron_length)
    }

    if(no_scaled_intron){
      upstream_intron <- rep(0,intron_length)
      downstream_intron <- rep(0,intron_length)
    }

    # central exon
    upstream_flank <- cov[s:(s+exon_flank-1)]
    downstream_flank <- cov[ (e-exon_flank+1):e ]
    exon <- scale_exon( cov[ (s+exon_flank):(e-exon_flank)], exon_length - (2*exon_flank) )

    exon_sum <- sum( cov[s:e] )

    if(no_scaled_exon){
      exon <- rep(0, exon_length - (2*exon_flank) )
    }
    # paste together all the parts
    total <- c(
      upstream_exon_with_flank,
      upstream_intron,
      threeSS_flank,
      upstream_flank,
      exon,
      downstream_flank,
      fiveSS_flank,
      downstream_intron,
      downstream_exon_with_flank
    )
    return( list(all_coverage = total, exon_sum = exon_sum) )
  } )

  #save(scaled_coverage, file = "/Users/Jack/google_drive/TDP_paper/RNA_maps/iCLIP/scripts/test.Rdata")

  all_coverage <- lapply(scaled_coverage, FUN = function(x){ return(x[[1]]) })

  scaled <- as.data.frame( do.call( what = rbind, args = all_coverage) )
  # how many exons have overlapping coverage?

  exon_count <- unlist(lapply(scaled_coverage, FUN = function(x){ return(x[[2]]) } ))
  # divide by total number of exons - not just the ones that have any coverage at all.
  exon_prop <- sum( exon_count > 0  ) / num.exons
  print(exon_prop)

  #return(scaled
  scaled.count <- colSums(scaled,na.rm = TRUE)
  scaled.norm <- scaled.count / num.exons
  #return(scaled.norm)
  smoothing <- 30
  scaled.smooth <- smth(scaled.norm, window = smoothing, method = "gaussian", tails = TRUE)

  return( list(smoothed = scaled.smooth, matrix = scaled, count = exon_count, proportion = exon_prop) )
}

# run loop three times
create_coverage <- function(coverage_file, exons_bed){
  coverage <- readLines(coverage_file)
  coverage <- sapply(coverage, FUN = function(x) str_split(x, ","))

  intron_df <- create_intron_df(coverage)

  exons <- read.table(exons_bed, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  matched_exons <- create_matched_exons(exons, intron_df)
  num.exons <- nrow(exons)

  smooth_coverage <- smooth_scaled_coverage(coverage, intron_df, matched_exons, num.exons )
  return( list( cov = smooth_coverage, exons = matched_exons) )
}

get_exon_number <- function(exons_bed){
  exons <- read.table(exons_bed, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  num_exons <- nrow(exons)
  return(num_exons)
}



# create numeric coverage vectors
#  coverage_file <- "/Users/Jack/Google Drive/TDP_paper/RNA_maps/results/F210I_included_100_coverage.csv"
#  exons_bed <- "/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data//noheader/F210I_embryonic_brain_se_intron_included.bed"
#  test <- create_coverage( coverage_file, exons_bed)
# #


# run each set through
cov_list <- list() # smooth coverage
matrix_list <- list() # matrix of counts
exon_list <- list() # a list of exons
exon_num_list <- list() # exon number
prop_exon_list <- list() # proportion of overlapping exons
count_exon_list <- list()

if( !is.null(included) & !is.null(included_exons) ){
  cov_object <- create_coverage(included,included_exons)
  cov_list[[1]] <- cov_object$cov$smoothed
  matrix_list[[1]] <- cov_object$cov$matrix
  exon_list[[1]] <- cov_object$exons
  prop_exon_list[[1]] <- cov_object$cov$proportion
  count_exon_list[[1]] <- cov_object$cov$count
  exon_num_list[[1]] <- get_exon_number(included_exons)
}else{
  cov_list[[1]] <- NA
  exon_num_list[[1]] <- 0
}

if( !is.null(skipped) & !is.null(skipped_exons) ){
  cov_object <- create_coverage(skipped,skipped_exons)
  cov_list[[2]] <- cov_object$cov$smoothed
  matrix_list[[2]] <- cov_object$cov$matrix
  exon_list[[2]] <- cov_object$exons
  prop_exon_list[[2]] <- cov_object$cov$proportion
  count_exon_list[[2]] <- cov_object$cov$count
  exon_num_list[[2]] <- get_exon_number(skipped_exons)
}else{
  cov_list[[2]] <- NA
  exon_num_list[[2]] <- 0
}

if( !is.null(control) & !is.null(control_exons) ){
  cov_object <- create_coverage(control,control_exons)
  cov_list[[3]] <- cov_object$cov$smoothed
  # always way too big - don't add
  #matrix_list[[3]] <- cov_object$cov$matrix
  count_exon_list[[3]] <- cov_object$cov$count
  prop_exon_list[[3]] <- cov_object$cov$proportion
  exon_list[[3]] <- cov_object$exons
  exon_num_list[[3]] <- get_exon_number(control_exons)
}else{
  cov_list[[3]] <- NA
  exon_num_list[[3]] <- 0
}




#results <- paste0(outFolder,"/", code, "_", "coverage_data.Rdata")
results <- paste0(outFolder, "/", gsub(" ", "_", code),"_", mode, "_", input,"_coverage_data.Rdata")
save.image(results)

# for debugging
coverage_file <- included
exons_bed <- included_exons


# quit()
#
#load("/Users/Jack/Google Drive/TDP_paper/RNA_maps/results/coverage_data.Rdata")
#

# plotting
my_ymax <- max(unlist(lapply(cov_list, max)),na.rm = TRUE )
my_ymax <- my_ymax + 0.2 * my_ymax

print(my_ymax)

# hardcode ymax
#my_ymax <- 0.2

x_breaks <- c( 1+(exon_length/2),
               1+exon_length,
               1+exon_length+flank,
               1+exon_length+flank+intron_length,
               1+exon_length+flank+intron_length+flank,
               1+exon_length+flank+intron_length+flank+exon_flank,
               1+exon_length+flank+intron_length+flank+(exon_length/2), # midpoint!
               1+exon_length+flank+intron_length+flank+exon_length-exon_flank,
               1+exon_length+flank+intron_length+flank+exon_length,
               1+exon_length+flank+intron_length+flank+exon_length+flank,
               1+exon_length+flank+intron_length+flank+exon_length+flank+intron_length,
               1+exon_length+flank+intron_length+flank+exon_length+flank+intron_length+flank,
               1+exon_length+flank+intron_length+flank+exon_length+flank+intron_length+flank + (exon_length/2)
)


x_labels <- c("upstream\nexon", "",flank, -flank, "", exon_flank,"\ncentral\nexon",-exon_flank,  "", flank, -flank, "", "downstream\nexon" )
total_length <- exon_length + flank + intron_length + flank + exon_length + flank + intron_length + flank + exon_length
centre_point <- 1+exon_length+flank+intron_length+flank+(exon_length/2)
# draw box on plot to represent exon

mySubTitle <- paste0( "Included exons: ", exon_num_list[[1]], "; Skipped exons: ", exon_num_list[[2]], "; Control exons: ", exon_num_list[[3]] )

exon_df <- data.frame(
  x = c(
    x_breaks[5],
    x_breaks[9],
    x_breaks[9],
    x_breaks[5],

    1,
    x_breaks[2],
    x_breaks[2],
    1,

    x_breaks[12],
    x_breaks[12] + exon_length -1,
    x_breaks[12] + exon_length -1,
    x_breaks[12]

  ),
  y = rep( c( exon_height * my_ymax, exon_height * my_ymax, -exon_height * my_ymax, -exon_height * my_ymax), 3 ),
  group = c(
    1,1,1,1,
    2,2,2,2,
    3,3,3,3
  )
)
included_df <- data.frame(
  x = c(
    x_breaks[2],
    x_breaks[9] ),
  xend = c(
    x_breaks[5],
    x_breaks[12]
  )
)
skipped_df <- data.frame(
  x = x_breaks[2],
  xend = x_breaks[12]
)

p <- ggplot() + theme_classic()
# included exons
if( !is.na(cov_list[[1]])){
  p <- p + geom_area( aes(1:length(cov_list[[1]]), cov_list[[1]] ), color = NA, fill = "blue", alpha = 0.9)
}

# skipped exons
if( !is.na( cov_list[[2]])){
  p <- p + geom_area( aes(1:length(cov_list[[2]]), -1*cov_list[[2]] ), color = NA, fill = "red", alpha = 0.9)
}

# control exons
if( !is.na(cov_list[[3]]) & !is.na( cov_list[[1]]) ){
  p <- p + geom_area( aes(1:length(cov_list[[3]]), cov_list[[3]] ), color = NA, fill = "gray", alpha = 0.75 )
}

# control exons
if( !is.na(cov_list[[3]]) & !is.na( cov_list[[2]] ) ){
  p <- p + geom_area( aes(1:length(cov_list[[3]]), -1*cov_list[[3]] ), color = NA, fill = "gray", alpha = 0.75 )
}

if( mode == "motif"){
  title <- paste0( gsub("_", " ", code), " motif RNA map" )
}else{
  title <- paste0(code, " iCLIP cluster RNA map")
}

p <- p +
  scale_x_continuous(
    "",
    breaks = x_breaks,
    label = x_labels,
    limits = c(0,total_length+1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    "Per-nucleotide normalised coverage",
    limits = c(-my_ymax, my_ymax ),
    labels = scales::percent,
    expand = c(0, 0),
    sec.axis = sec_axis(~./my_ymax*0.5,
                        name = "Percentage of exons covered",
                        labels = scales::percent)
  ) +
  labs(title = title, subtitle = mySubTitle ) +
  # # add exon
  geom_polygon(data=exon_df, aes(x,y, group = group), fill = NA, colour = "black", linetype="dashed", size = 1.1) +

  # add intron lines
  geom_curve( data=included_df, aes(x=x,xend=xend, y=0,yend=0), curvature = -0.2, size = 1.1 ) +
  geom_curve( data=skipped_df, aes(x=x, xend=xend, y=0,yend=0), curvature = 0.2, size = 1.1 ) +
  # add dashed lines for locations
  geom_vline( data = data.frame(xintercept = x_breaks[c(2:6,8:12)], slope=0), aes(xintercept=xintercept),colour = "black", linetype="dashed",alpha=0.4 )

if( !is.null( prop_exon_list[[1]]) ){
  p <- p + geom_bar( aes( x = centre_point, y = prop_exon_list[[1]] * my_ymax * 2), stat = "identity" , fill = "blue", alpha = 1.0, width = bar_width)
}

if( !is.null(prop_exon_list[[3]]) & !is.null( prop_exon_list[[1]] ) ){
  p <- p + geom_bar( aes( x = centre_point, y = prop_exon_list[[3]] * my_ymax * 2), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width )
}

if( !is.null(prop_exon_list[[2]])){
  p <- p + geom_bar( aes( x = centre_point, y = -prop_exon_list[[2]] * my_ymax * 2), stat = "identity" , fill = "red", alpha = 1.0, width = bar_width )
}

if( !is.null(prop_exon_list[[3]]) & !is.null(prop_exon_list[[2]])){
  p <- p + geom_bar( aes( x = centre_point, y = -prop_exon_list[[3]] * my_ymax * 2), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width )
}



fileName <- paste0(outFolder, "/", gsub(" ", "_", code),"_", mode, "_", input,"_RNAMap_whole_intron.pdf")

cat("saving plot as: ", fileName,'\n')

ggsave(p,
       units = "in",
       width = 20,
       height = 5,
       filename = fileName
)

# making heatmaps!

library(stringr)
library(reshape2)

options(echo=TRUE)



#data <- opt$data

n <- 20
# load in Rdata from iCLIP and UG RNA maps
annotations <- read.table(biomart,header=TRUE, stringsAsFactors = FALSE)
# remember - the two coverage data are not in matching row order - need to use exon list to match correctly.
# BUT! there are NA values and duplicates in the exon list - how stupid!

# make heatmaps of the most decorated genes

# for testing!
mat <- matrix_list[[1]]
exons <- exon_list[[1]]
matching <- NULL
n <- 20
counts <- count_exon_list[[1]]

prepareHeatMap <- function(mat, exons, counts, n, matching=NULL, centre_point ){
  # exon list is in same order as coverage matrix so a straightforward merge is fine
  row.names(exons) <- 1:nrow(exons)
  # remove duplicates - by gene name
  duplicates <- duplicated(exons$V4)
  exons <- exons[!duplicates,]
  mat <- mat[!duplicates,]


  counts <- counts[!duplicates]
  counts <- data.frame(x = counts)

  # replace EnsemblIDs with gene symbols
  if( any(grepl("ENS.*G.*0",exons$V4)) ){
    gene_symbol <- annotations$external_gene_name[ match(exons$V4, annotations$EnsemblID)]
    # fix compound names
    if( nrow( exons[is.na(gene_symbol),] ) > 0 ){
      fixed_names <- apply( exons[is.na(gene_symbol),], MAR = 1, FUN = function(x){
        genes <- str_split( x[4], "\\+", 2)[[1]]
        symbols <- annotations$external_gene_name[ match(genes, annotations$EnsemblID)]
        symbols <- paste(symbols, collapse="+")
        return(symbols)
      })
      gene_symbol[ is.na(gene_symbol)] <- as.character(fixed_names)
    }
    exons$V4 <- gene_symbol
  }

  # either find the top genes by CLIP cluster number or match a set of genes
  if( is.null(matching)){

    #exons$V4[ match( matching, exons$V4 ) ]

    # find the top genes and keep their exons
    mat <- mat[ order(rowSums(mat), decreasing = TRUE), ]
    mat <- mat[1:n,] # pick top

    counts <- counts$x[ match( row.names(mat), row.names(counts))]

    genes <- exons$V4[ match(row.names(mat), row.names(exons))]
    row.names(mat) <- genes

    counts <- data.frame( gene = genes, N = counts)


  }else{
    # get the rows of the matrix matching the row numbers of the list of genes to match
    mat <- mat[ match( matching, exons$V4 ) , ]
    x <- counts$x[ match( row.names(mat), row.names(counts))]

    row.names(mat) <- matching
    counts <- data.frame(gene = matching, N = x)
  }
  #pvalues <- signif(exons$V5[ match(row.names(mat), row.names(exons))],digits=2)
  names(mat) <- 1:ncol(mat)

  mat.melt <- melt(as.matrix(mat))
  small <- subset(mat.melt, value > 0)
  small$Var1 <- factor(small$Var1, levels = rev(row.names(mat)) )
  #print(exons)
  print(counts)
  counts$x <- centre_point
  counts <- subset(counts, N > 0)

  return(list(melted=small, genes = genes, counts = counts) )
}


makePlot <- function(data, motifs = NULL, colour, direction){
  title <- paste(code, mode, input, direction)
  p <- ggplot(data$melted, aes(x = Var2, y = Var1 ) ) +
    geom_tile(fill = colour, colour = colour, alpha = 1.0 ) +
    scale_y_discrete("")
  if( !is.null(motifs) ){
    p <- p +
      geom_tile(data=motifs$melted, aes(x = Var2, y = Var1 ), fill = "gray", colour = "gray", alpha = 1  ) +
      geom_point(data=motifs$counts, aes(x = x, y = gene ), colour = "gray", alpha = 1)
  }
  p <- p +

    scale_x_continuous(
      "",
      breaks = x_breaks,
      label = x_labels,
      limits = c(0,total_length+1),
      expand = c(0, 0)
    ) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(colour = "gray")
    ) +
    geom_vline( data = data.frame(xintercept = x_breaks[c(2:6,8:12)], slope=0),
                aes(xintercept=xintercept),colour = "black", linetype="dashed",alpha=0.5 ) +
    geom_polygon(data=blank_df, aes(x,y, group = group), fill = "white", colour = "NA") +
    geom_point(data=data$counts, aes(x = x, y = gene ), colour = colour, alpha = 0.8) +
    ggtitle(title)
  return(p)
}


#load("../iCLIP/results/F210I_coverage_data.Rdata")
#data <- "/Users/Jack/google_drive/TDP_paper/RNA_maps/iCLIP/results//M323K_iCLIP_F210I_WT_coverage_data.Rdata"
#load(data)
#print(centre_point)

#load("TGT/F210I_TGT_coverage_data.Rdata")
#included_motifs <- prepareHeatMap(matrix_list[[1]], exons = exon_list[[1]], matching = clip$genes, n = n, counts = count_exon_list[[1]], centre_point = centre_point)
#skipped_motifs <- prepareHeatMap(matrix_list[[2]], exons = exon_list[[2]], matching = clip$genes, n = n, counts = count_exon_list[[2]], centre_point = centre_point)
# fudge - x_breaks are loaded in with the data
blank_df <- data.frame(
  x = c(
    x_breaks[3],
    x_breaks[3],
    x_breaks[4],
    x_breaks[4],

    # x_breaks[6],
    # x_breaks[6],
    # x_breaks[8],
    # x_breaks[8],

    x_breaks[10],
    x_breaks[10],
    x_breaks[11],
    x_breaks[11]
  ),
  y = rep(c(n+1,0,0,n+1),2),
  group = c(
    1,1,1,1,
    # 2,2,2,2,
    3,3,3,3
  )
)


if( !is.null(exon_list[[1]])){
  included_iCLIP <- prepareHeatMap(matrix_list[[1]], exons = exon_list[[1]], n, counts = count_exon_list[[1]], centre_point = centre_point)
  p1 <- makePlot(included_iCLIP, "blue", motifs = NULL, direction = "included cassette exons")
}

if( !is.null(exon_list[[2]]) ){
  skipped_iCLIP <- prepareHeatMap(matrix_list[[2]], exons = exon_list[[2]], n, counts = count_exon_list[[2]], centre_point = centre_point)
  p2 <- makePlot(skipped_iCLIP, "red", motifs = NULL, direction = "skipped cassette exons")
}

#p3 <- makePlot(included_iCLIP, "red", motifs = included_motifs)
#p4 <- makePlot(skipped_iCLIP, "blue", motifs = skipped_motifs)
fileName <- paste0(outFolder, "/", gsub(" ", "_", code),"_", mode, "_", input,"_heatmap.pdf")

pdf(fileName,width = 20,height=5)

#grid.arrange(p1,p2, nrow=2)

if( !is.null(exon_list[[1]])){
  print(p1)
}
if( !is.null(exon_list[[2]]) ){
  print(p2)
}
#print(p3)
#print(p4)
dev.off()

# save new version if you get this far.
results <- paste0(outFolder, "/", gsub(" ", "_", code),"_", mode, "_", input,"_coverage_data.Rdata")
save.image(results)



