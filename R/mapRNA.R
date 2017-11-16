
#' ggplot smoothed coverage as density plot
#'
#' @param coverage
#' @param fill
#' @param alpha
#'
#' @import ggplot2
#' @return
#' @export
#'
#' @examples
plotRNAmap <- function(coverage, fill, alpha,
                       bar_width = NA,
                       bar_fill = NA,
                       bar_alpha = NA){
  # extract objects from list
  scaled <- coverage$scaled
  metadata <- coverage$metadata
  print(metadata)
  # density plot
  p <- ggplot() + theme_classic() # make this customisable?

  p <- p + ggplot2::geom_area( aes(1:length(scaled), scaled ), colour = NA, fill = fill, alpha = alpha)

  # PROPORTION
  # if centre == "proportion" then add a bar plot with a label - this will be a very different proportion so label it with the percentage overlap
  # but make the height?

  if( metadata$centre_mode == "proportion"){
    centre_point <- length(scaled) / 2
    my_ymax <- max(scaled)
    if( is.na(bar_alpha) ){
      bar_alpha <- alpha
    }
    if( is.na(bar_width) ){
      bar_width <- 0.5 * metadata$centre_seg
    }
    if( is.na(bar_fill) ){
      bar_fill <- fill
    }
    p <- p +
      #geom_bar( aes( x = centre_point, y = my_ymax * 0.90 ), stat = "identity" , fill = bar_fill, alpha = bar_alpha, width = bar_width) +
      annotate( geom = "text", x = centre_point, y = my_ymax * 0.10, label = paste0( signif(coverage$overlap * 100, digits = 3), "%" ), colour = fill )
  }

  # AXES
  # X axis marks should annotate the original flanking parameters AND the parameters set by formatting
  # example: an exon may be flanked by 50bp either side but the user also wants to see the unscaled coverage of 20bp inside at either end
  # so flank = 50 but left_seg and right_seg = 70
  # the centre is then hidden and the length is then arbitrary. There should be double tick marks to denote the change in scale
  # build up x breaks

  total_length <- 1 + metadata$left_seg + metadata$centre_seg + metadata$right_seg

  # special case when the segment lengths equal the flank on both sides
  if( metadata$A_flank == metadata$left_seg & metadata$A_flank == metadata$right_seg){
    x_breaks <- c(
      1,
      metadata$A_flank,
      metadata$A_flank + metadata$centre_seg,
      metadata$A_flank + metadata$centre_seg + metadata$A_flank
    )
    x_labels <- c(
      -metadata$A_flank,
      0,
      0,
      paste("+", metadata$A_flank)
    )
  }

  # flank is smaller than both segments (symmetric)
  if( metadata$A_flank < metadata$left_seg & metadata$A_flank < metadata$right_seg){
    x_breaks <- c(
      1,
      1 + metadata$A_flank,
      1 + metadata$left_seg,
      1 + metadata$left_seg + metadata$centre_seg,
      1 + metadata$left_seg + metadata$centre_seg + metadata$right_seg - metadata$A_flank,
      metadata$left_seg + metadata$centre_seg + metadata$right_seg
     )
    x_labels <- c(
      paste0("-", metadata$A_flank),
      "5\'",
      paste0("+", metadata$left_seg - metadata$A_flank),
      paste0("-", metadata$left_seg - metadata$A_flank),
      "3\'",
      paste0("+", metadata$A_flank)
    )
    total_length <- metadata$left_seg + metadata$centre_seg + metadata$right_seg
  }

  # flank is larger than both segments (symmetric)
  if( metadata$A_flank > metadata$left_seg & metadata$A_flank > metadata$right_seg){
    # TODO
    x_breaks <- c(
      1,
      1 + metadata$left_seg,
      1 + metadata$A_flank,
      1 + metadata$A_flank + metadata$centre_seg,
      1 + metadata$A_flank + metadata$centre_seg + metadata$A_flank - metadata$right_seg,
      metadata$A_flank + metadata$centre_seg + metadata$A_flank
    )
    x_labels <- c(
      paste0("-", metadata$A_flank),
      0,
      paste0("+", metadata$left_seg - metadata$A_flank),
      paste0("-", metadata$left_seg - metadata$A_flank),
      0,
      paste0("+", metadata$A_flank)
    )
  }

  if( !exists("x_breaks")){
    message("not supported yet")
    return(NULL)
  }

  p <- p +
    scale_x_continuous(
      "",
      breaks = x_breaks,
      label = x_labels,
      limits = c(1,total_length),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      "Normalised coverage",
      #limits = c(-my_ymax, my_ymax ),
      labels = scales::percent,
      expand = c(0, 0) ) +
    theme( axis.line.x = element_line(linetype = 3))





  return(p)
}

addCoverageTrack <- function(coverage2, fill, alpha=0.5,
                              bar_width = NA,
                              bar_fill = NA,
                              bar_alpha = NA){
  density <- ggplot2::geom_area(
    aes(1:length(coverage2$scaled), coverage2$scaled ),
    colour = NA, fill = fill, alpha = alpha
    )

  return(density)
}
