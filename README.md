# RNAmaps
An R package for creating sophisticated plots of genomic interval intersections.

## Dependencies
* R >= 3.3
  * dplyr
  * stringr
  * ggplot2
  * reshape2
  * smoother
  * GenomicRanges

## Installation

Install the package from GitHub using the devtools package.

```
library(devtools)
devtools::install_github("jackhump/RNAmaps")
```

# Creating beautiful intersection plots with *RNAmaps*

RNAmaps allows you to visualise the intersection of two sets of genomic intervals where one set of intervals (parent or *A*) overlaps the other (child or *B*). 

Examples of parent intervals include:
* introns
* exons
* entire genes
* transposable elements

Examples of child intervals include:
* iCLIP or CLIP peaks
* ChIP peaks
* Transposable elements
* locations of sequence motifs

RNAmaps first performs an intersection of the two sets of genomic intervals (stored as BED files) and records the position of each element in B that overlaps each element in A. As A intervals can vary in length, the intersected intervals are then scaled to a uniform size in a scheme set by the user which allows for unscaled sections. The relative density of at position is then calculated ( sum of the number of overlapping B features at a given position / total number of A features) to create a vector of densities.

## Worked example

In the /test directory there are two BED files, XX.bed and YY.bed. We can use RNAmaps to test for overlap and localisation of features in YY within each interval in XX

1. Reading in files and intersecting the two
``` 
intersect <- intersectRNA( A = XX.bed, B = YY.bed, A_flank = 50, B_flank = 0, stranded=TRUE)
```

2. Formatting and scaling the intersections

```
intersect <- formatRNA( left, right, centre_flank )
```

3. Plotting with ggplot2

The scaled density vector for the intersection of A and B is stored in intersect$scaled and so can be simply plot using base R:
```
plot(intersect$scaled)
```

RNAmaps was created with the ggplot2 plotting system in mind and so provides a interface to create good-looking plots.
The RNAmaps() function takes the list of arguments given to formatRNA() and creates x and y axes. The scaled coverage vector can then be plotted as an area or a line. 

```
RNAmaps( data = intersect, plotType = "area", fill = "blue", alpha = 0.5 )
```


## Visualising multiple intersections in the same plot

Often it is useful to compare the intersection of features in B with multiple sets of A intervals to look for a relative enrichment or depletion.

To save writing out lots of identical formatRNA() calls for each intersection you can set and apply a formatting scheme.

To add multiple intersections to a plot use the ggplot `+` syntax. Eg:

```
myscheme <- plotScheme( centre = "scaled", left = 0, right = 0, centre_length = 1000)
intersect1 <- intersectRNA(A = XX.bed, B = YY.bed) %>% formatRNA(scheme=myscheme)
intersect2 <- intersectRNA(A = ZZ.bed, B = YY.bed) %>% formatRNA(scheme=myscheme)

multiplot <- RNAmaps(scheme=myscheme) + 
  geom_area( makeCov(intersect1, scheme = myscheme), aes(x,y), fill = "blue", alpha = 0.5 ) +
  geom_area( makeCov(intersect2, scheme = myscheme), aes(x,y), fill = "red", alpha = 0.5 )
```

## TODO: heatmaps
## TODO: statistical tests between groups

```
a function that takes two intersect objects and computes Fisher's exact test P values at each position/range of positions
fisherP <- enrichIntersect( x = intersect1, y = intersect2, method = "fisher" )
# add the significance stars 
multiplot + plotSigTest( test = fisherP, display = "stars" ) 
```
