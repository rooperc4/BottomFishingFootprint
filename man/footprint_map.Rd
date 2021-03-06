% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{footprint_map}
\alias{footprint_map}
\title{Function to produce raster layers of effort by gear type and year group}
\usage{
footprint_map(footprint_data_object)
}
\arguments{
\item{footprint_data_object}{A list of data objects output from a footprint_data function}
}
\description{
This function takes a fishing effort data object (created by the footprint_data function) and makes a raster map of the effort data. The outputs are a raster 
stack with layers for each year group and gear type combination and a .png figure for each year group and gear type combination. The .png object
is overlaid on a static google map of the eastern or western Pacific, based on the input longitude values of the data. The footprint_data_object must 
be a list including a data frame of raw fishing effort data (Fishing_effort_data), spatially, temporally effort data compiled by gear type (fishing_
footprint_table), the number of grid cells where their are not enough vessels represented to plot the data due to confidentiality rules (not_enough_vessels),
gear type-year-year groups over which the data is summarized (years), gear types over which the data is summarized (gear_types), the spatial resolution 
of the data summary, #' the spatial resolution units (units) and the map projection of the data (map_projection). These 8 items are ouput from the 
footprint_data function.
#'
}
\examples{
effort_rasters<-footprint_data_object<-example_dd
}
\keyword{NPFC,}
\keyword{effort,}
\keyword{fishing}
\keyword{footprint}
