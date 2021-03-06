% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Metrics.R
\name{TCAI_140}
\alias{TCAI_140}
\title{Total Core Area Index with edge of 140 meters}
\usage{
TCAI_140(x)
}
\arguments{
\item{x}{Raster of data with patches identified as classes in RasterLayer format.}
}
\value{
The value of the total core area index in percentage, with subtracted edge of 140 meters.
}
\description{
Measurement of relative quantity of core area in the landscape, in percentage, when excluded the edge from 140 meters.
}
\references{
MCGARIGAL. K.; MARKS, B. J. Fragstats: Spatial pattern analysis program for quantifying landscape structure. Reference manual. Forest Science Department Oregon State University. Corvallis Oregon 1994. 60 p.
}
