% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dscan.R
\name{Dbscan.cluster}
\alias{Dbscan.cluster}
\title{Density clustering}
\usage{
Dbscan.cluster(type.pos, eps = NULL, MinPts)
}
\arguments{
\item{type.pos}{data frame contain corordinate of spots}

\item{eps}{distance measure that will be used to locate the points in the neighborhood of any point}

\item{MinPts}{the minimum number of points (a threshold) clustered together for a region to be considered dense}
}
\value{
a data frame containing the cell coordinates and cluster id of clustering
}
\description{
Density clustering using the coordinates of spots using dbscan method.
}
\examples{
fx <- Dbscan.cluster(type.pos = as.data.frame(matrix(rnorm(100), ncol = 2)), eps = 3, MinPts = 5)
}
