% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omics_network.R
\docType{methods}
\name{position-methods}
\alias{position-methods}
\alias{position}
\alias{position,omics_network-method}
\title{Returns the position of edges in the network}
\usage{
\S4method{position}{omics_network}(net, nv = 0)
}
\arguments{
\item{net}{a omics_network object}

\item{nv}{the level of cutoff at which the analysis should be done}
}
\value{
Matrix with as many rows as the number of edges of the network and three
columns (name, xcoord, ycoord).
}
\description{
Returns the position of edges in the network
Retrieve network position for consistent plotting.
Utility function to plot networks.
}
\section{Methods}{
 \describe{

\item{list("signature(net = \"omics_network\")")}{ Returns a matrix with the
position of the node. This matrix can then be used as an argument in the
plot function. } }
}

\examples{

data(network)
position(network)

}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}
\keyword{dplots}
\keyword{methods}
